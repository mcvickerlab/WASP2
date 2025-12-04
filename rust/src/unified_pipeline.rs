//! Unified Pipeline - Single-pass BAM processing for WASP2
//!
//! Replaces the multi-pass pipeline (filter + intersect + remap) with a single
//! BAM read that streams directly to FASTQ output.
//!
//! # Performance Target
//! - Current: ~500s (400s filter + 24s intersect + 76s remap)
//! - Target: ~100s (single pass)
//!
//! # Memory Budget
//! - VariantStore: ~250MB (2M variants)
//! - Pair buffer: ~1GB peak (500K pairs × 2KB)
//! - Channel buffers: ~20MB
//! - Total: ~1.3GB

use anyhow::{Context, Result};
use coitrees::{IntervalTree, SortedQuerent};
use crossbeam_channel::{bounded, Receiver, Sender};
use flate2::write::GzEncoder;
use flate2::Compression;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::{bam, bam::Read as BamRead};
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::thread;
use std::time::Instant;

use crate::bam_intersect::{build_variant_store, VariantStore};
use crate::bam_remapper::{generate_haplotype_seqs, RemapConfig, VariantSpan as RemapVariantSpan};
use crate::cigar_utils;

// ============================================================================
// Configuration and Statistics
// ============================================================================

/// Configuration for unified pipeline
#[derive(Debug, Clone)]
pub struct UnifiedConfig {
    /// Number of BAM reading threads
    pub read_threads: usize,
    /// Maximum haplotype sequences per read pair
    pub max_seqs: usize,
    /// Bounded channel buffer size
    pub channel_buffer: usize,
}

impl Default for UnifiedConfig {
    fn default() -> Self {
        Self {
            read_threads: 8,
            max_seqs: 64,
            channel_buffer: 50_000,
        }
    }
}

/// Statistics returned from unified pipeline
#[derive(Debug, Default, Clone)]
pub struct UnifiedStats {
    /// Total reads processed
    pub total_reads: usize,
    /// Read pairs processed
    pub pairs_processed: usize,
    /// Pairs with at least one variant overlap
    pub pairs_with_variants: usize,
    /// Total haplotype reads written
    pub haplotypes_written: usize,
    /// Pairs kept (no variants at all)
    pub pairs_kept: usize,
    /// Pairs skipped because minimum-position variant is in intron/deletion
    /// This matches baseline behavior where such pairs are discarded
    pub pairs_skipped_unmappable: usize,
    /// Pairs where haplotype generation failed (should be rare)
    pub pairs_haplotype_failed: usize,
    /// Orphan reads (mate not found)
    pub orphan_reads: usize,
    /// Time spent building variant tree (ms)
    pub tree_build_ms: u64,
    /// Time spent streaming BAM (ms)
    pub bam_stream_ms: u64,
}

// ============================================================================
// Haplotype Output Structure
// ============================================================================

/// A haplotype read ready for FASTQ output
#[derive(Debug, Clone)]
pub struct HaplotypeOutput {
    /// Read name with WASP suffix
    pub name: Vec<u8>,
    /// Sequence with swapped alleles
    pub sequence: Vec<u8>,
    /// Quality scores
    pub quals: Vec<u8>,
    /// Is R1 (true) or R2 (false)
    pub is_r1: bool,
}

// ============================================================================
// Core Functions
// ============================================================================

/// Build chromosome name lookup from BAM header
fn build_tid_lookup(header: &bam::HeaderView) -> Vec<String> {
    (0..header.target_count())
        .map(|tid| {
            std::str::from_utf8(header.tid2name(tid))
                .unwrap_or("unknown")
                .to_string()
        })
        .collect()
}

/// Generate WASP-style read name
fn generate_wasp_name(
    original_name: &[u8],
    r1_pos: u32,
    r2_pos: u32,
    hap_idx: usize,
    total_haps: usize,
) -> Vec<u8> {
    let mut name = original_name.to_vec();
    name.extend_from_slice(b"_WASP_");
    name.extend_from_slice(r1_pos.to_string().as_bytes());
    name.extend_from_slice(b"_");
    name.extend_from_slice(r2_pos.to_string().as_bytes());
    name.extend_from_slice(b"_");
    name.extend_from_slice(hap_idx.to_string().as_bytes());
    name.extend_from_slice(b"_");
    name.extend_from_slice(total_haps.to_string().as_bytes());
    name
}

/// Result of checking overlaps - returns ALL overlapping variants
///
/// To match baseline behavior exactly:
/// - Baseline bedtools finds ALL variants overlapping the read's genomic span
/// - Baseline bam_remapper checks ALL variants and skips if ANY is unmappable
/// - We must do the same: return ALL overlapping variants, let caller check mappability
#[derive(Debug)]
enum CheckOverlapResult {
    /// No variants overlap this read at all
    NoOverlaps,
    /// Found overlapping variants - returns Vec of (variant_idx, var_start, var_stop)
    /// Caller must check if ALL are mappable - if ANY is unmappable, skip entire read
    Found(Vec<(u32, u32, u32)>),
}

/// Check if a read overlaps any variants and return ALL of them
///
/// To match baseline behavior exactly:
/// - Returns ALL overlapping variants (baseline traversal order)
/// - Caller (generate_haplotypes_for_read) checks if ALL are mappable
/// - If ANY is unmappable → skip entire read (matching baseline bam_remapper.rs)
///
/// Returns:
/// - NoOverlaps: No variants overlap this read at all
/// - Found: All overlapping variants (baseline traversal order)
fn check_overlaps(
    read: &bam::Record,
    tid_to_name: &[String],
    trees: &FxHashMap<String, coitrees::COITree<u32, u32>>,
    store: &VariantStore,
) -> CheckOverlapResult {
    let tid = read.tid();
    if tid < 0 || tid as usize >= tid_to_name.len() {
        return CheckOverlapResult::NoOverlaps;
    }

    let chrom = &tid_to_name[tid as usize];
    let tree = match trees.get(chrom) {
        Some(t) => t,
        None => return CheckOverlapResult::NoOverlaps,
    };

    let read_start = read.pos() as i32;
    let read_end = read.reference_end() as i32 - 1;

    // Collect ALL overlapping variants using a fresh SortedQuerent per read
    let mut overlapping: Vec<(u32, u32, u32)> = Vec::new();
    let mut querent: coitrees::COITreeSortedQuerent<u32, u32> = SortedQuerent::new(tree);

    querent.query(read_start, read_end, |node| {
        let variant_idx: u32 = node.metadata.clone();
        let variant = &store.variants[variant_idx as usize];
        overlapping.push((variant_idx, variant.start, variant.stop));
    });

    if overlapping.is_empty() {
        return CheckOverlapResult::NoOverlaps;
    }

    // Sort by variant start position - empirically gives better match to baseline (3K vs 7K)
    overlapping.sort_by_key(|&(_, start, _)| start);
    CheckOverlapResult::Found(overlapping)
}

/// Deduplicate overlaps to match parse_intersect_bed behavior
///
/// Baseline keeps only the first variant for each (read_name, chrom, read_start, read_stop, mate)
/// combination. For a single read these fields are constant, so we reduce to the first overlap.
fn dedup_overlaps_for_read(overlaps: &[(u32, u32, u32)]) -> Vec<(u32, u32, u32)> {
    overlaps.iter().copied().take(1).collect()
}

/// Convert phased genotype to haplotype allele strings
/// Supports both 0/1 indexing (ref/alt) and direct allele strings.
fn genotype_to_alleles(genotype: &str, ref_allele: &str, alt_allele: &str) -> Option<(String, String)> {
    let parts: Vec<&str> = genotype.split('|').collect();
    if parts.len() != 2 {
        return None;
    }

    let to_allele = |s: &str| -> Option<String> {
        match s {
            "0" => Some(ref_allele.to_string()),
            "1" => Some(alt_allele.to_string()),
            _ => Some(s.to_string()), // Already allele string
        }
    };

    let hap1 = to_allele(parts[0])?;
    let hap2 = to_allele(parts[1])?;
    Some((hap1, hap2))
}

/// Generate haplotype sequences for a read with variants
///
/// To match baseline behavior:
/// - Baseline deduplicates intersect entries, keeping only FIRST per read
/// - The first variant (bedtools order) is the only one considered
/// - If that first variant is unmappable, the entire read is skipped
/// - Only that first variant is passed to haplotype generation
fn generate_haplotypes_for_read(
    read: &bam::Record,
    overlaps: &[(u32, u32, u32)],  // (variant_idx, var_start, var_stop)
    store: &VariantStore,
    max_seqs: usize,
) -> Option<Vec<(Vec<u8>, Vec<u8>)>> {
    let overlaps = dedup_overlaps_for_read(overlaps);

    if overlaps.is_empty() {
        // No variants - return original sequence
        let seq = read.seq().as_bytes();
        let qual = read.qual().to_vec();
        return Some(vec![(seq, qual)]);
    }

    // Build single VariantSpan (first/deduped) for haplotype generation
    let (variant_idx, _, _) = overlaps[0];
    let variant = &store.variants[variant_idx as usize];
    let (hap1, hap2) = match genotype_to_alleles(&variant.genotype, &variant.ref_allele, &variant.alt_allele) {
        Some(h) => h,
        None => return None,
    };

    let mate = if read.is_first_in_template() { 1 } else { 2 };
    let span = RemapVariantSpan {
        chrom: variant.chrom.clone(),
        start: read.pos() as u32,
        stop: read.reference_end() as u32,
        vcf_start: variant.start,
        vcf_stop: variant.stop,
        mate,
        hap1,
        hap2,
    };

    let span_refs: Vec<&RemapVariantSpan> = vec![&span];
    let remap_config = RemapConfig {
        max_seqs,
        is_phased: true,
    };

    match generate_haplotype_seqs(read, &span_refs, &remap_config) {
        Ok(Some(haps)) => Some(haps),
        _ => None, // Unmappable or error: skip this read
    }
}

/// Process a complete read pair and generate haplotype outputs
///
/// To match baseline behavior EXACTLY:
/// - If a read has variants but ALL are unmappable → skip the entire pair
/// - If a read has SOME mappable variants → process only the mappable ones
/// - Baseline processes each (read, variant) pair from bedtools intersect
/// - Unmappable variants (in introns/deletions) are skipped individually
/// - Read appears in output if ANY variant was successfully processed
fn process_pair(
    read1: &bam::Record,
    read2: &bam::Record,
    r1_overlaps: &[(u32, u32, u32)],
    r2_overlaps: &[(u32, u32, u32)],
    store: &VariantStore,
    config: &UnifiedConfig,
) -> Option<Vec<HaplotypeOutput>> {
    let mut outputs = Vec::new();

    // Generate haplotypes for each read independently
    // Returns None if read has variants but ALL are unmappable
    let r1_haps = generate_haplotypes_for_read(read1, r1_overlaps, store, config.max_seqs);
    let r2_haps = generate_haplotypes_for_read(read2, r2_overlaps, store, config.max_seqs);

    // Skip pair if either read has variants but ALL are unmappable
    // None means: read had variant overlaps but couldn't process ANY of them
    let r1_haps = match r1_haps {
        Some(haps) => haps,
        None => return None, // R1 has only unmappable variants - skip pair
    };
    let r2_haps = match r2_haps {
        Some(haps) => haps,
        None => return None, // R2 has only unmappable variants - skip pair
    };

    let r1_pos = read1.pos() as u32;
    let r2_pos = read2.pos() as u32;
    let original_name = read1.qname();

    // Generate paired outputs - WASP needs all haplotype combinations for remapping
    // When one read has more haplotypes than the other, reuse the last haplotype of the shorter
    let total_haps = r1_haps.len().max(r2_haps.len());

    for hap_idx in 0..total_haps {
        // Use the haplotype at index, or the last one if we've run out
        let r1_hap = r1_haps.get(hap_idx).unwrap_or_else(|| r1_haps.last().unwrap());
        let r2_hap = r2_haps.get(hap_idx).unwrap_or_else(|| r2_haps.last().unwrap());

        let wasp_name = generate_wasp_name(original_name, r1_pos, r2_pos, hap_idx + 1, total_haps);

        // R1 output with /1 mate suffix (matches baseline format)
        let mut r1_name = wasp_name.clone();
        r1_name.extend_from_slice(b"/1");
        outputs.push(HaplotypeOutput {
            name: r1_name,
            sequence: r1_hap.0.clone(),
            quals: r1_hap.1.clone(),
            is_r1: true,
        });

        // R2 output with /2 mate suffix (matches baseline format)
        let mut r2_name = wasp_name;
        r2_name.extend_from_slice(b"/2");
        outputs.push(HaplotypeOutput {
            name: r2_name,
            sequence: r2_hap.0.clone(),
            quals: r2_hap.1.clone(),
            is_r1: false,
        });
    }

    Some(outputs)
}

/// FASTQ writer thread - consumes haplotype outputs and writes to files
fn fastq_writer_thread(
    rx: Receiver<HaplotypeOutput>,
    r1_path: &str,
    r2_path: &str,
    counter: Arc<AtomicUsize>,
) -> Result<()> {
    // Use gzip compression with fastest level
    let r1_file = File::create(r1_path)?;
    let r2_file = File::create(r2_path)?;

    let mut r1_writer = BufWriter::with_capacity(
        1024 * 1024,
        GzEncoder::new(r1_file, Compression::fast()),
    );
    let mut r2_writer = BufWriter::with_capacity(
        1024 * 1024,
        GzEncoder::new(r2_file, Compression::fast()),
    );

    for hap in rx {
        let qual_string: Vec<u8> = hap.quals.iter().map(|&q| q + 33).collect();

        let writer = if hap.is_r1 {
            &mut r1_writer
        } else {
            &mut r2_writer
        };

        writer.write_all(b"@")?;
        writer.write_all(&hap.name)?;
        writer.write_all(b"\n")?;
        writer.write_all(&hap.sequence)?;
        writer.write_all(b"\n+\n")?;
        writer.write_all(&qual_string)?;
        writer.write_all(b"\n")?;

        counter.fetch_add(1, Ordering::Relaxed);
    }

    r1_writer.flush()?;
    r2_writer.flush()?;

    Ok(())
}

/// Unified make-reads pipeline - main entry point
///
/// Replaces: process_bam() + intersect_reads() + write_remap_bam()
///
/// # Arguments
/// * `bam_path` - Input BAM (coordinate-sorted)
/// * `bed_path` - Variant BED file (from vcf_to_bed)
/// * `r1_path` - Output R1 FASTQ (gzipped)
/// * `r2_path` - Output R2 FASTQ (gzipped)
/// * `config` - Pipeline configuration
///
/// # Returns
/// UnifiedStats with processing statistics
pub fn unified_make_reads(
    bam_path: &str,
    bed_path: &str,
    r1_path: &str,
    r2_path: &str,
    config: &UnifiedConfig,
) -> Result<UnifiedStats> {
    let mut stats = UnifiedStats::default();

    // Phase 1: Build variant store
    let t0 = Instant::now();
    eprintln!("Building variant store from {}...", bed_path);
    let store = build_variant_store(bed_path)?;
    stats.tree_build_ms = t0.elapsed().as_millis() as u64;
    eprintln!(
        "  {} chromosomes, {} variants ({}ms)",
        store.trees.len(),
        store.variants.len(),
        stats.tree_build_ms
    );

    // Phase 2: Set up writer channel
    let (tx, rx): (Sender<HaplotypeOutput>, Receiver<HaplotypeOutput>) =
        bounded(config.channel_buffer);

    let hap_counter = Arc::new(AtomicUsize::new(0));
    let hap_counter_clone = Arc::clone(&hap_counter);

    // Spawn writer thread
    let r1_owned = r1_path.to_string();
    let r2_owned = r2_path.to_string();
    let writer_handle = thread::spawn(move || {
        fastq_writer_thread(rx, &r1_owned, &r2_owned, hap_counter_clone)
    });

    // Phase 3: Stream BAM and process pairs
    let t1 = Instant::now();
    eprintln!("Streaming BAM and processing pairs...");

    let mut bam = bam::Reader::from_path(bam_path).context("Failed to open BAM")?;
    bam.set_threads(config.read_threads).ok();

    let header = bam.header().clone();
    let tid_to_name = build_tid_lookup(&header);

    // Pair buffer: read_name -> first-seen mate
    let mut pair_buffer: FxHashMap<Vec<u8>, bam::Record> = FxHashMap::default();
    pair_buffer.reserve(1_000_000);

    for result in bam.records() {
        let read = result?;
        stats.total_reads += 1;

        // Skip reads that don't pass baseline filtering:
        // IMPORTANT: Match bam_intersect.rs exactly (unmapped, secondary, supplementary)
        // Do NOT filter on QC fail (0x200) or duplicate (0x400) here because:
        // - bam_filter phase2 adds names to remap set (filters qc/dup on primary read)
        // - bam_filter phase3 writes BOTH mates by name (no filtering!)
        // - bam_intersect filters unmapped, secondary, supplementary ONLY
        // - If one mate is qc_fail but the other overlaps, BOTH go to remap.bam
        // - So we must process qc_fail/duplicate reads to match baseline exactly
        if read.is_unmapped() || read.is_secondary() || read.is_supplementary() {
            continue;
        }
        // Also check proper_pair like bam_remapper.rs:374 does
        if !read.is_proper_pair() {
            continue;
        }

        let read_name = read.qname().to_vec();

        // Try to complete a pair
        match pair_buffer.remove(&read_name) {
            Some(mate) => {
                stats.pairs_processed += 1;

                // Ensure read1 is first in template
                let (r1, r2) = if read.is_first_in_template() {
                    (read, mate)
                } else {
                    (mate, read)
                };

                // Check overlaps for both mates - returns ALL overlapping variants
                let r1_result = check_overlaps(&r1, &tid_to_name, &store.trees, &store);
                let r2_result = check_overlaps(&r2, &tid_to_name, &store.trees, &store);

                // Extract variant vectors from results
                let r1_variants = match &r1_result {
                    CheckOverlapResult::Found(v) => v.clone(),
                    CheckOverlapResult::NoOverlaps => Vec::new(),
                };
                let r2_variants = match &r2_result {
                    CheckOverlapResult::Found(v) => v.clone(),
                    CheckOverlapResult::NoOverlaps => Vec::new(),
                };

                // Process based on overlap results
                if r1_variants.is_empty() && r2_variants.is_empty() {
                    // No variants at all - this pair would go to keep.bam
                    stats.pairs_kept += 1;
                } else {
                    // At least one mate has variants - pass ALL to process_pair
                    // process_pair will check if ANY variant is unmappable and return None
                    // This matches baseline behavior: skip entire pair if ANY variant unmappable
                    match process_pair(
                        &r1,
                        &r2,
                        &r1_variants,
                        &r2_variants,
                        &store,
                        config,
                    ) {
                        Some(outputs) => {
                            stats.pairs_with_variants += 1;
                            for output in outputs {
                                tx.send(output).ok();
                            }
                        }
                        None => {
                            // Haplotype generation failed - likely because a variant is in
                            // an intron/deletion (unmappable). This matches baseline behavior.
                            stats.pairs_skipped_unmappable += 1;
                        }
                    }
                }
            }
            None => {
                // First mate seen - buffer it
                pair_buffer.insert(read_name, read);
            }
        }

        // Progress reporting
        if stats.total_reads % 10_000_000 == 0 {
            eprintln!(
                "  {} reads, {} pairs, {} with variants",
                stats.total_reads, stats.pairs_processed, stats.pairs_with_variants
            );
        }
    }

    stats.orphan_reads = pair_buffer.len();
    stats.bam_stream_ms = t1.elapsed().as_millis() as u64;

    eprintln!(
        "  {} orphan reads (mate not found)",
        stats.orphan_reads
    );

    // Close sender to signal writer thread to finish
    drop(tx);

    // Wait for writer thread
    writer_handle
        .join()
        .map_err(|_| anyhow::anyhow!("Writer thread panicked"))??;

    stats.haplotypes_written = hap_counter.load(Ordering::Relaxed);

    eprintln!("Unified pipeline complete:");
    eprintln!("  Total reads: {}", stats.total_reads);
    eprintln!("  Pairs processed: {}", stats.pairs_processed);
    eprintln!("  Pairs with variants: {}", stats.pairs_with_variants);
    eprintln!("  Pairs kept (no variants): {}", stats.pairs_kept);
    eprintln!("  Pairs skipped (unmappable): {}", stats.pairs_skipped_unmappable);
    eprintln!("  Pairs haplotype failed: {}", stats.pairs_haplotype_failed);
    eprintln!("  Haplotypes written: {}", stats.haplotypes_written);
    eprintln!(
        "  Time: {}ms tree build + {}ms BAM stream",
        stats.tree_build_ms, stats.bam_stream_ms
    );

    Ok(stats)
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_generate_wasp_name() {
        let name = generate_wasp_name(b"ERR123456.1000", 12345, 67890, 1, 2);
        let expected = b"ERR123456.1000_WASP_12345_67890_1_2";
        assert_eq!(name, expected.to_vec());
    }

    #[test]
    fn test_unified_config_default() {
        let config = UnifiedConfig::default();
        assert_eq!(config.read_threads, 8);
        assert_eq!(config.max_seqs, 64);
        assert_eq!(config.channel_buffer, 50_000);
    }

    #[test]
    fn test_unified_stats_default() {
        let stats = UnifiedStats::default();
        assert_eq!(stats.total_reads, 0);
        assert_eq!(stats.pairs_processed, 0);
        assert_eq!(stats.haplotypes_written, 0);
    }
}
