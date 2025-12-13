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
use coitrees::SortedQuerent;
use crossbeam_channel::{bounded, Receiver, Sender};
use flate2::Compression;
use gzp::{deflate::Gzip, ZBuilder};
use itoa::Buffer as ItoaBuffer;
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
use crate::bam_remapper::{
    apply_trim_combination, calculate_indel_delta, classify_variant_location,
    generate_haplotype_seqs, generate_trim_combinations, IndelConfig, RemapConfig, VariantLocation,
    VariantSpan as RemapVariantSpan,
};

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
    /// Number of compression threads per FASTQ file (0 = auto)
    pub compression_threads: usize,
    /// Compress output FASTQs (set to false for named pipe streaming)
    pub compress_output: bool,
    /// Enable INDEL mode with length-preserving trim combinations
    pub indel_mode: bool,
    /// Maximum INDEL size to handle (larger INDELs are skipped)
    pub max_indel_size: usize,
    /// Optional path to output read names of "keep-no-flip" pairs
    /// These are pairs that overlap variants but don't need remapping
    pub keep_no_flip_names_path: Option<String>,
    /// Optional path to output read names of pairs that were sent for remapping
    /// These are the pairs that have haplotypes written to FASTQ
    /// Use this to create the correct reference BAM for filter_bam_wasp
    pub remap_names_path: Option<String>,
}

impl Default for UnifiedConfig {
    fn default() -> Self {
        Self {
            read_threads: 8,
            max_seqs: 64,
            channel_buffer: 50_000,
            compression_threads: 4, // 4 threads per FASTQ file for parallel gzip
            compress_output: true,  // Default to compressed for disk storage
            indel_mode: false,      // Default to SNV-only mode for backward compatibility
            max_indel_size: 50,     // 50bp max INDEL (standard threshold)
            keep_no_flip_names_path: None, // Don't output keep-no-flip names by default
            remap_names_path: None, // Don't output remap names by default
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
    /// Pairs overlapping SNVs only (no indels)
    pub pairs_with_snvs_only: usize,
    /// Pairs overlapping indels only (no SNVs)
    pub pairs_with_indels_only: usize,
    /// Pairs overlapping both SNVs and indels
    pub pairs_with_snvs_and_indels: usize,
    /// Total haplotype reads written
    pub haplotypes_written: usize,
    /// Pairs kept (no variants at all)
    pub pairs_kept: usize,
    /// Pairs that overlap variants but don't need remapping (sequence unchanged)
    /// These should be KEPT in final output, not discarded!
    pub pairs_keep_no_flip: usize,
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

impl UnifiedStats {
    /// Merge stats from multiple threads into a single aggregate
    pub fn merge(self, other: Self) -> Self {
        Self {
            total_reads: self.total_reads + other.total_reads,
            pairs_processed: self.pairs_processed + other.pairs_processed,
            pairs_with_variants: self.pairs_with_variants + other.pairs_with_variants,
            pairs_with_snvs_only: self.pairs_with_snvs_only + other.pairs_with_snvs_only,
            pairs_with_indels_only: self.pairs_with_indels_only + other.pairs_with_indels_only,
            pairs_with_snvs_and_indels: self.pairs_with_snvs_and_indels
                + other.pairs_with_snvs_and_indels,
            haplotypes_written: self.haplotypes_written + other.haplotypes_written,
            pairs_kept: self.pairs_kept + other.pairs_kept,
            pairs_keep_no_flip: self.pairs_keep_no_flip + other.pairs_keep_no_flip,
            pairs_skipped_unmappable: self.pairs_skipped_unmappable
                + other.pairs_skipped_unmappable,
            pairs_haplotype_failed: self.pairs_haplotype_failed + other.pairs_haplotype_failed,
            orphan_reads: self.orphan_reads + other.orphan_reads,
            // Keep maximum time values (they represent wall clock for parallel execution)
            tree_build_ms: self.tree_build_ms.max(other.tree_build_ms),
            bam_stream_ms: self.bam_stream_ms.max(other.bam_stream_ms),
        }
    }
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
    /// Whether original BAM read was on reverse strand
    /// IMPORTANT: Used to reverse-complement before FASTQ output
    /// BAM stores reverse-strand reads as already rev-comped, but FASTQ needs original orientation
    pub is_reverse: bool,
}

/// A paired haplotype output (R1 + R2 together) for atomic writing
/// This ensures paired reads are written in the same order to both FASTQ files
#[derive(Debug, Clone)]
pub struct HaplotypePair {
    pub r1: HaplotypeOutput,
    pub r2: HaplotypeOutput,
    /// Shared trim combination ID (both mates use same combo for coordinated trimming)
    /// Encoded as: hap_idx * 1000 + combo_idx
    pub trim_combo_id: u16,
    /// Total number of trim combinations for this read pair (for filtering denominator)
    pub total_combos: u16,
    /// Expected positions for this haplotype+trim combo (variant-aware)
    pub exp_pos1: u32,
    pub exp_pos2: u32,
    /// Bitmask describing overlap types for the ORIGINAL read pair:
    /// 1 = SNV-only, 2 = INDEL-only, 3 = SNV+INDEL.
    pub overlap_mask: u8,
}

/// Result of processing a read pair with variants
/// This enum distinguishes between pairs that need remapping vs those that can be kept as-is
#[derive(Debug)]
pub enum ProcessPairResult {
    /// Pair needs remapping - contains haplotype pairs to write to FASTQ
    NeedsRemap(Vec<HaplotypePair>),
    /// Pair overlaps variants but sequences are unchanged - keep original reads
    /// (Both haplotypes match original sequence, so no allele flip needed)
    KeepAsIs,
    /// Pair is unmappable (variant in intron/deletion) - discard
    Unmappable,
}

// ============================================================================
// Core Functions
// ============================================================================

#[inline]
fn complement_base(b: u8) -> u8 {
    match b {
        b'A' | b'a' => b'T',
        b'T' | b't' => b'A',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'N' | b'n' => b'N',
        _ => b'N',
    }
}

/// Compute expected reference start for a read in a haplotype/trim combo.
///
/// CIGAR-AWARE: Uses `classify_variant_location` from bam_remapper to properly
/// classify variants relative to the read's CIGAR-derived reference span.
///
/// Only variants classified as:
/// - `Upstream`: entirely before read start → shift expected position
/// - `SpansStart`: deletion/insertion spanning read start → shift expected position
///
/// Variants classified as `WithinRead` or `Downstream` do NOT shift the anchor.
fn expected_start_upstream_only(
    read: &bam::Record,
    overlaps: &[(u32, u32, u32)],
    store: &VariantStore,
    hap_idx: usize,
) -> u32 {
    let read_start = read.pos() as i64;
    let mut shift: i64 = 0;

    for (idx, _s, _e) in overlaps {
        let v = &store.variants[*idx as usize];

        // Get variant's reference span
        let v_start = v.start;
        let v_stop = v.stop;

        // Use CIGAR-aware classification from bam_remapper
        let location = classify_variant_location(read, v_start, v_stop);

        // Get haplotype-specific allele for delta calculation
        let (hap1, hap2) = genotype_to_alleles(&v.genotype, &v.ref_allele, &v.alt_allele)
            .unwrap_or_else(|| (v.ref_allele.clone(), v.alt_allele.clone()));
        let ref_len = v.ref_allele.len() as i64;
        let alt_len = if hap_idx == 0 {
            hap1.len() as i64
        } else {
            hap2.len() as i64
        };
        let delta = alt_len - ref_len;

        match location {
            VariantLocation::Upstream => {
                // Fully upstream variant - shifts expected position
                shift += delta;
            }
            VariantLocation::SpansStart => {
                // Variant spans read start boundary
                // Deletions and insertions starting before read shift position
                if delta != 0 {
                    shift += delta;
                }
                // SNVs at boundary: no shift
            }
            VariantLocation::WithinRead | VariantLocation::Downstream => {
                // No shift for within-read or downstream variants
            }
        }
    }

    let exp = read_start + shift;
    if exp < 0 {
        0
    } else {
        exp as u32
    }
}

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
    let mut name = Vec::with_capacity(original_name.len() + 64);
    name.extend_from_slice(original_name);
    name.extend_from_slice(b"_WASP_");
    let mut tmp = ItoaBuffer::new();
    name.extend_from_slice(tmp.format(r1_pos).as_bytes());
    name.extend_from_slice(b"_");
    name.extend_from_slice(tmp.format(r2_pos).as_bytes());
    name.extend_from_slice(b"_");
    name.extend_from_slice(tmp.format(hap_idx).as_bytes());
    name.extend_from_slice(b"_");
    name.extend_from_slice(tmp.format(total_haps).as_bytes());
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
        let variant_idx: u32 = u32::from(node.metadata.clone());
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

/// Classify overlap types for a read pair.
///
/// Mask bits:
/// - 1: SNV overlap present (ref/alt same length)
/// - 2: INDEL overlap present (ref/alt different length)
fn overlap_mask_for_pair(
    r1_variants: &[(u32, u32, u32)],
    r2_variants: &[(u32, u32, u32)],
    store: &VariantStore,
) -> u8 {
    let mut has_snv = false;
    let mut has_indel = false;
    for (idx, _s, _e) in r1_variants.iter().chain(r2_variants.iter()) {
        let v = &store.variants[*idx as usize];
        if v.ref_allele.len() != v.alt_allele.len() {
            has_indel = true;
        } else {
            has_snv = true;
        }
        if has_snv && has_indel {
            break;
        }
    }
    match (has_snv, has_indel) {
        (true, false) => 1,
        (false, true) => 2,
        (true, true) => 3,
        _ => 0,
    }
}

fn increment_overlap_stats(stats: &mut UnifiedStats, mask: u8) {
    match mask {
        1 => stats.pairs_with_snvs_only += 1,
        2 => stats.pairs_with_indels_only += 1,
        3 => stats.pairs_with_snvs_and_indels += 1,
        _ => {}
    }
}

/// Convert phased genotype to haplotype allele strings
/// Supports both 0/1 indexing (ref/alt) and direct allele strings.
fn genotype_to_alleles(
    genotype: &str,
    ref_allele: &str,
    alt_allele: &str,
) -> Option<(String, String)> {
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
/// FIX: Process ALL overlapping variants (not just first) to match Python DEV behavior.
/// For phased data, this generates exactly 2 haplotype sequences with ALL alleles substituted.
///
/// # Algorithm (matching Python DEV make_remap_reads.py):
/// 1. Collect ALL variants overlapping the read
/// 2. Sort by genomic position for deterministic substitution order
/// 3. Build VariantSpan for each variant
/// 4. Call generate_haplotype_seqs which:
///    - Splits read sequence at all variant positions
///    - Substitutes hap1 alleles at odd indices -> haplotype 1
///    - Substitutes hap2 alleles at odd indices -> haplotype 2
/// 5. Return 2 haplotype sequences (for phased data)
fn generate_haplotypes_for_read(
    read: &bam::Record,
    overlaps: &[(u32, u32, u32)], // (variant_idx, var_start, var_stop)
    store: &VariantStore,
    max_seqs: usize,
) -> Option<Vec<(Vec<u8>, Vec<u8>)>> {
    if overlaps.is_empty() {
        // No variants - return original sequence TWICE (matches baseline bam_remapper.rs)
        // This is needed for correct zip pairing with the other read's haplotypes
        let seq = read.seq().as_bytes();
        let qual = read.qual().to_vec();
        return Some(vec![(seq.clone(), qual.clone()), (seq, qual)]);
    }

    // Sort overlaps by genomic position for deterministic substitution order
    let mut sorted_overlaps = overlaps.to_vec();
    sorted_overlaps.sort_by_key(|&(_, start, _)| start);

    // Build VariantSpan for EACH overlapping variant (not just first!)
    let mate = if read.is_first_in_template() { 1 } else { 2 };
    let mut spans: Vec<RemapVariantSpan> = Vec::with_capacity(sorted_overlaps.len());

    for (variant_idx, _, _) in &sorted_overlaps {
        let variant = &store.variants[*variant_idx as usize];
        let (hap1, hap2) = match genotype_to_alleles(
            &variant.genotype,
            &variant.ref_allele,
            &variant.alt_allele,
        ) {
            Some(h) => h,
            None => return None, // Invalid genotype format
        };

        spans.push(RemapVariantSpan {
            chrom: variant.chrom.clone(),
            start: read.pos() as u32,
            stop: read.reference_end() as u32,
            vcf_start: variant.start,
            vcf_stop: variant.stop,
            mate,
            hap1,
            hap2,
        });
    }

    // Pass ALL spans to generate_haplotype_seqs (which already supports multiple variants)
    let span_refs: Vec<&RemapVariantSpan> = spans.iter().collect();
    let remap_config = RemapConfig {
        max_seqs,
        is_phased: true,
    };

    match generate_haplotype_seqs(read, &span_refs, &remap_config) {
        Ok(Some(haps)) => Some(haps),
        _ => None, // Unmappable or error: skip this read
    }
}

/// Process a complete read pair and generate haplotype pair outputs
///
/// To match baseline behavior EXACTLY:
/// - If a read has variants but ALL are unmappable → skip the entire pair
/// - If a read has SOME mappable variants → process only the mappable ones
/// - Baseline processes each (read, variant) pair from bedtools intersect
/// - Unmappable variants (in introns/deletions) are skipped individually
/// - Read appears in output if ANY variant was successfully processed
///
/// Returns ProcessPairResult to distinguish between:
/// - NeedsRemap: pairs that need remapping (has sequence changes)
/// - KeepAsIs: pairs that overlap variants but have no sequence changes (keep original)
/// - Unmappable: pairs where variant is in intron/deletion (discard)
fn process_pair(
    read1: &bam::Record,
    read2: &bam::Record,
    r1_overlaps: &[(u32, u32, u32)],
    r2_overlaps: &[(u32, u32, u32)],
    store: &VariantStore,
    config: &UnifiedConfig,
    overlap_mask: u8,
) -> ProcessPairResult {
    let mut outputs = Vec::new();

    // Original sequences for unchanged check
    let r1_original = read1.seq().as_bytes();
    let r2_original = read2.seq().as_bytes();

    // Generate haplotypes for each read independently
    // Returns None if read has variants but ALL are unmappable
    // Returns exactly 2 haplotypes: either (orig, orig) for no variants, or (hap1, hap2) for variants
    let r1_haps = match generate_haplotypes_for_read(read1, r1_overlaps, store, config.max_seqs) {
        Some(h) => h,
        None => return ProcessPairResult::Unmappable,
    };
    let r2_haps = match generate_haplotypes_for_read(read2, r2_overlaps, store, config.max_seqs) {
        Some(h) => h,
        None => return ProcessPairResult::Unmappable,
    };

    let r1_pos = read1.pos() as u32;
    let r2_pos = read2.pos() as u32;
    let original_name = read1.qname();

    // First pass: filter to only pairs where at least one sequence changed
    // This matches Python DEV behavior where total = actual count of written pairs
    let changed_pairs: Vec<_> = r1_haps
        .iter()
        .zip(r2_haps.iter())
        .filter(|(r1_hap, r2_hap)| {
            // Keep if at least one read is changed (matches baseline bam_remapper.rs line 476-479)
            r1_hap.0 != r1_original || r2_hap.0 != r2_original
        })
        .collect();

    let total_seqs = changed_pairs.len();
    if total_seqs == 0 {
        // No sequence changes needed - the read already has reference alleles
        // This is NOT an error - the read should be KEPT, just not remapped
        return ProcessPairResult::KeepAsIs;
    }

    // Track reverse strand status for FASTQ output
    // IMPORTANT: BAM stores reverse-strand reads as already reverse-complemented
    // For FASTQ output (for remapping), we need to reverse-complement back to original orientation
    let r1_is_reverse = read1.is_reverse();
    let r2_is_reverse = read2.is_reverse();

    // Second pass: generate outputs with correct total count
    for (write_idx, (r1_hap, r2_hap)) in changed_pairs.into_iter().enumerate() {
        // Use actual count of changed pairs as total (matches Python DEV make_remap_reads.py)
        let wasp_name =
            generate_wasp_name(original_name, r1_pos, r2_pos, write_idx + 1, total_seqs);

        // R1 output
        let mut r1_name = wasp_name.clone();
        r1_name.extend_from_slice(b"/1");
        let r1_output = HaplotypeOutput {
            name: r1_name,
            sequence: r1_hap.0.clone(),
            quals: r1_hap.1.clone(),
            is_r1: true,
            is_reverse: r1_is_reverse,
        };

        // R2 output
        let mut r2_name = wasp_name;
        r2_name.extend_from_slice(b"/2");
        let r2_output = HaplotypeOutput {
            name: r2_name,
            sequence: r2_hap.0.clone(),
            quals: r2_hap.1.clone(),
            is_r1: false,
            is_reverse: r2_is_reverse,
        };

        // Bundle as pair for atomic writing
        // For SNV-only mode, use default trim combo values (no trimming)
        outputs.push(HaplotypePair {
            r1: r1_output,
            r2: r2_output,
            trim_combo_id: 0, // No trim combo in SNV mode
            total_combos: 1,  // Single combination (no trimming)
            exp_pos1: r1_pos,
            exp_pos2: r2_pos,
            overlap_mask,
        });
    }

    if outputs.is_empty() {
        ProcessPairResult::KeepAsIs
    } else {
        ProcessPairResult::NeedsRemap(outputs)
    }
}

/// Process a complete read pair with coordinated trim combinations for INDEL support
///
/// This is the INDEL-aware version that:
/// 1. Generates raw haplotype sequences (may be extended for insertions)
/// 2. Calculates the max INDEL delta across both reads
/// 3. Generates coordinated trim combinations (same for both R1 and R2)
/// 4. Applies the SAME trim to both mates, ensuring length preservation
///
/// Returns HaplotypePairs (R1+R2 together) with trim_combo_id for filtering
#[allow(dead_code)]
fn process_pair_with_trims(
    read1: &bam::Record,
    read2: &bam::Record,
    r1_overlaps: &[(u32, u32, u32)],
    r2_overlaps: &[(u32, u32, u32)],
    store: &VariantStore,
    config: &UnifiedConfig,
    indel_config: &IndelConfig,
    overlap_mask: u8,
) -> ProcessPairResult {
    let mut outputs = Vec::new();

    let r1_original_len = read1.seq().len();
    let r2_original_len = read2.seq().len();
    let r1_original = read1.seq().as_bytes();
    let r2_original = read2.seq().as_bytes();

    // Generate raw haplotypes for each read (may have different lengths due to INDELs)
    let r1_haps = match generate_haplotypes_for_read(read1, r1_overlaps, store, config.max_seqs)
    {
        Some(h) => h,
        None => return ProcessPairResult::Unmappable,
    };
    let r2_haps = match generate_haplotypes_for_read(read2, r2_overlaps, store, config.max_seqs)
    {
        Some(h) => h,
        None => return ProcessPairResult::Unmappable,
    };

    // ---------------------------------------------------------------------
    // New approach: trim combinations per read (pi guidance). We generate
    // combos independently for R1/R2 based on their own deltas and take
    // the cartesian product per haplotype pair. A small cap prevents
    // explosion on large deltas.
    // ---------------------------------------------------------------------
    const MAX_TRIM_COMBO_PRODUCT: usize = 256;
    let r1_pos = read1.pos() as u32;
    let r2_pos = read2.pos() as u32;
    let original_name = read1.qname();

    // Track reverse strand status for FASTQ output
    let r1_is_reverse = read1.is_reverse();
    let r2_is_reverse = read2.is_reverse();

    // Collect all outputs first to compute total_seqs accurately
    struct PendingOutput {
        hap_idx: usize,
        combo_idx_r1: usize,
        combo_idx_r2: usize,
        total_combos_pair: u16,
        r1_delta: i32,
        r2_delta: i32,
        r1_seq: Vec<u8>,
        r1_qual: Vec<u8>,
        r2_seq: Vec<u8>,
        r2_qual: Vec<u8>,
        exp_pos1: u32,
        exp_pos2: u32,
    }
    let mut pending: Vec<PendingOutput> = Vec::new();
    let mut any_non_skipped_hap = false;

    for (hap_idx, (r1_hap, r2_hap)) in r1_haps.iter().zip(r2_haps.iter()).enumerate() {
        let r1_delta = calculate_indel_delta(r1_hap.0.len(), r1_original_len);
        let r2_delta = calculate_indel_delta(r2_hap.0.len(), r2_original_len);
        // CIGAR-aware: only upstream variants shift the start anchor
        let exp_pos1 = expected_start_upstream_only(read1, r1_overlaps, store, hap_idx);
        let exp_pos2 = expected_start_upstream_only(read2, r2_overlaps, store, hap_idx);

        // Skip pairs with indels larger than threshold
        if (r1_delta.abs() as usize) > indel_config.max_indel_size
            || (r2_delta.abs() as usize) > indel_config.max_indel_size
        {
            if indel_config.skip_large_indels {
                continue;
            }
        }
        any_non_skipped_hap = true;

        let mut r1_combos = generate_trim_combinations(r1_delta, r1_original_len);
        let mut r2_combos = generate_trim_combinations(r2_delta, r2_original_len);

        // Cap combo explosion (sqrt of max product)
        let max_per_side = (MAX_TRIM_COMBO_PRODUCT as f64).sqrt().floor() as usize;
        if r1_combos.len() * r2_combos.len() > MAX_TRIM_COMBO_PRODUCT {
            r1_combos.truncate(max_per_side.max(1));
            r2_combos.truncate(max_per_side.max(1));
        }
        let total_combos_pair = (r1_combos.len() * r2_combos.len()) as u16;

        for (combo_idx_r1, trim_r1) in r1_combos.iter().enumerate() {
            let (r1_seq, r1_qual) =
                apply_trim_combination(&r1_hap.0, &r1_hap.1, r1_original_len, trim_r1);
            for (combo_idx_r2, trim_r2) in r2_combos.iter().enumerate() {
                let (r2_seq, r2_qual) =
                    apply_trim_combination(&r2_hap.0, &r2_hap.1, r2_original_len, trim_r2);

                // Skip if both unchanged from original
                if r1_seq == r1_original && r2_seq == r2_original {
                    continue;
                }

                pending.push(PendingOutput {
                    hap_idx,
                    combo_idx_r1,
                    combo_idx_r2,
                    total_combos_pair,
                    r1_delta,
                    r2_delta,
                    r1_seq: r1_seq.clone(),
                    r1_qual: r1_qual.clone(),
                    r2_seq: r2_seq.clone(),
                    r2_qual: r2_qual.clone(),
                    exp_pos1,
                    exp_pos2,
                });
            }
        }
    }

    let total_seqs = pending.len();
    if total_seqs == 0 {
        if any_non_skipped_hap {
            return ProcessPairResult::KeepAsIs;
        }
        return ProcessPairResult::Unmappable;
    }

    for (seq_idx, p) in pending.into_iter().enumerate() {
        let trim_combo_id = ((p.combo_idx_r1 as u16) << 8) | (p.combo_idx_r2 as u16);
        let wasp_name = generate_wasp_name_extended(
            original_name,
            r1_pos,
            r2_pos,
            seq_idx + 1, // 1-based sequence index
            total_seqs,  // total expected sequences
            trim_combo_id,
            p.total_combos_pair,
            p.r1_delta,
            p.r2_delta,
        );

        // R1 output
        let mut r1_name = wasp_name.clone();
        r1_name.extend_from_slice(b"/1");
        let r1_output = HaplotypeOutput {
            name: r1_name,
            sequence: p.r1_seq,
            quals: p.r1_qual,
            is_r1: true,
            is_reverse: r1_is_reverse,
        };

        // R2 output
        let mut r2_name = wasp_name;
        r2_name.extend_from_slice(b"/2");
        let r2_output = HaplotypeOutput {
            name: r2_name,
            sequence: p.r2_seq,
            quals: p.r2_qual,
            is_r1: false,
            is_reverse: r2_is_reverse,
        };

        outputs.push(HaplotypePair {
            r1: r1_output,
            r2: r2_output,
            trim_combo_id,
            total_combos: p.total_combos_pair,
            exp_pos1: p.exp_pos1,
            exp_pos2: p.exp_pos2,
            overlap_mask,
        });
    }

    ProcessPairResult::NeedsRemap(outputs)
}

/// Generate extended WASP-style read name including trim combo information
/// Format: {name}_WASP_{pos1}_{pos2}_{seq}_{total}_{trim_combo}_{total_combos}
fn generate_wasp_name_extended(
    original_name: &[u8],
    r1_pos: u32,
    r2_pos: u32,
    hap_idx: usize,
    total_haps: usize,
    trim_combo_id: u16,
    total_combos: u16,
    r1_delta: i32,
    r2_delta: i32,
) -> Vec<u8> {
    let mut name = Vec::with_capacity(original_name.len() + 128);
    name.extend_from_slice(original_name);
    name.extend_from_slice(b"_WASP_");
    let mut tmp = ItoaBuffer::new();
    name.extend_from_slice(tmp.format(r1_pos).as_bytes());
    name.extend_from_slice(b"_");
    name.extend_from_slice(tmp.format(r2_pos).as_bytes());
    name.extend_from_slice(b"_");
    name.extend_from_slice(tmp.format(hap_idx).as_bytes());
    name.extend_from_slice(b"_");
    name.extend_from_slice(tmp.format(total_haps).as_bytes());
    name.extend_from_slice(b"_");
    name.extend_from_slice(tmp.format(trim_combo_id).as_bytes());
    name.extend_from_slice(b"_");
    name.extend_from_slice(tmp.format(total_combos).as_bytes());
    name.extend_from_slice(b"_");
    name.extend_from_slice(tmp.format(r1_delta.abs()).as_bytes());
    name.extend_from_slice(b"_");
    name.extend_from_slice(tmp.format(r2_delta.abs()).as_bytes());
    name
}

/// Helper to write a single FASTQ record.
///
/// Uses caller-provided scratch buffers to avoid per-record allocations.
fn write_fastq_record<W: Write>(
    writer: &mut W,
    hap: &HaplotypeOutput,
    seq_buf: &mut Vec<u8>,
    qual_buf: &mut Vec<u8>,
) -> Result<()> {
    writer.write_all(b"@")?;
    writer.write_all(&hap.name)?;
    writer.write_all(b"\n")?;

    // Sequence
    if hap.is_reverse {
        seq_buf.clear();
        seq_buf.reserve(hap.sequence.len());
        seq_buf.extend(hap.sequence.iter().rev().map(|&b| complement_base(b)));
        writer.write_all(seq_buf)?;
    } else {
        writer.write_all(&hap.sequence)?;
    }
    writer.write_all(b"\n+\n")?;

    // Quals (+33, reverse if needed)
    qual_buf.clear();
    qual_buf.reserve(hap.quals.len());
    if hap.is_reverse {
        for &q in hap.quals.iter().rev() {
            qual_buf.push(q + 33);
        }
    } else {
        for &q in &hap.quals {
            qual_buf.push(q + 33);
        }
    }
    writer.write_all(qual_buf)?;
    writer.write_all(b"\n")?;
    Ok(())
}

/// FASTQ writer thread - consumes haplotype PAIRS and writes atomically to files
/// Uses gzp for parallel gzip compression (pigz-like) when compress=true
/// Uses plain buffered write when compress=false (faster for named pipes/streaming)
///
/// CRITICAL: Receives HaplotypePair to ensure R1 and R2 are written in the same order
/// This fixes the parallel pipeline bug where R1/R2 could get out of sync
fn fastq_writer_thread(
    rx: Receiver<HaplotypePair>,
    r1_path: &str,
    r2_path: &str,
    sidecar_path: &str,
    counter: Arc<AtomicUsize>,
    compression_threads: usize,
    compress: bool,
) -> Result<()> {
    let r1_file = File::create(r1_path)?;
    let r2_file = File::create(r2_path)?;
    let mut sidecar = File::create(sidecar_path)?;
    let mut seq_buf: Vec<u8> = Vec::new();
    let mut qual_buf: Vec<u8> = Vec::new();
    let mut itoa_buf = ItoaBuffer::new();

    if compress {
        // Use gzp for parallel gzip compression (similar to pigz)
        // This provides significant speedup for I/O-bound workloads
        let mut r1_writer = ZBuilder::<Gzip, _>::new()
            .num_threads(compression_threads)
            .compression_level(Compression::fast())
            .from_writer(BufWriter::with_capacity(1024 * 1024, r1_file));

        let mut r2_writer = ZBuilder::<Gzip, _>::new()
            .num_threads(compression_threads)
            .compression_level(Compression::fast())
            .from_writer(BufWriter::with_capacity(1024 * 1024, r2_file));

        for pair in rx {
            // Write R1 and R2 atomically - they arrive together and are written together
            write_fastq_record(&mut r1_writer, &pair.r1, &mut seq_buf, &mut qual_buf)?;
            write_fastq_record(&mut r2_writer, &pair.r2, &mut seq_buf, &mut qual_buf)?;
            // Sidecar: qname  exp_pos1  exp_pos2  trim_combo_id  total_combos  overlap_mask
            sidecar.write_all(&pair.r1.name)?;
            sidecar.write_all(b"\t")?;
            sidecar.write_all(itoa_buf.format(pair.exp_pos1).as_bytes())?;
            sidecar.write_all(b"\t")?;
            sidecar.write_all(itoa_buf.format(pair.exp_pos2).as_bytes())?;
            sidecar.write_all(b"\t")?;
            sidecar.write_all(itoa_buf.format(pair.trim_combo_id).as_bytes())?;
            sidecar.write_all(b"\t")?;
            sidecar.write_all(itoa_buf.format(pair.total_combos).as_bytes())?;
            sidecar.write_all(b"\t")?;
            sidecar.write_all(&[b'0' + pair.overlap_mask])?;
            sidecar.write_all(b"\n")?;
            counter.fetch_add(2, Ordering::Relaxed); // Count both reads
        }

        // Finish flushes and finalizes the gzip streams
        r1_writer.finish().context("Failed to finish R1 gzip")?;
        r2_writer.finish().context("Failed to finish R2 gzip")?;
    } else {
        // Uncompressed output - faster for named pipes and streaming to STAR
        // Use larger buffer (4MB) for better throughput
        let mut r1_writer = BufWriter::with_capacity(4 * 1024 * 1024, r1_file);
        let mut r2_writer = BufWriter::with_capacity(4 * 1024 * 1024, r2_file);

        for pair in rx {
            // Write R1 and R2 atomically - they arrive together and are written together
            write_fastq_record(&mut r1_writer, &pair.r1, &mut seq_buf, &mut qual_buf)?;
            write_fastq_record(&mut r2_writer, &pair.r2, &mut seq_buf, &mut qual_buf)?;
            // Sidecar: qname  exp_pos1  exp_pos2  trim_combo_id  total_combos  overlap_mask
            sidecar.write_all(&pair.r1.name)?;
            sidecar.write_all(b"\t")?;
            sidecar.write_all(itoa_buf.format(pair.exp_pos1).as_bytes())?;
            sidecar.write_all(b"\t")?;
            sidecar.write_all(itoa_buf.format(pair.exp_pos2).as_bytes())?;
            sidecar.write_all(b"\t")?;
            sidecar.write_all(itoa_buf.format(pair.trim_combo_id).as_bytes())?;
            sidecar.write_all(b"\t")?;
            sidecar.write_all(itoa_buf.format(pair.total_combos).as_bytes())?;
            sidecar.write_all(b"\t")?;
            sidecar.write_all(&[b'0' + pair.overlap_mask])?;
            sidecar.write_all(b"\n")?;
            counter.fetch_add(2, Ordering::Relaxed); // Count both reads
        }

        // Flush uncompressed writers
        r1_writer.flush().context("Failed to flush R1")?;
        r2_writer.flush().context("Failed to flush R2")?;
    }

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

    // Phase 2: Set up writer channel (sends pairs for atomic writing)
    let (tx, rx): (Sender<HaplotypePair>, Receiver<HaplotypePair>) = bounded(config.channel_buffer);

    let hap_counter = Arc::new(AtomicUsize::new(0));
    let hap_counter_clone = Arc::clone(&hap_counter);

    // Spawn writer thread (with optional compression)
    let r1_owned = r1_path.to_string();
    let r2_owned = r2_path.to_string();
    let sidecar_owned = format!("{}.expected_positions.tsv", r1_owned);
    let compression_threads = config.compression_threads;
    let compress = config.compress_output;
    let writer_handle = thread::spawn(move || {
        fastq_writer_thread(
            rx,
            &r1_owned,
            &r2_owned,
            &sidecar_owned,
            hap_counter_clone,
            compression_threads,
            compress,
        )
    });

    // Optional: Set up keep-no-flip names output
    let mut keep_no_flip_writer: Option<BufWriter<File>> =
        config.keep_no_flip_names_path.as_ref().map(|path| {
            let file = File::create(path).expect("Failed to create keep_no_flip_names file");
            BufWriter::with_capacity(1024 * 1024, file)
        });

    // Optional: Set up remap names output (for creating correct reference BAM for filter)
    let mut remap_names_writer: Option<BufWriter<File>> =
        config.remap_names_path.as_ref().map(|path| {
            let file = File::create(path).expect("Failed to create remap_names file");
            BufWriter::with_capacity(1024 * 1024, file)
        });

    // Phase 3: Stream BAM and process pairs
    // OPTIMIZATION: Use pre-allocated Record with bam.read() instead of .records() iterator
    // The docs say: "Using the iterator is about 10% slower than the read-based API"
    // We move the record into the buffer when buffering first mates, then allocate fresh
    let t1 = Instant::now();
    eprintln!("Streaming BAM and processing pairs...");

    let mut bam = bam::Reader::from_path(bam_path).context("Failed to open BAM")?;
    bam.set_threads(config.read_threads).ok();

    let header = bam.header().clone();
    let tid_to_name = build_tid_lookup(&header);

    // Pair buffer: read_name -> first-seen mate
    let mut pair_buffer: FxHashMap<Vec<u8>, bam::Record> = FxHashMap::default();
    pair_buffer.reserve(1_000_000);

    // Pre-allocate a single record for reading - avoids per-read allocation
    let mut record = bam::Record::new();

    // Use read() instead of records() iterator for ~10% speedup
    loop {
        match bam.read(&mut record) {
            Some(Ok(())) => {
                stats.total_reads += 1;

                // Skip reads that don't pass baseline filtering:
                // IMPORTANT: Match bam_intersect.rs exactly (unmapped, secondary, supplementary)
                // Do NOT filter on QC fail (0x200) or duplicate (0x400) here because:
                // - bam_filter phase2 adds names to remap set (filters qc/dup on primary read)
                // - bam_filter phase3 writes BOTH mates by name (no filtering!)
                // - bam_intersect filters unmapped, secondary, supplementary ONLY
                // - If one mate is qc_fail but the other overlaps, BOTH go to remap.bam
                // - So we must process qc_fail/duplicate reads to match baseline exactly
                if record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
                    continue;
                }
                // Also check proper_pair like bam_remapper.rs:374 does
                if !record.is_proper_pair() {
                    continue;
                }

                let read_name = record.qname().to_vec();

                // Try to complete a pair
                if let Some(mate) = pair_buffer.remove(&read_name) {
                    // Pair complete - process it
                    stats.pairs_processed += 1;

                    // Ensure read1 is first in template - use references to avoid moving record
                    let (r1, r2): (&bam::Record, &bam::Record) = if record.is_first_in_template() {
                        (&record, &mate)
                    } else {
                        (&mate, &record)
                    };

                    // Check overlaps for both mates - returns ALL overlapping variants
                    let r1_result = check_overlaps(r1, &tid_to_name, &store.trees, &store);
                    let r2_result = check_overlaps(r2, &tid_to_name, &store.trees, &store);

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
                        // process_pair returns ProcessPairResult to distinguish outcomes
                        let overlap_mask =
                            overlap_mask_for_pair(&r1_variants, &r2_variants, &store);
                        increment_overlap_stats(&mut stats, overlap_mask);
                        if config.indel_mode {
                            // INDEL mode: use trim combinations for length preservation
                            let indel_config = IndelConfig {
                                max_indel_size: config.max_indel_size,
                                skip_large_indels: true,
                            };
                            match process_pair_with_trims(
                                r1,
                                r2,
                                &r1_variants,
                                &r2_variants,
                                &store,
                                config,
                                &indel_config,
                                overlap_mask,
                            ) {
                                ProcessPairResult::NeedsRemap(pairs) => {
                                    stats.pairs_with_variants += 1;
                                    // Write read name to remap names file if configured
                                    if let Some(ref mut writer) = remap_names_writer {
                                        writer.write_all(r1.qname()).ok();
                                        writer.write_all(b"\n").ok();
                                    }
                                    for pair in pairs {
                                        tx.send(pair).ok();
                                    }
                                }
                                ProcessPairResult::KeepAsIs => {
                                    stats.pairs_keep_no_flip += 1;
                                    if let Some(ref mut writer) = keep_no_flip_writer {
                                        writer.write_all(r1.qname()).ok();
                                        writer.write_all(b"\n").ok();
                                    }
                                }
                                ProcessPairResult::Unmappable => {
                                    stats.pairs_skipped_unmappable += 1;
                                }
                            }
                        } else {
                            // SNV-only mode: use process_pair with ProcessPairResult
                            match process_pair(
                                r1,
                                r2,
                                &r1_variants,
                                &r2_variants,
                                &store,
                                config,
                                overlap_mask,
                            ) {
                                ProcessPairResult::NeedsRemap(pairs) => {
                                    stats.pairs_with_variants += 1;
                                    // Write read name to remap names file if configured
                                    if let Some(ref mut writer) = remap_names_writer {
                                        writer.write_all(r1.qname()).ok();
                                        writer.write_all(b"\n").ok();
                                    }
                                    for pair in pairs {
                                        tx.send(pair).ok();
                                    }
                                }
                                ProcessPairResult::KeepAsIs => {
                                    // Pair overlaps variants but no sequence change needed
                                    // These reads should be KEPT in final output!
                                    stats.pairs_keep_no_flip += 1;
                                    // Write read name to keep-no-flip file if configured
                                    if let Some(ref mut writer) = keep_no_flip_writer {
                                        writer.write_all(r1.qname()).ok();
                                        writer.write_all(b"\n").ok();
                                    }
                                }
                                ProcessPairResult::Unmappable => {
                                    // Variant in intron/deletion - discard this pair
                                    stats.pairs_skipped_unmappable += 1;
                                }
                            }
                        }
                    }
                    // `mate` is dropped here, `record` is reused for next iteration
                } else {
                    // First mate seen - move record into buffer and allocate new one
                    // This avoids cloning while still allowing record reuse for completed pairs
                    pair_buffer.insert(read_name, record);
                    record = bam::Record::new();
                }

                // Progress reporting
                if stats.total_reads % 10_000_000 == 0 {
                    eprintln!(
                        "  {} reads, {} pairs, {} with variants",
                        stats.total_reads, stats.pairs_processed, stats.pairs_with_variants
                    );
                }
            }
            Some(Err(e)) => return Err(e.into()),
            None => break, // End of file
        }
    }

    stats.orphan_reads = pair_buffer.len();
    stats.bam_stream_ms = t1.elapsed().as_millis() as u64;

    eprintln!("  {} orphan reads (mate not found)", stats.orphan_reads);

    // Flush keep-no-flip writer if configured
    if let Some(mut writer) = keep_no_flip_writer {
        writer
            .flush()
            .context("Failed to flush keep_no_flip_names file")?;
    }

    // Flush remap names writer if configured
    if let Some(mut writer) = remap_names_writer {
        writer.flush().context("Failed to flush remap_names file")?;
    }

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
    eprintln!(
        "  Pairs with variants (needs remap): {}",
        stats.pairs_with_variants
    );
    eprintln!("  Pairs kept (no variants): {}", stats.pairs_kept);
    eprintln!(
        "  Pairs keep-no-flip (variant overlap, no change): {}",
        stats.pairs_keep_no_flip
    );
    eprintln!(
        "  Pairs skipped (unmappable): {}",
        stats.pairs_skipped_unmappable
    );
    eprintln!("  Pairs haplotype failed: {}", stats.pairs_haplotype_failed);
    eprintln!("  Haplotypes written: {}", stats.haplotypes_written);

    eprintln!(
        "  Time: {}ms tree build + {}ms BAM stream",
        stats.tree_build_ms, stats.bam_stream_ms
    );

    Ok(stats)
}

// ============================================================================
// Parallel Chromosome Processing
// ============================================================================
//
// SAFETY NOTE: rust-htslib has a known thread safety issue (GitHub Issue #293):
// - bam::Record contains Rc<HeaderView> which is NOT thread-safe
// - Passing Records between threads causes random segfaults
//
// SAFE PATTERN (used here):
// - Each thread opens its OWN IndexedReader
// - Records are processed entirely within that thread
// - Only primitive data (HaplotypeOutput with Vec<u8>) crosses thread boundaries

/// Process a single chromosome using a per-thread IndexedReader
///
/// SAFETY: This function is designed to be called from rayon parallel iterator.
/// Each thread gets its own BAM reader instance to avoid rust-htslib thread safety issues.
fn process_chromosome(
    bam_path: &str,
    chrom: &str,
    store: &VariantStore,
    tx: &Sender<HaplotypePair>,
    config: &UnifiedConfig,
) -> Result<UnifiedStats> {
    use rust_htslib::bam::Read as BamRead;

    let mut stats = UnifiedStats::default();
    let t0 = Instant::now();

    // CRITICAL: Open a fresh IndexedReader for this thread
    // This avoids the Rc<HeaderView> thread safety bug in rust-htslib
    let mut bam = bam::IndexedReader::from_path(bam_path).context("Failed to open indexed BAM")?;

    // Fetch reads for this chromosome
    bam.fetch(chrom).context("Failed to fetch chromosome")?;

    // Use a few threads for BAM decompression within this worker
    bam.set_threads(2).ok();

    let header = bam.header().clone();
    let tid_to_name = build_tid_lookup(&header);

    // Per-chromosome pair buffer
    let mut pair_buffer: FxHashMap<Vec<u8>, bam::Record> = FxHashMap::default();
    pair_buffer.reserve(100_000); // Smaller per-chromosome

    // Pre-allocated record for reading
    let mut record = bam::Record::new();

    loop {
        match bam.read(&mut record) {
            Some(Ok(())) => {
                stats.total_reads += 1;

                // Apply same filters as sequential version
                if record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
                    continue;
                }
                if !record.is_proper_pair() {
                    continue;
                }

                let read_name = record.qname().to_vec();

                if let Some(mate) = pair_buffer.remove(&read_name) {
                    // Pair complete
                    stats.pairs_processed += 1;

                    let (r1, r2): (&bam::Record, &bam::Record) = if record.is_first_in_template() {
                        (&record, &mate)
                    } else {
                        (&mate, &record)
                    };

                    let r1_result = check_overlaps(r1, &tid_to_name, &store.trees, store);
                    let r2_result = check_overlaps(r2, &tid_to_name, &store.trees, store);

                    let r1_variants = match &r1_result {
                        CheckOverlapResult::Found(v) => v.clone(),
                        CheckOverlapResult::NoOverlaps => Vec::new(),
                    };
                    let r2_variants = match &r2_result {
                        CheckOverlapResult::Found(v) => v.clone(),
                        CheckOverlapResult::NoOverlaps => Vec::new(),
                    };

                    if r1_variants.is_empty() && r2_variants.is_empty() {
                        stats.pairs_kept += 1;
                    } else {
                        let overlap_mask =
                            overlap_mask_for_pair(&r1_variants, &r2_variants, store);
                        increment_overlap_stats(&mut stats, overlap_mask);
                        if config.indel_mode {
                            // INDEL mode: use trim combinations for length preservation
                            let indel_config = IndelConfig {
                                max_indel_size: config.max_indel_size,
                                skip_large_indels: true,
                            };
                            match process_pair_with_trims(
                                r1,
                                r2,
                                &r1_variants,
                                &r2_variants,
                                store,
                                config,
                                &indel_config,
                                overlap_mask,
                            ) {
                                ProcessPairResult::NeedsRemap(pairs) => {
                                    stats.pairs_with_variants += 1;
                                    for pair in pairs {
                                        tx.send(pair).ok();
                                    }
                                }
                                ProcessPairResult::KeepAsIs => {
                                    stats.pairs_keep_no_flip += 1;
                                }
                                ProcessPairResult::Unmappable => {
                                    stats.pairs_skipped_unmappable += 1;
                                }
                            }
                        } else {
                            // SNV-only mode: use process_pair with ProcessPairResult
                            match process_pair(
                                r1,
                                r2,
                                &r1_variants,
                                &r2_variants,
                                store,
                                config,
                                overlap_mask,
                            ) {
                                ProcessPairResult::NeedsRemap(pairs) => {
                                    stats.pairs_with_variants += 1;
                                    for pair in pairs {
                                        // Send pairs to writer thread - only Vec<u8> data crosses threads
                                        tx.send(pair).ok();
                                    }
                                }
                                ProcessPairResult::KeepAsIs => {
                                    // Pair overlaps variants but no sequence change needed
                                    stats.pairs_keep_no_flip += 1;
                                }
                                ProcessPairResult::Unmappable => {
                                    stats.pairs_skipped_unmappable += 1;
                                }
                            }
                        }
                    }
                } else {
                    // First mate - buffer it
                    pair_buffer.insert(read_name, record);
                    record = bam::Record::new();
                }
            }
            Some(Err(e)) => return Err(e.into()),
            None => break,
        }
    }

    stats.orphan_reads = pair_buffer.len();
    stats.bam_stream_ms = t0.elapsed().as_millis() as u64;

    Ok(stats)
}

/// Parallel unified pipeline - processes chromosomes in parallel for 3-8x speedup
///
/// REQUIREMENTS:
/// - BAM must be coordinate-sorted and indexed (.bai file must exist)
/// - Falls back to sequential if BAM index is missing
///
/// THREAD SAFETY:
/// - Each worker thread opens its own IndexedReader (avoids rust-htslib Issue #293)
/// - Records never cross thread boundaries
/// - Only HaplotypePair (paired Vec<u8>) is sent via channel for atomic writing
/// - VariantStore is shared read-only via Arc
pub fn unified_make_reads_parallel(
    bam_path: &str,
    bed_path: &str,
    r1_path: &str,
    r2_path: &str,
    config: &UnifiedConfig,
) -> Result<UnifiedStats> {
    use rayon::prelude::*;

    // Check BAM index exists - fall back to sequential if not
    let bai_path = format!("{}.bai", bam_path);
    if !std::path::Path::new(&bai_path).exists() {
        eprintln!(
            "BAM index not found ({}), falling back to sequential processing",
            bai_path
        );
        return unified_make_reads(bam_path, bed_path, r1_path, r2_path, config);
    }

    // If keep_no_flip_names_path is set, fall back to sequential
    // (parallel version would need thread-safe file writing)
    if config.keep_no_flip_names_path.is_some() {
        eprintln!(
            "keep_no_flip_names_path set, using sequential processing for thread-safe writes"
        );
        return unified_make_reads(bam_path, bed_path, r1_path, r2_path, config);
    }

    // Phase 1: Build variant store (shared, read-only)
    let t0 = Instant::now();
    eprintln!("Building variant store from {}...", bed_path);
    let store = Arc::new(build_variant_store(bed_path)?);
    let tree_build_ms = t0.elapsed().as_millis() as u64;
    eprintln!(
        "  {} chromosomes, {} variants ({}ms)",
        store.trees.len(),
        store.variants.len(),
        tree_build_ms
    );

    // Phase 2: Get chromosome list from BAM header
    let bam = bam::Reader::from_path(bam_path).context("Failed to open BAM")?;
    let chroms: Vec<String> = (0..bam.header().target_count())
        .map(|tid| String::from_utf8_lossy(bam.header().tid2name(tid)).to_string())
        .filter(|c| store.trees.contains_key(c)) // Only chromosomes with variants
        .collect();
    drop(bam);

    eprintln!(
        "Processing {} chromosomes with variants in parallel...",
        chroms.len()
    );

    // Phase 3: Set up output channel and writer thread (sends pairs for atomic writing)
    let (tx, rx): (Sender<HaplotypePair>, Receiver<HaplotypePair>) = bounded(config.channel_buffer);

    let hap_counter = Arc::new(AtomicUsize::new(0));
    let hap_counter_clone = Arc::clone(&hap_counter);

    let r1_owned = r1_path.to_string();
    let r2_owned = r2_path.to_string();
    let sidecar_owned = format!("{}.expected_positions.tsv", r1_owned);
    let compression_threads = config.compression_threads;
    let compress = config.compress_output;
    let writer_handle = thread::spawn(move || {
        fastq_writer_thread(
            rx,
            &r1_owned,
            &r2_owned,
            &sidecar_owned,
            hap_counter_clone,
            compression_threads,
            compress,
        )
    });

    // Phase 4: Process chromosomes in parallel
    // SAFE: Each thread opens its own IndexedReader
    let t1 = Instant::now();
    let bam_path_owned = bam_path.to_string();

    let results: Vec<Result<UnifiedStats>> = chroms
        .par_iter()
        .map(|chrom| {
            // Each thread processes one chromosome with its own reader
            process_chromosome(&bam_path_owned, chrom, &store, &tx, config)
        })
        .collect();

    // Close sender to signal writer thread
    drop(tx);

    // Wait for writer
    writer_handle
        .join()
        .map_err(|_| anyhow::anyhow!("Writer thread panicked"))??;

    // Phase 5: Aggregate stats from all chromosomes
    let mut final_stats = UnifiedStats::default();
    final_stats.tree_build_ms = tree_build_ms;

    for result in results {
        match result {
            Ok(stats) => {
                final_stats = final_stats.merge(stats);
            }
            Err(e) => {
                eprintln!("Warning: Chromosome processing failed: {}", e);
            }
        }
    }

    final_stats.haplotypes_written = hap_counter.load(Ordering::Relaxed);
    final_stats.bam_stream_ms = t1.elapsed().as_millis() as u64;

    eprintln!("Parallel unified pipeline complete:");
    eprintln!("  Total reads: {}", final_stats.total_reads);
    eprintln!("  Pairs processed: {}", final_stats.pairs_processed);
    eprintln!(
        "  Pairs with variants (needs remap): {}",
        final_stats.pairs_with_variants
    );
    eprintln!("  Pairs kept (no variants): {}", final_stats.pairs_kept);
    eprintln!(
        "  Pairs keep-no-flip (variant overlap, no change): {}",
        final_stats.pairs_keep_no_flip
    );
    eprintln!(
        "  Pairs skipped (unmappable): {}",
        final_stats.pairs_skipped_unmappable
    );
    eprintln!("  Haplotypes written: {}", final_stats.haplotypes_written);
    eprintln!(
        "  Time: {}ms tree build + {}ms parallel BAM ({}x potential speedup)",
        final_stats.tree_build_ms,
        final_stats.bam_stream_ms,
        chroms.len().min(rayon::current_num_threads())
    );

    Ok(final_stats)
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
    fn test_generate_wasp_name_extended() {
        let name = generate_wasp_name_extended(b"ERR123456.1000", 10, 20, 3, 5, 257, 16, -2, 4);
        let expected = b"ERR123456.1000_WASP_10_20_3_5_257_16_2_4";
        assert_eq!(name, expected.to_vec());
    }

    #[test]
    fn test_write_fastq_record_forward() {
        let hap = HaplotypeOutput {
            name: b"read/1".to_vec(),
            sequence: b"ACGTN".to_vec(),
            quals: vec![0, 1, 2, 3, 4],
            is_r1: true,
            is_reverse: false,
        };
        let mut out: Vec<u8> = Vec::new();
        let mut seq_buf: Vec<u8> = Vec::new();
        let mut qual_buf: Vec<u8> = Vec::new();
        write_fastq_record(&mut out, &hap, &mut seq_buf, &mut qual_buf).unwrap();
        assert_eq!(out, b"@read/1\nACGTN\n+\n!\"#$%\n".to_vec());
    }

    #[test]
    fn test_write_fastq_record_reverse() {
        let hap = HaplotypeOutput {
            name: b"read/1".to_vec(),
            sequence: b"ACGTN".to_vec(),
            quals: vec![0, 1, 2, 3, 4],
            is_r1: true,
            is_reverse: true,
        };
        let mut out: Vec<u8> = Vec::new();
        let mut seq_buf: Vec<u8> = Vec::new();
        let mut qual_buf: Vec<u8> = Vec::new();
        write_fastq_record(&mut out, &hap, &mut seq_buf, &mut qual_buf).unwrap();
        assert_eq!(out, b"@read/1\nNACGT\n+\n%$#\"!\n".to_vec());
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
