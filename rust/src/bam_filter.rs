//! BAM Variant Filter - Fast BAM splitting by variant overlap
//!
//! Replaces Python process_bam() with 4-5x faster Rust implementation.
//! Uses existing coitrees infrastructure from bam_intersect.rs.
//!
//! # Performance
//! - Current Python/samtools: ~450s for 56M reads
//! - Target Rust: ~100s (4-5x faster)
//!
//! # Algorithm
//! 1. Build variant interval tree from BED (reuse bam_intersect)
//! 2. Stream BAM, collect read names overlapping variants
//! 3. Stream BAM again, split to remap/keep based on name membership

use anyhow::{Context, Result};
use coitrees::{COITreeSortedQuerent, SortedQuerent};
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::{bam, bam::Read as BamRead};
use rustc_hash::{FxHashMap, FxHashSet};
use std::time::Instant;

use crate::bam_intersect::{build_variant_store, VariantStore};

// ============================================================================
// Data Structures
// ============================================================================

/// Statistics returned from filtering operation
#[derive(Debug, Clone, Default)]
pub struct FilterStats {
    /// Total reads processed
    pub total_reads: usize,
    /// Reads sent to remap BAM (overlapping variants or their mates)
    pub remap_reads: usize,
    /// Reads sent to keep BAM (no variant overlap)
    pub keep_reads: usize,
    /// Unique read names overlapping variants
    pub unique_remap_names: usize,
    /// Time spent in each phase (ms)
    pub phase1_ms: u64,
    pub phase2_ms: u64,
    pub phase3_ms: u64,
}

/// Configuration for BAM filtering
#[derive(Debug, Clone)]
pub struct FilterConfig {
    /// Number of threads for BAM reading
    pub read_threads: usize,
    /// Number of threads for BAM writing
    pub write_threads: usize,
    /// Whether input is paired-end
    pub is_paired: bool,
}

impl Default for FilterConfig {
    fn default() -> Self {
        Self {
            read_threads: 4,
            write_threads: 4,
            is_paired: true,
        }
    }
}

// ============================================================================
// Helper Functions
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

// ============================================================================
// Core Algorithm
// ============================================================================

/// Phase 2: Stream BAM, find reads overlapping variants, collect their names
///
/// # Key optimizations
/// - Parallel BAM decompression (rust-htslib thread pool)
/// - SortedQuerent for cache-efficient overlap queries on sorted BAM
/// - FxHashSet for O(1) membership (vs Python set)
fn phase2_collect_remap_names(
    bam_path: &str,
    store: &VariantStore,
    config: &FilterConfig,
) -> Result<FxHashSet<Vec<u8>>> {
    let mut bam = bam::Reader::from_path(bam_path).context("Failed to open BAM for phase 2")?;

    // Enable multi-threaded BAM decompression (use all available threads)
    let num_threads = config.read_threads.min(rayon::current_num_threads());
    bam.set_threads(num_threads).ok();

    let header = bam.header().clone();
    let tid_to_name = build_tid_lookup(&header);

    // Pre-allocate for expected ~10% overlap rate
    // For 56M reads with ~10% overlap, ~5.6M unique names
    let mut remap_names: FxHashSet<Vec<u8>> = FxHashSet::default();
    remap_names.reserve(2_000_000);

    // Create SortedQuerent per chromosome (2-5x faster for sorted BAM)
    let mut querents: FxHashMap<String, COITreeSortedQuerent<u32, u32>> = store
        .trees
        .iter()
        .map(|(k, v)| (k.clone(), SortedQuerent::new(v)))
        .collect();

    let mut processed = 0usize;
    let mut overlapping = 0usize;

    // Use read() with pre-allocated Record instead of records() iterator for better performance
    let mut read = bam::Record::new();
    while let Some(result) = bam.read(&mut read) {
        result?;
        processed += 1;

        // Skip unmapped, secondary, supplementary, QC fail, duplicate
        // Flags: 0x4=unmapped, 0x100=secondary, 0x800=supplementary, 0x200=QC fail, 0x400=duplicate
        if read.flags() & (0x4 | 0x100 | 0x800 | 0x200 | 0x400) != 0 {
            continue;
        }

        let tid = read.tid();
        if tid < 0 || tid as usize >= tid_to_name.len() {
            continue;
        }

        let chrom = &tid_to_name[tid as usize];

        // Skip if no variants on this chromosome
        let querent = match querents.get_mut(chrom) {
            Some(q) => q,
            None => continue,
        };

        // Read coordinates (0-based, half-open)
        let read_start = read.pos();
        let read_end = read.reference_end();

        // Check for overlap with any variant
        let mut has_overlap = false;
        querent.query(read_start as i32, read_end as i32 - 1, |_| {
            has_overlap = true;
        });

        if has_overlap {
            // Store read name (as bytes, no String allocation)
            remap_names.insert(read.qname().to_vec());
            overlapping += 1;
        }
    }

    eprintln!(
        "  Phase 2: {} reads processed, {} overlapping, {} unique names",
        processed,
        overlapping,
        remap_names.len()
    );

    Ok(remap_names)
}

/// Phase 3: Stream BAM, split to remap/keep based on read name membership
///
/// # Key optimizations
/// - Single pass through BAM
/// - FxHashSet O(1) membership check
/// - Parallel BGZF compression for both output files
fn phase3_split_bam(
    bam_path: &str,
    remap_names: &FxHashSet<Vec<u8>>,
    remap_bam_path: &str,
    keep_bam_path: &str,
    config: &FilterConfig,
) -> Result<(usize, usize)> {
    let mut bam = bam::Reader::from_path(bam_path).context("Failed to open BAM for phase 3")?;

    // Enable multi-threaded BAM reading (use all available threads)
    bam.set_threads(config.read_threads.min(rayon::current_num_threads())).ok();

    // Convert HeaderView to Header for writer
    let header = bam::Header::from_template(bam.header());

    // Create writers with parallel compression (use all available threads, fastest compression)
    let mut remap_writer =
        bam::Writer::from_path(remap_bam_path, &header, bam::Format::Bam)
            .context("Failed to create remap BAM writer")?;
    remap_writer.set_threads(config.write_threads.min(rayon::current_num_threads())).ok();
    remap_writer.set_compression_level(bam::CompressionLevel::Fastest).ok();

    let mut keep_writer =
        bam::Writer::from_path(keep_bam_path, &header, bam::Format::Bam)
            .context("Failed to create keep BAM writer")?;
    keep_writer.set_threads(config.write_threads.min(rayon::current_num_threads())).ok();
    keep_writer.set_compression_level(bam::CompressionLevel::Fastest).ok();

    let mut remap_count = 0usize;
    let mut keep_count = 0usize;

    // Use read() with pre-allocated Record instead of records() iterator for better performance
    let mut record = bam::Record::new();
    while let Some(result) = bam.read(&mut record) {
        result?;

        // For paired-end: if THIS read's name is in the set, BOTH mates go to remap
        // This ensures pairs stay together
        if remap_names.contains(record.qname()) {
            remap_writer.write(&record)?;
            remap_count += 1;
        } else {
            keep_writer.write(&record)?;
            keep_count += 1;
        }
    }

    eprintln!(
        "  Phase 3: {} remap, {} keep ({} total)",
        remap_count,
        keep_count,
        remap_count + keep_count
    );

    Ok((remap_count, keep_count))
}

/// Filter BAM by variant overlap - main entry point
///
/// Replaces process_bam() from intersect_variant_data.py
///
/// # Arguments
/// * `bam_path` - Input BAM file (should be coordinate-sorted)
/// * `bed_path` - Variant BED file (from vcf_to_bed)
/// * `remap_bam_path` - Output BAM for reads needing remapping
/// * `keep_bam_path` - Output BAM for reads not needing remapping
/// * `is_paired` - Whether reads are paired-end
/// * `threads` - Number of threads to use
///
/// # Returns
/// Tuple of (remap_count, keep_count, unique_names)
pub fn filter_bam_by_variants(
    bam_path: &str,
    bed_path: &str,
    remap_bam_path: &str,
    keep_bam_path: &str,
    is_paired: bool,
    threads: usize,
) -> Result<FilterStats> {
    let config = FilterConfig {
        read_threads: threads,
        write_threads: threads,
        is_paired,
    };

    let mut stats = FilterStats::default();

    // Phase 1: Build variant store (reuse from bam_intersect)
    let t0 = Instant::now();
    eprintln!("Phase 1: Building variant store from {}...", bed_path);
    let store = build_variant_store(bed_path)?;
    stats.phase1_ms = t0.elapsed().as_millis() as u64;
    eprintln!(
        "  {} chromosomes, {} variants ({}ms)",
        store.trees.len(),
        store.variants.len(),
        stats.phase1_ms
    );

    // Phase 2: Collect overlapping read names
    let t1 = Instant::now();
    eprintln!("Phase 2: Collecting overlapping read names...");
    let remap_names = phase2_collect_remap_names(bam_path, &store, &config)?;
    stats.phase2_ms = t1.elapsed().as_millis() as u64;
    stats.unique_remap_names = remap_names.len();
    eprintln!(
        "  {} unique read names to remap ({}ms)",
        remap_names.len(),
        stats.phase2_ms
    );

    // Phase 3: Split BAM
    let t2 = Instant::now();
    eprintln!("Phase 3: Splitting BAM into remap/keep...");
    let (remap_count, keep_count) =
        phase3_split_bam(bam_path, &remap_names, remap_bam_path, keep_bam_path, &config)?;
    stats.phase3_ms = t2.elapsed().as_millis() as u64;
    stats.remap_reads = remap_count;
    stats.keep_reads = keep_count;
    stats.total_reads = remap_count + keep_count;

    let total_ms = stats.phase1_ms + stats.phase2_ms + stats.phase3_ms;
    eprintln!(
        "✅ Filter complete: {} remap, {} keep, {} unique names",
        remap_count, keep_count, remap_names.len()
    );
    eprintln!(
        "   Total time: {}ms (phase1: {}ms, phase2: {}ms, phase3: {}ms)",
        total_ms, stats.phase1_ms, stats.phase2_ms, stats.phase3_ms
    );

    Ok(stats)
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write as IoWrite;
    use tempfile::{tempdir, NamedTempFile};

    /// Create a minimal BED file for testing
    fn create_test_bed() -> NamedTempFile {
        let mut bed = NamedTempFile::new().unwrap();
        writeln!(bed, "chr1\t100\t101\tA\tG\tA|G").unwrap();
        writeln!(bed, "chr1\t200\t201\tC\tT\tC|T").unwrap();
        writeln!(bed, "chr1\t300\t301\tG\tA\tG|A").unwrap();
        bed.flush().unwrap();
        bed
    }

    #[test]
    fn test_build_tid_lookup() {
        // This would need a real BAM file to test properly
        // For now, just verify the function signature works
    }

    #[test]
    fn test_filter_config_default() {
        let config = FilterConfig::default();
        assert_eq!(config.read_threads, 4);
        assert_eq!(config.write_threads, 4);
        assert!(config.is_paired);
    }

    #[test]
    fn test_filter_stats_default() {
        let stats = FilterStats::default();
        assert_eq!(stats.total_reads, 0);
        assert_eq!(stats.remap_reads, 0);
        assert_eq!(stats.keep_reads, 0);
        assert_eq!(stats.unique_remap_names, 0);
    }
}
