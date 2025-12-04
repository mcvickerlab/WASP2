//! BAM-BED Intersect - Fast read-variant intersection using coitrees
//!
//! Replaces pybedtools intersect with 50-100x faster Rust implementation.
//! Uses coitrees van Emde Boas layout for cache-efficient interval queries.
//!
//! # Performance Optimizations
//! - Index-based metadata: 12-byte tree nodes (vs 112 bytes) = 9x cache efficiency
//! - AVX2 SIMD: ~2x speedup on tree queries (when compiled with target-cpu=native)
//! - SortedQuerent: 2-5x speedup for sorted BAM files
//!
//! # Expected Speedup
//! - 20M reads: 152s (pybedtools) -> ~2-3s (coitrees+AVX2) = 50-75x faster

use anyhow::{Context, Result};
use coitrees::{COITree, COITreeSortedQuerent, IntervalNode, IntervalTree, SortedQuerent};
use rayon::prelude::*;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::{bam, bam::Read as BamRead};
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

// ============================================================================
// Data Structures
// ============================================================================

/// Variant metadata - stored separately from tree for cache efficiency
///
/// Contains all information needed to reconstruct pybedtools output format
#[derive(Clone, Debug)]
pub struct VariantInfo {
    /// Chromosome name (for output)
    pub chrom: String,
    /// Variant start position (0-based)
    pub start: u32,
    /// Variant end position (exclusive)
    pub stop: u32,
    /// Reference allele
    pub ref_allele: String,
    /// Alternate allele
    pub alt_allele: String,
    /// Phased genotype (e.g., "C|T")
    pub genotype: String,
}

/// Per-chromosome interval tree storing indices (not full data)
///
/// Using u32 indices instead of VariantInfo enables:
/// - AVX2 SIMD support (u32 is Copy + Default)
/// - 12-byte nodes vs 112-byte nodes = 9x better cache density
/// - Faster tree traversal for the 90% of reads with no overlaps
pub type VariantTree = COITree<u32, u32>;
pub type ChromTrees = FxHashMap<String, VariantTree>;

/// Combined storage: variants vector + per-chromosome interval trees
///
/// Trees store indices into the variants vector, enabling:
/// - Tiny tree nodes for fast traversal
/// - Full variant data only accessed on matches
pub struct VariantStore {
    /// All variants in a contiguous vector (cache-friendly for sequential access)
    pub variants: Vec<VariantInfo>,
    /// Per-chromosome interval trees with u32 indices as metadata
    pub trees: ChromTrees,
}

// ============================================================================
// Core Functions
// ============================================================================

/// Build variant store from BED file
///
/// # BED Format Expected (from vcf_to_bed output)
/// ```text
/// chrom  start  stop  ref  alt  GT
/// chr10  87400  87401  C   T    C|T
/// ```
///
/// # Arguments
/// * `bed_path` - Path to variant BED file
///
/// # Returns
/// VariantStore with variants vector and per-chromosome trees
///
/// # Performance
/// - Parsing: ~0.5s for 2M variants
/// - Tree construction: ~0.3s for 2M variants
/// - Memory: ~23MB for trees + ~200MB for variant data (2M variants)
pub fn build_variant_store(bed_path: &str) -> Result<VariantStore> {
    let file = File::open(bed_path).context("Failed to open BED file")?;
    let reader = BufReader::with_capacity(1024 * 1024, file); // 1MB buffer

    // Store all variants in a vector
    let mut variants: Vec<VariantInfo> = Vec::new();

    // Collect interval nodes per chromosome (storing indices)
    let mut chrom_intervals: FxHashMap<String, Vec<IntervalNode<u32, u32>>> =
        FxHashMap::default();

    for line in reader.lines() {
        let line = line?;

        // Skip comments and empty lines
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 6 {
            continue; // Skip malformed lines
        }

        let chrom = fields[0].to_string();
        let start = fields[1]
            .parse::<u32>()
            .context("Failed to parse start position")?;
        let stop = fields[2]
            .parse::<u32>()
            .context("Failed to parse stop position")?;

        // Store variant data
        let idx = variants.len() as u32;
        variants.push(VariantInfo {
            chrom: chrom.clone(),
            start,
            stop,
            ref_allele: fields[3].to_string(),
            alt_allele: fields[4].to_string(),
            genotype: fields[5].to_string(),
        });

        // coitrees uses end-inclusive intervals, BED is half-open [start, stop)
        // Store the INDEX as metadata (not the full VariantInfo)
        let node = IntervalNode::new(start as i32, (stop - 1) as i32, idx);

        chrom_intervals
            .entry(chrom)
            .or_insert_with(Vec::new)
            .push(node);
    }

    eprintln!("  Parsed {} variants from BED file", variants.len());

    // Build trees in parallel using rayon
    let chrom_list: Vec<_> = chrom_intervals.into_iter().collect();
    let trees_vec: Vec<_> = chrom_list
        .into_par_iter()
        .map(|(chrom, intervals)| {
            let interval_count = intervals.len();
            let tree = COITree::new(&intervals);
            eprintln!("    {}: {} variants", chrom, interval_count);
            (chrom, tree)
        })
        .collect();

    let trees: ChromTrees = trees_vec.into_iter().collect();

    Ok(VariantStore { variants, trees })
}

/// Intersect BAM reads with variant store, output bedtools-compatible format
///
/// Uses SortedQuerent for 2-5x speedup on sorted BAM files.
/// With AVX2 enabled, tree queries are ~2x faster.
///
/// # Arguments
/// * `bam_path` - Path to BAM file (should be sorted, indexed)
/// * `store` - VariantStore with trees and variant data
/// * `out_path` - Output file path
///
/// # Output Format (matches pybedtools wb=True, bed=True)
/// ```text
/// read_chrom  read_start  read_end  read_name/mate  mapq  strand  \
/// vcf_chrom   vcf_start   vcf_end   ref  alt  GT
/// ```
///
/// # Returns
/// Number of intersections written
///
/// # Performance
/// - Streams BAM: O(1) memory per read
/// - coitrees query: O(log n + k) per read
/// - Index lookup: O(1) per match
pub fn intersect_bam_with_store(
    bam_path: &str,
    store: &VariantStore,
    out_path: &str,
) -> Result<usize> {
    let mut bam = bam::Reader::from_path(bam_path).context("Failed to open BAM")?;

    // Enable multi-threaded BAM decompression (use all available threads)
    let num_threads = rayon::current_num_threads();
    bam.set_threads(num_threads).ok();

    let header = bam.header().clone();

    let out_file = File::create(out_path)?;
    let mut writer = BufWriter::with_capacity(1024 * 1024, out_file); // 1MB buffer

    let mut intersection_count = 0;
    let mut read_count = 0;
    let mut reads_with_overlaps = 0;

    // Build chromosome name lookup
    let mut tid_to_name: Vec<String> = Vec::new();
    for tid in 0..header.target_count() {
        let name = std::str::from_utf8(header.tid2name(tid))
            .unwrap_or("unknown")
            .to_string();
        tid_to_name.push(name);
    }

    // Create SortedQuerent for each chromosome (2-5x faster for sorted BAM)
    // Now works with AVX2 because u32 is Copy + Default!
    let mut querents: FxHashMap<String, COITreeSortedQuerent<u32, u32>> = store
        .trees
        .iter()
        .map(|(k, v)| (k.clone(), SortedQuerent::new(v)))
        .collect();

    for result in bam.records() {
        let read = result?;
        read_count += 1;

        // Skip unmapped, secondary, supplementary
        if read.is_unmapped() || read.is_secondary() || read.is_supplementary() {
            continue;
        }

        // Get chromosome name
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

        // Determine mate number and strand for output
        let mate = if read.is_first_in_template() { 1 } else { 2 };
        let strand = if read.is_reverse() { '-' } else { '+' };
        let mapq = read.mapq();
        let read_name = String::from_utf8_lossy(read.qname());

        let mut has_overlap = false;

        // Query overlapping variants using SortedQuerent + AVX2
        // coitrees uses inclusive intervals, so query [start, end-1]
        querent.query(read_start as i32, read_end as i32 - 1, |node| {
            // Lookup full variant data by index (only on matches!)
            let idx: usize = node.metadata.clone() as usize;
            let info = &store.variants[idx];
            has_overlap = true;

            // Write bedtools-compatible output format
            writeln!(
                writer,
                "{}\t{}\t{}\t{}/{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                chrom,
                read_start,
                read_end,
                read_name,
                mate,
                mapq,
                strand,
                info.chrom,
                info.start,
                info.stop,
                info.ref_allele,
                info.alt_allele,
                info.genotype,
            )
            .ok();

            intersection_count += 1;
        });

        if has_overlap {
            reads_with_overlaps += 1;
        }
    }

    writer.flush()?;

    eprintln!(
        "  Processed {} reads, {} with overlaps, {} total intersections",
        read_count, reads_with_overlaps, intersection_count
    );

    Ok(intersection_count)
}

/// Combined function: build store and intersect in one call
///
/// This is the main entry point from Python.
///
/// # Arguments
/// * `bam_path` - Path to sorted, indexed BAM file
/// * `bed_path` - Path to variant BED file
/// * `out_path` - Output path for intersections
///
/// # Returns
/// Number of intersections found
pub fn intersect_bam_with_variants(
    bam_path: &str,
    bed_path: &str,
    out_path: &str,
) -> Result<usize> {
    eprintln!("Building variant store from {}...", bed_path);
    let store = build_variant_store(bed_path)?;
    eprintln!(
        "  {} chromosomes, {} total variants",
        store.trees.len(),
        store.variants.len()
    );

    eprintln!("Intersecting reads with variants...");
    let count = intersect_bam_with_store(bam_path, &store, out_path)?;
    eprintln!("  {} intersections found", count);

    Ok(count)
}

// ============================================================================
// Multi-Sample Support
// ============================================================================

/// Variant metadata for multi-sample processing
#[derive(Clone, Debug)]
pub struct VariantInfoMulti {
    /// Chromosome name (for output)
    pub chrom: String,
    /// Variant start position (0-based)
    pub start: u32,
    /// Variant end position (exclusive)
    pub stop: u32,
    /// Reference allele
    pub ref_allele: String,
    /// Alternate allele
    pub alt_allele: String,
    /// Per-sample genotypes (e.g., ["A|G", "A|A", "G|T"])
    pub sample_genotypes: Vec<String>,
}

/// Multi-sample variant store
pub struct VariantStoreMulti {
    pub variants: Vec<VariantInfoMulti>,
    pub trees: ChromTrees,
    pub num_samples: usize,
}

/// Build multi-sample variant store from BED file
///
/// # BED Format Expected (multi-sample)
/// ```text
/// chrom  start  stop  ref  alt  GT_S1  GT_S2  GT_S3  ...
/// chr10  87400  87401  C   T    C|T    C|C    T|T
/// ```
pub fn build_variant_store_multi(bed_path: &str, num_samples: usize) -> Result<VariantStoreMulti> {
    let file = File::open(bed_path).context("Failed to open BED file")?;
    let reader = BufReader::with_capacity(1024 * 1024, file);

    let mut variants: Vec<VariantInfoMulti> = Vec::new();
    let mut chrom_intervals: FxHashMap<String, Vec<IntervalNode<u32, u32>>> =
        FxHashMap::default();

    let expected_cols = 5 + num_samples; // chrom, start, stop, ref, alt, GT1, GT2, ...

    for line in reader.lines() {
        let line = line?;

        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < expected_cols {
            continue;
        }

        let chrom = fields[0].to_string();
        let start = fields[1].parse::<u32>().context("Failed to parse start")?;
        let stop = fields[2].parse::<u32>().context("Failed to parse stop")?;

        // Collect sample genotypes
        let mut sample_genotypes = Vec::with_capacity(num_samples);
        for i in 0..num_samples {
            sample_genotypes.push(fields[5 + i].to_string());
        }

        let idx = variants.len() as u32;
        variants.push(VariantInfoMulti {
            chrom: chrom.clone(),
            start,
            stop,
            ref_allele: fields[3].to_string(),
            alt_allele: fields[4].to_string(),
            sample_genotypes,
        });

        let node = IntervalNode::new(start as i32, (stop - 1) as i32, idx);
        chrom_intervals
            .entry(chrom)
            .or_insert_with(Vec::new)
            .push(node);
    }

    eprintln!(
        "  Parsed {} multi-sample variants ({} samples)",
        variants.len(),
        num_samples
    );

    // Build trees in parallel
    let chrom_list: Vec<_> = chrom_intervals.into_iter().collect();
    let trees_vec: Vec<_> = chrom_list
        .into_par_iter()
        .map(|(chrom, intervals)| {
            let tree = COITree::new(&intervals);
            (chrom, tree)
        })
        .collect();

    let trees: ChromTrees = trees_vec.into_iter().collect();

    Ok(VariantStoreMulti {
        variants,
        trees,
        num_samples,
    })
}

/// Intersect BAM with multi-sample variant store
///
/// Output format includes all sample genotypes:
/// ```text
/// chrom  start  end  read/mate  mapq  strand  vcf_chrom  vcf_start  vcf_end  ref  alt  GT_S1  GT_S2  ...
/// ```
pub fn intersect_bam_with_store_multi(
    bam_path: &str,
    store: &VariantStoreMulti,
    out_path: &str,
) -> Result<usize> {
    let mut bam = bam::Reader::from_path(bam_path).context("Failed to open BAM")?;

    let num_threads = rayon::current_num_threads();
    bam.set_threads(num_threads).ok();

    let header = bam.header().clone();

    let out_file = File::create(out_path)?;
    let mut writer = BufWriter::with_capacity(1024 * 1024, out_file);

    let mut intersection_count = 0;
    let mut read_count = 0;

    // Build chromosome name lookup
    let mut tid_to_name: Vec<String> = Vec::new();
    for tid in 0..header.target_count() {
        let name = std::str::from_utf8(header.tid2name(tid))
            .unwrap_or("unknown")
            .to_string();
        tid_to_name.push(name);
    }

    // Create SortedQuerent for each chromosome
    let mut querents: FxHashMap<String, COITreeSortedQuerent<u32, u32>> = store
        .trees
        .iter()
        .map(|(k, v)| (k.clone(), SortedQuerent::new(v)))
        .collect();

    for result in bam.records() {
        let read = result?;
        read_count += 1;

        if read.is_unmapped() || read.is_secondary() || read.is_supplementary() {
            continue;
        }

        let tid = read.tid();
        if tid < 0 || tid as usize >= tid_to_name.len() {
            continue;
        }
        let chrom = &tid_to_name[tid as usize];

        let querent = match querents.get_mut(chrom) {
            Some(q) => q,
            None => continue,
        };

        let read_start = read.pos();
        let read_end = read.reference_end();
        let mate = if read.is_first_in_template() { 1 } else { 2 };
        let strand = if read.is_reverse() { '-' } else { '+' };
        let mapq = read.mapq();
        let read_name = String::from_utf8_lossy(read.qname());

        querent.query(read_start as i32, read_end as i32 - 1, |node| {
            let idx: usize = node.metadata.clone() as usize;
            let info = &store.variants[idx];

            // Write base columns
            write!(
                writer,
                "{}\t{}\t{}\t{}/{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                chrom,
                read_start,
                read_end,
                read_name,
                mate,
                mapq,
                strand,
                info.chrom,
                info.start,
                info.stop,
                info.ref_allele,
                info.alt_allele,
            )
            .ok();

            // Write all sample genotypes
            for gt in &info.sample_genotypes {
                write!(writer, "\t{}", gt).ok();
            }
            writeln!(writer).ok();

            intersection_count += 1;
        });
    }

    writer.flush()?;

    eprintln!(
        "  Processed {} reads, {} intersections ({} samples)",
        read_count, intersection_count, store.num_samples
    );

    Ok(intersection_count)
}

/// Combined multi-sample function: build store and intersect
pub fn intersect_bam_with_variants_multi(
    bam_path: &str,
    bed_path: &str,
    out_path: &str,
    num_samples: usize,
) -> Result<usize> {
    eprintln!(
        "Building multi-sample variant store from {} ({} samples)...",
        bed_path, num_samples
    );
    let store = build_variant_store_multi(bed_path, num_samples)?;
    eprintln!(
        "  {} chromosomes, {} total variants",
        store.trees.len(),
        store.variants.len()
    );

    eprintln!("Intersecting reads with variants (multi-sample)...");
    let count = intersect_bam_with_store_multi(bam_path, &store, out_path)?;
    eprintln!("  {} intersections found", count);

    Ok(count)
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write as IoWrite;
    use tempfile::NamedTempFile;

    #[test]
    fn test_build_variant_store() {
        let mut bed = NamedTempFile::new().unwrap();
        writeln!(bed, "chr1\t100\t101\tA\tG\tA|G").unwrap();
        writeln!(bed, "chr1\t200\t201\tC\tT\tC|T").unwrap();
        writeln!(bed, "chr2\t300\t301\tG\tA\tG|A").unwrap();
        bed.flush().unwrap();

        let store = build_variant_store(bed.path().to_str().unwrap()).unwrap();

        assert_eq!(store.variants.len(), 3, "Should have 3 variants");
        assert_eq!(store.trees.len(), 2, "Should have 2 chromosomes");
        assert!(store.trees.contains_key("chr1"), "Should have chr1");
        assert!(store.trees.contains_key("chr2"), "Should have chr2");
    }

    #[test]
    fn test_build_variant_store_with_comments() {
        let mut bed = NamedTempFile::new().unwrap();
        writeln!(bed, "# This is a comment").unwrap();
        writeln!(bed, "chr1\t100\t101\tA\tG\tA|G").unwrap();
        writeln!(bed, "").unwrap(); // Empty line
        writeln!(bed, "chr1\t200\t201\tC\tT\tC|T").unwrap();
        bed.flush().unwrap();

        let store = build_variant_store(bed.path().to_str().unwrap()).unwrap();

        assert_eq!(store.variants.len(), 2, "Should have 2 variants");
        assert_eq!(store.trees.len(), 1, "Should have 1 chromosome");
        assert!(store.trees.contains_key("chr1"), "Should have chr1");
    }

    #[test]
    fn test_index_based_tree_query() {
        // Build a simple tree with indices
        let variants = vec![
            VariantInfo {
                chrom: "chr1".to_string(),
                start: 100,
                stop: 101,
                ref_allele: "A".to_string(),
                alt_allele: "G".to_string(),
                genotype: "A|G".to_string(),
            },
            VariantInfo {
                chrom: "chr1".to_string(),
                start: 200,
                stop: 201,
                ref_allele: "C".to_string(),
                alt_allele: "T".to_string(),
                genotype: "C|T".to_string(),
            },
        ];

        let intervals: Vec<IntervalNode<u32, u32>> = vec![
            IntervalNode::new(100, 100, 0u32), // Index 0
            IntervalNode::new(200, 200, 1u32), // Index 1
        ];

        let tree: COITree<u32, u32> = COITree::new(&intervals);

        // Query that should hit first variant
        let mut found_indices: Vec<u32> = Vec::new();
        tree.query(50, 150, |node| {
            found_indices.push(*node.metadata);
        });
        assert_eq!(found_indices.len(), 1);
        assert_eq!(found_indices[0], 0);
        assert_eq!(variants[found_indices[0] as usize].ref_allele, "A");

        // Query that should hit both variants
        found_indices.clear();
        tree.query(50, 250, |node| {
            found_indices.push(*node.metadata);
        });
        assert_eq!(found_indices.len(), 2);

        // Query that should hit nothing
        found_indices.clear();
        tree.query(300, 400, |node| {
            found_indices.push(*node.metadata);
        });
        assert_eq!(found_indices.len(), 0);
    }

    #[test]
    fn test_sorted_querent_with_indices() {
        // Verify SortedQuerent works with u32 indices
        let intervals: Vec<IntervalNode<u32, u32>> = vec![
            IntervalNode::new(100, 100, 0u32),
            IntervalNode::new(200, 200, 1u32),
            IntervalNode::new(300, 300, 2u32),
        ];

        let tree: COITree<u32, u32> = COITree::new(&intervals);
        let mut querent: COITreeSortedQuerent<u32, u32> = SortedQuerent::new(&tree);

        // Sorted queries (simulating sorted BAM)
        let mut count = 0;
        querent.query(50, 150, |_| count += 1);
        assert_eq!(count, 1);

        count = 0;
        querent.query(150, 250, |_| count += 1);
        assert_eq!(count, 1);

        count = 0;
        querent.query(250, 350, |_| count += 1);
        assert_eq!(count, 1);
    }
}
