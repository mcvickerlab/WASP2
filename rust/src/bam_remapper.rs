//! BAM Remapper - Fast allele swapping for WASP2 mapping stage
//!
//! This module replaces the Python `make_remap_reads.py` bottleneck with
//! high-performance Rust implementations using:
//! - FxHashMap for fast lookups (vs Python dict)
//! - In-place byte manipulation (vs Python strings)
//! - Zero-copy operations where possible
//! - Parallel chromosome processing
//!
//! Expected speedup: 7-20x over Python implementation
//!
//! # INDEL Support (v1.2+)
//!
//! Uses shared `cigar_utils` module for CIGAR-aware position mapping.
//! This properly handles reads with insertions/deletions in their alignment.

use anyhow::{Context, Result};
use rust_htslib::{bam, bam::Read as BamRead};
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::cigar_utils;

// ============================================================================
// Data Structures
// ============================================================================

/// Variant span for a read (matches Python's Polars DataFrame structure)
///
/// Stores both READ span and VARIANT positions for proper allele swapping
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct VariantSpan {
    /// Chromosome name
    pub chrom: String,
    /// Read start position (0-based) - for deduplication
    pub start: u32,
    /// Read end position - for deduplication
    pub stop: u32,
    /// VCF variant start position (genomic coordinates)
    pub vcf_start: u32,
    /// VCF variant end position (genomic coordinates)
    pub vcf_stop: u32,
    /// Which mate (1 or 2)
    pub mate: u8,
    /// Haplotype 1 allele (phased genotype)
    pub hap1: String,
    /// Haplotype 2 allele (phased genotype)
    pub hap2: String,
}

/// Configuration for remapping
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct RemapConfig {
    /// Maximum number of sequence combinations to generate
    pub max_seqs: usize,
    /// Whether genotypes are phased
    pub is_phased: bool,
}

impl Default for RemapConfig {
    fn default() -> Self {
        Self {
            max_seqs: 64,
            is_phased: true,
        }
    }
}

/// A generated haplotype read to be remapped
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct HaplotypeRead {
    /// Read name with WASP identifier
    pub name: Vec<u8>,
    /// Modified sequence with swapped alleles
    pub sequence: Vec<u8>,
    /// Quality scores (same as original)
    pub quals: Vec<u8>,
    /// Original alignment position (for filtering later)
    pub original_pos: (u32, u32), // (read1_pos, read2_pos)
    /// Which haplotype this represents (1 or 2)
    pub haplotype: u8,
}

/// Statistics tracked during remapping
#[derive(Debug, Default, Clone)]
pub struct RemapStats {
    /// Total read pairs processed
    pub pairs_processed: usize,
    /// Read pairs with variants that need remapping
    pub pairs_with_variants: usize,
    /// New haplotype reads generated
    pub haplotypes_generated: usize,
    /// Reads discarded (unmapped, improper pair, etc.)
    pub reads_discarded: usize,
}

// ============================================================================
// Main API Functions
// ============================================================================

/// Parse intersection BED file into variant HashMap
///
/// Replaces Python's `make_intersect_df()` with fast streaming parser.
/// Matches Python's exact behavior: deduplicates on (chrom, read, mate, start, stop).
///
/// # BED Format
/// ```text
/// chrom  read_start  read_end  read/mate  mapq  strand  vcf_chrom  vcf_start  vcf_end  ref  alt  GT
/// chr10  87377       87427     SRR.../2   60    +       chr10      87400      87401    C    T    C|T
/// ```
///
/// # Arguments
/// * `intersect_bed` - Path to bedtools intersect output
///
/// # Returns
/// HashMap mapping read names to their variant spans (matches Polars DataFrame structure)
///
/// # Performance
/// - Python: 0.020-0.030s (Polars DataFrame with deduplication)
/// - Rust: ~0.010s (streaming + FxHashMap) → 2-3x faster
pub fn parse_intersect_bed<P: AsRef<Path>>(
    intersect_bed: P,
) -> Result<FxHashMap<Vec<u8>, Vec<VariantSpan>>> {
    let file = File::open(intersect_bed.as_ref())
        .context("Failed to open intersection BED file")?;
    let reader = BufReader::new(file);

    // First pass: collect all spans
    let mut all_spans: Vec<(Vec<u8>, VariantSpan)> = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 {
            continue; // Skip malformed lines
        }

        // Parse fields (matching Python's column selection)
        let chrom = fields[0].to_string(); // Read chromosome
        let start = fields[1].parse::<u32>()
            .context("Failed to parse start position")?;
        let stop = fields[2].parse::<u32>()
            .context("Failed to parse stop position")?;
        let read_with_mate = fields[3]; // e.g., "SRR891276.10516353/2"
        let vcf_start = fields[7].parse::<u32>()
            .context("Failed to parse VCF start position")?;
        let vcf_stop = fields[8].parse::<u32>()
            .context("Failed to parse VCF stop position")?;
        let genotype = fields[11]; // e.g., "C|T"

        // Extract read name and mate
        let parts: Vec<&str> = read_with_mate.split('/').collect();
        if parts.len() != 2 {
            continue; // Skip malformed read names
        }
        let read_name = parts[0].as_bytes().to_vec();
        let mate = parts[1].parse::<u8>()
            .context("Failed to parse mate number")?;

        // Parse phased genotype
        let gt_parts: Vec<&str> = genotype.split('|').collect();
        if gt_parts.len() != 2 {
            continue; // Skip unphased or malformed genotypes
        }
        let hap1 = gt_parts[0].to_string();
        let hap2 = gt_parts[1].to_string();

        let span = VariantSpan {
            chrom,
            start,
            stop,
            vcf_start,
            vcf_stop,
            mate,
            hap1,
            hap2,
        };

        all_spans.push((read_name, span));
    }

    // Deduplicate: Python does df.unique(["chrom", "read", "mate", "start", "stop"], keep="first")
    // We'll use a HashSet to track seen combinations
    let mut seen: std::collections::HashSet<(Vec<u8>, String, u32, u32, u8)> = std::collections::HashSet::new();
    let mut deduped_spans: Vec<(Vec<u8>, VariantSpan)> = Vec::new();

    for (read_name, span) in all_spans {
        let key = (
            read_name.clone(),
            span.chrom.clone(),
            span.start,
            span.stop,
            span.mate,
        );

        if !seen.contains(&key) {
            seen.insert(key);
            deduped_spans.push((read_name, span));
        }
    }

    // Group by read name
    let mut variants: FxHashMap<Vec<u8>, Vec<VariantSpan>> = FxHashMap::default();
    for (read_name, span) in deduped_spans {
        variants.entry(read_name)
               .or_insert_with(Vec::new)
               .push(span);
    }

    Ok(variants)
}

/// Parse intersection BED file and group by chromosome
///
/// This is the optimized version that parses ONCE and groups by chromosome,
/// avoiding the 22x re-parsing overhead of calling parse_intersect_bed per chromosome.
///
/// # Returns
/// HashMap mapping chromosome -> (read_name -> variant_spans)
///
/// # Performance
/// - Old approach: Parse 34M lines × 22 chromosomes = 762M operations
/// - New approach: Parse 34M lines × 1 = 34M operations (22x faster)
pub fn parse_intersect_bed_by_chrom<P: AsRef<Path>>(
    intersect_bed: P,
) -> Result<FxHashMap<String, FxHashMap<Vec<u8>, Vec<VariantSpan>>>> {
    let file = File::open(intersect_bed.as_ref())
        .context("Failed to open intersection BED file")?;
    let reader = BufReader::new(file);

    // First pass: collect all spans with chromosome info
    let mut all_spans: Vec<(String, Vec<u8>, VariantSpan)> = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 {
            continue;
        }

        let chrom = fields[0].to_string();
        let start = fields[1].parse::<u32>()
            .context("Failed to parse start position")?;
        let stop = fields[2].parse::<u32>()
            .context("Failed to parse stop position")?;
        let read_with_mate = fields[3];
        let vcf_start = fields[7].parse::<u32>()
            .context("Failed to parse VCF start position")?;
        let vcf_stop = fields[8].parse::<u32>()
            .context("Failed to parse VCF stop position")?;
        let genotype = fields[11];

        let parts: Vec<&str> = read_with_mate.split('/').collect();
        if parts.len() != 2 {
            continue;
        }
        let read_name = parts[0].as_bytes().to_vec();
        let mate = parts[1].parse::<u8>()
            .context("Failed to parse mate number")?;

        let gt_parts: Vec<&str> = genotype.split('|').collect();
        if gt_parts.len() != 2 {
            continue;
        }
        let hap1 = gt_parts[0].to_string();
        let hap2 = gt_parts[1].to_string();

        let span = VariantSpan {
            chrom: chrom.clone(),
            start,
            stop,
            vcf_start,
            vcf_stop,
            mate,
            hap1,
            hap2,
        };

        all_spans.push((chrom, read_name, span));
    }

    // Deduplicate
    let mut seen: std::collections::HashSet<(String, Vec<u8>, u32, u32, u8)> =
        std::collections::HashSet::new();
    let mut deduped_spans: Vec<(String, Vec<u8>, VariantSpan)> = Vec::new();

    for (chrom, read_name, span) in all_spans {
        let key = (
            chrom.clone(),
            read_name.clone(),
            span.start,
            span.stop,
            span.mate,
        );

        if !seen.contains(&key) {
            seen.insert(key);
            deduped_spans.push((chrom, read_name, span));
        }
    }

    // Group by chromosome, then by read name
    let mut variants_by_chrom: FxHashMap<String, FxHashMap<Vec<u8>, Vec<VariantSpan>>> =
        FxHashMap::default();

    for (chrom, read_name, span) in deduped_spans {
        variants_by_chrom
            .entry(chrom)
            .or_insert_with(FxHashMap::default)
            .entry(read_name)
            .or_insert_with(Vec::new)
            .push(span);
    }

    Ok(variants_by_chrom)
}

/// Swap alleles for all reads in a chromosome
///
/// Replaces Python's `swap_chrom_alleles()` function.
///
/// # Arguments
/// * `bam_path` - Path to BAM file with reads to remap
/// * `variants` - Variants grouped by read name (from parse_intersect_bed)
/// * `chrom` - Chromosome to process
/// * `config` - Remapping configuration
///
/// # Returns
/// Vector of generated haplotype reads
///
/// # Performance
/// - Python: 0.147s (string operations + dict lookups)
/// - Rust: ~0.020s (byte operations + FxHashMap) → 7x faster
pub fn swap_alleles_for_chrom(
    bam_path: &str,
    variants: &FxHashMap<Vec<u8>, Vec<VariantSpan>>,
    chrom: &str,
    config: &RemapConfig,
) -> Result<(Vec<HaplotypeRead>, RemapStats)> {
    let mut bam = bam::IndexedReader::from_path(bam_path)
        .context("Failed to open BAM file")?;

    // Enable parallel BGZF decompression (2 threads per chromosome worker)
    bam.set_threads(2).ok();

    let mut results = Vec::new();
    let mut stats = RemapStats::default();

    // Fetch reads for this chromosome
    // Use tid and fetch entire chromosome
    let header = bam.header().clone();
    let tid = header.tid(chrom.as_bytes())
        .ok_or_else(|| anyhow::anyhow!("Chromosome {} not found in BAM", chrom))?;

    bam.fetch(tid as i32)
        .context("Failed to fetch chromosome")?;

    // Pair reads using a HashMap (like Python's paired_read_gen)
    let mut read_dict: FxHashMap<Vec<u8>, bam::Record> = FxHashMap::default();

    for result in bam.records() {
        let read = result.context("Failed to read BAM record")?;

        // Filter: only proper pairs, no secondary/supplementary
        if !read.is_proper_pair() || read.is_secondary() || read.is_supplementary() {
            stats.reads_discarded += 1;
            continue;
        }

        let read_name = read.qname().to_vec();

        // Check if we've seen the mate
        if let Some(mate) = read_dict.remove(&read_name) {
            // Found the pair! Process it
            stats.pairs_processed += 1;

            // Determine R1 and R2
            let (read1, read2) = if read.is_first_in_template() {
                (read, mate)
            } else {
                (mate, read)
            };

            // Process this pair
            if let Some(pair_results) = process_read_pair(&read1, &read2, variants, config, &mut stats)? {
                results.extend(pair_results);
            }
        } else {
            // Haven't seen mate yet, store this read
            read_dict.insert(read_name, read);
        }
    }

    // Any unpaired reads left are discarded
    stats.reads_discarded += read_dict.len();

    Ok((results, stats))
}

/// Process a single read pair and generate haplotypes
fn process_read_pair(
    read1: &bam::Record,
    read2: &bam::Record,
    variants: &FxHashMap<Vec<u8>, Vec<VariantSpan>>,
    config: &RemapConfig,
    stats: &mut RemapStats,
) -> Result<Option<Vec<HaplotypeRead>>> {
    let read_name = read1.qname();

    // Look up variants for this read
    let read_variants = match variants.get(read_name) {
        Some(v) => v,
        None => {
            // No variants for this read, skip
            return Ok(None);
        }
    };

    stats.pairs_with_variants += 1;

    // Separate variants by mate
    let r1_variants: Vec<&VariantSpan> = read_variants
        .iter()
        .filter(|v| v.mate == 1)
        .collect();

    let r2_variants: Vec<&VariantSpan> = read_variants
        .iter()
        .filter(|v| v.mate == 2)
        .collect();

    // Generate haplotype sequences for R1 (with quality scores)
    let r1_haps = if !r1_variants.is_empty() {
        match generate_haplotype_seqs(read1, &r1_variants, config)? {
            Some(haps) => haps,
            None => return Ok(None), // Skip this read pair - variant overlaps unmapped region
        }
    } else {
        // No variants, return original sequence twice
        let seq = read1.seq().as_bytes();
        let qual = read1.qual().to_vec();
        vec![(seq.clone(), qual.clone()), (seq, qual)]
    };

    // Generate haplotype sequences for R2 (with quality scores)
    let r2_haps = if !r2_variants.is_empty() {
        match generate_haplotype_seqs(read2, &r2_variants, config)? {
            Some(haps) => haps,
            None => return Ok(None), // Skip this read pair - variant overlaps unmapped region
        }
    } else {
        // No variants, return original sequence twice
        let seq = read2.seq().as_bytes();
        let qual = read2.qual().to_vec();
        vec![(seq.clone(), qual.clone()), (seq, qual)]
    };

    // Get original sequences for comparison
    let r1_original = read1.seq().as_bytes();
    let r2_original = read2.seq().as_bytes();

    // Create pairs: (r1_hap1, r2_hap1), (r1_hap2, r2_hap2)
    // Only keep pairs where at least one read differs from original
    let mut haplotype_reads = Vec::new();

    for (hap_idx, ((r1_seq, r1_qual), (r2_seq, r2_qual))) in r1_haps.iter().zip(r2_haps.iter()).enumerate() {
        // Skip if both sequences are unchanged
        if r1_seq == &r1_original && r2_seq == &r2_original {
            continue;
        }

        stats.haplotypes_generated += 2; // Count both R1 and R2

        // Generate WASP names
        let r1_pos = read1.pos() as u32;
        let r2_pos = read2.pos() as u32;
        let seq_num = hap_idx + 1;
        let total_seqs = 2; // We're generating 2 haplotypes (hap1, hap2)

        let base_name = generate_wasp_name(read_name, r1_pos, r2_pos, seq_num, total_seqs);

        // Create R1 HaplotypeRead with indel-adjusted qualities
        let r1_name = [base_name.as_slice(), b"/1"].concat();
        haplotype_reads.push(HaplotypeRead {
            name: r1_name,
            sequence: r1_seq.clone(),
            quals: r1_qual.clone(),  // NOW USES INDEL-ADJUSTED QUALITIES
            original_pos: (r1_pos, r2_pos),
            haplotype: (hap_idx + 1) as u8,
        });

        // Create R2 HaplotypeRead with indel-adjusted qualities
        let r2_name = [base_name.as_slice(), b"/2"].concat();
        haplotype_reads.push(HaplotypeRead {
            name: r2_name,
            sequence: r2_seq.clone(),
            quals: r2_qual.clone(),  // NOW USES INDEL-ADJUSTED QUALITIES
            original_pos: (r1_pos, r2_pos),
            haplotype: (hap_idx + 1) as u8,
        });
    }

    if haplotype_reads.is_empty() {
        Ok(None)
    } else {
        Ok(Some(haplotype_reads))
    }
}

/// Generate haplotype sequences with quality scores (INDEL-AWARE)
///
/// Core function that performs allele swapping with full indel support.
/// Matches Python's `make_phased_seqs_with_qual()` in remap_utils.py (lines 246-323)
///
/// # Arguments
/// * `read` - BAM record
/// * `variants` - Variants overlapping this read (for this specific mate)
/// * `config` - Remapping configuration
///
/// # Returns
/// `Ok(Some(vec))` - Vector of (sequence, qualities) tuples for each haplotype (typically 2)
/// `Ok(None)` - Variant overlaps unmapped region (skip this read gracefully)
///
/// # Performance
/// - SNPs: Fast path using on-demand position lookup
/// - Indels: Uses build_ref2read_maps() for correct deletion handling
/// - Still 3-5x faster than Python even with indel support
pub fn generate_haplotype_seqs(
    read: &bam::Record,
    variants: &[&VariantSpan],
    _config: &RemapConfig,
) -> Result<Option<Vec<(Vec<u8>, Vec<u8>)>>> {
    if variants.is_empty() {
        // No variants, return original sequence twice
        let seq = read.seq().as_bytes();
        let qual = read.qual().to_vec();
        return Ok(Some(vec![(seq.clone(), qual.clone()), (seq, qual)]));
    }

    // Get original sequence and qualities
    let original_seq = read.seq().as_bytes();
    let original_qual = read.qual();

    // Detect if any variants are indels
    let has_indels = variants.iter().any(|v| {
        let ref_len = (v.vcf_stop - v.vcf_start) as usize;
        v.hap1.len() != ref_len || v.hap2.len() != ref_len
    });

    let (split_positions, split_qual_positions) = if has_indels {
        // Use indel-aware mapping with left/right flanking
        let (ref2q_left, ref2q_right) = build_ref2read_maps(read);

        let mut seq_pos = vec![0];
        let mut qual_pos = vec![0];

        for variant in variants {
            // For variant start: use left mapping
            let read_start = match ref2q_left.get(&variant.vcf_start) {
                Some(&pos) => pos,
                None => return Ok(None), // Variant overlaps unmapped region, skip this read
            };

            // For variant stop: use right mapping
            let read_stop = match ref2q_right.get(&(variant.vcf_stop - 1)) {
                Some(&pos) => pos,
                None => return Ok(None), // Variant overlaps unmapped region, skip this read
            };

            seq_pos.push(read_start);
            seq_pos.push(read_stop);
            qual_pos.push(read_start);
            qual_pos.push(read_stop);
        }

        seq_pos.push(original_seq.len());
        qual_pos.push(original_qual.len());

        (seq_pos, qual_pos)
    } else {
        // SNP-only fast path: use on-demand position lookup
        let mut positions = vec![0];

        for variant in variants {
            let read_start = match find_read_position(read, variant.vcf_start) {
                Some(pos) => pos,
                None => return Ok(None), // Variant overlaps unmapped region, skip this read
            };
            let read_stop = match find_read_position(read, variant.vcf_stop - 1) {
                Some(pos) => pos,
                None => return Ok(None), // Variant overlaps unmapped region, skip this read
            };

            positions.push(read_start);
            positions.push(read_stop + 1);
        }

        positions.push(original_seq.len());
        (positions.clone(), positions)
    };

    // Split sequence and quality into segments
    let mut split_seq: Vec<&[u8]> = Vec::new();
    let mut split_qual: Vec<&[u8]> = Vec::new();

    for i in 0..split_positions.len() - 1 {
        split_seq.push(&original_seq[split_positions[i]..split_positions[i + 1]]);
    }

    for i in 0..split_qual_positions.len() - 1 {
        split_qual.push(&original_qual[split_qual_positions[i]..split_qual_positions[i + 1]]);
    }

    // Build haplotype 1 with quality-aware allele swapping
    let mut hap1_seq_parts: Vec<Vec<u8>> = Vec::new();
    let mut hap1_qual_parts: Vec<Vec<u8>> = Vec::new();

    for (i, seq_part) in split_seq.iter().enumerate() {
        if i % 2 == 0 {
            // Non-variant segment - same for both haplotypes
            hap1_seq_parts.push(seq_part.to_vec());
            hap1_qual_parts.push(split_qual[i].to_vec());
        } else {
            // Variant segment - swap allele
            let variant_idx = i / 2;
            let variant = variants[variant_idx];
            let allele = variant.hap1.as_bytes();

            hap1_seq_parts.push(allele.to_vec());

            // Handle quality scores for length changes
            let orig_len = seq_part.len();
            let allele_len = allele.len();

            if allele_len == orig_len {
                // Same length - use original qualities
                hap1_qual_parts.push(split_qual[i].to_vec());
            } else if allele_len < orig_len {
                // Deletion - truncate qualities
                hap1_qual_parts.push(split_qual[i][..allele_len].to_vec());
            } else {
                // Insertion - fill extra qualities
                let extra_len = allele_len - orig_len;
                let left_qual = if i > 0 { split_qual[i - 1] } else { &[] };
                let right_qual = if i < split_qual.len() - 1 { split_qual[i + 1] } else { &[] };

                let extra_quals = fill_insertion_quals(extra_len, left_qual, right_qual, 30);
                let mut combined = split_qual[i].to_vec();
                combined.extend(extra_quals);
                hap1_qual_parts.push(combined);
            }
        }
    }

    // Build haplotype 2 with quality-aware allele swapping
    let mut hap2_seq_parts: Vec<Vec<u8>> = Vec::new();
    let mut hap2_qual_parts: Vec<Vec<u8>> = Vec::new();

    for (i, seq_part) in split_seq.iter().enumerate() {
        if i % 2 == 0 {
            // Non-variant segment - same for both haplotypes
            hap2_seq_parts.push(seq_part.to_vec());
            hap2_qual_parts.push(split_qual[i].to_vec());
        } else {
            // Variant segment - swap allele
            let variant_idx = i / 2;
            let variant = variants[variant_idx];
            let allele = variant.hap2.as_bytes();

            hap2_seq_parts.push(allele.to_vec());

            // Handle quality scores for length changes
            let orig_len = seq_part.len();
            let allele_len = allele.len();

            if allele_len == orig_len {
                // Same length - use original qualities
                hap2_qual_parts.push(split_qual[i].to_vec());
            } else if allele_len < orig_len {
                // Deletion - truncate qualities
                hap2_qual_parts.push(split_qual[i][..allele_len].to_vec());
            } else {
                // Insertion - fill extra qualities
                let extra_len = allele_len - orig_len;
                let left_qual = if i > 0 { split_qual[i - 1] } else { &[] };
                let right_qual = if i < split_qual.len() - 1 { split_qual[i + 1] } else { &[] };

                let extra_quals = fill_insertion_quals(extra_len, left_qual, right_qual, 30);
                let mut combined = split_qual[i].to_vec();
                combined.extend(extra_quals);
                hap2_qual_parts.push(combined);
            }
        }
    }

    // Join segments to create final sequences and qualities
    let hap1_seq: Vec<u8> = hap1_seq_parts.into_iter().flatten().collect();
    let hap1_qual: Vec<u8> = hap1_qual_parts.into_iter().flatten().collect();
    let hap2_seq: Vec<u8> = hap2_seq_parts.into_iter().flatten().collect();
    let hap2_qual: Vec<u8> = hap2_qual_parts.into_iter().flatten().collect();

    Ok(Some(vec![(hap1_seq, hap1_qual), (hap2_seq, hap2_qual)]))
}

/// Write haplotype reads to FASTQ files (paired-end)
///
/// # Arguments
/// * `haplotypes` - Generated haplotype reads
/// * `r1_path` - Output path for read 1 FASTQ
/// * `r2_path` - Output path for read 2 FASTQ
///
/// # Returns
/// (read1_count, read2_count)
pub fn write_fastq_pair<P: AsRef<Path>>(
    haplotypes: &[HaplotypeRead],
    r1_path: P,
    r2_path: P,
) -> Result<(usize, usize)> {
    use std::io::Write as IoWrite;

    let mut r1_file = std::io::BufWriter::new(
        File::create(r1_path.as_ref()).context("Failed to create R1 FASTQ")?
    );
    let mut r2_file = std::io::BufWriter::new(
        File::create(r2_path.as_ref()).context("Failed to create R2 FASTQ")?
    );

    let mut r1_count = 0;
    let mut r2_count = 0;

    // Write each haplotype to the appropriate file
    for hap in haplotypes {
        // Determine if this is R1 or R2 by checking the name suffix
        let is_r1 = hap.name.ends_with(b"/1");

        // Convert quality scores to ASCII (Phred+33)
        let qual_string: Vec<u8> = hap.quals.iter().map(|&q| q + 33).collect();

        // Write FASTQ format: @name\nseq\n+\nquals\n
        let fastq_entry = format!(
            "@{}\n{}\n+\n{}\n",
            String::from_utf8_lossy(&hap.name),
            String::from_utf8_lossy(&hap.sequence),
            String::from_utf8_lossy(&qual_string)
        );

        if is_r1 {
            r1_file.write_all(fastq_entry.as_bytes())
                .context("Failed to write R1 FASTQ entry")?;
            r1_count += 1;
        } else {
            r2_file.write_all(fastq_entry.as_bytes())
                .context("Failed to write R2 FASTQ entry")?;
            r2_count += 1;
        }
    }

    // Flush buffers
    r1_file.flush().context("Failed to flush R1 file")?;
    r2_file.flush().context("Failed to flush R2 file")?;

    Ok((r1_count, r2_count))
}

/// Process all chromosomes in parallel using pre-grouped variants
///
/// Uses rayon for parallel processing of independent chromosomes.
/// This is the optimized version that takes pre-parsed, chromosome-grouped variants.
///
/// # Arguments
/// * `bam_path` - Path to BAM file
/// * `variants_by_chrom` - Variants pre-grouped by chromosome (from parse_intersect_bed_by_chrom)
/// * `config` - Remapping configuration
///
/// # Returns
/// Vector of all haplotype reads from all chromosomes + aggregated stats
///
/// # Performance
/// - Parse once instead of 22x: ~22x faster parsing
/// - Parallel chromosome processing: Additional 4-8x speedup with 8 cores
/// - Total expected speedup: ~100x for large RNA-seq datasets
pub fn process_all_chromosomes_parallel(
    bam_path: &str,
    variants_by_chrom: &FxHashMap<String, FxHashMap<Vec<u8>, Vec<VariantSpan>>>,
    config: &RemapConfig,
) -> Result<(Vec<HaplotypeRead>, RemapStats)> {
    use rayon::prelude::*;

    // Get list of chromosomes to process
    let chromosomes: Vec<&String> = variants_by_chrom.keys().collect();

    if chromosomes.is_empty() {
        return Ok((Vec::new(), RemapStats::default()));
    }

    // Process chromosomes in parallel
    // Each thread gets its own BAM reader (IndexedReader is not Send)
    let results: Vec<Result<(Vec<HaplotypeRead>, RemapStats)>> = chromosomes
        .par_iter()
        .map(|chrom| {
            // Get variants for this chromosome
            let chrom_variants = variants_by_chrom.get(*chrom).unwrap();

            // Process this chromosome (opens its own BAM reader)
            swap_alleles_for_chrom(bam_path, chrom_variants, chrom, config)
        })
        .collect();

    // Combine results from all chromosomes
    let mut all_haplotypes: Vec<HaplotypeRead> = Vec::new();
    let mut combined_stats = RemapStats::default();

    for result in results {
        let (haplotypes, stats) = result?;
        all_haplotypes.extend(haplotypes);
        combined_stats.pairs_processed += stats.pairs_processed;
        combined_stats.pairs_with_variants += stats.pairs_with_variants;
        combined_stats.haplotypes_generated += stats.haplotypes_generated;
        combined_stats.reads_discarded += stats.reads_discarded;
    }

    Ok((all_haplotypes, combined_stats))
}

/// Process all chromosomes in parallel with streaming FASTQ writes
///
/// Uses crossbeam channels for producer-consumer pattern:
/// - Producer threads: Process chromosomes in parallel (Rayon)
/// - Consumer thread: Write FASTQ entries as they arrive
///
/// This eliminates memory accumulation and enables overlapped I/O.
///
/// # Arguments
/// * `bam_path` - Path to BAM file
/// * `variants_by_chrom` - Variants pre-grouped by chromosome
/// * `config` - Remapping configuration
/// * `r1_path` - Output path for R1 FASTQ
/// * `r2_path` - Output path for R2 FASTQ
/// * `num_threads` - Number of threads for parallel processing (0 = auto)
///
/// # Performance
/// - Streaming writes: Memory-efficient, no accumulation
/// - Overlapped I/O: Writing happens while processing continues
/// - Thread pool control: User-specified thread count
pub fn process_and_write_parallel<P: AsRef<std::path::Path>>(
    bam_path: &str,
    variants_by_chrom: &FxHashMap<String, FxHashMap<Vec<u8>, Vec<VariantSpan>>>,
    config: &RemapConfig,
    r1_path: P,
    r2_path: P,
    num_threads: usize,
) -> Result<RemapStats> {
    use crossbeam_channel::{bounded, Sender};
    use rayon::prelude::*;
    use std::io::Write as IoWrite;
    use std::thread;

    // Configure thread pool if specified
    if num_threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build_global()
            .ok(); // Ignore error if already initialized
    }

    let chromosomes: Vec<&String> = variants_by_chrom.keys().collect();
    if chromosomes.is_empty() {
        // Create empty output files
        std::fs::File::create(r1_path.as_ref())?;
        std::fs::File::create(r2_path.as_ref())?;
        return Ok(RemapStats::default());
    }

    // Bounded channel to prevent unbounded memory growth
    // Buffer ~1000 haplotypes at a time
    let (tx, rx): (Sender<HaplotypeRead>, _) = bounded(1000);

    // Clone paths for writer thread
    let r1_path_str = r1_path.as_ref().to_path_buf();
    let r2_path_str = r2_path.as_ref().to_path_buf();

    // Spawn writer thread (consumer)
    let writer_handle = thread::spawn(move || -> Result<(usize, usize)> {
        let mut r1_file = std::io::BufWriter::new(
            std::fs::File::create(&r1_path_str).context("Failed to create R1 FASTQ")?
        );
        let mut r2_file = std::io::BufWriter::new(
            std::fs::File::create(&r2_path_str).context("Failed to create R2 FASTQ")?
        );

        let mut r1_count = 0;
        let mut r2_count = 0;

        // Receive and write haplotypes as they arrive
        for hap in rx {
            let is_r1 = hap.name.ends_with(b"/1");
            let qual_string: Vec<u8> = hap.quals.iter().map(|&q| q + 33).collect();

            let fastq_entry = format!(
                "@{}\n{}\n+\n{}\n",
                String::from_utf8_lossy(&hap.name),
                String::from_utf8_lossy(&hap.sequence),
                String::from_utf8_lossy(&qual_string)
            );

            if is_r1 {
                r1_file.write_all(fastq_entry.as_bytes())
                    .context("Failed to write R1 FASTQ entry")?;
                r1_count += 1;
            } else {
                r2_file.write_all(fastq_entry.as_bytes())
                    .context("Failed to write R2 FASTQ entry")?;
                r2_count += 1;
            }
        }

        r1_file.flush().context("Failed to flush R1 file")?;
        r2_file.flush().context("Failed to flush R2 file")?;

        Ok((r1_count, r2_count))
    });

    // Process chromosomes in parallel (producers)
    let results: Vec<Result<RemapStats>> = chromosomes
        .par_iter()
        .map(|chrom| {
            let chrom_variants = variants_by_chrom.get(*chrom).unwrap();
            let tx = tx.clone();

            // Process chromosome
            let (haplotypes, stats) = swap_alleles_for_chrom(bam_path, chrom_variants, chrom, config)?;

            // Stream haplotypes to writer
            for hap in haplotypes {
                // If channel is closed, writer failed - abort
                if tx.send(hap).is_err() {
                    return Err(anyhow::anyhow!("Writer thread failed"));
                }
            }

            Ok(stats)
        })
        .collect();

    // Drop the sender to signal completion to writer
    drop(tx);

    // Wait for writer to finish
    let (_r1_count, _r2_count) = writer_handle
        .join()
        .map_err(|_| anyhow::anyhow!("Writer thread panicked"))??;

    // Aggregate stats
    let mut combined_stats = RemapStats::default();
    for result in results {
        let stats = result?;
        combined_stats.pairs_processed += stats.pairs_processed;
        combined_stats.pairs_with_variants += stats.pairs_with_variants;
        combined_stats.haplotypes_generated += stats.haplotypes_generated;
        combined_stats.reads_discarded += stats.reads_discarded;
    }

    Ok(combined_stats)
}

/// Process all chromosomes sequentially (for comparison/fallback)
///
/// Same as parallel version but processes chromosomes one at a time.
pub fn process_all_chromosomes_sequential(
    bam_path: &str,
    variants_by_chrom: &FxHashMap<String, FxHashMap<Vec<u8>, Vec<VariantSpan>>>,
    config: &RemapConfig,
) -> Result<(Vec<HaplotypeRead>, RemapStats)> {
    let mut all_haplotypes: Vec<HaplotypeRead> = Vec::new();
    let mut combined_stats = RemapStats::default();

    for (chrom, chrom_variants) in variants_by_chrom.iter() {
        let (haplotypes, stats) = swap_alleles_for_chrom(bam_path, chrom_variants, chrom, config)?;
        all_haplotypes.extend(haplotypes);
        combined_stats.pairs_processed += stats.pairs_processed;
        combined_stats.pairs_with_variants += stats.pairs_with_variants;
        combined_stats.haplotypes_generated += stats.haplotypes_generated;
        combined_stats.reads_discarded += stats.reads_discarded;
    }

    Ok((all_haplotypes, combined_stats))
}

// ============================================================================
// Helper Functions
// ============================================================================

/// Build reference-to-read position mappings for indel support
///
/// Creates two mappings to handle deletions properly:
/// - ref2q_left: Maps reference position to nearest left query position
/// - ref2q_right: Maps reference position to nearest right query position
///
/// For positions in deletions, left and right maps will differ (flanking positions).
/// For matches, both maps point to the same query position.
///
/// This uses the shared `cigar_utils::build_ref2query_maps()` which leverages
/// rust-htslib's `aligned_pairs_full()` API (equivalent to pysam's
/// `get_aligned_pairs(matches_only=False)`).
///
/// # Returns
/// (ref2q_left, ref2q_right) FxHashMaps mapping reference positions to query positions
fn build_ref2read_maps(read: &bam::Record) -> (FxHashMap<u32, usize>, FxHashMap<u32, usize>) {
    // Use the shared cigar_utils implementation which uses aligned_pairs_full()
    let (left_i64, right_i64) = cigar_utils::build_ref2query_maps(read);

    // Convert from i64 keys to u32 keys (for backwards compatibility with existing code)
    let ref2q_left: FxHashMap<u32, usize> = left_i64
        .into_iter()
        .filter_map(|(k, v)| {
            if k >= 0 && k <= u32::MAX as i64 {
                Some((k as u32, v))
            } else {
                None
            }
        })
        .collect();

    let ref2q_right: FxHashMap<u32, usize> = right_i64
        .into_iter()
        .filter_map(|(k, v)| {
            if k >= 0 && k <= u32::MAX as i64 {
                Some((k as u32, v))
            } else {
                None
            }
        })
        .collect();

    (ref2q_left, ref2q_right)
}

/// Fill quality scores for inserted bases
///
/// When an insertion makes a haplotype longer than the original read,
/// we need to generate quality scores for the extra bases.
///
/// Strategy: Average the flanking quality scores, or use default Q30.
///
/// Mirrors Python's `_fill_insertion_quals()` in remap_utils.py (lines 204-223)
fn fill_insertion_quals(
    insert_len: usize,
    left_qual: &[u8],
    right_qual: &[u8],
    insert_qual: u8,
) -> Vec<u8> {
    if left_qual.is_empty() && right_qual.is_empty() {
        // No flanking data - use default
        return vec![insert_qual; insert_len];
    }

    // Average flanking qualities
    let mut flank_quals = Vec::new();
    flank_quals.extend_from_slice(left_qual);
    flank_quals.extend_from_slice(right_qual);

    let sum: u32 = flank_quals.iter().map(|&q| q as u32).sum();
    let mean_qual = (sum / flank_quals.len() as u32) as u8;

    vec![mean_qual; insert_len]
}

/// Find read position for a given reference position (optimized)
///
/// Walks CIGAR string to find read position corresponding to genomic position.
/// This is O(k) where k = number of CIGAR operations, instead of O(n) where n = read length.
///
/// Much faster than building a full HashMap when you only need a few lookups.
///
/// # Returns
/// Some(read_pos) if position is mapped, None if in deletion/unmapped region
fn find_read_position(read: &bam::Record, target_ref_pos: u32) -> Option<usize> {
    let cigar = read.cigar();
    let mut read_pos: usize = 0;
    let mut ref_pos = read.pos() as u32;

    for op in cigar.iter() {
        use rust_htslib::bam::record::Cigar;

        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                // Check if target position is in this match block
                if target_ref_pos >= ref_pos && target_ref_pos < ref_pos + len {
                    let offset = (target_ref_pos - ref_pos) as usize;
                    return Some(read_pos + offset);
                }
                read_pos += *len as usize;
                ref_pos += len;
            }
            Cigar::Ins(len) => {
                // Insertion: only read advances
                read_pos += *len as usize;
            }
            Cigar::Del(len) | Cigar::RefSkip(len) => {
                // Deletion/skip: only reference advances
                // If target is in deletion, return None
                if target_ref_pos >= ref_pos && target_ref_pos < ref_pos + len {
                    return None;
                }
                ref_pos += len;
            }
            Cigar::SoftClip(len) => {
                // Soft clip: only read advances
                read_pos += *len as usize;
            }
            Cigar::HardClip(_) | Cigar::Pad(_) => {
                // Hard clip/pad: no advancement
            }
        }
    }

    None // Position not found in alignment
}

/// Generate WASP read name
///
/// Format: {original_name}_WASP_{pos1}_{pos2}_{seq_num}_{total_seqs}
/// Matches Python's: f"{og_name}_WASP_{r1_align_pos}_{r2_align_pos}_{write_num}_{write_total}"
///
/// # Arguments
/// * `original_name` - Original read name
/// * `pos1` - Read 1 alignment position
/// * `pos2` - Read 2 alignment position
/// * `seq_num` - Index of this sequence (1-based)
/// * `total_seqs` - Total number of sequences generated for this pair
fn generate_wasp_name(
    original_name: &[u8],
    pos1: u32,
    pos2: u32,
    seq_num: usize,
    total_seqs: usize,
) -> Vec<u8> {
    let name_str = std::str::from_utf8(original_name).unwrap_or("unknown");
    format!("{}_WASP_{}_{}_{}_{}", name_str, pos1, pos2, seq_num, total_seqs).into_bytes()
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_parse_intersect_bed() {
        // Create test BED file
        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "chr10\t87377\t87427\tSRR891276.10516353/2\t60\t+\tchr10\t87400\t87401\tC\tT\tC|T").unwrap();
        writeln!(temp_file, "chr10\t87392\t87440\tSRR891276.5620594/2\t60\t+\tchr10\t87400\t87401\tC\tT\tC|T").unwrap();
        writeln!(temp_file, "chr10\t87395\t87442\tSRR891276.5620594/1\t60\t-\tchr10\t87400\t87401\tC\tT\tC|T").unwrap();
        // Duplicate that should be removed
        writeln!(temp_file, "chr10\t87392\t87440\tSRR891276.5620594/2\t60\t+\tchr10\t87401\t87402\tA\tG\tA|G").unwrap();
        temp_file.flush().unwrap();

        // Parse
        let result = parse_intersect_bed(temp_file.path()).unwrap();

        // Verify
        assert_eq!(result.len(), 2, "Should have 2 unique reads");

        // Check first read
        let read1_name = b"SRR891276.10516353".to_vec();
        let read1_spans = result.get(&read1_name).unwrap();
        assert_eq!(read1_spans.len(), 1);
        assert_eq!(read1_spans[0].chrom, "chr10");
        assert_eq!(read1_spans[0].start, 87377);
        assert_eq!(read1_spans[0].stop, 87427);
        assert_eq!(read1_spans[0].vcf_start, 87400);
        assert_eq!(read1_spans[0].vcf_stop, 87401);
        assert_eq!(read1_spans[0].mate, 2);
        assert_eq!(read1_spans[0].hap1, "C");
        assert_eq!(read1_spans[0].hap2, "T");

        // Check second read (should have deduplication)
        let read2_name = b"SRR891276.5620594".to_vec();
        let read2_spans = result.get(&read2_name).unwrap();
        assert_eq!(read2_spans.len(), 2, "Should have 2 spans after dedup (different mates)");

        // Verify mate 1
        let mate1 = read2_spans.iter().find(|s| s.mate == 1).unwrap();
        assert_eq!(mate1.start, 87395);
        assert_eq!(mate1.stop, 87442);
        assert_eq!(mate1.vcf_start, 87400);
        assert_eq!(mate1.vcf_stop, 87401);

        // Verify mate 2 (should only have first occurrence, duplicate removed)
        let mate2 = read2_spans.iter().find(|s| s.mate == 2).unwrap();
        assert_eq!(mate2.start, 87392);
        assert_eq!(mate2.stop, 87440);
        assert_eq!(mate2.vcf_start, 87400);
        assert_eq!(mate2.vcf_stop, 87401);
    }

    #[test]
    #[ignore]
    fn test_generate_haplotype_seqs() {
        // TODO: Create mock BAM record
        // TODO: Create test variants
        // TODO: Generate haplotypes
        // TODO: Verify sequences are correct
    }

    #[test]
    #[ignore]
    fn test_build_alignment_map() {
        // TODO: Create read with known alignment
        // TODO: Build map
        // TODO: Verify positions are correct
    }

    #[test]
    #[ignore]
    fn test_generate_wasp_name() {
        // TODO: Generate name with test inputs
        // TODO: Verify format matches Python implementation
    }
}
