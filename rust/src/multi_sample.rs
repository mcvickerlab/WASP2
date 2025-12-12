//! Multi-sample support for BAM remapping
//!
//! Extends the single-sample Rust implementation to handle multiple samples.
//! This enables the full Rust acceleration path for multi-sample WASP2 runs.
//!
//! # Key Differences from Single-Sample
//!
//! Single-sample: Always generates 2 haplotypes (hap1, hap2)
//! Multi-sample: Generates all unique haplotype combinations across samples
//!
//! For example, with 2 samples at 1 variant:
//! - Sample1: A|G
//! - Sample2: A|T
//! - Unique combinations: [A], [G], [T] = 3 sequences (not 4, since A appears twice)
//!
//! # Data Flow
//! 1. VCF → BED with multi-sample genotypes
//! 2. BAM-BED intersection outputs all sample GTs per read-variant overlap
//! 3. parse_intersect_bed_multi() parses multi-sample genotypes
//! 4. generate_unique_combinations() finds unique allele sets
//! 5. Each unique combination generates one output sequence
//!
//! # INDEL Support (v1.2+)
//!
//! Uses CIGAR-aware position mapping via `cigar_utils::build_ref2query_maps()`.
//! This properly handles reads with insertions/deletions in their alignment.

use anyhow::{Context, Result};
use rustc_hash::FxHashMap;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::cigar_utils;

// ============================================================================
// Data Structures
// ============================================================================

/// Variant span for multi-sample processing
///
/// Unlike single-sample VariantSpan which stores just (hap1, hap2),
/// this stores alleles for ALL samples at this variant position.
#[derive(Debug, Clone)]
pub struct VariantSpanMulti {
    /// Chromosome name
    pub chrom: String,
    /// Read start position (from BAM)
    pub start: u32,
    /// Read stop position (from BAM)
    pub stop: u32,
    /// Variant start position (from VCF/BED)
    pub vcf_start: u32,
    /// Variant stop position (from VCF/BED)
    pub vcf_stop: u32,
    /// Mate number (1 or 2)
    pub mate: u8,
    /// Reference allele
    pub ref_allele: String,
    /// Alternate allele
    pub alt_allele: String,
    /// Per-sample alleles: [(hap1_s1, hap2_s1), (hap1_s2, hap2_s2), ...]
    pub sample_alleles: Vec<(String, String)>,
}

/// Multi-sample variant store for intersection output
pub type MultiSampleVariants = FxHashMap<Vec<u8>, Vec<VariantSpanMulti>>;

// ============================================================================
// BED Parsing
// ============================================================================

/// Parse multi-sample intersection BED file
///
/// Expected format (12 + N columns for N samples):
/// ```text
/// chrom  start  end  read/mate  mapq  strand  vcf_chrom  vcf_start  vcf_end  ref  alt  GT_S1  GT_S2  ...
/// chr10  100    200  readA/1    60    +       chr10      150        151      A    G    A|G    A|A    ...
/// ```
///
/// # Arguments
/// * `intersect_bed` - Path to bedtools intersect output
/// * `num_samples` - Number of samples (determines column count)
///
/// # Returns
/// HashMap mapping read names to their variant spans with all sample genotypes
pub fn parse_intersect_bed_multi<P: AsRef<Path>>(
    intersect_bed: P,
    num_samples: usize,
) -> Result<MultiSampleVariants> {
    let file =
        File::open(intersect_bed.as_ref()).context("Failed to open intersection BED file")?;
    let reader = BufReader::with_capacity(1024 * 1024, file);

    let mut variants: MultiSampleVariants = FxHashMap::default();
    let mut seen: HashSet<(Vec<u8>, String, u32, u32, u8)> = HashSet::default();

    let mut line_count = 0;
    let mut skipped_count = 0;

    for line in reader.lines() {
        let line = line?;
        line_count += 1;

        let fields: Vec<&str> = line.split('\t').collect();

        // Expected columns: 11 base columns + num_samples genotype columns
        let expected_cols = 11 + num_samples;
        if fields.len() < expected_cols {
            skipped_count += 1;
            continue;
        }

        // Parse basic fields
        let chrom = fields[0].to_string();
        let start = fields[1]
            .parse::<u32>()
            .context("Failed to parse read start")?;
        let stop = fields[2]
            .parse::<u32>()
            .context("Failed to parse read stop")?;
        let read_with_mate = fields[3];

        // Parse VCF fields
        let vcf_start = fields[7]
            .parse::<u32>()
            .context("Failed to parse vcf_start")?;
        let vcf_stop = fields[8]
            .parse::<u32>()
            .context("Failed to parse vcf_stop")?;
        let ref_allele = fields[9].to_string();
        let alt_allele = fields[10].to_string();

        // Parse read name and mate
        let parts: Vec<&str> = read_with_mate.split('/').collect();
        if parts.len() != 2 {
            skipped_count += 1;
            continue;
        }
        let read_name = parts[0].as_bytes().to_vec();
        let mate = parts[1]
            .parse::<u8>()
            .context("Failed to parse mate number")?;

        // Deduplication key (same as Python's unique(["chrom", "read", "mate", "start", "stop"]))
        let key = (read_name.clone(), chrom.clone(), start, stop, mate);
        if seen.contains(&key) {
            continue;
        }
        seen.insert(key);

        // Parse per-sample genotypes (columns 11, 12, 13, ...)
        let mut sample_alleles = Vec::with_capacity(num_samples);
        for i in 0..num_samples {
            let gt_col = 11 + i;
            let gt = fields[gt_col];

            // Try phased first (|), then unphased (/)
            let alleles: Vec<&str> = if gt.contains('|') {
                gt.split('|').collect()
            } else {
                gt.split('/').collect()
            };

            if alleles.len() == 2 {
                sample_alleles.push((alleles[0].to_string(), alleles[1].to_string()));
            } else {
                // Missing or malformed - use reference
                sample_alleles.push((".".to_string(), ".".to_string()));
            }
        }

        let span = VariantSpanMulti {
            chrom,
            start,
            stop,
            vcf_start,
            vcf_stop,
            mate,
            ref_allele,
            alt_allele,
            sample_alleles,
        };

        variants
            .entry(read_name)
            .or_insert_with(Vec::new)
            .push(span);
    }

    eprintln!(
        "  Parsed {} lines, {} unique read-variant pairs, {} skipped",
        line_count,
        variants.len(),
        skipped_count
    );

    Ok(variants)
}

// ============================================================================
// Unique Haplotype Column Generation (Matches Python Logic)
// ============================================================================

/// Generate unique haplotype columns across samples
///
/// This matches the Python logic in swap_chrom_alleles_multi:
/// 1. Each sample has 2 haplotype columns (hap1, hap2)
/// 2. Concatenate alleles in each column across all variants
/// 3. Find unique concatenated strings (columns with identical patterns)
/// 4. Return unique column indices to use for sequence generation
///
/// # Example
/// 2 samples, 2 variants:
/// - Sample1: pos100=A|G, pos200=C|T  → col0="AC", col1="GT"
/// - Sample2: pos100=A|A, pos200=C|C  → col2="AC", col3="CC"
/// Unique columns: ["AC", "GT", "CC"] → indices [0, 1, 3]
///
/// # Arguments
/// * `variants` - Slice of variant spans for a single read (must have same sample count)
///
/// # Returns
/// Vector of unique (column_index, alleles_vec) pairs
pub fn generate_unique_haplotype_columns(
    variants: &[&VariantSpanMulti],
) -> Vec<(usize, Vec<String>)> {
    if variants.is_empty() {
        return vec![];
    }

    // Determine number of haplotype columns (2 per sample)
    let num_samples = variants[0].sample_alleles.len();
    let num_columns = num_samples * 2;

    // Build concatenated string for each column across all variants
    let mut column_signatures: Vec<(usize, String, Vec<String>)> = Vec::with_capacity(num_columns);

    for col_idx in 0..num_columns {
        let sample_idx = col_idx / 2;
        let is_hap2 = col_idx % 2 == 1;

        let mut signature = String::new();
        let mut alleles = Vec::with_capacity(variants.len());

        for v in variants {
            if sample_idx < v.sample_alleles.len() {
                let (hap1, hap2) = &v.sample_alleles[sample_idx];
                let allele = if is_hap2 { hap2 } else { hap1 };
                signature.push_str(allele);
                alleles.push(allele.clone());
            }
        }

        column_signatures.push((col_idx, signature, alleles));
    }

    // Find unique signatures
    let mut seen_signatures: HashSet<String> = HashSet::new();
    let mut unique_columns: Vec<(usize, Vec<String>)> = Vec::new();

    for (col_idx, signature, alleles) in column_signatures {
        // Skip columns with missing data
        if signature.contains('.') {
            continue;
        }

        if !seen_signatures.contains(&signature) {
            seen_signatures.insert(signature);
            unique_columns.push((col_idx, alleles));
        }
    }

    unique_columns
}

/// Generate all unique allele combinations across variants
///
/// Wrapper that extracts just the allele vectors from unique columns.
///
/// # Arguments
/// * `variants` - Slice of variant spans for a single read
///
/// # Returns
/// Vector of allele combinations, where each inner vector has one allele per variant
pub fn generate_unique_combinations(variants: &[&VariantSpanMulti]) -> Vec<Vec<String>> {
    let unique_cols = generate_unique_haplotype_columns(variants);
    unique_cols
        .into_iter()
        .map(|(_, alleles)| alleles)
        .collect()
}

// ============================================================================
// Sequence Generation (CIGAR-Aware)
// ============================================================================

/// Apply allele substitutions using CIGAR-aware position mapping
///
/// This is the CORRECT implementation that handles reads with insertions/deletions
/// in their CIGAR string. The naive `offset = ref_pos - read_start` approach fails
/// when the read's alignment includes indels.
///
/// # Arguments
/// * `seq` - Original read sequence
/// * `qual` - Original quality scores
/// * `variants` - Variant spans overlapping this read
/// * `alleles` - Alleles to substitute (one per variant)
/// * `ref2query_left` - Left position mapping from cigar_utils
/// * `ref2query_right` - Right position mapping from cigar_utils
///
/// # Returns
/// (new_sequence, new_quality) with substitutions applied
pub fn apply_allele_substitutions_cigar_aware(
    seq: &[u8],
    qual: &[u8],
    variants: &[&VariantSpanMulti],
    alleles: &[String],
    ref2query_left: &FxHashMap<i64, usize>,
    ref2query_right: &FxHashMap<i64, usize>,
) -> Result<(Vec<u8>, Vec<u8>)> {
    if variants.is_empty() {
        return Ok((seq.to_vec(), qual.to_vec()));
    }

    // Convert variants to position tuples for segmentation
    let mut variant_positions: Vec<(usize, usize)> = Vec::with_capacity(variants.len());

    for variant in variants.iter() {
        let ref_start = variant.vcf_start as i64;
        let ref_end = variant.vcf_stop as i64;

        // Get query positions using CIGAR-aware mapping
        let query_start = ref2query_left.get(&ref_start).copied().ok_or_else(|| {
            anyhow::anyhow!(
                "Variant at ref {} not in left map (read may not cover variant)",
                ref_start
            )
        })?;

        // For end: use right mapping for ref_end - 1, then add 1
        let query_end = ref2query_right
            .get(&(ref_end - 1))
            .map(|&p| p + 1)
            .ok_or_else(|| anyhow::anyhow!("Variant at ref {} not in right map", ref_end - 1))?;

        variant_positions.push((query_start, query_end.min(seq.len())));
    }

    // Segment the sequence at variant positions
    let (seq_segments, qual_segments) =
        cigar_utils::segment_sequence(seq, qual, &variant_positions);

    // Build new sequence with allele substitutions
    let mut new_seq = Vec::with_capacity(seq.len());
    let mut new_qual = Vec::with_capacity(qual.len());

    for (i, (seq_seg, qual_seg)) in seq_segments.iter().zip(qual_segments.iter()).enumerate() {
        if i % 2 == 0 {
            // Non-variant segment: copy as-is
            new_seq.extend_from_slice(seq_seg);
            new_qual.extend_from_slice(qual_seg);
        } else {
            // Variant segment: substitute with allele
            let variant_idx = i / 2;
            if variant_idx < alleles.len() {
                let allele = &alleles[variant_idx];
                let allele_bytes = allele.as_bytes();

                new_seq.extend_from_slice(allele_bytes);

                // Handle quality scores for length changes
                let orig_len = seq_seg.len();
                let allele_len = allele_bytes.len();

                if allele_len == orig_len {
                    // Same length: use original qualities
                    new_qual.extend_from_slice(qual_seg);
                } else if allele_len < orig_len {
                    // Deletion: truncate qualities
                    new_qual.extend_from_slice(&qual_seg[..allele_len.min(qual_seg.len())]);
                } else {
                    // Insertion: use original + fill extra with Q30
                    new_qual.extend_from_slice(qual_seg);
                    let extra_needed = allele_len.saturating_sub(orig_len);
                    new_qual.extend(std::iter::repeat(30u8).take(extra_needed));
                }
            }
        }
    }

    Ok((new_seq, new_qual))
}

/// Legacy function for backwards compatibility (DEPRECATED)
///
/// WARNING: This function uses naive offset calculation that fails for reads
/// with insertions/deletions in their CIGAR string. Use
/// `apply_allele_substitutions_cigar_aware` or `generate_multi_sample_sequences_from_record`
/// instead.
#[deprecated(
    since = "1.2.0",
    note = "Use apply_allele_substitutions_cigar_aware instead"
)]
#[allow(dead_code)]
pub fn apply_allele_substitutions(
    seq: &[u8],
    qual: &[u8],
    variants: &[&VariantSpanMulti],
    alleles: &[String],
    read_start: u32,
) -> Result<(Vec<u8>, Vec<u8>)> {
    let mut new_seq = seq.to_vec();
    let mut new_qual = qual.to_vec();

    // Apply each substitution (naive offset - ONLY works for simple CIGAR like 150M)
    for (variant, allele) in variants.iter().zip(alleles.iter()) {
        let var_pos = variant.vcf_start;

        if var_pos >= read_start {
            let offset = (var_pos - read_start) as usize;

            if offset < new_seq.len() {
                let ref_len = variant.ref_allele.len();
                let alt_len = allele.len();

                if ref_len == 1 && alt_len == 1 {
                    new_seq[offset] = allele.as_bytes()[0];
                } else if ref_len > alt_len {
                    if offset + ref_len <= new_seq.len() {
                        for (i, b) in allele.bytes().enumerate() {
                            if offset + i < new_seq.len() {
                                new_seq[offset + i] = b;
                            }
                        }
                        let remove_start = offset + alt_len;
                        let remove_end = offset + ref_len;
                        if remove_end <= new_seq.len() {
                            new_seq.drain(remove_start..remove_end);
                            new_qual.drain(remove_start..remove_end);
                        }
                    }
                } else if alt_len > ref_len {
                    if offset + ref_len <= new_seq.len() {
                        for (i, b) in allele.bytes().take(ref_len).enumerate() {
                            new_seq[offset + i] = b;
                        }
                        let insert_pos = offset + ref_len;
                        let extra_bases: Vec<u8> = allele.bytes().skip(ref_len).collect();
                        let extra_qual: Vec<u8> = vec![30; extra_bases.len()];

                        for (i, (b, q)) in extra_bases.iter().zip(extra_qual.iter()).enumerate() {
                            new_seq.insert(insert_pos + i, *b);
                            new_qual.insert(insert_pos + i, *q);
                        }
                    }
                }
            }
        }
    }

    Ok((new_seq, new_qual))
}

/// Generate haplotype sequences from a BAM record with CIGAR awareness
///
/// This is the CORRECT entry point for multi-sample sequence generation.
/// It uses the BAM record's CIGAR string to properly map variant positions.
///
/// # Arguments
/// * `read` - BAM record with CIGAR information
/// * `variants` - Variant spans overlapping this read
///
/// # Returns
/// Vector of (sequence, quality) pairs, one per unique haplotype
pub fn generate_multi_sample_sequences_from_record(
    read: &rust_htslib::bam::Record,
    variants: &[&VariantSpanMulti],
) -> Result<Vec<(Vec<u8>, Vec<u8>)>> {
    if variants.is_empty() {
        let seq = read.seq().as_bytes();
        let qual = read.qual().to_vec();
        return Ok(vec![(seq, qual)]);
    }

    // Build CIGAR-aware position maps
    let (ref2query_left, ref2query_right) = cigar_utils::build_ref2query_maps(read);

    let seq = read.seq().as_bytes();
    let qual = read.qual().to_vec();

    // Generate unique allele combinations
    let combinations = generate_unique_combinations(variants);

    let mut results = Vec::with_capacity(combinations.len());

    for alleles in combinations {
        match apply_allele_substitutions_cigar_aware(
            &seq,
            &qual,
            variants,
            &alleles,
            &ref2query_left,
            &ref2query_right,
        ) {
            Ok((new_seq, new_qual)) => results.push((new_seq, new_qual)),
            Err(e) => {
                // Log error but continue - variant may not overlap read properly
                eprintln!("Warning: failed to apply substitution: {}", e);
                continue;
            }
        }
    }

    // If all combinations failed, return original
    if results.is_empty() {
        results.push((seq, qual));
    }

    Ok(results)
}

/// Legacy function - DEPRECATED
///
/// Use `generate_multi_sample_sequences_from_record` instead.
#[deprecated(
    since = "1.2.0",
    note = "Use generate_multi_sample_sequences_from_record instead"
)]
#[allow(dead_code)]
pub fn generate_multi_sample_sequences(
    seq: &[u8],
    qual: &[u8],
    variants: &[&VariantSpanMulti],
    read_start: u32,
) -> Result<Vec<(Vec<u8>, Vec<u8>)>> {
    let combinations = generate_unique_combinations(variants);

    let mut results = Vec::with_capacity(combinations.len());

    #[allow(deprecated)]
    for alleles in combinations {
        let (new_seq, new_qual) =
            apply_allele_substitutions(seq, qual, variants, &alleles, read_start)?;
        results.push((new_seq, new_qual));
    }

    Ok(results)
}

// ============================================================================
// Full Multi-Sample Remapping Pipeline
// ============================================================================

use rust_htslib::{bam, bam::Read as BamRead};
use std::io::{BufWriter, Write};

/// Statistics for multi-sample remapping
#[derive(Debug, Default, Clone)]
pub struct MultiSampleRemapStats {
    pub pairs_processed: usize,
    pub pairs_with_variants: usize,
    pub haplotypes_generated: usize,
    pub reads_discarded: usize,
}

/// Remap reads for a chromosome with multi-sample support
///
/// This is the multi-sample equivalent of `swap_alleles_for_chrom` in bam_remapper.rs.
/// Uses the unique haplotype column logic to match Python's `swap_chrom_alleles_multi`.
///
/// # Arguments
/// * `bam_path` - Path to BAM file
/// * `variants` - Multi-sample variants from `parse_intersect_bed_multi`
/// * `chrom` - Chromosome to process
/// * `out_r1` - Output FASTQ path for R1
/// * `out_r2` - Output FASTQ path for R2
/// * `max_seqs` - Maximum sequences to generate per read pair
///
/// # Returns
/// (pairs_processed, haplotypes_generated)
pub fn swap_alleles_for_chrom_multi(
    bam_path: &str,
    variants: &MultiSampleVariants,
    chrom: &str,
    out_r1: &str,
    out_r2: &str,
    max_seqs: usize,
) -> Result<MultiSampleRemapStats> {
    use rustc_hash::FxHashMap;

    let mut bam = bam::IndexedReader::from_path(bam_path).context("Failed to open BAM file")?;

    // Enable parallel BGZF decompression (2 threads per chromosome worker)
    bam.set_threads(2).ok();

    let mut stats = MultiSampleRemapStats::default();

    // Get chromosome tid
    let header = bam.header().clone();
    let tid = match header.tid(chrom.as_bytes()) {
        Some(t) => t,
        None => {
            eprintln!("  Chromosome {} not found in BAM, skipping", chrom);
            return Ok(stats);
        }
    };

    bam.fetch(tid as i32)
        .context("Failed to fetch chromosome")?;

    // Open output files
    let r1_file = std::fs::File::create(out_r1).context("Failed to create R1 output file")?;
    let r2_file = std::fs::File::create(out_r2).context("Failed to create R2 output file")?;
    let mut r1_writer = BufWriter::with_capacity(1024 * 1024, r1_file);
    let mut r2_writer = BufWriter::with_capacity(1024 * 1024, r2_file);

    // Pair reads using HashMap
    let mut read_dict: FxHashMap<Vec<u8>, bam::Record> = FxHashMap::default();

    for result in bam.records() {
        let read = result.context("Failed to read BAM record")?;

        // Filter: proper pairs only, no secondary/supplementary
        if !read.is_proper_pair() || read.is_secondary() || read.is_supplementary() {
            stats.reads_discarded += 1;
            continue;
        }

        let read_name = read.qname().to_vec();

        if let Some(mate) = read_dict.remove(&read_name) {
            stats.pairs_processed += 1;

            // Determine R1 and R2
            let (read1, read2) = if read.is_first_in_template() {
                (read, mate)
            } else {
                (mate, read)
            };

            // Process this pair
            process_read_pair_multi(
                &read1,
                &read2,
                variants,
                &mut r1_writer,
                &mut r2_writer,
                &mut stats,
                max_seqs,
            )?;
        } else {
            read_dict.insert(read_name, read);
        }
    }

    stats.reads_discarded += read_dict.len();

    r1_writer.flush()?;
    r2_writer.flush()?;

    Ok(stats)
}

/// Process a read pair for multi-sample remapping (CIGAR-aware)
///
/// Uses `generate_multi_sample_sequences_from_record` which properly handles
/// reads with insertions/deletions in their CIGAR string.
fn process_read_pair_multi<W: Write>(
    read1: &bam::Record,
    read2: &bam::Record,
    variants: &MultiSampleVariants,
    r1_writer: &mut BufWriter<W>,
    r2_writer: &mut BufWriter<W>,
    stats: &mut MultiSampleRemapStats,
    max_seqs: usize,
) -> Result<()> {
    let read_name = read1.qname();

    // Look up variants for this read
    let read_variants = match variants.get(read_name) {
        Some(v) => v,
        None => return Ok(()), // No variants for this read
    };

    stats.pairs_with_variants += 1;

    // Separate variants by mate
    let r1_variants: Vec<&VariantSpanMulti> =
        read_variants.iter().filter(|v| v.mate == 1).collect();

    let r2_variants: Vec<&VariantSpanMulti> =
        read_variants.iter().filter(|v| v.mate == 2).collect();

    // Get original sequences for comparison
    let r1_seq = read1.seq().as_bytes();
    let r1_qual = read1.qual().to_vec();
    let r2_seq = read2.seq().as_bytes();
    let r2_qual = read2.qual().to_vec();

    // Generate unique haplotype sequences for R1 using CIGAR-aware mapping
    let r1_haps = if !r1_variants.is_empty() {
        // Use the new CIGAR-aware function that takes the BAM record
        generate_multi_sample_sequences_from_record(read1, &r1_variants)?
    } else {
        // No variants - use original for all haplotypes
        let num_haps = if !r2_variants.is_empty() {
            generate_unique_combinations(&r2_variants).len()
        } else {
            1
        };
        vec![(r1_seq.clone(), r1_qual.clone()); num_haps]
    };

    // Generate unique haplotype sequences for R2 using CIGAR-aware mapping
    let r2_haps = if !r2_variants.is_empty() {
        // Use the new CIGAR-aware function that takes the BAM record
        generate_multi_sample_sequences_from_record(read2, &r2_variants)?
    } else {
        vec![(r2_seq.clone(), r2_qual.clone()); r1_haps.len()]
    };

    // Ensure same number of haplotypes (use minimum)
    let num_haps = r1_haps.len().min(r2_haps.len()).min(max_seqs);

    // Get positions for WASP naming
    let r1_pos = read1.pos() as u32;
    let r2_pos = read2.pos() as u32;

    // Write pairs where at least one sequence differs from original
    let mut write_num = 0;
    let mut pairs_to_write = Vec::new();

    for (idx, ((r1_hap_seq, r1_hap_qual), (r2_hap_seq, r2_hap_qual))) in r1_haps
        .iter()
        .zip(r2_haps.iter())
        .take(num_haps)
        .enumerate()
    {
        // Skip if both sequences are unchanged
        if r1_hap_seq == &r1_seq && r2_hap_seq == &r2_seq {
            continue;
        }
        pairs_to_write.push((idx, r1_hap_seq, r1_hap_qual, r2_hap_seq, r2_hap_qual));
    }

    let write_total = pairs_to_write.len();

    for (_, r1_hap_seq, r1_hap_qual, r2_hap_seq, r2_hap_qual) in pairs_to_write {
        write_num += 1;
        stats.haplotypes_generated += 2;

        // Generate WASP read name
        let new_name = format!(
            "{}_WASP_{}_{}_{}_{}",
            String::from_utf8_lossy(read_name),
            r1_pos,
            r2_pos,
            write_num,
            write_total
        );

        // Write R1 FASTQ
        write_fastq_record(r1_writer, &new_name, r1_hap_seq, r1_hap_qual)?;

        // Write R2 FASTQ
        write_fastq_record(r2_writer, &new_name, r2_hap_seq, r2_hap_qual)?;
    }

    Ok(())
}

/// Write a FASTQ record
fn write_fastq_record<W: Write>(
    writer: &mut BufWriter<W>,
    name: &str,
    seq: &[u8],
    qual: &[u8],
) -> Result<()> {
    writeln!(writer, "@{}", name)?;
    writer.write_all(seq)?;
    writeln!(writer)?;
    writeln!(writer, "+")?;
    // Convert quality scores to ASCII (Phred+33)
    let qual_ascii: Vec<u8> = qual.iter().map(|q| q + 33).collect();
    writer.write_all(&qual_ascii)?;
    writeln!(writer)?;
    Ok(())
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    fn make_test_variant(vcf_start: u32, sample_alleles: Vec<(&str, &str)>) -> VariantSpanMulti {
        VariantSpanMulti {
            chrom: "chr1".to_string(),
            start: 0,
            stop: 100,
            vcf_start,
            vcf_stop: vcf_start + 1,
            mate: 1,
            ref_allele: "A".to_string(),
            alt_allele: "G".to_string(),
            sample_alleles: sample_alleles
                .into_iter()
                .map(|(a, b)| (a.to_string(), b.to_string()))
                .collect(),
        }
    }

    #[test]
    fn test_generate_unique_haplotype_columns_single_variant() {
        // Two samples at one position: Sample1=A|G, Sample2=A|T
        // Columns: col0=A, col1=G, col2=A, col3=T
        // Unique signatures: "A" (col0, col2), "G" (col1), "T" (col3)
        // After dedup: col0=A, col1=G, col3=T (3 unique)
        let variant = make_test_variant(10, vec![("A", "G"), ("A", "T")]);
        let variants: Vec<&VariantSpanMulti> = vec![&variant];

        let unique_cols = generate_unique_haplotype_columns(&variants);

        // 4 columns (2 samples * 2), but "A" appears twice, so 3 unique
        assert_eq!(unique_cols.len(), 3);

        let allele_sets: HashSet<Vec<String>> = unique_cols.into_iter().map(|(_, a)| a).collect();
        assert!(allele_sets.contains(&vec!["A".to_string()]));
        assert!(allele_sets.contains(&vec!["G".to_string()]));
        assert!(allele_sets.contains(&vec!["T".to_string()]));
    }

    #[test]
    fn test_generate_unique_haplotype_columns_two_variants_same_pattern() {
        // Two samples, two variants
        // Sample1: pos10=A|G, pos20=C|T  → col0="AC", col1="GT"
        // Sample2: pos10=A|G, pos20=C|T  → col2="AC", col3="GT" (same as Sample1!)
        // Unique: only 2 patterns ("AC" and "GT")
        let v1 = make_test_variant(10, vec![("A", "G"), ("A", "G")]);
        let v2 = make_test_variant(20, vec![("C", "T"), ("C", "T")]);

        let variants: Vec<&VariantSpanMulti> = vec![&v1, &v2];

        let unique_cols = generate_unique_haplotype_columns(&variants);

        // Only 2 unique column patterns (not 4!)
        assert_eq!(unique_cols.len(), 2);

        let allele_sets: HashSet<Vec<String>> = unique_cols.into_iter().map(|(_, a)| a).collect();
        assert!(allele_sets.contains(&vec!["A".to_string(), "C".to_string()]));
        assert!(allele_sets.contains(&vec!["G".to_string(), "T".to_string()]));
    }

    #[test]
    fn test_generate_unique_haplotype_columns_different_patterns() {
        // Two samples, two variants with different patterns
        // Sample1: pos10=A|G, pos20=C|T  → col0="AC", col1="GT"
        // Sample2: pos10=A|A, pos20=C|C  → col2="AC", col3="AC"
        // Unique: "AC" (col0,2,3), "GT" (col1) = 2 unique
        let v1 = make_test_variant(10, vec![("A", "G"), ("A", "A")]);
        let v2 = make_test_variant(20, vec![("C", "T"), ("C", "C")]);

        let variants: Vec<&VariantSpanMulti> = vec![&v1, &v2];

        let unique_cols = generate_unique_haplotype_columns(&variants);

        // 2 unique patterns
        assert_eq!(unique_cols.len(), 2);

        let allele_sets: HashSet<Vec<String>> = unique_cols.into_iter().map(|(_, a)| a).collect();
        assert!(allele_sets.contains(&vec!["A".to_string(), "C".to_string()]));
        assert!(allele_sets.contains(&vec!["G".to_string(), "T".to_string()]));
    }

    #[test]
    fn test_generate_unique_combinations_wrapper() {
        // Same as test_generate_unique_haplotype_columns_single_variant
        let variant = make_test_variant(10, vec![("A", "G"), ("A", "T")]);
        let variants: Vec<&VariantSpanMulti> = vec![&variant];

        let combos = generate_unique_combinations(&variants);

        assert_eq!(combos.len(), 3);

        let combo_set: HashSet<Vec<String>> = combos.into_iter().collect();
        assert!(combo_set.contains(&vec!["A".to_string()]));
        assert!(combo_set.contains(&vec!["G".to_string()]));
        assert!(combo_set.contains(&vec!["T".to_string()]));
    }

    #[test]
    fn test_apply_snp_substitution() {
        let variant = make_test_variant(5, vec![("A", "G")]);
        let variants: Vec<&VariantSpanMulti> = vec![&variant];

        let seq = b"AAAAAAAAA".to_vec(); // Position 5 is 'A'
        let qual = vec![30; 9];
        let alleles = vec!["G".to_string()];

        let (new_seq, _new_qual) =
            apply_allele_substitutions(&seq, &qual, &variants, &alleles, 0).unwrap();

        assert_eq!(&new_seq, b"AAAAAGAAA"); // Position 5 changed to G
    }

    #[test]
    fn test_generate_multi_sample_sequences() {
        let variant = make_test_variant(2, vec![("A", "G"), ("A", "T")]);
        let variants: Vec<&VariantSpanMulti> = vec![&variant];

        let seq = b"AAAAAAA".to_vec();
        let qual = vec![30; 7];

        #[allow(deprecated)]
        let results = generate_multi_sample_sequences(&seq, &qual, &variants, 0).unwrap();

        // Should have 3 unique sequences (unique columns: A, G, T)
        assert_eq!(results.len(), 3);

        let seqs: HashSet<Vec<u8>> = results.into_iter().map(|(s, _)| s).collect();
        assert!(seqs.contains(&b"AAAAAAA".to_vec())); // A at pos 2
        assert!(seqs.contains(&b"AAGAAAA".to_vec())); // G at pos 2
        assert!(seqs.contains(&b"AATAAAA".to_vec())); // T at pos 2
    }

    // ========================================================================
    // CIGAR-Aware INDEL Tests
    // ========================================================================

    fn make_position_maps(
        positions: &[(i64, usize)],
    ) -> (FxHashMap<i64, usize>, FxHashMap<i64, usize>) {
        let left: FxHashMap<i64, usize> = positions.iter().cloned().collect();
        let right: FxHashMap<i64, usize> = positions.iter().cloned().collect();
        (left, right)
    }

    #[test]
    fn test_cigar_aware_snp_substitution() {
        // Test SNP substitution with CIGAR-aware function
        let mut variant = make_test_variant(5, vec![("A", "G")]);
        variant.ref_allele = "A".to_string();
        variant.alt_allele = "G".to_string();
        variant.vcf_stop = 6; // end = start + 1 for SNP
        let variants: Vec<&VariantSpanMulti> = vec![&variant];

        let seq = b"AAAAAAAAA".to_vec();
        let qual = vec![30; 9];
        let alleles = vec!["G".to_string()];

        // Create position maps: simple 1:1 mapping (no CIGAR complexity)
        let (ref2q_left, ref2q_right) = make_position_maps(&[
            (0, 0),
            (1, 1),
            (2, 2),
            (3, 3),
            (4, 4),
            (5, 5),
            (6, 6),
            (7, 7),
            (8, 8),
        ]);

        let (new_seq, new_qual) = apply_allele_substitutions_cigar_aware(
            &seq,
            &qual,
            &variants,
            &alleles,
            &ref2q_left,
            &ref2q_right,
        )
        .unwrap();

        assert_eq!(&new_seq, b"AAAAAGAAA"); // Position 5 changed to G
        assert_eq!(new_qual.len(), 9); // Same length
    }

    #[test]
    fn test_cigar_aware_deletion_substitution() {
        // Test deletion: ACG -> A (remove 2 bases)
        let mut variant = make_test_variant(3, vec![("ACG", "A")]);
        variant.ref_allele = "ACG".to_string();
        variant.alt_allele = "A".to_string();
        variant.vcf_stop = 6; // end = start + 3
        let variants: Vec<&VariantSpanMulti> = vec![&variant];

        // Sequence: AAACGAAAA (9 bases)
        //              ^^^ variant at positions 3-5
        let seq = b"AAACGAAAA".to_vec();
        let qual = vec![30; 9];
        let alleles = vec!["A".to_string()]; // Delete CG

        // Simple 1:1 position mapping
        let (ref2q_left, ref2q_right) = make_position_maps(&[
            (0, 0),
            (1, 1),
            (2, 2),
            (3, 3),
            (4, 4),
            (5, 5),
            (6, 6),
            (7, 7),
            (8, 8),
        ]);

        let (new_seq, new_qual) = apply_allele_substitutions_cigar_aware(
            &seq,
            &qual,
            &variants,
            &alleles,
            &ref2q_left,
            &ref2q_right,
        )
        .unwrap();

        // After deletion: AAA + A + AAAA = AAAAAAA (7 bases)
        assert_eq!(&new_seq, b"AAAAAAA");
        assert_eq!(new_qual.len(), 7);
    }

    #[test]
    fn test_cigar_aware_insertion_substitution() {
        // Test insertion: A -> ACGT (insert 3 bases)
        let mut variant = make_test_variant(3, vec![("A", "ACGT")]);
        variant.ref_allele = "A".to_string();
        variant.alt_allele = "ACGT".to_string();
        variant.vcf_stop = 4; // end = start + 1
        let variants: Vec<&VariantSpanMulti> = vec![&variant];

        // Sequence: AAAAAAA (7 bases, positions 0-6)
        let seq = b"AAAAAAA".to_vec();
        let qual = vec![30; 7];
        let alleles = vec!["ACGT".to_string()]; // Replace A with ACGT

        // Simple 1:1 position mapping
        let (ref2q_left, ref2q_right) =
            make_position_maps(&[(0, 0), (1, 1), (2, 2), (3, 3), (4, 4), (5, 5), (6, 6)]);

        let (new_seq, new_qual) = apply_allele_substitutions_cigar_aware(
            &seq,
            &qual,
            &variants,
            &alleles,
            &ref2q_left,
            &ref2q_right,
        )
        .unwrap();

        // Segmentation:
        // - Before (pos 0-2): "AAA" (3 chars)
        // - Variant (pos 3): "A" -> replaced with "ACGT" (4 chars)
        // - After (pos 4-6): "AAA" (3 chars)
        // Final: "AAA" + "ACGT" + "AAA" = "AAAACGTAAA" (10 chars)
        assert_eq!(&new_seq, b"AAAACGTAAA");
        assert_eq!(new_qual.len(), 10);

        // Check that quality scores for inserted bases are Q30 (default)
        // Original qual at pos 3 goes to new pos 3, extra bases at 4, 5, 6
        assert_eq!(new_qual[4], 30); // C quality (extra)
        assert_eq!(new_qual[5], 30); // G quality (extra)
        assert_eq!(new_qual[6], 30); // T quality (extra)
    }

    #[test]
    fn test_cigar_aware_with_deletion_in_cigar() {
        // Simulate a read with a 2bp deletion in CIGAR at position 5-6
        // Read sequence: AAAAABBBBB (10 bp)
        // Reference:     AAAAA--BBBBB (positions 0-4, skip 5-6, then 7-11)
        //
        // For a variant at ref position 7, the query position should be 5 (not 7!)

        let mut variant = make_test_variant(7, vec![("B", "X")]);
        variant.ref_allele = "B".to_string();
        variant.alt_allele = "X".to_string();
        variant.vcf_stop = 8;
        let variants: Vec<&VariantSpanMulti> = vec![&variant];

        // Read sequence (no gap - deletions are in reference, not read)
        let seq = b"AAAAABBBBB".to_vec();
        let qual = vec![30; 10];
        let alleles = vec!["X".to_string()];

        // Position mapping accounting for deletion at ref 5-6
        // ref 0-4 -> query 0-4 (1:1)
        // ref 5-6 -> deleted (mapped to flanking: 4 for left, 5 for right)
        // ref 7-11 -> query 5-9 (shifted by 2)
        let (ref2q_left, ref2q_right) = make_position_maps(&[
            (0, 0),
            (1, 1),
            (2, 2),
            (3, 3),
            (4, 4),
            // ref 5-6 would be deleted - but we need them for flanking
            (7, 5),
            (8, 6),
            (9, 7),
            (10, 8),
            (11, 9),
        ]);

        let (new_seq, new_qual) = apply_allele_substitutions_cigar_aware(
            &seq,
            &qual,
            &variants,
            &alleles,
            &ref2q_left,
            &ref2q_right,
        )
        .unwrap();

        // The variant at ref 7 should map to query position 5
        // So sequence should be AAAAAXBBBB
        assert_eq!(&new_seq, b"AAAAAXBBBB");
        assert_eq!(new_qual.len(), 10);
    }

    #[test]
    fn test_cigar_aware_multiple_variants() {
        // Two SNPs at ref positions 2 and 6
        let mut v1 = make_test_variant(2, vec![("A", "G")]);
        v1.ref_allele = "A".to_string();
        v1.alt_allele = "G".to_string();
        v1.vcf_stop = 3;

        let mut v2 = make_test_variant(6, vec![("A", "T")]);
        v2.ref_allele = "A".to_string();
        v2.alt_allele = "T".to_string();
        v2.vcf_stop = 7;

        let variants: Vec<&VariantSpanMulti> = vec![&v1, &v2];

        let seq = b"AAAAAAAAA".to_vec();
        let qual = vec![30; 9];
        let alleles = vec!["G".to_string(), "T".to_string()];

        let (ref2q_left, ref2q_right) = make_position_maps(&[
            (0, 0),
            (1, 1),
            (2, 2),
            (3, 3),
            (4, 4),
            (5, 5),
            (6, 6),
            (7, 7),
            (8, 8),
        ]);

        let (new_seq, new_qual) = apply_allele_substitutions_cigar_aware(
            &seq,
            &qual,
            &variants,
            &alleles,
            &ref2q_left,
            &ref2q_right,
        )
        .unwrap();

        // Positions 2 and 6 changed
        assert_eq!(&new_seq, b"AAGAAATAA");
        assert_eq!(new_qual.len(), 9);
    }
}
