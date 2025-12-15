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
//! Uses CIGAR-walk coordinate mapping (no per-base aligned-pairs expansion),
//! properly handling reads with insertions/deletions in their alignment.

use anyhow::{Context, Result};
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::{bam, bam::Read as BamRead};
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

// ============================================================================
// Data Structures
// ============================================================================

fn complement_base(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'a' => b't',
        b'c' => b'g',
        b'g' => b'c',
        b't' => b'a',
        _ => b'N',
    }
}

fn reverse_complement_in_place(seq: &mut [u8]) {
    seq.reverse();
    for b in seq.iter_mut() {
        *b = complement_base(*b);
    }
}

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

/// Lightweight view of a variant span for allele swapping.
///
/// `generate_haplotype_seqs()` only needs the VCF coordinates and haplotype alleles,
/// so the unified pipeline can avoid per-read `String` allocations by using this
/// borrowed form.
#[derive(Debug, Clone, Copy)]
pub struct VariantSpanView<'a> {
    /// VCF variant start position (genomic coordinates)
    pub vcf_start: u32,
    /// VCF variant end position (genomic coordinates, exclusive)
    pub vcf_stop: u32,
    /// Haplotype 1 allele (phased genotype)
    pub hap1: &'a str,
    /// Haplotype 2 allele (phased genotype)
    pub hap2: &'a str,
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
// INDEL Length-Preserving Trim Structures (Phase 1 of INDEL fix)
// ============================================================================

/// Represents a single trim combination for length-preserving INDEL handling
///
/// When processing INDELs, the swapped allele may change the read length.
/// For an N-bp insertion, we need to trim N bases to restore original length.
/// This struct represents one way to distribute the trim between left and right ends.
///
/// # Example
/// For a 2bp insertion, we generate 3 combinations:
/// - TrimCombination { trim_left: 0, trim_right: 2 }  // All from right
/// - TrimCombination { trim_left: 1, trim_right: 1 }  // Split evenly
/// - TrimCombination { trim_left: 2, trim_right: 0 }  // All from left
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct TrimCombination {
    /// Bases to trim from left (5') end of the read
    pub trim_left: usize,
    /// Bases to trim from right (3') end of the read
    pub trim_right: usize,
}

impl TrimCombination {
    /// Create a new trim combination
    pub fn new(trim_left: usize, trim_right: usize) -> Self {
        Self {
            trim_left,
            trim_right,
        }
    }

    /// Total bases trimmed (should equal the INDEL delta)
    pub fn total_trim(&self) -> usize {
        self.trim_left + self.trim_right
    }

    /// Check if this is an identity (no-op) trim
    pub fn is_identity(&self) -> bool {
        self.trim_left == 0 && self.trim_right == 0
    }
}

/// Configuration for INDEL-aware remapping
#[derive(Debug, Clone)]
pub struct IndelConfig {
    /// Maximum INDEL size to process (default: 50bp)
    /// INDELs larger than this are skipped to avoid combinatorial explosion
    pub max_indel_size: usize,
    /// Whether to skip reads with large INDELs (vs failing)
    pub skip_large_indels: bool,
}

impl Default for IndelConfig {
    fn default() -> Self {
        Self {
            max_indel_size: 50,
            skip_large_indels: true,
        }
    }
}

// ============================================================================
// Main API Functions
// ============================================================================

/// Parse intersection BED file into variant HashMap
///
/// Replaces Python's `make_intersect_df()` with fast streaming parser.
/// Deduplicates exact duplicate overlaps on (chrom, read, mate, vcf_start, vcf_stop).
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
    let file =
        File::open(intersect_bed.as_ref()).context("Failed to open intersection BED file")?;
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
        let start = fields[1]
            .parse::<u32>()
            .context("Failed to parse start position")?;
        let stop = fields[2]
            .parse::<u32>()
            .context("Failed to parse stop position")?;
        let read_with_mate = fields[3]; // e.g., "SRR891276.10516353/2"
        let vcf_start = fields[7]
            .parse::<u32>()
            .context("Failed to parse VCF start position")?;
        let vcf_stop = fields[8]
            .parse::<u32>()
            .context("Failed to parse VCF stop position")?;
        let genotype = fields[11]; // e.g., "C|T"

        // Extract read name and mate
        let parts: Vec<&str> = read_with_mate.split('/').collect();
        if parts.len() != 2 {
            continue; // Skip malformed read names
        }
        let read_name = parts[0].as_bytes().to_vec();
        let mate = parts[1]
            .parse::<u8>()
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

    // Deduplicate exact duplicates on the variant span for each read/mate.
    // We'll use a HashSet to track seen combinations
    let mut seen: std::collections::HashSet<(Vec<u8>, String, u32, u32, u8)> =
        std::collections::HashSet::new();
    let mut deduped_spans: Vec<(Vec<u8>, VariantSpan)> = Vec::new();

    for (read_name, span) in all_spans {
        let key = (
            read_name.clone(),
            span.chrom.clone(),
            span.vcf_start,
            span.vcf_stop,
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
        variants
            .entry(read_name)
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
    let file =
        File::open(intersect_bed.as_ref()).context("Failed to open intersection BED file")?;
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
        let start = fields[1]
            .parse::<u32>()
            .context("Failed to parse start position")?;
        let stop = fields[2]
            .parse::<u32>()
            .context("Failed to parse stop position")?;
        let read_with_mate = fields[3];
        let vcf_start = fields[7]
            .parse::<u32>()
            .context("Failed to parse VCF start position")?;
        let vcf_stop = fields[8]
            .parse::<u32>()
            .context("Failed to parse VCF stop position")?;
        let genotype = fields[11];

        let parts: Vec<&str> = read_with_mate.split('/').collect();
        if parts.len() != 2 {
            continue;
        }
        let read_name = parts[0].as_bytes().to_vec();
        let mate = parts[1]
            .parse::<u8>()
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

    // Deduplicate exact duplicates on the variant span for each read/mate.
    let mut seen: std::collections::HashSet<(String, Vec<u8>, u32, u32, u8)> =
        std::collections::HashSet::new();
    let mut deduped_spans: Vec<(String, Vec<u8>, VariantSpan)> = Vec::new();

    for (chrom, read_name, span) in all_spans {
        let key = (
            chrom.clone(),
            read_name.clone(),
            span.vcf_start,
            span.vcf_stop,
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
    let mut bam = bam::IndexedReader::from_path(bam_path).context("Failed to open BAM file")?;

    // Enable parallel BGZF decompression (2 threads per chromosome worker)
    bam.set_threads(2).ok();

    let mut results = Vec::new();
    let mut stats = RemapStats::default();

    // Fetch reads for this chromosome
    // Use tid and fetch entire chromosome
    let header = bam.header().clone();
    let tid = header
        .tid(chrom.as_bytes())
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
            if let Some(pair_results) =
                process_read_pair(&read1, &read2, variants, config, &mut stats)?
            {
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
    let r1_variants: Vec<&VariantSpan> = read_variants.iter().filter(|v| v.mate == 1).collect();

    let r2_variants: Vec<&VariantSpan> = read_variants.iter().filter(|v| v.mate == 2).collect();

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

    for (hap_idx, ((r1_seq, r1_qual), (r2_seq, r2_qual))) in
        r1_haps.iter().zip(r2_haps.iter()).enumerate()
    {
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
        let mut r1_seq_out = r1_seq.clone();
        let mut r1_qual_out = r1_qual.clone();
        if read1.is_reverse() {
            reverse_complement_in_place(&mut r1_seq_out);
            r1_qual_out.reverse();
        }
        haplotype_reads.push(HaplotypeRead {
            name: r1_name,
            sequence: r1_seq_out,
            quals: r1_qual_out, // NOW USES INDEL-ADJUSTED QUALITIES
            original_pos: (r1_pos, r2_pos),
            haplotype: (hap_idx + 1) as u8,
        });

        // Create R2 HaplotypeRead with indel-adjusted qualities
        let r2_name = [base_name.as_slice(), b"/2"].concat();
        let mut r2_seq_out = r2_seq.clone();
        let mut r2_qual_out = r2_qual.clone();
        if read2.is_reverse() {
            reverse_complement_in_place(&mut r2_seq_out);
            r2_qual_out.reverse();
        }
        haplotype_reads.push(HaplotypeRead {
            name: r2_name,
            sequence: r2_seq_out,
            quals: r2_qual_out, // NOW USES INDEL-ADJUSTED QUALITIES
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
/// - Indels: CIGAR-walk boundary mapping (no aligned_pairs_full)
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
        // Indel-aware mapping: map BED half-open coordinates [start, stop) to query positions.
        // This matches Python’s remap_utils.py behavior:
        //   query_start = ref2q_left[start]
        //   query_stop  = ref2q_right[stop]
        let mut seq_pos = vec![0];
        let mut qual_pos = vec![0];

        for variant in variants {
            let read_start = match find_query_boundary(read, variant.vcf_start) {
                Some(pos) => pos,
                None => return Ok(None), // Variant overlaps unmapped region (e.g. splice), skip
            };
            let read_stop = match find_query_boundary(read, variant.vcf_stop) {
                Some(pos) => pos,
                None => return Ok(None),
            };

            // Skip reads where variant positions are inverted (complex CIGAR or overlapping variants)
            if read_start > read_stop {
                return Ok(None);
            }

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

            // Skip reads where variant positions are inverted (complex CIGAR or overlapping variants)
            if read_start > read_stop {
                return Ok(None);
            }

            positions.push(read_start);
            positions.push(read_stop + 1);
        }

        positions.push(original_seq.len());
        (positions.clone(), positions)
    };

    // Validate positions are monotonically increasing (overlapping variants or complex CIGARs can cause issues)
    for i in 1..split_positions.len() {
        if split_positions[i] < split_positions[i - 1] {
            return Ok(None); // Skip reads with overlapping or out-of-order variant positions
        }
    }
    for i in 1..split_qual_positions.len() {
        if split_qual_positions[i] < split_qual_positions[i - 1] {
            return Ok(None);
        }
    }

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
                let right_qual = if i < split_qual.len() - 1 {
                    split_qual[i + 1]
                } else {
                    &[]
                };

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
                let right_qual = if i < split_qual.len() - 1 {
                    split_qual[i + 1]
                } else {
                    &[]
                };

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

pub fn generate_haplotype_seqs_view(
    read: &bam::Record,
    variants: &[VariantSpanView<'_>],
    _config: &RemapConfig,
) -> Result<Option<Vec<(Vec<u8>, Vec<u8>)>>> {
    if variants.is_empty() {
        let seq = read.seq().as_bytes();
        let qual = read.qual().to_vec();
        return Ok(Some(vec![(seq.clone(), qual.clone()), (seq, qual)]));
    }

    let original_seq = read.seq().as_bytes();
    let original_qual = read.qual();

    let has_indels = variants.iter().any(|v| {
        let ref_len = (v.vcf_stop - v.vcf_start) as usize;
        v.hap1.len() != ref_len || v.hap2.len() != ref_len
    });

    // Fast path (common case): no INDEL variants AND the mapped query slice length matches allele length.
    // This avoids splitting/allocating segment vectors for SNVs/MNPs.
    if !has_indels {
        // Precompute all query ranges; fall back to slow path if any mapping is odd (e.g., read CIGAR indel
        // within the variant span causing query_len != ref_len).
        let mut edits: Vec<(usize, usize, &[u8], &[u8])> = Vec::with_capacity(variants.len());
        let mut prev_end: usize = 0;

        let mut can_fast = true;
        for v in variants {
            if v.vcf_stop <= v.vcf_start {
                can_fast = false;
                break;
            }
            let start = match find_read_position(read, v.vcf_start) {
                Some(pos) => pos,
                None => return Ok(None),
            };
            let stop_inclusive = match find_read_position(read, v.vcf_stop - 1) {
                Some(pos) => pos,
                None => return Ok(None),
            };
            let stop = stop_inclusive + 1;

            if start >= stop || stop > original_seq.len() {
                return Ok(None);
            }
            if start < prev_end {
                can_fast = false;
                break;
            }

            let a1 = v.hap1.as_bytes();
            let a2 = v.hap2.as_bytes();
            let span_len = stop - start;
            if a1.len() != span_len || a2.len() != span_len {
                can_fast = false;
                break;
            }

            edits.push((start, stop, a1, a2));
            prev_end = stop;
        }

        if can_fast {
            let mut hap1_seq = original_seq.to_vec();
            let mut hap2_seq = original_seq.to_vec();
            for (start, stop, a1, a2) in edits {
                hap1_seq[start..stop].copy_from_slice(a1);
                hap2_seq[start..stop].copy_from_slice(a2);
            }
            let qual = original_qual.to_vec();
            return Ok(Some(vec![(hap1_seq, qual.clone()), (hap2_seq, qual)]));
        }
    }

    let (split_positions, split_qual_positions) = if has_indels {
        let mut seq_pos = vec![0];
        let mut qual_pos = vec![0];

        for variant in variants {
            let read_start = match find_query_boundary(read, variant.vcf_start) {
                Some(pos) => pos,
                None => return Ok(None),
            };
            let read_stop = match find_query_boundary(read, variant.vcf_stop) {
                Some(pos) => pos,
                None => return Ok(None),
            };

            if read_start > read_stop {
                return Ok(None);
            }

            seq_pos.push(read_start);
            seq_pos.push(read_stop);
            qual_pos.push(read_start);
            qual_pos.push(read_stop);
        }

        seq_pos.push(original_seq.len());
        qual_pos.push(original_qual.len());

        (seq_pos, qual_pos)
    } else {
        let mut positions = vec![0];
        for variant in variants {
            let read_start = match find_read_position(read, variant.vcf_start) {
                Some(pos) => pos,
                None => return Ok(None),
            };
            let read_stop = match find_read_position(read, variant.vcf_stop - 1) {
                Some(pos) => pos,
                None => return Ok(None),
            };

            if read_start > read_stop {
                return Ok(None);
            }

            positions.push(read_start);
            positions.push(read_stop + 1);
        }

        positions.push(original_seq.len());
        (positions.clone(), positions)
    };

    for i in 1..split_positions.len() {
        if split_positions[i] < split_positions[i - 1] {
            return Ok(None);
        }
    }
    for i in 1..split_qual_positions.len() {
        if split_qual_positions[i] < split_qual_positions[i - 1] {
            return Ok(None);
        }
    }

    let mut split_seq: Vec<&[u8]> = Vec::new();
    let mut split_qual: Vec<&[u8]> = Vec::new();

    for i in 0..split_positions.len() - 1 {
        split_seq.push(&original_seq[split_positions[i]..split_positions[i + 1]]);
    }
    for i in 0..split_qual_positions.len() - 1 {
        split_qual.push(&original_qual[split_qual_positions[i]..split_qual_positions[i + 1]]);
    }

    let mut hap1_seq_parts: Vec<Vec<u8>> = Vec::new();
    let mut hap1_qual_parts: Vec<Vec<u8>> = Vec::new();

    for (i, seq_part) in split_seq.iter().enumerate() {
        if i % 2 == 0 {
            hap1_seq_parts.push(seq_part.to_vec());
            hap1_qual_parts.push(split_qual[i].to_vec());
        } else {
            let variant_idx = i / 2;
            let variant = &variants[variant_idx];
            let allele = variant.hap1.as_bytes();

            hap1_seq_parts.push(allele.to_vec());

            let orig_len = seq_part.len();
            let allele_len = allele.len();

            if allele_len == orig_len {
                hap1_qual_parts.push(split_qual[i].to_vec());
            } else if allele_len < orig_len {
                hap1_qual_parts.push(split_qual[i][..allele_len].to_vec());
            } else {
                let extra_len = allele_len - orig_len;
                let left_qual = if i > 0 { split_qual[i - 1] } else { &[] };
                let right_qual = if i < split_qual.len() - 1 {
                    split_qual[i + 1]
                } else {
                    &[]
                };

                let extra_quals = fill_insertion_quals(extra_len, left_qual, right_qual, 30);
                let mut combined = split_qual[i].to_vec();
                combined.extend(extra_quals);
                hap1_qual_parts.push(combined);
            }
        }
    }

    let mut hap2_seq_parts: Vec<Vec<u8>> = Vec::new();
    let mut hap2_qual_parts: Vec<Vec<u8>> = Vec::new();

    for (i, seq_part) in split_seq.iter().enumerate() {
        if i % 2 == 0 {
            hap2_seq_parts.push(seq_part.to_vec());
            hap2_qual_parts.push(split_qual[i].to_vec());
        } else {
            let variant_idx = i / 2;
            let variant = &variants[variant_idx];
            let allele = variant.hap2.as_bytes();

            hap2_seq_parts.push(allele.to_vec());

            let orig_len = seq_part.len();
            let allele_len = allele.len();

            if allele_len == orig_len {
                hap2_qual_parts.push(split_qual[i].to_vec());
            } else if allele_len < orig_len {
                hap2_qual_parts.push(split_qual[i][..allele_len].to_vec());
            } else {
                let extra_len = allele_len - orig_len;
                let left_qual = if i > 0 { split_qual[i - 1] } else { &[] };
                let right_qual = if i < split_qual.len() - 1 {
                    split_qual[i + 1]
                } else {
                    &[]
                };

                let extra_quals = fill_insertion_quals(extra_len, left_qual, right_qual, 30);
                let mut combined = split_qual[i].to_vec();
                combined.extend(extra_quals);
                hap2_qual_parts.push(combined);
            }
        }
    }

    let hap1_seq: Vec<u8> = hap1_seq_parts.into_iter().flatten().collect();
    let hap1_qual: Vec<u8> = hap1_qual_parts.into_iter().flatten().collect();
    let hap2_seq: Vec<u8> = hap2_seq_parts.into_iter().flatten().collect();
    let hap2_qual: Vec<u8> = hap2_qual_parts.into_iter().flatten().collect();

    Ok(Some(vec![(hap1_seq, hap1_qual), (hap2_seq, hap2_qual)]))
}

// ============================================================================
// INDEL Length-Preserving Trim Functions (Phase 2 of INDEL fix)
// ============================================================================

/// Generate all valid trim combinations for a given net length change
///
/// For an N-bp insertion (delta > 0), we need to trim N bases total.
/// Generates N+1 combinations: (0,N), (1,N-1), ..., (N,0)
///
/// # Arguments
/// * `indel_delta` - Net length change (positive = insertion bytes to trim)
/// * `read_len` - Original read length (to validate trim doesn't exceed)
///
/// # Returns
/// Vector of TrimCombination structs
///
/// # Examples
/// ```
/// let combos = generate_trim_combinations(2, 51);
/// assert_eq!(combos.len(), 3);  // (0,2), (1,1), (2,0)
/// ```
pub fn generate_trim_combinations(indel_delta: i32, read_len: usize) -> Vec<TrimCombination> {
    if indel_delta <= 0 {
        // Deletion or SNP: no trim needed, single "identity" combination
        return vec![TrimCombination::new(0, 0)];
    }

    let trim_needed = indel_delta as usize;

    // Safety: don't trim more than half the read from either side
    let max_trim_per_side = read_len / 2;

    let mut combinations = Vec::with_capacity(trim_needed + 1);

    for left_trim in 0..=trim_needed {
        let right_trim = trim_needed - left_trim;

        // Validate this combination is feasible (don't trim too much from either side)
        if left_trim <= max_trim_per_side && right_trim <= max_trim_per_side {
            combinations.push(TrimCombination::new(left_trim, right_trim));
        }
    }

    // Fallback for very large indels where no combination works
    if combinations.is_empty() {
        // Fall back to splitting evenly
        let half = trim_needed / 2;
        let remainder = trim_needed % 2;
        combinations.push(TrimCombination::new(half, half + remainder));
    }

    combinations
}

/// Apply trim combination to sequence and quality scores
///
/// Trims the extended sequence back to original length for insertions,
/// or pads with N's for deletions (to maintain consistent length).
///
/// # Arguments
/// * `seq` - The (possibly extended) sequence after allele swapping
/// * `qual` - The quality scores corresponding to seq
/// * `original_len` - The original read length we want to restore
/// * `trim` - Which trim combination to apply
///
/// # Returns
/// Tuple of (trimmed_sequence, trimmed_qualities) both with length = original_len
pub fn apply_trim_combination(
    seq: &[u8],
    qual: &[u8],
    original_len: usize,
    trim: &TrimCombination,
) -> (Vec<u8>, Vec<u8>) {
    let seq_len = seq.len();

    if seq_len <= original_len {
        // Deletion case: sequence is shorter or equal to original
        // Pad with N's to restore original length
        let mut padded_seq = seq.to_vec();
        let mut padded_qual = qual.to_vec();

        while padded_seq.len() < original_len {
            padded_seq.push(b'N');
            padded_qual.push(0); // Quality 0 for padded bases
        }
        return (padded_seq, padded_qual);
    }

    // Insertion case: sequence is longer than original, need to trim
    // Calculate start and end indices after trimming
    let start = trim.trim_left.min(seq_len);
    let end = seq_len.saturating_sub(trim.trim_right);
    let end = end.max(start); // Ensure end >= start

    // Extract the trimmed region
    let trimmed_seq: Vec<u8> = seq[start..end].to_vec();
    let trimmed_qual: Vec<u8> = qual[start..end.min(qual.len())].to_vec();

    // Ensure exact length (should already be correct, but safety check)
    let mut final_seq = trimmed_seq;
    let mut final_qual = trimmed_qual;

    final_seq.truncate(original_len);
    final_qual.truncate(original_len);

    // Pad if somehow still short (shouldn't happen with correct trim values)
    while final_seq.len() < original_len {
        final_seq.push(b'N');
    }
    while final_qual.len() < original_len {
        final_qual.push(0);
    }

    (final_seq, final_qual)
}

/// Calculate the INDEL delta (length change) for a haplotype sequence
///
/// # Arguments
/// * `hap_seq_len` - Length of the generated haplotype sequence
/// * `original_len` - Original read length
///
/// # Returns
/// Positive value for insertions (need to trim), negative for deletions, 0 for SNPs
#[inline]
pub fn calculate_indel_delta(hap_seq_len: usize, original_len: usize) -> i32 {
    hap_seq_len as i32 - original_len as i32
}

/// Generate haplotype sequences with trim combinations for length preservation
///
/// This is the INDEL-aware version that maintains original read length.
/// For each raw haplotype, generates multiple trimmed versions if the sequence
/// was extended by an insertion.
///
/// # Arguments
/// * `read` - BAM record
/// * `variants` - Variants overlapping this read
/// * `config` - Remapping configuration
/// * `indel_config` - INDEL handling configuration
///
/// # Returns
/// `Ok(Some(vec))` - Vector of (sequence, qualities, trim_combo_id) tuples
/// `Ok(None)` - Read should be skipped (unmappable variant position or too large INDEL)
pub fn generate_haplotype_seqs_with_trims(
    read: &bam::Record,
    variants: &[&VariantSpan],
    config: &RemapConfig,
    indel_config: &IndelConfig,
) -> Result<Option<Vec<(Vec<u8>, Vec<u8>, u16)>>> {
    let original_len = read.seq().len();

    // Check for oversized INDELs
    for variant in variants {
        let ref_len = (variant.vcf_stop - variant.vcf_start) as usize;
        let max_allele_len = variant.hap1.len().max(variant.hap2.len());
        let indel_size = (max_allele_len as i32 - ref_len as i32).unsigned_abs() as usize;

        if indel_size > indel_config.max_indel_size {
            if indel_config.skip_large_indels {
                return Ok(None); // Skip this read
            }
        }
    }

    // First, generate raw (potentially extended) haplotype sequences
    let raw_haps = match generate_haplotype_seqs(read, variants, config)? {
        Some(h) => h,
        None => return Ok(None),
    };

    let mut result: Vec<(Vec<u8>, Vec<u8>, u16)> = Vec::new();

    for (hap_idx, (raw_seq, raw_qual)) in raw_haps.iter().enumerate() {
        let indel_delta = calculate_indel_delta(raw_seq.len(), original_len);

        let trim_combos = generate_trim_combinations(indel_delta, original_len);

        for (combo_idx, trim) in trim_combos.iter().enumerate() {
            let (trimmed_seq, trimmed_qual) =
                apply_trim_combination(raw_seq, raw_qual, original_len, trim);

            // Encode: hap_idx * 1000 + combo_idx (allows up to 1000 combos per haplotype)
            let trim_combo_id = (hap_idx as u16) * 1000 + (combo_idx as u16);

            result.push((trimmed_seq, trimmed_qual, trim_combo_id));
        }
    }

    if result.is_empty() {
        Ok(None)
    } else {
        Ok(Some(result))
    }
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
        File::create(r1_path.as_ref()).context("Failed to create R1 FASTQ")?,
    );
    let mut r2_file = std::io::BufWriter::new(
        File::create(r2_path.as_ref()).context("Failed to create R2 FASTQ")?,
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
            r1_file
                .write_all(fastq_entry.as_bytes())
                .context("Failed to write R1 FASTQ entry")?;
            r1_count += 1;
        } else {
            r2_file
                .write_all(fastq_entry.as_bytes())
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
            std::fs::File::create(&r1_path_str).context("Failed to create R1 FASTQ")?,
        );
        let mut r2_file = std::io::BufWriter::new(
            std::fs::File::create(&r2_path_str).context("Failed to create R2 FASTQ")?,
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
                r1_file
                    .write_all(fastq_entry.as_bytes())
                    .context("Failed to write R1 FASTQ entry")?;
                r1_count += 1;
            } else {
                r2_file
                    .write_all(fastq_entry.as_bytes())
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
            let (haplotypes, stats) =
                swap_alleles_for_chrom(bam_path, chrom_variants, chrom, config)?;

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

/// Map a reference coordinate to a query (read) coordinate using CIGAR.
///
/// Returns the query position corresponding to the *boundary before* `target_ref_pos`
/// in the reference coordinate system, which matches the semantics used by WASP2’s
/// Python implementation for indel-aware splitting:
/// - query_start = ref2q_left[start]
/// - query_stop  = ref2q_right[stop]
///
/// We treat:
/// - `D` (deletion) as mappable using the current query position (flank)
/// - `N` (ref-skip / splice) as NOT mappable (returns None)
fn find_query_boundary(read: &bam::Record, target_ref_pos: u32) -> Option<usize> {
    use rust_htslib::bam::record::Cigar;

    let mut query_pos: usize = 0;
    let mut ref_pos: u32 = read.pos() as u32;

    for op in read.cigar().iter() {
        match op {
            Cigar::Ins(len) | Cigar::SoftClip(len) => {
                // Query advances, reference stays. This must be applied before mapping the
                // next reference-consuming operation at the same ref_pos.
                query_pos += *len as usize;
            }
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                let op_ref_len = *len;
                if target_ref_pos < ref_pos {
                    return None;
                }
                if target_ref_pos < ref_pos + op_ref_len {
                    let offset = (target_ref_pos - ref_pos) as usize;
                    return Some(query_pos + offset);
                }
                // target is at or after end of this op
                query_pos += op_ref_len as usize;
                ref_pos += op_ref_len;
            }
            Cigar::Del(len) => {
                let op_ref_len = *len;
                if target_ref_pos < ref_pos {
                    return None;
                }
                if target_ref_pos < ref_pos + op_ref_len {
                    // Inside a deletion: return flank (query doesn't advance)
                    return Some(query_pos);
                }
                ref_pos += op_ref_len;
            }
            Cigar::RefSkip(len) => {
                let op_ref_len = *len;
                if target_ref_pos < ref_pos {
                    return None;
                }
                if target_ref_pos < ref_pos + op_ref_len {
                    // Splice/intron skip: treat as unmappable
                    return None;
                }
                ref_pos += op_ref_len;
            }
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }

    // If target is exactly at the end of the reference span, return boundary at end of read.
    if target_ref_pos == ref_pos {
        Some(query_pos)
    } else {
        None
    }
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

// ============================================================================
// CIGAR-Aware Expected Position Calculation
// ============================================================================

/// Classification of a variant relative to a read's CIGAR alignment
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantLocation {
    /// Variant ends strictly before the read's reference start - shifts expected position
    Upstream,
    /// Variant overlaps the read's aligned region - no shift
    WithinRead,
    /// Variant starts after the read's reference end - no shift
    Downstream,
    /// Variant spans the read start boundary - treated as within-read (no shift)
    SpansStart,
}

/// Classify a variant's location relative to a read using CIGAR information
///
/// This uses the read's CIGAR-derived reference span to determine if a variant
/// is truly upstream (before alignment start), within the read's aligned region,
/// or downstream (after alignment end).
///
/// # Arguments
/// * `read` - BAM record with CIGAR information
/// * `variant_start` - Variant start position (0-based, reference coordinates)
/// * `variant_end` - Variant end position (0-based, exclusive, reference coordinates)
///
/// # Returns
/// `VariantLocation` classification
pub fn classify_variant_location(
    read: &bam::Record,
    variant_start: u32,
    variant_end: u32,
) -> VariantLocation {
    // Get read's reference span from alignment
    let read_ref_start = read.pos() as u32;
    let read_ref_end = read.reference_end() as u32;

    // Variant ends before read starts on reference
    if variant_end <= read_ref_start {
        return VariantLocation::Upstream;
    }

    // Variant starts after read ends on reference
    if variant_start >= read_ref_end {
        return VariantLocation::Downstream;
    }

    // Variant spans the read start boundary
    if variant_start < read_ref_start && variant_end > read_ref_start {
        return VariantLocation::SpansStart;
    }

    // Otherwise, variant is within the read's aligned region
    VariantLocation::WithinRead
}

/// Compute expected alignment position for a read after applying haplotype variants
///
/// This is CIGAR-aware: it uses the read's CIGAR-derived reference span to
/// classify variants as upstream vs within-read. Only **upstream** variants
/// (those ending strictly before the read's reference start) shift the expected
/// alignment position.
///
/// Within-read variants change the read sequence but don't change where it
/// should align on the reference.
///
/// # Arguments
/// * `read` - BAM record with CIGAR information
/// * `variants` - Iterator of (variant_start, variant_end, delta) tuples where:
///   - variant_start: 0-based reference position
///   - variant_end: 0-based exclusive end position
///   - delta: len(alt) - len(ref), positive for insertions, negative for deletions
///
/// # Returns
/// Expected alignment position (0-based) after applying upstream variant shifts
///
/// # Example
/// ```ignore
/// // Read at pos=100, upstream 5bp insertion at pos=50
/// // Expected position = 100 + 5 = 105
/// let expected = compute_expected_position_cigar_aware(&read, &[(50, 51, 5)]);
/// assert_eq!(expected, 105);
/// ```
pub fn compute_expected_position_cigar_aware<'a, I>(read: &bam::Record, variants: I) -> i64
where
    I: IntoIterator<Item = &'a (u32, u32, i32)>,
{
    let read_start = read.pos();
    let mut cumulative_shift: i64 = 0;

    for &(var_start, var_end, delta) in variants {
        let location = classify_variant_location(read, var_start, var_end);

        match location {
            VariantLocation::Upstream => {
                // Variant is fully upstream - shifts expected position
                cumulative_shift += delta as i64;
            }
            VariantLocation::SpansStart => {
                // Variant spans read start - complex case
                // For deletions spanning into the read: the read start moves
                // For insertions at boundary: treat as upstream shift
                if delta < 0 {
                    // Deletion spanning into read - shifts position
                    cumulative_shift += delta as i64;
                } else if delta > 0 && var_start < read_start as u32 {
                    // Insertion before read start - shifts position
                    cumulative_shift += delta as i64;
                }
                // SNVs at boundary: no shift
            }
            VariantLocation::WithinRead | VariantLocation::Downstream => {
                // No shift for within-read or downstream variants
            }
        }
    }

    read_start + cumulative_shift
}

/// Simplified interface for compute_expected_position_cigar_aware
///
/// Takes variants as (position, delta) pairs where position is the variant start
/// and delta is len(alt) - len(ref). Computes variant end as:
/// - For deletions (delta < 0): end = start + |delta|
/// - For insertions (delta > 0): end = start + 1 (point insertion)
/// - For SNVs (delta == 0): end = start + 1
///
/// # Arguments
/// * `read` - BAM record
/// * `variants` - Iterator of (position, delta) pairs
///
/// # Returns
/// Expected alignment position after upstream shifts
pub fn compute_expected_position<'a, I>(read: &bam::Record, variants: I) -> i64
where
    I: IntoIterator<Item = &'a (u32, i32)>,
{
    let read_start = read.pos();
    let read_ref_start = read_start as u32;
    let mut cumulative_shift: i64 = 0;

    for &(var_pos, delta) in variants {
        // Compute variant end based on delta
        let var_end = if delta < 0 {
            // Deletion: spans |delta| reference bases
            var_pos + ((-delta) as u32)
        } else {
            // Insertion or SNV: point position
            var_pos + 1
        };

        // Check if variant is upstream
        if var_end <= read_ref_start {
            // Fully upstream - shift expected position
            cumulative_shift += delta as i64;
        } else if var_pos < read_ref_start && delta < 0 {
            // Deletion spanning into read start - still shifts
            cumulative_shift += delta as i64;
        } else if var_pos < read_ref_start && delta > 0 {
            // Insertion before read start - shifts
            cumulative_shift += delta as i64;
        }
        // Within-read or downstream: no shift
    }

    read_start + cumulative_shift
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
    format!(
        "{}_WASP_{}_{}_{}_{}",
        name_str, pos1, pos2, seq_num, total_seqs
    )
    .into_bytes()
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
        writeln!(
            temp_file,
            "chr10\t87377\t87427\tSRR891276.10516353/2\t60\t+\tchr10\t87400\t87401\tC\tT\tC|T"
        )
        .unwrap();
        writeln!(
            temp_file,
            "chr10\t87392\t87440\tSRR891276.5620594/2\t60\t+\tchr10\t87400\t87401\tC\tT\tC|T"
        )
        .unwrap();
        // Second distinct variant overlap for the same read/mate (should be preserved)
        writeln!(
            temp_file,
            "chr10\t87392\t87440\tSRR891276.5620594/2\t60\t+\tchr10\t87401\t87402\tA\tG\tA|G"
        )
        .unwrap();
        writeln!(
            temp_file,
            "chr10\t87395\t87442\tSRR891276.5620594/1\t60\t-\tchr10\t87400\t87401\tC\tT\tC|T"
        )
        .unwrap();
        // Duplicate that should be removed (same read/mate + same variant span)
        writeln!(
            temp_file,
            "chr10\t87392\t87440\tSRR891276.5620594/2\t60\t+\tchr10\t87401\t87402\tA\tG\tA|G"
        )
        .unwrap();
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
        assert_eq!(
            read2_spans.len(),
            3,
            "Should keep both variant overlaps for mate 2, plus mate 1"
        );

        // Verify mate 1
        let mate1 = read2_spans.iter().find(|s| s.mate == 1).unwrap();
        assert_eq!(mate1.start, 87395);
        assert_eq!(mate1.stop, 87442);
        assert_eq!(mate1.vcf_start, 87400);
        assert_eq!(mate1.vcf_stop, 87401);

        // Verify mate 2 (should have two distinct variant overlaps; duplicate removed)
        let mate2: Vec<_> = read2_spans.iter().filter(|s| s.mate == 2).collect();
        assert_eq!(mate2.len(), 2);
        assert!(mate2.iter().any(|s| s.vcf_start == 87400 && s.vcf_stop == 87401));
        assert!(mate2.iter().any(|s| s.vcf_start == 87401 && s.vcf_stop == 87402));
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

    // ============================================================================
    // INDEL Trim Combination Tests
    // ============================================================================

    #[test]
    fn test_trim_combination_struct() {
        let trim = TrimCombination::new(2, 3);
        assert_eq!(trim.trim_left, 2);
        assert_eq!(trim.trim_right, 3);
        assert_eq!(trim.total_trim(), 5);
        assert!(!trim.is_identity());

        let identity = TrimCombination::new(0, 0);
        assert!(identity.is_identity());
    }

    #[test]
    fn test_generate_trim_combinations_2bp_insertion() {
        // 2bp insertion → need to trim 2 bases total
        // Should generate 3 combinations: (0,2), (1,1), (2,0)
        let combos = generate_trim_combinations(2, 51);
        assert_eq!(combos.len(), 3, "2bp insertion should give 3 combos");
        assert_eq!(combos[0], TrimCombination::new(0, 2));
        assert_eq!(combos[1], TrimCombination::new(1, 1));
        assert_eq!(combos[2], TrimCombination::new(2, 0));
    }

    #[test]
    fn test_generate_trim_combinations_snv() {
        // SNV (delta=0) → no trimming needed
        let combos = generate_trim_combinations(0, 51);
        assert_eq!(combos.len(), 1);
        assert_eq!(combos[0], TrimCombination::new(0, 0));
        assert!(combos[0].is_identity());
    }

    #[test]
    fn test_generate_trim_combinations_deletion() {
        // Deletion (delta=-2) → no trimming needed (padding is separate)
        let combos = generate_trim_combinations(-2, 51);
        assert_eq!(combos.len(), 1);
        assert_eq!(combos[0], TrimCombination::new(0, 0));
    }

    #[test]
    fn test_generate_trim_combinations_5bp_insertion() {
        // 5bp insertion → 6 combinations
        let combos = generate_trim_combinations(5, 51);
        assert_eq!(combos.len(), 6, "5bp insertion should give 6 combos");
        // Check all combinations sum to 5
        for combo in &combos {
            assert_eq!(combo.total_trim(), 5);
        }
    }

    #[test]
    fn test_apply_trim_combination_insertion() {
        // Original: 10bp, Extended: 12bp (2bp insertion)
        let seq = b"ACGTACGTACGT".to_vec(); // 12bp
        let qual = vec![30; 12];
        let original_len = 10;

        // Trim 1 from left, 1 from right → should get middle 10bp
        let trim = TrimCombination::new(1, 1);
        let (trimmed_seq, trimmed_qual) = apply_trim_combination(&seq, &qual, original_len, &trim);

        assert_eq!(
            trimmed_seq.len(),
            original_len,
            "Trimmed seq should match original length"
        );
        assert_eq!(
            trimmed_qual.len(),
            original_len,
            "Trimmed qual should match original length"
        );
        assert_eq!(trimmed_seq, b"CGTACGTACG".to_vec());
    }

    #[test]
    fn test_apply_trim_combination_trim_all_left() {
        // Trim all from left
        let seq = b"ACGTACGTACGT".to_vec(); // 12bp
        let qual = vec![30; 12];
        let original_len = 10;

        let trim = TrimCombination::new(2, 0);
        let (trimmed_seq, _) = apply_trim_combination(&seq, &qual, original_len, &trim);

        assert_eq!(trimmed_seq.len(), original_len);
        assert_eq!(trimmed_seq, b"GTACGTACGT".to_vec());
    }

    #[test]
    fn test_apply_trim_combination_trim_all_right() {
        // Trim all from right
        let seq = b"ACGTACGTACGT".to_vec(); // 12bp
        let qual = vec![30; 12];
        let original_len = 10;

        let trim = TrimCombination::new(0, 2);
        let (trimmed_seq, _) = apply_trim_combination(&seq, &qual, original_len, &trim);

        assert_eq!(trimmed_seq.len(), original_len);
        assert_eq!(trimmed_seq, b"ACGTACGTAC".to_vec());
    }

    #[test]
    fn test_apply_trim_combination_deletion_pads() {
        // Deletion case: seq shorter than original → should pad with N's
        let seq = b"ACGTACGT".to_vec(); // 8bp
        let qual = vec![30; 8];
        let original_len = 10;

        let trim = TrimCombination::new(0, 0); // No trim for deletions
        let (trimmed_seq, trimmed_qual) = apply_trim_combination(&seq, &qual, original_len, &trim);

        assert_eq!(trimmed_seq.len(), original_len);
        assert_eq!(trimmed_qual.len(), original_len);
        // Should be padded with N's
        assert_eq!(&trimmed_seq[8..], b"NN");
        assert_eq!(&trimmed_qual[8..], &[0, 0]);
    }

    #[test]
    fn test_calculate_indel_delta() {
        // Insertion: hap_len > original
        assert_eq!(calculate_indel_delta(53, 51), 2);

        // Deletion: hap_len < original
        assert_eq!(calculate_indel_delta(49, 51), -2);

        // SNV: hap_len == original
        assert_eq!(calculate_indel_delta(51, 51), 0);
    }

    #[test]
    fn test_indel_config_default() {
        let config = IndelConfig::default();
        assert_eq!(config.max_indel_size, 50);
        assert!(config.skip_large_indels);
    }

    // ========================================================================
    // CIGAR-Aware Expected Position Tests
    // ========================================================================

    /// Helper to create a minimal BAM record with specified pos and CIGAR
    fn create_test_record(pos: i64, cigar_str: &str) -> bam::Record {
        use rust_htslib::bam::record::{Cigar, CigarString};

        let mut rec = bam::Record::new();
        rec.set_pos(pos);

        // Parse simple CIGAR string (e.g., "50M", "10M5D10M", "5S45M")
        let mut cigars = Vec::new();
        let mut num_str = String::new();

        for c in cigar_str.chars() {
            if c.is_ascii_digit() {
                num_str.push(c);
            } else {
                let len: u32 = num_str.parse().unwrap_or(1);
                num_str.clear();
                let op = match c {
                    'M' => Cigar::Match(len),
                    'I' => Cigar::Ins(len),
                    'D' => Cigar::Del(len),
                    'S' => Cigar::SoftClip(len),
                    'N' => Cigar::RefSkip(len),
                    '=' => Cigar::Equal(len),
                    'X' => Cigar::Diff(len),
                    'H' => Cigar::HardClip(len),
                    _ => Cigar::Match(len),
                };
                cigars.push(op);
            }
        }

        let query_len: usize = cigars
            .iter()
            .map(|op| match op {
                Cigar::Match(len)
                | Cigar::Ins(len)
                | Cigar::SoftClip(len)
                | Cigar::Equal(len)
                | Cigar::Diff(len) => *len as usize,
                Cigar::Del(_) | Cigar::RefSkip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => 0,
            })
            .sum();

        let cigar_string = CigarString(cigars);
        let seq = vec![b'A'; query_len];
        let qual = vec![30u8; query_len];
        rec.set(
            b"test_read",
            Some(&cigar_string),
            &seq,  // Dummy sequence
            &qual, // Dummy qualities
        );
        rec.set_pos(pos);

        rec
    }

    #[test]
    fn test_find_query_boundary_simple_match() {
        let rec = create_test_record(100, "50M");

        assert_eq!(find_query_boundary(&rec, 100), Some(0));
        assert_eq!(find_query_boundary(&rec, 101), Some(1));
        assert_eq!(find_query_boundary(&rec, 150), Some(50)); // end boundary
        assert_eq!(find_query_boundary(&rec, 99), None);
    }

    #[test]
    fn test_find_query_boundary_softclip() {
        // 5S45M: aligned portion starts at query offset 5
        let rec = create_test_record(100, "5S45M");
        assert_eq!(find_query_boundary(&rec, 100), Some(5));
        assert_eq!(find_query_boundary(&rec, 101), Some(6));
        assert_eq!(find_query_boundary(&rec, 145), Some(50)); // 5 + 45
    }

    #[test]
    fn test_find_query_boundary_insertion_shifts_downstream() {
        // 10M2I40M: insertion occurs at ref_pos=110, pushing downstream query coords by +2
        let rec = create_test_record(100, "10M2I40M");
        assert_eq!(find_query_boundary(&rec, 109), Some(9));
        assert_eq!(find_query_boundary(&rec, 110), Some(12));
        assert_eq!(find_query_boundary(&rec, 111), Some(13));
    }

    #[test]
    fn test_find_query_boundary_deletion_keeps_query_constant() {
        // 10M2D40M: deletion consumes ref 110-111 with no query advance
        let rec = create_test_record(100, "10M2D40M");
        assert_eq!(find_query_boundary(&rec, 110), Some(10));
        assert_eq!(find_query_boundary(&rec, 111), Some(10));
        assert_eq!(find_query_boundary(&rec, 112), Some(10));
    }

    #[test]
    fn test_find_query_boundary_refskip_is_unmappable() {
        // 10M100N40M: positions within N are unmappable
        let rec = create_test_record(100, "10M100N40M");
        assert_eq!(find_query_boundary(&rec, 110), None);
        assert_eq!(find_query_boundary(&rec, 150), None);
        assert_eq!(find_query_boundary(&rec, 210), Some(10));
    }

    #[test]
    fn test_generate_haplotype_seqs_view_insertion_uses_stop_boundary() {
        // Insertion at [125,126): should replace 1 ref base with 3 bases, net +2 length
        let rec = create_test_record(100, "50M");
        let view = vec![VariantSpanView {
            vcf_start: 125,
            vcf_stop: 126,
            hap1: "A",
            hap2: "ATG",
        }];
        let cfg = RemapConfig::default();
        let out = generate_haplotype_seqs_view(&rec, &view, &cfg).unwrap().unwrap();

        assert_eq!(out[0].0.len(), 50); // hap1: ref allele
        assert_eq!(out[1].0.len(), 52); // hap2: insertion allele, replaces 1 base with 3
        assert_eq!(&out[1].0[25..28], b"ATG");
    }

    #[test]
    fn test_generate_haplotype_seqs_view_deletion_contracts_sequence() {
        // Deletion at [120,122): replaces 2 ref bases with 1 base, net -1 length
        let rec = create_test_record(100, "50M");
        let view = vec![VariantSpanView {
            vcf_start: 120,
            vcf_stop: 122,
            hap1: "AA",
            hap2: "A",
        }];
        let cfg = RemapConfig::default();
        let out = generate_haplotype_seqs_view(&rec, &view, &cfg).unwrap().unwrap();

        assert_eq!(out[0].0.len(), 50); // hap1 matches ref length
        assert_eq!(out[1].0.len(), 49); // hap2 shorter by 1
    }

    #[test]
    fn test_generate_haplotype_seqs_view_matches_owned_snp() {
        let rec = create_test_record(100, "50M");
        let owned = vec![VariantSpan {
            chrom: "chr1".to_string(),
            start: 100,
            stop: 150,
            vcf_start: 120,
            vcf_stop: 121,
            mate: 1,
            hap1: "A".to_string(),
            hap2: "G".to_string(),
        }];
        let owned_refs: Vec<&VariantSpan> = owned.iter().collect();

        let view = vec![VariantSpanView {
            vcf_start: 120,
            vcf_stop: 121,
            hap1: "A",
            hap2: "G",
        }];

        let cfg = RemapConfig::default();
        let out_owned = generate_haplotype_seqs(&rec, &owned_refs, &cfg).unwrap();
        let out_view = generate_haplotype_seqs_view(&rec, &view, &cfg).unwrap();
        assert_eq!(out_owned, out_view);
    }

    #[test]
    fn test_generate_haplotype_seqs_view_matches_owned_insertion() {
        let rec = create_test_record(100, "50M");
        let owned = vec![VariantSpan {
            chrom: "chr1".to_string(),
            start: 100,
            stop: 150,
            vcf_start: 125,
            vcf_stop: 126,
            mate: 1,
            hap1: "A".to_string(),
            hap2: "ATG".to_string(), // 2bp insertion relative to ref len=1
        }];
        let owned_refs: Vec<&VariantSpan> = owned.iter().collect();

        let view = vec![VariantSpanView {
            vcf_start: 125,
            vcf_stop: 126,
            hap1: "A",
            hap2: "ATG",
        }];

        let cfg = RemapConfig::default();
        let out_owned = generate_haplotype_seqs(&rec, &owned_refs, &cfg).unwrap();
        let out_view = generate_haplotype_seqs_view(&rec, &view, &cfg).unwrap();
        assert_eq!(out_owned, out_view);
    }

    #[test]
    fn test_classify_variant_upstream() {
        // Read at pos=100 with 50M CIGAR (covers ref 100-149)
        let rec = create_test_record(100, "50M");

        // Variant at 50-51 is upstream (ends before read starts)
        let loc = classify_variant_location(&rec, 50, 51);
        assert_eq!(loc, VariantLocation::Upstream);

        // Variant at 90-99 is upstream (ends at 99, before read start at 100)
        let loc = classify_variant_location(&rec, 90, 99);
        assert_eq!(loc, VariantLocation::Upstream);

        // Variant at 90-100 is upstream (ends exactly at read start)
        let loc = classify_variant_location(&rec, 90, 100);
        assert_eq!(loc, VariantLocation::Upstream);
    }

    #[test]
    fn test_classify_variant_within_read() {
        // Read at pos=100 with 50M CIGAR (covers ref 100-149)
        let rec = create_test_record(100, "50M");

        // Variant at 110-111 is within read
        let loc = classify_variant_location(&rec, 110, 111);
        assert_eq!(loc, VariantLocation::WithinRead);

        // Variant at 100-101 is within read (at read start)
        let loc = classify_variant_location(&rec, 100, 101);
        assert_eq!(loc, VariantLocation::WithinRead);

        // Variant at 148-150 overlaps read end - still within
        let loc = classify_variant_location(&rec, 148, 150);
        assert_eq!(loc, VariantLocation::WithinRead);
    }

    #[test]
    fn test_classify_variant_downstream() {
        // Read at pos=100 with 50M CIGAR (covers ref 100-149)
        let rec = create_test_record(100, "50M");

        // Variant at 150-151 is downstream (starts at read end)
        let loc = classify_variant_location(&rec, 150, 151);
        assert_eq!(loc, VariantLocation::Downstream);

        // Variant at 200-201 is downstream
        let loc = classify_variant_location(&rec, 200, 201);
        assert_eq!(loc, VariantLocation::Downstream);
    }

    #[test]
    fn test_classify_variant_spans_start() {
        // Read at pos=100 with 50M CIGAR (covers ref 100-149)
        let rec = create_test_record(100, "50M");

        // Variant at 95-105 spans read start (starts before, ends after)
        let loc = classify_variant_location(&rec, 95, 105);
        assert_eq!(loc, VariantLocation::SpansStart);

        // Deletion from 98-102 spans read start
        let loc = classify_variant_location(&rec, 98, 102);
        assert_eq!(loc, VariantLocation::SpansStart);
    }

    #[test]
    fn test_compute_expected_position_no_variants() {
        let rec = create_test_record(100, "50M");
        let variants: Vec<(u32, i32)> = vec![];
        let expected = compute_expected_position(&rec, &variants);
        assert_eq!(expected, 100);
    }

    #[test]
    fn test_compute_expected_position_upstream_insertion() {
        // Read at pos=100, upstream 5bp insertion at pos=50
        let rec = create_test_record(100, "50M");
        let variants = vec![(50u32, 5i32)]; // 5bp insertion
        let expected = compute_expected_position(&rec, &variants);
        // Upstream insertion shifts expected position right
        assert_eq!(expected, 105);
    }

    #[test]
    fn test_compute_expected_position_upstream_deletion() {
        // Read at pos=100, upstream 3bp deletion at pos=50
        let rec = create_test_record(100, "50M");
        let variants = vec![(50u32, -3i32)]; // 3bp deletion (spans 50-52)
        let expected = compute_expected_position(&rec, &variants);
        // Upstream deletion shifts expected position left
        assert_eq!(expected, 97);
    }

    #[test]
    fn test_compute_expected_position_upstream_snv() {
        // Read at pos=100, upstream SNV at pos=50
        let rec = create_test_record(100, "50M");
        let variants = vec![(50u32, 0i32)]; // SNV (delta=0)
        let expected = compute_expected_position(&rec, &variants);
        // SNV doesn't shift position
        assert_eq!(expected, 100);
    }

    #[test]
    fn test_compute_expected_position_within_read_variants() {
        // Read at pos=100, within-read variants shouldn't shift
        let rec = create_test_record(100, "50M");

        // Insertion within read
        let variants = vec![(120u32, 5i32)];
        let expected = compute_expected_position(&rec, &variants);
        assert_eq!(expected, 100); // No shift

        // Deletion within read
        let variants = vec![(120u32, -3i32)];
        let expected = compute_expected_position(&rec, &variants);
        assert_eq!(expected, 100); // No shift
    }

    #[test]
    fn test_compute_expected_position_downstream_variants() {
        // Read at pos=100 with 50M (ends at 149), downstream variant at 200
        let rec = create_test_record(100, "50M");
        let variants = vec![(200u32, 10i32)]; // Far downstream insertion
        let expected = compute_expected_position(&rec, &variants);
        assert_eq!(expected, 100); // No shift
    }

    #[test]
    fn test_compute_expected_position_multiple_upstream() {
        // Read at pos=100, multiple upstream variants
        let rec = create_test_record(100, "50M");
        let variants = vec![
            (30u32, 5i32),  // +5bp insertion
            (50u32, -2i32), // -2bp deletion
            (70u32, 3i32),  // +3bp insertion
        ];
        let expected = compute_expected_position(&rec, &variants);
        // Net shift: +5 - 2 + 3 = +6
        assert_eq!(expected, 106);
    }

    #[test]
    fn test_compute_expected_position_mixed_locations() {
        // Read at pos=100, variants at different locations
        let rec = create_test_record(100, "50M");
        let variants = vec![
            (30u32, 5i32),   // Upstream insertion: +5
            (120u32, 10i32), // Within-read: no shift
            (200u32, -3i32), // Downstream: no shift
        ];
        let expected = compute_expected_position(&rec, &variants);
        // Only upstream counts: +5
        assert_eq!(expected, 105);
    }

    #[test]
    fn test_compute_expected_position_deletion_spanning_start() {
        // Read at pos=100, deletion from 95-105 spans read start
        let rec = create_test_record(100, "50M");
        let variants = vec![(95u32, -10i32)]; // 10bp deletion spanning 95-104
        let expected = compute_expected_position(&rec, &variants);
        // Spanning deletion still shifts (it started upstream)
        assert_eq!(expected, 90);
    }

    #[test]
    fn test_compute_expected_position_insertion_at_boundary() {
        // Read at pos=100, insertion right before read start (at pos=99)
        let rec = create_test_record(100, "50M");
        let variants = vec![(99u32, 5i32)]; // 5bp insertion at 99
        let expected = compute_expected_position(&rec, &variants);
        // Insertion before read start shifts position
        assert_eq!(expected, 105);
    }

    #[test]
    fn test_compute_expected_position_cigar_with_deletion() {
        // Read at pos=100 with deletion in CIGAR: 20M5D30M
        // This covers ref 100-154 (20 + 5 + 30 - 1 = 54 bases)
        let rec = create_test_record(100, "20M5D30M");

        // Upstream variant should still work
        let variants = vec![(50u32, 3i32)];
        let expected = compute_expected_position(&rec, &variants);
        assert_eq!(expected, 103);

        // Within-read variant (in CIGAR deletion region)
        let variants = vec![(120u32, 5i32)]; // pos 120 is in CIGAR deletion
        let expected = compute_expected_position(&rec, &variants);
        assert_eq!(expected, 100); // No shift - within read's ref span
    }

    #[test]
    fn test_compute_expected_position_cigar_with_softclip() {
        // Read at pos=100 with soft clip: 5S45M
        // Soft clip doesn't affect reference span
        let rec = create_test_record(100, "5S45M");

        // Upstream variant
        let variants = vec![(50u32, 5i32)];
        let expected = compute_expected_position(&rec, &variants);
        assert_eq!(expected, 105);

        // Within-read variant
        let variants = vec![(110u32, 5i32)];
        let expected = compute_expected_position(&rec, &variants);
        assert_eq!(expected, 100); // No shift
    }

    #[test]
    fn test_compute_expected_position_large_indels() {
        // Test with larger indels (50bp)
        let rec = create_test_record(1000, "100M");

        // Large upstream insertion
        let variants = vec![(500u32, 50i32)];
        let expected = compute_expected_position(&rec, &variants);
        assert_eq!(expected, 1050);

        // Large upstream deletion
        let variants = vec![(500u32, -50i32)];
        let expected = compute_expected_position(&rec, &variants);
        assert_eq!(expected, 950);
    }

    #[test]
    fn test_compute_expected_position_cigar_aware_full_api() {
        // Test the full API with (start, end, delta) tuples
        let rec = create_test_record(100, "50M");

        // Upstream insertion
        let variants = vec![(50u32, 51u32, 5i32)];
        let expected = compute_expected_position_cigar_aware(&rec, &variants);
        assert_eq!(expected, 105);

        // Within-read deletion
        let variants = vec![(110u32, 115u32, -5i32)];
        let expected = compute_expected_position_cigar_aware(&rec, &variants);
        assert_eq!(expected, 100); // No shift
    }
}
