#![allow(non_local_definitions)]

use pyo3::prelude::*;
use pyo3::exceptions::PyRuntimeError;

// Modules
mod bam_counter;
mod bam_remapper;
mod bam_intersect;
mod bam_filter;  // Fast BAM filtering by variant overlap (replaces samtools process_bam)
mod cigar_utils;  // Shared CIGAR-aware position mapping utilities
mod read_pairer;
mod analysis;
mod mapping_filter;
mod vcf_to_bed;
mod multi_sample;
mod unified_pipeline;  // Single-pass unified make-reads (5x faster)

use bam_counter::BamCounter;
use mapping_filter::filter_bam_wasp;

// ============================================================================
// PyO3 Bindings for BAM Remapping
// ============================================================================

/// Parse intersection BED file (Rust implementation)
///
/// Fast streaming parser that replaces Python's `make_intersect_df()`.
/// Expected speedup: 3.7-6.1x over Polars implementation.
///
/// # Arguments
/// * `intersect_bed` - Path to bedtools intersect output
///
/// # Returns
/// Dictionary mapping read names (bytes) to list of variant spans
///
/// # Example (Python)
/// ```python
/// import wasp2_rust
/// variants = wasp2_rust.parse_intersect_bed("intersect.bed")
/// print(f"Parsed {len(variants)} reads")
/// ```
#[pyfunction]
fn parse_intersect_bed(py: Python, intersect_bed: &str) -> PyResult<PyObject> {
    use pyo3::types::{PyDict, PyList};

    // Call Rust parser
    let variants = bam_remapper::parse_intersect_bed(intersect_bed)
        .map_err(|e| PyRuntimeError::new_err(format!("Failed to parse BED: {}", e)))?;

    // Convert to Python dict
    let py_dict = PyDict::new(py);

    for (read_name, spans) in variants.iter() {
        let py_list = PyList::empty(py);

        for span in spans {
            let span_dict = PyDict::new(py);
            span_dict.set_item("chrom", &span.chrom)?;
            span_dict.set_item("start", span.start)?;
            span_dict.set_item("stop", span.stop)?;
            span_dict.set_item("vcf_start", span.vcf_start)?;
            span_dict.set_item("vcf_stop", span.vcf_stop)?;
            span_dict.set_item("mate", span.mate)?;
            span_dict.set_item("hap1", &span.hap1)?;
            span_dict.set_item("hap2", &span.hap2)?;
            py_list.append(span_dict)?;
        }

        py_dict.set_item(pyo3::types::PyBytes::new(py, read_name), py_list)?;
    }

    Ok(py_dict.into())
}

/// Remap reads for a single chromosome (Rust implementation)
///
/// Replaces Python's `swap_chrom_alleles()` function.
///
/// # Arguments
/// * `bam_path` - Path to BAM file with reads to remap
/// * `intersect_bed` - Path to bedtools intersect output
/// * `chrom` - Chromosome to process (e.g., "chr10")
/// * `out_r1` - Output path for read 1 FASTQ
/// * `out_r2` - Output path for read 2 FASTQ
///
/// # Returns
/// (pairs_processed, haplotypes_generated)
///
/// # Example (Python)
/// ```python
/// import wasp2_rust
/// pairs, haps = wasp2_rust.remap_chromosome(
///     "input.bam",
///     "intersect.bed",
///     "chr10",
///     "remap_r1.fq",
///     "remap_r2.fq"
/// )
/// print(f"Processed {pairs} pairs, generated {haps} haplotypes")
/// ```
#[pyfunction]
#[pyo3(signature = (bam_path, intersect_bed, chrom, out_r1, out_r2, max_seqs=64))]
fn remap_chromosome(
    bam_path: &str,
    intersect_bed: &str,
    chrom: &str,
    out_r1: &str,
    out_r2: &str,
    max_seqs: usize,
) -> PyResult<(usize, usize)> {
    let config = bam_remapper::RemapConfig {
        max_seqs,
        is_phased: true,
    };

    // Parse intersection file
    let variants = bam_remapper::parse_intersect_bed(intersect_bed)
        .map_err(|e| PyRuntimeError::new_err(format!("Failed to parse BED: {}", e)))?;

    // Process chromosome
    let (haplotypes, stats) = bam_remapper::swap_alleles_for_chrom(
        bam_path, &variants, chrom, &config
    ).map_err(|e| PyRuntimeError::new_err(format!("Failed to swap alleles: {}", e)))?;

    // Write FASTQ files
    let (_r1_count, _r2_count) = bam_remapper::write_fastq_pair(
        &haplotypes, out_r1, out_r2
    ).map_err(|e| PyRuntimeError::new_err(format!("Failed to write FASTQ: {}", e)))?;

    Ok((stats.pairs_processed, stats.haplotypes_generated))
}

/// Remap all chromosomes in parallel (Rust implementation)
///
/// High-performance parallel processing of all chromosomes with streaming FASTQ writes.
/// Uses crossbeam channels for producer-consumer pattern - writes happen as processing continues.
///
/// # Arguments
/// * `bam_path` - Path to BAM file
/// * `intersect_bed` - Path to bedtools intersect output
/// * `out_r1` - Output path for read 1 FASTQ
/// * `out_r2` - Output path for read 2 FASTQ
/// * `max_seqs` - Maximum haplotype sequences per read pair (default 64)
/// * `parallel` - Use parallel processing (default true)
/// * `num_threads` - Number of threads (0 = auto-detect, default 0)
///
/// # Returns
/// (pairs_processed, haplotypes_generated)
#[pyfunction]
#[pyo3(signature = (bam_path, intersect_bed, out_r1, out_r2, max_seqs=64, parallel=true, num_threads=0))]
fn remap_all_chromosomes(
    bam_path: &str,
    intersect_bed: &str,
    out_r1: &str,
    out_r2: &str,
    max_seqs: usize,
    parallel: bool,
    num_threads: usize,
) -> PyResult<(usize, usize)> {
    let config = bam_remapper::RemapConfig {
        max_seqs,
        is_phased: true,
    };

    // Parse intersect file ONCE, grouped by chromosome
    // This is the key optimization: 22x fewer parse operations for RNA-seq
    let variants_by_chrom = bam_remapper::parse_intersect_bed_by_chrom(intersect_bed)
        .map_err(|e| PyRuntimeError::new_err(format!("Failed to parse intersect BED: {}", e)))?;

    // Report chromosome count
    let num_chroms = variants_by_chrom.len();
    let total_reads: usize = variants_by_chrom.values().map(|v| v.len()).sum();
    eprintln!("Parsed {} chromosomes with {} reads from intersect file", num_chroms, total_reads);

    let stats = if parallel {
        // Use streaming parallel version with crossbeam channels
        let effective_threads = if num_threads > 0 { num_threads } else { rayon::current_num_threads() };
        eprintln!("Processing {} chromosomes in parallel ({} threads) with streaming writes...", num_chroms, effective_threads);

        bam_remapper::process_and_write_parallel(
            bam_path,
            &variants_by_chrom,
            &config,
            out_r1,
            out_r2,
            num_threads
        ).map_err(|e| PyRuntimeError::new_err(format!("Failed to process chromosomes: {}", e)))?
    } else {
        eprintln!("Processing {} chromosomes sequentially...", num_chroms);
        let (haplotypes, stats) = bam_remapper::process_all_chromosomes_sequential(bam_path, &variants_by_chrom, &config)
            .map_err(|e| PyRuntimeError::new_err(format!("Failed to process chromosomes: {}", e)))?;

        // Write FASTQ output files
        bam_remapper::write_fastq_pair(&haplotypes, out_r1, out_r2)
            .map_err(|e| PyRuntimeError::new_err(format!("Failed to write FASTQ: {}", e)))?;

        stats
    };

    eprintln!(
        "✅ Processed {} pairs → {} haplotypes",
        stats.pairs_processed, stats.haplotypes_generated
    );

    Ok((stats.pairs_processed, stats.haplotypes_generated))
}

// ============================================================================
// PyO3 Bindings for Analysis
// ============================================================================

/// Analyze allelic imbalance (Rust implementation)
///
/// Replaces Python's `get_imbalance()` function from as_analysis.py.
///
/// # Arguments
/// * `tsv_path` - Path to TSV file with allele counts
/// * `min_count` - Minimum total count threshold
/// * `pseudocount` - Pseudocount to add to allele counts
/// * `method` - Analysis method ("single" or "linear")
///
/// # Returns
/// List of dictionaries with imbalance results
///
/// # Example (Python)
/// ```python
/// import wasp2_rust
/// results = wasp2_rust.analyze_imbalance(
///     "counts.tsv",
///     min_count=10,
///     pseudocount=1,
///     method="single"
/// )
/// for r in results:
///     print(f"{r['region']}: pval={r['pval']:.4f}")
/// ```
#[pyfunction]
#[pyo3(signature = (tsv_path, min_count=10, pseudocount=1, method="single"))]
fn analyze_imbalance(
    py: Python,
    tsv_path: &str,
    min_count: u32,
    pseudocount: u32,
    method: &str,
) -> PyResult<PyObject> {
    use pyo3::types::{PyDict, PyList};
    use std::fs::File;
    use std::io::{BufRead, BufReader};

    // Parse method
    let analysis_method = match method {
        "single" => analysis::AnalysisMethod::Single,
        "linear" => analysis::AnalysisMethod::Linear,
        _ => return Err(PyRuntimeError::new_err(format!("Unknown method: {}", method))),
    };

    let config = analysis::AnalysisConfig {
        min_count,
        pseudocount,
        method: analysis_method,
    };

    // Read TSV file
    let file = File::open(tsv_path)
        .map_err(|e| PyRuntimeError::new_err(format!("Failed to open TSV: {}", e)))?;
    let reader = BufReader::new(file);

    let mut variants = Vec::new();
    let mut header_seen = false;

    for line in reader.lines() {
        let line = line.map_err(|e| PyRuntimeError::new_err(format!("Failed to read line: {}", e)))?;

        if !header_seen {
            header_seen = true;
            continue; // Skip header
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 7 {
            continue;
        }

        // Parse fields: chrom, pos, ref, alt, region, ref_count, alt_count, other_count
        let chrom = fields[0].to_string();
        let pos = fields[1].parse::<u32>()
            .map_err(|e| PyRuntimeError::new_err(format!("Invalid pos: {}", e)))?;
        let ref_count = fields[5].parse::<u32>()
            .map_err(|e| PyRuntimeError::new_err(format!("Invalid ref_count: {}", e)))?;
        let alt_count = fields[6].parse::<u32>()
            .map_err(|e| PyRuntimeError::new_err(format!("Invalid alt_count: {}", e)))?;

        // Create region identifier (chrom_pos_pos+1 format to match Python)
        let region = format!("{}_{}_{}",  chrom, pos, pos + 1);

        variants.push(analysis::VariantCounts {
            chrom,
            pos,
            ref_count,
            alt_count,
            region,
        });
    }

    // Run analysis
    let results = analysis::analyze_imbalance(variants, &config)
        .map_err(|e| PyRuntimeError::new_err(format!("Analysis failed: {}", e)))?;

    // Convert to Python list of dicts
    let py_list = PyList::empty(py);

    for result in results {
        let py_dict = PyDict::new(py);
        py_dict.set_item("region", result.region)?;
        py_dict.set_item("ref_count", result.ref_count)?;
        py_dict.set_item("alt_count", result.alt_count)?;
        py_dict.set_item("N", result.n)?;
        py_dict.set_item("snp_count", result.snp_count)?;
        py_dict.set_item("null_ll", result.null_ll)?;
        py_dict.set_item("alt_ll", result.alt_ll)?;
        py_dict.set_item("mu", result.mu)?;
        py_dict.set_item("lrt", result.lrt)?;
        py_dict.set_item("pval", result.pval)?;
        py_dict.set_item("fdr_pval", result.fdr_pval)?;
        py_list.append(py_dict)?;
    }

    Ok(py_list.into())
}

// ============================================================================
// PyO3 Bindings for BAM-BED Intersection (coitrees)
// ============================================================================

/// Intersect BAM reads with variant BED file (Rust/coitrees implementation)
///
/// Replaces pybedtools intersect with 15-30x faster Rust implementation
/// using coitrees van Emde Boas layout interval trees.
///
/// # Arguments
/// * `bam_path` - Path to sorted BAM file
/// * `bed_path` - Path to variant BED file (chrom, start, stop, ref, alt, GT)
/// * `out_path` - Output path for intersections
///
/// # Returns
/// Number of intersections found
///
/// # Example (Python)
/// ```python
/// import wasp2_rust
/// count = wasp2_rust.intersect_bam_bed("reads.bam", "variants.bed", "out.bed")
/// print(f"Found {count} read-variant overlaps")
/// ```
#[pyfunction]
fn intersect_bam_bed(bam_path: &str, bed_path: &str, out_path: &str) -> PyResult<usize> {
    bam_intersect::intersect_bam_with_variants(bam_path, bed_path, out_path)
        .map_err(|e| PyRuntimeError::new_err(format!("Intersect failed: {}", e)))
}

/// Intersect BAM reads with multi-sample variant BED file
///
/// # Arguments
/// * `bam_path` - Path to sorted BAM file
/// * `bed_path` - Path to variant BED file with multiple GT columns
/// * `out_path` - Output path for intersections
/// * `num_samples` - Number of sample genotype columns in BED
///
/// # Returns
/// Number of intersections found
#[pyfunction]
fn intersect_bam_bed_multi(
    bam_path: &str,
    bed_path: &str,
    out_path: &str,
    num_samples: usize,
) -> PyResult<usize> {
    bam_intersect::intersect_bam_with_variants_multi(bam_path, bed_path, out_path, num_samples)
        .map_err(|e| PyRuntimeError::new_err(format!("Multi-sample intersect failed: {}", e)))
}

// ============================================================================
// PyO3 Bindings for BAM Filtering (replaces samtools process_bam)
// ============================================================================

/// Filter BAM by variant overlap (Rust implementation)
///
/// Replaces Python's process_bam() function which uses samtools subprocess calls.
/// Expected speedup: 4-5x (from ~450s to ~100s for 56M reads).
///
/// # Algorithm
/// 1. Build coitrees interval tree from variant BED file
/// 2. Stream BAM, collect read names overlapping variants
/// 3. Stream BAM again, split to remap/keep based on name membership
///
/// # Arguments
/// * `bam_path` - Input BAM file (should be coordinate-sorted)
/// * `bed_path` - Variant BED file (chrom, start, stop, ref, alt, GT)
/// * `remap_bam_path` - Output BAM for reads needing remapping
/// * `keep_bam_path` - Output BAM for reads not needing remapping
/// * `is_paired` - Whether reads are paired-end (default: true)
/// * `threads` - Number of threads to use (default: 4)
///
/// # Returns
/// Tuple of (remap_count, keep_count, unique_names)
///
/// # Example (Python)
/// ```python
/// import wasp2_rust
/// remap, keep, names = wasp2_rust.filter_bam_by_variants(
///     "input.bam",
///     "variants.bed",
///     "remap.bam",
///     "keep.bam",
///     is_paired=True,
///     threads=4
/// )
/// print(f"Split: {remap} remap, {keep} keep ({names} unique names)")
/// ```
#[pyfunction]
#[pyo3(signature = (bam_path, bed_path, remap_bam_path, keep_bam_path, is_paired=true, threads=4))]
fn filter_bam_by_variants_py(
    bam_path: &str,
    bed_path: &str,
    remap_bam_path: &str,
    keep_bam_path: &str,
    is_paired: bool,
    threads: usize,
) -> PyResult<(usize, usize, usize)> {
    let stats = bam_filter::filter_bam_by_variants(
        bam_path,
        bed_path,
        remap_bam_path,
        keep_bam_path,
        is_paired,
        threads,
    )
    .map_err(|e| PyRuntimeError::new_err(format!("BAM filter failed: {}", e)))?;

    Ok((stats.remap_reads, stats.keep_reads, stats.unique_remap_names))
}

// ============================================================================
// PyO3 Bindings for Unified Pipeline (Single-pass make-reads)
// ============================================================================

/// Unified single-pass make-reads pipeline (Rust implementation)
///
/// Replaces the multi-step Python pipeline (filter + intersect + remap) with a
/// single-pass Rust implementation that streams directly from BAM to FASTQ.
/// Expected speedup: 5x (from ~500s to ~100s for 56M reads).
///
/// # Algorithm
/// 1. Build coitrees interval tree from variant BED file
/// 2. Stream BAM ONCE, buffer pairs, check variant overlap
/// 3. For overlapping pairs: generate haplotypes, write to FASTQ
/// 4. Track stats: pairs processed, haplotypes generated
///
/// # Arguments
/// * `bam_path` - Input BAM file (should be coordinate-sorted)
/// * `bed_path` - Variant BED file (chrom, start, stop, ref, alt, GT)
/// * `out_r1` - Output path for read 1 FASTQ
/// * `out_r2` - Output path for read 2 FASTQ
/// * `max_seqs` - Maximum haplotype sequences per read pair (default: 64)
/// * `threads` - Number of threads to use (default: 8)
/// * `channel_buffer` - Channel buffer size for streaming (default: 50000)
///
/// # Returns
/// Dictionary with stats: pairs_processed, pairs_with_variants, haplotypes_written, etc.
///
/// # Example (Python)
/// ```python
/// import wasp2_rust
/// stats = wasp2_rust.unified_make_reads(
///     "input.bam",
///     "variants.bed",
///     "remap_r1.fq",
///     "remap_r2.fq",
///     max_seqs=64,
///     threads=8
/// )
/// print(f"Processed {stats['pairs_processed']} pairs -> {stats['haplotypes_written']} haplotypes")
/// ```
#[pyfunction]
#[pyo3(signature = (bam_path, bed_path, out_r1, out_r2, max_seqs=64, threads=8, channel_buffer=50000, compression_threads=4))]
fn unified_make_reads_py(
    py: Python,
    bam_path: &str,
    bed_path: &str,
    out_r1: &str,
    out_r2: &str,
    max_seqs: usize,
    threads: usize,
    channel_buffer: usize,
    compression_threads: usize,
) -> PyResult<PyObject> {
    use pyo3::types::PyDict;

    let config = unified_pipeline::UnifiedConfig {
        read_threads: threads,
        max_seqs,
        channel_buffer,
        compression_threads,
    };

    let stats = unified_pipeline::unified_make_reads(
        bam_path,
        bed_path,
        out_r1,
        out_r2,
        &config,
    )
    .map_err(|e| PyRuntimeError::new_err(format!("Unified pipeline failed: {}", e)))?;

    // Return stats as Python dict
    let py_dict = PyDict::new(py);
    py_dict.set_item("total_reads", stats.total_reads)?;
    py_dict.set_item("pairs_processed", stats.pairs_processed)?;
    py_dict.set_item("pairs_with_variants", stats.pairs_with_variants)?;
    py_dict.set_item("haplotypes_written", stats.haplotypes_written)?;
    py_dict.set_item("pairs_kept", stats.pairs_kept)?;
    py_dict.set_item("pairs_skipped_unmappable", stats.pairs_skipped_unmappable)?;
    py_dict.set_item("pairs_haplotype_failed", stats.pairs_haplotype_failed)?;
    py_dict.set_item("orphan_reads", stats.orphan_reads)?;
    py_dict.set_item("tree_build_ms", stats.tree_build_ms)?;
    py_dict.set_item("bam_stream_ms", stats.bam_stream_ms)?;

    Ok(py_dict.into())
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
/// - Only HaplotypeOutput (Vec<u8>) is sent via channel
///
/// # Arguments
/// * `bam_path` - Input BAM file (must be coordinate-sorted and indexed)
/// * `bed_path` - Variant BED file (chrom, start, stop, ref, alt, GT)
/// * `out_r1` - Output path for read 1 FASTQ
/// * `out_r2` - Output path for read 2 FASTQ
/// * `max_seqs` - Maximum haplotype sequences per read pair (default: 64)
/// * `threads` - Number of threads to use (default: 8)
/// * `channel_buffer` - Channel buffer size for streaming (default: 50000)
/// * `compression_threads` - Threads per FASTQ file for gzip (default: 4)
///
/// # Returns
/// Dictionary with stats: pairs_processed, pairs_with_variants, haplotypes_written, etc.
///
/// # Example (Python)
/// ```python
/// import wasp2_rust
/// stats = wasp2_rust.unified_make_reads_parallel(
///     "input.bam",  # Must have .bai index
///     "variants.bed",
///     "remap_r1.fq.gz",
///     "remap_r2.fq.gz",
///     max_seqs=64,
///     threads=8
/// )
/// print(f"Processed {stats['pairs_processed']} pairs -> {stats['haplotypes_written']} haplotypes")
/// ```
#[pyfunction]
#[pyo3(signature = (bam_path, bed_path, out_r1, out_r2, max_seqs=64, threads=8, channel_buffer=50000, compression_threads=4))]
fn unified_make_reads_parallel_py(
    py: Python,
    bam_path: &str,
    bed_path: &str,
    out_r1: &str,
    out_r2: &str,
    max_seqs: usize,
    threads: usize,
    channel_buffer: usize,
    compression_threads: usize,
) -> PyResult<PyObject> {
    use pyo3::types::PyDict;

    // Configure rayon thread pool
    if threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .ok(); // Ignore if already set
    }

    let config = unified_pipeline::UnifiedConfig {
        read_threads: threads,
        max_seqs,
        channel_buffer,
        compression_threads,
    };

    let stats = unified_pipeline::unified_make_reads_parallel(
        bam_path,
        bed_path,
        out_r1,
        out_r2,
        &config,
    )
    .map_err(|e| PyRuntimeError::new_err(format!("Parallel unified pipeline failed: {}", e)))?;

    // Return stats as Python dict
    let py_dict = PyDict::new(py);
    py_dict.set_item("total_reads", stats.total_reads)?;
    py_dict.set_item("pairs_processed", stats.pairs_processed)?;
    py_dict.set_item("pairs_with_variants", stats.pairs_with_variants)?;
    py_dict.set_item("haplotypes_written", stats.haplotypes_written)?;
    py_dict.set_item("pairs_kept", stats.pairs_kept)?;
    py_dict.set_item("pairs_skipped_unmappable", stats.pairs_skipped_unmappable)?;
    py_dict.set_item("pairs_haplotype_failed", stats.pairs_haplotype_failed)?;
    py_dict.set_item("orphan_reads", stats.orphan_reads)?;
    py_dict.set_item("tree_build_ms", stats.tree_build_ms)?;
    py_dict.set_item("bam_stream_ms", stats.bam_stream_ms)?;

    Ok(py_dict.into())
}

// ============================================================================
// PyO3 Bindings for VCF/BCF to BED Conversion
// ============================================================================

/// Convert VCF/BCF to BED format (Rust/noodles implementation)
///
/// Replaces bcftools subprocess with 5-6x faster pure Rust implementation.
/// Supports VCF, VCF.gz, and BCF formats.
///
/// # Arguments
/// * `vcf_path` - Path to VCF/BCF file
/// * `bed_path` - Output BED file path
/// * `samples` - Optional list of sample names to extract (None = all)
/// * `het_only` - Only output heterozygous sites (default: true)
/// * `include_indels` - Include indels, not just SNPs (default: false)
/// * `max_indel_len` - Maximum indel length to include (default: 10)
/// * `include_genotypes` - Include genotype column in output (default: true)
///
/// # Returns
/// Number of variants written to BED file
///
/// # Example (Python)
/// ```python
/// import wasp2_rust
/// count = wasp2_rust.vcf_to_bed(
///     "variants.vcf.gz",
///     "variants.bed",
///     samples=["NA12878"],
///     het_only=True
/// )
/// print(f"Wrote {count} het variants")
/// ```
#[pyfunction]
#[pyo3(signature = (vcf_path, bed_path, samples=None, het_only=true, include_indels=false, max_indel_len=10, include_genotypes=true))]
fn vcf_to_bed_py(
    vcf_path: &str,
    bed_path: &str,
    samples: Option<Vec<String>>,
    het_only: bool,
    include_indels: bool,
    max_indel_len: usize,
    include_genotypes: bool,
) -> PyResult<usize> {
    let config = vcf_to_bed::VcfToBedConfig {
        samples,
        het_only,
        include_indels,
        max_indel_len,
        include_genotypes,
    };

    vcf_to_bed::vcf_to_bed(vcf_path, bed_path, &config)
        .map_err(|e| PyRuntimeError::new_err(format!("VCF to BED failed: {}", e)))
}

// ============================================================================
// PyO3 Bindings for Multi-Sample Processing
// ============================================================================

/// Parse multi-sample intersection BED file (Rust implementation)
///
/// Parses BED file with multiple sample genotype columns.
/// Used for multi-sample WASP2 processing.
///
/// # Arguments
/// * `intersect_bed` - Path to intersection BED file
/// * `num_samples` - Number of sample genotype columns
///
/// # Returns
/// Dictionary mapping read names to variant spans with all sample genotypes
///
/// # Example (Python)
/// ```python
/// import wasp2_rust
/// variants = wasp2_rust.parse_intersect_bed_multi("intersect.bed", num_samples=3)
/// ```
#[pyfunction]
fn parse_intersect_bed_multi(
    py: Python,
    intersect_bed: &str,
    num_samples: usize,
) -> PyResult<PyObject> {
    use pyo3::types::{PyDict, PyList};

    let variants = multi_sample::parse_intersect_bed_multi(intersect_bed, num_samples)
        .map_err(|e| PyRuntimeError::new_err(format!("Failed to parse multi-sample BED: {}", e)))?;

    // Convert to Python dict
    let py_dict = PyDict::new(py);

    for (read_name, spans) in variants.iter() {
        let py_list = PyList::empty(py);

        for span in spans {
            let span_dict = PyDict::new(py);
            span_dict.set_item("chrom", &span.chrom)?;
            span_dict.set_item("start", span.start)?;
            span_dict.set_item("stop", span.stop)?;
            span_dict.set_item("vcf_start", span.vcf_start)?;
            span_dict.set_item("vcf_stop", span.vcf_stop)?;
            span_dict.set_item("mate", span.mate)?;
            span_dict.set_item("ref_allele", &span.ref_allele)?;
            span_dict.set_item("alt_allele", &span.alt_allele)?;

            // Convert sample_alleles to list of tuples
            let alleles_list = PyList::empty(py);
            for (h1, h2) in &span.sample_alleles {
                let tuple = pyo3::types::PyTuple::new(py, &[h1.as_str(), h2.as_str()]);
                alleles_list.append(tuple)?;
            }
            span_dict.set_item("sample_alleles", alleles_list)?;

            py_list.append(span_dict)?;
        }

        py_dict.set_item(pyo3::types::PyBytes::new(py, read_name), py_list)?;
    }

    Ok(py_dict.into())
}

/// Remap reads for a single chromosome - multi-sample version (Rust implementation)
///
/// Replaces Python's `swap_chrom_alleles_multi()` function.
/// Generates unique haplotype sequences across all samples.
///
/// # Arguments
/// * `bam_path` - Path to BAM file with reads to remap
/// * `intersect_bed` - Path to bedtools intersect output (multi-sample format)
/// * `chrom` - Chromosome to process (e.g., "chr10")
/// * `out_r1` - Output path for read 1 FASTQ
/// * `out_r2` - Output path for read 2 FASTQ
/// * `num_samples` - Number of samples in the intersection BED
/// * `max_seqs` - Maximum haplotype sequences per read pair (default 64)
///
/// # Returns
/// (pairs_processed, haplotypes_generated)
///
/// # Example (Python)
/// ```python
/// import wasp2_rust
/// pairs, haps = wasp2_rust.remap_chromosome_multi(
///     "input.bam",
///     "intersect.bed",
///     "chr10",
///     "remap_r1.fq",
///     "remap_r2.fq",
///     num_samples=3,
///     max_seqs=64
/// )
/// print(f"Processed {pairs} pairs, generated {haps} haplotypes")
/// ```
#[pyfunction]
#[pyo3(signature = (bam_path, intersect_bed, chrom, out_r1, out_r2, num_samples, max_seqs=64))]
fn remap_chromosome_multi(
    bam_path: &str,
    intersect_bed: &str,
    chrom: &str,
    out_r1: &str,
    out_r2: &str,
    num_samples: usize,
    max_seqs: usize,
) -> PyResult<(usize, usize)> {
    // Parse multi-sample intersection file
    let variants = multi_sample::parse_intersect_bed_multi(intersect_bed, num_samples)
        .map_err(|e| PyRuntimeError::new_err(format!("Failed to parse multi-sample BED: {}", e)))?;

    // Process chromosome
    let stats = multi_sample::swap_alleles_for_chrom_multi(
        bam_path, &variants, chrom, out_r1, out_r2, max_seqs
    ).map_err(|e| PyRuntimeError::new_err(format!("Failed to swap alleles: {}", e)))?;

    Ok((stats.pairs_processed, stats.haplotypes_generated))
}

// ============================================================================
// Legacy Functions (keep for compatibility)
// ============================================================================

/// Simple test function to verify PyO3 is working
#[pyfunction]
fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
    Ok((a + b).to_string())
}

// ============================================================================
// Module Definition
// ============================================================================

/// WASP2 Rust acceleration module
///
/// Provides high-performance implementations of bottleneck functions:
/// - BamCounter: Fast allele counting (IMPLEMENTED)
/// - intersect_bam_bed: Fast BAM-BED intersection using coitrees (41x faster)
/// - intersect_bam_bed_multi: Multi-sample BAM-BED intersection (41x faster)
/// - vcf_to_bed: Fast VCF/BCF to BED conversion using noodles (5-6x faster)
/// - remap_chromosome: Fast allele swapping for mapping stage (IMPLEMENTED)
/// - remap_chromosome_multi: Multi-sample allele swapping (IMPLEMENTED)
/// - remap_all_chromosomes: Parallel processing of all chromosomes (skeleton)
/// - parse_intersect_bed_multi: Multi-sample intersection parsing (IMPLEMENTED)
/// - analyze_imbalance: Fast beta-binomial analysis for AI detection (IMPLEMENTED)
#[pymodule]
fn wasp2_rust(_py: Python, m: &PyModule) -> PyResult<()> {
    // Legacy test function
    m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;

    // Counting module (IMPLEMENTED)
    m.add_class::<BamCounter>()?;

    // BAM-BED intersection using coitrees (41x faster than pybedtools)
    m.add_function(wrap_pyfunction!(intersect_bam_bed, m)?)?;
    m.add_function(wrap_pyfunction!(intersect_bam_bed_multi, m)?)?;

    // VCF/BCF to BED conversion using noodles (5-6x faster than bcftools)
    m.add_function(wrap_pyfunction!(vcf_to_bed_py, m)?)?;

    // Remapping module - parser (IMPLEMENTED)
    m.add_function(wrap_pyfunction!(parse_intersect_bed, m)?)?;

    // Multi-sample intersection parsing (NEW)
    m.add_function(wrap_pyfunction!(parse_intersect_bed_multi, m)?)?;

    // Remapping module - full pipeline (IMPLEMENTED)
    m.add_function(wrap_pyfunction!(remap_chromosome, m)?)?;
    m.add_function(wrap_pyfunction!(remap_chromosome_multi, m)?)?;
    m.add_function(wrap_pyfunction!(remap_all_chromosomes, m)?)?;

    // Mapping filter (WASP remap filter)
    m.add_function(wrap_pyfunction!(filter_bam_wasp, m)?)?;

    // BAM filtering by variant overlap (replaces samtools process_bam, 4-5x faster)
    m.add_function(wrap_pyfunction!(filter_bam_by_variants_py, m)?)?;

    // Unified single-pass pipeline (replaces filter + intersect + remap, 5x faster)
    m.add_function(wrap_pyfunction!(unified_make_reads_py, m)?)?;

    // Parallel unified pipeline (3-8x speedup over sequential, requires BAM index)
    m.add_function(wrap_pyfunction!(unified_make_reads_parallel_py, m)?)?;

    // Analysis module (beta-binomial allelic imbalance detection)
    m.add_function(wrap_pyfunction!(analyze_imbalance, m)?)?;

    Ok(())
}
