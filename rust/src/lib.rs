#![allow(non_local_definitions)]

use pyo3::prelude::*;
use pyo3::exceptions::PyRuntimeError;

// Modules
mod bam_counter;
mod bam_remapper;
mod read_pairer;

use bam_counter::BamCounter;

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
/// High-performance parallel processing of all chromosomes.
///
/// # Arguments
/// * `bam_path` - Path to BAM file
/// * `intersect_bed` - Path to bedtools intersect output
/// * `out_r1` - Output path for read 1 FASTQ
/// * `out_r2` - Output path for read 2 FASTQ
///
/// # Returns
/// (pairs_processed, haplotypes_generated)
#[pyfunction]
#[pyo3(signature = (bam_path, intersect_bed, out_r1, out_r2, max_seqs=64))]
fn remap_all_chromosomes(
    bam_path: &str,
    intersect_bed: &str,
    out_r1: &str,
    out_r2: &str,
    max_seqs: usize,
) -> PyResult<(usize, usize)> {
    // TODO: Implement when bam_remapper functions are ready
    let _ = (bam_path, intersect_bed, out_r1, out_r2, max_seqs);

    Err(PyRuntimeError::new_err(
        "remap_all_chromosomes not yet implemented - skeleton only"
    ))

    // Future implementation with rayon parallelism:
    // let config = bam_remapper::RemapConfig { max_seqs, is_phased: true };
    // let variants = bam_remapper::parse_intersect_bed(intersect_bed)
    //     .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    // let (haplotypes, stats) = bam_remapper::process_all_chromosomes_parallel(
    //     bam_path, &variants, &config
    // ).map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    // let (r1_count, r2_count) = bam_remapper::write_fastq_pair(
    //     &haplotypes, out_r1, out_r2
    // ).map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    // Ok((stats.pairs_processed, stats.haplotypes_generated))
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
/// - BamCounter: Fast allele counting (already implemented)
/// - remap_chromosome: Fast allele swapping for mapping stage (NEW)
/// - remap_all_chromosomes: Parallel processing of all chromosomes (NEW)
#[pymodule]
fn wasp2_rust(_py: Python, m: &PyModule) -> PyResult<()> {
    // Legacy test function
    m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;

    // Counting module (existing)
    m.add_class::<BamCounter>()?;

    // Remapping module - parser (IMPLEMENTED)
    m.add_function(wrap_pyfunction!(parse_intersect_bed, m)?)?;

    // Remapping module - full pipeline (skeleton)
    m.add_function(wrap_pyfunction!(remap_chromosome, m)?)?;
    m.add_function(wrap_pyfunction!(remap_all_chromosomes, m)?)?;

    Ok(())
}
