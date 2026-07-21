#![allow(non_local_definitions)]

use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use pyo3::types::PyModule;

// Modules
mod analysis;
mod bam_counter;
mod bam_counter_sc; // Per-cell (single-cell) allele counter, mirrors bam_counter
mod bam_filter; // Fast BAM filtering by variant overlap (replaces samtools process_bam)
mod bam_intersect;
mod bam_remapper;
mod cigar_utils; // Shared CIGAR-aware position mapping utilities
mod mapping_filter;
mod multi_sample;
mod read_pairer;
mod seq_decode;
mod unified_pipeline;
mod vcf_to_bed;

pub use unified_pipeline::{
    unified_make_reads, unified_make_reads_parallel, UnifiedConfig, UnifiedStats,
};

use bam_counter::BamCounter;
use bam_counter_sc::BamCounterSC;
use mapping_filter::filter_bam_wasp;

fn parse_numeric_phased_gt(gt: &str) -> Option<u8> {
    match gt {
        "0|1" => Some(0),
        "1|0" => Some(1),
        _ => None,
    }
}

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
fn parse_intersect_bed(py: Python, intersect_bed: &str) -> PyResult<Py<PyAny>> {
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

    Ok(py_dict.unbind().into_any())
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
    let (haplotypes, stats) =
        bam_remapper::swap_alleles_for_chrom(bam_path, &variants, chrom, &config)
            .map_err(|e| PyRuntimeError::new_err(format!("Failed to swap alleles: {}", e)))?;

    // Write FASTQ files
    let (_r1_count, _r2_count) = bam_remapper::write_fastq_pair(&haplotypes, out_r1, out_r2)
        .map_err(|e| PyRuntimeError::new_err(format!("Failed to write FASTQ: {}", e)))?;

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
    eprintln!(
        "Parsed {} chromosomes with {} reads from intersect file",
        num_chroms, total_reads
    );

    let stats = if parallel {
        // Use streaming parallel version with crossbeam channels
        let effective_threads = if num_threads > 0 {
            num_threads
        } else {
            rayon::current_num_threads()
        };
        eprintln!(
            "Processing {} chromosomes in parallel ({} threads) with streaming writes...",
            num_chroms, effective_threads
        );

        bam_remapper::process_and_write_parallel(
            bam_path,
            &variants_by_chrom,
            &config,
            out_r1,
            out_r2,
            num_threads,
        )
        .map_err(|e| PyRuntimeError::new_err(format!("Failed to process chromosomes: {}", e)))?
    } else {
        eprintln!("Processing {} chromosomes sequentially...", num_chroms);
        let (haplotypes, stats) =
            bam_remapper::process_all_chromosomes_sequential(bam_path, &variants_by_chrom, &config)
                .map_err(|e| {
                    PyRuntimeError::new_err(format!("Failed to process chromosomes: {}", e))
                })?;

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
fn parse_analysis_method(method: &str) -> PyResult<analysis::AnalysisMethod> {
    match method {
        "single" => Ok(analysis::AnalysisMethod::Single),
        "linear" => Ok(analysis::AnalysisMethod::Linear),
        "per_donor" => Err(PyRuntimeError::new_err(
            "method='per_donor' has been removed because it pooled effect likelihoods across donors; split the count table by donor, fit each donor with fit_imbalance_dispersion(), and pass those fixed nuisance parameters to analyze_imbalance_run()",
        )),
        _ => Err(PyRuntimeError::new_err(format!(
            "Unknown method '{method}'; expected 'single' or 'linear'"
        ))),
    }
}

fn fixed_dispersion_parameters(
    method: analysis::AnalysisMethod,
    rho: Option<f64>,
    linear_d1: Option<f64>,
    linear_d2: Option<f64>,
    linear_depth_center: Option<f64>,
    linear_depth_scale: Option<f64>,
) -> PyResult<Option<analysis::DispersionParameters>> {
    match method {
        analysis::AnalysisMethod::Single => {
            if linear_d1.is_some()
                || linear_d2.is_some()
                || linear_depth_center.is_some()
                || linear_depth_scale.is_some()
            {
                return Err(PyRuntimeError::new_err(
                    "method='single' accepts rho, not linear dispersion parameters",
                ));
            }
            Ok(rho.map(|rho| analysis::DispersionParameters::Single { rho }))
        }
        analysis::AnalysisMethod::Linear => {
            if rho.is_some() {
                return Err(PyRuntimeError::new_err(
                    "method='linear' accepts linear_d1/linear_d2, not rho",
                ));
            }
            match (
                linear_d1,
                linear_d2,
                linear_depth_center,
                linear_depth_scale,
            ) {
                (None, None, None, None) => Ok(None),
                (Some(d1), Some(d2), Some(depth_center), Some(depth_scale)) => {
                    Ok(Some(analysis::DispersionParameters::Linear {
                        d1,
                        d2,
                        depth_center,
                        depth_scale,
                    }))
                }
                _ => Err(PyRuntimeError::new_err(
                    "fixed method='linear' requires linear_d1, linear_d2, \
                     linear_depth_center, and linear_depth_scale",
                )),
            }
        }
    }
}

fn read_analysis_variants(
    tsv_path: &str,
    region_col: Option<&str>,
) -> PyResult<Vec<analysis::VariantCounts>> {
    use flate2::read::MultiGzDecoder;
    use std::collections::HashSet;
    use std::fs::File;
    use std::io::{BufRead, BufReader};

    let file = File::open(tsv_path)
        .map_err(|e| PyRuntimeError::new_err(format!("Failed to open TSV: {e}")))?;
    let mut reader: Box<dyn BufRead> = if tsv_path.ends_with(".gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut header_line = String::new();
    if reader
        .read_line(&mut header_line)
        .map_err(|e| PyRuntimeError::new_err(format!("Failed to read TSV header: {e}")))?
        == 0
    {
        return Err(PyRuntimeError::new_err("Count TSV is empty"));
    }
    let headers: Vec<&str> = header_line
        .trim_end_matches(['\r', '\n'])
        .split('\t')
        .collect();
    let unique: HashSet<&str> = headers.iter().copied().collect();
    if unique.len() != headers.len() {
        return Err(PyRuntimeError::new_err(
            "Count TSV contains duplicate column names",
        ));
    }
    for name in ["chrom", "pos", "ref", "alt", "ref_count", "alt_count"] {
        if !unique.contains(name) {
            return Err(PyRuntimeError::new_err(format!(
                "Required column '{name}' is missing from count TSV"
            )));
        }
    }
    let index = |name: &str| {
        headers
            .iter()
            .position(|header| *header == name)
            .expect("required header checked")
    };
    let chrom_idx = index("chrom");
    let pos_idx = index("pos");
    let ref_idx = index("ref");
    let alt_idx = index("alt");
    let ref_count_idx = index("ref_count");
    let alt_count_idx = index("alt_count");
    let gt_idx = headers.iter().position(|header| *header == "GT");
    let region_idx = match region_col {
        Some(name) => Some(
            headers
                .iter()
                .position(|header| *header == name)
                .ok_or_else(|| {
                    PyRuntimeError::new_err(format!(
                        "region_col '{name}' not found in header columns: {headers:?}"
                    ))
                })?,
        ),
        None => None,
    };

    let mut variants = Vec::new();
    for (offset, line) in reader.lines().enumerate() {
        let line_number = offset + 2;
        let line = line.map_err(|e| {
            PyRuntimeError::new_err(format!("Failed to read line {line_number}: {e}"))
        })?;
        if line.is_empty() {
            return Err(PyRuntimeError::new_err(format!(
                "Blank row at line {line_number}"
            )));
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != headers.len() {
            return Err(PyRuntimeError::new_err(format!(
                "Line {line_number} has {} fields; expected {}",
                fields.len(),
                headers.len()
            )));
        }

        let chrom = fields[chrom_idx].to_string();
        let ref_allele = fields[ref_idx];
        let alt_allele = fields[alt_idx];
        if chrom.is_empty() || ref_allele.is_empty() || alt_allele.is_empty() {
            return Err(PyRuntimeError::new_err(format!(
                "Empty chromosome, REF, or ALT value at line {line_number}"
            )));
        }
        if ref_allele == alt_allele || alt_allele.contains(',') {
            return Err(PyRuntimeError::new_err(format!(
                "Line {line_number} is not a biallelic REF/ALT observation"
            )));
        }
        let pos = fields[pos_idx].parse::<u32>().map_err(|e| {
            PyRuntimeError::new_err(format!("Invalid pos at line {line_number}: {e}"))
        })?;
        if pos == 0 {
            return Err(PyRuntimeError::new_err(format!(
                "Position must be one-based at line {line_number}"
            )));
        }
        let ref_count = fields[ref_count_idx].parse::<u32>().map_err(|e| {
            PyRuntimeError::new_err(format!("Invalid ref_count at line {line_number}: {e}"))
        })?;
        let alt_count = fields[alt_count_idx].parse::<u32>().map_err(|e| {
            PyRuntimeError::new_err(format!("Invalid alt_count at line {line_number}: {e}"))
        })?;
        let region = region_idx
            .map(|idx| fields[idx].to_string())
            .unwrap_or_else(|| format!("{chrom}_{pos}"));
        if region.is_empty() {
            return Err(PyRuntimeError::new_err(format!(
                "Empty region at line {line_number}"
            )));
        }
        let gt = gt_idx.and_then(|idx| parse_numeric_phased_gt(fields[idx]));
        variants.push(analysis::VariantCounts {
            chrom,
            pos,
            ref_count,
            alt_count,
            region,
            gt,
        });
    }
    Ok(variants)
}

fn results_to_py<'py>(
    py: Python<'py>,
    results: Vec<analysis::ImbalanceResult>,
) -> PyResult<Bound<'py, pyo3::types::PyList>> {
    use pyo3::types::{PyDict, PyList};

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
    Ok(py_list)
}

#[pyfunction]
#[pyo3(signature = (tsv_path, min_count=10, pseudocount=1, method="single", phased=false, region_col=None))]
fn analyze_imbalance(
    py: Python,
    tsv_path: &str,
    min_count: u32,
    pseudocount: u32,
    method: &str,
    phased: bool,
    region_col: Option<String>,
) -> PyResult<Py<PyAny>> {
    let config = analysis::AnalysisConfig {
        min_count,
        pseudocount,
        method: parse_analysis_method(method)?,
        phased,
    };
    let variants = read_analysis_variants(tsv_path, region_col.as_deref())?;
    let output = analysis::analyze_imbalance_detailed(variants, &config)
        .map_err(|e| PyRuntimeError::new_err(format!("Analysis failed: {e}")))?;
    Ok(results_to_py(py, output.results)?.unbind().into_any())
}

/// Fit balanced-null dispersion without testing allelic effects.
#[pyfunction]
#[pyo3(signature = (tsv_path, min_count=10, pseudocount=1, method="single", phased=false, region_col=None))]
fn fit_imbalance_dispersion(
    py: Python,
    tsv_path: &str,
    min_count: u32,
    pseudocount: u32,
    method: &str,
    phased: bool,
    region_col: Option<String>,
) -> PyResult<Py<PyAny>> {
    use pyo3::types::PyDict;

    let analysis_method = parse_analysis_method(method)?;
    let config = analysis::AnalysisConfig {
        min_count,
        pseudocount,
        method: analysis_method,
        phased,
    };
    let variants = read_analysis_variants(tsv_path, region_col.as_deref())?;
    let fit = analysis::fit_imbalance_dispersion(variants, &config)
        .map_err(|e| PyRuntimeError::new_err(format!("Dispersion fitting failed: {e}")))?;

    let result = PyDict::new(py);
    result.set_item("method", method)?;
    result.set_item("rho", fit.parameters.rho())?;
    result.set_item(
        "linear_d1",
        fit.parameters.linear_params().map(|params| params.0),
    )?;
    result.set_item(
        "linear_d2",
        fit.parameters.linear_params().map(|params| params.1),
    )?;
    result.set_item(
        "linear_depth_center",
        fit.parameters
            .linear_depth_standardization()
            .map(|standardization| standardization.0),
    )?;
    result.set_item(
        "linear_depth_scale",
        fit.parameters
            .linear_depth_standardization()
            .map(|standardization| standardization.1),
    )?;
    result.set_item("n_observations", fit.n_observations)?;
    Ok(result.unbind().into_any())
}

/// Analyze one independent count table and return dispersion metadata.
///
/// Supplying `rho` (single) or all four linear parameters runs a fixed-nuisance
/// test. Omitting them retains the legacy fit-and-test behavior.
#[pyfunction]
#[pyo3(signature = (tsv_path, min_count=10, pseudocount=1, phased=false, region_col=None, method="single", rho=None, linear_d1=None, linear_d2=None, linear_depth_center=None, linear_depth_scale=None, exact_snv=None))]
#[allow(clippy::too_many_arguments)]
fn analyze_imbalance_run(
    py: Python,
    tsv_path: &str,
    min_count: u32,
    pseudocount: u32,
    phased: bool,
    region_col: Option<String>,
    method: &str,
    rho: Option<f64>,
    linear_d1: Option<f64>,
    linear_d2: Option<f64>,
    linear_depth_center: Option<f64>,
    linear_depth_scale: Option<f64>,
    exact_snv: Option<bool>,
) -> PyResult<Py<PyAny>> {
    use pyo3::types::PyDict;

    let analysis_method = parse_analysis_method(method)?;
    let fixed_parameters = fixed_dispersion_parameters(
        analysis_method,
        rho,
        linear_d1,
        linear_d2,
        linear_depth_center,
        linear_depth_scale,
    )?;
    let fixed_nuisance = fixed_parameters.is_some();
    let exact_snv = exact_snv.unwrap_or(fixed_nuisance && region_col.is_none());
    if exact_snv && !fixed_nuisance {
        return Err(PyRuntimeError::new_err(
            "exact_snv=True requires nuisance parameters from fit_imbalance_dispersion()",
        ));
    }
    let config = analysis::AnalysisConfig {
        min_count,
        pseudocount,
        method: analysis_method,
        phased,
    };
    let variants = read_analysis_variants(tsv_path, region_col.as_deref())?;
    let output = match fixed_parameters {
        Some(parameters) => analysis::analyze_imbalance_with_fixed_dispersion(
            variants, &config, parameters, exact_snv,
        ),
        None => analysis::analyze_imbalance_detailed(variants, &config),
    }
    .map_err(|e| PyRuntimeError::new_err(format!("Analysis failed: {e}")))?;
    if output.results.is_empty() {
        return Err(PyRuntimeError::new_err(
            "No variants passed the count threshold; no results were produced",
        ));
    }

    let result = PyDict::new(py);
    result.set_item("results", results_to_py(py, output.results)?)?;
    result.set_item("method", method)?;
    result.set_item("rho", output.rho)?;
    result.set_item("linear_d1", output.linear_params.map(|params| params.0))?;
    result.set_item("linear_d2", output.linear_params.map(|params| params.1))?;
    result.set_item("linear_depth_center", output.linear_depth_center)?;
    result.set_item("linear_depth_scale", output.linear_depth_scale)?;
    result.set_item("n_observations", output.n_observations)?;
    result.set_item("requested_phased", phased)?;
    result.set_item("effective_phased", output.effective_phased)?;
    result.set_item(
        "nuisance_source",
        if fixed_nuisance { "fixed" } else { "fitted" },
    )?;
    result.set_item("exact_snv", exact_snv)?;
    result.set_item(
        "pvalue_method",
        if exact_snv {
            "beta_binomial_exact_two_sided"
        } else {
            "chi_square_lrt"
        },
    )?;
    Ok(result.unbind().into_any())
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

    Ok((
        stats.remap_reads,
        stats.keep_reads,
        stats.unique_remap_names,
    ))
}

// ============================================================================
// Helper: UnifiedStats → PyDict
// NOTE: Update this function when adding fields to UnifiedStats.
// ============================================================================

fn stats_to_pydict(py: Python, stats: &unified_pipeline::UnifiedStats) -> PyResult<Py<PyAny>> {
    use pyo3::types::PyDict;
    let d = PyDict::new(py);
    d.set_item("total_reads", stats.total_reads)?;
    d.set_item("pairs_processed", stats.pairs_processed)?;
    d.set_item("pairs_with_variants", stats.pairs_with_variants)?;
    d.set_item("pairs_with_snvs_only", stats.pairs_with_snvs_only)?;
    d.set_item("pairs_with_indels_only", stats.pairs_with_indels_only)?;
    d.set_item(
        "pairs_with_snvs_and_indels",
        stats.pairs_with_snvs_and_indels,
    )?;
    d.set_item("haplotypes_written", stats.haplotypes_written)?;
    d.set_item("pairs_kept", stats.pairs_kept)?;
    d.set_item("pairs_keep_no_flip", stats.pairs_keep_no_flip)?;
    d.set_item("pairs_skipped_unmappable", stats.pairs_skipped_unmappable)?;
    d.set_item("pairs_haplotype_failed", stats.pairs_haplotype_failed)?;
    d.set_item("orphan_reads", stats.orphan_reads)?;
    d.set_item("tree_build_ms", stats.tree_build_ms)?;
    d.set_item("bam_stream_ms", stats.bam_stream_ms)?;
    d.set_item("overlap_query_ms", stats.overlap_query_ms)?;
    d.set_item("pair_process_ms", stats.pair_process_ms)?;
    d.set_item("send_ms", stats.send_ms)?;
    d.set_item("writer_thread_ms", stats.writer_thread_ms)?;
    Ok(d.unbind().into_any())
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
#[pyo3(signature = (bam_path, bed_path, out_r1, out_r2, max_seqs=64, threads=8, channel_buffer=50000, compression_threads=1, compress_output=true, indel_mode=false, max_indel_size=50, keep_no_flip_names_path=None, remap_names_path=None, pair_buffer_reserve=100000))]
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
    compress_output: bool,
    indel_mode: bool,
    max_indel_size: usize,
    keep_no_flip_names_path: Option<String>,
    remap_names_path: Option<String>,
    pair_buffer_reserve: usize,
) -> PyResult<Py<PyAny>> {
    let config = unified_pipeline::UnifiedConfig {
        read_threads: threads,
        max_seqs,
        pair_buffer_reserve,
        channel_buffer,
        compression_threads,
        compress_output,
        indel_mode,
        max_indel_size,
        keep_no_flip_names_path,
        remap_names_path,
    };

    let stats = unified_pipeline::unified_make_reads(bam_path, bed_path, out_r1, out_r2, &config)
        .map_err(|e| PyRuntimeError::new_err(format!("Unified pipeline failed: {}", e)))?;

    stats_to_pydict(py, &stats)
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
#[pyo3(signature = (bam_path, bed_path, out_r1, out_r2, max_seqs=64, threads=8, channel_buffer=50000, compression_threads=1, compress_output=true, indel_mode=false, max_indel_size=50, keep_no_flip_names_path=None, remap_names_path=None, pair_buffer_reserve=100000))]
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
    compress_output: bool,
    indel_mode: bool,
    max_indel_size: usize,
    keep_no_flip_names_path: Option<String>,
    remap_names_path: Option<String>,
    pair_buffer_reserve: usize,
) -> PyResult<Py<PyAny>> {
    let config = unified_pipeline::UnifiedConfig {
        read_threads: threads,
        max_seqs,
        pair_buffer_reserve,
        channel_buffer,
        compression_threads,
        compress_output,
        indel_mode,
        max_indel_size,
        keep_no_flip_names_path,
        remap_names_path,
    };

    let run = || {
        unified_pipeline::unified_make_reads_parallel(bam_path, bed_path, out_r1, out_r2, &config)
    };

    // Use a per-call Rayon thread pool so repeated calls can use different thread counts.
    let stats = if threads > 0 {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .map_err(|e| {
                PyRuntimeError::new_err(format!("Failed to build Rayon thread pool: {}", e))
            })?;
        pool.install(run)
    } else {
        run()
    }
    .map_err(|e| PyRuntimeError::new_err(format!("Parallel unified pipeline failed: {}", e)))?;

    stats_to_pydict(py, &stats)
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
#[pyo3(signature = (vcf_path, bed_path, samples=None, het_only=true, include_indels=false, max_indel_len=10, include_genotypes=true, biallelic_only=true))]
fn vcf_to_bed_py(
    vcf_path: &str,
    bed_path: &str,
    samples: Option<Vec<String>>,
    het_only: bool,
    include_indels: bool,
    max_indel_len: usize,
    include_genotypes: bool,
    biallelic_only: bool,
) -> PyResult<usize> {
    let config = vcf_to_bed::VcfToBedConfig {
        samples,
        het_only,
        include_indels,
        max_indel_len,
        include_genotypes,
        biallelic_only,
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
) -> PyResult<Py<PyAny>> {
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
                let tuple = pyo3::types::PyTuple::new(py, &[h1.as_str(), h2.as_str()])?;
                alleles_list.append(&tuple)?;
            }
            span_dict.set_item("sample_alleles", alleles_list)?;

            py_list.append(span_dict)?;
        }

        py_dict.set_item(pyo3::types::PyBytes::new(py, read_name), py_list)?;
    }

    Ok(py_dict.unbind().into_any())
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
        bam_path, &variants, chrom, out_r1, out_r2, max_seqs,
    )
    .map_err(|e| PyRuntimeError::new_err(format!("Failed to swap alleles: {}", e)))?;

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
fn wasp2_rust(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;

    // Legacy test function
    m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;

    // Counting module (IMPLEMENTED)
    m.add_class::<BamCounter>()?;

    // Single-cell per-barcode counting module (IMPLEMENTED)
    m.add_class::<BamCounterSC>()?;

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
    // Mapping filter with explicit sidecar argument (CIGAR-aware expected positions)
    m.add_function(wrap_pyfunction!(filter_bam_wasp_with_sidecar, m)?)?;

    // BAM filtering by variant overlap (replaces samtools process_bam, 4-5x faster)
    m.add_function(wrap_pyfunction!(filter_bam_by_variants_py, m)?)?;

    // Unified single-pass pipeline (replaces filter + intersect + remap, 5x faster)
    m.add_function(wrap_pyfunction!(unified_make_reads_py, m)?)?;

    // Parallel unified pipeline (3-8x speedup over sequential, requires BAM index)
    m.add_function(wrap_pyfunction!(unified_make_reads_parallel_py, m)?)?;

    // Analysis module (beta-binomial allelic imbalance detection)
    m.add_function(wrap_pyfunction!(analyze_imbalance, m)?)?;
    m.add_function(wrap_pyfunction!(fit_imbalance_dispersion, m)?)?;
    m.add_function(wrap_pyfunction!(analyze_imbalance_run, m)?)?;

    Ok(())
}

/// Explicit binding exposing expected_sidecar argument (CIGAR-aware expected positions)
#[pyfunction]
#[pyo3(signature = (to_remap_bam, remapped_bam, remap_keep_bam, keep_read_file=None, threads=1, same_locus_slop=0, expected_sidecar=None))]
fn filter_bam_wasp_with_sidecar(
    to_remap_bam: String,
    remapped_bam: String,
    remap_keep_bam: String,
    keep_read_file: Option<String>,
    threads: usize,
    same_locus_slop: i64,
    expected_sidecar: Option<String>,
) -> PyResult<(u64, u64, u64)> {
    mapping_filter::filter_bam_wasp(
        to_remap_bam,
        remapped_bam,
        remap_keep_bam,
        keep_read_file,
        threads,
        same_locus_slop,
        expected_sidecar,
    )
}

#[cfg(test)]
mod tests {
    use super::parse_numeric_phased_gt;

    #[test]
    fn test_parse_numeric_phased_gt_accepts_only_aho_numeric_hets() {
        assert_eq!(parse_numeric_phased_gt("0|1"), Some(0));
        assert_eq!(parse_numeric_phased_gt("1|0"), Some(1));

        for gt in ["A|G", "G|A", "0/1", "0|0", "1|1", "./.", ".|."] {
            assert_eq!(
                parse_numeric_phased_gt(gt),
                None,
                "{} should not be treated as phased GT",
                gt
            );
        }
    }
}
