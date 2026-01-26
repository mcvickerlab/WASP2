"""Type stubs for wasp2_rust - PyO3 Rust acceleration module.

This module provides high-performance implementations of bottleneck functions
for allele-specific analysis in WASP2.
"""

from typing import TypedDict

class UnifiedStats(TypedDict):
    """Statistics from unified pipeline execution."""

    total_reads: int
    pairs_processed: int
    pairs_with_variants: int
    pairs_with_snvs_only: int
    pairs_with_indels_only: int
    pairs_with_snvs_and_indels: int
    haplotypes_written: int
    pairs_kept: int
    pairs_keep_no_flip: int
    pairs_skipped_unmappable: int
    pairs_haplotype_failed: int
    orphan_reads: int
    tree_build_ms: float
    bam_stream_ms: float
    overlap_query_ms: float
    pair_process_ms: float
    send_ms: float
    writer_thread_ms: float


class ImbalanceResult(TypedDict):
    """Result from allelic imbalance analysis."""

    region: str
    ref_count: int
    alt_count: int
    N: int
    snp_count: int
    null_ll: float
    alt_ll: float
    mu: float
    lrt: float
    pval: float
    fdr_pval: float


class VariantSpan(TypedDict):
    """Variant span information from intersection parsing."""

    chrom: str
    start: int
    stop: int
    vcf_start: int
    vcf_stop: int
    mate: int
    hap1: str
    hap2: str


class MultiSampleVariantSpan(TypedDict):
    """Multi-sample variant span information."""

    chrom: str
    start: int
    stop: int
    vcf_start: int
    vcf_stop: int
    mate: int
    ref_allele: str
    alt_allele: str
    sample_alleles: list[tuple[str, str]]


class BamCounter:
    """Fast BAM allele counter using Rust/htslib.

    Parameters
    ----------
    bam_path : str
        Path to BAM file (must be indexed).

    Examples
    --------
    >>> counter = BamCounter("sample.bam")
    >>> regions = [("chr1", 12345, "A", "G"), ("chr1", 12400, "C", "T")]
    >>> counts = counter.count_alleles(regions, min_qual=0, threads=4)
    >>> for ref, alt, other in counts:
    ...     print(f"ref={ref}, alt={alt}, other={other}")
    """

    def __init__(self, bam_path: str) -> None: ...
    def count_alleles(
        self,
        regions: list[tuple[str, int, str, str]],
        min_qual: int = 0,
        threads: int = 1,
    ) -> list[tuple[int, int, int]]:
        """Count alleles at specified positions.

        Parameters
        ----------
        regions : list[tuple[str, int, str, str]]
            List of (chrom, pos, ref, alt) tuples.
        min_qual : int, optional
            Minimum base quality threshold, by default 0.
        threads : int, optional
            Number of threads to use, by default 1.

        Returns
        -------
        list[tuple[int, int, int]]
            List of (ref_count, alt_count, other_count) tuples.
        """
        ...


# Test function
def sum_as_string(a: int, b: int) -> str:
    """Simple test function to verify PyO3 is working."""
    ...


# BAM-BED intersection functions
def intersect_bam_bed(bam_path: str, bed_path: str, out_path: str) -> int:
    """Intersect BAM reads with variant BED file using coitrees.

    Parameters
    ----------
    bam_path : str
        Path to sorted BAM file.
    bed_path : str
        Path to variant BED file (chrom, start, stop, ref, alt, GT).
    out_path : str
        Output path for intersections.

    Returns
    -------
    int
        Number of intersections found.
    """
    ...


def intersect_bam_bed_multi(
    bam_path: str,
    bed_path: str,
    out_path: str,
    num_samples: int,
) -> int:
    """Intersect BAM reads with multi-sample variant BED file.

    Parameters
    ----------
    bam_path : str
        Path to sorted BAM file.
    bed_path : str
        Path to variant BED file with multiple GT columns.
    out_path : str
        Output path for intersections.
    num_samples : int
        Number of sample genotype columns in BED.

    Returns
    -------
    int
        Number of intersections found.
    """
    ...


# VCF/BCF to BED conversion
def vcf_to_bed(
    vcf_path: str,
    bed_path: str,
    samples: list[str] | None = None,
    het_only: bool = True,
    include_indels: bool = False,
    max_indel_len: int = 10,
    include_genotypes: bool = True,
) -> int:
    """Convert VCF/BCF to BED format using noodles.

    Parameters
    ----------
    vcf_path : str
        Path to VCF/BCF file.
    bed_path : str
        Output BED file path.
    samples : list[str] | None, optional
        List of sample names to extract (None = all), by default None.
    het_only : bool, optional
        Only output heterozygous sites, by default True.
    include_indels : bool, optional
        Include indels, not just SNPs, by default False.
    max_indel_len : int, optional
        Maximum indel length to include, by default 10.
    include_genotypes : bool, optional
        Include genotype column in output, by default True.

    Returns
    -------
    int
        Number of variants written to BED file.
    """
    ...


# Intersection parsing functions
def parse_intersect_bed(intersect_bed: str) -> dict[bytes, list[VariantSpan]]:
    """Parse intersection BED file.

    Parameters
    ----------
    intersect_bed : str
        Path to bedtools intersect output.

    Returns
    -------
    dict[bytes, list[VariantSpan]]
        Dictionary mapping read names (bytes) to list of variant spans.
    """
    ...


def parse_intersect_bed_multi(
    intersect_bed: str,
    num_samples: int,
) -> dict[bytes, list[MultiSampleVariantSpan]]:
    """Parse multi-sample intersection BED file.

    Parameters
    ----------
    intersect_bed : str
        Path to intersection BED file.
    num_samples : int
        Number of sample genotype columns.

    Returns
    -------
    dict[bytes, list[MultiSampleVariantSpan]]
        Dictionary mapping read names to variant spans with all sample genotypes.
    """
    ...


# Remapping functions
def remap_chromosome(
    bam_path: str,
    intersect_bed: str,
    chrom: str,
    out_r1: str,
    out_r2: str,
    max_seqs: int = 64,
) -> tuple[int, int]:
    """Remap reads for a single chromosome.

    Parameters
    ----------
    bam_path : str
        Path to BAM file with reads to remap.
    intersect_bed : str
        Path to bedtools intersect output.
    chrom : str
        Chromosome to process (e.g., "chr10").
    out_r1 : str
        Output path for read 1 FASTQ.
    out_r2 : str
        Output path for read 2 FASTQ.
    max_seqs : int, optional
        Maximum haplotype sequences per read pair, by default 64.

    Returns
    -------
    tuple[int, int]
        (pairs_processed, haplotypes_generated).
    """
    ...


def remap_chromosome_multi(
    bam_path: str,
    intersect_bed: str,
    chrom: str,
    out_r1: str,
    out_r2: str,
    num_samples: int,
    max_seqs: int = 64,
) -> tuple[int, int]:
    """Remap reads for a single chromosome - multi-sample version.

    Parameters
    ----------
    bam_path : str
        Path to BAM file with reads to remap.
    intersect_bed : str
        Path to bedtools intersect output (multi-sample format).
    chrom : str
        Chromosome to process (e.g., "chr10").
    out_r1 : str
        Output path for read 1 FASTQ.
    out_r2 : str
        Output path for read 2 FASTQ.
    num_samples : int
        Number of samples in the intersection BED.
    max_seqs : int, optional
        Maximum haplotype sequences per read pair, by default 64.

    Returns
    -------
    tuple[int, int]
        (pairs_processed, haplotypes_generated).
    """
    ...


def remap_all_chromosomes(
    bam_path: str,
    intersect_bed: str,
    out_r1: str,
    out_r2: str,
    max_seqs: int = 64,
    parallel: bool = True,
    num_threads: int = 0,
) -> tuple[int, int]:
    """Remap all chromosomes in parallel.

    Parameters
    ----------
    bam_path : str
        Path to BAM file.
    intersect_bed : str
        Path to bedtools intersect output.
    out_r1 : str
        Output path for read 1 FASTQ.
    out_r2 : str
        Output path for read 2 FASTQ.
    max_seqs : int, optional
        Maximum haplotype sequences per read pair, by default 64.
    parallel : bool, optional
        Use parallel processing, by default True.
    num_threads : int, optional
        Number of threads (0 = auto-detect), by default 0.

    Returns
    -------
    tuple[int, int]
        (pairs_processed, haplotypes_generated).
    """
    ...


# BAM filtering functions
def filter_bam_wasp(
    to_remap_bam: str,
    remapped_bam: str,
    remap_keep_bam: str,
    keep_read_file: str | None = None,
    threads: int = 1,
    same_locus_slop: int = 0,
    expected_sidecar: str | None = None,
) -> tuple[int, int, int]:
    """Filter BAM reads using WASP mapping filter.

    Parameters
    ----------
    to_remap_bam : str
        Path to BAM file with reads to remap.
    remapped_bam : str
        Path to remapped BAM file.
    remap_keep_bam : str
        Output path for kept reads.
    keep_read_file : str | None, optional
        Output path for read names, by default None.
    threads : int, optional
        Number of threads, by default 1.
    same_locus_slop : int, optional
        Slop for same-locus detection, by default 0.
    expected_sidecar : str | None, optional
        Path to expected positions sidecar file, by default None.

    Returns
    -------
    tuple[int, int, int]
        (kept_count, filtered_count, total_count).
    """
    ...


def filter_bam_wasp_with_sidecar(
    to_remap_bam: str,
    remapped_bam: str,
    remap_keep_bam: str,
    keep_read_file: str | None = None,
    threads: int = 1,
    same_locus_slop: int = 0,
    expected_sidecar: str | None = None,
) -> tuple[int, int, int]:
    """Filter BAM reads using WASP mapping filter with explicit sidecar argument.

    Parameters
    ----------
    to_remap_bam : str
        Path to BAM file with reads to remap.
    remapped_bam : str
        Path to remapped BAM file.
    remap_keep_bam : str
        Output path for kept reads.
    keep_read_file : str | None, optional
        Output path for read names, by default None.
    threads : int, optional
        Number of threads, by default 1.
    same_locus_slop : int, optional
        Slop for same-locus detection, by default 0.
    expected_sidecar : str | None, optional
        Path to expected positions sidecar file (CIGAR-aware), by default None.

    Returns
    -------
    tuple[int, int, int]
        (kept_count, filtered_count, total_count).
    """
    ...


def filter_bam_by_variants(
    bam_path: str,
    bed_path: str,
    remap_bam_path: str,
    keep_bam_path: str,
    is_paired: bool = True,
    threads: int = 4,
) -> tuple[int, int, int]:
    """Filter BAM by variant overlap.

    Parameters
    ----------
    bam_path : str
        Input BAM file (should be coordinate-sorted).
    bed_path : str
        Variant BED file (chrom, start, stop, ref, alt, GT).
    remap_bam_path : str
        Output BAM for reads needing remapping.
    keep_bam_path : str
        Output BAM for reads not needing remapping.
    is_paired : bool, optional
        Whether reads are paired-end, by default True.
    threads : int, optional
        Number of threads to use, by default 4.

    Returns
    -------
    tuple[int, int, int]
        (remap_count, keep_count, unique_names).
    """
    ...


# Unified pipeline functions
def unified_make_reads(
    bam_path: str,
    bed_path: str,
    out_r1: str,
    out_r2: str,
    max_seqs: int = 64,
    threads: int = 8,
    channel_buffer: int = 50000,
    compression_threads: int = 1,
    compress_output: bool = True,
    indel_mode: bool = False,
    max_indel_size: int = 50,
    keep_no_flip_names_path: str | None = None,
    remap_names_path: str | None = None,
    pair_buffer_reserve: int = 100000,
) -> UnifiedStats:
    """Unified single-pass make-reads pipeline.

    Parameters
    ----------
    bam_path : str
        Input BAM file (should be coordinate-sorted).
    bed_path : str
        Variant BED file (chrom, start, stop, ref, alt, GT).
    out_r1 : str
        Output path for read 1 FASTQ.
    out_r2 : str
        Output path for read 2 FASTQ.
    max_seqs : int, optional
        Maximum haplotype sequences per read pair, by default 64.
    threads : int, optional
        Number of threads to use, by default 8.
    channel_buffer : int, optional
        Channel buffer size for streaming, by default 50000.
    compression_threads : int, optional
        Threads per FASTQ file for gzip, by default 1.
    compress_output : bool, optional
        Whether to gzip output, by default True.
    indel_mode : bool, optional
        Enable indel processing, by default False.
    max_indel_size : int, optional
        Maximum indel size to process, by default 50.
    keep_no_flip_names_path : str | None, optional
        Path to write names of kept reads without flip, by default None.
    remap_names_path : str | None, optional
        Path to write names of remapped reads, by default None.
    pair_buffer_reserve : int, optional
        Pair buffer reserve size, by default 100000.

    Returns
    -------
    UnifiedStats
        Dictionary with processing statistics.
    """
    ...


def unified_make_reads_parallel(
    bam_path: str,
    bed_path: str,
    out_r1: str,
    out_r2: str,
    max_seqs: int = 64,
    threads: int = 8,
    channel_buffer: int = 50000,
    compression_threads: int = 1,
    compress_output: bool = True,
    indel_mode: bool = False,
    max_indel_size: int = 50,
    keep_no_flip_names_path: str | None = None,
    remap_names_path: str | None = None,
    pair_buffer_reserve: int = 100000,
) -> UnifiedStats:
    """Parallel unified pipeline - processes chromosomes in parallel.

    Requires BAM to be coordinate-sorted and indexed (.bai file must exist).
    Falls back to sequential if BAM index is missing.

    Parameters
    ----------
    bam_path : str
        Input BAM file (must be coordinate-sorted and indexed).
    bed_path : str
        Variant BED file (chrom, start, stop, ref, alt, GT).
    out_r1 : str
        Output path for read 1 FASTQ.
    out_r2 : str
        Output path for read 2 FASTQ.
    max_seqs : int, optional
        Maximum haplotype sequences per read pair, by default 64.
    threads : int, optional
        Number of threads to use, by default 8.
    channel_buffer : int, optional
        Channel buffer size for streaming, by default 50000.
    compression_threads : int, optional
        Threads per FASTQ file for gzip, by default 1.
    compress_output : bool, optional
        Whether to gzip output, by default True.
    indel_mode : bool, optional
        Enable indel processing, by default False.
    max_indel_size : int, optional
        Maximum indel size to process, by default 50.
    keep_no_flip_names_path : str | None, optional
        Path to write names of kept reads without flip, by default None.
    remap_names_path : str | None, optional
        Path to write names of remapped reads, by default None.
    pair_buffer_reserve : int, optional
        Pair buffer reserve size, by default 100000.

    Returns
    -------
    UnifiedStats
        Dictionary with processing statistics.
    """
    ...


# Analysis functions
def analyze_imbalance(
    tsv_path: str,
    min_count: int = 10,
    pseudocount: int = 1,
    method: str = "single",
) -> list[ImbalanceResult]:
    """Analyze allelic imbalance using beta-binomial model.

    Parameters
    ----------
    tsv_path : str
        Path to TSV file with allele counts.
    min_count : int, optional
        Minimum total count threshold, by default 10.
    pseudocount : int, optional
        Pseudocount to add to allele counts, by default 1.
    method : str, optional
        Analysis method ("single" or "linear"), by default "single".

    Returns
    -------
    list[ImbalanceResult]
        List of imbalance results per region.
    """
    ...
