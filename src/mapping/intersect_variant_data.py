"""Variant intersection and BAM filtering utilities.

Provides functions for converting variants to BED format, filtering BAM files
by variant overlap, and creating intersection files for the WASP pipeline.
"""

from __future__ import annotations

import logging
import os
import subprocess
from pathlib import Path

import numpy as np
import polars as pl
import pysam

# Multi-format variant support
from wasp2.io import variants_to_bed as _variants_to_bed

# Rust acceleration (required; no fallback)
from wasp2_rust import filter_bam_by_variants_py as _rust_filter_bam
from wasp2_rust import intersect_bam_bed as _rust_intersect
from wasp2_rust import intersect_bam_bed_multi as _rust_intersect_multi

logger = logging.getLogger(__name__)


def vcf_to_bed(
    vcf_file: str | Path,
    out_bed: str | Path,
    samples: list[str] | None = None,
    include_indels: bool = False,
    max_indel_len: int = 10,
) -> str:
    """Convert variant file to BED format.

    Supports VCF, VCF.GZ, BCF, and PGEN formats via the VariantSource API.

    Note: Parameter name 'vcf_file' is kept for backward compatibility,
    but accepts any supported variant format (VCF, BCF, PGEN).

    Args:
        vcf_file: Path to variant file (VCF, VCF.GZ, BCF, or PGEN)
        out_bed: Output BED file path
        samples: Optional list of sample IDs. If provided, filters to het sites.
        include_indels: Include indels in addition to SNPs
        max_indel_len: Maximum indel length (bp) to include

    Returns:
        Path to output BED file as string
    """
    # Use new unified interface with Rust VCF parser (5-6x faster than bcftools)
    # include_gt=True for mapping (needs genotypes for allele assignment)
    result = _variants_to_bed(
        variant_file=vcf_file,
        out_bed=out_bed,
        samples=samples,
        include_gt=True,
        het_only=bool(samples),
        include_indels=include_indels,
        max_indel_len=max_indel_len,
        use_legacy=False,  # Use Rust VCF parser (5-6x faster than bcftools)
    )
    return str(result)


def process_bam(
    bam_file: str,
    vcf_bed: str,
    remap_bam: str,
    remap_reads: str,
    keep_bam: str,
    is_paired: bool = True,
    threads: int = 1,
) -> str:
    """Filter BAM by variant overlap, splitting into remap/keep BAMs.

    Uses Rust acceleration (~2x faster than samtools).

    Args:
        bam_file: Input BAM file (coordinate-sorted)
        vcf_bed: Variant BED file from vcf_to_bed
        remap_bam: Output BAM for reads needing remapping
        remap_reads: Output file for unique read names
        keep_bam: Output BAM for reads not needing remapping
        is_paired: Whether reads are paired-end
        threads: Number of threads

    Returns:
        Path to remap BAM file
    """
    logger.info("Using Rust acceleration for BAM filtering...")
    remap_count, keep_count, unique_names = _rust_filter_bam(
        bam_file, vcf_bed, remap_bam, keep_bam, is_paired, threads
    )
    logger.info(
        "Rust filter: %s remap, %s keep, %s unique names",
        f"{remap_count:,}", f"{keep_count:,}", f"{unique_names:,}",
    )

    # Write read names file for compatibility
    with pysam.AlignmentFile(remap_bam, "rb") as bam, open(remap_reads, "w") as f:
        names = {read.query_name for read in bam.fetch(until_eof=True) if read.query_name is not None}
        f.write("\n".join(names))

    # Sort the remap BAM (Rust outputs unsorted)
    remap_bam_tmp = remap_bam + ".sorting.tmp"
    subprocess.run(
        ["samtools", "sort", "-@", str(threads), "-o", remap_bam_tmp, remap_bam], check=True
    )
    os.rename(remap_bam_tmp, remap_bam)

    subprocess.run(["samtools", "index", "-@", str(threads), str(remap_bam)], check=True)

    return remap_bam


def intersect_reads(remap_bam: str, vcf_bed: str, out_bed: str, num_samples: int = 1) -> str:
    """Intersect BAM reads with variant BED file.

    Uses Rust/coitrees (15-30x faster than pybedtools).

    Args:
        remap_bam: Path to BAM file with reads overlapping variants
        vcf_bed: Path to BED file with variant positions
        out_bed: Output path for intersection results
        num_samples: Number of sample genotype columns in BED file (default 1)

    Returns:
        Path to output BED file
    """
    if num_samples == 1:
        logger.info("Using Rust acceleration for intersection...")
        count = _rust_intersect(remap_bam, vcf_bed, out_bed)
    else:
        logger.info("Using Rust multi-sample intersection (%d samples)...", num_samples)
        count = _rust_intersect_multi(remap_bam, vcf_bed, out_bed, num_samples)
    logger.info("Rust intersect: %d overlaps found", count)
    return out_bed


def make_intersect_df(
    intersect_file: str,
    samples: list[str],
    is_paired: bool = True,
) -> pl.DataFrame:
    """Parse intersection file into a typed polars DataFrame.

    Parameters
    ----------
    intersect_file : str
        Path to intersection BED file.
    samples : list[str]
        List of sample column names.
    is_paired : bool, optional
        Whether reads are paired-end, by default True.

    Returns
    -------
    pl.DataFrame
        Parsed intersection data with alleles split by sample.
    """
    # Create Dataframe
    df = pl.scan_csv(intersect_file, separator="\t", has_header=False, infer_schema_length=0)

    # Parse sample data
    num_samps = len(samples)

    subset_cols = [df.columns[i] for i in np.r_[0, 3, 1, 2, -num_samps:0]]
    new_cols = ["chrom", "read", "start", "stop", *samples]

    rename_cols = dict(zip(subset_cols, new_cols))

    base_schema = [
        pl.col("chrom").cast(pl.Categorical),
        pl.col("read").cast(pl.Utf8),
        pl.col("start").cast(pl.UInt32),
        pl.col("stop").cast(pl.UInt32),
    ]

    sample_schema = [pl.col(samp).cast(pl.Utf8) for samp in samples]
    col_schema = [*base_schema, *sample_schema]

    # Make sure types are correct
    df = df.select(subset_cols).rename(rename_cols).with_columns(col_schema)

    expr_list = []
    cast_list = []

    for s in samples:
        a1 = f"{s}_a1"
        a2 = f"{s}_a2"

        # Add split per sample
        expr_list.append(pl.col(s).str.split_exact(by="|", n=1).struct.rename_fields([a1, a2]))

        # cast new gt cols
        cast_list.append(pl.col(a1).cast(pl.Categorical))
        cast_list.append(pl.col(a2).cast(pl.Categorical))

    # Split mate expr
    expr_list.append(
        pl.col("read").str.split_exact(by="/", n=1).struct.rename_fields(["read", "mate"])
    )

    cast_list.append(pl.col("mate").cast(pl.UInt8))

    df = df.with_columns(expr_list).unnest([*samples, "read"]).with_columns(cast_list)

    # should i remove instead of keep first?
    df = df.unique(
        ["chrom", "read", "mate", "start", "stop"], keep="first"
    )  # Doesnt remove dup snp in pair?

    return df.collect()
