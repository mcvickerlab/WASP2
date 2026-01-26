"""Generate allele-swapped reads for remapping.

Provides functions for creating FASTQ files with haplotype-swapped reads
that need to be remapped to check for mapping bias.
"""

from __future__ import annotations

import shutil
import tempfile
from pathlib import Path

import pysam

# Rust acceleration (required; no fallback)
from wasp2_rust import remap_all_chromosomes, remap_chromosome, remap_chromosome_multi


def _write_remap_bam_rust_optimized(
    bam_file: str,
    intersect_file: str,
    r1_out: str,
    r2_out: str,
    max_seqs: int = 64,
    parallel: bool = True,
) -> None:
    """
    Optimized Rust remapping - parses intersect file ONCE, processes chromosomes in parallel.

    This is the fastest implementation:
    - Parses intersect file once (22x fewer parse operations for RNA-seq)
    - Uses rayon for parallel chromosome processing (4-8x speedup with 8 cores)
    - Total expected speedup: ~100x for large RNA-seq datasets
    """
    import inspect

    print(
        f"Using optimized Rust remapper (parse-once, {'parallel' if parallel else 'sequential'})..."
    )

    # Check if the Rust function accepts 'parallel' parameter (backward compatibility)
    sig = inspect.signature(remap_all_chromosomes)
    has_parallel_param = "parallel" in sig.parameters

    if has_parallel_param:
        # New version with parallel parameter
        pairs, haps = remap_all_chromosomes(
            bam_file, intersect_file, r1_out, r2_out, max_seqs=max_seqs, parallel=parallel
        )
    else:
        # Old version without parallel parameter (always runs in parallel)
        print("  Note: Using Rust version without 'parallel' parameter (parallel by default)")
        pairs, haps = remap_all_chromosomes(
            bam_file, intersect_file, r1_out, r2_out, max_seqs=max_seqs
        )

    print(f"\n✅ Rust remapper (optimized): {pairs} pairs → {haps} haplotypes")
    print(f"Reads to remap written to:\n{r1_out}\n{r2_out}")


def _write_remap_bam_rust(
    bam_file: str, intersect_file: str, r1_out: str, r2_out: str, max_seqs: int = 64
) -> None:
    """Rust-accelerated remapping implementation (5-7x faster than Python) - LEGACY per-chromosome version"""
    # Get chromosomes that have variants in the intersect file
    # This avoids processing ~170 empty chromosomes (major speedup!)
    intersect_chroms = set()
    with open(intersect_file) as f:
        for line in f:
            chrom = line.split("\t")[0]
            intersect_chroms.add(chrom)

    # Filter BAM chromosomes to only those with variants
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        chromosomes = [c for c in bam.header.references if c in intersect_chroms]

    print(
        f"Processing {len(chromosomes)} chromosomes with variants (filtered from {len(intersect_chroms)} in intersect)"
    )

    # Create temp directory for per-chromosome outputs
    with tempfile.TemporaryDirectory() as tmpdir:
        total_pairs = 0
        total_haps = 0

        # Process each chromosome with Rust
        for chrom in chromosomes:
            chrom_r1 = f"{tmpdir}/{chrom}_r1.fq"
            chrom_r2 = f"{tmpdir}/{chrom}_r2.fq"

            try:
                pairs, haps = remap_chromosome(
                    bam_file, intersect_file, chrom, chrom_r1, chrom_r2, max_seqs=max_seqs
                )
                total_pairs += pairs
                total_haps += haps
                if pairs > 0:
                    print(f"  {chrom}: {pairs} pairs → {haps} haplotypes")
            except Exception as e:
                print(f"  {chrom}: Error - {e}")
                continue

        # Concatenate all R1 files
        r1_files = sorted(Path(tmpdir).glob("*_r1.fq"))
        with open(r1_out, "wb") as outfile:
            for f in r1_files:
                with open(f, "rb") as infile:
                    shutil.copyfileobj(infile, outfile)

        # Concatenate all R2 files
        r2_files = sorted(Path(tmpdir).glob("*_r2.fq"))
        with open(r2_out, "wb") as outfile:
            for f in r2_files:
                with open(f, "rb") as infile:
                    shutil.copyfileobj(infile, outfile)

        print(f"\n✅ Rust remapper: {total_pairs} pairs → {total_haps} haplotypes")
        print(f"Reads to remapped written to \n{r1_out}\n{r2_out}")


def _write_remap_bam_rust_multi(
    bam_file: str,
    intersect_file: str,
    r1_out: str,
    r2_out: str,
    num_samples: int,
    max_seqs: int = 64,
) -> None:
    """Rust-accelerated multi-sample remapping implementation"""
    # Get chromosomes that have variants in the intersect file
    intersect_chroms = set()
    with open(intersect_file) as f:
        for line in f:
            chrom = line.split("\t")[0]
            intersect_chroms.add(chrom)

    # Filter BAM chromosomes to only those with variants
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        chromosomes = [c for c in bam.header.references if c in intersect_chroms]

    print(f"Processing {len(chromosomes)} chromosomes with variants ({num_samples} samples)")

    # Create temp directory for per-chromosome outputs
    with tempfile.TemporaryDirectory() as tmpdir:
        total_pairs = 0
        total_haps = 0

        # Process each chromosome with Rust multi-sample
        for chrom in chromosomes:
            chrom_r1 = f"{tmpdir}/{chrom}_r1.fq"
            chrom_r2 = f"{tmpdir}/{chrom}_r2.fq"

            try:
                pairs, haps = remap_chromosome_multi(
                    bam_file,
                    intersect_file,
                    chrom,
                    chrom_r1,
                    chrom_r2,
                    num_samples=num_samples,
                    max_seqs=max_seqs,
                )
                total_pairs += pairs
                total_haps += haps
                if pairs > 0:
                    print(f"  {chrom}: {pairs} pairs → {haps} haplotypes")
            except Exception as e:
                print(f"  {chrom}: Error - {e}")
                continue

        # Concatenate all R1 files
        r1_files = sorted(Path(tmpdir).glob("*_r1.fq"))
        with open(r1_out, "wb") as outfile:
            for f in r1_files:
                with open(f, "rb") as infile:
                    shutil.copyfileobj(infile, outfile)

        # Concatenate all R2 files
        r2_files = sorted(Path(tmpdir).glob("*_r2.fq"))
        with open(r2_out, "wb") as outfile:
            for f in r2_files:
                with open(f, "rb") as infile:
                    shutil.copyfileobj(infile, outfile)

        print(f"\n✅ Rust multi-sample remapper: {total_pairs} pairs → {total_haps} haplotypes")
        print(f"Reads to remapped written to \n{r1_out}\n{r2_out}")


def write_remap_bam(
    bam_file: str,
    intersect_file: str,
    r1_out: str,
    r2_out: str,
    samples: list[str],
    max_seqs: int = 64,
    include_indels: bool = False,
    insert_qual: int = 30,
) -> None:
    """Rust-accelerated remapping - parses intersect file once, processes chromosomes in parallel.

    Uses Rust acceleration (required; no fallback).

    Args:
        bam_file: Input BAM file
        intersect_file: Intersect BED file
        r1_out: Output FASTQ for read 1
        r2_out: Output FASTQ for read 2
        samples: List of sample IDs
        max_seqs: Maximum haplotype sequences per read pair
        include_indels: Include indels in remapping (not yet supported in Rust)
        insert_qual: Quality score for inserted bases (not yet supported in Rust)
    """
    num_samples = len(samples)

    if num_samples == 1:
        # Single sample: use optimized all-chromosome Rust
        _write_remap_bam_rust_optimized(
            bam_file, intersect_file, r1_out, r2_out, max_seqs, parallel=True
        )
    else:
        # Multi-sample: use per-chromosome Rust
        _write_remap_bam_rust_multi(bam_file, intersect_file, r1_out, r2_out, num_samples, max_seqs)
