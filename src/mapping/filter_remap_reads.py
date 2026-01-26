"""Filter and merge remapped BAM reads using WASP algorithm.

Provides functions for filtering reads that remap to the same locus
after allele swapping and merging with non-remapped reads.
"""

from __future__ import annotations

import subprocess
import timeit

# Rust acceleration (required; no fallback)
from wasp2_rust import filter_bam_wasp


def filt_remapped_reads(
    to_remap_bam: str,
    remapped_bam: str,
    filt_out_bam: str,
    keep_read_file: str | None = None,
    threads: int = 1,
    same_locus_slop: int = 0,
) -> None:
    """Filter remapped reads using WASP algorithm.

    Uses Rust acceleration.

    Args:
        to_remap_bam: Original BAM with reads to remap
        remapped_bam: Remapped BAM with swapped alleles
        filt_out_bam: Output filtered BAM
        keep_read_file: Optional file to write kept read names
        threads: Number of threads for BAM I/O
        same_locus_slop: Tolerance (bp) for same locus test (for indels)
    """
    filter_bam_wasp(
        to_remap_bam,
        remapped_bam,
        filt_out_bam,
        keep_read_file=keep_read_file,
        threads=threads,
        same_locus_slop=same_locus_slop,
    )


def merge_filt_bam(keep_bam: str, remapped_filt_bam: str, out_bam: str, threads: int = 1) -> None:
    """Merge filtered BAM files using samtools (faster than pysam).

    Both input BAMs are already coordinate-sorted, so samtools merge
    produces sorted output without needing an explicit sort step.

    Args:
        keep_bam: BAM with reads that didn't need remapping
        remapped_filt_bam: BAM with filtered remapped reads
        out_bam: Output merged BAM
        threads: Number of threads for samtools
    """
    start_time = timeit.default_timer()

    # Merge using samtools (faster than pysam, inputs are already sorted)
    subprocess.run(
        ["samtools", "merge", "-@", str(threads), "-f", "-o", out_bam, keep_bam, remapped_filt_bam],
        check=True,
    )
    print(f"Merged BAM in {timeit.default_timer() - start_time:.2f} seconds")

    # Index the merged BAM (no sort needed - inputs were already sorted)
    start_index = timeit.default_timer()
    subprocess.run(["samtools", "index", "-@", str(threads), out_bam], check=True)
    print(f"Indexed BAM in {timeit.default_timer() - start_index:.2f} seconds")
