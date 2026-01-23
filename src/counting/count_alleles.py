import os
import timeit
from pathlib import Path
from typing import Optional

import polars as pl

# Try to import Rust acceleration (required; no Python fallback)
try:
    from wasp2_rust import BamCounter as RustBamCounter
    RUST_AVAILABLE = True
except ImportError:
    RUST_AVAILABLE = False

def count_snp_alleles_rust(bam_file, chrom, snp_list, threads: Optional[int] = None):
    """
    Rust-accelerated version of count_snp_alleles

    :param str bam_file: Path to BAM file
    :param str chrom: Chromosome name
    :param snp_list: Iterator of (pos, ref, alt) tuples
    :param int threads: Optional number of threads (default 1 or WASP2_RUST_THREADS env)
    :return list: List of (chrom, pos, ref_count, alt_count, other_count) tuples
    """
    rust_threads_env = os.environ.get("WASP2_RUST_THREADS") if threads is None else None
    try:
        rust_threads = threads if threads is not None else (int(rust_threads_env) if rust_threads_env else 1)
    except ValueError:
        rust_threads = 1
    rust_threads = max(1, rust_threads)

    # Convert snp_list to list of regions for Rust
    regions = [(chrom, pos, ref, alt) for pos, ref, alt in snp_list]

    # Create Rust BAM counter
    counter = RustBamCounter(bam_file)

    # Count alleles (returns list of (ref_count, alt_count, other_count))
    # min_qual=0 matches WASP2 behavior (no quality filtering)
    counts = counter.count_alleles(regions, min_qual=0, threads=rust_threads)

    # Combine with chromosome and position info
    allele_counts = [
        (chrom, pos, ref_count, alt_count, other_count)
        for (_, pos, _, _), (ref_count, alt_count, other_count) in zip(regions, counts)
    ]

    return allele_counts


def make_count_df(bam_file, df, use_rust=True):
    """
    Make DF containing all intersections and allele counts

    :param str bam_file: Path to BAM file
    :param DataFrame df: Dataframe of intersections, output from
        parse_(intersect/gene)_df()
    :param bool use_rust: Use Rust acceleration if available (default: True)
    :return DataFrame: DataFrame of counts
    """
    count_list = []

    chrom_list = df.get_column("chrom").unique(
        maintain_order=True)

    # Require Rust path (no Python fallback)
    if not (use_rust and RUST_AVAILABLE):
        raise RuntimeError(
            "Rust BAM counter not available. Build the extension with "
            "`maturin develop --release` in the WASP2 env."
        )

    rust_threads_env = os.environ.get("WASP2_RUST_THREADS")
    try:
        rust_threads = int(rust_threads_env) if rust_threads_env else 1
    except ValueError:
        rust_threads = 1
    rust_threads = max(1, rust_threads)
    print(f"Using Rust acceleration for BAM counting 🦀 (threads={rust_threads})")

    total_start = timeit.default_timer()

    for chrom in chrom_list:
        chrom_df = df.filter(pl.col("chrom") == chrom)

        snp_list = chrom_df.select(
            ["pos", "ref", "alt"]).unique(
            subset=["pos"], maintain_order=True).iter_rows()

        start = timeit.default_timer()

        try:
            count_list.extend(count_snp_alleles_rust(bam_file, chrom, snp_list, threads=rust_threads))
        except Exception as e:
            print(f"Skipping {chrom}: {e}\n")
        else:
            print(f"{chrom}: Counted {chrom_df.height} SNP's in {timeit.default_timer() - start:.2f} seconds!")

    total_end = timeit.default_timer()
    print(f"Counted all SNP's in {total_end - total_start:.2f} seconds!")

    # Previously used str as chrom instead of cat
    chrom_enum = pl.Enum(df.get_column("chrom").cat.get_categories())

    count_df = pl.DataFrame(
        count_list,
        schema={"chrom": chrom_enum,
                "pos": pl.UInt32,
                "ref_count": pl.UInt16,
                "alt_count": pl.UInt16,
                "other_count": pl.UInt16
               },
        orient="row"
    )

    # possibly find better solution
    df = df.with_columns([pl.col("chrom").cast(chrom_enum)]
                         ).join(count_df, on=["chrom", "pos"], how="left")

    # df = df.join(count_df, on=["chrom", "pos"], how="left")

    return df

# Legacy helper retained for imports in counting/count_alleles_sc.py
def find_read_aln_pos(read, pos):
    """
    Binary search over aligned pairs to find query position for a given reference pos.
    """
    aln_list = read.get_aligned_pairs(True)
    # bisect_left using manual loop to avoid Python <3.10 key support
    lo, hi = 0, len(aln_list)
    while lo < hi:
        mid = (lo + hi) // 2
        if aln_list[mid][1] < pos:
            lo = mid + 1
        else:
            hi = mid
    if lo != len(aln_list) and aln_list[lo][1] == pos:
        return aln_list[lo][0]
    return None
