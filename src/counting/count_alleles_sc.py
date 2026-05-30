"""Single-cell allele counting functions."""

from __future__ import annotations

import logging
import os
import timeit
from collections import defaultdict
from collections.abc import Iterator

import anndata as ad
import numpy as np
import pandas as pd
import polars as pl
from pysam.libcalignmentfile import AlignmentFile
from scipy.sparse import csr_matrix

# Local imports
from .count_alleles import find_read_aln_pos

logger = logging.getLogger(__name__)

# Try to import the Rust single-cell counter (optional; Python fallback always available)
try:
    from wasp2_rust import BamCounterSC as _RustSC

    RUST_SC_AVAILABLE = True
except ImportError:
    RUST_SC_AVAILABLE = False


def _sparse_from_counts(
    counts: defaultdict[tuple[int, int], int],
    shape: tuple[int, int],
) -> csr_matrix:
    if not counts:
        return csr_matrix(shape, dtype=np.uint16)
    return csr_matrix(
        (list(counts.values()), list(zip(*counts.keys()))),
        shape=shape,
        dtype=np.uint16,
    )


class CountStatsSC:
    """Container for mutable single-cell counting statistics.

    Tracks allele counts and metadata per chromosome during counting.
    """

    def __init__(self) -> None:
        self.ref_count: defaultdict[tuple[int, int], int] = defaultdict(int)
        self.alt_count: defaultdict[tuple[int, int], int] = defaultdict(int)
        self.other_count: defaultdict[tuple[int, int], int] = defaultdict(int)

        # Keep track of metadata

        # Number
        self.num_snps: defaultdict[str, int] = defaultdict(int)
        self.num_barcodes: defaultdict[str, int] = defaultdict(int)
        self.reads_counted: defaultdict[str, int] = defaultdict(int)

        # Reads that were not counted
        self.reads_skipped_no_barcode: defaultdict[str, int] = defaultdict(int)
        self.reads_skipped_barcode_no_index: defaultdict[str, int] = defaultdict(int)
        self.reads_skipped_prev_counted: defaultdict[str, int] = defaultdict(int)
        self.reads_skipped_no_sequence: defaultdict[str, int] = defaultdict(int)
        self.reads_skipped_no_aln_pos: defaultdict[str, int] = defaultdict(int)
        self.reads_skipped_seq_error: defaultdict[str, int] = defaultdict(int)

    def stats_to_df(self) -> pd.DataFrame:
        """Convert statistics to a pandas DataFrame."""
        stat_attributes = [
            "num_snps",
            "num_barcodes",
            "reads_counted",
            "reads_skipped_no_barcode",
            "reads_skipped_barcode_no_index",
            "reads_skipped_prev_counted",
            "reads_skipped_no_sequence",
            "reads_skipped_no_aln_pos",
            "reads_skipped_seq_error",
        ]

        stat_df = pd.DataFrame({key: getattr(self, key) for key in stat_attributes}).reset_index(
            names="chrom"
        )

        return stat_df


def make_count_matrix(
    bam_file: str,
    df: pl.DataFrame,
    bc_dict: dict[str, int],
    include_samples: list[str] | None = None,
    include_features: list[str] | None = None,
    use_rust: bool = True,
) -> ad.AnnData:
    """Create sparse count matrix from BAM and variant data.

    Parameters
    ----------
    bam_file : str
        Path to BAM file with cell barcodes.
    df : pl.DataFrame
        DataFrame with variant positions from intersection.
    bc_dict : dict[str, int]
        Mapping of cell barcodes to integer indices.
    include_samples : list[str] | None, optional
        Sample columns to include from variant data, by default None.
    include_features : list[str] | None, optional
        Feature columns to include, by default None.
    use_rust : bool, optional
        Use the Rust per-cell counter when available, by default True. Falls back
        to the pure-Python ``count_bc_snp_alleles`` if the extension is missing.

    Returns
    -------
    ad.AnnData
        AnnData object with count matrices in layers (ref, alt, other).
    """
    chrom_list = df.get_column("chrom").unique(maintain_order=True)

    # Add genotypes annotations
    snp_df_cols = ["chrom", "pos", "ref", "alt"]
    if include_samples is not None:
        sample_cols = list(include_samples)
        missing_sample_cols = [col for col in sample_cols if col not in df.columns]
        if missing_sample_cols and len(sample_cols) == 1 and "GT" in df.columns:
            sample_name = sample_cols[0]
            df = df.with_columns(pl.col("GT").alias(sample_name))
        snp_df_cols.extend(sample_cols)

    snp_df = df.select(snp_df_cols).unique(maintain_order=True).with_row_index()

    sc_counts = CountStatsSC()  # Class that holds total count data

    rust_path = use_rust and RUST_SC_AVAILABLE

    if rust_path:
        start = timeit.default_timer()
        count_bc_snp_alleles_rust(
            bam_file=bam_file,
            bc_dict=bc_dict,
            snp_df=snp_df,
            sc_counts=sc_counts,
        )
        logger.info(
            "Rust SC counter: %d SNPs across %d chromosomes in %.2f s",
            snp_df.shape[0],
            chrom_list.len(),
            timeit.default_timer() - start,
        )
    else:
        with AlignmentFile(bam_file, "rb") as bam:
            for chrom in chrom_list:
                chrom_df = snp_df.filter(pl.col("chrom") == chrom)

                start = timeit.default_timer()

                try:
                    count_bc_snp_alleles(
                        bam=bam,
                        bc_dict=bc_dict,
                        chrom=chrom,
                        snp_list=chrom_df.select(["index", "pos", "ref", "alt"]).iter_rows(),
                        sc_counts=sc_counts,
                    )
                except ValueError:
                    logger.warning("Skipping %s: Contig not found!", chrom)
                else:
                    logger.info(
                        "%s: Counted %d SNPs in %.2f seconds",
                        chrom,
                        chrom_df.height,
                        timeit.default_timer() - start,
                    )

    # Create sparse matrices
    matrix_shape = (snp_df.shape[0], len(bc_dict))
    sparse_ref = _sparse_from_counts(sc_counts.ref_count, matrix_shape)
    sparse_alt = _sparse_from_counts(sc_counts.alt_count, matrix_shape)
    sparse_other = _sparse_from_counts(sc_counts.other_count, matrix_shape)

    # Create anndata With total as X
    adata = ad.AnnData(
        X=sparse_ref + sparse_alt + sparse_other,
        layers={"ref": sparse_ref, "alt": sparse_alt, "other": sparse_other},
    )

    adata.obs = snp_df.to_pandas()
    adata.obs["ref_count"] = adata.layers["ref"].sum(axis=1).T.A1
    adata.obs["alt_count"] = adata.layers["alt"].sum(axis=1).T.A1

    # Add barcode names
    adata.var_names = bc_dict.keys()

    # Add genotypes to anndata
    if include_samples is not None:
        adata.uns["samples"] = include_samples

    feature_cols = include_features or []
    if not feature_cols and "region" in df.columns:
        feature_cols = ["region"]

    if feature_cols:
        feature_df = (
            df.join(snp_df, on=["chrom", "pos", "ref", "alt"], how="left")
            .select([*feature_cols, "index"])
            .unique(maintain_order=True)
        )
        if "region" not in feature_df.columns:
            feature_df = feature_df.with_columns(pl.col(feature_cols[0]).alias("region"))
        adata.uns["feature"] = feature_df.to_pandas()

        # region_snp_dict = dict(
        #     df.join(snp_df, on=["chrom", "pos", "ref", "alt"], how="left"
        #            ).group_by("region", maintain_order=True
        #                      ).agg("index").iter_rows()
        # )
        # adata.uns["region_snps"] = region_snp_dict

    # Write out count stats
    adata.uns["count_stats"] = sc_counts.stats_to_df()

    return adata


def count_bc_snp_alleles_rust(
    bam_file: str,
    bc_dict: dict[str, int],
    snp_df: pl.DataFrame,
    sc_counts: CountStatsSC,
    threads: int | None = None,
) -> None:
    """Rust-accelerated per-cell allele counting for ALL chromosomes in one call.

    Mirrors :func:`count_bc_snp_alleles`, filling the same
    ``sc_counts.ref_count`` / ``alt_count`` / ``other_count`` defaultdicts keyed
    by ``(snp_idx, bc_idx)`` so downstream sparse-matrix assembly is unchanged.

    Passing every region in a single call lets the Rust counter parallelize
    across chromosomes via Rayon (the previous per-chromosome calls never did).

    Parameters
    ----------
    bam_file : str
        Path to BAM file with cell barcodes.
    bc_dict : dict[str, int]
        Mapping of cell barcodes to indices.
    snp_df : pl.DataFrame
        DataFrame with columns ``index``, ``chrom``, ``pos``, ``ref``, ``alt``
        for every SNP across all chromosomes.
    sc_counts : CountStatsSC
        Statistics container to update with counts (COO-filled in place).
    threads : int | None, optional
        Number of Rayon worker threads (parallelizes across chromosomes). When
        ``None``, reads ``WASP2_RUST_THREADS``; if that is unset, defaults to
        ``min(8, os.cpu_count())`` so the accelerated path is fast out-of-the-box
        (single-thread is ~1.5x slower than pysam; cross-chromosome parallelism is
        the speedup). Set ``WASP2_RUST_THREADS=1`` to force the serial path.
    """
    cpu_default = min(8, os.cpu_count() or 1)
    rust_threads_env = os.environ.get("WASP2_RUST_THREADS") if threads is None else None
    try:
        rust_threads = (
            threads
            if threads is not None
            else (int(rust_threads_env) if rust_threads_env else cpu_default)
        )
    except ValueError:
        rust_threads = cpu_default
    rust_threads = max(1, rust_threads)

    # Build region tuples: (snp_idx, chrom, pos, ref, alt) for ALL chromosomes.
    regions = list(snp_df.select(["index", "chrom", "pos", "ref", "alt"]).iter_rows())

    # Call Rust counter; min_qual=0 matches WASP2 behavior (no quality filtering)
    counter = _RustSC(bam_file)
    coo = counter.count_bc_alleles(regions, bc_dict, 0, rust_threads)

    # Fill the COO results into the shared accumulators keyed by (snp_idx, bc_idx)
    for snp_idx, bc_idx, ref_count, alt_count, other_count in coo:
        key = (snp_idx, bc_idx)
        if ref_count:
            sc_counts.ref_count[key] += ref_count
        if alt_count:
            sc_counts.alt_count[key] += alt_count
        if other_count:
            sc_counts.other_count[key] += other_count

    # Bookkeeping consistent with count_bc_snp_alleles (per-chrom num_snps).
    for chrom, n in snp_df.group_by("chrom").len().iter_rows():
        sc_counts.num_snps[chrom] += n


def count_bc_snp_alleles(
    bam: AlignmentFile,
    bc_dict: dict[str, int],
    chrom: str,
    snp_list: Iterator[tuple[int, int, str, str]],
    sc_counts: CountStatsSC,
) -> None:
    """Count alleles at SNP positions for each cell barcode.

    Parameters
    ----------
    bam : AlignmentFile
        Open BAM file handle.
    bc_dict : dict[str, int]
        Mapping of cell barcodes to indices.
    chrom : str
        Chromosome to process.
    snp_list : Iterator[tuple[int, int, str, str]]
        Iterator of (index, pos, ref, alt) tuples.
    sc_counts : CountStatsSC
        Statistics container to update with counts.
    """
    read_set = set()  # Keep track of reads seen
    bc_set = set()

    for idx, pos, ref, alt in snp_list:
        for read in bam.fetch(chrom, pos - 1, pos):
            # If already counted allele or pair in read
            if read.query_name in read_set:
                sc_counts.reads_skipped_prev_counted[chrom] += 1
                continue

            # Check if there is a read barcode
            try:
                read_bc = str(read.get_tag("CB"))
            except KeyError:
                sc_counts.reads_skipped_no_barcode[chrom] += 1
                continue

            # If barcode not in dict
            if read_bc not in bc_dict:
                sc_counts.reads_skipped_barcode_no_index[chrom] += 1
                continue

            seq = read.query_sequence
            if seq is None:
                sc_counts.reads_skipped_no_sequence[chrom] += 1
                continue

            # Binary search for alignment position
            qpos = find_read_aln_pos(read, pos - 1)
            if qpos is None:
                sc_counts.reads_skipped_no_aln_pos[chrom] += 1
                continue

            try:
                if seq[qpos] == ref:
                    sc_counts.ref_count[(idx, bc_dict[read_bc])] += 1
                elif seq[qpos] == alt:
                    sc_counts.alt_count[(idx, bc_dict[read_bc])] += 1
                else:
                    sc_counts.other_count[(idx, bc_dict[read_bc])] += 1

            except (TypeError, IndexError) as e:
                # Narrow exception handling: only catch sequence access errors
                # Log the actual exception for debugging unexpected errors
                sc_counts.reads_skipped_seq_error[chrom] += 1
                logger.debug(
                    "Skipping read %s: sequence access error at %s:%d (qpos=%s): %s",
                    read.query_name,
                    chrom,
                    pos,
                    qpos,
                    e,
                )
                continue
            else:
                read_set.add(read.query_name)
                bc_set.add(read_bc)
                sc_counts.reads_counted[chrom] += 1

        sc_counts.num_snps[chrom] += 1  # Put here in case of error
    sc_counts.num_barcodes[chrom] = len(bc_set)  # Add unique barcodes
