"""Allelic imbalance analysis pipeline.

Main entry point for running the beta-binomial allelic imbalance analysis
using the Rust-accelerated backend.
"""

from __future__ import annotations

import logging
from csv import reader
from pathlib import Path
from typing import Literal

import pandas as pd

# Rust analysis (required; no Python fallback)
try:
    from wasp2_rust import analyze_imbalance as rust_analyze_imbalance
except ImportError:
    rust_analyze_imbalance = None

logger = logging.getLogger(__name__)

_COUNT_METADATA_COLUMNS = {
    "chrom",
    "pos0",
    "start",
    "pos",
    "end",
    "stop",
    "ref",
    "alt",
    "GT",
    "genotype",
    "ref_count",
    "alt_count",
    "other_count",
    "total_count",
    "N",
    "sample",
    "donor_id",
    "vcf_sample",
}


class WaspAnalysisData:
    """Container for allelic imbalance analysis configuration.

    Attributes
    ----------
    count_file : str | Path
        Path to the count TSV file.
    region_col : str | None
        Column name for grouping variants by region.
    groupby : str | None
        Column name for additional grouping (e.g., parent gene).
    out_file : str
        Output file path for results.
    phased : bool
        Whether to use phased genotype information.
    model : Literal["single", "linear"]
        Dispersion model type.
    min_count : int
        Minimum total allele count threshold.
    pseudocount : int
        Pseudocount to add to allele counts.
    """

    def __init__(
        self,
        count_file: str | Path,
        min_count: int | None = None,
        pseudocount: int | None = None,
        phased: bool | None = None,
        model: str | None = None,
        out_file: str | None = None,
        region_col: str | None = None,
        groupby: str | None = None,
        per_variant: bool = False,
    ) -> None:
        # Per-variant (SNV-solo) mode cannot also name a region column: per-variant tests each
        # SNV independently, whereas region_col groups variants by that column.
        if per_variant and region_col is not None:
            raise ValueError(
                "per_variant cannot be combined with region_col: per-variant tests each SNV "
                "independently, while region_col groups by that column. Choose one."
            )
        self.per_variant: bool = per_variant

        # User input data
        self.count_file = count_file
        self.region_col = region_col
        self.groupby = groupby
        self.out_file = out_file

        self.phased: bool = bool(phased)

        # Default to single dispersion model
        if model == "linear":
            self.model: Literal["single", "linear"] = "linear"
        else:
            self.model = "single"

        # Default min count of 10, pseudocount of 1
        self.min_count: int = 10 if min_count is None else min_count
        self.pseudocount: int = 1 if pseudocount is None else pseudocount

        # Read header only for validation
        with open(self.count_file) as f:
            count_cols = next(reader(f, delimiter="\t"))

        # Feature columns are the fields not owned by the variant/count schema. Identifying
        # them by name keeps pos0, sample, and donor metadata from becoming grouping keys.
        feature_cols = [col for col in count_cols if col not in _COUNT_METADATA_COLUMNS]

        # Skip auto-detection in per-variant mode so the Rust backend groups by chrom/pos.
        if self.region_col is None and not self.per_variant:
            if feature_cols:
                self.region_col = feature_cols[0]

        # By default group by feature rather than parent?
        if self.groupby is not None:
            # If denoting to group by feature
            if (self.region_col is None) or (self.groupby == self.region_col):
                self.groupby = None

            elif self.groupby not in feature_cols:
                logger.warning("%s not found in columns %s", self.groupby, count_cols)
                self.groupby = None

        # Create default outfile
        if self.out_file is None:
            self.out_file = str(Path.cwd() / "ai_results.tsv")


def run_ai_analysis(
    count_file: str | Path,
    min_count: int | None = None,
    pseudocount: int | None = None,
    phased: bool | None = None,
    model: str | None = None,
    out_file: str | None = None,
    region_col: str | None = None,
    groupby: str | None = None,
    per_variant: bool = False,
) -> None:
    """Run allelic imbalance analysis pipeline.

    Parameters
    ----------
    count_file : str | Path
        Path to TSV file with allele counts.
    min_count : int | None, optional
        Minimum total count threshold, by default 10.
    pseudocount : int | None, optional
        Pseudocount to add, by default 1.
    phased : bool | None, optional
        Use phased genotype information, by default False.
    model : str | None, optional
        Dispersion model ('single' or 'linear'), by default 'single'.
    out_file : str | None, optional
        Output file path, by default 'ai_results.tsv'.
    region_col : str | None, optional
        Column name for grouping variants.
    groupby : str | None, optional
        Additional grouping column (not supported by the Rust backend; raises if set).
    per_variant : bool, optional
        Test each SNV independently (per-variant) instead of grouping by a region column.
        Forces per-variant even when a region column is present. Default False.

    Raises
    ------
    RuntimeError
        If Rust analysis extension is not available.
    """
    # Fail closed: --groupby is not supported by the Rust backend. It only re-keys the grouping
    # column (region_col = groupby), so the identical result is obtained with --region_col <parent>.
    # Erroring prevents silently returning feature-level results when parent-level was requested.
    if groupby is not None:
        raise RuntimeError(
            "--groupby (parent-level grouping) is not supported by the Rust analysis backend. "
            "Since groupby only re-keys the grouping column, group by the parent column directly "
            "with --region_col <parent_column> instead."
        )

    # Store analysis data and params
    ai_files = WaspAnalysisData(
        count_file,
        min_count=min_count,
        pseudocount=pseudocount,
        phased=phased,
        model=model,
        out_file=out_file,
        region_col=region_col,
        groupby=groupby,
        per_variant=per_variant,
    )

    # Run analysis pipeline (Rust only)
    if rust_analyze_imbalance is None:
        raise RuntimeError(
            "Rust analysis extension not available. Build it with "
            "`maturin develop --release` in the WASP2 env."
        )

    results = rust_analyze_imbalance(
        str(ai_files.count_file),
        min_count=ai_files.min_count,
        pseudocount=ai_files.pseudocount,
        method=ai_files.model,
        phased=ai_files.phased,
        region_col=ai_files.region_col,
    )
    ai_df = pd.DataFrame(results)

    if "fdr_pval" in ai_df.columns:
        ai_df = ai_df.sort_values(by="fdr_pval", ascending=True)

    # Write results
    ai_df.to_csv(ai_files.out_file, sep="\t", header=True, index=False)
