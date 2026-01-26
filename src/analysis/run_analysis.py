"""Allelic imbalance analysis pipeline.

Main entry point for running the beta-binomial allelic imbalance analysis
using the Rust-accelerated backend.
"""

from __future__ import annotations

from csv import reader
from pathlib import Path
from typing import Literal

import pandas as pd

# Rust analysis (required; no Python fallback)
try:
    from wasp2_rust import analyze_imbalance as rust_analyze_imbalance
except ImportError:
    rust_analyze_imbalance = None


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
    ) -> None:
        # User input data
        self.count_file = count_file
        self.region_col = region_col
        self.groupby = groupby  # group by region or parent?
        self.out_file = out_file

        # TODO parse vcf for phased instead of default unphased
        self.phased: bool = bool(phased)

        # Default to single dispersion model
        if model in {"single", "linear"}:
            self.model: Literal["single", "linear"] = model  # type: ignore[assignment]
        else:
            self.model = "single"

        # Default min count of 10, pseudocount of 1
        self.min_count: int = 10 if min_count is None else min_count
        self.pseudocount: int = 1 if pseudocount is None else pseudocount

        # Read header only for validation
        with open(self.count_file) as f:
            count_cols = next(reader(f, delimiter="\t"))

        # 7 columns at minimum, 10 at maximum
        # 3required : chr, pos, ref, alt
        # 3 optional: <GT>, <region>, <parent>
        # 3 required: ref_count, alt_count, other_count
        # [chr, pos, ref, alt, <GT>, <region>, <parent>, ref_c, alt_c, other_c]

        if "GT" in count_cols:
            min_cols = 8
            region_idx = 5
        else:
            min_cols = 7
            region_idx = 4

        # Check regions
        if self.region_col is None:
            if len(count_cols) > min_cols:
                self.region_col = count_cols[region_idx]

        # By default group by feature rather than parent?
        if self.groupby is not None:
            # If denoting to group by feature
            if (self.region_col is None) or (self.groupby == self.region_col):
                self.groupby = None

            elif (len(count_cols) > (min_cols + 1)) and self.groupby in {
                count_cols[region_idx + 1],
                "Parent",
                "parent",
            }:
                self.groupby = count_cols[region_idx + 1]  # Set group
            else:
                # Maybe throw error instead
                print(f"{self.groupby} not found in columns \n{count_cols}")
                self.groupby = None

        # Create default outfile
        if self.out_file is None:
            self.out_file = str(Path.cwd() / "ai_results.tsv")  # do this after


def run_ai_analysis(
    count_file: str | Path,
    min_count: int | None = None,
    pseudocount: int | None = None,
    phased: bool | None = None,
    model: str | None = None,
    out_file: str | None = None,
    region_col: str | None = None,
    groupby: str | None = None,
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
        Additional grouping column.

    Raises
    ------
    RuntimeError
        If Rust analysis extension is not available.
    """
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
    )
    ai_df = pd.DataFrame(results)

    # Maybe give option to sort or not sort by pval
    if "fdr_pval" in ai_df.columns:
        ai_df = ai_df.sort_values(by="fdr_pval", ascending=True)

    # Write results
    ai_df.to_csv(ai_files.out_file, sep="\t", header=True, index=False)
