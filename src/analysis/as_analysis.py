"""
Author: Aaron Ho
Python Version: 3.9
"""

# Default Python package Imports
import inspect
import time
import timeit
from collections.abc import Callable
from pathlib import Path
from typing import Any, Literal, cast

import numpy as np

# External package imports
import pandas as pd
from numpy.typing import NDArray
from scipy.optimize import OptimizeResult, minimize, minimize_scalar
from scipy.special import expit
from scipy.stats import betabinom, chi2, false_discovery_control

from wasp2.cli import create_spinner_progress, error, info, success

# =============================================================================
# BETA-BINOMIAL RHO PARAMETER BOUNDS (Issue #228)
# =============================================================================
# The beta-binomial parameterization uses alpha = mu * (1-rho) / rho, which
# causes division by zero when rho=0 and produces zero alpha/beta when rho=1.
# We clamp rho to (epsilon, 1-epsilon) to prevent numerical instability.
# =============================================================================

RHO_EPSILON: float = 1e-10


def clamp_rho(rho: float | NDArray[np.float64]) -> float | NDArray[np.float64]:
    """
    Clamp dispersion parameter rho to safe range (epsilon, 1-epsilon).

    The beta-binomial parameterization uses alpha = mu * (1-rho) / rho, which
    causes division by zero when rho=0 and produces zero alpha/beta when rho=1.
    This function prevents these boundary issues.

    Args:
        rho: Dispersion parameter (scalar or array), expected in [0, 1]

    Returns:
        Clamped rho in range (RHO_EPSILON, 1 - RHO_EPSILON)
    """
    return np.clip(rho, RHO_EPSILON, 1.0 - RHO_EPSILON)


def opt_linear(
    disp_params: NDArray[np.float64],
    ref_counts: NDArray[np.integer[Any]],
    n_array: NDArray[np.integer[Any]],
) -> float:
    """
    Optimize dispersion parameter weighted by N
    (Function called by optimizer)

    :param disp_params: Array of dispersion parameters [disp1, disp2]
    :param ref_counts: Array of reference allele counts
    :param n_array: Array of total counts (N)
    :return: Negative log-likelihood value
    """
    disp1, disp2 = disp_params

    exp_in = disp1 + (n_array * disp2)
    exp_in = np.select([exp_in > 10, exp_in < -10], [10, -10], default=exp_in)

    rho = clamp_rho(expit(exp_in))

    ll = -np.sum(
        betabinom.logpmf(ref_counts, n_array, (0.5 * (1 - rho) / rho), (0.5 * (1 - rho) / rho))
    )  # If alpha is beta

    return float(ll)


def opt_prob(
    in_prob: float | NDArray[np.float64],
    in_rho: float | NDArray[np.float64],
    k: int | NDArray[np.integer[Any]],
    n: int | NDArray[np.integer[Any]],
    log: bool = True,
) -> float | NDArray[np.float64]:
    """
    Optimize Probability value that maximizes imbalance likelihood.
    (Function called by optimizer)

    **CRITICAL FUNCTION** - Used by as_analysis_sc.py and compare_ai.py

    :param in_prob: Probability parameter (scalar or array)
    :param in_rho: Dispersion parameter (scalar or array)
    :param k: Reference allele count(s)
    :param n: Total count(s)
    :param log: If True, return negative log-likelihood; if False, return pmf
    :return: Negative log-likelihood (if log=True) or probability mass (if log=False)
    """
    prob = in_prob
    rho = clamp_rho(in_rho)  # Prevent division by zero at boundaries

    alpha = prob * (1 - rho) / rho
    beta = (1 - prob) * (1 - rho) / rho

    if log is True:
        ll = -1 * betabinom.logpmf(k, n, alpha, beta)
    else:
        ll = betabinom.pmf(k, n, alpha, beta)

    return cast(float | NDArray[np.float64], ll)


# updated phasing optimizer: currently used in single-cell analysis
# This version modifies prob arr outside of func
# GT phase should be with respect to first snp on first chrom
def opt_phased_new(
    prob: float,
    disp: float | NDArray[np.float64],
    ref_data: NDArray[np.integer[Any]],
    n_data: NDArray[np.integer[Any]],
    gt_data: NDArray[np.integer[Any]],
) -> float:
    """
    Optimize likelihood for phased data (updated version for single-cell analysis).

    **CRITICAL FUNCTION** - Used by as_analysis_sc.py and compare_ai.py

    :param prob: Probability parameter to optimize
    :param disp: Dispersion parameter (scalar or array)
    :param ref_data: Array of reference allele counts
    :param n_data: Array of total counts
    :param gt_data: Array of genotype phase information
    :return: Negative log-likelihood value
    """
    # phase and prob with respect to snp1 as ref
    phased_ll = opt_prob(np.abs(prob - gt_data), disp, ref_data, n_data)

    return float(np.sum(phased_ll))


# Updated unphasing optimizer using DP
def opt_unphased_dp(
    prob: float,
    disp: float | NDArray[np.float64],
    first_ref: NDArray[np.integer[Any]],
    first_n: NDArray[np.integer[Any]],
    phase_ref: NDArray[np.integer[Any]],
    phase_n: NDArray[np.integer[Any]],
) -> float:
    """
    Optimize likelihood while taking phase into account using dynamic programming.

    **CRITICAL FUNCTION** - Used by as_analysis_sc.py and compare_ai.py

    :param prob: Probability parameter to optimize
    :param disp: Dispersion parameter (scalar or array)
    :param first_ref: Reference count for first position (length 1 array)
    :param first_n: Total count for first position (length 1 array)
    :param phase_ref: Array of reference counts for subsequent positions
    :param phase_n: Array of total counts for subsequent positions
    :return: Negative log-likelihood value
    """
    # Get likelihood of first pos
    first_ll = opt_prob(prob, disp, first_ref[0], first_n[0])

    # Get likelihood witth regard to phasing of first pos
    phase1_like = opt_prob(prob, disp, phase_ref, phase_n, log=False)
    phase2_like = opt_prob(1 - prob, disp, phase_ref, phase_n, log=False)

    prev_like: float = 1.0
    # phase1_like and phase2_like are arrays when phase_ref/phase_n are arrays
    phase1_arr = cast(NDArray[np.float64], phase1_like)
    phase2_arr = cast(NDArray[np.float64], phase2_like)
    for p1, p2 in zip(phase1_arr, phase2_arr):
        p1_combined_like = prev_like * p1
        p2_combined_like = prev_like * p2
        prev_like = float((0.5 * p1_combined_like) + (0.5 * p2_combined_like))

    return float(first_ll + -np.log(prev_like))


def parse_opt(
    df: pd.DataFrame, disp: float | NDArray[np.float64] | None = None, phased: bool = False
) -> tuple[float, float]:
    """
    Optimize necessary data when running model

    :param df: Dataframe with allele counts
    :param disp: pre-computed dispersion parameter, defaults to None
    :param phased: Whether data is phased
    :return: Tuple of (alt_ll, mu) - likelihood of alternate model and imbalance proportion
    """
    snp_count = df.shape[0]

    # Create array used for AI analysis
    ref_array = df["ref_count"].to_numpy()
    n_array = df["N"].to_numpy()

    # In the case that we do use linear model with disp per N
    if disp is None:
        disp = df["disp"].to_numpy()

    res: OptimizeResult
    if snp_count > 1:
        # If data is phased
        if phased:
            # Use known phasing info
            gt_array = df["GT"].to_numpy()

            # First pos with respect to ref
            if gt_array[0] > 0:
                gt_array = 1 - gt_array

            res = minimize_scalar(
                opt_phased_new,
                args=(disp, ref_array, n_array, gt_array),
                method="bounded",
                bounds=(0, 1),
            )

        else:
            # Use unphased algorithm for subsequent phases
            first_ref = ref_array[:1]
            first_n = n_array[:1]

            phase_ref = ref_array[1:]
            phase_n = n_array[1:]

            res = minimize_scalar(
                opt_unphased_dp,
                args=(disp, first_ref, first_n, phase_ref, phase_n),
                method="bounded",
                bounds=(0, 1),
            )

    else:
        # Single site optimize
        res = minimize_scalar(
            opt_prob, args=(disp, ref_array[0], n_array[0]), method="bounded", bounds=(0, 1)
        )

    # Get res data
    mu: float = res["x"]
    alt_ll: float = -1 * res["fun"]

    return alt_ll, mu


def single_model(df: pd.DataFrame, region_col: str, phased: bool = False) -> pd.DataFrame:
    """
    Find allelic imbalance using normal beta-binomial model

    :param df: Dataframe with allele counts
    :param region_col: Name of column to group by
    :param phased: Whether data is phased
    :return: Dataframe with imbalance likelihood
    """
    info("Running analysis with single dispersion model")

    def opt_disp(rho: float, ref_data: NDArray[Any], n_data: NDArray[Any]) -> float:
        """Negative log-likelihood for dispersion optimization (rho clamped)."""
        rho_safe = float(clamp_rho(rho))
        return float(
            -np.sum(
                betabinom.logpmf(
                    ref_data,
                    n_data,
                    (0.5 * (1 - rho_safe) / rho_safe),
                    (0.5 * (1 - rho_safe) / rho_safe),
                )
            )
        )

    ref_array = df["ref_count"].to_numpy()
    n_array = df["N"].to_numpy()

    disp_start = timeit.default_timer()

    with create_spinner_progress() as progress:
        progress.add_task("Optimizing dispersion parameter", total=None)
        disp: float = float(
            clamp_rho(
                minimize_scalar(
                    opt_disp, args=(ref_array, n_array), method="bounded", bounds=(0, 1)
                )["x"]
            )
        )

    success(f"Optimized dispersion parameter ({timeit.default_timer() - disp_start:.2f}s)")

    group_df = df.groupby(region_col, sort=False)
    include_groups_supported = "include_groups" in inspect.signature(group_df.apply).parameters
    apply_kwargs = {"include_groups": False} if include_groups_supported else {}

    ll_start = timeit.default_timer()

    with create_spinner_progress() as progress:
        progress.add_task("Optimizing imbalance likelihood", total=None)
        null_test = group_df.apply(
            lambda x: np.sum(
                betabinom.logpmf(
                    x["ref_count"].to_numpy(),
                    x["N"].to_numpy(),
                    (0.5 * (1 - disp) / disp),
                    (0.5 * (1 - disp) / disp),
                )
            ),
            **apply_kwargs,
        )

        # Optimize Alt
        alt_test = group_df.apply(lambda x: parse_opt(x, disp, phased=phased), **apply_kwargs)
        alt_df = pd.DataFrame(alt_test.tolist(), columns=["alt_ll", "mu"], index=alt_test.index)

    success(f"Optimized imbalance likelihood ({timeit.default_timer() - ll_start:.2f}s)")

    ll_df = pd.concat([null_test, alt_df], axis=1).reset_index()
    ll_df.columns = [region_col, "null_ll", "alt_ll", "mu"]

    ll_df["lrt"] = -2 * (ll_df["null_ll"] - ll_df["alt_ll"])
    ll_df["pval"] = chi2.sf(ll_df["lrt"], 1)

    return ll_df


def linear_model(df: pd.DataFrame, region_col: str, phased: bool = False) -> pd.DataFrame:
    """
    Find allelic imbalance using linear allelic imbalance model,
    weighting imbalance linear with N counts

    :param df: Dataframe with allele counts
    :param region_col: Name of column to group by
    :param phased: Whether data is phased
    :return: Dataframe with imbalance likelihood
    """
    info("Running analysis with linear dispersion model")
    in_data = df[["ref_count", "N"]].to_numpy().T

    disp_start = time.time()

    with create_spinner_progress() as progress:
        progress.add_task("Optimizing dispersion parameters", total=None)
        res: OptimizeResult = minimize(
            opt_linear, x0=(0, 0), method="Nelder-Mead", args=(in_data[0], in_data[1])
        )
    disp1: float
    disp2: float
    disp1, disp2 = res["x"]
    df["disp"] = clamp_rho(expit(disp1 + (in_data[1] * disp2)))

    success(f"Optimized dispersion parameters ({time.time() - disp_start:.2f}s)")

    # Group by region
    group_df = df.groupby(region_col, sort=False)
    include_groups_supported = "include_groups" in inspect.signature(group_df.apply).parameters
    apply_kwargs = {"include_groups": False} if include_groups_supported else {}

    # Get null test
    ll_start = time.time()

    with create_spinner_progress() as progress:
        progress.add_task("Optimizing imbalance likelihood", total=None)
        null_test = group_df.apply(
            lambda x: np.sum(
                betabinom.logpmf(
                    x["ref_count"].to_numpy(),
                    x["N"].to_numpy(),
                    (0.5 * (1 - x["disp"].to_numpy()) / x["disp"].to_numpy()),
                    (0.5 * (1 - x["disp"].to_numpy()) / x["disp"].to_numpy()),
                )
            ),
            **apply_kwargs,
        )

        # Optimize Alt
        alt_test = group_df.apply(lambda x: parse_opt(x), **apply_kwargs)
        alt_df = pd.DataFrame(alt_test.tolist(), columns=["alt_ll", "mu"], index=alt_test.index)

    success(f"Optimized imbalance likelihood ({time.time() - ll_start:.2f}s)")

    ll_df = pd.concat([null_test, alt_df], axis=1).reset_index()
    ll_df.columns = [region_col, "null_ll", "alt_ll", "mu"]

    ll_df["lrt"] = -2 * (ll_df["null_ll"] - ll_df["alt_ll"])
    ll_df["pval"] = chi2.sf(ll_df["lrt"], 1)

    return ll_df


def get_imbalance(
    in_data: pd.DataFrame | str | Path,
    min_count: int = 10,
    pseudocount: int = 1,
    method: Literal["single", "linear"] = "single",
    phased: bool = False,
    region_col: str | None = None,
    groupby: str | None = None,
) -> pd.DataFrame:
    """
    Process input data and method for finding allelic imbalance.

    **CRITICAL FUNCTION** - Main analysis entry point used by run_analysis.py

    :param in_data: Dataframe with allele counts or filepath to TSV file
    :param min_count: minimum allele count for analysis
    :param pseudocount: pseudocount to add to allele counts
    :param method: analysis method ("single" or "linear")
    :param phased: whether to use phased genotype information
    :param region_col: column name to group variants by (e.g., gene, peak)
    :param groupby: alternative grouping column (overrides region_col if provided)
    :return: DataFrame with imbalance statistics per region
    """
    model_dict: dict[str, Callable[[pd.DataFrame, str, bool], pd.DataFrame]] = {
        "single": single_model,
        "linear": linear_model,
    }

    # If preparsed dataframe or filepath
    if isinstance(in_data, pd.DataFrame):
        df = in_data
    else:
        df = pd.read_csv(
            in_data,
            sep="\t",
            dtype={
                "chrom": "category",
                "pos": np.uint32,
                "ref": "category",
                "alt": "category",
                "ref_count": np.uint16,
                "alt_count": np.uint16,
                "other_count": np.uint16,
            },
        )

    # If no region_col measure imbalance per variant
    if region_col is None:
        region_col = "variant"
        groupby = None  # no parent

        df[region_col] = df["chrom"].astype("string") + "_" + df["pos"].astype("string")

    # Process pseudocount values and filter data by min
    df[["ref_count", "alt_count"]] += pseudocount
    df["N"] = df["ref_count"] + df["alt_count"]
    df = df.loc[df["N"].ge(min_count + (2 * pseudocount)), :]

    # Get unique values based on group
    if groupby is not None:
        region_col = groupby

    keep_cols = ["chrom", "pos", "ref_count", "alt_count", "N", region_col]

    # Check validity of phasing info
    if phased:
        # Check if GT are actually phased - use error() so always shown (even in quiet mode)
        if "GT" not in df.columns:
            error("Genotypes not found: Switching to unphased model")
            phased = False
        elif len(df["GT"].unique()) <= 1:
            error(f"All genotypes {df['GT'].unique()}: Switching to unphased model")
            phased = False
        elif not any(i in ["1|0", "0|1"] for i in df["GT"].unique()):
            error(f"Expected GT as 0|1 and 1|0 but found: {df['GT'].unique()}")
            error("Switching to unphased model")
            phased = False
        else:
            # GT is indeed phased
            df["GT"] = df["GT"].str.split("|", n=1).str[0].astype(dtype=np.uint8)
            keep_cols.append("GT")

    df = df[keep_cols].drop_duplicates()

    p_df = model_dict[method](df, region_col, phased)  # Perform analysis

    # remove pseudocount
    df[["ref_count", "alt_count"]] -= pseudocount
    df["N"] -= pseudocount * 2

    snp_counts = pd.DataFrame(df[region_col].value_counts(sort=False)).reset_index()
    snp_counts.columns = [region_col, "snp_count"]

    count_alleles = (
        df[[region_col, "ref_count", "alt_count", "N"]].groupby(region_col, sort=False).sum()
    )

    merge_df = pd.merge(snp_counts, p_df, how="left", on=region_col)

    as_df = pd.merge(count_alleles, merge_df, how="left", on=region_col)
    as_df["fdr_pval"] = false_discovery_control(as_df["pval"], method="bh")

    return as_df
