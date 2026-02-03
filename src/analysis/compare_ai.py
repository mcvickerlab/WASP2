import logging
from collections import namedtuple
from collections.abc import Callable
from itertools import combinations
from typing import Any

import numpy as np
import pandas as pd

# AnnData for single-cell analysis
from anndata import AnnData
from numpy.typing import NDArray
from scipy.optimize import OptimizeResult, minimize_scalar
from scipy.stats import betabinom, chi2, false_discovery_control

# Local imports
from .as_analysis import clamp_rho, opt_phased_new, opt_prob, opt_unphased_dp

logger = logging.getLogger(__name__)


# Use these functions to figure out how to optimize per group
def get_imbalance_func(
    ref_count: NDArray[np.integer[Any]],
    n_count: NDArray[np.integer[Any]],
    phase_array: NDArray[np.integer[Any]] | None = None,
) -> tuple[Callable[..., Any], tuple[Any, ...]]:
    """
    Determine which imbalance function to use based on data characteristics.

    :param ref_count: Array of reference allele counts
    :param n_count: Array of total counts
    :param phase_array: Optional phasing information array
    :return: Tuple of (likelihood function, function arguments)
    """
    like_func: Callable[..., Any]
    like_func_args: tuple[Any, ...]

    if len(ref_count) == 1:
        # Parse single opt
        like_func = opt_prob

        # This excludes disp since we always use disp
        like_func_args = (ref_count[0], n_count[0])
    elif phase_array is None:
        # Do unphased
        like_func = opt_unphased_dp
        like_func_args = (
            ref_count[:1],
            n_count[:1],
            ref_count[1:],
            n_count[1:],
        )
    else:
        # Do phased
        like_func = opt_phased_new
        like_func_args = (ref_count, n_count, phase_array)

    return like_func, like_func_args


def opt_combined_imbalance(
    prob: float,
    disp: float,
    like_func1: Callable[..., float],
    like_func1_args: tuple[Any, ...],
    like_func2: Callable[..., float],
    like_func2_args: tuple[Any, ...],
) -> float:
    """
    Optimize combined imbalance likelihood for two groups.

    :param prob: Probability parameter
    :param disp: Dispersion parameter
    :param like_func1: Likelihood function for group 1
    :param like_func1_args: Arguments for group 1 likelihood function
    :param like_func2: Likelihood function for group 2
    :param like_func2_args: Arguments for group 2 likelihood function
    :return: Combined negative log-likelihood
    """
    return like_func1(prob, disp, *like_func1_args) + like_func2(prob, disp, *like_func2_args)


# Current version that uses shared snps
def get_compared_imbalance(
    adata: AnnData,
    min_count: int = 10,
    pseudocount: int = 1,
    phased: bool = False,
    sample: str | None = None,
    groups: list[str] | None = None,
) -> dict[tuple[str, str], pd.DataFrame]:
    """
    Compare allelic imbalance between groups using shared SNPs.

    :param adata: AnnData object containing SNP count data
    :param min_count: Minimum allele count threshold
    :param pseudocount: Pseudocount to add to avoid zero counts
    :param phased: Whether to use phased analysis
    :param sample: Sample column name for phasing information
    :param groups: List of groups to compare (if None, compare all)
    :return: Dict mapping (group1, group2) tuples to comparison DataFrames
    """
    # Failsafe in case preparse somehow misses these
    if sample is None:
        phased = False

    # Should I be comparing all combos by default??? Seems like a lot
    if groups is None:
        groups = list(adata.var["group"].dropna().unique())
        logger.info("Comparing all combinations of available groups")
    elif len(groups) == 1:
        raise ValueError("Please provide 2 or more groups to compare.")

    # Process initial minimums for whole data dispersion
    min_count + (2 * pseudocount)
    snp_cutoff: int = 2 * pseudocount

    ref_counts: NDArray[np.uint16] = (
        adata.layers["ref"].sum(axis=1, dtype=np.uint16).T.A1 + pseudocount
    )
    alt_counts: NDArray[np.uint16] = (
        adata.layers["alt"].sum(axis=1, dtype=np.uint16).T.A1 + pseudocount
    )
    n_counts: NDArray[np.uint16] = ref_counts + alt_counts

    # Calculate dispersion across dataset
    def opt_disp(
        rho: float, ref_data: NDArray[np.uint16], n_data: NDArray[np.uint16]
    ) -> float:
        rho_safe = float(clamp_rho(rho))  # Prevent division by zero (Issue #228)
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

    disp: float = float(
        clamp_rho(
            minimize_scalar(
                opt_disp, args=(ref_counts, n_counts), method="bounded", bounds=(0, 1)
            )["x"]
        )
    )

    if phased:
        gt_array: NDArray[np.uint8] | None = (
            adata.obs[sample].str.split("|", n=1).str[0].to_numpy(dtype=np.uint8)
        )
    else:
        gt_array = None

    # process counts on a per group basis to avoid recalculating
    group_dict: dict[str, Any] = {}
    # group_data = namedtuple("group_data", ["ref_counts", "n_counts", "phase_data", "region_snp_dict"]) # Maybe include the gt_array instead of min_idx
    group_data = namedtuple("group_data", ["ref_counts", "n_counts", "region_snp_df"])

    for group_name in groups:
        # Subset by group
        adata_sub = adata[:, adata.var["group"] == group_name]

        # Create count data per group, should i do pseudocount now or later?
        ref_counts_group = adata_sub.layers["ref"].sum(axis=1, dtype=np.uint16).T.A1 + pseudocount
        alt_counts_group = adata_sub.layers["alt"].sum(axis=1, dtype=np.uint16).T.A1 + pseudocount
        n_counts_group = ref_counts_group + alt_counts_group

        nonzero_idx = np.where(
            n_counts_group > snp_cutoff
        )  # Get indices where no counts were found

        if nonzero_idx[0].size == 0:
            logger.warning("Skipping %s: no SNP counts found", group_name)
            continue

        # Remove snps with 0 counts from regions
        idx_df = pd.DataFrame({"index": nonzero_idx[0]}, dtype=np.uint32).reset_index(
            names="filt_index"
        )
        region_idx_df = adata.uns["feature"].merge(idx_df, on="index")

        # Check total allele counts/N per region
        region_n_df = region_idx_df.merge(
            pd.DataFrame(n_counts_group, columns=["N"]).reset_index(names="index"), on="index"
        )

        group_dict[group_name] = group_data(ref_counts_group, n_counts_group, region_n_df)

    # Create group combinations and process shared snps
    group_combos: list[tuple[str, str]] = list(combinations(group_dict.keys(), r=2))

    df_dict: dict[tuple[str, str], pd.DataFrame] = {}
    for group1, group2 in group_combos:
        # Get relevant counts and nonzero snps
        ref_counts1, n_counts1, region_snp_df1 = group_dict[group1]
        ref_counts2, n_counts2, region_snp_df2 = group_dict[group2]

        # Get shared snps -> get regions that meet cutoff
        shared_df = region_snp_df1[["region", "index", "N"]].merge(
            region_snp_df2[["index", "N"]], on="index", suffixes=("1", "2")
        )

        # Take into account pseudocounts added to total N
        region_agg_df = shared_df.groupby("region", sort=False).agg(
            snp_idx=("index", tuple),
            num_snps=("index", "size"),
            N1=("N1", np.sum),
            N2=("N2", np.sum),
        )

        region_agg_df["region_cutoff"] = (region_agg_df["num_snps"] * snp_cutoff) + min_count

        # Find regions where N is satisfied for both
        # region_agg_df = shared_df.groupby("region", sort=False).agg(
        #     snp_idx=("index", tuple), N1=("N1", np.sum), N2=("N2", np.sum)
        # )

        # Per group snp_dict
        region_snp_dict = region_agg_df.loc[
            (
                (region_agg_df["N1"] >= region_agg_df["region_cutoff"])
                & (region_agg_df["N2"] >= region_agg_df["region_cutoff"])
            ),
            "snp_idx",
        ].to_dict()

        # region_snp_dict = region_agg_df.loc[
        #     (region_agg_df["N1"] >= region_cutoff) & (region_agg_df["N2"] >= region_cutoff),
        #     "snp_idx"].to_dict()

        if not region_snp_dict:
            logger.warning(
                "Skipping %s-%s comparison: no shared regions with allele counts >= %d",
                group1, group2, min_count,
            )

            continue

        # This sub function name kinda long...find better name maybe?
        df = compare_imbalance_between_groups(
            disp, ref_counts1, n_counts1, ref_counts2, n_counts2, region_snp_dict, gt_array
        )

        # Using a tuple as key
        df_dict[(group1, group2)] = df

    return df_dict


def compare_imbalance_between_groups(
    disp: float,
    ref_counts1: NDArray[np.uint16],
    n_counts1: NDArray[np.uint16],
    ref_counts2: NDArray[np.uint16],
    n_counts2: NDArray[np.uint16],
    region_snp_dict: dict[str, tuple[int, ...]],
    gt_array: NDArray[np.uint8] | None = None,
) -> pd.DataFrame:
    """
    Compare allelic imbalance between two groups for shared regions.

    :param disp: Dispersion parameter
    :param ref_counts1: Reference allele counts for group 1
    :param n_counts1: Total counts for group 1
    :param ref_counts2: Reference allele counts for group 2
    :param n_counts2: Total counts for group 2
    :param region_snp_dict: Dict mapping region names to SNP index tuples
    :param gt_array: Optional genotype/phasing array
    :return: DataFrame with comparison statistics and p-values
    """
    # Helper func called by get_compared_imbalance()

    group_results: list[
        tuple[str, int, float, float, float, float, float, float]
    ] = []  # Store imbalance results

    # Compare allelic imbalance difference per region
    for region, snp_list in region_snp_dict.items():
        # Get per region snps and counts
        region_ref1 = ref_counts1[snp_list,]
        region_n1 = n_counts1[snp_list,]

        region_ref2 = ref_counts2[snp_list,]
        region_n2 = n_counts2[snp_list,]

        # Process which model we'll use to process likelihood per group
        like_func: Callable[..., Any]
        like_func_args1: tuple[Any, ...]
        like_func_args2: tuple[Any, ...]

        if len(snp_list) == 1:
            # Parse single opt
            like_func = opt_prob

            # This excludes disp since we always use disp
            like_func_args1 = (region_ref1[0], region_n1[0])
            like_func_args2 = (region_ref2[0], region_n2[0])

        elif gt_array is None:
            # Do unphased
            like_func = opt_unphased_dp

            like_func_args1 = (
                region_ref1[:1],
                region_n1[:1],
                region_ref1[1:],
                region_n1[1:],
            )

            like_func_args2 = (
                region_ref2[:1],
                region_n2[:1],
                region_ref2[1:],
                region_n2[1:],
            )

        else:
            # Do phased

            # Get phasing info
            region_gt = gt_array[snp_list,]

            # Make sure phase with respect to first snp ref
            if region_gt[0] > 0:
                region_gt = 1 - region_gt

            like_func = opt_phased_new

            like_func_args1 = (region_ref1, region_n1, region_gt)
            like_func_args2 = (region_ref2, region_n2, region_gt)

        # Null Hypothesis: Imbalance is the same
        null_res: OptimizeResult = minimize_scalar(
            opt_combined_imbalance,
            args=(disp, like_func, like_func_args1, like_func, like_func_args2),
            method="bounded",
            bounds=(0, 1),
        )

        combined_mu: float = null_res["x"]
        null_ll: float = -1 * null_res["fun"]

        # Alt Hypothesis: Imbalance is different between groups
        alt_res1: OptimizeResult = minimize_scalar(
            like_func, args=(disp, *like_func_args1), method="bounded", bounds=(0, 1)
        )

        alt_res2: OptimizeResult = minimize_scalar(
            like_func, args=(disp, *like_func_args2), method="bounded", bounds=(0, 1)
        )

        # Get separate mu
        alt_mu1: float = alt_res1["x"]
        alt_mu2: float = alt_res2["x"]

        # get Alternative likelihood
        alt_ll1: float = alt_res1["fun"]
        alt_ll2: float = alt_res2["fun"]

        alt_ll: float = -1 * (alt_ll1 + alt_ll2)

        # Log ratio ttest
        lrt: float = -2 * (null_ll - alt_ll)
        pval: float = chi2.sf(lrt, 1)

        # Add data to output list

        # How should i format this, lots of possible outputs
        group_results.append(
            (region, len(snp_list), combined_mu, alt_mu1, alt_mu2, null_ll, alt_ll, pval)
        )

    # Create allelic imbalance df

    # Polars implementation might be more performant
    df: pd.DataFrame = pd.DataFrame(
        group_results,
        columns=["region", "num_snps", "combined_mu", "mu1", "mu2", "null_ll", "alt_ll", "pval"],
    )

    # fdr correction
    df["fdr_pval"] = false_discovery_control(df["pval"], method="bh")

    return df


# THIS IS A V0 VERSION THAT DIDN'T USE SHARED SNPS BETWEEN REGIONS
# COULD BE USEFUL AS AN OPTION POSSIBLY
def get_compared_imbalance_diff_snps(
    adata: AnnData,
    min_count: int = 10,
    pseudocount: int = 1,
    phased: bool = False,
    sample: str | None = None,
    groups: list[str] | None = None,
) -> dict[tuple[str, str], pd.DataFrame]:
    """
    Compare allelic imbalance between groups (V0 version without shared SNPs).

    :param adata: AnnData object containing SNP count data
    :param min_count: Minimum allele count threshold
    :param pseudocount: Pseudocount to add to avoid zero counts
    :param phased: Whether to use phased analysis
    :param sample: Sample column name for phasing information
    :param groups: List of groups to compare (if None, compare all)
    :return: Dict mapping (group1, group2) tuples to comparison DataFrames
    """
    # Failsafe in case preparse somehow misses these
    if sample is None:
        phased = False

    # Should I be comparing all combos by default??? Seems like a lot
    if groups is None:
        groups = list(adata.var["group"].dropna().unique())
        logger.info("Comparing all combinations of available groups")
    elif len(groups) == 1:
        raise ValueError("Please provide 2 or more groups to compare.")

    # Process initial minimums for whole data dispersion
    cutoff: int = min_count + (2 * pseudocount)

    ref_counts: NDArray[np.uint16] = (
        adata.layers["ref"].sum(axis=1, dtype=np.uint16).T.A1 + pseudocount
    )
    alt_counts: NDArray[np.uint16] = (
        adata.layers["alt"].sum(axis=1, dtype=np.uint16).T.A1 + pseudocount
    )

    n_counts: NDArray[np.uint16] = ref_counts + alt_counts
    min_idx: tuple[NDArray[np.intp], ...] = np.where(
        n_counts >= cutoff
    )  # Get indices for min_count

    ref_counts_filt: NDArray[np.uint16]
    n_counts_filt: NDArray[np.uint16]
    ref_counts_filt, n_counts_filt = ref_counts[min_idx], n_counts[min_idx]

    # Calculate dispersion across dataset
    def opt_disp_filt(
        rho: float, ref_data: NDArray[np.uint16], n_data: NDArray[np.uint16]
    ) -> float:
        rho_safe = float(clamp_rho(rho))  # Prevent division by zero (Issue #228)
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

    disp: float = float(
        clamp_rho(
            minimize_scalar(
                opt_disp_filt,
                args=(ref_counts_filt, n_counts_filt),
                method="bounded",
                bounds=(0, 1),
            )["x"]
        )
    )

    # process counts on a per group basis to avoid recalculating
    group_dict: dict[str, Any] = {}
    group_data = namedtuple(
        "group_data", ["ref_counts", "n_counts", "phase_data", "region_snp_dict"]
    )  # Maybe include the gt_array instead of min_idx

    for group_name in groups:
        # Subset by group
        adata_sub = adata[:, adata.var["group"] == group_name]

        # Create count data per group, should i do pseudocount now or later?
        ref_counts_group = adata_sub.layers["ref"].sum(axis=1, dtype=np.uint16).T.A1 + pseudocount
        alt_counts_group = adata_sub.layers["alt"].sum(axis=1, dtype=np.uint16).T.A1 + pseudocount
        n_counts_group = ref_counts_group + alt_counts_group

        min_idx_group = np.where(n_counts_group >= cutoff)  # Get indices for min_count

        ref_counts_group_filt, n_counts_group_filt = (
            ref_counts_group[min_idx_group],
            n_counts_group[min_idx_group],
        )

        if phased:
            phase_array = (
                adata.obs.iloc[min_idx_group][sample]
                .str.split("|", n=1)
                .str[0]
                .to_numpy(dtype=np.uint8)
            )
        else:
            phase_array = None

        # Create region_snp_dict but for each group
        idx_df = pd.DataFrame({"index": min_idx_group[0]}, dtype=np.uint32).reset_index(
            names="filt_index"
        )

        region_snp_dict = (
            adata.uns["feature"]
            .merge(idx_df, on="index")[["region", "filt_index"]]
            .groupby("region", sort=False)
            .agg(tuple)["filt_index"]
            .to_dict()
        )

        group_dict[group_name] = group_data(
            ref_counts_group_filt, n_counts_group_filt, phase_array, region_snp_dict
        )

    # Create group combinations and process shared snps
    group_combos: list[tuple[str, str]] = list(combinations(group_dict.keys(), r=2))

    df_dict: dict[tuple[str, str], pd.DataFrame] = {}
    for group1, group2 in group_combos:
        # Might be smart to create a cache to prevent repeating calculations
        # This sub function name kinda long...find better name maybe?
        df = compare_imbalance_between_groups_diff_snps(
            disp, *group_dict[group1], *group_dict[group2]
        )

        if df.empty:
            logger.warning("Skipping %s - %s comparison: no shared regions", group1, group2)
        else:
            # Using a tuple as key
            df_dict[(group1, group2)] = df

    return df_dict


def compare_imbalance_between_groups_diff_snps(
    disp: float,
    ref_counts1: NDArray[np.uint16],
    n_counts1: NDArray[np.uint16],
    phase_array1: NDArray[np.uint8] | None,
    region_snp_dict1: dict[str, tuple[int, ...]],
    ref_counts2: NDArray[np.uint16],
    n_counts2: NDArray[np.uint16],
    phase_array2: NDArray[np.uint8] | None,
    region_snp_dict2: dict[str, tuple[int, ...]],
) -> pd.DataFrame:
    """
    Compare allelic imbalance between two groups with different SNPs per region.

    :param disp: Dispersion parameter
    :param ref_counts1: Reference allele counts for group 1
    :param n_counts1: Total counts for group 1
    :param phase_array1: Optional phasing array for group 1
    :param region_snp_dict1: Dict mapping region names to SNP index tuples for group 1
    :param ref_counts2: Reference allele counts for group 2
    :param n_counts2: Total counts for group 2
    :param phase_array2: Optional phasing array for group 2
    :param region_snp_dict2: Dict mapping region names to SNP index tuples for group 2
    :return: DataFrame with comparison statistics and p-values
    """
    # These values are unpacked versions of named tuple
    # Helper func called by get_compared_imbalance()

    # Check if phasing info available
    phased: bool = (phase_array1 is not None) and (phase_array2 is not None)

    # Get shared regions
    shared_regions: list[str] = [i for i in region_snp_dict1 if i in region_snp_dict2]

    group_results: list[
        tuple[str, int, int, float, float, float, float, float, float]
    ] = []  # Store imbalance results

    # Compare allelic imbalance difference per region
    for region in shared_regions:
        # Get per region snps and counts
        snp_list1 = region_snp_dict1[region]
        region_ref1 = ref_counts1[snp_list1,]
        region_n1 = n_counts1[snp_list1,]

        snp_list2 = region_snp_dict2[region]
        region_ref2 = ref_counts2[snp_list2,]
        region_n2 = n_counts2[snp_list2,]

        if phased:
            assert phase_array1 is not None and phase_array2 is not None
            region_phasing1 = phase_array1[snp_list1,]
            region_phasing2 = phase_array2[snp_list2,]
        else:
            region_phasing1, region_phasing2 = None, None

        # Process which model we'll use to process likelihood per group
        like_func1, like_func_inputs1 = get_imbalance_func(
            region_ref1, region_n1, phase_array=region_phasing1
        )

        like_func2, like_func_inputs2 = get_imbalance_func(
            region_ref2, region_n2, phase_array=region_phasing2
        )

        # Null Hypothesis: Imbalance is the same
        null_res: OptimizeResult = minimize_scalar(
            opt_combined_imbalance,
            args=(disp, like_func1, like_func_inputs1, like_func2, like_func_inputs2),
            method="bounded",
            bounds=(0, 1),
        )

        combined_mu: float = null_res["x"]
        null_ll: float = -1 * null_res["fun"]

        # Alt Hypothesis: Imbalance is different between groups
        alt_res1: OptimizeResult = minimize_scalar(
            like_func1, args=(disp, *like_func_inputs1), method="bounded", bounds=(0, 1)
        )

        alt_res2: OptimizeResult = minimize_scalar(
            like_func2, args=(disp, *like_func_inputs2), method="bounded", bounds=(0, 1)
        )

        # Get separate mu
        alt_mu1: float = alt_res1["x"]
        alt_mu2: float = alt_res2["x"]

        # get Alternative likelihood
        alt_ll: float = -1 * (alt_res1["fun"] + alt_res2["fun"])

        # Log ratio ttest
        lrt: float = -2 * (null_ll - alt_ll)
        pval: float = chi2.sf(lrt, 1)

        # Add data to output list

        # How should i format this, lots of possible outputs
        group_results.append(
            (
                region,
                len(snp_list1),
                len(snp_list2),
                combined_mu,
                alt_mu1,
                alt_mu2,
                null_ll,
                alt_ll,
                pval,
            )
        )

    # Create allelic imbalance df

    # Polars implementation might be more performant
    df: pd.DataFrame = pd.DataFrame(
        group_results,
        columns=[
            "region",
            "num_snps_group1",
            "num_snps_group2",
            "combined_mu",
            "mu1",
            "mu2",
            "null_ll",
            "alt_ll",
            "pval",
        ],
    )

    # fdr correction
    df["fdr_pval"] = false_discovery_control(df["pval"], method="bh")

    return df
