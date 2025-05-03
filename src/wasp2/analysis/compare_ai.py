"""
Module for comparing allelic imbalance between groups.

This module provides functions to compute allelic imbalance statistics for both single-group 
and pairwise comparisons using shared SNPs (or alternative SNP selections) from an AnnData object.
It includes functions to choose the appropriate likelihood function, combine likelihoods, and perform 
likelihood ratio tests for comparing imbalance between groups.

Note:
    Several commented-out code lines are preserved for potential alternative implementations or future extensions.
"""

import sys
import warnings
from pathlib import Path
from collections import namedtuple
from itertools import combinations
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

# import polars as pl
# import anndata as ad

from scipy.stats import betabinom, chi2, false_discovery_control
from scipy.optimize import minimize_scalar

# Local imports
from wasp2.analysis.as_analysis import opt_prob, opt_unphased_dp, opt_phased_new, bh_correction
# from run_analysis_sc import WaspAnalysisSC, process_adata_inputs


def get_imbalance_func(ref_count: np.ndarray, 
                       n_count: np.ndarray, 
                       phase_array: Optional[np.ndarray] = None
                      ) -> Tuple[Callable[..., float], Tuple[Any, ...]]:
    """
    Determine the appropriate likelihood function and its arguments for imbalance estimation.

    Based on the number of SNPs and the availability of phasing information, this function selects one of:
    
    - A single-SNP optimization using `opt_prob`
    - An unphased multi-SNP optimization using `opt_unphased_dp`
    - A phased multi-SNP optimization using `opt_phased_new`
    
    Parameters
    ----------
    ref_count : np.ndarray
        Array of reference allele counts.
    n_count : np.ndarray
        Array of total allele counts.
    phase_array : np.ndarray, optional
        Array of phasing information. If provided, indicates that the data is phased.
    
    Returns
    -------
    tuple
        A tuple containing:
            - like_func: The selected likelihood function.
            - like_func_args: A tuple of arguments for the likelihood function.
    """
    if len(ref_count) == 1:
        # Parse single opt
        like_func = opt_prob
        # This excludes disp since we always use disp
        like_func_args = (ref_count[0], n_count[0])
    elif phase_array is None:
        # Do unphased
        like_func = opt_unphased_dp
        like_func_args = (ref_count[:1], n_count[:1],
                          ref_count[1:], n_count[1:])
    else:
        # Do phased
        like_func = opt_phased_new
        like_func_args = (ref_count, n_count, phase_array)
    
    return like_func, like_func_args


def opt_combined_imbalance(prob: float, 
                           disp: float,
                           like_func1: Callable[..., float], 
                           like_func1_args: Tuple[Any, ...],
                           like_func2: Callable[..., float], 
                           like_func2_args: Tuple[Any, ...]
                          ) -> float:
    """
    Compute the combined likelihood for imbalance estimation across two groups.

    This function returns the sum of the likelihoods from two groups by calling the respective 
    likelihood functions with the provided arguments.
    
    Parameters
    ----------
    prob : float
        The imbalance probability parameter.
    disp : float
        Dispersion parameter.
    like_func1 : callable
        Likelihood function for the first group.
    like_func1_args : tuple
        Arguments for the first group's likelihood function.
    like_func2 : callable
        Likelihood function for the second group.
    like_func2_args : tuple
        Arguments for the second group's likelihood function.
    
    Returns
    -------
    float
        The combined likelihood value.
    """
    return (like_func1(prob, disp, *like_func1_args) +
            like_func2(prob, disp, *like_func2_args))


def get_compared_imbalance(adata: Any,  # ideally: ad.AnnData
                           min_count: int = 10,
                           pseudocount: int = 1,
                           phased: bool = False,
                           sample: Optional[str] = None,
                           groups: Optional[List[str]] = None
                          ) -> Dict[Tuple[str, str], pd.DataFrame]:
    """
    Compare allelic imbalance between groups using shared SNPs.

    This function computes pairwise comparisons of allelic imbalance between groups present in the AnnData 
    object. It first calculates count data and dispersion for the entire dataset, then for each pair of groups 
    it identifies shared SNPs, aggregates allele counts per region, and performs likelihood ratio tests to compare 
    imbalance estimates.
    
    Parameters
    ----------
    adata : anndata.AnnData
        AnnData object containing layers "ref" and "alt", along with metadata in `adata.var` and `adata.obs`.
    min_count : int, optional
        Minimum total allele count required for a region (default is 10).
    pseudocount : int, optional
        Pseudocount added to allele counts to avoid division by zero (default is 1).
    phased : bool, optional
        Whether genotype data is phased. If True, specialized processing is applied.
    sample : str, optional
        Column name in `adata.obs` containing phasing information, used if `phased` is True.
    groups : list, optional
        List of groups to compare. If None, all unique groups from `adata.var["group"]` are used.
    
    Returns
    -------
    dict
        Dictionary where keys are tuples (group1, group2) and values are pandas DataFrames containing 
        comparison results with columns such as region, number of SNPs, combined imbalance estimates, null 
        and alternative log-likelihoods, p-value, and FDR-corrected p-value.
    """
    # Failsafe in case preparse somehow misses these
    if sample is None:
        phased = False

    # Should I be comparing all combinations by default? Seems like a lot.
    if groups is None:
        groups = list(adata.var["group"].dropna().unique())
        print("Comparing all combinations of available groups")
    elif len(groups) == 1:
        raise ValueError("Please provide 2 or more groups to compare.")
    
    # Process initial minimums for whole data dispersion
    region_cutoff = min_count + (2 * pseudocount)
    snp_cutoff = (2 * pseudocount)
    
    ref_counts = adata.layers["ref"].sum(axis=1, dtype=np.uint16).T.A1 + pseudocount
    alt_counts = adata.layers["alt"].sum(axis=1, dtype=np.uint16).T.A1 + pseudocount
    n_counts = ref_counts + alt_counts

    # Calculate dispersion across dataset
    opt_disp = lambda rho, ref_data, n_data: -np.sum(
        betabinom.logpmf(ref_data, n_data, (0.5 * (1 - rho) / rho), (0.5 * (1 - rho) / rho))
    )
    
    disp = minimize_scalar(opt_disp, args=(ref_counts, n_counts), method="bounded", bounds=(0, 1))["x"]
    
    if phased:
        gt_array = adata.obs[sample].str.split("|", n=1).str[0].to_numpy(dtype=np.uint8)
    else:
        gt_array = None

    # Process counts on a per-group basis to avoid recalculating
    group_dict: Dict[str, Any] = {}
    # group_data = namedtuple("group_data", ["ref_counts", "n_counts", "phase_data", "region_snp_dict"]) 
    group_data = namedtuple("group_data", ["ref_counts", "n_counts", "region_snp_df"])
    
    for group_name in groups:
        # Subset by group
        adata_sub = adata[:, adata.var["group"] == group_name]

        # Create count data per group; should I do pseudocount now or later?
        ref_counts_group = adata_sub.layers["ref"].sum(axis=1, dtype=np.uint16).T.A1 + pseudocount
        alt_counts_group = adata_sub.layers["alt"].sum(axis=1, dtype=np.uint16).T.A1 + pseudocount
        n_counts_group = ref_counts_group + alt_counts_group

        min_idx_group = np.where(n_counts_group >= snp_cutoff)  # Get indices where counts were found
    
        ref_counts_group_filt, n_counts_group_filt = ref_counts_group[min_idx_group], n_counts_group[min_idx_group]

        if phased:
            phase_array = adata.obs.iloc[min_idx_group][sample].str.split("|", n=1).str[0].to_numpy(dtype=np.uint8)
        else:
            phase_array = None

        # Create region_snp_dict for each group
        idx_df = pd.DataFrame({"index": min_idx_group[0]}, dtype=np.uint32).reset_index(names="filt_index")

        region_snp_dict = adata.uns["feature"].merge(
            idx_df, on="index")[["region", "filt_index"]].groupby("region", sort=False).agg(tuple)["filt_index"].to_dict()

        group_dict[group_name] = group_data(ref_counts_group_filt, n_counts_group_filt, region_snp_dict)
    
    # Create group combinations and process shared SNPs
    group_combos = list(combinations(group_dict.keys(), r=2))
    
    df_dict: Dict[Tuple[str, str], pd.DataFrame] = {}
    for group1, group2 in group_combos:
        df = compare_imbalance_between_groups_diff_snps(
            disp,
            *group_dict[group1],
            *group_dict[group2]
        )
        if df.empty:
            print(f"Skipping {group1}-{group2} comparison. No shared regions.")
        else:
            df_dict[(group1, group2)] = df

    return df_dict


def compare_imbalance_between_groups_diff_snps(
    disp: float,
    ref_counts1: np.ndarray,
    n_counts1: np.ndarray,
    phase_array1: Optional[np.ndarray],
    region_snp_dict1: Dict[str, Tuple[int, ...]],
    ref_counts2: np.ndarray,
    n_counts2: np.ndarray,
    phase_array2: Optional[np.ndarray],
    region_snp_dict2: Dict[str, Tuple[int, ...]]
) -> pd.DataFrame:
    """
    Compare allelic imbalance between two groups using different SNPs (non-shared regions).

    This function performs pairwise comparisons of imbalance estimates for shared regions between two groups,
    using separate SNP sets defined in `region_snp_dict1` and `region_snp_dict2`. It selects the appropriate 
    likelihood functions based on the available data and performs optimization under both the null and alternative 
    hypotheses to compute a likelihood ratio test.

    Parameters
    ----------
    disp : float
        Dispersion parameter estimated from the dataset.
    ref_counts1 : np.ndarray
        Array of reference allele counts for group 1.
    n_counts1 : np.ndarray
        Array of total allele counts for group 1.
    phase_array1 : np.ndarray, optional
        Array of phasing information for group 1.
    region_snp_dict1 : dict
        Dictionary mapping region identifiers to tuples of shared SNP indices for group 1.
    ref_counts2 : np.ndarray
        Array of reference allele counts for group 2.
    n_counts2 : np.ndarray
        Array of total allele counts for group 2.
    phase_array2 : np.ndarray, optional
        Array of phasing information for group 2.
    region_snp_dict2 : dict
        Dictionary mapping region identifiers to tuples of shared SNP indices for group 2.
    
    Returns
    -------
    pd.DataFrame
        DataFrame containing comparison results for each region with columns:
            - region: The region identifier.
            - num_snps_group1: Number of SNPs in group 1 for the region.
            - num_snps_group2: Number of SNPs in group 2 for the region.
            - combined_mu: Combined imbalance estimate under the null hypothesis.
            - mu1: Imbalance estimate for group 1.
            - mu2: Imbalance estimate for group 2.
            - null_ll: Log-likelihood under the null hypothesis.
            - alt_ll: Log-likelihood under the alternative hypothesis.
            - pval: p-value from the likelihood ratio test.
        
        FDR correction is applied to the p-values.
    """
    group_results: List[Tuple[Any, ...]] = []  # Store imbalance results
    
    # Compare allelic imbalance difference per region
    for region, snp_list in region_snp_dict1.items():
        # Get per-region SNPs and counts for both groups
        region_ref1 = ref_counts1[snp_list,]
        region_n1 = n_counts1[snp_list,]

        region_ref2 = ref_counts2[snp_list,]
        region_n2 = n_counts2[snp_list,]

        if phase_array1 is not None and phase_array2 is not None:
            phased = True
            region_phasing1 = phase_array1[snp_list,]
            region_phasing2 = phase_array2[snp_list,]
        else:
            phased = False
            region_phasing1, region_phasing2 = None, None
        
        # Process which model to use for likelihood estimation for each group
        like_func1, like_func_inputs1 = get_imbalance_func(region_ref1, region_n1, phase_array=region_phasing1)
        like_func2, like_func_inputs2 = get_imbalance_func(region_ref2, region_n2, phase_array=region_phasing2)

        # Null Hypothesis: Imbalance is the same between groups
        null_res = minimize_scalar(
            opt_combined_imbalance,
            args=(disp, like_func1, like_func_inputs1, like_func2, like_func_inputs2),
            method="bounded",
            bounds=(0, 1)
        )
        combined_mu = null_res["x"]
        null_ll = -1 * null_res["fun"]

        # Alternative Hypothesis: Imbalance is different between groups
        alt_res1 = minimize_scalar(
            like_func1,
            args=(disp, *like_func_inputs1),
            method="bounded",
            bounds=(0, 1)
        )
        alt_res2 = minimize_scalar(
            like_func2,
            args=(disp, *like_func_inputs2),
            method="bounded",
            bounds=(0, 1)
        )
        alt_mu1 = alt_res1["x"]
        alt_mu2 = alt_res2["x"]
        alt_ll = -1 * (alt_res1["fun"] + alt_res2["fun"])

        # Likelihood ratio test
        lrt = -2 * (null_ll - alt_ll)
        pval = chi2.sf(lrt, 1)

        # Add data to output list
        group_results.append(
            (region, len(snp_list), len(snp_list), combined_mu, alt_mu1, alt_mu2, null_ll, alt_ll, pval)
        )
    
    # Create allelic imbalance DataFrame
    df = pd.DataFrame(
        group_results,
        columns=["region", "num_snps_group1", "num_snps_group2", "combined_mu",
                 "mu1", "mu2", "null_ll", "alt_ll", "pval"]
    )
    
    # fdr correction
    df["fdr_pval"] = false_discovery_control(df["pval"], method="bh")
    
    return df


# The following function is missing from the original code.
def get_compared_imbalance_diff_snps(
    adata: Any,  # ideally: ad.AnnData
    min_count: int = 10,
    pseudocount: int = 1,
    phased: bool = False,
    sample: Optional[str] = None,
    groups: Optional[List[str]] = None
) -> Dict[Tuple[str, str], pd.DataFrame]:
    """
    Compare allelic imbalance between groups using different SNPs (non-shared regions).

    This alternative version compares imbalance estimates between groups without requiring that SNPs be shared 
    between regions. It filters data based on a minimum count cutoff, computes dispersion, and then for each group 
    calculates allele counts and creates region dictionaries. Pairwise comparisons are then performed on the shared 
    regions between groups.
    
    Parameters
    ----------
    adata : anndata.AnnData
        AnnData object with layers "ref" and "alt", and group information in `adata.var["group"]`.
    min_count : int, optional
        Minimum total allele count required for a region (default is 10).
    pseudocount : int, optional
        Pseudocount added to allele counts (default is 1).
    phased : bool, optional
        Whether genotype data is phased (default is False).
    sample : str, optional
        Column name in `adata.obs` with phasing information, used if `phased` is True.
    groups : list, optional
        List of groups to compare. If None, all unique groups in `adata.var["group"]` are used.
    
    Returns
    -------
    dict
        Dictionary where keys are tuples (group1, group2) and values are DataFrames containing imbalance 
        comparison results. If no shared regions are found for a pair, that comparison is skipped.
    """
    # Failsafe in case preparse somehow misses these
    if sample is None:
        phased = False

    # If groups not provided, use all unique groups from adata.var["group"]
    if groups is None:
        groups = list(adata.var["group"].dropna().unique())
        print("Comparing all combinations of available groups")
    elif len(groups) == 1:
        raise ValueError("Please provide 2 or more groups to compare.")

    # Process counts for the entire dataset
    cutoff = min_count + (2 * pseudocount)
    ref_counts = adata.layers["ref"].sum(axis=1, dtype=np.uint16).T.A1 + pseudocount
    alt_counts = adata.layers["alt"].sum(axis=1, dtype=np.uint16).T.A1 + pseudocount
    n_counts = ref_counts + alt_counts
    min_idx = np.where(n_counts >= cutoff)

    ref_counts_filt, n_counts_filt = ref_counts[min_idx], n_counts[min_idx]

    # Calculate dispersion across dataset
    opt_disp = lambda rho, ref_data, n_data: -np.sum(
        betabinom.logpmf(ref_data, n_data, (0.5 * (1 - rho) / rho), (0.5 * (1 - rho) / rho))
    )
    disp = minimize_scalar(opt_disp, args=(ref_counts_filt, n_counts_filt), method="bounded", bounds=(0, 1))["x"]

    # Process counts on a per-group basis to avoid recalculating
    group_dict: Dict[str, Any] = {}
    group_data = namedtuple("group_data", ["ref_counts", "n_counts", "phase_data", "region_snp_dict"])

    for group_name in groups:
        adata_sub = adata[:, adata.var["group"] == group_name]

        ref_counts_group = adata_sub.layers["ref"].sum(axis=1, dtype=np.uint16).T.A1 + pseudocount
        alt_counts_group = adata_sub.layers["alt"].sum(axis=1, dtype=np.uint16).T.A1 + pseudocount
        n_counts_group = ref_counts_group + alt_counts_group

        min_idx_group = np.where(n_counts_group >= cutoff)
        ref_counts_group_filt, n_counts_group_filt = ref_counts_group[min_idx_group], n_counts_group[min_idx_group]

        if phased:
            phase_array = adata.obs.iloc[min_idx_group][sample].str.split("|", n=1).str[0].to_numpy(dtype=np.uint8)
        else:
            phase_array = None

        # Create region_snp_dict for each group
        idx_df = pd.DataFrame({"index": min_idx_group[0]}, dtype=np.uint32).reset_index(names="filt_index")
        region_snp_dict = adata.uns["feature"].merge(
            idx_df, on="index")[["region", "filt_index"]].groupby("region", sort=False).agg(tuple)["filt_index"].to_dict()

        group_dict[group_name] = group_data(ref_counts_group_filt, n_counts_group_filt, phase_array, region_snp_dict)
    
    # Create group combinations and process shared SNPs
    group_combos = list(combinations(group_dict.keys(), r=2))
    df_dict: Dict[Tuple[str, str], pd.DataFrame] = {}
    
    for group1, group2 in group_combos:
        df = compare_imbalance_between_groups_diff_snps(
            disp,
            *group_dict[group1],
            *group_dict[group2]
        )
        if df.empty:
            print(f"Skipping {group1}-{group2} comparison. No shared regions.")
        else:
            df_dict[(group1, group2)] = df

    return df_dict


def compare_imbalance_between_groups_diff_snps(
    disp: float,
    ref_counts1: np.ndarray,
    n_counts1: np.ndarray,
    phase_array1: Optional[np.ndarray],
    region_snp_dict1: Dict[str, Tuple[int, ...]],
    ref_counts2: np.ndarray,
    n_counts2: np.ndarray,
    phase_array2: Optional[np.ndarray],
    region_snp_dict2: Dict[str, Tuple[int, ...]]
) -> pd.DataFrame:
    """
    Compare allelic imbalance between two groups using different SNPs (non-shared regions).

    This function performs pairwise comparisons of imbalance estimates for shared regions between two groups,
    using separate SNP sets defined in `region_snp_dict1` and `region_snp_dict2`. It selects the appropriate 
    likelihood functions based on the available data and performs optimization under both the null and alternative 
    hypotheses to compute a likelihood ratio test.

    Parameters
    ----------
    disp : float
        Dispersion parameter estimated from the dataset.
    ref_counts1 : np.ndarray
        Array of reference allele counts for group 1.
    n_counts1 : np.ndarray
        Array of total allele counts for group 1.
    phase_array1 : np.ndarray, optional
        Array of phasing information for group 1.
    region_snp_dict1 : dict
        Dictionary mapping region identifiers to tuples of shared SNP indices for group 1.
    ref_counts2 : np.ndarray
        Array of reference allele counts for group 2.
    n_counts2 : np.ndarray
        Array of total allele counts for group 2.
    phase_array2 : np.ndarray, optional
        Array of phasing information for group 2.
    region_snp_dict2 : dict
        Dictionary mapping region identifiers to tuples of shared SNP indices for group 2.
    
    Returns
    -------
    pd.DataFrame
        DataFrame containing comparison results for each region with columns:
            - region: The region identifier.
            - num_snps_group1: Number of SNPs in group 1 for the region.
            - num_snps_group2: Number of SNPs in group 2 for the region.
            - combined_mu: Combined imbalance estimate under the null hypothesis.
            - mu1: Imbalance estimate for group 1.
            - mu2: Imbalance estimate for group 2.
            - null_ll: Log-likelihood under the null hypothesis.
            - alt_ll: Log-likelihood under the alternative hypothesis.
            - pval: p-value from the likelihood ratio test.
        
        FDR correction is applied to the p-values.
    """
    group_results: List[Tuple[Any, ...]] = []  # Store imbalance results
    
    # Compare allelic imbalance difference per region
    for region, snp_list in region_snp_dict1.items():
        # Get per-region SNPs and counts for both groups
        region_ref1 = ref_counts1[snp_list,]
        region_n1 = n_counts1[snp_list,]

        region_ref2 = ref_counts2[snp_list,]
        region_n2 = n_counts2[snp_list,]

        if phase_array1 is not None and phase_array2 is not None:
            phased = True
            region_phasing1 = phase_array1[snp_list,]
            region_phasing2 = phase_array2[snp_list,]
        else:
            phased = False
            region_phasing1, region_phasing2 = None, None
        
        # Process which model to use for likelihood estimation for each group
        like_func1, like_func_inputs1 = get_imbalance_func(region_ref1, region_n1, phase_array=region_phasing1)
        like_func2, like_func_inputs2 = get_imbalance_func(region_ref2, region_n2, phase_array=region_phasing2)

        # Null Hypothesis: Imbalance is the same between groups
        null_res = minimize_scalar(
            opt_combined_imbalance,
            args=(disp, like_func1, like_func_inputs1, like_func2, like_func_inputs2),
            method="bounded",
            bounds=(0, 1)
        )
        combined_mu = null_res["x"]
        null_ll = -1 * null_res["fun"]

        # Alternative Hypothesis: Imbalance is different between groups
        alt_res1 = minimize_scalar(
            like_func1,
            args=(disp, *like_func_inputs1),
            method="bounded",
            bounds=(0, 1)
        )
        alt_res2 = minimize_scalar(
            like_func2,
            args=(disp, *like_func_inputs2),
            method="bounded",
            bounds=(0, 1)
        )
        alt_mu1 = alt_res1["x"]
        alt_mu2 = alt_res2["x"]
        alt_ll = -1 * (alt_res1["fun"] + alt_res2["fun"])

        # Likelihood ratio test
        lrt = -2 * (null_ll - alt_ll)
        pval = chi2.sf(lrt, 1)

        # Add data to output list
        group_results.append(
            (region, len(snp_list), len(snp_list), combined_mu, alt_mu1, alt_mu2, null_ll, alt_ll, pval)
        )
    
    # Create allelic imbalance DataFrame
    df = pd.DataFrame(
        group_results,
        columns=["region", "num_snps_group1", "num_snps_group2", "combined_mu",
                 "mu1", "mu2", "null_ll", "alt_ll", "pval"]
    )
    
    # FDR correction
    df["fdr_pval"] = false_discovery_control(df["pval"], method="bh")
    
    return df
