"""
This module provides functions for quality control and allelic imbalance analysis 
on AnnData count datasets. The functions include filtering out SNP outliers 
(adata_count_qc) and computing imbalance statistics across regions and groups 
(get_imbalance_sc and get_imbalance_per_group).

The analysis leverages beta‐binomial models, likelihood ratio tests, and (optionally) 
phased genotype data.
"""

import sys
import warnings
from pathlib import Path
from typing import Optional, Union, Dict, List, Any

import numpy as np
import pandas as pd

import anndata as ad

from scipy.stats import betabinom, chi2, zscore, false_discovery_control
from scipy.optimize import minimize_scalar

# Local imports
from wasp2.analysis.as_analysis import opt_prob, opt_unphased_dp, opt_phased_new, bh_correction


def adata_count_qc(adata: ad.AnnData,
                   z_cutoff: Optional[float] = None,
                   gt_error: Optional[Any] = None) -> ad.AnnData:
    """
    Perform quality control and prefiltering on an AnnData count dataset.
    
    This function filters out SNP outliers based on a z-score cutoff and, optionally,
    prepares the data for genotype error filtering. If neither `z_cutoff` nor 
    `gt_error` is provided, the original AnnData object is returned unmodified.
    
    Parameters
    ----------
    adata : anndata.AnnData
        An AnnData object containing count data, with observations (obs) and features (uns).
    z_cutoff : float, optional
        Z-score threshold to identify and filter out SNP outliers based on the sum of 
        reference and alternative allele counts.
    gt_error : any, optional
        Placeholder for future genotype error filtering options.
    
    Returns
    -------
    anndata.AnnData
        The filtered AnnData object after outlier removal and reindexing.
    
    Notes
    -----
    - The function computes the total allele count ("N") for each SNP.
    - SNP outliers are identified where the absolute z-score of "N" exceeds the provided cutoff.
    - Outlier regions (regions with one or more outlier SNPs) are removed from the dataset.
    - Genotype error filtering (gt_error) is not yet implemented.
    
    Examples
    --------
    >>> adata_filtered = adata_count_qc(adata, z_cutoff=3)
    """
    
    # No need to prefilt
    if z_cutoff is None and gt_error is None:
        return adata
    
    # Filt outliers
    if z_cutoff is not None:
        snp_outliers = adata.obs[["index", "ref_count", "alt_count"]].copy()
        snp_outliers["N"] = snp_outliers["ref_count"] + snp_outliers["alt_count"]
        snp_outliers = snp_outliers[np.abs(zscore(snp_outliers["N"])) > z_cutoff]  # At least 3
        
        # Todo: add option if there aren't any features
        # Get regions containing 1 or more outlier snps
        snp_outliers = snp_outliers.merge(adata.uns["feature"], on="index", how="left")
        
        outlier_regions = adata.uns["feature"].loc[
            adata.uns["feature"]["region"].isin(snp_outliers["region"].unique()), :
        ]
        
        # Remove outlier regions from adata
        adata = adata[~adata.obs["index"].isin(outlier_regions["index"]), :].copy()
        adata.obs = adata.obs.reset_index(drop=True)  # update index
        
        # Update valid regions and snps
        adata.uns["feature"] = adata.uns["feature"].merge(
            adata.obs[["index"]].reset_index(names="filt_index"),
            on="index"
        )[["region", "filt_index"]].rename(columns={"filt_index": "index"})
        
        adata.obs["index"] = adata.obs.index  # Replace index column
    
    # TODO add options to identify and filter GT errors
    if gt_error is not None:
        pass
    
    return adata


def get_imbalance_sc(adata: ad.AnnData,
                     min_count: int = 10,
                     pseudocount: int = 1,
                     phased: bool = False,
                     sample: Optional[str] = None,
                     groups: Optional[List[str]] = None) -> Dict[str, pd.DataFrame]:
    """
    Compute allelic imbalance statistics for each group in the dataset.
    
    This function processes an AnnData object to calculate allelic imbalance per region 
    across specified groups. It computes total allele counts, estimates a dispersion 
    parameter using a beta-binomial model, and calculates imbalance statistics using 
    likelihood ratio tests.
    
    Parameters
    ----------
    adata : anndata.AnnData
        An AnnData object containing count data with layers "ref" and "alt" and 
        associated metadata.
    min_count : int, optional
        Minimum total allele count required for a region to be considered (default is 10).
    pseudocount : int, optional
        Pseudocount added to allele counts to avoid division by zero (default is 1).
    phased : bool, optional
        Indicates whether genotype data is phased. If True, specialized processing is applied.
    sample : str, optional
        Column name in `adata.obs` that contains genotype information when data is phased.
    groups : list, optional
        List of group names to process. If None, all unique groups from `adata.var["group"]`
        (excluding NaN) are used.
    
    Returns
    -------
    dict
        A dictionary where keys are group names and values are pandas DataFrames containing 
        imbalance statistics for each region. Each DataFrame includes columns for region, number 
        of SNPs, estimated imbalance (mu), null log-likelihood, alternative log-likelihood, and p-value.
    
    Notes
    -----
    - A beta-binomial dispersion parameter is estimated across the dataset using a minimization 
      of the negative log-likelihood.
    - For each group, SNPs with insufficient counts are filtered out.
    - When phased genotype data is provided, the function uses a specialized optimization routine.
    - If no regions meet the minimum allele count criteria for a group, that group is skipped.
    
    See Also
    --------
    get_imbalance_per_group : Computes imbalance per region for a single group.
    
    Examples
    --------
    >>> imbalance_results = get_imbalance_sc(adata, min_count=10, pseudocount=1)
    """
    
    # Need to preparse input using process_adata_inputs()
    
    # Failsafe in case preparse somehow misses these
    if sample is None:
        phased = False
    
    if groups is None:
        groups = list(adata.var["group"].dropna().unique())

    # Process initial minimums for whole data dispersion
    # region_cutoff = min_count + (2*pseudocount)
    snp_cutoff = (2 * pseudocount)
    
    ref_counts = adata.layers["ref"].sum(axis=1, dtype=np.uint16).T.A1 + pseudocount
    alt_counts = adata.layers["alt"].sum(axis=1, dtype=np.uint16).T.A1 + pseudocount
    n_counts = ref_counts + alt_counts

    # Calculate dispersion across dataset
    opt_disp = lambda rho, ref_data, n_data: -np.sum(
        betabinom.logpmf(ref_data, n_data, (0.5 * (1 - rho) / rho), (0.5 * (1 - rho) / rho))
    )
    
    disp = minimize_scalar(opt_disp, args=(ref_counts, n_counts), method="bounded", bounds=(0, 1))["x"]
    
    print(disp)  # DEEBUG BY SHOWING DISP
    
    df_dict: Dict[str, pd.DataFrame] = {}
    
    # Loop through groups
    for group_name in groups:
        
        # Subset by group
        adata_sub = adata[:, adata.var["group"] == group_name]
        
        # Create count data per group
        ref_counts_group = adata_sub.layers["ref"].sum(axis=1, dtype=np.uint16).T.A1 + pseudocount
        alt_counts_group = adata_sub.layers["alt"].sum(axis=1, dtype=np.uint16).T.A1 + pseudocount
        n_counts_group = ref_counts_group + alt_counts_group
        
        nonzero_idx = np.where(n_counts_group > snp_cutoff)  # Get indices where counts were found
        
        if nonzero_idx[0].size == 0:
            print(f"Skipping {group_name}: No SNP counts found")
            continue
        
        # Remove snps with 0 counts from regions
        idx_df = pd.DataFrame({"index": nonzero_idx[0]}, dtype=np.uint32).reset_index(names="filt_index")
        region_idx_df = adata.uns["feature"].merge(idx_df, on="index")
        
        # Check total allele counts/N per region
        region_n_df = region_idx_df.merge(
            pd.DataFrame(n_counts_group, columns=["N"]).reset_index(names="index"),
            on="index"
        )
        
        # region_n_df = adata.uns["feature"].merge(
        #     pd.DataFrame(n_counts_group, columns=["N"]).reset_index(names="index"),
        #     on="index")
        
        
        # Take into account pseudocounts added to total N
        region_agg_df = region_n_df.groupby("region", sort=False).agg(
            snp_idx=("index", tuple), num_snps=("index", "size"), N=("N", np.sum)
        )
        
        # Take into account pseudocounts added to total N
        region_agg_df["region_cutoff"] = (region_agg_df["num_snps"] * snp_cutoff) + min_count
        
        # Per group snp_dict
        region_snp_dict: Dict[str, tuple] = region_agg_df.loc[
            region_agg_df["N"] >= region_agg_df["region_cutoff"], "snp_idx"
        ].to_dict()
        # region_snp_dict = region_agg_df.loc[region_agg_df["N"] >= region_cutoff, "snp_idx"].to_dict()
        
        if not region_snp_dict:
            print(f"Skipping {group_name}: No regions with total allele counts >= {min_count}")
            continue

        if phased:
            gt_array = adata.obs[sample].str.split("|", n=1).str[0].to_numpy(dtype=np.uint8)
        else:
            gt_array = None

        # CREATE sub function that processes subgroup
        df = get_imbalance_per_group(
            ref_counts_group,
            n_counts_group,
            region_snp_dict,
            disp,
            gt_array=gt_array
        )
        
        df_dict[group_name] = df
        
    # Should I return something?
    # Maybe compile all of the dataframes?
    
    return df_dict


def get_imbalance_per_group(ref_counts: np.ndarray,
                            n_counts: np.ndarray,
                            region_snp_dict: Dict[str, tuple],
                            disp: float,
                            gt_array: Optional[np.ndarray] = None) -> pd.DataFrame:
    """
    Compute allelic imbalance for a single group across multiple regions.
    
    This function iterates over regions defined in `region_snp_dict` and computes the likelihood 
    of allelic imbalance using a beta-binomial model. Different optimization routines are applied 
    depending on whether genotype phasing data is available.
    
    Parameters
    ----------
    ref_counts : numpy.ndarray
        Array of reference allele counts for SNPs.
    n_counts : numpy.ndarray
        Array of total allele counts (reference + alternative) for SNPs.
    region_snp_dict : dict
        Dictionary mapping region identifiers to tuples of SNP indices.
    disp : float
        Dispersion parameter estimated from the dataset using a beta-binomial model.
    gt_array : numpy.ndarray, optional
        Array of genotype phasing information. If provided, indicates that data is phased.
    
    Returns
    -------
    pandas.DataFrame
        DataFrame containing imbalance statistics for each region with columns:
        - 'region': Region identifier.
        - 'num_snps': Number of SNPs in the region.
        - 'mu': Estimated imbalance parameter.
        - 'null_ll': Log-likelihood under the null hypothesis.
        - 'alt_ll': Log-likelihood under the alternative hypothesis.
        - 'pval': p-value from the likelihood ratio test.
        - 'fdr_pval': False discovery rate (FDR) corrected p-value.
    
    Notes
    -----
    - For regions with more than one SNP:
      - If phased data is available, the optimization is performed with respect to phased counts.
      - Otherwise, an unphased optimization routine is used.
    - For regions with a single SNP, the imbalance is computed directly.
    - The likelihood ratio test (LRT) is used to compare the null and alternative models.
    
    Examples
    --------
    >>> df = get_imbalance_per_group(ref_counts, n_counts, region_snp_dict, disp)
    """
    
    # Check if genotype phasing info available
    if gt_array is None:
        phased = False
    else:
        phased = True
    
    group_results: List[tuple] = []  # Store imbalance results
    
    # Would the old method of grouped dataframe work better?
    for region, snp_list in region_snp_dict.items():

        region_ref = ref_counts[list(snp_list)]
        region_n = n_counts[list(snp_list)]

        # Null test
        null_ll = np.sum(betabinom.logpmf(
            region_ref, region_n, (0.5 * (1 - disp) / disp), (0.5 * (1 - disp) / disp)
        ))

        # Handle phasing stuff
        snp_count = region_ref.shape[0]

        if snp_count > 1:

            if phased:
                
                region_gt = gt_array[list(snp_list)]
                
                # Make sure phase with respect to first snp ref
                if region_gt[0] > 0:
                    region_gt = 1 - region_gt

                res = minimize_scalar(
                    opt_phased_new,
                    args=(disp, region_ref, region_n, region_gt),
                    method="bounded",
                    bounds=(0, 1)
                )
                mu = res["x"]
                opt_ll = res["fun"]

            else:
                first_ref = region_ref[:1]
                first_n = region_n[:1]

                phase_ref = region_ref[1:]
                phase_n = region_n[1:]

                # Using some minimize scalar
                res = minimize_scalar(
                    opt_unphased_dp,
                    args=(disp, first_ref, first_n, phase_ref, phase_n),
                    method="bounded",
                    bounds=(0, 1)
                )

                mu = res["x"]
                opt_ll = res["fun"]

        else:
            # If only one snp
            if 0 < region_ref[0] < region_n[0]:
                mu = region_ref[0] / region_n[0]
                opt_ll = opt_prob(mu, disp, int(region_ref[0]), int(region_n[0]))
            else:
                res = minimize_scalar(
                    opt_prob,
                    args=(disp, int(region_ref[0]), int(region_n[0])),
                    method="bounded",
                    bounds=(0, 1)
                )
                # Get res data
                mu = res["x"]
                opt_ll = res["fun"]

        # Process LRT
        alt_ll = -1 * opt_ll

        # OUTSIDE OF FUNCTION
        lrt = -2 * (null_ll - alt_ll)
        pval = chi2.sf(lrt, 1)

        # Add data to output list
        group_results.append(
            (region, snp_count, mu, null_ll, alt_ll, pval)
        )
    
    # Create allelic imbalance df
    df = pd.DataFrame(group_results,
                      columns=["region", "num_snps", "mu",
                               "null_ll", "alt_ll", "pval"]
                     )

    # fdr correction
    df["fdr_pval"] = false_discovery_control(df["pval"], method="bh")
    
    return df
