import sys
import warnings
from pathlib import Path


import numpy as np
import pandas as pd

import anndata as ad

from scipy.stats import betabinom, chi2
from scipy.optimize import minimize_scalar

# Local imports
from as_analysis import opt_prob, opt_phased_new, opt_unphased_dp, bh_correction


def get_imbalance_sc(adata,
                     min_count=10,
                     pseudocount=1,
                     phased=False,
                     sample=None,
                     groups=None):
    
    # Need to preparse input using process_adata_inputs()
    
    # Failsafe in case preparse somehow misses these
    if sample is None:
        phased = False

    if groups is None:
        groups = list(adata.var["group"].dropna().unique())

    
    # Process initial minimums for whole data dispersion
    cutoff = min_count + (2*pseudocount)
    
    ref_counts = adata.layers["ref"].sum(axis=1, dtype=np.uint16).T.A1 + pseudocount
    alt_counts = adata.layers["alt"].sum(axis=1, dtype=np.uint16).T.A1 + pseudocount
    
    n_counts = ref_counts + alt_counts
    min_idx = np.where(n_counts >= cutoff) # Get indices for min_count

    ref_counts_filt, n_counts_filt = ref_counts[min_idx], n_counts[min_idx]
    
    # Calculate dispersion across dataset
    opt_disp = lambda rho, ref_data, n_data: -np.sum(
        betabinom.logpmf(ref_data, n_data, (0.5 * (1 - rho) / rho), (0.5 * (1 - rho) / rho))
    )
    
    disp = minimize_scalar(opt_disp, args=(ref_counts_filt, n_counts_filt), method="bounded", bounds=(0,1))["x"]

    df_dict = {}
    
    # Loop through groups
    for group_name in groups:
        
        # Subset by group
        adata_sub = adata[:, adata.var["group"] == group_name]
        
        # Create count data per group
        ref_counts_group = adata_sub.layers["ref"].sum(axis=1, dtype=np.uint16).T.A1 + pseudocount
        alt_counts_group = adata_sub.layers["alt"].sum(axis=1, dtype=np.uint16).T.A1 + pseudocount
        n_counts_group = ref_counts_group + alt_counts_group

        min_idx_group = np.where(n_counts_group >= cutoff) # Get indices for min_count

        if min_idx_group[0].size == 0:
            print(f"Skipping {group_name}: No SNP's with allele counts >= {min_count}")
            continue

        ref_counts_group_filt, n_counts_group_filt = ref_counts_group[min_idx_group], n_counts_group[min_idx_group]
        
        # Per group snp_dict
        idx_df = pd.DataFrame({"index": min_idx_group[0]}, dtype=np.uint32).reset_index(names="filt_index")
        
        region_snp_dict = adata.uns["feature"].merge(
            idx_df, on="index")[["region", "filt_index"]].groupby(
            "region", sort=False).agg(tuple)["filt_index"].to_dict()


        if phased:
            gt_array = adata.obs.iloc[min_idx_group][sample].str.split("|", n=1).str[0].to_numpy(dtype=np.uint8)
        else:
            gt_array = None

        
        # CREATE sub function that processes subgroup
        df = get_imbalance_per_group(ref_counts_group_filt, n_counts_group_filt,
                                     region_snp_dict, disp,
                                     gt_array=gt_array
                                    )
        
        df_dict[group_name] = df
        
    # Should I return something?
    # Maybe compile all of the dataframes?
    
    return df_dict


def get_imbalance_per_group(ref_counts,
                            n_counts,
                            region_snp_dict,
                            disp,
                            gt_array=None
                           ):
    
    # Check if genotype phasing info available
    if gt_array is None:
        phased = False
    else:
        phased = True
    
    group_results = [] # Store imbalance results
    
    # Would the old method of grouped dataframe work better?
    for region, snp_list in region_snp_dict.items():

        region_ref = ref_counts[snp_list,]
        region_n = n_counts[snp_list,]

        # Null test
        null_ll = np.sum(betabinom.logpmf(
            region_ref, region_n, (0.5 * (1 - disp) / disp), (0.5 * (1 - disp) / disp)))


        # Handle phasing stuff
        snp_count = region_ref.shape[0]

        if snp_count > 1:

            if phased:

                res = minimize_scalar(opt_phased_new,
                                      args=(disp, region_ref, region_n, gt_array[snp_list,]),
                                      method="bounded", bounds=(0, 1))
                mu = res["x"]
                opt_ll = res["fun"]

            else:
                first_ref = region_ref[:1]
                first_n = region_n[:1]

                phase_ref = region_ref[1:]
                phase_n = region_n[1:]


                # Using some minimize scalar
                res = minimize_scalar(opt_unphased_dp,
                                      args=(disp, first_ref, first_n, phase_ref, phase_n),
                                      method="bounded", bounds=(0, 1))

                mu = res["x"]
                opt_ll = res["fun"]

        else:

            # If only one snp
            if 0 < region_ref[0] < region_n[0]:
                mu = region_ref[0]/region_n[0]
                opt_ll = opt_prob(mu, disp, region_ref[0], region_n[0])
            else:
                res = minimize_scalar(opt_prob, args=(disp, region_ref[0], region_n[0]), 
                                      method="bounded", bounds=(0, 1))
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
    # Polars vs pandas??
    df = pd.DataFrame(group_results,
                      columns=["region", "num_snps", "mu",
                               "null_ll", "alt_ll", "pval"]
                     )

    # fdr correction
    df = bh_correction(df)
    
    return df
