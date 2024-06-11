import sys
import warnings
from pathlib import Path


import numpy as np
import pandas as pd

import anndata as ad

from scipy.stats import betabinom, chi2, zscore, false_discovery_control
from scipy.optimize import minimize_scalar

# Local imports
from as_analysis import opt_prob, opt_phased_new, opt_unphased_dp, bh_correction


# Performs qc and prefilters anndata count data
# Should this be a decorator instead?
def adata_count_qc(adata, z_cutoff=None, gt_error=None):
    
    # No need to prefilt
    if z_cutoff is None and gt_error is None:
        return adata
    
    # Filt outliers
    if z_cutoff is not None:
        snp_outliers = adata.obs[["index", "ref_count", "alt_count"]].copy()
        snp_outliers["N"] = snp_outliers["ref_count"] + snp_outliers["alt_count"]
        snp_outliers = snp_outliers[np.abs(zscore(snp_outliers["N"])) > z_cutoff] # At least 3
        
        # Todo: add option if there aren't any features
        # Get regions containing 1 or more outlier snps
        snp_outliers = snp_outliers.merge(adata.uns["feature"], on="index", how="left")
        
        outlier_regions = adata.uns["feature"].loc[adata.uns["feature"]["region"].isin(
            snp_outliers["region"].unique()), :]
        
        # Remove outlier regions from adata
        adata = adata[~adata.obs["index"].isin(outlier_regions["index"]), :].copy()
        adata.obs = adata.obs.reset_index(drop=True) # update index
        
        # Update valid regions and snps
        adata.uns["feature"] = adata.uns["feature"].merge(
            adata.obs[["index"]].reset_index(names="filt_index"),
            on="index")[["region", "filt_index"]].rename(
            columns={"filt_index": "index"})
        
        adata.obs["index"] = adata.obs.index # Replace index column
    
    # TODO add options to identify and filter GT errors
    if gt_error is not None:
        pass
    
    return adata


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
    # region_cutoff = min_count + (2*pseudocount)
    snp_cutoff = (2*pseudocount)
    
    ref_counts = adata.layers["ref"].sum(axis=1, dtype=np.uint16).T.A1 + pseudocount
    alt_counts = adata.layers["alt"].sum(axis=1, dtype=np.uint16).T.A1 + pseudocount
    n_counts = ref_counts + alt_counts

    # Calculate dispersion across dataset
    opt_disp = lambda rho, ref_data, n_data: -np.sum(
        betabinom.logpmf(ref_data, n_data, (0.5 * (1 - rho) / rho), (0.5 * (1 - rho) / rho))
    )
    
    disp = minimize_scalar(opt_disp, args=(ref_counts, n_counts), method="bounded", bounds=(0,1))["x"]
    
    print(disp) # DEEBUG BY SHOWING DISP
    
    df_dict = {}
    
    # Loop through groups
    for group_name in groups:
        
        # Subset by group
        adata_sub = adata[:, adata.var["group"] == group_name]
        
        # Create count data per group
        ref_counts_group = adata_sub.layers["ref"].sum(axis=1, dtype=np.uint16).T.A1 + pseudocount
        alt_counts_group = adata_sub.layers["alt"].sum(axis=1, dtype=np.uint16).T.A1 + pseudocount
        n_counts_group = ref_counts_group + alt_counts_group
        
        nonzero_idx = np.where(n_counts_group > snp_cutoff) # Get indices where counts were found
        
        if nonzero_idx[0].size == 0:
            print(f"Skipping {group_name}: No SNP counts found")
            continue
        
        # Remove snps with 0 counts from regions
        idx_df = pd.DataFrame({"index": nonzero_idx[0]}, dtype=np.uint32).reset_index(names="filt_index")
        region_idx_df = adata.uns["feature"].merge(idx_df, on="index")

        # Check total allele counts/N per region
        region_n_df = region_idx_df.merge(
            pd.DataFrame(n_counts_group, columns=["N"]).reset_index(names="index"),
            on="index")
        
        # region_n_df = adata.uns["feature"].merge(
        #     pd.DataFrame(n_counts_group, columns=["N"]).reset_index(names="index"),
        #     on="index")
        
        
        # Take into account pseudocounts added to total N
        region_agg_df = region_n_df.groupby("region", sort=False).agg(
            snp_idx=("index", tuple), num_snps=("index", "size"), N=("N", np.sum))
        
        # Take into account pseudocounts added to total N
        region_agg_df["region_cutoff"] = (region_agg_df["num_snps"] * snp_cutoff) + min_count
        
        # Per group snp_dict
        region_snp_dict = region_agg_df.loc[region_agg_df["N"] >= region_agg_df["region_cutoff"], "snp_idx"].to_dict()
        # region_snp_dict = region_agg_df.loc[region_agg_df["N"] >= region_cutoff, "snp_idx"].to_dict()
        
        if not region_snp_dict:
            print(f"Skipping {group_name}: No regions with total allele counts >= {min_count}")
            continue

        if phased:
            gt_array = adata.obs[sample].str.split("|", n=1).str[0].to_numpy(dtype=np.uint8)
        else:
            gt_array = None

        # CREATE sub function that processes subgroup
        df = get_imbalance_per_group(ref_counts_group,
                                     n_counts_group,
                                     region_snp_dict,
                                     disp,
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
                
                region_gt = gt_array[snp_list,]
                
                # Make sure phase with respect to first snp ref
                if region_gt[0] > 0:
                    region_gt = 1 - region_gt

                res = minimize_scalar(opt_phased_new,
                                      args=(disp, region_ref, region_n, region_gt),
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
    df["fdr_pval"] = false_discovery_control(df["pval"], method="bh")
    
    return df
