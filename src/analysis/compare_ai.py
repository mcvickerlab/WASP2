import sys
import warnings
from pathlib import Path

from collections import namedtuple
from itertools import combinations

import numpy as np
import pandas as pd

from scipy.stats import betabinom, chi2, false_discovery_control
from scipy.optimize import minimize_scalar


# Local imports
from as_analysis import opt_prob, opt_unphased_dp, opt_phased_new, bh_correction


# Use these functions to figure out how to optimize per group
def get_imbalance_func(ref_count, n_count, phase_array=None):
    
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


def opt_combined_imbalance(prob, disp,
                           like_func1, like_func1_args,
                           like_func2, like_func2_args):
    
    return (like_func1(prob, disp, *like_func1_args) +
            like_func2(prob, disp, *like_func2_args))


# Current version that uses shared snps
def get_compared_imbalance(adata,
                           min_count=10,
                           pseudocount=1,
                           phased=False,
                           sample=None,
                           groups=None):
    
    # Failsafe in case preparse somehow misses these
    if sample is None:
        phased = False

    # Should I be comparing all combos by default??? Seems like a lot
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
    
    disp = minimize_scalar(opt_disp, args=(ref_counts, n_counts), method="bounded", bounds=(0,1))["x"]
    
    if phased:
        gt_array = adata.obs[sample].str.split("|", n=1).str[0].to_numpy(dtype=np.uint8)
    else:
        gt_array = None

    
    # process counts on a per group basis to avoid recalculating
    group_dict = {}
    # group_data = namedtuple("group_data", ["ref_counts", "n_counts", "phase_data", "region_snp_dict"]) # Maybe include the gt_array instead of min_idx
    group_data = namedtuple("group_data", ["ref_counts", "n_counts", "region_snp_df"])
    
    for group_name in groups:

        # Subset by group
        adata_sub = adata[:, adata.var["group"] == group_name]

        # Create count data per group, should i do pseudocount now or later?
        ref_counts_group = adata_sub.layers["ref"].sum(axis=1, dtype=np.uint16).T.A1 + pseudocount
        alt_counts_group = adata_sub.layers["alt"].sum(axis=1, dtype=np.uint16).T.A1 + pseudocount
        n_counts_group = ref_counts_group + alt_counts_group
        
        nonzero_idx = np.where(n_counts_group > snp_cutoff) # Get indices where no counts were found
    
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

        group_dict[group_name] = group_data(ref_counts_group, n_counts_group, region_n_df)
    
    
    # Create group combinations and process shared snps
    group_combos = list(combinations(group_dict.keys(), r=2))
    
    df_dict = {}
    for group1, group2 in group_combos:
        
        # Get relevant counts and nonzero snps
        ref_counts1, n_counts1, region_snp_df1 = group_dict[group1]
        ref_counts2, n_counts2, region_snp_df2 = group_dict[group2]
        
        
        # Get shared snps -> get regions that meet cutoff
        shared_df = region_snp_df1[["region", "index", "N"]].merge(
            region_snp_df2[["index", "N"]], on="index", suffixes=("1", "2"))
        
        
        # Take into account pseudocounts added to total N
        region_agg_df = shared_df.groupby("region", sort=False).agg(
            snp_idx=("index", tuple), num_snps=("index", "size"),
            N1=("N1", np.sum), N2=("N2", np.sum)
        )
        
        region_agg_df["region_cutoff"] = (region_agg_df["num_snps"] * snp_cutoff) + min_count


        # Find regions where N is satisfied for both 
        # region_agg_df = shared_df.groupby("region", sort=False).agg(
        #     snp_idx=("index", tuple), N1=("N1", np.sum), N2=("N2", np.sum)
        # )

        # Per group snp_dict
        region_snp_dict = region_agg_df.loc[
            (
                (region_agg_df["N1"] >= region_agg_df["region_cutoff"]) & 
                (region_agg_df["N2"] >= region_agg_df["region_cutoff"])
                ),
            "snp_idx"].to_dict()
        
        # region_snp_dict = region_agg_df.loc[
        #     (region_agg_df["N1"] >= region_cutoff) & (region_agg_df["N2"] >= region_cutoff),
        #     "snp_idx"].to_dict()

        if not region_snp_dict:
            print(
                (f"Skipping {group1}-{group2} Comparison: "
                 f"No shared regions with allele counts >= {min_count}"
                )
            )

            continue


        # This sub function name kinda long...find better name maybe?
        df = compare_imbalance_between_groups(disp,
                                              ref_counts1,
                                              n_counts1,
                                              ref_counts2,
                                              n_counts2,
                                              region_snp_dict,
                                              gt_array
                                              )
        
        # Using a tuple as key
        df_dict[(group1, group2)] = df

    return df_dict


def compare_imbalance_between_groups(disp,
                                     ref_counts1,
                                     n_counts1,
                                     ref_counts2,
                                     n_counts2,
                                     region_snp_dict,
                                     gt_array=None
                                     ):
    
    # Helper func called by get_compared_imbalance()
    
    group_results = [] # Store imbalance results
    
    # Compare allelic imbalance difference per region
    for region, snp_list in region_snp_dict.items():
        
        # Get per region snps and counts
        region_ref1 = ref_counts1[snp_list,]
        region_n1 = n_counts1[snp_list,]

        region_ref2 = ref_counts2[snp_list,]
        region_n2 = n_counts2[snp_list,]

        
        # Process which model we'll use to process likelihood per group
        if len(snp_list) == 1:
            # Parse single opt
            like_func = opt_prob

            # This excludes disp since we always use disp
            like_func_args1 = (region_ref1[0], region_n1[0])
            like_func_args2 = (region_ref2[0], region_n2[0])

        elif gt_array is None:
            # Do unphased
            like_func = opt_unphased_dp

            like_func_args1 = (region_ref1[:1], region_n1[:1],
                               region_ref1[1:], region_n1[1:])

            like_func_args2 = (region_ref2[:1], region_n2[:1],
                               region_ref2[1:], region_n2[1:])

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
        null_res = minimize_scalar(opt_combined_imbalance,
                                   args=(disp,
                                         like_func, like_func_args1,
                                         like_func, like_func_args2), 
                                   method="bounded", bounds=(0, 1))

        combined_mu = null_res["x"]
        null_ll = -1 * null_res["fun"]


        # Alt Hypothesis: Imbalance is different between groups
        alt_res1 = minimize_scalar(like_func,
                                   args=(disp, *like_func_args1),
                                   method="bounded", bounds=(0, 1))

        alt_res2 = minimize_scalar(like_func,
                                   args=(disp, *like_func_args2),
                                   method="bounded", bounds=(0, 1))

        # Get separate mu
        alt_mu1 = alt_res1["x"]
        alt_mu2 = alt_res2["x"]

        # get Alternative likelihood
        alt_ll1 = alt_res1["fun"]
        alt_ll2 = alt_res2["fun"]

        alt_ll = -1 * (alt_ll1 + alt_ll2)

        # Log ratio ttest
        lrt = -2 * (null_ll - alt_ll)
        pval = chi2.sf(lrt, 1)

        # Add data to output list
        
        # How should i format this, lots of possible outputs
        group_results.append(
            (region, len(snp_list), combined_mu, alt_mu1, alt_mu2, null_ll, alt_ll, pval)
        )
        
    # Create allelic imbalance df
    
    # Polars implementation might be more performant
    df = pd.DataFrame(group_results,
                      columns=["region",
                               "num_snps", 
                               "combined_mu",
                               "mu1", "mu2",
                               "null_ll",
                               "alt_ll",
                               "pval"]
                     )
    
    # fdr correction
    df["fdr_pval"] = false_discovery_control(df["pval"], method="bh")
    
    return df


# THIS IS A V0 VERSION THAT DIDN'T USE SHARED SNPS BETWEEN REGIONS
# COULD BE USEFUL AS AN OPTION POSSIBLY
def get_compared_imbalance_diff_snps(adata,
                           min_count=10,
                           pseudocount=1,
                           phased=False,
                           sample=None,
                           groups=None):
    
    # Failsafe in case preparse somehow misses these
    if sample is None:
        phased = False

    # Should I be comparing all combos by default??? Seems like a lot
    if groups is None:
        groups = list(adata.var["group"].dropna().unique())
        print("Comparing all combinations of available groups")
    elif len(groups) == 1:
        raise ValueError("Please provide 2 or more groups to compare.")
    

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

    # process counts on a per group basis to avoid recalculating
    group_dict = {}
    group_data = namedtuple("group_data", ["ref_counts", "n_counts", "phase_data", "region_snp_dict"]) # Maybe include the gt_array instead of min_idx

    for group_name in groups:

        # Subset by group
        adata_sub = adata[:, adata.var["group"] == group_name]

        # Create count data per group, should i do pseudocount now or later?
        ref_counts_group = adata_sub.layers["ref"].sum(axis=1, dtype=np.uint16).T.A1 + pseudocount
        alt_counts_group = adata_sub.layers["alt"].sum(axis=1, dtype=np.uint16).T.A1 + pseudocount
        n_counts_group = ref_counts_group + alt_counts_group

        min_idx_group = np.where(n_counts_group >= cutoff) # Get indices for min_count

        ref_counts_group_filt, n_counts_group_filt = ref_counts_group[min_idx_group], n_counts_group[min_idx_group]

        if phased:
            phase_array = adata.obs.iloc[min_idx_group][sample].str.split("|", n=1).str[0].to_numpy(dtype=np.uint8)
        else:
            phase_array = None

        # Create region_snp_dict but for each group
        idx_df = pd.DataFrame({"index": min_idx_group[0]}, dtype=np.uint32).reset_index(names="filt_index")

        region_snp_dict = adata.uns["feature"].merge(
            idx_df, on="index")[["region", "filt_index"]].groupby(
            "region", sort=False).agg(tuple)["filt_index"].to_dict()

        group_dict[group_name] = group_data(ref_counts_group_filt, n_counts_group_filt,
                                            phase_array, region_snp_dict)
    
    # Create group combinations and process shared snps
    group_combos = list(combinations(group_dict.keys(), r=2))
    
    df_dict = {}
    for group1, group2 in group_combos:

        # Might be smart to create a cache to prevent repeating calculations
        # This sub function name kinda long...find better name maybe?
        df = compare_imbalance_between_groups_diff_snps(disp,
                                              *group_dict[group1],
                                              *group_dict[group2]
                                             )
        
        if df.empty:
            print(f"Skipping {group1} - {group2} comparison. No shared regions.")
        else:
            # Using a tuple as key
            df_dict[(group1, group2)] = df

        
    return df_dict


def compare_imbalance_between_groups_diff_snps(disp,
                                     ref_counts1,
                                     n_counts1,
                                     phase_array1,
                                     region_snp_dict1,
                                     ref_counts2,
                                     n_counts2,
                                     phase_array2,
                                     region_snp_dict2):
    
    # These values are unpacked versions of named tuple
    # Helper func called by get_compared_imbalance()
    
    # Check if phasing info available
    phased = ((phase_array1 is not None) and
              (phase_array2 is not None))
    
    # Get shared regions
    shared_regions = [i for i in region_snp_dict1.keys()
                      if i in region_snp_dict2]
    
    
    group_results = [] # Store imbalance results
    
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
            region_phasing1 = phase_array1[snp_list1,]
            region_phasing2 = phase_array2[snp_list2,]
        else:
            region_phasing1, region_phasing2 = None, None
        
        # Process which model we'll use to process likelihood per group
        like_func1, like_func_inputs1 = get_imbalance_func(
            region_ref1, region_n1, phase_array=region_phasing1)
        
        like_func2, like_func_inputs2 = get_imbalance_func(
            region_ref2, region_n2, phase_array=region_phasing2)


        # Null Hypothesis: Imbalance is the same
        null_res = minimize_scalar(opt_combined_imbalance,
                                   args=(disp,
                                         like_func1, like_func_inputs1,
                                         like_func2, like_func_inputs2), 
                                   method="bounded", bounds=(0, 1))

        combined_mu = null_res["x"]
        null_ll = -1 * null_res["fun"]


        # Alt Hypothesis: Imbalance is different between groups
        alt_res1 = minimize_scalar(like_func1,
                                   args=(disp, *like_func_inputs1),
                                   method="bounded", bounds=(0, 1))


        alt_res2 = minimize_scalar(like_func2,
                                   args=(disp, *like_func_inputs2),
                                   method="bounded", bounds=(0, 1))


        # Get separate mu
        alt_mu1 = alt_res1["x"]
        alt_mu2 = alt_res2["x"]

        # get Alternative likelihood
        alt_ll = -1 * (alt_res1["fun"] + alt_res2["fun"])


        # Log ratio ttest
        lrt = -2 * (null_ll - alt_ll)
        pval = chi2.sf(lrt, 1)

        # Add data to output list
        
        # How should i format this, lots of possible outputs
        group_results.append(
            (region, len(snp_list1), len(snp_list2), combined_mu, alt_mu1, alt_mu2, null_ll, alt_ll, pval)
        )
        
    # Create allelic imbalance df
    
    # Polars implementation might be more performant
    df = pd.DataFrame(group_results,
                      columns=["region",
                               "num_snps_group1", "num_snps_group2",
                               "combined_mu", "mu1", "mu2",
                               "null_ll", "alt_ll", "pval"]
                     )
    
    # fdr correction
    df = bh_correction(df)
    
    return df

