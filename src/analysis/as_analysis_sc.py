"""Single-cell allelic imbalance analysis functions.

Provides functions for analyzing allelic imbalance in single-cell data
stored in AnnData format with SNP counts in layers.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from anndata import AnnData
from numpy.typing import NDArray
from scipy.optimize import OptimizeResult, minimize_scalar
from scipy.stats import betabinom, chi2, false_discovery_control, zscore

# Local imports
from .as_analysis import opt_phased_new, opt_prob, opt_unphased_dp


def adata_count_qc(
    adata: AnnData, z_cutoff: float | None = None, gt_error: Any | None = None
) -> AnnData:
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
        adata.uns["feature"] = (
            adata.uns["feature"]
            .merge(adata.obs[["index"]].reset_index(names="filt_index"), on="index")[
                ["region", "filt_index"]
            ]
            .rename(columns={"filt_index": "index"})
        )

        adata.obs["index"] = adata.obs.index  # Replace index column

    # TODO add options to identify and filter GT errors
    if gt_error is not None:
        pass

    return adata


def get_imbalance_sc(
    adata: AnnData,
    min_count: int = 10,
    pseudocount: int = 1,
    phased: bool = False,
    sample: str | None = None,
    groups: list[str] | None = None,
) -> dict[str, pd.DataFrame]:
    # Need to preparse input using process_adata_inputs()

    # Failsafe in case preparse somehow misses these
    if sample is None:
        phased = False

    if groups is None:
        groups = list(adata.var["group"].dropna().unique())

    # Process initial minimums for whole data dispersion
    # region_cutoff = min_count + (2*pseudocount)
    snp_cutoff = 2 * pseudocount

    ref_counts = adata.layers["ref"].sum(axis=1, dtype=np.uint16).T.A1 + pseudocount
    alt_counts = adata.layers["alt"].sum(axis=1, dtype=np.uint16).T.A1 + pseudocount
    n_counts = ref_counts + alt_counts

    # Calculate dispersion across dataset
    def opt_disp(rho: float, ref_data: NDArray[np.uint16], n_data: NDArray[np.uint16]) -> float:
        return float(
            -np.sum(
                betabinom.logpmf(ref_data, n_data, (0.5 * (1 - rho) / rho), (0.5 * (1 - rho) / rho))
            )
        )

    disp_result: OptimizeResult = minimize_scalar(
        opt_disp, args=(ref_counts, n_counts), method="bounded", bounds=(0, 1)
    )
    disp: float = float(disp_result["x"])

    df_dict: dict[str, pd.DataFrame] = {}

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
        idx_df = pd.DataFrame({"index": nonzero_idx[0]}, dtype=np.uint32).reset_index(
            names="filt_index"
        )
        region_idx_df = adata.uns["feature"].merge(idx_df, on="index")

        # Check total allele counts/N per region
        region_n_df = region_idx_df.merge(
            pd.DataFrame(n_counts_group, columns=["N"]).reset_index(names="index"), on="index"
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
        region_snp_dict = region_agg_df.loc[
            region_agg_df["N"] >= region_agg_df["region_cutoff"], "snp_idx"
        ].to_dict()
        # region_snp_dict = region_agg_df.loc[region_agg_df["N"] >= region_cutoff, "snp_idx"].to_dict()

        if not region_snp_dict:
            print(f"Skipping {group_name}: No regions with total allele counts >= {min_count}")
            continue

        gt_array_typed: NDArray[np.uint8] | None
        if phased:
            gt_array_typed = adata.obs[sample].str.split("|", n=1).str[0].to_numpy(dtype=np.uint8)
        else:
            gt_array_typed = None

        # CREATE sub function that processes subgroup
        df: pd.DataFrame = get_imbalance_per_group(
            ref_counts_group, n_counts_group, region_snp_dict, disp, gt_array=gt_array_typed
        )

        df_dict[group_name] = df

    # Should I return something?
    # Maybe compile all of the dataframes?

    return df_dict


def get_imbalance_per_group(
    ref_counts: NDArray[np.integer[Any]],
    n_counts: NDArray[np.integer[Any]],
    region_snp_dict: dict[int, tuple[int, ...]],
    disp: float,
    gt_array: NDArray[np.uint8] | None = None,
) -> pd.DataFrame:
    # Check if genotype phasing info available
    phased: bool
    if gt_array is None:
        phased = False
    else:
        phased = True

    group_results: list[tuple[int, int, float, float, float, float]] = []  # Store imbalance results

    # Would the old method of grouped dataframe work better?
    for region, snp_list in region_snp_dict.items():
        region_ref: NDArray[np.integer[Any]] = ref_counts[snp_list,]
        region_n: NDArray[np.integer[Any]] = n_counts[snp_list,]

        # Null test
        null_ll: float = float(
            np.sum(
                betabinom.logpmf(
                    region_ref, region_n, (0.5 * (1 - disp) / disp), (0.5 * (1 - disp) / disp)
                )
            )
        )

        # Handle phasing stuff
        snp_count: int = region_ref.shape[0]

        if snp_count > 1:
            if phased:
                assert gt_array is not None  # Type guard for mypy
                region_gt: NDArray[np.uint8] = gt_array[snp_list,]

                # Make sure phase with respect to first snp ref
                if region_gt[0] > 0:
                    region_gt = 1 - region_gt

                res: OptimizeResult = minimize_scalar(
                    opt_phased_new,
                    args=(disp, region_ref, region_n, region_gt),
                    method="bounded",
                    bounds=(0, 1),
                )
                mu: float = float(res["x"])
                opt_ll: float = float(res["fun"])

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
                    bounds=(0, 1),
                )

                mu = float(res["x"])
                opt_ll = float(res["fun"])

        else:
            # If only one snp
            if 0 < region_ref[0] < region_n[0]:
                mu = float(region_ref[0]) / float(region_n[0])
                opt_ll_result = opt_prob(mu, disp, region_ref[0], region_n[0])
                opt_ll = float(opt_ll_result)
            else:
                res = minimize_scalar(
                    opt_prob,
                    args=(disp, region_ref[0], region_n[0]),
                    method="bounded",
                    bounds=(0, 1),
                )
                # Get res data
                mu = float(res["x"])
                opt_ll = float(res["fun"])

        # Process LRT
        alt_ll: float = -1 * opt_ll

        # OUTSIDE OF FUNCTION
        lrt: float = -2 * (null_ll - alt_ll)
        pval: float = float(chi2.sf(lrt, 1))

        # Add data to output list
        group_results.append((region, snp_count, mu, null_ll, alt_ll, pval))

    # Create allelic imbalance df
    # Polars vs pandas??
    df: pd.DataFrame = pd.DataFrame(
        group_results, columns=["region", "num_snps", "mu", "null_ll", "alt_ll", "pval"]
    )

    # fdr correction
    df["fdr_pval"] = false_discovery_control(df["pval"], method="bh")

    return df
