"""
as_analysis.py
Author: Aaron Ho
Python Version: 3.9
"""

# Default Python package Imports
from pathlib import Path
import time
import timeit
from typing import Tuple, List, Union, Optional

# External package imports
import pandas as pd
import numpy as np
from scipy.stats import betabinom, chi2, binom, rankdata, false_discovery_control
from scipy.optimize import minimize_scalar, minimize
from scipy.special import expit


def opt_linear(disp_params: Tuple[float, float],
               ref_counts: np.ndarray,
               n_array: np.ndarray) -> float:
    """
    Optimize dispersion parameter weighted by N.
    
    This function is called by an optimizer to compute the negative log-likelihood 
    of a beta-binomial model where the dispersion parameter is modeled as a linear 
    function of the total count (N).

    :param disp_params: Tuple of dispersion parameters (disp1, disp2).
    :type disp_params: tuple of floats
    :param ref_counts: Array of reference allele counts.
    :type ref_counts: numpy.ndarray
    :param n_array: Array of total counts (ref + alt) for each SNP.
    :type n_array: numpy.ndarray
    :return: Negative log-likelihood.
    :rtype: float
    """
    disp1, disp2 = disp_params

    exp_in = disp1 + (n_array * disp2)
    exp_in = np.select([exp_in > 10, exp_in < -10], [10, -10], default=exp_in)

    rho = expit(exp_in)

    ll = -np.sum(betabinom.logpmf(ref_counts, n_array, (0.5 * (1 - rho) / rho),
                                  (0.5 * (1 - rho) / rho)))  # If alpha is beta
    
    return ll


def opt_prob(in_prob: float,
             in_rho: float,
             k: int,
             n: int,
             log: bool = True) -> float:
    """
    Optimize the probability value to maximize imbalance likelihood.
    
    This function calculates the (negative) log-likelihood (or likelihood) of the 
    beta-binomial distribution given a probability parameter.

    :param in_prob: The input probability.
    :type in_prob: float
    :param in_rho: Dispersion or rate parameter.
    :type in_rho: float
    :param k: Observed reference allele count.
    :type k: int
    :param n: Total count (ref + alt).
    :type n: int
    :param log: If True, returns negative log-likelihood; otherwise returns likelihood.
    :type log: bool, optional
    :return: The (negative) likelihood value.
    :rtype: float
    """
    prob = in_prob

    alpha = prob * (1 - in_rho) / in_rho
    beta = (1 - prob) * (1 - in_rho) / in_rho
     
    if log:
        ll = -1 * betabinom.logpmf(k, n, alpha, beta)
    else:
        ll = betabinom.pmf(k, n, alpha, beta)

    return ll


def opt_phased(prob: float,
               first_data: np.ndarray,
               phase_data: np.ndarray) -> float:
    """
    Optimize likelihood while considering phasing information.

    This function computes the total (negative) log-likelihood by combining 
    the likelihood from the first SNP (used as reference) with the likelihoods 
    calculated for the phased data.

    :param prob: Probability parameter.
    :type prob: float
    :param first_data: Tuple/list containing (dispersion, ref_count, total count) for the first SNP.
    :type first_data: list or tuple of numbers (or numpy.ndarray of length 3)
    :param phase_data: Tuple/list containing (dispersion, ref_count, total count) for phased SNPs.
    :type phase_data: list or tuple of numbers (or numpy.ndarray with three rows)
    :return: Combined (negative) log-likelihood.
    :rtype: float
    """
    first_ll: float = opt_prob(prob, first_data[0], int(first_data[1]), int(first_data[2]))
    
    # Compute likelihoods for phased data
    phase1_lls: float = opt_prob(prob, phase_data[0], int(phase_data[1]), int(phase_data[2]), log=False)
    phase2_lls: float = opt_prob(1 - prob, phase_data[0], int(phase_data[1]), int(phase_data[2]), log=False)

    combined_lls = (0.5 * phase1_lls) + (0.5 * phase2_lls)
    return first_ll - np.sum(np.log(combined_lls))


def opt_phased_new(prob: float,
                   disp: float,
                   ref_data: np.ndarray,
                   n_data: np.ndarray,
                   gt_data: np.ndarray) -> float:
    """
    Updated phasing optimizer for single-cell analysis.

    This function adjusts the probability with respect to the first SNP 
    (used as reference) and then computes the negative log-likelihood using the 
    beta-binomial model.

    :param prob: Input probability.
    :type prob: float
    :param disp: Dispersion parameter.
    :type disp: float
    :param ref_data: Array of reference allele counts.
    :type ref_data: numpy.ndarray
    :param n_data: Array of total counts.
    :type n_data: numpy.ndarray
    :param gt_data: Array of genotype information for phasing.
    :type gt_data: numpy.ndarray
    :return: Sum of negative log-likelihood values.
    :rtype: float
    """
    # Compute likelihood using absolute difference between prob and genotype data
    phased_ll = opt_prob(np.abs(prob - gt_data), disp, int(ref_data.sum()), int(n_data.sum()))
    return np.sum(phased_ll)


def opt_unphased(prob: float,
                 first_data: np.ndarray,
                 phase_data: np.ndarray) -> float:
    """
    Optimize likelihood for unphased data.

    This function computes the negative log-likelihood without using phasing information.
    
    :param prob: Probability parameter.
    :type prob: float
    :param first_data: Data for the first SNP (dispersion, ref_count, total count).
    :type first_data: list or tuple of numbers (or numpy.ndarray of length 3)
    :param phase_data: Data for subsequent SNPs.
    :type phase_data: list or tuple of numbers (or numpy.ndarray with three rows)
    :return: Combined negative log-likelihood.
    :rtype: float
    """
    first_ll: float = opt_prob(prob, first_data[0], int(first_data[1]), int(first_data[2]))
    
    phase1_lls: float = opt_prob(prob, phase_data[0], int(phase_data[1]), int(phase_data[2]), log=False)
    phase2_lls: float = opt_prob(1 - prob, phase_data[0], int(phase_data[1]), int(phase_data[2]), log=False)

    combined_lls = (0.5 * phase1_lls) + (0.5 * phase2_lls)
    return first_ll - np.sum(np.log(combined_lls))


def opt_unphased_dp(prob: float,
                    disp: float,
                    first_ref: Union[np.ndarray, List[float]],
                    first_n: Union[np.ndarray, List[float]],
                    phase_ref: Union[np.ndarray, List[float]],
                    phase_n: Union[np.ndarray, List[float]]) -> float:
    """
    Optimize likelihood for unphased data using dynamic programming (DP).

    This function computes the negative log-likelihood by iteratively combining 
    likelihoods from the first SNP and the phased data for subsequent SNPs.

    :param prob: Probability parameter.
    :type prob: float
    :param disp: Dispersion parameter.
    :type disp: float
    :param first_ref: Reference allele count for the first SNP.
    :type first_ref: array-like (usually a one-element array)
    :param first_n: Total count for the first SNP.
    :type first_n: array-like (usually a one-element array)
    :param phase_ref: Array of reference allele counts for phased SNPs.
    :type phase_ref: array-like
    :param phase_n: Array of total counts for phased SNPs.
    :type phase_n: array-like
    :return: Combined negative log-likelihood.
    :rtype: float
    """
    first_ll: float = opt_prob(prob, disp, int(first_ref[0]), int(first_n[0]))
    phase1_like = opt_prob(prob, disp, int(phase_ref.sum()), int(phase_n.sum()), log=False)
    phase2_like = opt_prob(1 - prob, disp, int(phase_ref.sum()), int(phase_n.sum()), log=False)
    
    prev_like = 1.0
    for p1, p2 in zip(np.atleast_1d(phase1_like), np.atleast_1d(phase2_like)):
        p1_combined_like = prev_like * p1
        p2_combined_like = prev_like * p2
        prev_like = (0.5 * p1_combined_like) + (0.5 * p2_combined_like)

    return first_ll - np.log(prev_like)


def parse_opt(df: pd.DataFrame,
              disp: Optional[Union[float, np.ndarray]] = None,
              phased: bool = False) -> Tuple[float, float]:
    """
    Optimize model parameters and compute the imbalance likelihood.

    This function prepares the data from a dataframe with allele counts, computes
    the likelihood of the alternative model, and returns the estimated imbalance proportion.

    :param df: DataFrame containing allele counts and dispersion information.
    :type df: pandas.DataFrame
    :param disp: Pre-computed dispersion parameter; if None, the dispersion is taken from df["disp"].
    :type disp: float, optional
    :param phased: Boolean indicating if the data is phased.
    :type phased: bool, optional
    :return: Tuple containing the likelihood of the alternative model and the estimated imbalance proportion.
    :rtype: tuple (float, float)
    """
    snp_count: int = df.shape[0]
    ref_array: np.ndarray = df["ref_count"].to_numpy()
    n_array: np.ndarray = df["N"].to_numpy()

    if disp is None:
        disp = df["disp"].to_numpy()

    if snp_count > 1:
        if phased:
            gt_array: np.ndarray = df["GT"].to_numpy()
            if gt_array[0] > 0:
                gt_array = 1 - gt_array
            res = minimize_scalar(opt_phased_new,
                                  args=(disp, ref_array, n_array, gt_array),
                                  method="bounded", bounds=(0, 1))
        else:
            first_ref: np.ndarray = ref_array[:1]
            first_n: np.ndarray = n_array[:1]
            phase_ref: np.ndarray = ref_array[1:]
            phase_n: np.ndarray = n_array[1:]
            res = minimize_scalar(opt_unphased_dp,
                                  args=(disp, first_ref, first_n, phase_ref, phase_n),
                                  method="bounded", bounds=(0, 1))
    else:
        res = minimize_scalar(opt_prob,
                              args=(disp, int(ref_array[0]), int(n_array[0])),
                              method="bounded", bounds=(0, 1))

    mu: float = res["x"]
    alt_ll: float = -1 * res["fun"]

    return alt_ll, mu


def single_model(df: pd.DataFrame,
                 region_col: str,
                 phased: bool = False) -> pd.DataFrame:
    """
    Find allelic imbalance using a standard beta-binomial model.

    This function computes the imbalance likelihood for each group (or region)
    by comparing the observed reference counts to the expected counts under a 
    beta-binomial model with a single dispersion parameter.

    :param df: DataFrame with allele counts, including columns 'ref_count', 'N', and 'disp'.
    :type df: pandas.DataFrame
    :param region_col: The column name in df that defines the grouping (e.g., gene or peak).
    :type region_col: str
    :param phased: Boolean flag indicating whether the genotype data is phased.
    :type phased: bool, optional
    :return: A DataFrame with the following columns:
             - region (as defined by region_col)
             - null_ll: Likelihood under the null model
             - alt_ll: Likelihood under the alternative model
             - mu: Estimated imbalance proportion
             - lrt: Likelihood ratio test statistic
             - pval: p-value from the chi-square test
    :rtype: pandas.DataFrame
    """
    print("Running analysis with single dispersion model")
    opt_disp = lambda rho, ref_data, n_data: -np.sum(
        betabinom.logpmf(ref_data, n_data, (0.5 * (1 - rho) / rho),
                         (0.5 * (1 - rho) / rho))
    )
    
    ref_array: np.ndarray = df["ref_count"].to_numpy()
    n_array: np.ndarray = df["N"].to_numpy()

    disp_start = timeit.default_timer()
    disp: float = minimize_scalar(opt_disp, args=(ref_array, n_array),
                                  method="bounded", bounds=(0, 1))["x"]
    print(f"Optimized dispersion parameter in {timeit.default_timer() - disp_start:.2f} seconds")

    group_df = df.groupby(region_col, sort=False)

    print("Optimizing imbalance likelihood")
    ll_start = timeit.default_timer()
    null_test = group_df.apply(lambda x: np.sum(
        betabinom.logpmf(x["ref_count"].to_numpy(),
                          x["N"].to_numpy(),
                          (0.5 * (1 - disp) / disp),
                          (0.5 * (1 - disp) / disp)
                         )
    ))
    alt_test = group_df.apply(lambda x: parse_opt(x, disp, phased=phased))
    alt_df = pd.DataFrame(alt_test.to_list(), columns=["alt_ll", "mu"], index=alt_test.index)

    print(f"Optimized imbalance likelihood in {timeit.default_timer() - ll_start:.2f} seconds")

    ll_df = pd.concat([null_test, alt_df], axis=1).reset_index()
    ll_df.columns = [region_col, "null_ll", "alt_ll", "mu"]

    ll_df["lrt"] = -2 * (ll_df["null_ll"] - ll_df["alt_ll"])
    ll_df["pval"] = chi2.sf(ll_df["lrt"], 1)

    return ll_df


def linear_model(df: pd.DataFrame,
                 region_col: str,
                 phased: bool = False) -> pd.DataFrame:
    """
    Find allelic imbalance using a linear beta-binomial model.

    This function models the dispersion parameter as a linear function of total counts (N),
    and then computes the imbalance likelihood for each group (or region).

    :param df: DataFrame with allele counts, which must include 'ref_count', 'N', and a preliminary 'disp'.
    :type df: pandas.DataFrame
    :param region_col: The column name that defines grouping (e.g., gene or peak).
    :type region_col: str
    :param phased: Boolean flag indicating whether genotype data is phased.
    :type phased: bool, optional
    :return: A DataFrame with the following columns:
             - region (as defined by region_col)
             - null_ll: Log-likelihood under the null model
             - alt_ll: Log-likelihood under the alternative model
             - mu: Estimated imbalance proportion
             - lrt: Likelihood ratio test statistic
             - pval: p-value from the chi-square test
    :rtype: pandas.DataFrame
    """
    print("Running analysis with linear dispersion model")
    in_data: np.ndarray = df[["ref_count", "N"]].to_numpy().T
    
    print("Optimizing dispersion parameters...")
    disp_start = time.time()
    res = minimize(opt_linear, x0=(0, 0), method="Nelder-Mead", args=(in_data[0], in_data[1]))
    disp1, disp2 = res["x"]
    df["disp"] = expit(disp1 + (in_data[1] * disp2))
    print(f"Optimized dispersion parameters in {time.time() - disp_start} seconds")

    group_df = df.groupby(region_col, sort=False)
    
    print("Optimizing imbalance likelihood")
    ll_start = time.time()
    null_test = group_df.apply(lambda x: np.sum(
        betabinom.logpmf(x["ref_count"].to_numpy(),
                          x["N"].to_numpy(),
                          (0.5 * (1 - x["disp"].to_numpy()) / x["disp"].to_numpy()),
                          (0.5 * (1 - x["disp"].to_numpy()) / x["disp"].to_numpy())
                         )
    ))
    alt_test = group_df.apply(lambda x: parse_opt(x))
    alt_df = pd.DataFrame(alt_test.to_list(), columns=["alt_ll", "mu"], index=alt_test.index)
    print(f"Optimized imbalance likelihood in {time.time() - ll_start} seconds")
    
    ll_df = pd.concat([null_test, alt_df], axis=1).reset_index()
    ll_df.columns = [region_col, "null_ll", "alt_ll", "mu"]
    ll_df["lrt"] = -2 * (ll_df["null_ll"] - ll_df["alt_ll"])
    ll_df["pval"] = chi2.sf(ll_df["lrt"], 1)
    
    return ll_df


# def binom_model(df):
#     """
#     Find allelic imbalance using a standard binomial model
#
#     :param df: Dataframe with allele counts
#     :type df: DataFrame
#     :return: Dataframe with imbalance likelihood
#     :rtype: DataFrame
#     """
#
#     print("Running analysis with binomial model")
#     group_df = df.groupby("peak", sort=False)
#     
#     print(f"Calculating imbalance likelihood")
#     ll_start = time.time()
#     
#     # Get null test
#     null_test = group_df.apply(lambda x: np.sum(binom.logpmf(x["ref_count"].to_numpy(), x["N"].to_numpy(), 0.5)))
#     
#     # Optimize Alt
#     alt_test = group_df.apply(lambda x: binom_phase(x))
#
#     print(f"Calculated imbalance likelihood in {time.time() - ll_start} seconds")
#
#     ll_df = pd.concat([null_test, alt_test], axis=1).reset_index()
#     ll_df.columns = ["peak", "null_ll", "alt_ll"]
#     
#     ll_df["lrt"] = -2 * (ll_df["null_ll"] - ll_df["alt_ll"])
#     ll_df["pval"] = chi2.sf(ll_df["lrt"], 1)
#     
#     return ll_df


def bh_correction(df: pd.DataFrame) -> pd.DataFrame:
    """
    Apply Benjamini-Hochberg correction to control the false discovery rate.

    This function adjusts the p-values in the DataFrame using the Benjamini-Hochberg
    method, returning the DataFrame with an additional column 'fdr_pval'.

    :param df: DataFrame containing a 'pval' column of p-values.
    :type df: pandas.DataFrame
    :return: DataFrame with added 'fdr_pval' column.
    :rtype: pandas.DataFrame
    """
    if "pval" in df.columns:
        pcol = "pval"
    elif "pval" in df.columns[-1]:
        pcol = str(df.columns[-1])
    else:
        print("Pvalues not found! Returning Original Data")
        return df
    
    num_test = df.shape[0]

    if num_test == 1:
        df["fdr_pval"] = df[pcol]
        return df
    
    df["rank"] = rankdata(df[pcol], method="max").astype(int)
    df["adj_pval"] = df[pcol] * (num_test / df["rank"])
    
    rank_df = df[["rank", "adj_pval"]].drop_duplicates()
    rank_df = rank_df.sort_values(by=["rank"], ascending=False)

    rank_p = rank_df.set_index("rank").squeeze()
    rank_p = rank_p.rename("fdr_pval")
    rank_p[rank_p > 1] = 1
    
    prev = None
    for index, value in rank_p.items():
        if prev is None:
            prev = value
        elif value > prev:
            rank_p.at[index] = prev
        else:
            prev = value

    return_df = pd.merge(df, rank_p, left_on="rank", right_index=True).sort_index()
    return_df = return_df.drop(columns=["rank", "adj_pval"])

    return return_df


def get_imbalance(in_data: Union[pd.DataFrame, str],
                  min_count: int = 10,
                  pseudocount: int = 1,
                  method: str = "single",
                  phased: bool = False,
                  region_col: Optional[str] = None,
                  groupby: Optional[str] = None) -> pd.DataFrame:
    """
    Compute allelic imbalance for bulk data.

    This function reads allele count data (either from a DataFrame or a file), 
    applies a pseudocount and filters based on a minimum total count. It then 
    computes the likelihood of the alternative model using either the 'single' or 
    'linear' beta-binomial model, and returns a DataFrame with imbalance statistics.

    :param in_data: Either a DataFrame or a file path to a tab-separated file with allele counts.
    :type in_data: pandas.DataFrame or str
    :param min_count: Minimum total count (N) required for a SNP to be included in the analysis.
                      Defaults to 10.
    :type min_count: int, optional
    :param pseudocount: Value added to allele counts to avoid zero counts. Defaults to 1.
    :type pseudocount: int, optional
    :param method: Analysis method to use. Options are "single" or "linear". Defaults to "single".
    :type method: str, optional
    :param phased: Boolean flag indicating whether the genotype data is phased.
    :type phased: bool, optional
    :param region_col: Column name that specifies the grouping (e.g., gene or peak). If None, a unique variant identifier is created.
    :type region_col: str, optional
    :param groupby: Column name to use for grouping instead of the automatically generated variant identifier.
    :type groupby: str, optional
    :return: DataFrame with columns for region, null log-likelihood, alternative log-likelihood,
             estimated imbalance proportion (mu), likelihood ratio test (lrt), and p-values.
    :rtype: pandas.DataFrame
    """
    model_dict = {"single": single_model, "linear": linear_model}
    
    # Load data from file if necessary
    if isinstance(in_data, pd.DataFrame):
        df = in_data.copy()
    else:
        df = pd.read_csv(in_data,
                         sep="\t",
                         dtype={
                             "chrom": "category",
                             "pos": "uint32",
                             "ref": "category",
                             "alt": "category",
                             "ref_count": "uint16",
                             "alt_count": "uint16",
                             "other_count": "uint16"}
                        )

    if region_col is None:
        region_col = "variant"
        groupby = None
        df[region_col] = df["chrom"].astype("string") + "_" + df["pos"].astype("string")
    
    df[["ref_count", "alt_count"]] += pseudocount
    df["N"] = df["ref_count"] + df["alt_count"]
    df = df.loc[df["N"].ge(min_count + (2 * pseudocount)), :]

    if groupby is not None:
        region_col = groupby

    keep_cols = ["chrom", "pos", "ref_count", "alt_count", "N", region_col]
    
    if phased:
        if "GT" not in df.columns:
            print("Genotypes not found: Switching to unphased model")
            phased = False
        elif len(df["GT"].unique()) <= 1:
            print(f"All genotypes {df['GT'].unique()}: Switching to unphased model")
            phased = False
        elif not any(i in ['1|0', '0|1'] for i in df["GT"].unique()):
            print(f"Expected GT as 0|1 and 1|0 but found: {df['GT'].unique()}")
            print("Switching to unphased model")
            phased = False
        else:
            df["GT"] = df["GT"].str.split("|", n=1).str[0].astype("uint8")
            keep_cols.append("GT")
    
    df = df[keep_cols].drop_duplicates()
    p_df = model_dict[method](df, region_col, phased=phased)
    
    df[["ref_count", "alt_count"]] -= pseudocount
    df["N"] -= pseudocount * 2
    
    snp_counts = pd.DataFrame(df[region_col].value_counts(sort=False)).reset_index()
    snp_counts.columns = [region_col, "snp_count"]
    
    count_alleles = df[[region_col, "ref_count", "alt_count", "N"]].groupby(region_col, sort=False).sum()
    merge_df = pd.merge(snp_counts, p_df, how="left", on=region_col)
    
    as_df = pd.merge(count_alleles, merge_df, how="left", on=region_col)
    as_df["fdr_pval"] = false_discovery_control(as_df["pval"], method="bh")

    return as_df


def get_imbalance_sc(in_data: Union[pd.DataFrame, str],
                     min_count: int = 10,
                     method: str = "single",
                     out_dir: Optional[str] = None,
                     is_gene: bool = False,
                     feature: Optional[str] = None) -> pd.DataFrame:
    """
    Process input data and perform single-cell allelic imbalance analysis.

    This function processes a DataFrame (or reads from a file) containing single-cell 
    allele counts. It splits the data by cell, applies the specified imbalance model, 
    and performs Benjamini-Hochberg correction for multiple testing. Optionally, it 
    writes the results to disk.

    :param in_data: A pandas DataFrame or file path (TSV) containing allele count data.
    :type in_data: pandas.DataFrame or str
    :param min_count: Minimum total count (N) required for a SNP to be included (default: 10).
    :type min_count: int, optional
    :param method: Analysis method to use ("single" or "linear"). Defaults to "single".
    :type method: str, optional
    :param out_dir: Directory where output files should be written. If None, results are not saved to disk.
    :type out_dir: str, optional
    :param is_gene: Boolean flag. If True, treat gene names as regions.
    :type is_gene: bool, optional
    :param feature: Feature label to use in output file names. Defaults to "peak" if not provided.
    :type feature: str, optional
    :return: DataFrame with single-cell allelic imbalance results.
    :rtype: pandas.DataFrame
    """
    model_dict = {"single": single_model, "linear": linear_model}
    # model_dict = {"single": single_model, "linear": linear_model, "binomial": binom_model}

    if method not in model_dict:
        print("Please input a valid method (single, linear, binomial)")
        return pd.DataFrame()

    if isinstance(in_data, pd.DataFrame):
        df = in_data.copy()
    else:
        df = pd.read_csv(in_data, sep="\t")
    
    if is_gene:
        df = df.rename(columns={"genes": "peak"})

    default_df = df.iloc[:, :5]
    df_dict = {}
    start_index = min([df.columns.get_loc(c) for c in df.columns if "_ref" in c])
    for i in range(start_index, len(df.columns), 2):
        df_key = df.columns[i].split("_ref")[0]
        cell_df = pd.merge(default_df, df.iloc[:, [i, i+1]], left_index=True, right_index=True)
        cell_df.columns = ["chrom", "pos", "ref", "alt", "peak", "ref_count", "alt_count"]
        cell_df["N"] = cell_df["ref_count"] + cell_df["alt_count"]
        df_dict[df_key] = cell_df
    
    as_dict = {}
    return_df = df["peak"].drop_duplicates().reset_index(drop=True)
    fdr_df = df["peak"].drop_duplicates().reset_index(drop=True)
    
    for key, cell_df in df_dict.items():
        print(f"Analyzing imbalance for {key}")
        cell_df = cell_df.loc[cell_df["N"] >= min_count]
        if not cell_df.empty:
            p_df = model_dict[method](cell_df, "peak", phased=False)
            p_df = bh_correction(p_df)
            return_df = pd.merge(return_df, p_df[["peak", "pval"]], on="peak", how="left")
            return_df = return_df.rename(columns={"pval": f"{key}_pval"})
            fdr_df = pd.merge(fdr_df, p_df[["peak", "fdr_pval"]], on="peak", how="left")
            fdr_df = fdr_df.rename(columns={"fdr_pval": f"{key}_fdr"})
            snp_counts = pd.DataFrame(cell_df["peak"].value_counts(sort=False)).reset_index()
            snp_counts.columns = ["peak", "snp_count"]
            count_alleles = cell_df[["peak", "ref_count", "alt_count", "N"]].groupby("peak", sort=False).sum()
            merge_df = pd.merge(snp_counts, p_df, how="left", on="peak")
            as_df = pd.merge(count_alleles, merge_df, how="left", on="peak")
            as_dict[key] = as_df
        else:
            print(f"Not enough data to perform analysis on {key}")

    return_df = return_df.set_index("peak")
    return_df = return_df.dropna(axis=0, how="all").reset_index()
    fdr_df = fdr_df.set_index("peak")
    fdr_df = fdr_df.dropna(axis=0, how="all").reset_index()

    if is_gene:
        return_df = return_df.rename(columns={"peak": "genes"})
        fdr_df = fdr_df.rename(columns={"peak": "genes"})

    if feature is None:
        feature = "peak"

    if out_dir is not None:
        Path(out_dir).mkdir(parents=True, exist_ok=True)
        out_file = str(Path(out_dir) / f"as_results_{feature}_{method}_singlecell.tsv")
        return_df.to_csv(out_file, sep="\t", index=False)
        fdr_file = str(Path(out_dir) / f"as_results_{feature}_{method}_singlecell_fdr.tsv")
        fdr_df.to_csv(fdr_file, sep="\t", index=False)
        feat_dir = Path(out_dir) / f"cell_results_{feature}"
        feat_dir.mkdir(parents=True, exist_ok=True)
        for key, as_df in as_dict.items():
            if is_gene:
                as_df = as_df.rename(columns={"peak": "genes"})
            as_df.to_csv(str(feat_dir / f"{key}_results_{feature}_{method}.tsv"), sep="\t", index=False)
        print(f"Results written to {out_file}")

    return return_df


# LEGACY, NOT REALLY USED
def get_imbalance_sc(in_data: Union[pd.DataFrame, str],
                     min_count: int = 10,
                     method: str = "single",
                     out_dir: Optional[str] = None,
                     is_gene: bool = False,
                     feature: Optional[str] = None) -> pd.DataFrame:
    """
    Process input data and perform single-cell allelic imbalance analysis.

    :param in_data: DataFrame or file path (TSV) with allele counts.
    :type in_data: pandas.DataFrame or str
    :param min_count: Minimum total count required per SNP for inclusion (default: 10).
    :type min_count: int, optional
    :param method: Analysis method; options are "single" or "linear" (default: "single").
    :type method: str, optional
    :param out_dir: Directory where output files will be saved (if provided).
    :type out_dir: str, optional
    :param is_gene: If True, treat gene names as regions (default: False).
    :type is_gene: bool, optional
    :param feature: Feature label for output file names; defaults to "peak" if not provided.
    :type feature: str, optional
    :return: DataFrame with single-cell allelic imbalance results.
    :rtype: pandas.DataFrame
    """
    model_dict = {"single": single_model, "linear": linear_model}
    # model_dict = {"single": single_model, "linear": linear_model, "binomial": binom_model}

    if method not in model_dict:
        print("Please input a valid method (single, linear, binomial)")
        return pd.DataFrame()

    if isinstance(in_data, pd.DataFrame):
        df = in_data.copy()
    else:
        df = pd.read_csv(in_data, sep="\t")
    
    if is_gene:
        df = df.rename(columns={"genes": "peak"})

    default_df = df.iloc[:, :5]
    df_dict = {}
    start_index = min([df.columns.get_loc(c) for c in df.columns if "_ref" in c])
    for i in range(start_index, len(df.columns), 2):
        df_key = df.columns[i].split("_ref")[0]
        cell_df = pd.merge(default_df, df.iloc[:, [i, i+1]], left_index=True, right_index=True)
        cell_df.columns = ["chrom", "pos", "ref", "alt", "peak", "ref_count", "alt_count"]
        cell_df["N"] = cell_df["ref_count"] + cell_df["alt_count"]
        df_dict[df_key] = cell_df
    
    as_dict = {}
    return_df = df["peak"].drop_duplicates().reset_index(drop=True)
    fdr_df = df["peak"].drop_duplicates().reset_index(drop=True)
    
    for key, cell_df in df_dict.items():
        print(f"Analyzing imbalance for {key}")
        cell_df = cell_df.loc[cell_df["N"] >= min_count]
        if not cell_df.empty:
            p_df = model_dict[method](cell_df, "peak", phased=False)
            p_df = bh_correction(p_df)
            return_df = pd.merge(return_df, p_df[["peak", "pval"]], on="peak", how="left")
            return_df = return_df.rename(columns={"pval": f"{key}_pval"})
            fdr_df = pd.merge(fdr_df, p_df[["peak", "fdr_pval"]], on="peak", how="left")
            fdr_df = fdr_df.rename(columns={"fdr_pval": f"{key}_fdr"})
            snp_counts = pd.DataFrame(cell_df["peak"].value_counts(sort=False)).reset_index()
            snp_counts.columns = ["peak", "snp_count"]
            count_alleles = cell_df[["peak", "ref_count", "alt_count", "N"]].groupby("peak", sort=False).sum()
            merge_df = pd.merge(snp_counts, p_df, how="left", on="peak")
            as_df = pd.merge(count_alleles, merge_df, how="left", on="peak")
            as_dict[key] = as_df
        else:
            print(f"Not enough data to perform analysis on {key}")

    return_df = return_df.set_index("peak")
    return_df = return_df.dropna(axis=0, how="all").reset_index()
    fdr_df = fdr_df.set_index("peak")
    fdr_df = fdr_df.dropna(axis=0, how="all").reset_index()

    if is_gene:
        return_df = return_df.rename(columns={"peak": "genes"})
        fdr_df = fdr_df.rename(columns={"peak": "genes"})

    if feature is None:
        feature = "peak"

    if out_dir is not None:
        Path(out_dir).mkdir(parents=True, exist_ok=True)
        out_file = str(Path(out_dir) / f"as_results_{feature}_{method}_singlecell.tsv")
        return_df.to_csv(out_file, sep="\t", index=False)
        fdr_file = str(Path(out_dir) / f"as_results_{feature}_{method}_singlecell_fdr.tsv")
        fdr_df.to_csv(fdr_file, sep="\t", index=False)
        feat_dir = Path(out_dir) / f"cell_results_{feature}"
        feat_dir.mkdir(parents=True, exist_ok=True)
        for key, as_df in as_dict.items():
            if is_gene:
                as_df = as_df.rename(columns={"peak": "genes"})
            as_df.to_csv(str(feat_dir / f"{key}_results_{feature}_{method}.tsv"), sep="\t", index=False)
        print(f"Results written to {out_file}")

    return return_df
