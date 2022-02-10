"""
Author: Aaron Ho
Python Version: 3.9
"""

# Default Python package Imports
from pathlib import Path
import time

# External package imports
import pandas as pd
import numpy as np
from scipy.stats import betabinom, chi2, binom, rankdata
from scipy.optimize import minimize_scalar, minimize
from scipy.special import expit


def opt_linear(disp_params, ref_counts, n_array):
    """
    Optimize dispersion parameter weighted by N
    (Function called by optimizer)
    """
    disp1, disp2 = disp_params

    exp_in = (disp1 + (n_array * disp2))
    exp_in = np.select([exp_in > 10, exp_in < -10], [10, -10], default=exp_in)

    rho = expit(exp_in)

    ll = -np.sum(betabinom.logpmf(ref_counts, n_array, (0.5 * (1 - rho) / rho), (0.5 * (1 - rho) / rho))) # If alpha is beta
    
    return ll


def opt_prob(in_prob, in_rho, k, n, log=True):
    """
    Optimize Probability value that maximizes imbalance likelihood.
    (Function called by optimizer)
    """
    prob = in_prob

    alpha = (prob * (1 - in_rho) / in_rho)
    beta = ((1 - prob) * (1 - in_rho) / in_rho)
    
    if log is True:
        ll = -1 * betabinom.logpmf(k, n, alpha, beta)
    else:
        ll = betabinom.pmf(k, n, alpha, beta)

    return ll


def opt_phase(prob, first_data, phase_data):
    """
    Optimize likelihood while taking phase into account
    (Function called by optimizer)
    """
    
    first_ll = opt_prob(prob, first_data[0], first_data[1], first_data[2])
    
    # Sum opts given prob
    phase1_lls = opt_prob(prob, phase_data[0], phase_data[1], phase_data[2], log=False)
    phase2_lls = opt_prob(1 - prob, phase_data[0], phase_data[1], phase_data[2], log=False)


    combined_lls = (0.5 * phase1_lls) + (0.5 * phase2_lls)
    return first_ll + -np.sum(np.log(combined_lls))


def parse_opt(df, in_disp=None):
    """
    Optimize necessary data when running model

    :param df: Dataframe with allele counts
    :type df: DataFrame
    :param in_disp: pre-computed dispersion parameter, defaults to None
    :type in_disp: float, optional
    :return: Liklihood of alternate model, and imbalance proportion
    :rtype: array, array
    """

    snp_count = df.shape[0]
    
    if in_disp is not None:
        df["disp"] = in_disp

    if snp_count > 1:
        first_data = df[:1][["disp", "ref_count", "N"]].to_numpy()[0]
        phase_data = df[1:][["disp", "ref_count", "N"]].to_numpy().T
        
        res = minimize_scalar(opt_phase, args=(first_data, phase_data), method="bounded", bounds=(0, 1))
        mu = res["x"]
        alt_ll = -1 * res["fun"]

    else:
        snp_data = df[["disp", "ref_count", "N"]].to_numpy()[0]

        res = minimize_scalar(opt_prob, args=(snp_data[0], snp_data[1], snp_data[2]), method="bounded", bounds=(0, 1))
        mu = res["x"]
        alt_ll = -1 * res["fun"]

    return alt_ll, mu


def binom_phase(df):
    """
    Process phasing data for binomial model

    :param df: Dataframe with allele counts
    :type df: DataFrame
    :return: Array with alternate hypothesis likelihoods
    :rtype: Array
    """
    
    snp_count = df.shape[0]
    
    if snp_count > 1:
        first_data = df[:1][["ref_count", "N"]].to_numpy()[0]
        first_ll = binom.logpmf(first_data[0], first_data[1], (first_data[0]/first_data[1]))


        phase_data = df[1:][["ref_count", "N"]].to_numpy().T

        phase1_lls = binom.pmf(phase_data[0], phase_data[1], (phase_data[0]/phase_data[1]))
        phase2_lls = binom.pmf(phase_data[0], phase_data[1], 1 - (phase_data[0]/phase_data[1]))

        combined_lls = (0.5 * phase1_lls) + (0.5 * phase2_lls)

        alt_ll = first_ll + np.sum(np.log(combined_lls))

    else:
        snp_data = df[["ref_count", "N"]].to_numpy()[0]
        alt_ll = binom.logpmf(snp_data[0], snp_data[1], (snp_data[0]/snp_data[1]))

    return alt_ll


def single_model(df):
    """
    Find allelic imbalance using normal beta-binomial model

    :param df: Dataframe with allele counts
    :type df: DataFrame
    :return: Dataframe with imbalance likelihood
    :rtype: DataFrame
    """

    print("Running analysis with single dispersion model")

    opt_disp = lambda rho, data: -np.sum(betabinom.logpmf(data[0], data[1], (0.5 * (1 - rho) / rho), (0.5 * (1 - rho) / rho)))
    in_data = df[["ref_count", "N"]].to_numpy().T


    print("Optimizing dispersion parameter...")
    disp_start = time.time()

    disp = minimize_scalar(opt_disp, args=(in_data),
                           method="bounded", bounds=(0,1))["x"]

    print(f"Optimized dispersion parameter in {time.time() - disp_start} seconds")


    group_df = df.groupby("peak", sort=False)

    print("Optimizing imbalance likelihood")
    ll_start = time.time()
    null_test = group_df.apply(lambda x: np.sum(betabinom.logpmf(x["ref_count"].to_numpy(), x["N"].to_numpy(),
                                                                 (0.5 * (1 - disp) / disp), (0.5 * (1 - disp) / disp))))

    # Optimize Alt
    alt_test = group_df.apply(lambda x: parse_opt(x, disp))
    alt_df = pd.DataFrame(alt_test.to_list(), columns=["alt_ll", "mu"], index=alt_test.index)

    print(f"Optimized imbalance likelihood in {time.time() - ll_start} seconds")

    ll_df = pd.concat([null_test, alt_df], axis=1).reset_index()
    ll_df.columns = ["peak", "null_ll", "alt_ll", "mu"]

    ll_df["lrt"] = -2 * (ll_df["null_ll"] - ll_df["alt_ll"])
    ll_df["pval"] = chi2.sf(ll_df["lrt"], 1)

    return ll_df


def linear_model(df):
    """
    Find allelic imbalance using linear allelic imbalance model,
    weighting imbalance linear with N counts

    :param df: Dataframe with allele counts
    :type df: DataFrame
    :return: Dataframe with imbalance likelihood
    :rtype: DataFrame
    """

    print("Running analysis with linear dispersion model")
    in_data = df[["ref_count", "N"]].to_numpy().T
    
    print("Optimizing dispersion parameters...")
    disp_start = time.time()

    res = minimize(opt_linear, x0=(0, 0), method="Nelder-Mead", args=(in_data[0], in_data[1]))
    disp1, disp2 = res["x"]
    df["disp"] = expit((disp1 + (in_data[1] * disp2)))

    print(f"Optimized dispersion parameters in {time.time() - disp_start} seconds")

    # Group by region
    group_df = df.groupby("peak", sort=False)

    # Get null test
    print("Optimizing imbalance likelihood")
    ll_start = time.time()
    null_test = group_df.apply(lambda x: np.sum(betabinom.logpmf(
        x["ref_count"].to_numpy(), x["N"].to_numpy(),
        (0.5 * (1 - x["disp"].to_numpy()) / x["disp"].to_numpy()),
        (0.5 * (1 - x["disp"].to_numpy()) / x["disp"].to_numpy()))))

    # Optimize Alt
    alt_test = group_df.apply(lambda x: parse_opt(x))
    alt_df = pd.DataFrame(alt_test.to_list(), columns=["alt_ll", "mu"], index=alt_test.index)
    
    print(f"Optimized imbalance likelihood in {time.time() - ll_start} seconds")
    
    ll_df = pd.concat([null_test, alt_df], axis=1).reset_index()
    ll_df.columns = ["peak", "null_ll", "alt_ll", "mu"]

    ll_df["lrt"] = -2 * (ll_df["null_ll"] - ll_df["alt_ll"])
    ll_df["pval"] = chi2.sf(ll_df["lrt"], 1)

    return ll_df


def binom_model(df):
    """
    Find allelic imbalance using a standard binomial model

    :param df: Dataframe with allele counts
    :type df: DataFrame
    :return: Dataframe with imbalance likelihood
    :rtype: DataFrame
    """

    print("Running analysis with binomial model")
    group_df = df.groupby("peak", sort=False)
    
    print(f"Calculating imbalance likelihood")
    ll_start = time.time()
    
    # Get null test
    null_test = group_df.apply(lambda x: np.sum(binom.logpmf(x["ref_count"].to_numpy(), x["N"].to_numpy(), 0.5)))
    
    # Optimize Alt
    alt_test = group_df.apply(lambda x: binom_phase(x))

    print(f"Calculated imbalance likelihood in {time.time() - ll_start} seconds")

    ll_df = pd.concat([null_test, alt_test], axis=1).reset_index()
    ll_df.columns = ["peak", "null_ll", "alt_ll"]
    
    ll_df["lrt"] = -2 * (ll_df["null_ll"] - ll_df["alt_ll"])
    ll_df["pval"] = chi2.sf(ll_df["lrt"], 1)
    
    return ll_df


def bh_correction(df):
    if "pval" in df.columns:
        pcol = "pval"
    elif "pval" in df.columns[-1]:
        pcol = str(df.columns[-1])
    else:
        print("Pvalues not found! Returning Original Data")
        return df
    
    df["rank"] = rankdata(df[pcol], method="max").astype(int)

    num_test = df.shape[0]
    df["adj_pval"] = df[pcol] * (num_test / df["rank"])
    
    rank_df = df[["rank", "adj_pval"]].drop_duplicates()
    rank_df = rank_df.sort_values(by=["rank"], ascending=False)

    rank_p = rank_df.set_index("rank").squeeze()
    rank_p = rank_p.rename("fdr_pval")
    rank_p[rank_p > 1] = 1
    
    # test_adj
    prev = None
    for index, value in rank_p.items():
        if prev is None:
            prev = value
        elif value > prev:
            rank_p.at[index] = prev
        else:
            prev = value

    # Combine back into df
    return_df = pd.merge(df, rank_p, left_on="rank", right_index=True).sort_index()
    return_df = return_df.drop(columns=["rank", "adj_pval"])

    return return_df


def get_imbalance(in_data, min_count=10, method="single", out_dir=None, is_gene=False, feature=None):
    """
    Process input data and method for finding allelic imbalance

    :param in_data: Dataframe with allele counts
    :type in_data: DataFrame
    :param min_count: minimum allele count for analysis, defaults to 10
    :type min_count: int, optional
    :param method: analysis method, defaults to "single"
    :type method: str, optional
    :param out_dir: output directory, defaults to None
    :type out_dir: str, optional
    :return: DataFrame with imbalance Pvals per region
    :rtype: DataFrame
    """

    model_dict = {"single": single_model, "linear": linear_model, "binomial": binom_model}
    
    if method not in model_dict:
        print("Please input a valid method (single, linear, binomial)")
        return -1
    
    if isinstance(in_data, pd.DataFrame):
        df = in_data
    else:
        df = pd.read_csv(in_data, sep="\t")

    # Change label for gene to peak temporarily
    if is_gene is True:
        df = df.rename(columns={"genes": "peak"})

    
    df["N"] = df["ref_count"] + df["alt_count"]
    df = df.loc[df["N"] >= min_count]
    
    p_df = model_dict[method](df) # Perform analysis
    
    snp_counts = pd.DataFrame(df["peak"].value_counts(sort=False)).reset_index()
    snp_counts.columns = ["peak", "snp_count"]
    
    
    count_alleles = df[["peak", "ref_count", "alt_count", "N"]].groupby("peak", sort=False).sum()
    
    merge_df = pd.merge(snp_counts, p_df, how="left", on="peak")
    
    as_df = pd.merge(count_alleles, merge_df, how="left", on="peak")
    as_df = bh_correction(as_df)

    # Change label for gene to peak temporarily
    if is_gene is True:
        as_df = as_df.rename(columns={"peak": "genes"})
    
    if feature is None:
        feature = "peak"

    if out_dir is not None:
        Path(out_dir).mkdir(parents=True, exist_ok=True)
        out_file = str(Path(out_dir) / f"as_results_{feature}_{method}.tsv")
        as_df.to_csv(out_file, sep="\t", index=False)
        print(f"Results written to {out_file}")

    return as_df


def get_imbalance_sc(in_data, min_count=10, method="single", out_dir=None, is_gene=False, feature=None):
    """
    Process input data and method for finding single-cell allelic imbalance

    :param in_data: Dataframe with allele counts
    :type in_data: DataFrame
    :param min_count: minimum allele count for analysis, defaults to 10
    :type min_count: int, optional
    :param method: analysis method, defaults to "single"
    :type method: str, optional
    :param out: output directory, defaults to None
    :type out: str, optional
    :return: DataFrame with imbalance Pvals per region and per cell type
    :rtype: DataFrame
    """

    model_dict = {"single": single_model, "linear": linear_model, "binomial": binom_model}

    if method not in model_dict:
        print("Please input a valid method (single, linear, binomial)")
        return -1

    if isinstance(in_data, pd.DataFrame):
        df = in_data
    else:
        df = pd.read_csv(in_data, sep="\t")
    
    # Change label for gene to peak temporarily
    if is_gene is True:
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
        
        cell_df = cell_df.loc[cell_df["N"] >= min_count] # Filter by N
        
        if not cell_df.empty:
            p_df = model_dict[method](cell_df)
            p_df = bh_correction(p_df)

            return_df = pd.merge(return_df, p_df[["peak", "pval"]], on="peak", how="left")
            return_df = return_df.rename(columns={"pval": f"{key}_pval"})

            fdr_df = pd.merge(fdr_df, p_df[["peak", "fdr_pval"]], on="peak", how="left")
            fdr_df = fdr_df.rename(columns={"fdr_pval": f"{key}_fdr"})
            
            snp_counts = pd.DataFrame(cell_df["peak"].value_counts(sort=False)).reset_index() # get individual counts
            snp_counts.columns = ["peak", "snp_count"]
            
            count_alleles = cell_df[["peak", "ref_count", "alt_count", "N"]].groupby("peak", sort=False).sum()
            merge_df = pd.merge(snp_counts, p_df, how="left", on="peak")
            
            as_df = pd.merge(count_alleles, merge_df, how="left", on="peak")
            as_dict[key] = as_df

        else:
            print(f"Not enough data to perform analysis on {key}")

    # Remove empty columns
    return_df = return_df.set_index("peak")
    return_df = return_df.dropna(axis=0, how="all").reset_index()

    fdr_df = fdr_df.set_index("peak")
    fdr_df = fdr_df.dropna(axis=0, how="all").reset_index()

    if is_gene is True:
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
            
            if is_gene is True:
                as_df = as_df.rename(columns={"peak": "genes"})

            as_df.to_csv(str(feat_dir / f"{key}_results_{feature}_{method}.tsv"), sep="\t", index=False)
        
        print(f"Results written to {out_file}")

    return return_df
