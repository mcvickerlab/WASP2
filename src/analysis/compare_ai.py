
# Default package imports
import timeit
from pathlib import Path
from itertools import combinations
from csv import DictReader

# External package imports
import pandas as pd
import numpy as np
from scipy.stats import betabinom, chi2
from scipy.optimize import minimize_scalar

# Local package imports
from as_analysis import opt_prob, opt_phase


def opt_compare_null(df, cell_a, cell_b):

    snp_count = df.shape[0]

    if snp_count > 1:
        phase_func = lambda prob, af, ap, bf, bp: opt_phase(prob, af, ap) + opt_phase(prob, bf, bp)

        # Make data for first snp and phased positions
        a_snp = df[:1][[f"{cell_a}_disp", f"{cell_a}_ref", f"{cell_a}_N"]].to_numpy()[0]
        a_phase = df[1:][[f"{cell_a}_disp", f"{cell_a}_ref", f"{cell_a}_N"]].to_numpy().T

        b_snp = df[:1][[f"{cell_b}_disp", f"{cell_b}_ref", f"{cell_b}_N"]].to_numpy()[0]
        b_phase = df[1:][[f"{cell_b}_disp", f"{cell_b}_ref", f"{cell_b}_N"]].to_numpy().T

        # Optimize joint likelihoods
        res = minimize_scalar(phase_func, args=(a_snp, a_phase, b_snp, b_phase), method="bounded", bounds=(0, 1))
        mu = res["x"]
        alt_ll = -1 * res["fun"]

    else:
        prob_func = lambda prob, a, b: opt_prob(prob, a[0], a[1], a[2]) + opt_prob(prob, b[0], b[1], b[2])

        # Make data for snp counts
        a_snp = df[[f"{cell_a}_disp", f"{cell_a}_ref", f"{cell_a}_N"]].to_numpy()[0]
        b_snp = df[[f"{cell_b}_disp", f"{cell_b}_ref", f"{cell_b}_N"]].to_numpy()[0]

        # Optimize joint likelihoods
        res = minimize_scalar(prob_func, args=(a_snp, b_snp), method="bounded", bounds=(0, 1))
        mu = res["x"]
        alt_ll = -1 * res["fun"]

    return alt_ll, mu


def opt_compare_alt(df, cell_a, cell_b):

    snp_count = df.shape[0]

    if snp_count > 1:
        # Make data for first snp and phased positions
        a_snp = df[:1][[f"{cell_a}_disp", f"{cell_a}_ref", f"{cell_a}_N"]].to_numpy()[0]
        a_phase = df[1:][[f"{cell_a}_disp", f"{cell_a}_ref", f"{cell_a}_N"]].to_numpy().T

        b_snp = df[:1][[f"{cell_b}_disp", f"{cell_b}_ref", f"{cell_b}_N"]].to_numpy()[0]
        b_phase = df[1:][[f"{cell_b}_disp", f"{cell_b}_ref", f"{cell_b}_N"]].to_numpy().T
        
        # Optimize first likelihood
        res_a = minimize_scalar(opt_phase, args=(a_snp, a_phase), method="bounded", bounds=(0, 1))
        mu_a = res_a["x"]

        # Optimize second likelihood
        res_b = minimize_scalar(opt_phase, args=(b_snp, b_phase), method="bounded", bounds=(0, 1))
        mu_b = res_b["x"]

        # Get joint log-likelihood
        alt_ll = -1 * (res_a["fun"] + res_b["fun"])

    else:
        # Make data for snp counts
        a_snp = df[[f"{cell_a}_disp", f"{cell_a}_ref", f"{cell_a}_N"]].to_numpy()[0]
        b_snp = df[[f"{cell_b}_disp", f"{cell_b}_ref", f"{cell_b}_N"]].to_numpy()[0]

        # Optimize first likelihood
        res_a = minimize_scalar(opt_prob, args=(a_snp[0], a_snp[1], a_snp[2]), method="bounded", bounds=(0, 1))
        mu_a = res_a["x"]

        # Optimize second likelihood
        res_b = minimize_scalar(opt_prob, args=(b_snp[0], b_snp[1], b_snp[2]), method="bounded", bounds=(0, 1))
        mu_b = res_b["x"]

        # Get joint log-likelihood
        alt_ll = -1 * (res_a["fun"] + res_b["fun"])

    return alt_ll, mu_a, mu_b


def compare_imbalance(df, cell_a, cell_b, show_counts=True):
    total_start = timeit.default_timer()

    # Optimize dispersion parameter using Single Model
    opt_disp = lambda rho, data: -np.sum(betabinom.logpmf(data[0], data[1], (0.5 * (1 - rho) / rho), (0.5 * (1 - rho) / rho)))
    
    a_data = df[[f"{cell_a}_ref", f"{cell_a}_N"]].to_numpy().T
    b_data = df[[f"{cell_b}_ref", f"{cell_b}_N"]].to_numpy().T

    df[f"{cell_a}_disp"] = minimize_scalar(opt_disp, args=(a_data), method="bounded", bounds=(0,1))["x"]
    df[f"{cell_b}_disp"] = minimize_scalar(opt_disp, args=(b_data), method="bounded", bounds=(0,1))["x"]

    # Group dataframe by regions
    group_df = df.groupby("peak", sort=False)

    # Optimize likelihoods assuming AI proportion is the same
    null_start = timeit.default_timer()

    null_test = group_df.apply(lambda x: opt_compare_null(x, cell_a, cell_b))
    null_df = pd.DataFrame(null_test.to_list(), columns=["null_ll", "shared_mu"], index=null_test.index)

    print(f"Optimized shared likelihood in {timeit.default_timer() - null_start} seconds")

    # Optimize likelihoods assuming AI proportion is different
    alt_start = timeit.default_timer()

    alt_test = group_df.apply(lambda x: opt_compare_alt(x, cell_a, cell_b))
    alt_df = pd.DataFrame(alt_test.to_list(),
                          columns=["alt_ll", f"{cell_a}_mu", f"{cell_b}_mu"], index=alt_test.index)

    print(f"Optimized different likelihoods in {timeit.default_timer() - alt_start} seconds")

    ll_df = pd.concat([null_df, alt_df], axis=1)

    ll_df["lrt"] = -2 * (ll_df["null_ll"] - ll_df["alt_ll"])
    ll_df["pval"] = chi2.sf(ll_df["lrt"], 1)
    
    ll_df = ll_df[["shared_mu", f"{cell_a}_mu", f"{cell_b}_mu",
                   "null_ll", "alt_ll", "lrt", "pval"]]

    runtime = timeit.default_timer() - total_start
    print(f"Finished in {runtime} seconds!\n")

    if show_counts:
        snp_counts = df["peak"].value_counts(sort=False) # get individual counts
        snp_counts.name = "snp_count"

        col_order = ["peak", f"{cell_a}_ref", f"{cell_a}_alt", f"{cell_a}_N",
                     f"{cell_b}_ref", f"{cell_b}_alt", f"{cell_b}_N"]

        count_alleles = df[col_order].groupby("peak", sort=False).sum()

        col_order[0] = "snp_count"

        merge_df = pd.concat([count_alleles, snp_counts], axis=1)
        merge_df = merge_df[col_order]

        ll_df = pd.concat([merge_df, ll_df], axis=1)

    return ll_df


def merge_samp_counts(s1_df, s2_df, cell, min_count=10):
    comp_cols = ["chrom", "pos", "ref", "alt", "peak", f"{cell}_ref", f"{cell}_alt"]

    s1_counts = s1_df.loc[:, comp_cols]
    s1_counts = s1_counts.rename(columns={f"{cell}_ref": "S1_ref", f"{cell}_alt": "S1_alt"})
    s1_counts["S1_N"] = s1_counts["S1_ref"] + s1_counts["S1_alt"]

    s2_counts = s2_df.loc[:, comp_cols]
    s2_counts = s2_counts.rename(columns={f"{cell}_ref": "S2_ref", f"{cell}_alt": "S2_alt"})
    s2_counts["S2_N"] = s2_counts["S2_ref"] + s2_counts["S2_alt"]

    # TODO: ALLOW FOR DIFFERENT SNP'S IF PEAKS MEET MINS
    c_df = pd.merge(s1_counts, s2_counts, how="inner", on=["chrom", "pos", "ref", "alt", "peak"])

    c_df = c_df.loc[(c_df["S1_N"] >= min_count) & (c_df["S2_N"] >= min_count)].reset_index(drop=True)
    
    return c_df


def sample_comp(count_files, cells, min_count=10, is_gene=False):

    with open(count_files[1], "r") as file:
        header = DictReader(file, delimiter="\t").fieldnames

    other_cells = [col.split("_ref")[0] for col in header if col.endswith("_ref")]
    cell_list = list(set(other_cells).intersection(cells))

    s1_df = pd.read_csv(count_files[0], sep="\t")
    s2_df = pd.read_csv(count_files[1], sep="\t")

    # Change label for gene to peak temporarily
    if is_gene:
        s1_df = s1_df.rename(columns={"genes": "peak"})
        s2_df = s2_df.rename(columns={"genes": "peak"})

    s1_idx = pd.Index(s1_df["peak"].drop_duplicates())
    s2_idx = pd.Index(s2_df["peak"].drop_duplicates())

    return_df = pd.DataFrame(index=s1_idx.union(s2_idx))
    as_dict = {}

    for cell in cell_list:
        c_df = merge_samp_counts(s1_df, s2_df, cell, min_count=min_count)

        if not c_df.empty:
            print(f"Comparing allelic imbalance in {cell} between samples")
            diff_df = compare_imbalance(c_df, "S1", "S2")

            return_df = pd.concat([return_df, diff_df["pval"]], axis=1)
            return_df = return_df.rename(columns={"pval": f"{cell}_pval"})

            diff_df = diff_df.reset_index()
            if is_gene:
                diff_df = diff_df.rename(columns={"index": "genes"})
            else:
                diff_df = diff_df.rename(columns={"index": "peak"})

            as_dict[cell] = diff_df
        else:
            print(f"Not enough {cell} data to compare samples\n")
    
    return return_df, as_dict


def group_comp(count_file, cells, min_count=10, is_gene=False):
    df = pd.read_csv(count_file, sep="\t")

    # Change label for gene to peak temporarily
    if is_gene:
        df = df.rename(columns={"genes": "peak"})

    as_dict = {}
    return_df = pd.DataFrame(index=df["peak"].drop_duplicates())

    for cell_a, cell_b in list(combinations(cells, 2)):
        comp_cols = ["peak", f"{cell_a}_ref", f"{cell_a}_alt", f"{cell_b}_ref", f"{cell_b}_alt"]
        c_df = df.loc[:, comp_cols]

        # Parse allele minimums
        c_df[f"{cell_a}_N"] = c_df[f"{cell_a}_ref"] + c_df[f"{cell_a}_alt"]
        c_df[f"{cell_b}_N"] = c_df[f"{cell_b}_ref"] + c_df[f"{cell_b}_alt"]

        c_df = c_df.loc[(c_df[f"{cell_a}_N"] >= min_count)
                        & (c_df[f"{cell_b}_N"] >= min_count)].reset_index(drop=True)

        if not c_df.empty:
            key=f"{cell_a}_{cell_b}"

            print(f"Comparing allelic imbalance in {cell_a} and {cell_b}")
            diff_df = compare_imbalance(c_df, cell_a, cell_b)

            return_df = pd.concat([return_df, diff_df["pval"]], axis=1)
            return_df = return_df.rename(columns={"pval": f"{key}_pval"})

            diff_df = diff_df.reset_index()
            if is_gene:
                diff_df = diff_df.rename(columns={"index": "genes"})
            else:
                diff_df = diff_df.rename(columns={"index": "peak"})

            as_dict[key] = diff_df
        else:
            print(f"Not enough data to compare {cell_a} and {cell_b}\n")

    return return_df, as_dict


def get_imbalance(count_files, cells, min_count=10, out_dir=None, is_gene=False):
    multi = True if len(count_files) > 1 else False

    # if cells not provided, get all groups:
    if not cells:
        print("Comparing all cells/clusters")
        with open(count_files[0], "r") as file:
            header = DictReader(file, delimiter="\t").fieldnames

        cells = [col.split("_ref")[0] for col in header if col.endswith("_ref")]

    if multi:
        return_df, as_dict = sample_comp(count_files, cells,
                                         min_count=min_count,is_gene=is_gene)
    else:
        return_df, as_dict = group_comp(count_files[0], cells,
                                        min_count=min_count,is_gene=is_gene)

    return_df = return_df.dropna(axis=0, how="all")
    return_df = return_df.reset_index()

    if is_gene:
        return_df = return_df.rename(columns={"index": "genes"})
    else:
        return_df = return_df.rename(columns={"index": "peak"})

    if out_dir is not None:
        if multi:
            out_file = Path(out_dir) / "ai_compare_samples.tsv"
            cell_dir = Path(out_dir) / "sample_results"
        else:
            out_file = Path(out_dir) / "ai_compare_groups.tsv"
            cell_dir = Path(out_dir) / "group_results"

        return_df.to_csv(str(out_file), sep="\t", index=False)
        cell_dir.mkdir(parents=True, exist_ok=True)

        for key, as_df in as_dict.items():
            cell_file = cell_dir / f"{key}_results.tsv"
            as_df.to_csv(str(cell_file), sep="\t", index=False)

        print(f"Results written to {out_file}")
    
    return return_df

