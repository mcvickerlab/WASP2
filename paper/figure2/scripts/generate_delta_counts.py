#!/usr/bin/env python3
"""
Generate Figure 2 Panel C delta-counts table: per-site Δref and Δalt counts
(processed − original) for each method (WASP2, GATK, phASER).

Key properties:
- Uses the same site set across methods (intersection).
- Applies a consistent per-site coverage cutoff based on orig_total and filt_total.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def compute_delta_table(before_after: pd.DataFrame, min_total: int) -> tuple[pd.DataFrame, dict]:
    required_cols = {
        "chrom",
        "pos",
        "method",
        "orig_ref",
        "orig_alt",
        "filt_ref",
        "filt_alt",
        "orig_total",
        "filt_total",
    }
    missing = required_cols - set(before_after.columns)
    if missing:
        raise ValueError(f"before_after_counts.tsv missing columns: {sorted(missing)}")

    methods = ["WASP2", "GATK", "phASER"]
    df = before_after[before_after["method"].isin(methods)].copy()

    if min_total > 0:
        df = df[(df["orig_total"] >= min_total) & (df["filt_total"] >= min_total)]

    # Find intersection of sites present for all methods after filtering.
    present = (
        df.groupby(["chrom", "pos"])["method"]
        .nunique()
        .reset_index(name="n_methods")
    )
    present = present[present["n_methods"] == len(methods)][["chrom", "pos"]]
    df = df.merge(present, on=["chrom", "pos"], how="inner")

    df["delta_ref"] = df["filt_ref"] - df["orig_ref"]
    df["delta_alt"] = df["filt_alt"] - df["orig_alt"]

    out = df[
        [
            "chrom",
            "pos",
            "method",
            "delta_ref",
            "delta_alt",
            "orig_total",
            "filt_total",
        ]
    ].copy()

    stats: dict[str, object] = {
        "min_total": int(min_total),
        "n_sites_intersection": int(len(present)),
    }
    for m in methods:
        d = out[out["method"] == m]
        if len(d) == 0:
            stats[f"{m}_n"] = 0
            continue
        stats[f"{m}_n"] = int(len(d))
        stats[f"{m}_delta_ref_mean"] = float(d["delta_ref"].mean())
        stats[f"{m}_delta_ref_median"] = float(d["delta_ref"].median())
        stats[f"{m}_delta_alt_mean"] = float(d["delta_alt"].mean())
        stats[f"{m}_delta_alt_median"] = float(d["delta_alt"].median())
        for col in ["delta_ref", "delta_alt"]:
            for q in [0.01, 0.05, 0.95, 0.99]:
                stats[f"{m}_{col}_p{int(q*100):02d}"] = float(np.quantile(d[col].to_numpy(), q))

    return out, stats


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", choices=["hg00731"], default="hg00731")
    ap.add_argument("--min-total", type=int, default=10)
    args = ap.parse_args()

    repo_root = Path(
        "/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp"
    )
    data_dir = repo_root / "paper" / "figure2" / "data" / args.dataset
    in_tsv = data_dir / "before_after_counts.tsv"
    out_tsv = data_dir / "delta_counts.tsv"
    out_stats = data_dir / "delta_counts_stats.txt"

    if not in_tsv.exists():
        raise SystemExit(f"Missing input: {in_tsv} (run generate_before_after_counts.py first)")

    df = pd.read_csv(in_tsv, sep="\t")
    out, stats = compute_delta_table(df, min_total=args.min_total)

    out.to_csv(out_tsv, sep="\t", index=False)
    with out_stats.open("w") as f:
        for k, v in stats.items():
            f.write(f"{k}: {v}\n")

    print(f"Wrote: {out_tsv}")
    print(f"Wrote: {out_stats}")


if __name__ == "__main__":
    main()

