#!/usr/bin/env python3
"""
Supplement for Figure 2: compare WASP1 vs WASP2 read filtering effects using a single counter.

We use the WASP2-rust BamCounter outputs on:
  - original BAM (pre-WASP)
  - WASP2-filtered BAM (post-WASP2)
  - WASP1-filtered BAM (post-WASP1)  [optional; requires wasp1_counts.filtered.tsv]

This avoids conflating "counting tool differences" (GATK/phASER) with "filtering pipeline differences".
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import sys

_this_file = globals().get("__file__")
if _this_file:
    _paper_dir = Path(_this_file).resolve().parents[2]
else:
    _paper_dir = Path.cwd().resolve() / "paper"
sys.path.insert(0, str(_paper_dir))

from config import COLORS, PLOT_SETTINGS, get_plot_path, get_data_path


def setup_style():
    plt.rcParams.update(
        {
            "font.family": PLOT_SETTINGS["font_family"],
            "font.sans-serif": PLOT_SETTINGS["font_sans_serif"],
            "font.size": PLOT_SETTINGS["font_size"],
            "axes.labelsize": PLOT_SETTINGS["axes_labelsize"],
            "axes.titlesize": PLOT_SETTINGS["axes_titlesize"],
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.dpi": PLOT_SETTINGS["figure_dpi"],
            "savefig.dpi": PLOT_SETTINGS["savefig_dpi"],
        }
    )


def _load_counts(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    needed = {"chrom", "pos", "ref_count", "alt_count"}
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"{path} missing columns: {sorted(missing)}")
    out = df[["chrom", "pos", "ref_count", "alt_count"]].copy()
    out["total"] = out["ref_count"] + out["alt_count"]
    out["ref_ratio"] = out["ref_count"] / out["total"].replace(0, np.nan)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", choices=["hg00731"], default="hg00731")
    ap.add_argument("--min-total", type=int, default=10)
    args = ap.parse_args()

    setup_style()

    data_dir = Path(get_data_path(2, f"{args.dataset}"))
    orig_path = data_dir / "wasp2_counts.original.tsv"
    wasp2_path = data_dir / "wasp2_counts.filtered.tsv"
    wasp1_path = data_dir / "wasp1_counts.filtered.tsv"

    if not orig_path.exists() or not wasp2_path.exists():
        raise SystemExit(f"Missing required inputs: {orig_path} and/or {wasp2_path}")

    orig = _load_counts(orig_path)
    wasp2 = _load_counts(wasp2_path)

    # Align on shared SNP IDs (same VCF query set); this keeps the comparison strict.
    key = ["chrom", "pos"]
    merged = orig.merge(
        wasp2[key + ["ref_count", "alt_count", "total", "ref_ratio"]],
        on=key,
        how="inner",
        suffixes=("_orig", "_wasp2"),
    )

    labels = ["WASP2"]
    filt_cols = [("ref_count_wasp2", "alt_count_wasp2", "total_wasp2", "ref_ratio_wasp2")]

    if wasp1_path.exists():
        wasp1 = _load_counts(wasp1_path)
        merged = merged.merge(
            wasp1[key + ["ref_count", "alt_count", "total", "ref_ratio"]],
            on=key,
            how="inner",
            suffixes=("", "_wasp1"),
        )
        labels.append("WASP1")
        filt_cols.append(("ref_count", "alt_count", "total", "ref_ratio"))  # wasp1 merge has no suffix

    # Apply coverage filter on both orig and filtered totals (per comparison).
    rows = []
    for label, (ref_c, alt_c, tot_c, rr_c) in zip(labels, filt_cols):
        df = merged.copy()
        df = df[(df["total_orig"] >= args.min_total) & (df[tot_c] >= args.min_total)]
        orig_tot = df["total_orig"].to_numpy()
        filt_tot = df[tot_c].to_numpy()
        dec = float((filt_tot < orig_tot).mean() * 100.0) if len(df) else 0.0
        same = float((filt_tot == orig_tot).mean() * 100.0) if len(df) else 0.0
        inc = float((filt_tot > orig_tot).mean() * 100.0) if len(df) else 0.0

        delta_ref_ratio = (df[rr_c] - df["ref_ratio_orig"]).to_numpy(dtype=float)
        delta_ref_ratio = delta_ref_ratio[np.isfinite(delta_ref_ratio)]
        mean_delta = float(delta_ref_ratio.mean()) if delta_ref_ratio.size else 0.0
        rows.append((label, len(df), dec, same, inc, mean_delta))

    summary = pd.DataFrame(
        rows,
        columns=[
            "filter",
            "n_sites",
            "pct_decrease",
            "pct_same",
            "pct_increase",
            "mean_delta_ref_ratio",
        ],
    )

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(5.2, 3.6), gridspec_kw={"hspace": 0.45})

    # Top: stacked %
    y = np.arange(len(summary))
    dec = summary["pct_decrease"].to_numpy()
    same = summary["pct_same"].to_numpy()
    inc = summary["pct_increase"].to_numpy()
    ax1.barh(y, dec, color=COLORS["blue"], alpha=0.55, edgecolor="black", linewidth=0.4, label="decrease")
    ax1.barh(y, same, left=dec, color=COLORS["black"], alpha=0.18, edgecolor="black", linewidth=0.4, label="same")
    ax1.barh(y, inc, left=dec + same, color=COLORS["vermillion"], alpha=0.55, edgecolor="black", linewidth=0.4, label="increase")
    ax1.set_yticks(y)
    ax1.set_yticklabels(summary["filter"], fontsize=8)
    ax1.set_xlim(0, 100)
    ax1.set_xlabel("% SNPs (Δ total counts)", fontsize=8)
    ax1.legend(fontsize=6.5, framealpha=0.9, loc="lower right")
    for yi, n in zip(y, summary["n_sites"].to_numpy()):
        ax1.text(0.01, yi, f"n={int(n):,}", transform=ax1.get_yaxis_transform(), ha="left", va="center", fontsize=7)

    # Bottom: mean Δ ref ratio
    ax2.bar(
        y,
        summary["mean_delta_ref_ratio"].to_numpy(),
        color=[COLORS["sky_blue"] if x == "WASP2" else COLORS["orange"] for x in summary["filter"]],
        edgecolor="black",
        linewidth=0.5,
        alpha=0.85,
    )
    ax2.axhline(0.0, color="black", linestyle="--", linewidth=1, alpha=0.6)
    ax2.set_xticks(y)
    ax2.set_xticklabels(summary["filter"], fontsize=8)
    ax2.set_ylabel("Mean Δ ref ratio\n(filt − orig)", fontsize=8)

    fig.suptitle("Supplement: WASP1 vs WASP2 filtering effect", fontsize=11, fontweight="bold", y=0.99)
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    out_dir = get_plot_path(2, "figure2").parent
    out_dir.mkdir(parents=True, exist_ok=True)
    png = out_dir / "figure2_c_wasp1_vs_wasp2.png"
    pdf = out_dir / "figure2_c_wasp1_vs_wasp2.pdf"
    fig.savefig(png, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    print(f"Saved: {png}")
    print(f"Saved: {pdf}")

    # Also save the summary table for debugging/captioning.
    out_tsv = data_dir / "wasp1_vs_wasp2_filtering_summary.tsv"
    summary.to_csv(out_tsv, sep="\t", index=False)
    print(f"Wrote: {out_tsv}")


if __name__ == "__main__":
    main()

