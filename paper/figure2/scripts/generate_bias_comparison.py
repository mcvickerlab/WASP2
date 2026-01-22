#!/usr/bin/env python3
"""
Generate bias comparison data for Figure 2 Panel C.

This compares ref/alt ratios between:
- an original alignment (pre-WASP2 filtering)
- a WASP2-processed alignment

To keep the comparison apples-to-apples and avoid very slow pileups over
millions of sites, this script consumes GATK ASEReadCounter output tables for
both alignments and computes bias metrics from those per-variant counts.
"""

from __future__ import annotations

import argparse
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd


def parse_gatk_ase_table(file_path: Path) -> pd.DataFrame:
    with file_path.open("r") as f:
        lines = [line for line in f if not line.startswith("#")]
    if not lines:
        return pd.DataFrame(columns=["chrom", "pos", "ref_count", "alt_count", "total_count"])

    df = pd.read_csv(StringIO("".join(lines)), sep="\t")
    df = df.rename(
        columns={
            "contig": "chrom",
            "position": "pos",
            "refCount": "ref_count",
            "altCount": "alt_count",
            "totalCount": "total_count",
        }
    )
    return df[["chrom", "pos", "ref_count", "alt_count", "total_count"]]


def merge_original_remapped(orig_df: pd.DataFrame, remap_df: pd.DataFrame) -> pd.DataFrame:
    merged = orig_df[["chrom", "pos"]].copy()
    merged["orig_ref"] = orig_df["ref_count"]
    merged["orig_alt"] = orig_df["alt_count"]
    merged["orig_total"] = orig_df["total_count"]

    remap_subset = remap_df[["chrom", "pos", "ref_count", "alt_count", "total_count"]].copy()
    remap_subset.columns = ["chrom", "pos", "remap_ref", "remap_alt", "remap_total"]
    merged = merged.merge(remap_subset, on=["chrom", "pos"], how="inner")
    return merged


def calculate_bias_metrics(df: pd.DataFrame):
    df = df.copy()
    df = df[(df["orig_total"] > 0) & (df["remap_total"] > 0)]
    df["orig_ref_ratio"] = df["orig_ref"] / df["orig_total"]
    df["remap_ref_ratio"] = df["remap_ref"] / df["remap_total"]

    df["orig_bias"] = np.abs(df["orig_ref_ratio"] - 0.5)
    df["remap_bias"] = np.abs(df["remap_ref_ratio"] - 0.5)

    stats = {
        "n_variants": int(len(df)),
        "orig_mean_bias": float(df["orig_bias"].mean()),
        "orig_median_bias": float(df["orig_bias"].median()),
        "remap_mean_bias": float(df["remap_bias"].mean()),
        "remap_median_bias": float(df["remap_bias"].median()),
    }
    stats["mean_bias_reduction_pct"] = (
        (stats["orig_mean_bias"] - stats["remap_mean_bias"]) / stats["orig_mean_bias"] * 100
        if stats["orig_mean_bias"] > 0
        else 0.0
    )
    stats["median_bias_reduction_pct"] = (
        (stats["orig_median_bias"] - stats["remap_median_bias"]) / stats["orig_median_bias"] * 100
        if stats["orig_median_bias"] > 0
        else 0.0
    )
    stats["n_improved"] = int((df["remap_bias"] < df["orig_bias"]).sum())
    stats["pct_improved"] = float(stats["n_improved"] / stats["n_variants"] * 100 if stats["n_variants"] else 0.0)
    return df, stats


def process_dataset(
    dataset_name: str, orig_gatk: Path, remap_gatk: Path, output_file: Path, min_total: int
):
    print(f"\n{'='*60}")
    print(f"Processing: {dataset_name}")
    print(f"{'='*60}")

    print("\nParsing GATK ASEReadCounter outputs...")
    print(f"  Original: {orig_gatk}")
    orig_df = parse_gatk_ase_table(orig_gatk)
    print(f"  Parsed {len(orig_df)} rows")

    print(f"  Remapped: {remap_gatk}")
    remap_df = parse_gatk_ase_table(remap_gatk)
    print(f"  Parsed {len(remap_df)} rows")

    print("\nMerging counts...")
    merged_df = merge_original_remapped(orig_df, remap_df)
    print(f"  Merged: {len(merged_df)} variants")

    if min_total > 0:
        before = len(merged_df)
        merged_df = merged_df[(merged_df["orig_total"] >= min_total) & (merged_df["remap_total"] >= min_total)]
        after = len(merged_df)
        print(f"\nCoverage filter: total_count >= {min_total} (kept {after:,} / {before:,})")

    print("\nCalculating bias metrics...")
    merged_df, stats = calculate_bias_metrics(merged_df)
    stats["min_total"] = int(min_total)

    print("\nBias Analysis Results:")
    print(f"  Variants analyzed: {stats['n_variants']}")
    print("\n  Original:")
    print(f"    Mean bias:   {stats['orig_mean_bias']:.4f}")
    print(f"    Median bias: {stats['orig_median_bias']:.4f}")
    print("\n  WASP2-processed:")
    print(f"    Mean bias:   {stats['remap_mean_bias']:.4f}")
    print(f"    Median bias: {stats['remap_median_bias']:.4f}")
    print("\n  Bias Reduction:")
    print(f"    Mean:   {stats['mean_bias_reduction_pct']:.1f}%")
    print(f"    Median: {stats['median_bias_reduction_pct']:.1f}%")
    print(f"\n  Variants with improved bias: {stats['n_improved']} ({stats['pct_improved']:.1f}%)")

    print(f"\nSaving to: {output_file}")
    output_file.parent.mkdir(parents=True, exist_ok=True)
    merged_df.to_csv(output_file, sep="\t", index=False)

    stats_file = output_file.parent / f"{output_file.stem}_stats.txt"
    with stats_file.open("w") as f:
        f.write(f"Bias Analysis: {dataset_name}\n")
        f.write("=" * 60 + "\n\n")
        for key, value in stats.items():
            f.write(f"{key}: {value}\n")
    print(f"Saved stats to: {stats_file}")


def main():
    ap = argparse.ArgumentParser(description="Generate bias comparison data for Figure 2 Panel C")
    ap.add_argument("--dataset", choices=["hg00731"], default="hg00731")
    ap.add_argument(
        "--min-total",
        type=int,
        default=10,
        help="Minimum total allele count required in both original and processed tables (default: 10)",
    )
    args = ap.parse_args()

    repo_root = Path(
        "/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp"
    )
    data_dir = repo_root / "paper" / "figure2" / "data"

    print("=" * 60)
    print("Bias Comparison Data Generation")
    print("=" * 60)

    if args.dataset == "hg00731":
        orig_gatk = data_dir / "hg00731" / "gatk_counts.original.table"
        remap_gatk = data_dir / "hg00731" / "gatk_counts.filtered.table"
        out_tsv = data_dir / "hg00731" / "original_vs_remapped.tsv"

        if not orig_gatk.exists():
            raise SystemExit(f"Missing original GATK output: {orig_gatk}")
        if not remap_gatk.exists():
            raise SystemExit(f"Missing remapped GATK output: {remap_gatk}")

        process_dataset("HG00731 RNA-seq", orig_gatk, remap_gatk, out_tsv, min_total=args.min_total)

    print("\n" + "=" * 60)
    print("Bias comparison data generation complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
