#!/usr/bin/env python3
"""
Generate Figure 3 Panel data from existing analyses.

Transforms raw experiment results into the format expected by generate_figure3.py:
  - Panel B: qq_plot_data.tsv (from stage1 experiments)
  - Panel C: volcano_data.tsv (from QTL-AI replication)
  - Panel E: qtl_stratification.tsv (from QTL-AI replication)
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from pathlib import Path

# Paths
REPO_ROOT = Path(__file__).resolve().parents[3]
CVPC_ROOT = REPO_ROOT.parent  # cvpc directory (parent of WASP2-exp)

# Source data locations
STAGE1_RESULTS = CVPC_ROOT / "experiments" / "stage1_global_dispersion" / "results"
QTL_AI_DIR = CVPC_ROOT / "results" / "analysis" / "peak_ai_qtl" / "from_genome_counts"

# Output directory
OUT_DIR = REPO_ROOT / "paper" / "figure3" / "data"


def generate_qq_data() -> None:
    """
    Generate Panel B QQ plot data from stage1 dispersion experiments.

    Aggregates p-values from multiple methods across 137 samples.
    """
    print("Generating QQ plot data (Panel B)...")

    # Find all experiment TSV files
    tsv_files = list(STAGE1_RESULTS.glob("exp*_*.tsv"))

    if not tsv_files:
        print(f"  WARNING: No experiment TSV files found in {STAGE1_RESULTS}")
        return

    print(f"  Found {len(tsv_files)} TSV files")

    # Aggregate p-values by method
    all_pvals = []

    for tsv in tsv_files[:500]:  # Limit to avoid memory issues
        try:
            df = pd.read_csv(tsv, sep='\t', usecols=['pvalue', 'experiment_id'], nrows=10000)
            if 'pvalue' in df.columns and 'experiment_id' in df.columns:
                method = df['experiment_id'].iloc[0] if len(df) > 0 else "unknown"
                for pval in df['pvalue'].dropna().sample(min(1000, len(df))):
                    all_pvals.append({'method': method, 'pvalue': pval})
        except Exception:
            continue

    if all_pvals:
        out_df = pd.DataFrame(all_pvals)
        out_path = OUT_DIR / "qq_plot_data.tsv"
        out_df.to_csv(out_path, sep='\t', index=False)
        print(f"  Saved: {out_path} ({len(out_df)} rows)")
    else:
        print("  WARNING: No p-values extracted")


def generate_volcano_data() -> None:
    """
    Generate Panel C volcano plot data from QTL-AI replication.

    Computes log2(alt/ref) for X-axis and -log10(pvalue) for Y-axis.
    """
    print("Generating volcano plot data (Panel C)...")

    atac_file = QTL_AI_DIR / "caqtl_ai_replication_atac.tsv"

    if not atac_file.exists():
        print(f"  WARNING: File not found: {atac_file}")
        return

    # Read ATAC QTL-AI data
    df = pd.read_csv(atac_file, sep='\t', usecols=[
        'ref_count', 'alt_count', 'total_reads', 'pvalue', 'qvalue', 'significant'
    ])

    # Filter for sites with reads
    df = df[(df['total_reads'] > 0) & (df['ref_count'] + df['alt_count'] > 0)].copy()

    # Compute log2 ratio (add pseudocount to avoid log(0))
    df['log2_ratio'] = np.log2((df['alt_count'] + 0.5) / (df['ref_count'] + 0.5))

    # Compute -log10(pvalue)
    df['neglog10p'] = -np.log10(df['pvalue'].clip(lower=1e-300))

    # Rename columns for figure script
    df = df.rename(columns={'qvalue': 'fdr'})

    # Select output columns
    out_df = df[['log2_ratio', 'pvalue', 'fdr']].dropna()

    # Sample if too large
    if len(out_df) > 50000:
        out_df = out_df.sample(50000, random_state=42)

    out_path = OUT_DIR / "volcano_data.tsv"
    out_df.to_csv(out_path, sep='\t', index=False)
    print(f"  Saved: {out_path} ({len(out_df)} rows)")


def generate_qtl_stratification() -> None:
    """
    Generate Panel E QTL stratification data.

    Stratifies allelic ratio by genotype (Het vs Hom lead SNP).
    """
    print("Generating QTL stratification data (Panel E)...")

    atac_file = QTL_AI_DIR / "caqtl_ai_replication_atac.tsv"

    if not atac_file.exists():
        print(f"  WARNING: File not found: {atac_file}")
        return

    # Read ATAC QTL-AI data
    df = pd.read_csv(atac_file, sep='\t', usecols=[
        'GT', 'ref_count', 'alt_count', 'total_reads', 'significant'
    ])

    # Filter for sites with reads
    df = df[(df['total_reads'] >= 10)].copy()

    # Parse genotype to het/hom
    def classify_gt(gt: str) -> str:
        if pd.isna(gt):
            return 'unknown'
        parts = gt.replace('|', '/').split('/')
        if len(parts) != 2:
            return 'unknown'
        return 'het' if parts[0] != parts[1] else 'hom'

    df['genotype'] = df['GT'].apply(classify_gt)

    # Filter to het/hom only
    df = df[df['genotype'].isin(['het', 'hom'])].copy()

    # Compute allelic ratio
    df['allelic_ratio'] = df['alt_count'] / (df['ref_count'] + df['alt_count'])

    # Select output columns
    out_df = df[['genotype', 'allelic_ratio']].dropna()

    # Sample if too large
    if len(out_df) > 50000:
        out_df = out_df.sample(50000, random_state=42)

    out_path = OUT_DIR / "qtl_stratification.tsv"
    out_df.to_csv(out_path, sep='\t', index=False)
    print(f"  Saved: {out_path} ({len(out_df)} rows)")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("Figure 3 Data Generation")
    print("=" * 60)
    print(f"Source: {CVPC_ROOT}")
    print(f"Output: {OUT_DIR}")
    print()

    generate_qq_data()
    print()
    generate_volcano_data()
    print()
    generate_qtl_stratification()

    print()
    print("=" * 60)
    print("Done! Now run: python paper/figure3/scripts/generate_figure3.py")
    print("=" * 60)


if __name__ == "__main__":
    main()
