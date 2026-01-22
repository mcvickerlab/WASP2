#!/usr/bin/env python3
"""
Panel F: Genomic Feature AI Enrichment Analysis

Compares allelic imbalance rates across genomic features:
- Promoter regions (is_promoter=True)
- Non-promoter regions (is_promoter=False)
- ATAC peaks (is_capeak=True for ATAC data)
- Background regions

Uses Fisher's exact test to compute enrichment vs background.

Data sources:
- RNA: stage1_with_promoters/ (promoter annotations)
- ATAC: stage1_with_peaks/ (chromatin accessibility peaks)

Author: WASP2 R&D Team
"""
import sys
from pathlib import Path
from typing import Dict, List, Tuple
import numpy as np
import pandas as pd
from scipy import stats
import warnings

warnings.filterwarnings('ignore', category=pd.errors.SettingWithCopyWarning)

# Add paper directory to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from config import COLORS, PLOT_SETTINGS, REPO_ROOT

# Data paths
CVPC_ROOT = Path("/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc")
RNA_PROMOTERS_DIR = CVPC_ROOT / "results/analysis/peak_ai_rna/stage1_with_promoters"
ATAC_PEAKS_DIR = CVPC_ROOT / "results/analysis/peak_ai/stage1_with_peaks"


def compute_wilson_ci(successes: int, n: int, z: float = 1.96) -> Tuple[float, float]:
    """Compute Wilson score confidence interval."""
    if n == 0:
        return (0.0, 0.0)
    p = successes / n
    denominator = 1 + z**2 / n
    center = (p + z**2 / (2 * n)) / denominator
    margin = z * np.sqrt((p * (1 - p) + z**2 / (4 * n)) / n) / denominator
    return (max(0, center - margin), min(1, center + margin))


def fisher_enrichment(a: int, b: int, c: int, d: int) -> Tuple[float, float]:
    """
    Compute Fisher's exact test for enrichment.

    Table:
              | Sig  | Not-Sig |
    Feature   |  a   |    b    |
    Background|  c   |    d    |

    Returns:
        (odds_ratio, p_value)
    """
    odds_ratio, pvalue = stats.fisher_exact([[a, b], [c, d]], alternative='greater')
    return odds_ratio, pvalue


def load_rna_promoter_data(experiment: str = 'exp09', max_files: int = None) -> pd.DataFrame:
    """
    Load RNA-seq data with promoter annotations.

    Args:
        experiment: Experiment ID (e.g., 'exp09')
        max_files: Maximum number of sample files to load

    Returns:
        Combined DataFrame with promoter annotations
    """
    pattern = f"{experiment}_*_with_promoters.tsv"
    files = sorted(RNA_PROMOTERS_DIR.glob(pattern))

    if not files:
        print(f"WARNING: No files found for {experiment} in {RNA_PROMOTERS_DIR}")
        return pd.DataFrame()

    if max_files:
        files = files[:max_files]

    dfs = []
    for f in files:
        try:
            df = pd.read_csv(f, sep='\t', usecols=[
                'chrom', 'pos0', 'qvalue', 'significant', 'is_promoter', 'total_reads'
            ])
            dfs.append(df)
        except Exception as e:
            print(f"Warning: Error reading {f}: {e}")

    if not dfs:
        return pd.DataFrame()

    combined = pd.concat(dfs, ignore_index=True)
    print(f"Loaded RNA promoter data: {len(combined):,} rows")
    return combined


def load_atac_peak_data(experiment: str = 'exp04', max_files: int = None) -> pd.DataFrame:
    """
    Load ATAC-seq data with peak annotations.

    Args:
        experiment: Experiment ID (e.g., 'exp04')
        max_files: Maximum number of sample files to load

    Returns:
        Combined DataFrame with peak annotations
    """
    pattern = f"{experiment}_*_with_peaks.tsv"
    files = sorted(ATAC_PEAKS_DIR.glob(pattern))

    if not files:
        print(f"WARNING: No files found for {experiment} in {ATAC_PEAKS_DIR}")
        return pd.DataFrame()

    if max_files:
        files = files[:max_files]

    dfs = []
    for f in files:
        try:
            df = pd.read_csv(f, sep='\t', usecols=[
                'chrom', 'pos0', 'qvalue', 'significant', 'is_capeak', 'total_reads'
            ])
            dfs.append(df)
        except Exception as e:
            print(f"Warning: Error reading {f}: {e}")

    if not dfs:
        return pd.DataFrame()

    combined = pd.concat(dfs, ignore_index=True)
    print(f"Loaded ATAC peak data: {len(combined):,} rows")
    return combined


def compute_feature_stats(
    df: pd.DataFrame,
    feature_col: str,
    fdr_threshold: float = 0.10,
    min_coverage: int = 10
) -> pd.DataFrame:
    """
    Compute AI rates by genomic feature.

    Args:
        df: DataFrame with qvalue and feature column
        feature_col: Column name for feature classification (e.g., 'is_promoter')
        fdr_threshold: FDR threshold for significance
        min_coverage: Minimum read coverage to include

    Returns:
        Summary DataFrame with AI rates per feature
    """
    if df.empty or feature_col not in df.columns:
        return pd.DataFrame()

    # Filter by coverage
    df = df[df['total_reads'] >= min_coverage].copy()
    print(f"After coverage filter (>={min_coverage}): {len(df):,} rows")

    results = []

    for feature_val in [True, False]:
        subset = df[df[feature_col] == feature_val]
        n_total = len(subset)

        if n_total == 0:
            continue

        n_sig = (subset['qvalue'] <= fdr_threshold).sum()
        ai_rate = n_sig / n_total
        ci_low, ci_high = compute_wilson_ci(n_sig, n_total)

        feature_name = 'In feature' if feature_val else 'Background'

        results.append({
            'feature': feature_name,
            'feature_value': feature_val,
            'n_total': n_total,
            'n_significant': n_sig,
            'ai_rate': ai_rate,
            'ci_low': ci_low,
            'ci_high': ci_high
        })

    if len(results) == 2:
        # Compute enrichment (feature vs background)
        feat_row = [r for r in results if r['feature_value']][0]
        bg_row = [r for r in results if not r['feature_value']][0]

        a = feat_row['n_significant']
        b = feat_row['n_total'] - a
        c = bg_row['n_significant']
        d = bg_row['n_total'] - c

        odds_ratio, pvalue = fisher_enrichment(a, b, c, d)

        for r in results:
            r['enrichment_or'] = odds_ratio if r['feature_value'] else 1.0
            r['enrichment_pval'] = pvalue if r['feature_value'] else 1.0

    return pd.DataFrame(results)


def generate_panel_f_data(max_files: int = None) -> Dict[str, pd.DataFrame]:
    """
    Generate Panel F data for genomic feature enrichment.

    Returns:
        Dictionary with 'rna_promoter' and 'atac_peak' DataFrames
    """
    print("=" * 60)
    print("PANEL F: GENOMIC FEATURE AI ENRICHMENT")
    print("=" * 60)

    results = {}

    # RNA promoter analysis
    print("\n" + "-" * 40)
    print("RNA-seq: Promoter vs Non-promoter")
    print("-" * 40)
    rna_df = load_rna_promoter_data('exp09', max_files)
    if not rna_df.empty:
        rna_stats = compute_feature_stats(rna_df, 'is_promoter')
        results['rna_promoter'] = rna_stats

        if not rna_stats.empty:
            print("\nRNA Promoter Results:")
            for _, row in rna_stats.iterrows():
                or_str = f" (OR={row['enrichment_or']:.2f}, p={row['enrichment_pval']:.2e})" if row['feature_value'] else ""
                print(f"  {row['feature']}: {row['ai_rate']*100:.2f}% "
                      f"({row['n_significant']:,}/{row['n_total']:,}){or_str}")

    # ATAC peak analysis
    print("\n" + "-" * 40)
    print("ATAC-seq: Peak vs Non-peak")
    print("-" * 40)
    atac_df = load_atac_peak_data('exp04', max_files)
    if not atac_df.empty:
        atac_stats = compute_feature_stats(atac_df, 'is_capeak')
        results['atac_peak'] = atac_stats

        if not atac_stats.empty:
            print("\nATAC Peak Results:")
            for _, row in atac_stats.iterrows():
                or_str = f" (OR={row['enrichment_or']:.2f}, p={row['enrichment_pval']:.2e})" if row['feature_value'] else ""
                print(f"  {row['feature']}: {row['ai_rate']*100:.2f}% "
                      f"({row['n_significant']:,}/{row['n_total']:,}){or_str}")

    return results


def generate_panel_f_plot(results: Dict[str, pd.DataFrame], output_dir: Path = None):
    """Generate publication-quality plot for Panel F."""
    import matplotlib.pyplot as plt

    if output_dir is None:
        output_dir = Path(__file__).parent.parent / "plots"
    output_dir.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update({
        'font.family': PLOT_SETTINGS['font_family'],
        'font.sans-serif': PLOT_SETTINGS['font_sans_serif'],
        'font.size': PLOT_SETTINGS['font_size'],
        'axes.labelsize': PLOT_SETTINGS['axes_labelsize'],
        'axes.titlesize': PLOT_SETTINGS['axes_titlesize'],
        'axes.spines.top': False,
        'axes.spines.right': False,
    })

    fig, axes = plt.subplots(1, 2, figsize=(8, 4))

    titles = [
        ('rna_promoter', 'RNA-seq: Promoter Enrichment'),
        ('atac_peak', 'ATAC-seq: Peak Enrichment')
    ]

    for ax, (data_key, title) in zip(axes, titles):
        if data_key not in results or results[data_key].empty:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(title)
            continue

        df = results[data_key]

        x = range(len(df))
        colors = [COLORS['blue'] if fv else COLORS['orange'] for fv in df['feature_value']]
        labels = ['In feature' if fv else 'Background' for fv in df['feature_value']]

        bars = ax.bar(x, df['ai_rate'] * 100, color=colors, edgecolor='black', linewidth=0.5)

        # Error bars
        yerr_low = (df['ai_rate'] - df['ci_low']) * 100
        yerr_high = (df['ci_high'] - df['ai_rate']) * 100
        ax.errorbar(x, df['ai_rate'] * 100, yerr=[yerr_low, yerr_high],
                    fmt='none', color='black', capsize=3)

        ax.set_xticks(x)
        ax.set_xticklabels([f"{lab}\n(n={n:,})" for lab, n in zip(labels, df['n_total'])], fontsize=7)
        ax.set_ylabel('AI Rate (%)')
        ax.set_title(title, fontweight='bold', fontsize=9)

        # Add enrichment annotation
        feat_df = df[df['feature_value']]
        if not feat_df.empty:
            or_val = feat_df['enrichment_or'].values[0]
            pval = feat_df['enrichment_pval'].values[0]
            sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "ns"
            ax.annotate(f'OR={or_val:.1f} {sig}',
                       xy=(0.5, max(df['ai_rate'] * 100) + 1),
                       ha='center', fontsize=8, fontweight='bold')

    fig.suptitle('Panel F: Genomic Feature AI Enrichment', fontsize=10, fontweight='bold')
    plt.tight_layout()

    output_path = output_dir / "panel_f_genomic_features.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(output_path.with_suffix('.pdf'), bbox_inches='tight', facecolor='white')
    print(f"\nSaved: {output_path}")
    print(f"Saved: {output_path.with_suffix('.pdf')}")

    return fig


def save_panel_f_data(results: Dict[str, pd.DataFrame], output_dir: Path = None):
    """Save Panel F data to TSV files."""
    if output_dir is None:
        output_dir = Path(__file__).parent.parent / "data" / "genomic_features"
    output_dir.mkdir(parents=True, exist_ok=True)

    for name, df in results.items():
        if not df.empty:
            output_path = output_dir / f"panel_f_{name}.tsv"
            df.to_csv(output_path, sep='\t', index=False)
            print(f"Saved: {output_path}")


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Generate Panel F genomic feature enrichment data')
    parser.add_argument('--max-files', type=int, default=None,
                        help='Max sample files to load (for testing)')

    args = parser.parse_args()

    results = generate_panel_f_data(args.max_files)

    if any(not df.empty for df in results.values()):
        save_panel_f_data(results)
        generate_panel_f_plot(results)
    else:
        print("ERROR: No data loaded. Check file paths.")
