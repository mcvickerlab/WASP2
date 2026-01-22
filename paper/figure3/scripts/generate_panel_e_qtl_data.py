#!/usr/bin/env python3
"""
Panel E: QTL Replication Data Generation

Generates data-driven QTL replication rates with:
- Wilson confidence intervals
- Breakdown by data type (RNA-seq eQTL, ATAC-seq caQTL)
- Proper sample size reporting

Data sources:
- caqtl_ai_replication_atac.tsv (ATAC-seq caQTL replication)
- eqtl_ai_replication_rna.tsv (RNA-seq eQTL replication)

Verified statistics (from 5-agent investigation):
- RNA-seq: 45.4% (109/240 variants with AI q≤0.10)
- ATAC-seq: 41.9% (6,698/15,975 variants with AI q≤0.10)

Author: WASP2 R&D Team
"""
import sys
from pathlib import Path
from typing import Dict, Tuple
import numpy as np
import pandas as pd

# Add paper directory to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from config import COLORS, PLOT_SETTINGS, REPO_ROOT

# Data paths
CVPC_ROOT = Path("/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc")
QTL_DATA_DIR = CVPC_ROOT / "results/analysis/peak_ai_qtl/from_genome_counts"

ATAC_REPLICATION = QTL_DATA_DIR / "caqtl_ai_replication_atac.tsv"
RNA_REPLICATION = QTL_DATA_DIR / "eqtl_ai_replication_rna.tsv"


def compute_wilson_ci(successes: int, n: int, z: float = 1.96) -> Tuple[float, float]:
    """Compute Wilson score confidence interval for binomial proportion."""
    if n == 0:
        return (0.0, 0.0)
    p = successes / n
    denominator = 1 + z**2 / n
    center = (p + z**2 / (2 * n)) / denominator
    margin = z * np.sqrt((p * (1 - p) + z**2 / (4 * n)) / n) / denominator
    return (max(0, center - margin), min(1, center + margin))


def load_replication_data(file_path: Path, data_type: str) -> pd.DataFrame:
    """
    Load QTL-AI replication data.

    Args:
        file_path: Path to replication TSV file
        data_type: 'atac' or 'rna' for logging

    Returns:
        DataFrame with replication data
    """
    if not file_path.exists():
        print(f"WARNING: {data_type} file not found: {file_path}")
        return pd.DataFrame()

    df = pd.read_csv(file_path, sep='\t')
    print(f"Loaded {data_type}: {len(df):,} rows")
    return df


def compute_replication_stats(
    df: pd.DataFrame,
    fdr_threshold: float = 0.10
) -> Dict:
    """
    Compute QTL replication statistics.

    The key insight: we count UNIQUE variants that show AI,
    not total rows (which can have multiple samples per variant).

    Args:
        df: DataFrame with qvalue and variant_key columns
        fdr_threshold: FDR threshold for significance

    Returns:
        Dictionary with replication statistics
    """
    if df.empty:
        return {}

    # Get unique variants
    if 'variant_key' in df.columns:
        variant_col = 'variant_key'
    else:
        # Construct variant key from chrom, pos, ref, alt
        df = df.copy()
        df['variant_key'] = df['chrom'] + ':' + df['pos0'].astype(str) + ':' + df['ref'] + ':' + df['alt']
        variant_col = 'variant_key'

    # Count unique testable variants (those with q-values)
    unique_variants = df[variant_col].unique()
    n_testable = len(unique_variants)

    # For each variant, check if ANY sample shows significant AI
    sig_variants = df[df['qvalue'] <= fdr_threshold][variant_col].unique()
    n_replicating = len(sig_variants)

    replication_rate = n_replicating / n_testable if n_testable > 0 else 0
    ci_low, ci_high = compute_wilson_ci(n_replicating, n_testable)

    return {
        'n_testable': n_testable,
        'n_replicating': n_replicating,
        'replication_rate': replication_rate,
        'ci_low': ci_low,
        'ci_high': ci_high,
        'fdr_threshold': fdr_threshold
    }


def generate_panel_e_data() -> pd.DataFrame:
    """
    Generate Panel E data for QTL replication rates.

    Returns:
        DataFrame with replication statistics for RNA and ATAC
    """
    print("=" * 60)
    print("PANEL E: QTL REPLICATION RATES")
    print("=" * 60)

    results = []

    # RNA-seq eQTL analysis
    print("\n" + "-" * 40)
    print("RNA-seq (eQTL) Replication")
    print("-" * 40)
    rna_df = load_replication_data(RNA_REPLICATION, 'RNA')
    if not rna_df.empty:
        rna_stats = compute_replication_stats(rna_df)
        rna_stats['data_type'] = 'RNA-seq (eQTL)'
        rna_stats['samples'] = 137
        results.append(rna_stats)

        print(f"Testable variants: {rna_stats['n_testable']:,}")
        print(f"Replicating variants (AI q≤0.10): {rna_stats['n_replicating']:,}")
        print(f"Replication rate: {rna_stats['replication_rate']*100:.1f}%")
        print(f"95% CI: [{rna_stats['ci_low']*100:.1f}%, {rna_stats['ci_high']*100:.1f}%]")

    # ATAC-seq caQTL analysis
    print("\n" + "-" * 40)
    print("ATAC-seq (caQTL) Replication")
    print("-" * 40)
    atac_df = load_replication_data(ATAC_REPLICATION, 'ATAC')
    if not atac_df.empty:
        atac_stats = compute_replication_stats(atac_df)
        atac_stats['data_type'] = 'ATAC-seq (caQTL)'
        atac_stats['samples'] = 138
        results.append(atac_stats)

        print(f"Testable variants: {atac_stats['n_testable']:,}")
        print(f"Replicating variants (AI q≤0.10): {atac_stats['n_replicating']:,}")
        print(f"Replication rate: {atac_stats['replication_rate']*100:.1f}%")
        print(f"95% CI: [{atac_stats['ci_low']*100:.1f}%, {atac_stats['ci_high']*100:.1f}%]")

    print()
    print("=" * 60)
    print("INTERPRETATION")
    print("=" * 60)
    print("""
High replication rates (>40%) indicate that:
1. QTL variants discovered in population-level analysis
   show significant allelic imbalance at the individual level
2. WASP2's AI detection captures biologically meaningful signal
3. The statistical framework is appropriately calibrated
""")

    return pd.DataFrame(results)


def generate_panel_e_plot(df: pd.DataFrame, output_dir: Path = None):
    """Generate publication-quality bar plot for Panel E."""
    import matplotlib.pyplot as plt

    if output_dir is None:
        output_dir = Path(__file__).parent.parent / "plots"
    output_dir.mkdir(parents=True, exist_ok=True)

    if df.empty:
        print("No data for plot")
        return None

    plt.rcParams.update({
        'font.family': PLOT_SETTINGS['font_family'],
        'font.sans-serif': PLOT_SETTINGS['font_sans_serif'],
        'font.size': PLOT_SETTINGS['font_size'],
        'axes.labelsize': PLOT_SETTINGS['axes_labelsize'],
        'axes.titlesize': PLOT_SETTINGS['axes_titlesize'],
        'axes.spines.top': False,
        'axes.spines.right': False,
    })

    fig, ax = plt.subplots(figsize=(5, 4))

    x = range(len(df))
    colors = [COLORS['blue'], COLORS['orange']][:len(df)]

    bars = ax.bar(x, df['replication_rate'] * 100, color=colors, edgecolor='black', linewidth=0.5)

    # Error bars
    yerr_low = (df['replication_rate'] - df['ci_low']) * 100
    yerr_high = (df['ci_high'] - df['replication_rate']) * 100
    ax.errorbar(x, df['replication_rate'] * 100, yerr=[yerr_low, yerr_high],
                fmt='none', color='black', capsize=5, linewidth=1.5)

    ax.set_xticks(x)
    labels = [f"{row['data_type']}\n{row['samples']} samples\n({row['n_replicating']:,}/{row['n_testable']:,})"
              for _, row in df.iterrows()]
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel('QTL Replication Rate (%)')
    ax.set_ylim(0, 60)

    # Add percentage labels on bars
    for i, (_, row) in enumerate(df.iterrows()):
        ax.annotate(f"{row['replication_rate']*100:.1f}%",
                   xy=(i, row['replication_rate'] * 100 + 3),
                   ha='center', va='bottom', fontsize=10, fontweight='bold')

    ax.set_title('Panel E: QTL Replication via Allelic Imbalance', fontweight='bold', fontsize=10)

    # Add horizontal line at 50% for reference
    ax.axhline(50, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax.text(len(df) - 0.5, 51, '50%', fontsize=7, color='gray', ha='right')

    plt.tight_layout()

    output_path = output_dir / "panel_e_qtl_replication.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(output_path.with_suffix('.pdf'), bbox_inches='tight', facecolor='white')
    print(f"\nSaved: {output_path}")
    print(f"Saved: {output_path.with_suffix('.pdf')}")

    return fig


def save_panel_e_data(df: pd.DataFrame, output_dir: Path = None):
    """Save Panel E data to TSV file."""
    if output_dir is None:
        output_dir = Path(__file__).parent.parent / "data" / "qtl_replication"
    output_dir.mkdir(parents=True, exist_ok=True)

    if not df.empty:
        output_path = output_dir / "panel_e_qtl_replication.tsv"
        df.to_csv(output_path, sep='\t', index=False)
        print(f"Saved: {output_path}")


if __name__ == '__main__':
    df = generate_panel_e_data()

    if not df.empty:
        save_panel_e_data(df)
        generate_panel_e_plot(df)
    else:
        print("ERROR: No data loaded. Check file paths.")
