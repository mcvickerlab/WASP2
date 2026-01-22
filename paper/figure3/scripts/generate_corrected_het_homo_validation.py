#!/usr/bin/env python3
"""
Corrected Het/Homo Allelic Imbalance Validation

ROOT CAUSE FIX: The original lead-het/homo analysis had a fundamental design flaw:
- Lead variant position != Test SNV position
- 0 of 4,753 lead variants found in Stage1 data
- No LD linkage considered between positions

This script implements Option B: Stratify SNVs by their OWN genotype (GT column)
from the QTL replication files, not the lead variant genotype.

Expected result (if WASP2 works correctly):
- Heterozygous sites should show higher AI rates (both alleles expressed → can differ)
- Homozygous sites should show lower AI rates (same allele → no allelic difference)

Author: WASP2 R&D Team
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from typing import Tuple, Dict

# Add paper directory to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from config import COLORS, PLOT_SETTINGS, REPO_ROOT

# Data paths
CVPC_ROOT = Path("/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc")
QTL_DATA_DIR = CVPC_ROOT / "results/analysis/peak_ai_qtl/from_genome_counts"

ATAC_REPLICATION = QTL_DATA_DIR / "caqtl_ai_replication_atac.tsv"
RNA_REPLICATION = QTL_DATA_DIR / "eqtl_ai_replication_rna.tsv"


def classify_genotype(gt_str: str) -> str:
    """
    Classify VCF genotype string as het, homo-ref, homo-alt, or unknown.

    Args:
        gt_str: Genotype string like "C|G", "A|A", "G|G", "0|1", etc.

    Returns:
        Classification: 'het', 'homo', or 'unknown'
    """
    if pd.isna(gt_str) or gt_str == '' or gt_str == '.':
        return 'unknown'

    # Handle phased (|) and unphased (/) separators
    if '|' in gt_str:
        parts = gt_str.split('|')
    elif '/' in gt_str:
        parts = gt_str.split('/')
    else:
        return 'unknown'

    if len(parts) != 2:
        return 'unknown'

    allele1, allele2 = parts[0].strip(), parts[1].strip()

    # Handle both numeric (0/1) and nucleotide (A/G) formats
    if allele1 == allele2:
        return 'homo'
    else:
        return 'het'


def compute_wilson_ci(successes: int, n: int, z: float = 1.96) -> Tuple[float, float]:
    """
    Compute Wilson score confidence interval for binomial proportion.

    More accurate than normal approximation, especially for extreme proportions.

    Args:
        successes: Number of successes (e.g., significant AI calls)
        n: Total trials
        z: Z-score for confidence level (1.96 = 95% CI)

    Returns:
        (lower, upper) bounds of confidence interval
    """
    if n == 0:
        return (0.0, 0.0)

    p = successes / n
    denominator = 1 + z**2 / n
    center = (p + z**2 / (2 * n)) / denominator
    margin = z * np.sqrt((p * (1 - p) + z**2 / (4 * n)) / n) / denominator

    return (max(0, center - margin), min(1, center + margin))


def load_and_classify_data(file_path: Path, data_type: str) -> pd.DataFrame:
    """
    Load QTL replication data and classify by SNV's own genotype.

    Args:
        file_path: Path to replication TSV file
        data_type: 'atac' or 'rna' for logging

    Returns:
        DataFrame with added 'gt_class' column
    """
    if not file_path.exists():
        print(f"WARNING: {data_type} file not found: {file_path}")
        return pd.DataFrame()

    df = pd.read_csv(file_path, sep='\t')
    print(f"Loaded {data_type}: {len(df):,} rows")

    # Classify genotypes using the GT column (SNV's OWN genotype)
    df['gt_class'] = df['GT'].apply(classify_genotype)

    # Print classification summary
    class_counts = df['gt_class'].value_counts()
    print(f"  Genotype classification:")
    for gt_class, count in class_counts.items():
        print(f"    {gt_class}: {count:,} ({100*count/len(df):.1f}%)")

    return df


def compute_validation_statistics(
    df: pd.DataFrame,
    fdr_threshold: float = 0.10,
    min_coverage: int = 10
) -> pd.DataFrame:
    """
    Compute corrected het/homo validation statistics.

    This is the KEY FIX: stratify by SNV's own genotype, not lead genotype.

    Args:
        df: DataFrame with gt_class and qvalue columns
        fdr_threshold: FDR threshold for significance
        min_coverage: Minimum read coverage to include

    Returns:
        Summary DataFrame with AI rates by genotype class
    """
    if df.empty:
        return pd.DataFrame()

    # Filter by coverage
    if 'total_reads' in df.columns:
        df = df[df['total_reads'] >= min_coverage].copy()
        print(f"After coverage filter (>={min_coverage}): {len(df):,} rows")

    results = []

    for gt_class in ['het', 'homo', 'unknown']:
        subset = df[df['gt_class'] == gt_class]
        n_total = len(subset)

        if n_total == 0:
            continue

        # Count significant AI (using qvalue from the data)
        n_sig = (subset['qvalue'] <= fdr_threshold).sum()
        ai_rate = n_sig / n_total

        # Wilson confidence interval
        ci_low, ci_high = compute_wilson_ci(n_sig, n_total)

        results.append({
            'genotype_class': gt_class,
            'n_total': n_total,
            'n_significant': n_sig,
            'ai_rate': ai_rate,
            'ci_low': ci_low,
            'ci_high': ci_high
        })

    return pd.DataFrame(results)


def run_validation(verbose: bool = True) -> Dict[str, pd.DataFrame]:
    """
    Run the corrected het/homo validation on both ATAC and RNA data.

    Returns:
        Dictionary with 'atac' and 'rna' DataFrames
    """
    results = {}

    print("=" * 60)
    print("CORRECTED HET/HOMO ALLELIC IMBALANCE VALIDATION")
    print("FIX: Using SNV's own genotype (GT column), not lead genotype")
    print("=" * 60)
    print()

    # Process ATAC data
    print("-" * 40)
    print("ATAC-seq (caQTL) Analysis")
    print("-" * 40)
    atac_df = load_and_classify_data(ATAC_REPLICATION, 'ATAC')
    if not atac_df.empty:
        atac_stats = compute_validation_statistics(atac_df)
        results['atac'] = atac_stats

        if verbose and not atac_stats.empty:
            print("\nATAC Results:")
            for _, row in atac_stats.iterrows():
                print(f"  {row['genotype_class']}: {row['ai_rate']*100:.2f}% "
                      f"({row['n_significant']:,}/{row['n_total']:,}) "
                      f"[95% CI: {row['ci_low']*100:.2f}-{row['ci_high']*100:.2f}%]")

    print()

    # Process RNA data
    print("-" * 40)
    print("RNA-seq (eQTL) Analysis")
    print("-" * 40)
    rna_df = load_and_classify_data(RNA_REPLICATION, 'RNA')
    if not rna_df.empty:
        # RNA file uses 'pvalue_ai' column, rename for consistency
        if 'pvalue_ai' in rna_df.columns and 'pvalue' not in rna_df.columns:
            rna_df = rna_df.rename(columns={'pvalue_ai': 'pvalue'})

        rna_stats = compute_validation_statistics(rna_df)
        results['rna'] = rna_stats

        if verbose and not rna_stats.empty:
            print("\nRNA Results:")
            for _, row in rna_stats.iterrows():
                print(f"  {row['genotype_class']}: {row['ai_rate']*100:.2f}% "
                      f"({row['n_significant']:,}/{row['n_total']:,}) "
                      f"[95% CI: {row['ci_low']*100:.2f}-{row['ci_high']*100:.2f}%]")

    print()
    print("=" * 60)
    print("INTERPRETATION")
    print("=" * 60)
    print("""
Expected pattern (if biological signal is real):
  - Het SNVs should show HIGHER AI rates
    (heterozygous sites can show allelic imbalance)
  - Homo SNVs should show LOWER AI rates
    (homozygous sites have same allele on both chromosomes)

If this pattern is observed, it validates that WASP2's AI detection
is capturing real biological signal, not statistical artifacts.
""")

    return results


def generate_validation_plot(results: Dict[str, pd.DataFrame], output_dir: Path = None):
    """
    Generate publication-quality bar plot comparing het vs homo AI rates.

    Args:
        results: Dictionary from run_validation()
        output_dir: Directory for output files (default: paper/figure3/plots)
    """
    import matplotlib.pyplot as plt

    if output_dir is None:
        output_dir = Path(__file__).parent.parent / "plots"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Set up style
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

    for ax, (data_type, title) in zip(axes, [('atac', 'ATAC-seq (caQTL)'), ('rna', 'RNA-seq (eQTL)')]):
        if data_type not in results or results[data_type].empty:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(title)
            continue

        df = results[data_type]
        # Only show het and homo (exclude unknown)
        df = df[df['genotype_class'].isin(['het', 'homo'])].copy()

        if df.empty:
            ax.text(0.5, 0.5, 'No het/homo data', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(title)
            continue

        x = range(len(df))
        colors = [COLORS['blue'] if gc == 'het' else COLORS['orange'] for gc in df['genotype_class']]

        bars = ax.bar(x, df['ai_rate'] * 100, color=colors, edgecolor='black', linewidth=0.5)

        # Add error bars (95% CI)
        yerr_low = (df['ai_rate'] - df['ci_low']) * 100
        yerr_high = (df['ci_high'] - df['ai_rate']) * 100
        ax.errorbar(x, df['ai_rate'] * 100, yerr=[yerr_low, yerr_high],
                    fmt='none', color='black', capsize=3)

        ax.set_xticks(x)
        ax.set_xticklabels([f"{row['genotype_class'].title()}\n(n={row['n_total']:,})"
                           for _, row in df.iterrows()], fontsize=8)
        ax.set_ylabel('AI Rate (%)')
        ax.set_title(title, fontweight='bold')

        # Add significance annotation if both categories exist
        if len(df) == 2:
            het_rate = df[df['genotype_class'] == 'het']['ai_rate'].values[0]
            homo_rate = df[df['genotype_class'] == 'homo']['ai_rate'].values[0]
            fold_change = het_rate / homo_rate if homo_rate > 0 else float('inf')
            ax.annotate(f'{fold_change:.1f}x', xy=(0.5, max(df['ai_rate'] * 100) + 2),
                       ha='center', fontsize=8, fontweight='bold')

    fig.suptitle('Corrected Het/Homo Validation\n(Using SNV own genotype, not lead genotype)',
                 fontsize=10, fontweight='bold')
    plt.tight_layout()

    # Save
    output_path = output_dir / "het_homo_validation_corrected.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(output_path.with_suffix('.pdf'), bbox_inches='tight', facecolor='white')
    print(f"\nSaved: {output_path}")
    print(f"Saved: {output_path.with_suffix('.pdf')}")

    return fig


def save_validation_data(results: Dict[str, pd.DataFrame], output_dir: Path = None):
    """Save validation results to TSV files."""
    if output_dir is None:
        output_dir = Path(__file__).parent.parent / "data" / "het_homo_validation"
    output_dir.mkdir(parents=True, exist_ok=True)

    for data_type, df in results.items():
        if not df.empty:
            output_path = output_dir / f"het_homo_validation_{data_type}.tsv"
            df.to_csv(output_path, sep='\t', index=False)
            print(f"Saved: {output_path}")


if __name__ == '__main__':
    results = run_validation()

    if any(not df.empty for df in results.values()):
        save_validation_data(results)
        generate_validation_plot(results)
    else:
        print("ERROR: No data loaded. Check file paths.")
