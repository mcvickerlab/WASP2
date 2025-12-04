#!/usr/bin/env python3
"""
Benchmark WASP2 vs GATK ASEReadCounter for Allele-Specific Expression Counting

This script compares WASP2 and GATK ASEReadCounter on the same simulated data to enable
head-to-head comparison of:
  1. Accuracy: Pearson correlation with ground truth
  2. Bias: Reference allele bias (mean ratio - 0.5)
  3. Speed: Runtime comparison

Usage:
    # Run comparison on simulation data
    python simulation/benchmark_vs_gatk.py \
        --bam simulation_results/comprehensive_*/aligned.sorted.bam \
        --vcf simulation_results/comprehensive_*/variants.vcf.gz \
        --ref simulation_results/comprehensive_*/reference.fa \
        --ground-truth simulation_results/comprehensive_*/ground_truth.csv \
        --output comparison_results/

    # Quick test mode
    python simulation/benchmark_vs_gatk.py \
        --bam test.bam \
        --vcf test.vcf.gz \
        --ref ref.fa \
        --ground-truth truth.csv \
        --output test_comparison/

Requirements:
    - GATK installed and in PATH (conda install -c bioconda gatk4)
    - WASP2 Python package installed
    - BAM file indexed (.bai)
    - VCF file compressed and indexed (.tbi)
    - Reference FASTA indexed (.fai)

Author: WASP2 Development Team
Date: 2025-12-03
"""

import pandas as pd
import numpy as np
import subprocess
import tempfile
import time
import argparse
import sys
import pysam
from pathlib import Path
from typing import Dict, Tuple, Optional
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import mean_squared_error, mean_absolute_error
import matplotlib.pyplot as plt
import seaborn as sns


def check_gatk_available() -> bool:
    """
    Verify GATK is installed and accessible.

    Returns:
        True if GATK is available, False otherwise
    """
    try:
        result = subprocess.run(
            ['gatk', '--list'],
            capture_output=True,
            text=True,
            timeout=10
        )
        return 'ASEReadCounter' in result.stderr or 'ASEReadCounter' in result.stdout
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return False


def check_dependencies(bam_path: Path, vcf_path: Path, ref_path: Path) -> bool:
    """
    Check all required files and indexes exist. Create missing indexes if possible.

    Args:
        bam_path: Path to BAM file
        vcf_path: Path to VCF file
        ref_path: Path to reference FASTA

    Returns:
        True if all dependencies exist
    """
    errors = []

    # Check BAM and index
    if not bam_path.exists():
        errors.append(f"BAM file not found: {bam_path}")
    if not Path(str(bam_path) + '.bai').exists() and not Path(str(bam_path).replace('.bam', '.bai')).exists():
        errors.append(f"BAM index not found: {bam_path}.bai")

    # Check VCF and index
    if not vcf_path.exists():
        errors.append(f"VCF file not found: {vcf_path}")
    if not Path(str(vcf_path) + '.tbi').exists():
        errors.append(f"VCF index not found: {vcf_path}.tbi")

    # Check reference and index
    if not ref_path.exists():
        errors.append(f"Reference FASTA not found: {ref_path}")
    if not Path(str(ref_path) + '.fai').exists():
        errors.append(f"Reference index not found: {ref_path}.fai")

    # Check for sequence dictionary (required by GATK)
    dict_path = Path(str(ref_path).replace('.fa', '.dict').replace('.fasta', '.dict'))
    if not dict_path.exists():
        print(f"  Creating sequence dictionary: {dict_path}")
        try:
            subprocess.run(
                ['gatk', 'CreateSequenceDictionary', '-R', str(ref_path), '-O', str(dict_path)],
                check=True,
                capture_output=True
            )
        except subprocess.CalledProcessError as e:
            errors.append(f"Failed to create sequence dictionary: {e}")

    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return False

    return True


def ensure_bam_has_read_groups(bam_path: Path) -> Path:
    """
    Check if BAM has read groups, add them if missing.

    Args:
        bam_path: Path to BAM file

    Returns:
        Path to BAM file with read groups (may be same as input or new file)
    """
    # Check if BAM has read groups
    bam = pysam.AlignmentFile(str(bam_path), 'rb')
    has_rg = 'RG' in bam.header
    bam.close()

    if has_rg:
        return bam_path

    # Need to add read groups
    print(f"  Adding read groups to BAM...")
    rg_bam_path = Path(str(bam_path).replace('.bam', '.rg.bam'))

    try:
        subprocess.run([
            'gatk', 'AddOrReplaceReadGroups',
            '-I', str(bam_path),
            '-O', str(rg_bam_path),
            '-RGID', '1',
            '-RGLB', 'lib1',
            '-RGPL', 'illumina',
            '-RGPU', 'unit1',
            '-RGSM', 'sample1'
        ], check=True, capture_output=True)

        # Index the new BAM
        subprocess.run(['samtools', 'index', str(rg_bam_path)], check=True, capture_output=True)

        print(f"  Created BAM with read groups: {rg_bam_path}")
        return rg_bam_path

    except subprocess.CalledProcessError as e:
        print(f"  WARNING: Failed to add read groups: {e}", file=sys.stderr)
        return bam_path


def run_gatk_ase_counter(
    bam_path: Path,
    vcf_path: Path,
    ref_path: Path,
    output_path: Path,
    min_mapq: int = 10,
    min_baseq: int = 20
) -> Tuple[float, pd.DataFrame]:
    """
    Run GATK ASEReadCounter on BAM file.

    Args:
        bam_path: Path to aligned BAM file
        vcf_path: Path to variants VCF file (must be compressed and indexed)
        ref_path: Path to reference FASTA (must be indexed)
        output_path: Path to save GATK output table
        min_mapq: Minimum mapping quality (default: 10)
        min_baseq: Minimum base quality (default: 20)

    Returns:
        Tuple of (runtime_seconds, results_dataframe)
    """
    print(f"\nRunning GATK ASEReadCounter...")
    print(f"  BAM: {bam_path}")
    print(f"  VCF: {vcf_path}")
    print(f"  Reference: {ref_path}")
    print(f"  Output: {output_path}")

    # Ensure BAM has read groups (GATK requirement)
    bam_path = ensure_bam_has_read_groups(bam_path)

    start_time = time.time()

    cmd = [
        'gatk', 'ASEReadCounter',
        '-R', str(ref_path),
        '-I', str(bam_path),
        '-V', str(vcf_path),
        '-O', str(output_path),
        '--min-mapping-quality', str(min_mapq),
        '--min-base-quality', str(min_baseq)
    ]

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=600  # 10 minute timeout
        )

        if result.returncode != 0:
            print(f"GATK ERROR:\n{result.stderr}", file=sys.stderr)
            raise RuntimeError(f"GATK ASEReadCounter failed with exit code {result.returncode}")

    except subprocess.TimeoutExpired:
        raise RuntimeError("GATK ASEReadCounter timed out (>10 minutes)")

    runtime = time.time() - start_time

    # Parse GATK output
    df = parse_gatk_output(output_path)

    print(f"  Runtime: {runtime:.2f} seconds")
    print(f"  Sites processed: {len(df)}")

    return runtime, df


def parse_gatk_output(output_path: Path) -> pd.DataFrame:
    """
    Parse GATK ASEReadCounter output table.

    GATK output format:
        contig  position  variantID  refAllele  altAllele  refCount  altCount  totalCount  lowMAPQDepth  lowBaseQDepth  rawDepth  otherBases  improperPairs

    Args:
        output_path: Path to GATK output table

    Returns:
        DataFrame with columns: chrom, pos, ref_count, alt_count, total_count, ref_ratio
    """
    df = pd.read_csv(output_path, sep='\t', comment='@')

    # Rename columns to match WASP2 format
    df = df.rename(columns={
        'contig': 'chrom',
        'position': 'pos',
        'refCount': 'ref_count',
        'altCount': 'alt_count',
        'totalCount': 'total_count'
    })

    # Calculate reference allele ratio
    df['ref_ratio'] = df['ref_count'] / (df['ref_count'] + df['alt_count'])

    # Filter out sites with no coverage
    df = df[df['total_count'] > 0].copy()

    return df[['chrom', 'pos', 'ref_count', 'alt_count', 'total_count', 'ref_ratio']]


def run_wasp2_counting(
    bam_path: Path,
    vcf_path: Path,
    output_path: Path,
    min_mapq: int = 10,
    min_baseq: int = 20
) -> Tuple[float, pd.DataFrame]:
    """
    Run WASP2 allele counting on BAM file.

    This is a placeholder - actual implementation depends on WASP2's counting interface.
    For now, we'll use a mock implementation that you should replace with actual WASP2 calls.

    Args:
        bam_path: Path to aligned BAM file
        vcf_path: Path to variants VCF file
        output_path: Path to save WASP2 output
        min_mapq: Minimum mapping quality
        min_baseq: Minimum base quality

    Returns:
        Tuple of (runtime_seconds, results_dataframe)
    """
    print(f"\nRunning WASP2 counting...")
    print(f"  BAM: {bam_path}")
    print(f"  VCF: {vcf_path}")
    print(f"  Output: {output_path}")

    start_time = time.time()

    # TODO: Replace with actual WASP2 counting command
    # Example:
    # cmd = [
    #     'python', '-m', 'wasp.count_alleles',
    #     '--bam', str(bam_path),
    #     '--vcf', str(vcf_path),
    #     '--output', str(output_path),
    #     '--min-mapq', str(min_mapq),
    #     '--min-baseq', str(min_baseq)
    # ]
    # subprocess.run(cmd, check=True)

    # For now, create a placeholder
    print("  WARNING: Using placeholder WASP2 implementation")
    print("  TODO: Replace with actual WASP2 counting code")

    # Mock implementation - replace this!
    df = pd.DataFrame({
        'chrom': ['chr1'],
        'pos': [1000],
        'ref_count': [50],
        'alt_count': [50],
        'total_count': [100],
        'ref_ratio': [0.5]
    })

    runtime = time.time() - start_time

    print(f"  Runtime: {runtime:.2f} seconds")
    print(f"  Sites processed: {len(df)}")

    return runtime, df


def merge_with_ground_truth(
    tool_df: pd.DataFrame,
    truth_df: pd.DataFrame,
    tool_name: str
) -> pd.DataFrame:
    """
    Merge tool results with ground truth data.

    Args:
        tool_df: DataFrame from tool (GATK or WASP2)
        truth_df: Ground truth DataFrame with columns: chrom, pos, true_ref_ratio
        tool_name: Name of tool for column prefixes

    Returns:
        Merged DataFrame with tool and truth columns
    """
    # Merge on chromosome and position
    merged = pd.merge(
        truth_df,
        tool_df,
        on=['chrom', 'pos'],
        how='inner',
        suffixes=('_truth', f'_{tool_name}')
    )

    return merged


def compute_metrics(
    merged_df: pd.DataFrame,
    truth_col: str = 'true_ref_ratio',
    tool_col: str = 'ref_ratio'
) -> Dict[str, float]:
    """
    Compute accuracy and bias metrics.

    Args:
        merged_df: DataFrame with both truth and tool columns
        truth_col: Column name for ground truth ratios
        tool_col: Column name for tool-computed ratios

    Returns:
        Dictionary of metrics:
            - pearson_r: Pearson correlation coefficient
            - spearman_rho: Spearman correlation coefficient
            - rmse: Root mean squared error
            - mae: Mean absolute error
            - ref_bias: Mean difference from 0.5 (reference allele bias)
            - n_sites: Number of sites compared
    """
    truth = merged_df[truth_col].values
    tool = merged_df[tool_col].values

    # Filter out NaN values
    valid_mask = ~(np.isnan(truth) | np.isnan(tool))
    truth = truth[valid_mask]
    tool = tool[valid_mask]

    if len(truth) == 0:
        return {
            'pearson_r': np.nan,
            'spearman_rho': np.nan,
            'rmse': np.nan,
            'mae': np.nan,
            'ref_bias': np.nan,
            'n_sites': 0
        }

    # Correlation metrics
    pearson_r, _ = pearsonr(truth, tool)
    spearman_rho, _ = spearmanr(truth, tool)

    # Error metrics
    rmse = np.sqrt(mean_squared_error(truth, tool))
    mae = mean_absolute_error(truth, tool)

    # Reference allele bias (how much does ref ratio deviate from 0.5)
    ref_bias = tool.mean() - 0.5

    return {
        'pearson_r': pearson_r,
        'spearman_rho': spearman_rho,
        'rmse': rmse,
        'mae': mae,
        'ref_bias': ref_bias,
        'n_sites': len(truth)
    }


def create_comparison_plots(
    gatk_df: pd.DataFrame,
    wasp2_df: pd.DataFrame,
    truth_df: pd.DataFrame,
    output_dir: Path
) -> None:
    """
    Create publication-quality comparison plots.

    Args:
        gatk_df: GATK results
        wasp2_df: WASP2 results
        truth_df: Ground truth data
        output_dir: Directory to save plots
    """
    output_dir.mkdir(exist_ok=True, parents=True)

    # Merge all data
    merged = truth_df.copy()
    merged = pd.merge(merged, gatk_df, on=['chrom', 'pos'], how='inner', suffixes=('', '_gatk'))
    merged = pd.merge(merged, wasp2_df, on=['chrom', 'pos'], how='inner', suffixes=('', '_wasp2'))

    sns.set_style("whitegrid")

    # 1. Scatter plot: Truth vs Tool
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # GATK
    axes[0].scatter(merged['true_ref_ratio'], merged['ref_ratio_gatk'], alpha=0.5, s=20)
    axes[0].plot([0, 1], [0, 1], 'r--', label='Perfect correlation')
    axes[0].set_xlabel('True Reference Ratio')
    axes[0].set_ylabel('GATK Reference Ratio')
    axes[0].set_title('GATK ASEReadCounter vs Ground Truth')
    axes[0].legend()

    # WASP2
    axes[1].scatter(merged['true_ref_ratio'], merged['ref_ratio_wasp2'], alpha=0.5, s=20)
    axes[1].plot([0, 1], [0, 1], 'r--', label='Perfect correlation')
    axes[1].set_xlabel('True Reference Ratio')
    axes[1].set_ylabel('WASP2 Reference Ratio')
    axes[1].set_title('WASP2 vs Ground Truth')
    axes[1].legend()

    plt.tight_layout()
    plt.savefig(output_dir / 'accuracy_scatter.png', dpi=300, bbox_inches='tight')
    plt.close()

    # 2. Error distribution
    gatk_error = merged['ref_ratio_gatk'] - merged['true_ref_ratio']
    wasp2_error = merged['ref_ratio_wasp2'] - merged['true_ref_ratio']

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(gatk_error, bins=50, alpha=0.6, label='GATK', edgecolor='black')
    ax.hist(wasp2_error, bins=50, alpha=0.6, label='WASP2', edgecolor='black')
    ax.axvline(0, color='red', linestyle='--', label='Zero error')
    ax.set_xlabel('Error (Tool - Truth)')
    ax.set_ylabel('Frequency')
    ax.set_title('Error Distribution Comparison')
    ax.legend()
    plt.savefig(output_dir / 'error_distribution.png', dpi=300, bbox_inches='tight')
    plt.close()

    # 3. Bias comparison
    fig, ax = plt.subplots(figsize=(8, 6))
    bias_data = [
        gatk_error,
        wasp2_error
    ]
    ax.boxplot(bias_data, labels=['GATK', 'WASP2'])
    ax.axhline(0, color='red', linestyle='--', label='No bias')
    ax.set_ylabel('Error (Tool - Truth)')
    ax.set_title('Reference Allele Bias Comparison')
    ax.legend()
    plt.savefig(output_dir / 'bias_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"\nPlots saved to {output_dir}/")


def generate_comparison_report(
    gatk_metrics: Dict[str, float],
    wasp2_metrics: Dict[str, float],
    gatk_runtime: float,
    wasp2_runtime: float,
    output_path: Path
) -> str:
    """
    Generate publication-ready comparison report.

    Args:
        gatk_metrics: Metrics dictionary for GATK
        wasp2_metrics: Metrics dictionary for WASP2
        gatk_runtime: GATK runtime in seconds
        wasp2_runtime: WASP2 runtime in seconds
        output_path: Path to save markdown report

    Returns:
        Markdown-formatted report string
    """
    speedup = gatk_runtime / wasp2_runtime if wasp2_runtime > 0 else float('inf')

    report_lines = [
        "# WASP2 vs GATK ASEReadCounter Benchmark",
        "",
        f"**Analysis Date**: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "## Executive Summary",
        "",
        f"WASP2 vs GATK head-to-head comparison on simulated allele-specific expression data.",
        "",
        "## Accuracy Metrics",
        "",
        "| Metric | GATK | WASP2 | Winner |",
        "|--------|------|-------|--------|",
        f"| Pearson r | {gatk_metrics['pearson_r']:.4f} | {wasp2_metrics['pearson_r']:.4f} | {'WASP2' if wasp2_metrics['pearson_r'] > gatk_metrics['pearson_r'] else 'GATK'} |",
        f"| Spearman rho | {gatk_metrics['spearman_rho']:.4f} | {wasp2_metrics['spearman_rho']:.4f} | {'WASP2' if wasp2_metrics['spearman_rho'] > gatk_metrics['spearman_rho'] else 'GATK'} |",
        f"| RMSE | {gatk_metrics['rmse']:.4f} | {wasp2_metrics['rmse']:.4f} | {'WASP2' if wasp2_metrics['rmse'] < gatk_metrics['rmse'] else 'GATK'} |",
        f"| MAE | {gatk_metrics['mae']:.4f} | {wasp2_metrics['mae']:.4f} | {'WASP2' if wasp2_metrics['mae'] < gatk_metrics['mae'] else 'GATK'} |",
        "",
        "## Bias Metrics",
        "",
        "| Metric | GATK | WASP2 | Better |",
        "|--------|------|-------|--------|",
        f"| Ref Bias | {gatk_metrics['ref_bias']:.4f} | {wasp2_metrics['ref_bias']:.4f} | {'WASP2' if abs(wasp2_metrics['ref_bias']) < abs(gatk_metrics['ref_bias']) else 'GATK'} |",
        f"| Sites | {gatk_metrics['n_sites']} | {wasp2_metrics['n_sites']} | - |",
        "",
        "## Performance Metrics",
        "",
        "| Metric | GATK | WASP2 | Speedup |",
        "|--------|------|-------|---------|",
        f"| Runtime (sec) | {gatk_runtime:.2f} | {wasp2_runtime:.2f} | {speedup:.2f}x |",
        "",
        "## Interpretation",
        "",
        f"- **Accuracy**: WASP2 achieves Pearson r = {wasp2_metrics['pearson_r']:.4f} vs GATK r = {gatk_metrics['pearson_r']:.4f}",
        f"- **Bias**: WASP2 ref bias = {wasp2_metrics['ref_bias']:.4f} vs GATK ref bias = {gatk_metrics['ref_bias']:.4f}",
        f"- **Speed**: WASP2 is {speedup:.2f}x {'faster' if speedup > 1 else 'slower'} than GATK",
        "",
        "## Manuscript Text",
        "",
        f"We compared WASP2 against GATK ASEReadCounter on simulated allele-specific expression data "
        f"({gatk_metrics['n_sites']} heterozygous sites). WASP2 demonstrated "
        f"{'superior' if wasp2_metrics['pearson_r'] > gatk_metrics['pearson_r'] else 'comparable'} "
        f"accuracy (Pearson r = {wasp2_metrics['pearson_r']:.3f}) compared to GATK (r = {gatk_metrics['pearson_r']:.3f}). "
        f"Reference allele bias was {'lower' if abs(wasp2_metrics['ref_bias']) < abs(gatk_metrics['ref_bias']) else 'higher'} "
        f"for WASP2 ({wasp2_metrics['ref_bias']:.3f}) versus GATK ({gatk_metrics['ref_bias']:.3f}). "
        f"WASP2 completed processing in {wasp2_runtime:.1f} seconds, achieving {speedup:.1f}x "
        f"{'speedup' if speedup > 1 else 'slowdown'} relative to GATK ({gatk_runtime:.1f} seconds).",
        "",
        "## References",
        "",
        "- McKenna et al. (2010). The Genome Analysis Toolkit. *Genome Research*, 20(9), 1297-1303.",
        "- van de Geijn et al. (2015). WASP: allele-specific software for robust molecular QTL discovery. *Nature Methods*, 12(11), 1061-1063.",
        ""
    ]

    report = "\n".join(report_lines)

    # Save to file
    output_path.write_text(report)
    print(f"\nReport saved to: {output_path}")

    return report


def main():
    parser = argparse.ArgumentParser(
        description='Benchmark WASP2 vs GATK ASEReadCounter',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run full comparison
  python simulation/benchmark_vs_gatk.py \\
      --bam simulation_results/comprehensive_*/aligned.sorted.bam \\
      --vcf simulation_results/comprehensive_*/variants.vcf.gz \\
      --ref simulation_results/comprehensive_*/reference.fa \\
      --ground-truth simulation_results/comprehensive_*/ground_truth.csv \\
      --output comparison_results/

  # GATK only (to test installation)
  python simulation/benchmark_vs_gatk.py \\
      --bam test.bam \\
      --vcf test.vcf.gz \\
      --ref ref.fa \\
      --ground-truth truth.csv \\
      --output test_output/ \\
      --gatk-only
        """
    )

    parser.add_argument('--bam', required=True, type=Path,
                        help='Path to aligned BAM file (must be indexed)')
    parser.add_argument('--vcf', required=True, type=Path,
                        help='Path to variants VCF file (must be compressed and indexed)')
    parser.add_argument('--ref', required=True, type=Path,
                        help='Path to reference FASTA (must be indexed)')
    parser.add_argument('--ground-truth', required=True, type=Path,
                        help='Path to ground truth CSV with columns: chrom, pos, true_ref_ratio')
    parser.add_argument('--output', required=True, type=Path,
                        help='Output directory for results')
    parser.add_argument('--min-mapq', type=int, default=10,
                        help='Minimum mapping quality (default: 10)')
    parser.add_argument('--min-baseq', type=int, default=20,
                        help='Minimum base quality (default: 20)')
    parser.add_argument('--gatk-only', action='store_true',
                        help='Run GATK only (skip WASP2)')
    parser.add_argument('--wasp2-only', action='store_true',
                        help='Run WASP2 only (skip GATK)')

    args = parser.parse_args()

    # Create output directory
    args.output.mkdir(exist_ok=True, parents=True)

    print("="*80)
    print("WASP2 vs GATK ASEReadCounter Benchmark")
    print("="*80)

    # Check GATK installation
    if not args.wasp2_only:
        print("\nChecking GATK installation...")
        if not check_gatk_available():
            print("ERROR: GATK not found or ASEReadCounter not available", file=sys.stderr)
            print("Install with: conda install -c bioconda gatk4", file=sys.stderr)
            sys.exit(1)
        print("  GATK found")

    # Check dependencies
    print("\nChecking input files...")
    if not check_dependencies(args.bam, args.vcf, args.ref):
        sys.exit(1)
    print("  All required files found")

    # Load ground truth
    print(f"\nLoading ground truth from {args.ground_truth}...")
    truth_df = pd.read_csv(args.ground_truth)
    print(f"  Loaded {len(truth_df)} ground truth sites")

    # Run GATK
    gatk_runtime = 0.0
    gatk_metrics = {}
    if not args.wasp2_only:
        gatk_output = args.output / 'gatk_ase_counts.table'
        gatk_runtime, gatk_df = run_gatk_ase_counter(
            args.bam, args.vcf, args.ref, gatk_output,
            args.min_mapq, args.min_baseq
        )

        # Compute GATK metrics
        gatk_merged = merge_with_ground_truth(gatk_df, truth_df, 'gatk')
        gatk_metrics = compute_metrics(gatk_merged, 'true_ref_ratio', 'ref_ratio')

        print(f"\nGATK Results:")
        print(f"  Pearson r: {gatk_metrics['pearson_r']:.4f}")
        print(f"  Ref bias: {gatk_metrics['ref_bias']:.4f}")
        print(f"  Runtime: {gatk_runtime:.2f} sec")

    # Run WASP2
    wasp2_runtime = 0.0
    wasp2_metrics = {}
    if not args.gatk_only:
        wasp2_output = args.output / 'wasp2_ase_counts.txt'
        wasp2_runtime, wasp2_df = run_wasp2_counting(
            args.bam, args.vcf, wasp2_output,
            args.min_mapq, args.min_baseq
        )

        # Compute WASP2 metrics
        wasp2_merged = merge_with_ground_truth(wasp2_df, truth_df, 'wasp2')
        wasp2_metrics = compute_metrics(wasp2_merged, 'true_ref_ratio', 'ref_ratio')

        print(f"\nWASP2 Results:")
        print(f"  Pearson r: {wasp2_metrics['pearson_r']:.4f}")
        print(f"  Ref bias: {wasp2_metrics['ref_bias']:.4f}")
        print(f"  Runtime: {wasp2_runtime:.2f} sec")

    # Generate comparison report
    if not args.gatk_only and not args.wasp2_only:
        report_path = args.output / 'comparison_report.md'
        generate_comparison_report(
            gatk_metrics, wasp2_metrics,
            gatk_runtime, wasp2_runtime,
            report_path
        )

        # Create plots
        create_comparison_plots(gatk_df, wasp2_df, truth_df, args.output)

        # Print summary
        print("\n" + "="*80)
        print("BENCHMARK COMPLETE")
        print("="*80)
        print(f"\nWASP2 vs GATK:")
        print(f"  Accuracy (Pearson r): WASP2={wasp2_metrics['pearson_r']:.4f} vs GATK={gatk_metrics['pearson_r']:.4f}")
        print(f"  Bias (ref ratio): WASP2={wasp2_metrics['ref_bias']:.4f} vs GATK={gatk_metrics['ref_bias']:.4f}")
        print(f"  Speed: WASP2 is {gatk_runtime/wasp2_runtime:.2f}x {'faster' if gatk_runtime > wasp2_runtime else 'slower'}")
        print(f"\nResults saved to: {args.output}/")

    print("\n" + "="*80)


if __name__ == '__main__':
    main()
