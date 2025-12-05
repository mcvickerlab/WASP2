#!/usr/bin/env python3
"""
Benchmark WASP2 vs GATK ASEReadCounter on simulation results.

This script:
1. Runs GATK ASEReadCounter (or uses existing results)
2. Runs WASP2 allele counting
3. Compares both against ground truth
4. Reports metrics: Pearson r, RMSE, MAE, REF bias, and runtime
"""

import argparse
import pandas as pd
import numpy as np
import subprocess
import time
import os
from pathlib import Path
from scipy import stats


def run_gatk_asereadcounter(bam_path, vcf_path, ref_path, output_path):
    """Run GATK ASEReadCounter."""
    print("Running GATK ASEReadCounter...")
    start_time = time.time()

    cmd = [
        "gatk", "ASEReadCounter",
        "-R", ref_path,
        "-I", bam_path,
        "-V", vcf_path,
        "-O", output_path,
        "--min-depth-of-non-filtered-base", "0",
        "--min-mapping-quality", "10",
        "--min-base-quality", "2"
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        runtime = time.time() - start_time
        print(f"GATK completed in {runtime:.2f} seconds")
        return runtime
    except subprocess.CalledProcessError as e:
        print(f"GATK error: {e.stderr}")
        raise


def run_wasp2_counting(bam_path, vcf_path, output_file):
    """Run WASP2 allele counting."""
    print("Running WASP2 allele counting...")
    start_time = time.time()

    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "python", "-m", "src.counting",
        "count-variants",
        bam_path,
        vcf_path,
        "--out", str(output_file)
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        runtime = time.time() - start_time
        print(f"WASP2 completed in {runtime:.2f} seconds")
        return runtime
    except subprocess.CalledProcessError as e:
        print(f"WASP2 error: {e.stderr}")
        raise


def parse_gatk_output(gatk_file):
    """Parse GATK ASEReadCounter output."""
    print(f"Parsing GATK output: {gatk_file}")

    # GATK output is tab-delimited
    df = pd.read_csv(gatk_file, sep='\t', comment='#')

    # Rename columns to match our format
    df = df.rename(columns={
        'contig': 'chrom',
        'position': 'pos',
        'variantID': 'variant_id',
        'refAllele': 'ref_allele',
        'altAllele': 'alt_allele',
        'refCount': 'ref_count',
        'altCount': 'alt_count',
        'totalCount': 'total_count'
    })

    # Calculate ref ratio
    df['ref_ratio'] = df['ref_count'] / (df['ref_count'] + df['alt_count'])
    df['ref_ratio'] = df['ref_ratio'].fillna(0.5)

    print(f"Parsed {len(df)} entries from GATK")
    return df[['chrom', 'pos', 'variant_id', 'ref_count', 'alt_count', 'total_count', 'ref_ratio']]


def parse_wasp2_output(wasp2_file):
    """Parse WASP2 output."""
    print(f"Parsing WASP2 output: {wasp2_file}")

    # Parse WASP2 output format (TSV file)
    df = pd.read_csv(wasp2_file, sep='\t')

    print(f"Available columns: {df.columns.tolist()}")

    # Ensure we have the right columns
    required_cols = ['chrom', 'pos', 'ref_count', 'alt_count']
    if not all(col in df.columns for col in required_cols):
        raise ValueError(f"Missing required columns in WASP2 output. Expected {required_cols}, got {df.columns.tolist()}")

    # Calculate total and ratio
    df['total_count'] = df['ref_count'] + df['alt_count']
    df['ref_ratio'] = df['ref_count'] / df['total_count']
    df['ref_ratio'] = df['ref_ratio'].fillna(0.5)

    print(f"Parsed {len(df)} entries from WASP2")
    return df[['chrom', 'pos', 'ref_count', 'alt_count', 'total_count', 'ref_ratio']]


def parse_ground_truth(ground_truth_file):
    """Parse ground truth CSV."""
    print(f"Parsing ground truth: {ground_truth_file}")

    df = pd.read_csv(ground_truth_file)

    # Ensure we have required columns
    required_cols = ['chrom', 'pos', 'true_ref_ratio']
    if not all(col in df.columns for col in required_cols):
        print(f"Available columns: {df.columns.tolist()}")
        raise ValueError(f"Missing required columns in ground truth")

    print(f"Parsed {len(df)} ground truth entries")
    return df[['chrom', 'pos', 'variant_id', 'true_ref_ratio']]


def calculate_metrics(predicted, true, label="Tool"):
    """Calculate accuracy metrics."""
    # Remove NaN values
    mask = ~(np.isnan(predicted) | np.isnan(true))
    predicted = predicted[mask]
    true = true[mask]

    if len(predicted) == 0:
        print(f"Warning: No valid data points for {label}")
        return {
            'pearson_r': np.nan,
            'pearson_p': np.nan,
            'rmse': np.nan,
            'mae': np.nan,
            'ref_bias': np.nan,
            'n_variants': 0
        }

    # Pearson correlation
    pearson_r, pearson_p = stats.pearsonr(predicted, true)

    # RMSE
    rmse = np.sqrt(np.mean((predicted - true) ** 2))

    # MAE
    mae = np.mean(np.abs(predicted - true))

    # REF bias (mean difference)
    ref_bias = np.mean(predicted - true)

    return {
        'pearson_r': pearson_r,
        'pearson_p': pearson_p,
        'rmse': rmse,
        'mae': mae,
        'ref_bias': ref_bias,
        'n_variants': len(predicted)
    }


def compare_tools(gatk_df, wasp2_df, ground_truth_df, gatk_runtime, wasp2_runtime):
    """Compare GATK and WASP2 against ground truth."""
    print("\n" + "="*80)
    print("COMPARISON RESULTS")
    print("="*80)

    # Merge with ground truth on chrom and pos
    gatk_merged = pd.merge(
        ground_truth_df,
        gatk_df,
        on=['chrom', 'pos'],
        how='inner',
        suffixes=('_truth', '_gatk')
    )

    wasp2_merged = pd.merge(
        ground_truth_df,
        wasp2_df,
        on=['chrom', 'pos'],
        how='inner',
        suffixes=('_truth', '_wasp2')
    )

    print(f"\nGATK: Matched {len(gatk_merged)} / {len(ground_truth_df)} ground truth variants")
    print(f"WASP2: Matched {len(wasp2_merged)} / {len(ground_truth_df)} ground truth variants")

    # Calculate metrics for GATK
    gatk_metrics = calculate_metrics(
        gatk_merged['ref_ratio'].values,
        gatk_merged['true_ref_ratio'].values,
        "GATK"
    )

    # Calculate metrics for WASP2
    wasp2_metrics = calculate_metrics(
        wasp2_merged['ref_ratio'].values,
        wasp2_merged['true_ref_ratio'].values,
        "WASP2"
    )

    # Print results
    print("\n" + "-"*80)
    print("GATK ASEReadCounter Results:")
    print("-"*80)
    print(f"  Pearson r:      {gatk_metrics['pearson_r']:.4f} (p={gatk_metrics['pearson_p']:.2e})")
    print(f"  RMSE:           {gatk_metrics['rmse']:.4f}")
    print(f"  MAE:            {gatk_metrics['mae']:.4f}")
    print(f"  REF bias:       {gatk_metrics['ref_bias']:+.4f}")
    print(f"  Runtime:        {gatk_runtime:.2f} seconds")
    print(f"  Variants:       {gatk_metrics['n_variants']}")

    print("\n" + "-"*80)
    print("WASP2 Results:")
    print("-"*80)
    print(f"  Pearson r:      {wasp2_metrics['pearson_r']:.4f} (p={wasp2_metrics['pearson_p']:.2e})")
    print(f"  RMSE:           {wasp2_metrics['rmse']:.4f}")
    print(f"  MAE:            {wasp2_metrics['mae']:.4f}")
    print(f"  REF bias:       {wasp2_metrics['ref_bias']:+.4f}")
    print(f"  Runtime:        {wasp2_runtime:.2f} seconds")
    print(f"  Variants:       {wasp2_metrics['n_variants']}")

    print("\n" + "-"*80)
    print("Winner Summary:")
    print("-"*80)

    # Determine winners
    winners = {
        'Pearson r': 'WASP2' if wasp2_metrics['pearson_r'] > gatk_metrics['pearson_r'] else 'GATK',
        'RMSE': 'WASP2' if wasp2_metrics['rmse'] < gatk_metrics['rmse'] else 'GATK',
        'MAE': 'WASP2' if wasp2_metrics['mae'] < gatk_metrics['mae'] else 'GATK',
        'REF bias (abs)': 'WASP2' if abs(wasp2_metrics['ref_bias']) < abs(gatk_metrics['ref_bias']) else 'GATK',
        'Runtime': 'WASP2' if wasp2_runtime < gatk_runtime else 'GATK'
    }

    for metric, winner in winners.items():
        symbol = ">" if winner == "WASP2" else "<"
        print(f"  {metric:20s}: {'WASP2':6s} {symbol} {'GATK':6s}")

    print("="*80)

    return {
        'gatk': gatk_metrics,
        'wasp2': wasp2_metrics,
        'gatk_runtime': gatk_runtime,
        'wasp2_runtime': wasp2_runtime,
        'winners': winners
    }


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark WASP2 vs GATK ASEReadCounter on simulation results"
    )
    parser.add_argument("--bam", required=True, help="Input BAM file")
    parser.add_argument("--vcf", required=True, help="Input VCF file")
    parser.add_argument("--ref", required=True, help="Reference FASTA file")
    parser.add_argument("--ground-truth", required=True, help="Ground truth CSV file")
    parser.add_argument("--output-dir", default="comparison_results", help="Output directory")
    parser.add_argument("--gatk-output", help="Existing GATK output file (skip running GATK)")
    parser.add_argument("--wasp2-output", help="Existing WASP2 output directory (skip running WASP2)")

    args = parser.parse_args()

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Run or load GATK results
    if args.gatk_output and os.path.exists(args.gatk_output):
        print(f"Using existing GATK output: {args.gatk_output}")
        gatk_output = args.gatk_output
        gatk_runtime = 0.0  # Unknown runtime
    else:
        gatk_output = output_dir / "gatk_ase_counts.table"
        gatk_runtime = run_gatk_asereadcounter(
            args.bam, args.vcf, args.ref, str(gatk_output)
        )

    # Run or load WASP2 results
    if args.wasp2_output and os.path.exists(args.wasp2_output):
        print(f"Using existing WASP2 output: {args.wasp2_output}")
        wasp2_output_file = Path(args.wasp2_output)
        wasp2_runtime = 0.0  # Unknown runtime
    else:
        wasp2_output_file = output_dir / "wasp2_counts.tsv"
        wasp2_runtime = run_wasp2_counting(
            args.bam, args.vcf, str(wasp2_output_file)
        )

    # Parse results
    gatk_df = parse_gatk_output(gatk_output)
    wasp2_df = parse_wasp2_output(wasp2_output_file)
    ground_truth_df = parse_ground_truth(args.ground_truth)

    # Compare
    results = compare_tools(gatk_df, wasp2_df, ground_truth_df, gatk_runtime, wasp2_runtime)

    # Save detailed results
    results_file = output_dir / "comparison_results.txt"
    with open(results_file, 'w') as f:
        f.write("GATK vs WASP2 Benchmark Results\n")
        f.write("="*80 + "\n\n")
        f.write(f"GATK Pearson r: {results['gatk']['pearson_r']:.4f}\n")
        f.write(f"GATK RMSE: {results['gatk']['rmse']:.4f}\n")
        f.write(f"GATK MAE: {results['gatk']['mae']:.4f}\n")
        f.write(f"GATK REF bias: {results['gatk']['ref_bias']:+.4f}\n")
        f.write(f"GATK Runtime: {results['gatk_runtime']:.2f}s\n\n")
        f.write(f"WASP2 Pearson r: {results['wasp2']['pearson_r']:.4f}\n")
        f.write(f"WASP2 RMSE: {results['wasp2']['rmse']:.4f}\n")
        f.write(f"WASP2 MAE: {results['wasp2']['mae']:.4f}\n")
        f.write(f"WASP2 REF bias: {results['wasp2']['ref_bias']:+.4f}\n")
        f.write(f"WASP2 Runtime: {results['wasp2_runtime']:.2f}s\n")

    print(f"\nResults saved to: {results_file}")


if __name__ == "__main__":
    main()
