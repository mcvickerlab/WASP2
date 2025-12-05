#!/usr/bin/env python3
"""
Simple GATK vs WASP2 comparison using direct BAM counting.

This script:
1. Reads ground truth and GATK results
2. Runs WASP2 Rust BAM counter directly on variants
3. Compares both against ground truth
4. Reports metrics: Pearson r, RMSE, MAE, REF bias
"""

import argparse
import pandas as pd
import numpy as np
import time
from pathlib import Path
from scipy import stats
import pysam


# Import Rust acceleration
try:
    from wasp2_rust import BamCounter as RustBamCounter
    RUST_AVAILABLE = True
except ImportError:
    RUST_AVAILABLE = False
    print("Warning: Rust extension not available")


def parse_gatk_output(gatk_file):
    """Parse GATK ASEReadCounter output."""
    print(f"Parsing GATK output: {gatk_file}")
    df = pd.read_csv(gatk_file, sep='\t', comment='#')

    df = df.rename(columns={
        'contig': 'chrom',
        'position': 'pos',
        'variantID': 'variant_id',
        'refCount': 'ref_count',
        'altCount': 'alt_count',
        'totalCount': 'total_count'
    })

    df['ref_ratio'] = df['ref_count'] / (df['ref_count'] + df['alt_count'])
    df['ref_ratio'] = df['ref_ratio'].fillna(0.5)

    print(f"Parsed {len(df)} entries from GATK")
    return df


def run_wasp2_rust_counting(bam_path, vcf_path):
    """Run WASP2 Rust BAM counter on VCF variants."""
    if not RUST_AVAILABLE:
        raise RuntimeError("Rust extension not available")

    print(f"Running WASP2 Rust counting on {bam_path}...")
    start_time = time.time()

    # Read VCF to get variants
    vcf = pysam.VariantFile(vcf_path)
    variants = []

    for record in vcf:
        # Get reference and first alt allele
        ref = record.ref
        alt = record.alts[0] if record.alts else ref

        # Store variant info
        variants.append({
            'chrom': record.chrom,
            'pos': record.pos,
            'ref': ref,
            'alt': alt,
            'variant_id': record.id
        })

    vcf.close()
    print(f"Found {len(variants)} variants in VCF")

    # Create Rust BAM counter
    counter = RustBamCounter(bam_path)

    # Prepare regions for Rust
    regions = [(v['chrom'], v['pos'], v['ref'], v['alt']) for v in variants]

    # Count alleles
    print("Counting alleles with Rust...")
    counts = counter.count_alleles(regions, min_qual=0, threads=1)

    # Combine results
    results = []
    for variant, (ref_count, alt_count, other_count) in zip(variants, counts):
        total = ref_count + alt_count
        ref_ratio = ref_count / total if total > 0 else 0.5

        results.append({
            'chrom': variant['chrom'],
            'pos': variant['pos'],
            'variant_id': variant['variant_id'],
            'ref_count': ref_count,
            'alt_count': alt_count,
            'other_count': other_count,
            'total_count': total,
            'ref_ratio': ref_ratio
        })

    runtime = time.time() - start_time
    print(f"WASP2 counting completed in {runtime:.2f} seconds")

    return pd.DataFrame(results), runtime


def parse_ground_truth(ground_truth_file):
    """Parse ground truth CSV."""
    print(f"Parsing ground truth: {ground_truth_file}")
    df = pd.read_csv(ground_truth_file)
    print(f"Parsed {len(df)} ground truth entries")
    return df


def calculate_metrics(predicted, true, label="Tool"):
    """Calculate accuracy metrics."""
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

    pearson_r, pearson_p = stats.pearsonr(predicted, true)
    rmse = np.sqrt(np.mean((predicted - true) ** 2))
    mae = np.mean(np.abs(predicted - true))
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

    # Merge with ground truth
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

    # Show some sample data
    print("\nSample GATK results:")
    print(gatk_merged[['chrom', 'pos', 'variant_id_truth', 'true_ref_ratio', 'ref_ratio']].head())

    print("\nSample WASP2 results:")
    print(wasp2_merged[['chrom', 'pos', 'variant_id_truth', 'true_ref_ratio', 'ref_ratio']].head())

    # Calculate metrics
    gatk_metrics = calculate_metrics(
        gatk_merged['ref_ratio'].values,
        gatk_merged['true_ref_ratio'].values,
        "GATK"
    )

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
    print(f"  Runtime:        {gatk_runtime:.2f} seconds (from existing run)")
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
        description="Simple GATK vs WASP2 comparison"
    )
    parser.add_argument("--bam", required=True, help="Input BAM file")
    parser.add_argument("--vcf", required=True, help="Input VCF file")
    parser.add_argument("--ground-truth", required=True, help="Ground truth CSV file")
    parser.add_argument("--gatk-output", required=True, help="GATK output file")
    parser.add_argument("--output-dir", default="comparison_results", help="Output directory")

    args = parser.parse_args()

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Parse GATK results (assume already run)
    gatk_df = parse_gatk_output(args.gatk_output)
    gatk_runtime = 0.0  # Unknown

    # Run WASP2 counting
    wasp2_df, wasp2_runtime = run_wasp2_rust_counting(args.bam, args.vcf)

    # Save WASP2 results
    wasp2_output = output_dir / "wasp2_counts.tsv"
    wasp2_df.to_csv(wasp2_output, sep='\t', index=False)
    print(f"\nWASP2 results saved to: {wasp2_output}")

    # Parse ground truth
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
