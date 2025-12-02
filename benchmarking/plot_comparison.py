#!/usr/bin/env python3
"""
Plot WASP2 Rust vs Python performance comparison.
"""
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path


def load_rust_data(filepath):
    """Load Rust benchmark results."""
    df = pd.read_csv(filepath, sep='\t',
                     names=['n_reads', 'seed', 'total', 'intersect', 'remap', 'filter'])
    return df.groupby('n_reads').agg({
        'total': 'mean',
        'intersect': 'mean',
        'remap': 'mean',
        'filter': 'mean'
    }).reset_index()


def load_aho_python_data(filepath):
    """Load aho's Python benchmark results."""
    df = pd.read_csv(filepath, sep='\t',
                     names=['n_reads', 'seed', 'total', 'intersect', 'remap', 'filter'])
    return df.groupby('n_reads').agg({
        'total': 'mean',
        'intersect': 'mean',
        'remap': 'mean',
        'filter': 'mean'
    }).reset_index()


def load_wasp1_data(filepath):
    """Load WASP1 benchmark results."""
    df = pd.read_csv(filepath, sep='\t',
                     names=['n_reads', 'seed', 'total', 'snp2h5', 'intersect', 'remap', 'filter', 'extra'])
    return df.groupby('n_reads').agg({
        'total': 'mean',
        'intersect': 'mean',
        'remap': 'mean',
        'filter': 'mean'
    }).reset_index()


def plot_scaling_comparison(rust_df, python_df, output_dir='plots', wasp1_df=None):
    """Plot scaling comparison."""
    Path(output_dir).mkdir(exist_ok=True)

    # Merge dataframes
    merged = rust_df.merge(python_df, on='n_reads', suffixes=('_rust', '_python'))
    if wasp1_df is not None:
        merged = merged.merge(wasp1_df[['n_reads', 'total']], on='n_reads', how='left')
        merged = merged.rename(columns={'total': 'total_wasp1'})

    # Convert to millions
    merged['n_reads_m'] = merged['n_reads'] / 1e6

    # Plot 1: Total time comparison
    fig, ax = plt.subplots(figsize=(10, 6))

    if wasp1_df is not None and 'total_wasp1' in merged.columns:
        ax.plot(merged['n_reads_m'], merged['total_wasp1'], 'o-', color='gray',
                label='WASP1 (original)', alpha=0.7, markersize=3)

    ax.plot(merged['n_reads_m'], merged['total_python'], 'o-', color='blue',
            label='WASP2 Python', alpha=0.7, markersize=3)
    ax.plot(merged['n_reads_m'], merged['total_rust'], 'o-', color='red',
            label='WASP2 Rust', alpha=0.7, markersize=3)

    ax.set_xlabel('Number of Reads (millions)', fontsize=12)
    ax.set_ylabel('Total Time (seconds)', fontsize=12)
    ax.set_title('WASP2 Performance: Rust vs Python', fontsize=14)
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f'{output_dir}/total_time_comparison.png', dpi=150)
    print(f"Saved: {output_dir}/total_time_comparison.png")
    plt.close()

    # Plot 2: Speedup ratio
    fig, ax = plt.subplots(figsize=(10, 6))

    merged['speedup'] = merged['total_python'] / merged['total_rust']
    ax.plot(merged['n_reads_m'], merged['speedup'], 'o-', color='green', markersize=4)
    ax.axhline(y=merged['speedup'].mean(), color='red', linestyle='--',
               label=f'Average: {merged["speedup"].mean():.1f}x')

    ax.set_xlabel('Number of Reads (millions)', fontsize=12)
    ax.set_ylabel('Speedup (Python time / Rust time)', fontsize=12)
    ax.set_title('WASP2 Rust Speedup over Python', fontsize=14)
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f'{output_dir}/speedup_ratio.png', dpi=150)
    print(f"Saved: {output_dir}/speedup_ratio.png")
    plt.close()

    # Plot 3: Stage breakdown
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    stages = ['intersect', 'remap', 'filter']
    titles = ['Intersect Stage', 'Remap Stage', 'Filter Stage']

    for ax, stage, title in zip(axes, stages, titles):
        ax.plot(merged['n_reads_m'], merged[f'{stage}_python'], 'o-',
                color='blue', label='Python', markersize=3)
        ax.plot(merged['n_reads_m'], merged[f'{stage}_rust'], 'o-',
                color='red', label='Rust', markersize=3)
        ax.set_xlabel('Reads (millions)')
        ax.set_ylabel('Time (seconds)')
        ax.set_title(title)
        ax.legend()
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f'{output_dir}/stage_breakdown_comparison.png', dpi=150)
    print(f"Saved: {output_dir}/stage_breakdown_comparison.png")
    plt.close()

    return merged


def main():
    parser = argparse.ArgumentParser(description='Plot WASP2 benchmark comparison')
    parser.add_argument('--rust', required=True, help='Rust benchmark results TSV')
    parser.add_argument('--python', help='Python benchmark results TSV')
    parser.add_argument('--compare-aho', action='store_true',
                        help='Use aho\'s benchmark data for Python')
    parser.add_argument('--include-wasp1', action='store_true',
                        help='Include WASP1 (original) in comparison')
    parser.add_argument('--output', default='plots', help='Output directory')

    args = parser.parse_args()

    # Load Rust data
    rust_df = load_rust_data(args.rust)

    # Load Python data
    # Use the CORRECT files from aho's preprint (plot_performance_v1.ipynb):
    # - WASP2: wasp2_perf_logs_7908428-Copy1.txt (multithreaded version)
    # - WASP1: wasp1_perf_logs_snp2h5_7881649-Copy1.txt
    if args.compare_aho:
        aho_python = '/iblm/netapp/home/aho/projects/wasp/testing/performance/outputs/test_logs_v1/wasp2_perf_logs_7908428-Copy1.txt'
        python_df = load_aho_python_data(aho_python)
    elif args.python:
        python_df = load_aho_python_data(args.python)
    else:
        print("Error: Need either --compare-aho or --python argument")
        return

    # Load WASP1 data if requested
    wasp1_df = None
    if args.include_wasp1:
        wasp1_path = '/iblm/netapp/home/aho/projects/wasp/testing/performance/outputs/test_logs_v1/wasp1_perf_logs_snp2h5_7881649-Copy1.txt'
        wasp1_df = load_wasp1_data(wasp1_path)

    # Generate plots
    merged = plot_scaling_comparison(rust_df, python_df, args.output, wasp1_df)

    # Print summary
    print("\n" + "=" * 60)
    print("PERFORMANCE SUMMARY")
    print("=" * 60)

    avg_speedup = (merged['total_python'] / merged['total_rust']).mean()
    print(f"\nAverage speedup: {avg_speedup:.2f}x")

    print("\nBy read count:")
    for _, row in merged.iterrows():
        n = int(row['n_reads'] / 1e6)
        speedup = row['total_python'] / row['total_rust']
        print(f"  {n:3d}M reads: Python {row['total_python']:7.1f}s vs Rust {row['total_rust']:7.1f}s ({speedup:.1f}x speedup)")


if __name__ == '__main__':
    main()
