#!/usr/bin/env python3
"""
Plot WASP performance comparison: WASP1 vs WASP2-Python vs WASP2-Rust (SNPs) vs WASP2-Rust (+INDELs)
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# Data sources (local copies for reproducibility)
WASP1_FILE = 'benchmarking/data/aho_perf_logs/wasp1_perf_logs_snp2h5_7881649-Copy1.txt'
WASP2_PYTHON_FILE = 'benchmarking/data/aho_perf_logs/wasp2_perf_logs_7908428-Copy1.txt'
RUST_SNP_FILE = 'benchmarking/results/wasp2_rust_scale_8512868.tsv'
RUST_INDEL_FILE = 'benchmarking/results/wasp2_rust_indels_scale_8558675.tsv'  # Full 150 tasks (1M increments)

def load_wasp1(filepath):
    """WASP1: n_reads, seed, total, snp2h5, intersect, remap, filter, extra"""
    df = pd.read_csv(filepath, sep='\t',
                     names=['n_reads', 'seed', 'total', 'snp2h5', 'intersect', 'remap', 'filter', 'extra'])
    return df.groupby('n_reads')['total'].mean().reset_index()

def load_wasp2_python(filepath):
    """WASP2 Python: n_reads, seed, total, intersect, remap, filter"""
    df = pd.read_csv(filepath, sep='\t',
                     names=['n_reads', 'seed', 'total', 'intersect', 'remap', 'filter'])
    return df.groupby('n_reads')['total'].mean().reset_index()

def load_rust(filepath):
    """Rust benchmark: n_reads, seed, total, make_reads, bwa, filter"""
    df = pd.read_csv(filepath, sep='\t',
                     names=['n_reads', 'seed', 'total', 'make_reads', 'bwa', 'filter'])
    return df.groupby('n_reads')['total'].mean().reset_index()

def main():
    output_dir = Path('plots')
    output_dir.mkdir(exist_ok=True)

    # Load all data
    wasp1 = load_wasp1(WASP1_FILE)
    wasp2_py = load_wasp2_python(WASP2_PYTHON_FILE)
    rust_snp = load_rust(RUST_SNP_FILE)

    # Try to load indel data (may not exist yet)
    try:
        rust_indel = load_rust(RUST_INDEL_FILE)
        has_indel_data = True
    except FileNotFoundError:
        rust_indel = None
        has_indel_data = False
        print(f"Note: Indel data file not found ({RUST_INDEL_FILE})")

    # Convert to millions
    wasp1['n_reads_m'] = wasp1['n_reads'] / 1e6
    wasp2_py['n_reads_m'] = wasp2_py['n_reads'] / 1e6
    rust_snp['n_reads_m'] = rust_snp['n_reads'] / 1e6
    if has_indel_data:
        rust_indel['n_reads_m'] = rust_indel['n_reads'] / 1e6

    # ============================================
    # Plot 1: Full comparison (WASP1, WASP2-Python, WASP2-Rust SNP, +INDELs)
    # ============================================
    fig, ax = plt.subplots(figsize=(12, 7))

    ax.plot(wasp1['n_reads_m'], wasp1['total'], 'o-', color='gray',
            label='WASP1', linewidth=2, markersize=5, alpha=0.8)
    ax.plot(wasp2_py['n_reads_m'], wasp2_py['total'], 's-', color='blue',
            label='WASP2-Python', linewidth=2, markersize=5, alpha=0.8)
    ax.plot(rust_snp['n_reads_m'], rust_snp['total'], 'D-', color='red',
            label='WASP2-Rust (SNPs)', linewidth=2, markersize=5, alpha=0.8)
    if has_indel_data:
        ax.plot(rust_indel['n_reads_m'], rust_indel['total'], '^-', color='forestgreen',
                label='WASP2-Rust (+INDELs)', linewidth=2, markersize=5, alpha=0.8)

    ax.set_xlabel('Number of Reads (millions)', fontsize=14)
    ax.set_ylabel('Total Time (seconds)', fontsize=14)
    ax.set_title('WASP Performance Comparison', fontsize=16)
    ax.legend(loc='upper left', fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, max(rust_snp['n_reads_m'].max(), 150) + 5)

    plt.tight_layout()
    plt.savefig(output_dir / 'wasp_all_versions_comparison.png', dpi=150)
    print(f"Saved: {output_dir}/wasp_all_versions_comparison.png")
    plt.close()

    # ============================================
    # Plot 2: Zoomed in on WASP2 versions
    # ============================================
    fig, ax = plt.subplots(figsize=(12, 7))

    max_reads = rust_snp['n_reads_m'].max()
    py_subset = wasp2_py[wasp2_py['n_reads_m'] <= max_reads]

    ax.plot(py_subset['n_reads_m'], py_subset['total'], 's-', color='blue',
            label='WASP2-Python', linewidth=2, markersize=6, alpha=0.8)
    ax.plot(rust_snp['n_reads_m'], rust_snp['total'], 'D-', color='red',
            label='WASP2-Rust (SNPs)', linewidth=2, markersize=6, alpha=0.8)
    if has_indel_data:
        ax.plot(rust_indel['n_reads_m'], rust_indel['total'], '^-', color='forestgreen',
                label='WASP2-Rust (+INDELs)', linewidth=2, markersize=6, alpha=0.8)

    ax.set_xlabel('Number of Reads (millions)', fontsize=14)
    ax.set_ylabel('Total Time (seconds)', fontsize=14)
    ax.set_title('WASP2 Performance: Python vs Rust', fontsize=16)
    ax.legend(loc='upper left', fontsize=12)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_dir / 'rust_versions_comparison.png', dpi=150)
    print(f"Saved: {output_dir}/rust_versions_comparison.png")
    plt.close()

    # ============================================
    # Plot 3: Speedup factors
    # ============================================
    fig, ax = plt.subplots(figsize=(12, 7))

    # Merge on common read counts
    merged = rust_snp.merge(wasp2_py, on='n_reads', suffixes=('_rust', '_py'))
    merged = merged.merge(wasp1, on='n_reads')
    merged.rename(columns={'total': 'total_wasp1'}, inplace=True)

    merged['speedup_vs_python'] = merged['total_py'] / merged['total_rust']
    merged['speedup_vs_wasp1'] = merged['total_wasp1'] / merged['total_rust']
    merged['n_reads_m'] = merged['n_reads'] / 1e6

    ax.plot(merged['n_reads_m'], merged['speedup_vs_wasp1'], 'o-', color='gray',
            label=f'vs WASP1 (avg {merged["speedup_vs_wasp1"].mean():.1f}x)', linewidth=2, markersize=5)
    ax.plot(merged['n_reads_m'], merged['speedup_vs_python'], 's-', color='blue',
            label=f'vs WASP2-Python (avg {merged["speedup_vs_python"].mean():.1f}x)', linewidth=2, markersize=5)

    ax.axhline(y=1, color='black', linestyle='--', alpha=0.5)
    ax.set_xlabel('Number of Reads (millions)', fontsize=14)
    ax.set_ylabel('Speedup Factor (X times faster)', fontsize=14)
    ax.set_title('WASP2-Rust Speedup vs Other Versions', fontsize=16)
    ax.legend(loc='upper right', fontsize=11)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_dir / 'speedup_factors.png', dpi=150)
    print(f"Saved: {output_dir}/speedup_factors.png")
    plt.close()

    # ============================================
    # Panel-style figure (A/B) for a "Nature" layout
    # A: total time (all versions)
    # B: head-to-head comparison
    # ============================================
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Panel A: total time comparison
    axA = axes[0]
    axA.plot(wasp1['n_reads_m'], wasp1['total'], 'o-', color='gray',
             label='WASP1', linewidth=2, markersize=4, alpha=0.8)
    axA.plot(wasp2_py['n_reads_m'], wasp2_py['total'], 's-', color='blue',
             label='WASP2-Python', linewidth=2, markersize=4, alpha=0.8)
    axA.plot(rust_snp['n_reads_m'], rust_snp['total'], 'D-', color='red',
             label='WASP2-Rust (SNPs)', linewidth=2, markersize=4, alpha=0.8)
    if has_indel_data:
        axA.plot(rust_indel['n_reads_m'], rust_indel['total'], '^-', color='forestgreen',
                 label='WASP2-Rust (+INDELs)', linewidth=2, markersize=4, alpha=0.8)
    axA.set_xlabel('Reads (millions)', fontsize=12)
    axA.set_ylabel('Total time (seconds)', fontsize=12)
    axA.set_title('A  Total runtime', loc='left', fontsize=13, fontweight='bold')
    axA.legend(fontsize=10)
    axA.grid(alpha=0.3)

    # Panel B: Head-to-head comparison STAR+WASP vs WASP2-Rust (SNP-only and +INDELS)
    # Same-cluster benchmark (2025-12-04): fair apples-to-apples comparison
    # Indel overhead from 15M benchmark: 20% total, applied to estimate 56M indel time
    axB = axes[1]
    methods = ['STAR+WASP\n(C++)', 'WASP2\n(Rust, SNPs)', 'WASP2\n(Rust, +INDELs)']
    times = [532.61, 560.19, 672.23]  # 672 = 560 * 1.20 (20% indel overhead from 15M bench)
    colors = ['steelblue', 'crimson', 'forestgreen']
    bars = axB.bar(methods, times, color=colors, edgecolor='black', width=0.5)
    for bar, t in zip(bars, times):
        axB.text(bar.get_x() + bar.get_width() / 2, t + max(times) * 0.02,
                 f'{t:.0f}s', ha='center', va='bottom', fontsize=10, fontweight='bold')
    axB.set_ylabel('Time (seconds)', fontsize=12)
    axB.set_title('B  HG00731 (56M reads, 8 threads)', loc='left', fontsize=13, fontweight='bold')
    axB.set_ylim(0, max(times) * 1.25)
    axB.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_dir / 'figure_benchmark_v3.png', dpi=300)
    plt.savefig(output_dir / 'figure_benchmark_v3.pdf', dpi=300)
    print(f"Saved: {output_dir}/figure_benchmark_v3.png/.pdf")
    plt.close()

    # ============================================
    # Print summary statistics
    # ============================================
    print("\n" + "=" * 70)
    print("PERFORMANCE SUMMARY - WASP2-Rust (SNPs)")
    print("=" * 70)

    print(f"\nSpeedup vs WASP1:        {merged['speedup_vs_wasp1'].mean():.1f}x average")
    print(f"Speedup vs WASP2-Python: {merged['speedup_vs_python'].mean():.1f}x average")

    print("\n" + "-" * 70)
    print("Detailed comparison (first 20 common data points):")
    print("-" * 70)
    print(f"{'Reads':>8} | {'WASP1':>8} | {'Python':>8} | {'Rust':>8} | {'vs Py':>6}")
    print("-" * 70)

    for _, row in merged.head(20).iterrows():
        print(f"{int(row['n_reads']/1e6):>7}M | {row['total_wasp1']:>7.0f}s | {row['total_py']:>7.0f}s | "
              f"{row['total_rust']:>7.0f}s | {row['speedup_vs_python']:>5.1f}x")

if __name__ == '__main__':
    main()
