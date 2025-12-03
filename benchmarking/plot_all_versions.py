#!/usr/bin/env python3
"""
Plot WASP performance comparison: WASP1 vs WASP2 Python vs Rust v1.2.0 vs Rust v1.3.0
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# Data sources
WASP1_FILE = '/iblm/netapp/home/aho/projects/wasp/testing/performance/outputs/test_logs_v1/wasp1_perf_logs_snp2h5_7881649-Copy1.txt'
WASP2_PYTHON_FILE = '/iblm/netapp/home/aho/projects/wasp/testing/performance/outputs/test_logs_v1/wasp2_perf_logs_7908428-Copy1.txt'
RUST_V120_FILE = 'results/wasp2_rust_indels_8498628.tsv'
RUST_V130_FILE = 'results/wasp2_rust_scale_8512868.tsv'

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
    rust_v120 = load_rust(RUST_V120_FILE)
    rust_v130 = load_rust(RUST_V130_FILE)

    # Convert to millions
    wasp1['n_reads_m'] = wasp1['n_reads'] / 1e6
    wasp2_py['n_reads_m'] = wasp2_py['n_reads'] / 1e6
    rust_v120['n_reads_m'] = rust_v120['n_reads'] / 1e6
    rust_v130['n_reads_m'] = rust_v130['n_reads'] / 1e6

    # ============================================
    # Plot 1: Full comparison (all 4 versions)
    # ============================================
    fig, ax = plt.subplots(figsize=(12, 7))

    ax.plot(wasp1['n_reads_m'], wasp1['total'], 'o-', color='gray',
            label='WASP1 (original)', linewidth=2, markersize=5, alpha=0.8)
    ax.plot(wasp2_py['n_reads_m'], wasp2_py['total'], 's-', color='blue',
            label='WASP2 Python', linewidth=2, markersize=5, alpha=0.8)
    ax.plot(rust_v120['n_reads_m'], rust_v120['total'], '^-', color='orange',
            label='WASP2 Rust v1.2.0', linewidth=2, markersize=5, alpha=0.8)
    ax.plot(rust_v130['n_reads_m'], rust_v130['total'], 'D-', color='red',
            label='WASP2 Rust v1.3.0 (parallel streaming)', linewidth=2, markersize=5, alpha=0.8)

    ax.set_xlabel('Number of Reads (millions)', fontsize=14)
    ax.set_ylabel('Total Time (seconds)', fontsize=14)
    ax.set_title('WASP Performance Evolution\nWASP1 → WASP2 Python → Rust v1.2.0 → Rust v1.3.0', fontsize=16)
    ax.legend(loc='upper left', fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, max(rust_v130['n_reads_m'].max(), 150) + 5)

    plt.tight_layout()
    plt.savefig(output_dir / 'wasp_all_versions_comparison.png', dpi=150)
    print(f"Saved: {output_dir}/wasp_all_versions_comparison.png")
    plt.close()

    # ============================================
    # Plot 2: Zoomed in on Rust versions (up to current data)
    # ============================================
    fig, ax = plt.subplots(figsize=(12, 7))

    # Only show data we have
    max_reads = min(rust_v120['n_reads_m'].max(), rust_v130['n_reads_m'].max())

    v120_subset = rust_v120[rust_v120['n_reads_m'] <= max_reads]
    v130_subset = rust_v130[rust_v130['n_reads_m'] <= max_reads]
    py_subset = wasp2_py[wasp2_py['n_reads_m'] <= max_reads]

    ax.plot(py_subset['n_reads_m'], py_subset['total'], 's-', color='blue',
            label='WASP2 Python', linewidth=2, markersize=6, alpha=0.8)
    ax.plot(v120_subset['n_reads_m'], v120_subset['total'], '^-', color='orange',
            label='Rust v1.2.0', linewidth=2, markersize=6, alpha=0.8)
    ax.plot(v130_subset['n_reads_m'], v130_subset['total'], 'D-', color='red',
            label='Rust v1.3.0 (parallel streaming)', linewidth=2, markersize=6, alpha=0.8)

    ax.set_xlabel('Number of Reads (millions)', fontsize=14)
    ax.set_ylabel('Total Time (seconds)', fontsize=14)
    ax.set_title('WASP2 Rust Optimization Progress', fontsize=16)
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
    merged = rust_v130.merge(rust_v120, on='n_reads', suffixes=('_v130', '_v120'))
    merged = merged.merge(wasp2_py, on='n_reads')
    merged = merged.merge(wasp1, on='n_reads', suffixes=('_py', '_wasp1'))

    merged['speedup_vs_python'] = merged['total_py'] / merged['total_v130']
    merged['speedup_vs_v120'] = merged['total_v120'] / merged['total_v130']
    merged['speedup_vs_wasp1'] = merged['total_wasp1'] / merged['total_v130']
    merged['n_reads_m'] = merged['n_reads'] / 1e6

    ax.plot(merged['n_reads_m'], merged['speedup_vs_wasp1'], 'o-', color='gray',
            label=f'vs WASP1 (avg {merged["speedup_vs_wasp1"].mean():.1f}x)', linewidth=2, markersize=5)
    ax.plot(merged['n_reads_m'], merged['speedup_vs_python'], 's-', color='blue',
            label=f'vs WASP2 Python (avg {merged["speedup_vs_python"].mean():.1f}x)', linewidth=2, markersize=5)
    ax.plot(merged['n_reads_m'], merged['speedup_vs_v120'], '^-', color='orange',
            label=f'vs Rust v1.2.0 (avg {merged["speedup_vs_v120"].mean():.1f}x)', linewidth=2, markersize=5)

    ax.axhline(y=1, color='black', linestyle='--', alpha=0.5)
    ax.set_xlabel('Number of Reads (millions)', fontsize=14)
    ax.set_ylabel('Speedup Factor (X times faster)', fontsize=14)
    ax.set_title('Rust v1.3.0 Speedup vs Other Versions', fontsize=16)
    ax.legend(loc='upper right', fontsize=11)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_dir / 'speedup_factors.png', dpi=150)
    print(f"Saved: {output_dir}/speedup_factors.png")
    plt.close()

    # ============================================
    # Print summary statistics
    # ============================================
    print("\n" + "=" * 70)
    print("PERFORMANCE SUMMARY - Rust v1.3.0 (Parallel Streaming FASTQ Writes)")
    print("=" * 70)

    print(f"\nSpeedup vs WASP1:        {merged['speedup_vs_wasp1'].mean():.1f}x average")
    print(f"Speedup vs WASP2 Python: {merged['speedup_vs_python'].mean():.1f}x average")
    print(f"Speedup vs Rust v1.2.0:  {merged['speedup_vs_v120'].mean():.1f}x average")

    print("\n" + "-" * 70)
    print("Detailed comparison (first 20 common data points):")
    print("-" * 70)
    print(f"{'Reads':>8} | {'WASP1':>8} | {'Python':>8} | {'v1.2.0':>8} | {'v1.3.0':>8} | {'vs Py':>6} | {'vs 1.2':>6}")
    print("-" * 70)

    for _, row in merged.head(20).iterrows():
        print(f"{int(row['n_reads']/1e6):>7}M | {row['total_wasp1']:>7.0f}s | {row['total_py']:>7.0f}s | "
              f"{row['total_v120']:>7.0f}s | {row['total_v130']:>7.0f}s | "
              f"{row['speedup_vs_python']:>5.1f}x | {row['speedup_vs_v120']:>5.1f}x")

if __name__ == '__main__':
    main()
