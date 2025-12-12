#!/usr/bin/env python3
"""
Plot Figure 1B (timing) and Figure 1C (SNV read retention) for WASP2 paper.
GM12878 ATAC-seq 159M reads benchmark.
"""

import json
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# IBM colorblind-safe palette
COLORS = {
    'original': '#648FFF',   # Blue
    'retained': '#DC267F',   # Magenta
    'wasp1': '#785EF0',      # Purple
    'wasp2py': '#FFB000',    # Gold
    'wasp2rust_snp': '#FE6100',  # Orange
    'wasp2rust_indel': '#DC267F',  # Magenta
}

# Results directory
RESULTS_BASE = Path("/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/atac_gm12878/results")

# Load benchmark results
def load_results():
    results = {}

    # WASP1
    wasp1_dir = list(RESULTS_BASE.glob("wasp1_snp_FIXED_*"))
    if wasp1_dir:
        with open(wasp1_dir[0] / "benchmark_results.json") as f:
            results['wasp1'] = json.load(f)

    # WASP2-Python Dev MT
    wasp2py_dir = list(RESULTS_BASE.glob("wasp2python_snp_DEV_MT_*"))
    if wasp2py_dir:
        with open(wasp2py_dir[0] / "benchmark_results.json") as f:
            results['wasp2py'] = json.load(f)

    # WASP2-Rust SNP
    wasp2rust_snp_dir = list(RESULTS_BASE.glob("wasp2rust_snp_fixed_*"))
    if wasp2rust_snp_dir:
        with open(wasp2rust_snp_dir[0] / "benchmark_results.json") as f:
            results['wasp2rust_snp'] = json.load(f)

    # WASP2-Rust Indel
    wasp2rust_indel_dir = list(RESULTS_BASE.glob("wasp2rust_indel_fixed_*"))
    if wasp2rust_indel_dir:
        with open(wasp2rust_indel_dir[0] / "benchmark_results.json") as f:
            results['wasp2rust_indel'] = json.load(f)

    return results

def plot_figure1b(results, output_path):
    """Plot Figure 1B: WASP-only timing comparison."""

    fig, ax = plt.subplots(figsize=(8, 5))

    # Data
    pipelines = ['WASP1', 'WASP2-Python\n(8 threads)', 'WASP2-Rust\n(SNP)', 'WASP2-Rust\n(SNP+Indel)']
    times = [
        results['wasp1']['wasp_only_s'],
        results['wasp2py']['wasp_only_s'],
        results['wasp2rust_snp']['wasp_only_s'],
        results['wasp2rust_indel']['wasp_only_s'],
    ]
    colors = ['#785EF0', '#FFB000', '#FE6100', '#DC267F']

    x = np.arange(len(pipelines))
    bars = ax.bar(x, times, color=colors, edgecolor='black', linewidth=1)

    # Add speedup labels
    wasp1_time = times[0]
    for i, (bar, t) in enumerate(zip(bars, times)):
        speedup = wasp1_time / t
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 100,
                f'{speedup:.1f}x', ha='center', va='bottom', fontsize=11, fontweight='bold')
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height()/2,
                f'{t:.0f}s', ha='center', va='center', fontsize=10, color='white', fontweight='bold')

    ax.set_ylabel('WASP-only Time (seconds)', fontsize=12)
    ax.set_xlabel('Pipeline', fontsize=12)
    ax.set_title('Figure 1B: WASP Pipeline Timing Comparison\n(GM12878 ATAC-seq, 159M reads, het-only SNPs)', fontsize=12)
    ax.set_xticks(x)
    ax.set_xticklabels(pipelines, fontsize=10)
    ax.set_ylim(0, max(times) * 1.15)

    # Add reference line at WASP1 time
    ax.axhline(y=wasp1_time, color='gray', linestyle='--', alpha=0.5, label='WASP1 baseline')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved Figure 1B to {output_path}")
    plt.close()

def plot_figure1c(results, output_path):
    """Plot Figure 1C: SNV-overlapping read retention."""

    fig, ax = plt.subplots(figsize=(10, 6))

    # Data
    pipelines = ['WASP1', 'WASP2-Python', 'WASP2-Rust\n(SNP)', 'WASP2-Rust\n(SNP+Indel)']

    pre_wasp = [
        results['wasp1']['snv_reads_pre'],
        results['wasp2py']['snv_reads_pre'],
        results['wasp2rust_snp']['snv_reads_pre'],
        results['wasp2rust_indel']['snv_reads_pre'],
    ]

    post_wasp = [
        results['wasp1']['snv_reads_post'],
        results['wasp2py']['snv_reads_post'],
        results['wasp2rust_snp']['snv_reads_post'],
        results['wasp2rust_indel']['snv_reads_post'],
    ]

    pass_rates = [post/pre * 100 for pre, post in zip(pre_wasp, post_wasp)]

    x = np.arange(len(pipelines))
    width = 0.35

    bars1 = ax.bar(x - width/2, np.array(pre_wasp)/1e6, width, label='Original (pre-WASP)',
                   color=COLORS['original'], edgecolor='black', linewidth=1)
    bars2 = ax.bar(x + width/2, np.array(post_wasp)/1e6, width, label='Retained (post-WASP)',
                   color=COLORS['retained'], edgecolor='black', linewidth=1)

    # Add pass rate labels
    for i, (bar, rate) in enumerate(zip(bars2, pass_rates)):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.05,
                f'{rate:.1f}%', ha='center', va='bottom', fontsize=10, fontweight='bold')

    ax.set_ylabel('SNV-overlapping Read Pairs (millions)', fontsize=12)
    ax.set_xlabel('Pipeline', fontsize=12)
    ax.set_title('Figure 1C: SNV-overlapping Read Retention\n(GM12878 ATAC-seq, 159M reads)', fontsize=12)
    ax.set_xticks(x)
    ax.set_xticklabels(pipelines, fontsize=10)
    ax.legend(loc='upper right', fontsize=10)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved Figure 1C to {output_path}")
    plt.close()

def plot_combined(results, output_path):
    """Plot combined Figure 1B and 1C side by side."""

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # --- Figure 1B: Timing ---
    pipelines_b = ['WASP1', 'WASP2-Py\n(8T)', 'WASP2-Rust\n(SNP)', 'WASP2-Rust\n(Indel)']
    times = [
        results['wasp1']['wasp_only_s'],
        results['wasp2py']['wasp_only_s'],
        results['wasp2rust_snp']['wasp_only_s'],
        results['wasp2rust_indel']['wasp_only_s'],
    ]
    colors_b = ['#785EF0', '#FFB000', '#FE6100', '#DC267F']

    x = np.arange(len(pipelines_b))
    bars = ax1.bar(x, times, color=colors_b, edgecolor='black', linewidth=1)

    wasp1_time = times[0]
    for i, (bar, t) in enumerate(zip(bars, times)):
        speedup = wasp1_time / t
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 100,
                f'{speedup:.1f}x', ha='center', va='bottom', fontsize=10, fontweight='bold')
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height()/2,
                f'{t:.0f}s', ha='center', va='center', fontsize=9, color='white', fontweight='bold')

    ax1.set_ylabel('WASP-only Time (seconds)', fontsize=11)
    ax1.set_title('B. Pipeline Timing', fontsize=12, fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(pipelines_b, fontsize=9)
    ax1.set_ylim(0, max(times) * 1.2)
    ax1.axhline(y=wasp1_time, color='gray', linestyle='--', alpha=0.5)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # --- Figure 1C: SNV retention ---
    pipelines_c = ['WASP1', 'WASP2-Py', 'WASP2-Rust\n(SNP)', 'WASP2-Rust\n(Indel)']

    pre_wasp = [
        results['wasp1']['snv_reads_pre'],
        results['wasp2py']['snv_reads_pre'],
        results['wasp2rust_snp']['snv_reads_pre'],
        results['wasp2rust_indel']['snv_reads_pre'],
    ]

    post_wasp = [
        results['wasp1']['snv_reads_post'],
        results['wasp2py']['snv_reads_post'],
        results['wasp2rust_snp']['snv_reads_post'],
        results['wasp2rust_indel']['snv_reads_post'],
    ]

    pass_rates = [post/pre * 100 for pre, post in zip(pre_wasp, post_wasp)]

    x = np.arange(len(pipelines_c))
    width = 0.35

    bars1 = ax2.bar(x - width/2, np.array(pre_wasp)/1e6, width, label='Original',
                   color=COLORS['original'], edgecolor='black', linewidth=1)
    bars2 = ax2.bar(x + width/2, np.array(post_wasp)/1e6, width, label='Retained',
                   color=COLORS['retained'], edgecolor='black', linewidth=1)

    for i, (bar, rate) in enumerate(zip(bars2, pass_rates)):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.05,
                f'{rate:.1f}%', ha='center', va='bottom', fontsize=9, fontweight='bold')

    ax2.set_ylabel('SNV-overlapping Reads (millions)', fontsize=11)
    ax2.set_title('C. SNV Read Retention', fontsize=12, fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels(pipelines_c, fontsize=9)
    ax2.legend(loc='upper right', fontsize=9)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    plt.suptitle('Figure 1: WASP2 Performance on GM12878 ATAC-seq (159M reads)',
                 fontsize=13, fontweight='bold', y=1.02)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved combined figure to {output_path}")
    plt.close()

def print_summary(results):
    """Print summary table."""
    print("\n" + "="*70)
    print("BENCHMARK SUMMARY: GM12878 ATAC-seq (159M reads)")
    print("="*70)

    print("\nFigure 1B - WASP-only Timing:")
    print("-"*50)
    print(f"{'Pipeline':<25} {'Time (s)':<12} {'Speedup':<10}")
    print("-"*50)

    wasp1_time = results['wasp1']['wasp_only_s']

    for name, key in [('WASP1', 'wasp1'), ('WASP2-Python (8T)', 'wasp2py'),
                       ('WASP2-Rust SNP', 'wasp2rust_snp'), ('WASP2-Rust Indel', 'wasp2rust_indel')]:
        t = results[key]['wasp_only_s']
        speedup = wasp1_time / t
        print(f"{name:<25} {t:>10.1f}s {speedup:>8.1f}x")

    print("\nFigure 1C - SNV Read Retention:")
    print("-"*70)
    print(f"{'Pipeline':<25} {'Pre-WASP':<12} {'Post-WASP':<12} {'Pass Rate':<10}")
    print("-"*70)

    for name, key in [('WASP1', 'wasp1'), ('WASP2-Python', 'wasp2py'),
                       ('WASP2-Rust SNP', 'wasp2rust_snp'), ('WASP2-Rust Indel', 'wasp2rust_indel')]:
        pre = results[key]['snv_reads_pre']
        post = results[key]['snv_reads_post']
        rate = post / pre * 100
        print(f"{name:<25} {pre:>10,} {post:>10,} {rate:>8.2f}%")

    print("="*70 + "\n")

def main():
    # Load results
    results = load_results()

    if not results:
        print("ERROR: No benchmark results found!")
        return

    # Print summary
    print_summary(results)

    # Output directory
    output_dir = RESULTS_BASE.parent / "figures"
    output_dir.mkdir(exist_ok=True)

    # Plot individual figures
    plot_figure1b(results, output_dir / "figure1b_timing.png")
    plot_figure1c(results, output_dir / "figure1c_snv_retention.png")

    # Plot combined figure
    plot_combined(results, output_dir / "figure1bc_combined.png")

    print(f"\nAll figures saved to: {output_dir}")

if __name__ == "__main__":
    main()
