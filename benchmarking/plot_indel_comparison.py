#!/usr/bin/env python3
"""
Plot WASP2-Rust benchmark: SNP-only vs WITH INDELS
15M reads comparison using unified pipeline
"""
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# Benchmark results from 2025-12-04 (unified pipeline, fixed)
# 15M reads, 8 threads

RESULTS = {
    'SNP-only': {
        'make_reads': 10.39,
        'bwa_remap': 28.66,
        'filter': 3.69,
        'wasp_only': 14.08,  # make_reads + filter
        'total': 42.78,
    },
    'WITH INDELS': {
        'make_reads': 11.13,
        'bwa_remap': 36.09,
        'filter': 4.25,
        'wasp_only': 15.37,  # make_reads + filter
        'total': 51.50,
    },
}

def create_comparison_figure(output_dir: str = 'plots'):
    """Create bar chart comparing SNP-only vs INDELS."""
    Path(output_dir).mkdir(exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Panel A: Total time comparison
    ax = axes[0]
    methods = ['SNP-only', 'WITH\nINDELS']
    totals = [RESULTS['SNP-only']['total'], RESULTS['WITH INDELS']['total']]
    colors = ['steelblue', 'forestgreen']

    bars = ax.bar(methods, totals, color=colors, edgecolor='black', width=0.6)

    for bar, t in zip(bars, totals):
        ax.text(bar.get_x() + bar.get_width()/2, t + 2,
                f'{t:.1f}s', ha='center', va='bottom', fontsize=14, fontweight='bold')

    # Add overhead annotation
    overhead = (RESULTS['WITH INDELS']['total'] / RESULTS['SNP-only']['total'] - 1) * 100
    ax.annotate(f'+{overhead:.0f}% overhead',
                xy=(1, RESULTS['WITH INDELS']['total']),
                xytext=(1.3, RESULTS['WITH INDELS']['total'] + 8),
                fontsize=12, color='forestgreen', fontweight='bold',
                arrowprops=dict(arrowstyle='->', color='forestgreen', lw=1.5))

    ax.set_ylabel('Total Time (seconds)', fontsize=14)
    ax.set_title('A  WASP2-Rust Total Runtime (15M reads)', loc='left', fontsize=14, fontweight='bold')
    ax.set_ylim(0, max(totals) * 1.35)
    ax.grid(axis='y', alpha=0.3)

    # Panel B: Stage breakdown comparison
    ax2 = axes[1]

    stages = ['Make-reads', 'BWA remap', 'Filter']
    snp_times = [RESULTS['SNP-only']['make_reads'],
                 RESULTS['SNP-only']['bwa_remap'],
                 RESULTS['SNP-only']['filter']]
    indel_times = [RESULTS['WITH INDELS']['make_reads'],
                   RESULTS['WITH INDELS']['bwa_remap'],
                   RESULTS['WITH INDELS']['filter']]

    x = np.arange(len(stages))
    width = 0.35

    bars1 = ax2.bar(x - width/2, snp_times, width, label='SNP-only',
                    color='steelblue', edgecolor='black')
    bars2 = ax2.bar(x + width/2, indel_times, width, label='WITH INDELS',
                    color='forestgreen', edgecolor='black')

    # Add value labels
    for bars, times in [(bars1, snp_times), (bars2, indel_times)]:
        for bar, t in zip(bars, times):
            ax2.text(bar.get_x() + bar.get_width()/2, t + 1,
                     f'{t:.1f}s', ha='center', va='bottom', fontsize=10)

    ax2.set_xticks(x)
    ax2.set_xticklabels(stages, fontsize=12)
    ax2.set_ylabel('Time (seconds)', fontsize=14)
    ax2.set_title('B  Stage Breakdown', loc='left', fontsize=14, fontweight='bold')
    ax2.legend(loc='upper right', fontsize=11)
    ax2.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(f'{output_dir}/indel_comparison.png', dpi=150)
    plt.savefig(f'{output_dir}/indel_comparison.pdf', dpi=150)
    print(f"Saved: {output_dir}/indel_comparison.png/.pdf")
    plt.close()

    # Print summary
    print("\n" + "=" * 60)
    print("WASP2-Rust Indel Benchmark Summary (15M reads, 8 threads)")
    print("=" * 60)
    print(f"\n{'Stage':<15} {'SNP-only':>10} {'INDELS':>10} {'Overhead':>12}")
    print("-" * 60)
    print(f"{'Make-reads':<15} {RESULTS['SNP-only']['make_reads']:>9.1f}s {RESULTS['WITH INDELS']['make_reads']:>9.1f}s {(RESULTS['WITH INDELS']['make_reads']/RESULTS['SNP-only']['make_reads']-1)*100:>+10.0f}%")
    print(f"{'BWA remap':<15} {RESULTS['SNP-only']['bwa_remap']:>9.1f}s {RESULTS['WITH INDELS']['bwa_remap']:>9.1f}s {(RESULTS['WITH INDELS']['bwa_remap']/RESULTS['SNP-only']['bwa_remap']-1)*100:>+10.0f}%")
    print(f"{'Filter':<15} {RESULTS['SNP-only']['filter']:>9.1f}s {RESULTS['WITH INDELS']['filter']:>9.1f}s {(RESULTS['WITH INDELS']['filter']/RESULTS['SNP-only']['filter']-1)*100:>+10.0f}%")
    print("-" * 60)
    print(f"{'WASP-only':<15} {RESULTS['SNP-only']['wasp_only']:>9.1f}s {RESULTS['WITH INDELS']['wasp_only']:>9.1f}s {(RESULTS['WITH INDELS']['wasp_only']/RESULTS['SNP-only']['wasp_only']-1)*100:>+10.0f}%")
    print(f"{'TOTAL':<15} {RESULTS['SNP-only']['total']:>9.1f}s {RESULTS['WITH INDELS']['total']:>9.1f}s {(RESULTS['WITH INDELS']['total']/RESULTS['SNP-only']['total']-1)*100:>+10.0f}%")
    print("=" * 60)
    print("\nConclusion: INDEL support adds only ~9% overhead to WASP processing time")
    print("            (20% total when including BWA remapping)")


if __name__ == '__main__':
    create_comparison_figure()
