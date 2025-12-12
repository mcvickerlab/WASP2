#!/usr/bin/env python3
"""
Figure 2: RNA-seq Read Filtering

Panel A: RNA-seq read filtering speed (HG00731, 56M reads, 8 threads)
         Comparing WASP1, STAR-WASP, and WASP2-Rust
Panel B: Ref/Alt count comparison across methods
Panel C: Original vs remapped BAM ref/alt comparison
"""
import sys
from pathlib import Path

# Add paper directory to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from config import (
    BENCHMARK_DATA, COLORS, TOOL_COLORS, PLOT_SETTINGS,
    get_plot_path, get_data_path
)


def setup_style():
    """Set up Nature Methods plotting style."""
    plt.rcParams.update({
        'font.family': PLOT_SETTINGS['font_family'],
        'font.sans-serif': PLOT_SETTINGS['font_sans_serif'],
        'font.size': PLOT_SETTINGS['font_size'],
        'axes.labelsize': PLOT_SETTINGS['axes_labelsize'],
        'axes.titlesize': PLOT_SETTINGS['axes_titlesize'],
        'axes.spines.top': False,
        'axes.spines.right': False,
        'figure.dpi': PLOT_SETTINGS['figure_dpi'],
        'savefig.dpi': PLOT_SETTINGS['savefig_dpi'],
    })


def panel_a_tool_comparison(ax):
    """Panel A: RNA-seq read filtering speed comparison."""

    # Get benchmark data
    rna = BENCHMARK_DATA['rnaseq_filtering']

    # RNA-seq data (HG00731, 56M reads, 8 threads)
    tools = ['WASP1', 'STAR-WASP', 'WASP2-Rust']
    times = [
        rna['wasp1']['total_s'],
        rna['star_wasp']['total_s'],
        rna['wasp2_rust']['total_s'],
    ]

    # Color scheme
    colors = [COLORS['black'], COLORS['sky_blue'], COLORS['green']]

    x = np.arange(len(tools))
    width = 0.6

    # Plot bars
    bars = ax.bar(x, times, width, color=colors, edgecolor='black', linewidth=0.5)

    # Add time labels on bars
    def format_time(val):
        if val >= 3600:
            return f'{val/3600:.1f}h'
        elif val >= 60:
            return f'{val/60:.0f}m'
        else:
            return f'{val:.0f}s'

    for bar, val in zip(bars, times):
        ax.text(bar.get_x() + bar.get_width()/2, val * 1.1, format_time(val),
                ha='center', va='bottom', fontsize=7)

    # X-axis labels
    ax.set_xticks(x)
    ax.set_xticklabels(tools, fontsize=8)

    # Add dataset label
    ax.text(0.5, -0.15, 'RNA-seq (HG00731, 56M reads, 8 threads)',
            transform=ax.transAxes, ha='center', va='top', fontsize=7, style='italic')

    ax.set_ylabel('Time (seconds)')
    ax.set_yscale('log')
    ax.set_ylim(100, 3000)

    # Speedup annotation - WASP2-Rust vs WASP1
    speedup = rna['wasp1']['total_s'] / rna['wasp2_rust']['total_s']
    ax.annotate(f'{speedup:.1f}×', xy=(2, rna['wasp2_rust']['total_s']),
                xytext=(1.5, rna['wasp1']['total_s'] * 0.6),
                fontsize=8, fontweight='bold', color=COLORS['green'],
                ha='center', va='bottom',
                arrowprops=dict(arrowstyle='->', color=COLORS['green'], lw=1.5))

    ax.set_title('A', fontsize=11, fontweight='bold', loc='left', x=-0.15)


def panel_b_count_comparison(ax):
    """Panel B: Ref/Alt count comparison across methods."""

    data_file = get_data_path(2, 'count_comparison.tsv')

    if data_file.exists():
        df = pd.read_csv(data_file, sep='\t')

        # Scatter: WASP2 vs GATK
        ax.scatter(df['wasp2_ref'], df['gatk_ref'], alpha=0.3, s=5,
                  color=COLORS['blue'], label='Ref counts')
        ax.scatter(df['wasp2_alt'], df['gatk_alt'], alpha=0.3, s=5,
                  color=COLORS['vermillion'], label='Alt counts')

        # Identity line
        max_val = max(df['wasp2_ref'].max(), df['gatk_ref'].max())
        ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.5, linewidth=1)

        ax.set_xlabel('WASP2-Rust counts')
        ax.set_ylabel('GATK counts')
        ax.legend(fontsize=7)

        # Correlation
        from scipy.stats import pearsonr
        r, _ = pearsonr(df['wasp2_ref'], df['gatk_ref'])
        ax.text(0.05, 0.95, f'r = {r:.3f}', transform=ax.transAxes,
                fontsize=8, va='top')
    else:
        # Placeholder
        ax.text(0.5, 0.5, 'Ref/Alt Count Comparison\n(WASP2 vs GATK vs phASER)',
                transform=ax.transAxes, ha='center', va='center',
                fontsize=9, style='italic', color='gray')
        ax.text(0.5, 0.25, 'Data needed:\nRun compare_counts_3tools.py',
                transform=ax.transAxes, ha='center', va='center',
                fontsize=7, color='gray')

    ax.set_title('B', fontsize=11, fontweight='bold', loc='left', x=-0.15)


def panel_c_original_vs_remapped(ax):
    """Panel C: Original vs remapped BAM ref/alt comparison."""

    data_file = get_data_path(2, 'original_vs_remapped.tsv')

    if data_file.exists():
        df = pd.read_csv(data_file, sep='\t')

        # Calculate ref bias (ref / (ref + alt))
        df['orig_bias'] = df['orig_ref'] / (df['orig_ref'] + df['orig_alt'])
        df['remap_bias'] = df['remap_ref'] / (df['remap_ref'] + df['remap_alt'])

        # Scatter plot
        ax.scatter(df['orig_bias'], df['remap_bias'], alpha=0.3, s=5, color=COLORS['blue'])

        # Identity line
        ax.plot([0, 1], [0, 1], 'k--', alpha=0.5, linewidth=1)
        ax.axhline(0.5, color='gray', linestyle=':', alpha=0.5)
        ax.axvline(0.5, color='gray', linestyle=':', alpha=0.5)

        ax.set_xlabel('Original BAM (ref ratio)')
        ax.set_ylabel('Remapped BAM (ref ratio)')
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)

        # Annotation about bias reduction
        orig_mean = abs(df['orig_bias'].mean() - 0.5)
        remap_mean = abs(df['remap_bias'].mean() - 0.5)
        reduction = (orig_mean - remap_mean) / orig_mean * 100 if orig_mean > 0 else 0
        ax.text(0.95, 0.05, f'{reduction:.0f}% bias\nreduction',
                transform=ax.transAxes, ha='right', va='bottom', fontsize=7,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    else:
        # Placeholder
        ax.text(0.5, 0.5, 'Original vs Remapped BAM\nRef/Alt Comparison',
                transform=ax.transAxes, ha='center', va='center',
                fontsize=9, style='italic', color='gray')
        ax.text(0.5, 0.25, 'Data needed:\nRun compare_original_vs_remapped.py',
                transform=ax.transAxes, ha='center', va='center',
                fontsize=7, color='gray')

    ax.set_title('C', fontsize=11, fontweight='bold', loc='left', x=-0.15)


def generate_figure2():
    """Generate complete Figure 2."""

    setup_style()

    fig, axes = plt.subplots(1, 3, figsize=(11, 3.5))

    panel_a_tool_comparison(axes[0])
    panel_b_count_comparison(axes[1])
    panel_c_original_vs_remapped(axes[2])

    fig.suptitle('Figure 2: RNA-seq Read Filtering', fontsize=12, fontweight='bold', y=1.02)

    plt.tight_layout()

    # Save
    out_dir = get_plot_path(2, 'figure2').parent
    out_dir.mkdir(parents=True, exist_ok=True)

    plt.savefig(get_plot_path(2, 'figure2', 'png'), bbox_inches='tight', facecolor='white')
    plt.savefig(get_plot_path(2, 'figure2', 'pdf'), bbox_inches='tight', facecolor='white')

    print(f"Saved: {get_plot_path(2, 'figure2', 'png')}")
    print(f"Saved: {get_plot_path(2, 'figure2', 'pdf')}")

    return fig


if __name__ == '__main__':
    generate_figure2()
