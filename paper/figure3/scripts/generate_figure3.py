#!/usr/bin/env python3
"""
Figure 3: Statistical Analysis

Panel A: Statistical methods schematic
Panel B: QQ plots for different methods
Panel C: SNV-based volcano plot
Panel D: Gene imprinting analysis (SNVs + INDELs) - HAVE DATA
Panel E: Het/homo QTL stratification
"""
import sys
from pathlib import Path

# Add paper directory to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from config import (
    COLORS, PLOT_SETTINGS, IMPRINTING_RESULTS,
    get_plot_path, get_data_path, REPO_ROOT,
    GENE_SIGNIFICANCE_TSV
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


def panel_a_methods_schematic(ax):
    """Panel A: Statistical methods schematic."""
    ax.text(0.5, 0.85, 'Statistical Methods',
            transform=ax.transAxes, ha='center', va='center',
            fontsize=10, fontweight='bold')

    methods_text = """
Overdispersion Models:
  • Binomial (naive)
  • Beta-binomial
  • Quasi-binomial

Aggregation:
  • SNV-level tests
  • Feature-level (Fisher)
  • FDR correction
"""
    ax.text(0.5, 0.45, methods_text,
            transform=ax.transAxes, ha='center', va='center',
            fontsize=7, color='gray', family='monospace')

    ax.text(0.5, 0.05, '(Create detailed schematic)',
            transform=ax.transAxes, ha='center', va='center',
            fontsize=6, style='italic', color='gray')

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')
    ax.set_title('A', fontsize=11, fontweight='bold', loc='left', x=-0.05)


def panel_b_qq_plots(ax):
    """Panel B: QQ plots for different methods."""
    data_file = get_data_path(3, 'qq_plot_data.tsv')

    if data_file.exists():
        df = pd.read_csv(data_file, sep='\t')

        # Expected vs observed p-values
        for method in df['method'].unique():
            subset = df[df['method'] == method]
            expected = -np.log10(np.linspace(1/len(subset), 1, len(subset)))
            observed = -np.log10(np.sort(subset['pvalue']))
            ax.scatter(expected, observed, alpha=0.5, s=10, label=method)

        # Identity line
        max_val = max(ax.get_xlim()[1], ax.get_ylim()[1])
        ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.5)

        ax.set_xlabel('Expected -log10(p)')
        ax.set_ylabel('Observed -log10(p)')
        ax.legend(fontsize=6)
    else:
        ax.text(0.5, 0.5, 'QQ Plots by Method',
                transform=ax.transAxes, ha='center', va='center',
                fontsize=9, style='italic', color='gray')
        ax.text(0.5, 0.3, 'iPSCORE data needed',
                transform=ax.transAxes, ha='center', va='center',
                fontsize=7, color='gray')

    ax.set_title('B', fontsize=11, fontweight='bold', loc='left', x=-0.15)


def panel_c_volcano(ax):
    """Panel C: SNV-based volcano plot."""
    data_file = get_data_path(3, 'volcano_data.tsv')

    if data_file.exists():
        df = pd.read_csv(data_file, sep='\t')

        # Volcano plot
        colors = np.where(df['fdr'] < 0.1, COLORS['vermillion'], 'gray')
        ax.scatter(df['log2_ratio'], -np.log10(df['pvalue']),
                  c=colors, alpha=0.5, s=10)

        ax.axhline(-np.log10(0.05), color='red', linestyle='--', alpha=0.5, linewidth=0.5)
        ax.axvline(0, color='gray', linestyle='-', alpha=0.3, linewidth=0.5)

        ax.set_xlabel('log2(Alt/Ref)')
        ax.set_ylabel('-log10(p-value)')

        # Count significant
        n_sig = (df['fdr'] < 0.1).sum()
        ax.text(0.95, 0.95, f'n={n_sig} sig',
                transform=ax.transAxes, ha='right', va='top', fontsize=7)
    else:
        ax.text(0.5, 0.5, 'Volcano Plot\n(Allelic Imbalance)',
                transform=ax.transAxes, ha='center', va='center',
                fontsize=9, style='italic', color='gray')
        ax.text(0.5, 0.3, 'Data needed',
                transform=ax.transAxes, ha='center', va='center',
                fontsize=7, color='gray')

    ax.set_title('C', fontsize=11, fontweight='bold', loc='left', x=-0.15)


def panel_d_imprinting(ax):
    """Panel D: Gene imprinting analysis (SNVs + INDELs)."""

    # Load gene_significance.tsv which has both SNP-only and SNP+INDEL FDRs
    if GENE_SIGNIFICANCE_TSV.exists():
        df = pd.read_csv(GENE_SIGNIFICANCE_TSV, sep='\t')

        # Exclude XIST (X-inactivation, not imprinting)
        df = df[df['gene'] != 'XIST']

        # Count significant genes (FDR < 0.1)
        snp_sig = df[df['snp_fdr'] < 0.1]['gene'].tolist()
        all_sig = df[df['all_fdr'] < 0.1]['gene'].tolist()

        # Bar chart
        categories = ['SNP-only', 'SNP+INDEL']
        counts = [len(snp_sig), len(all_sig)]

        bars = ax.bar(categories, counts,
                     color=[COLORS['blue'], COLORS['green']],
                     edgecolor='black', linewidth=0.5)

        for bar, val in zip(bars, counts):
            ax.text(bar.get_x() + bar.get_width()/2, val + 0.2, str(val),
                    ha='center', va='bottom', fontsize=9, fontweight='bold')

        ax.set_ylabel('Significant Genes\n(FDR < 0.1)')
        ax.set_ylim(0, max(counts) * 1.4)

        # New genes
        new_genes = [g for g in all_sig if g not in snp_sig]
        improvement = counts[1] / counts[0] if counts[0] > 0 else float('inf')

        annotation = f'{improvement:.1f}× more\nwith INDELs'
        ax.text(0.95, 0.95, annotation,
                transform=ax.transAxes, ha='right', va='top', fontsize=7,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    else:
        # Fallback
        categories = ['SNP-only', 'SNP+INDEL']
        counts = [IMPRINTING_RESULTS['snp_only_significant'],
                  IMPRINTING_RESULTS['snp_indel_significant']]

        bars = ax.bar(categories, counts,
                     color=[COLORS['blue'], COLORS['green']],
                     edgecolor='black', linewidth=0.5)

        for bar, val in zip(bars, counts):
            ax.text(bar.get_x() + bar.get_width()/2, val + 0.1, str(val),
                    ha='center', va='bottom', fontsize=9, fontweight='bold')

        ax.set_ylabel('Significant Genes\n(FDR < 0.1)')
        ax.set_ylim(0, 7)

    ax.set_title('D', fontsize=11, fontweight='bold', loc='left', x=-0.15)


def panel_e_qtl_stratification(ax):
    """Panel E: Het/homo QTL stratification."""
    data_file = get_data_path(3, 'qtl_stratification.tsv')

    if data_file.exists():
        df = pd.read_csv(data_file, sep='\t')

        # Box plots by genotype
        het_data = df[df['genotype'] == 'het']['allelic_ratio']
        hom_data = df[df['genotype'] == 'hom']['allelic_ratio']

        bp = ax.boxplot([het_data, hom_data], labels=['Het', 'Hom'],
                       patch_artist=True)
        bp['boxes'][0].set_facecolor(COLORS['blue'])
        bp['boxes'][1].set_facecolor(COLORS['green'])

        ax.axhline(0.5, color='gray', linestyle='--', alpha=0.5)
        ax.set_ylabel('Allelic ratio')
        ax.set_xlabel('Lead SNP genotype')
    else:
        ax.text(0.5, 0.5, 'QTL Stratification\n(Het vs Hom)',
                transform=ax.transAxes, ha='center', va='center',
                fontsize=9, style='italic', color='gray')
        ax.text(0.5, 0.3, 'iPSCORE QTL data needed',
                transform=ax.transAxes, ha='center', va='center',
                fontsize=7, color='gray')

    ax.set_title('E', fontsize=11, fontweight='bold', loc='left', x=-0.15)


def generate_figure3():
    """Generate complete Figure 3."""

    setup_style()

    # 2 rows: top row has 3 panels, bottom row has 2
    fig = plt.figure(figsize=(11, 7))

    # Top row: A, B, C
    ax_a = fig.add_subplot(2, 3, 1)
    ax_b = fig.add_subplot(2, 3, 2)
    ax_c = fig.add_subplot(2, 3, 3)

    # Bottom row: D, E (centered)
    ax_d = fig.add_subplot(2, 3, 4)
    ax_e = fig.add_subplot(2, 3, 5)
    # Leave 2,3,6 empty for centering effect

    panel_a_methods_schematic(ax_a)
    panel_b_qq_plots(ax_b)
    panel_c_volcano(ax_c)
    panel_d_imprinting(ax_d)
    panel_e_qtl_stratification(ax_e)

    fig.suptitle('Figure 3: Statistical Analysis', fontsize=12, fontweight='bold', y=1.02)

    plt.tight_layout()

    # Save
    out_dir = get_plot_path(3, 'figure3').parent
    out_dir.mkdir(parents=True, exist_ok=True)

    plt.savefig(get_plot_path(3, 'figure3', 'png'), bbox_inches='tight', facecolor='white')
    plt.savefig(get_plot_path(3, 'figure3', 'pdf'), bbox_inches='tight', facecolor='white')

    print(f"Saved: {get_plot_path(3, 'figure3', 'png')}")
    print(f"Saved: {get_plot_path(3, 'figure3', 'pdf')}")

    return fig


if __name__ == '__main__':
    generate_figure3()
