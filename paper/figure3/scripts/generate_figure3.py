#!/usr/bin/env python3
"""
Figure 3: Statistical Analysis

Panel A: Statistical methods schematic (binomial → FDR pipeline)
Panel B: QQ plots comparing 6 representative methods (of 14+ tested)
         - Shows lambda values for each method
         - Best method: BB-cov (exp15_b10, λ=0.28)
Panel C: SNV-based volcano plot showing allelic imbalance
Panel D: Gene imprinting analysis (SNV vs SNV+INDEL comparison)
Panel E: QTL replication rates (RNA-seq 45.4%, ATAC-seq 41.9%)
         - 137 RNA samples, 137 ATAC samples from iPSCORE CVPC
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
    get_plot_path, get_data_path, REPO_ROOT
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
    """Panel A: Statistical methods flowchart.

    K-Dense scientific-schematics compliant:
    - Clean layout with proper spacing
    - Colorblind-safe colors (Bang Wong palette)
    - Clear labels and arrows
    - Professional rounded boxes with subtle fill colors
    """
    from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
    from matplotlib.colors import to_rgba

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    # Colors from Bang Wong palette - lighter fills for boxes
    blue = COLORS['blue']
    green = COLORS['green']
    orange = COLORS['orange']
    sky_blue = COLORS['sky_blue']

    # Helper to create light fill color
    def light_fill(color, alpha=0.15):
        rgba = to_rgba(color)
        return (rgba[0], rgba[1], rgba[2], alpha)

    # Professional box dimensions
    box_width = 0.55
    box_height = 0.12
    center_x = 0.5

    # Vertical positions with better spacing
    y_positions = [0.88, 0.66, 0.44, 0.22]
    arrow_offset = 0.065  # Gap for arrows

    # Box styling - more professional with rounded corners and subtle fills
    def draw_box(ax, x, y, w, h, text, edge_color, fill_alpha=0.12):
        box = FancyBboxPatch(
            (x - w/2, y - h/2), w, h,
            boxstyle="round,pad=0.02,rounding_size=0.02",
            facecolor=light_fill(edge_color, fill_alpha),
            edgecolor=edge_color,
            linewidth=1.5,
            transform=ax.transData
        )
        ax.add_patch(box)
        ax.text(x, y, text, ha='center', va='center',
                fontsize=7, fontweight='medium', linespacing=1.2)

    # Draw boxes with improved styling
    draw_box(ax, center_x, y_positions[0], box_width, box_height,
             'Allele Counts\n(ref, alt per site)', blue)

    draw_box(ax, center_x, y_positions[1], box_width, box_height,
             'Binomial Test\nH\u2080: ratio = 0.5', green)

    draw_box(ax, center_x, y_positions[2], box_width, box_height,
             'FDR Correction\n(Benjamini-Hochberg)', orange)

    draw_box(ax, center_x, y_positions[3], box_width, box_height,
             'Significant Sites\n(FDR < 0.1)', sky_blue)

    # Draw clean arrows between boxes
    arrow_style = dict(arrowstyle='-|>', color='#444444', lw=1.2,
                       mutation_scale=10, shrinkA=0, shrinkB=0)

    for i in range(len(y_positions) - 1):
        y_start = y_positions[i] - box_height/2 - 0.01
        y_end = y_positions[i+1] + box_height/2 + 0.01
        ax.annotate('', xy=(center_x, y_end), xytext=(center_x, y_start),
                    arrowprops=arrow_style)

    # Side annotations - positioned better
    ax.text(0.88, y_positions[0], 'SNVs\n+\nINDELs', ha='center', va='center',
            fontsize=6, color='#666666', style='italic')

    # Connecting line to annotation
    ax.plot([center_x + box_width/2 + 0.02, 0.82], [y_positions[0], y_positions[0]],
            color='#999999', linestyle=':', linewidth=0.8)

    ax.text(0.12, y_positions[2], 'Per-gene\naggregation', ha='center', va='center',
            fontsize=6, color='#666666', style='italic')

    # Connecting line to annotation
    ax.plot([center_x - box_width/2 - 0.02, 0.18], [y_positions[2], y_positions[2]],
            color='#999999', linestyle=':', linewidth=0.8)

    ax.set_title('A', fontsize=9, fontweight='bold', loc='left')


def panel_b_qq_plots(ax):
    """Panel B: QQ plots for different methods.

    Shows comparison across 6 representative statistical approaches selected from
    14+ comprehensive experiments comparing:
    - Naive binomial vs beta-binomial models
    - Global vs coverage-binned dispersion estimation
    - Different FDR correction methods (BH, discrete, mid-p)
    """
    from scipy import stats

    data_file = get_data_path(3, 'qq_plot_data.tsv')

    # Human-readable method names with descriptions
    # These 6 methods represent key points in the method space from 14+ experiments
    METHOD_LABELS = {
        'exp12': 'Binomial',           # Naive binomial (baseline)
        'exp13_b10': 'BB-global',       # Beta-binomial, global dispersion
        'exp15_b10': 'BB-cov',          # Beta-binomial, coverage-binned (BEST)
        'exp16_b10': 'BB-cov-disc',     # Coverage-binned, discrete FDR
        'exp17': 'BB-pooled',           # Cross-sample pooled
        'exp20': 'BB-cross',            # Cross-sample aggregation
    }

    # Colors for each method - use colorblind-safe palette
    METHOD_COLORS = {
        'exp12': COLORS['blue'],
        'exp13_b10': COLORS['green'],
        'exp15_b10': COLORS['orange'],      # Highlight best method
        'exp16_b10': COLORS['vermillion'],
        'exp17': COLORS['purple'],
        'exp20': COLORS['sky_blue'],
    }

    # Track lambdas per method for legend
    method_lambdas = {}

    if data_file.exists():
        df = pd.read_csv(data_file, sep='\t')

        # Expected vs observed p-values
        methods = df['method'].unique()

        for method in methods:
            subset = df[df['method'] == method]
            expected = -np.log10(np.linspace(1/len(subset), 1, len(subset)))
            observed = -np.log10(np.sort(subset['pvalue']))

            # Calculate genomic inflation factor (lambda_gc) for this method
            pvals = subset['pvalue'].values
            pvals = pvals[pvals > 0]  # Exclude zero p-values
            lambda_gc = None
            if len(pvals) > 0:
                chi2_obs = stats.chi2.ppf(1 - pvals, df=1)
                lambda_gc = np.median(chi2_obs) / stats.chi2.ppf(0.5, df=1)
                method_lambdas[method] = lambda_gc

            # Create label with lambda value
            base_label = METHOD_LABELS.get(method, method)
            if lambda_gc is not None:
                label = f'{base_label} (\u03bb={lambda_gc:.2f})'
            else:
                label = base_label

            color = METHOD_COLORS.get(method, 'gray')

            ax.scatter(expected, observed, alpha=0.6, s=8, label=label,
                      color=color, edgecolors='none')

        # Identity line
        max_val = max(ax.get_xlim()[1], ax.get_ylim()[1])
        ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.5, linewidth=1)

        ax.set_xlabel('Expected -log\u2081\u2080(p)')
        ax.set_ylabel('Observed -log\u2081\u2080(p)')

        # Legend inside plot area - compact styling with lambdas
        ax.legend(fontsize=5, frameon=True, loc='lower right',
                  handletextpad=0.2, framealpha=0.95, edgecolor='none',
                  labelspacing=0.25, markerscale=0.7)

        # Add note about comprehensive comparison
        ax.text(0.03, 0.97, '6 of 14+ methods',
                transform=ax.transAxes, ha='left', va='top',
                fontsize=5, color='gray', style='italic')
    else:
        ax.text(0.5, 0.5, 'QQ Plots by Method',
                transform=ax.transAxes, ha='center', va='center',
                fontsize=9, style='italic', color='gray')
        ax.text(0.5, 0.3, 'iPSCORE data needed',
                transform=ax.transAxes, ha='center', va='center',
                fontsize=7, color='gray')

    # Clean panel label only
    ax.set_title('B', fontsize=9, fontweight='bold', loc='left')


def panel_c_volcano(ax):
    """Panel C: SNV-based volcano plot showing allelic imbalance."""
    data_file = get_data_path(3, 'volcano_data.tsv')

    # Track n_sig for annotation
    n_sig = None
    n_total = None

    if data_file.exists():
        df = pd.read_csv(data_file, sep='\t')
        n_total = len(df)

        # Separate significant and non-significant for layered plotting
        sig_mask = df['fdr'] < 0.1
        df_nonsig = df[~sig_mask]
        df_sig = df[sig_mask]

        # Plot non-significant first (background)
        ax.scatter(df_nonsig['log2_ratio'], -np.log10(df_nonsig['pvalue']),
                  c='#CCCCCC', alpha=0.4, s=6, edgecolors='none', rasterized=True)

        # Plot significant on top (foreground)
        ax.scatter(df_sig['log2_ratio'], -np.log10(df_sig['pvalue']),
                  c=COLORS['vermillion'], alpha=0.7, s=8, edgecolors='none', rasterized=True)

        # Reference line at ratio=1 (log2=0)
        ax.axvline(0, color='#666666', linestyle='-', alpha=0.4, linewidth=0.8)

        ax.set_xlabel('log\u2082(Alt/Ref)')
        ax.set_ylabel('-log\u2081\u2080(p-value)')

        # Count significant
        n_sig = sig_mask.sum()
    else:
        ax.text(0.5, 0.5, 'Volcano Plot\n(Allelic Imbalance)',
                transform=ax.transAxes, ha='center', va='center',
                fontsize=9, style='italic', color='gray')
        ax.text(0.5, 0.3, 'Data needed',
                transform=ax.transAxes, ha='center', va='center',
                fontsize=7, color='gray')

    # Clean panel label only - stats as small annotation inside plot
    ax.set_title('C', fontsize=9, fontweight='bold', loc='left')
    if n_sig is not None:
        ax.text(0.95, 0.95, f'n={n_sig:,} sig\n({n_sig/n_total*100:.1f}%)',
                transform=ax.transAxes, ha='right', va='top', fontsize=6, color='gray')


def panel_d_imprinting(ax):
    """Panel D: Gene imprinting analysis (SNVs + INDELs).

    Shows improvement in detection of imprinted genes when including INDEL variants.
    """

    data_file = get_data_path(3, "gene_significance.tsv")

    # Load gene_significance.tsv which has both SNP-only and SNP+INDEL FDRs
    if data_file.exists():
        df = pd.read_csv(data_file, sep="\t")

        # Restrict to the intended positive-control set used in the imprinting validation:
        # "Imprinted" + "Predicted" (Aaron's list), and genes with coverage.
        df = df[df["status"].isin(["Imprinted", "Predicted"])].copy()
        df = df[df["all_total"] > 0].copy()
        # Exclude XIST (X-inactivation, not canonical imprinting positive-control)
        df = df[df["gene"] != "XIST"].copy()

        # Count significant genes (FDR < 0.1)
        snp_sig = df[(df["snp_total"] > 0) & (df["snp_fdr"] < 0.1)]["gene"].tolist()
        all_sig = df[df["all_fdr"] < 0.1]["gene"].tolist()

        # Bar chart with improved styling
        categories = ['SNV-only', 'SNV+INDEL']
        counts = [len(snp_sig), len(all_sig)]

        bars = ax.bar(categories, counts,
                     color=[COLORS['blue'], COLORS['green']],
                     edgecolor='white', linewidth=1.5, width=0.6)

        # Value labels on top of bars
        for bar, val in zip(bars, counts):
            ax.text(bar.get_x() + bar.get_width()/2, val + max(counts) * 0.03, str(val),
                    ha='center', va='bottom', fontsize=8, fontweight='bold')

        ax.set_ylabel('Significant Imprinted\nGenes (FDR < 0.1)')
        ax.set_ylim(0, max(counts) * 1.35)

        # New genes detected only with INDELs included
        new_genes = [g for g in all_sig if g not in snp_sig]
        new_genes_count = len(new_genes)

        # Show net gain as annotation above the bars
        if new_genes_count > 0:
            # Place +N annotation above the second bar with a small bracket/arrow indicator
            gain_y = counts[1] + max(counts) * 0.15
            ax.annotate(f'+{new_genes_count} new',
                       xy=(0.5, gain_y), xycoords='data',
                       fontsize=7, ha='center', va='bottom',
                       color=COLORS['green'], fontweight='bold')
            # Draw a small bracket/indicator line
            ax.plot([0, 1], [counts[1] + max(counts)*0.08, counts[1] + max(counts)*0.08],
                   color=COLORS['green'], linewidth=1, alpha=0.6)
    else:
        ax.text(
            0.5,
            0.55,
            "Imprinting data missing\n(run generate_panel_d_imprinting_data.py)",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=8,
            style="italic",
            color="gray",
        )
        ax.text(
            0.5,
            0.35,
            str(get_data_path(3, "gene_significance.tsv")),
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=6,
            color="gray",
        )
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")

    ax.set_title('D', fontsize=9, fontweight='bold', loc='left')


def panel_e_qtl_replication(ax):
    """Panel E: QTL Replication Analysis - RNA vs ATAC comparison.

    Shows the percentage of QTL variants that replicate in allelic imbalance
    analysis, comparing RNA-seq (eQTL) and ATAC-seq (caQTL) from iPSCORE CVPC.

    Key results:
    - RNA-seq: 45.4% eQTL replication (109/240 variants)
    - ATAC-seq: 41.9% caQTL replication (6,698/15,975 variants)
    """
    # QTL replication data from comprehensive analysis
    # Source: cvpc/results/analysis/peak_ai_qtl/from_genome_counts/
    qtl_data = {
        'RNA-seq\n(eQTL)': {
            'replication_rate': 45.4,
            'n_replicated': 109,
            'n_tested': 240,
            'n_samples': 137,
        },
        'ATAC-seq\n(caQTL)': {
            'replication_rate': 41.9,
            'n_replicated': 6698,
            'n_tested': 15975,
            'n_samples': 138,
        }
    }

    categories = list(qtl_data.keys())
    rates = [qtl_data[c]['replication_rate'] for c in categories]

    # Create bars
    bars = ax.bar(categories, rates,
                  color=[COLORS['vermillion'], COLORS['blue']],
                  edgecolor='white', linewidth=1.5, width=0.6)

    # Value labels on top of bars
    for bar, cat in zip(bars, categories):
        val = qtl_data[cat]['replication_rate']
        n_rep = qtl_data[cat]['n_replicated']
        n_test = qtl_data[cat]['n_tested']
        ax.text(bar.get_x() + bar.get_width()/2, val + 2,
                f'{val:.1f}%', ha='center', va='bottom',
                fontsize=8, fontweight='bold')
        # Add n below the rate
        ax.text(bar.get_x() + bar.get_width()/2, val - 5,
                f'{n_rep:,}/{n_test:,}', ha='center', va='top',
                fontsize=5, color='white')

    ax.set_ylabel('QTL replication rate (%)')
    ax.set_ylim(0, 60)

    # Add reference line at 50% (chance level would be ~5-10% with FDR)
    ax.axhline(50, color='#888888', linestyle='--', alpha=0.4, linewidth=0.8)

    # Sample count annotation
    n_rna = qtl_data['RNA-seq\n(eQTL)']['n_samples']
    n_atac = qtl_data['ATAC-seq\n(caQTL)']['n_samples']
    ax.text(0.5, 0.02, f'iPSCORE CVPC: {n_rna} RNA, {n_atac} ATAC samples',
            transform=ax.transAxes, ha='center', va='bottom',
            fontsize=5, color='gray', style='italic')

    ax.set_title('E', fontsize=9, fontweight='bold', loc='left')


def generate_figure3():
    """Generate complete Figure 3.

    Compliant with Nature Methods figure specs:
    - Max width: 180mm (7.09 inches) for double column
    - Min font: 5pt, recommended 7pt
    - Colorblind-safe palette (Bang Wong)

    Layout: 5 panels showing statistical analysis workflow and results
    - Panel A: Methods schematic (statistical pipeline)
    - Panel B: QQ plots comparing 6 statistical methods
    - Panel C: Volcano plot showing allelic imbalance
    - Panel D: Imprinting gene detection (SNV vs SNV+INDEL)
    - Panel E: QTL het/hom stratification
    """

    setup_style()

    # Nature Methods compliant: 180mm max width = 7.09 inches
    # Using 7.09 x 6.0 inches for better vertical spacing
    fig = plt.figure(figsize=(7.09, 6.0))

    # Use GridSpec for flexible layout
    from matplotlib.gridspec import GridSpec

    # Create grid with explicit spacing for better control
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1.1, 1.1], height_ratios=[1, 0.95],
                         wspace=0.35, hspace=0.40,
                         left=0.08, right=0.97, top=0.92, bottom=0.08)

    # Top row: A, B, C
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[0, 2])

    # Bottom row: D, E (D and E span, leaving space on right)
    ax_d = fig.add_subplot(gs[1, 0])
    ax_e = fig.add_subplot(gs[1, 1])
    # gs[1,2] is intentionally empty for visual balance

    # Generate all panels
    panel_a_methods_schematic(ax_a)
    panel_b_qq_plots(ax_b)
    panel_c_volcano(ax_c)
    panel_d_imprinting(ax_d)
    panel_e_qtl_replication(ax_e)

    # Cleaner title without "Figure 3:" prefix (journals add their own numbering)
    fig.suptitle('Statistical Analysis Pipeline and Results',
                fontsize=10, fontweight='bold', y=0.98)

    # Save outputs
    out_dir = get_plot_path(3, 'figure3').parent
    out_dir.mkdir(parents=True, exist_ok=True)

    # Save in multiple formats
    for ext in ['png', 'pdf', 'tiff']:
        path = get_plot_path(3, 'figure3', ext)
        plt.savefig(path, dpi=300, facecolor='white')
        print(f"Saved: {path}")

    return fig


if __name__ == '__main__':
    generate_figure3()
