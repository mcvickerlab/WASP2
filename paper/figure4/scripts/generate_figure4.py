#!/usr/bin/env python3
"""
Figure 4: Single Cell ATAC-seq Allelic Analysis

Shows WASP2's capability for single-cell allelic imbalance analysis.
Data: GM12878 scATAC-seq from 10X Genomics (pre-processed by Aaron Ho)

K-Dense scanpy skill compliant:
- Colorblind-safe visualization
- Publication-quality figures (Nature Methods specs)
"""
import sys
from pathlib import Path

# Add paper directory to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from config import (
    COLORS, PLOT_SETTINGS,
    get_plot_path, get_data_path
)

# Paths to Aaron's pre-computed data
SC_DATA_DIR = Path("/iblm/netapp/data3/aho/project_data/wasp2/10x_cellranger_atac/gm12878_el4/gm12878_filt")
CELL_QC_PATH = SC_DATA_DIR / "gm12878_singlecell.csv"
ALLELIC_COUNTS_PATH = SC_DATA_DIR / "counts" / "gm12878_filt_all_peaks.h5ad"


def setup_style():
    """Set up Nature Methods plotting style (k-dense compliant)."""
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


def load_cell_qc():
    """Load cell-level QC metrics from 10X CellRanger output."""
    if CELL_QC_PATH.exists():
        df = pd.read_csv(CELL_QC_PATH)
        # Filter to valid cells (passed QC)
        df = df[df['is_GRCh38_cell_barcode'] == 1].copy()
        print(f"Loaded cell QC: {len(df):,} cells")
        return df
    return None


def load_allelic_counts():
    """Load allelic count data from WASP2 output."""
    if ALLELIC_COUNTS_PATH.exists():
        try:
            import scanpy as sc
            adata = sc.read_h5ad(ALLELIC_COUNTS_PATH)
            # This contains variant-level allelic counts
            df = adata.obs.copy()
            print(f"Loaded allelic counts: {len(df):,} variants")
            return df
        except Exception as e:
            print(f"Error loading allelic counts: {e}")
    return None


def panel_a_cell_quality(ax, cell_df):
    """Panel A: Cell quality - passed_filters distribution."""
    n_cells = 0
    if cell_df is not None and 'passed_filters' in cell_df.columns:
        counts = cell_df['passed_filters'].values
        n_cells = len(cell_df)
        ax.hist(np.log10(counts + 1), bins=50, color=COLORS['blue'],
                edgecolor='white', linewidth=0.5, alpha=0.95)
        ax.set_xlabel('log₁₀(Passed filter reads + 1)')
        ax.set_ylabel('Cells')

        median_val = np.median(counts)
        ax.axvline(np.log10(median_val + 1), color=COLORS['vermillion'],
                   linestyle='--', linewidth=1.5, label=f'Median: {int(median_val):,}')
        ax.legend(fontsize=6, frameon=False, loc='upper left')
    else:
        ax.text(0.5, 0.5, 'Cell Quality\n(QC metrics)',
                transform=ax.transAxes, ha='center', va='center',
                fontsize=9, style='italic', color='gray')
    # Clean panel label
    ax.set_title('A', fontsize=9, fontweight='bold', loc='left')
    # Keep y-axis label simple, add n as text annotation
    ax.set_ylabel('Cells')
    ax.text(0.98, 0.98, f'n = {n_cells:,}', transform=ax.transAxes,
            ha='right', va='top', fontsize=6, color='gray')


def panel_b_tss_enrichment(ax, cell_df):
    """Panel B: TSS enrichment as quality indicator."""
    if cell_df is not None and 'TSS_fragments' in cell_df.columns:
        tss = cell_df['TSS_fragments'].values
        total = cell_df['passed_filters'].values
        # TSS enrichment ratio (TSS fragments / total fragments)
        tss_ratio = tss / (total + 1)

        ax.hist(tss_ratio, bins=50, color=COLORS['green'],
                edgecolor='white', linewidth=0.5, alpha=0.95)
        ax.set_xlabel('TSS fragment fraction')
        ax.set_ylabel('Cells')

        median_ratio = np.median(tss_ratio)
        ax.axvline(median_ratio, color=COLORS['vermillion'],
                   linestyle='--', linewidth=1.5, label=f'Median: {median_ratio:.2f}')
        ax.legend(fontsize=7, frameon=False, loc='upper right')
    else:
        ax.text(0.5, 0.5, 'TSS Enrichment\nDistribution',
                transform=ax.transAxes, ha='center', va='center',
                fontsize=9, style='italic', color='gray')
    # Clean panel label
    ax.set_title('B', fontsize=9, fontweight='bold', loc='left')


def panel_c_allelic_ratio(ax, allelic_df):
    """Panel C: Allelic ratio distribution at heterozygous sites."""
    n_sites = 0
    median_ratio = 0.0
    if allelic_df is not None and 'ref_count' in allelic_df.columns:
        ref = allelic_df['ref_count'].values
        alt = allelic_df['alt_count'].values
        total = ref + alt

        # Filter sites with sufficient coverage
        mask = total >= 10
        ref_filt = ref[mask]
        alt_filt = alt[mask]
        total_filt = total[mask]

        if len(ref_filt) > 0:
            ratio = ref_filt / total_filt
            n_sites = len(ratio)
            median_ratio = np.median(ratio)

            ax.hist(ratio, bins=50, color=COLORS['orange'],
                    edgecolor='white', linewidth=0.5, alpha=0.95)
            ax.set_xlabel('Ref allele ratio')
            ax.set_ylabel('Variant sites')
            ax.axvline(0.5, color='gray', linestyle='--', linewidth=1, alpha=0.7,
                       label='Expected (0.5)')
            ax.legend(fontsize=6, frameon=False, loc='upper left')
        else:
            ax.text(0.5, 0.5, 'Insufficient coverage',
                    transform=ax.transAxes, ha='center', va='center',
                    fontsize=9, style='italic', color='gray')
    else:
        ax.text(0.5, 0.5, 'Allelic Ratio\nDistribution',
                transform=ax.transAxes, ha='center', va='center',
                fontsize=9, style='italic', color='gray')
    # Clean panel label
    ax.set_title('C', fontsize=9, fontweight='bold', loc='left')
    # Keep y-axis label simple, add stats as text annotation
    ax.set_ylabel('Variant sites')
    ax.text(0.98, 0.98, f'n = {n_sites:,}\nmedian = {median_ratio:.3f}',
            transform=ax.transAxes, ha='right', va='top', fontsize=6, color='gray')


def panel_d_coverage_vs_power(ax, allelic_df):
    """Panel D: Coverage vs statistical power for allelic imbalance."""
    pct = 0.0
    if allelic_df is not None and 'ref_count' in allelic_df.columns:
        ref = allelic_df['ref_count'].values
        alt = allelic_df['alt_count'].values
        total = ref + alt

        # Bin by coverage and show fraction of sites
        coverage_bins = [0, 10, 20, 50, 100, 200, 500, 1000, np.inf]
        bin_labels = ['0-10', '10-20', '20-50', '50-100', '100-200', '200-500', '500-1K', '>1K']

        counts = []
        for i in range(len(coverage_bins) - 1):
            mask = (total >= coverage_bins[i]) & (total < coverage_bins[i+1])
            counts.append(mask.sum())

        bars = ax.bar(range(len(bin_labels)), counts, color=COLORS['blue'],
                      edgecolor='black', linewidth=0.5, alpha=0.95)
        ax.set_xlabel('Read coverage')
        ax.set_ylabel('Variant sites')
        ax.set_xticks(range(len(bin_labels)))
        ax.set_xticklabels(bin_labels, rotation=45, ha='right', fontsize=6)

        # Annotate sufficient coverage threshold (removed confusing horizontal line)
        total_sites = len(total)
        high_cov = sum(counts[2:])  # >=20 reads
        pct = 100 * high_cov / total_sites if total_sites > 0 else 0
    else:
        ax.text(0.5, 0.5, 'Coverage\nDistribution',
                transform=ax.transAxes, ha='center', va='center',
                fontsize=9, style='italic', color='gray')
    # Clean panel label
    ax.set_title('D', fontsize=9, fontweight='bold', loc='left')
    # Keep y-axis label simple, add stats as text annotation
    ax.set_ylabel('Variant sites')
    ax.text(0.98, 0.98, f'{pct:.1f}% sites\n\u226520x coverage',
            transform=ax.transAxes, ha='right', va='top', fontsize=6, color='gray')


def generate_figure4():
    """Generate Figure 4: Single Cell Analysis.

    K-Dense compliant:
    - Nature Methods dimensions (8x7 inches for 2x2 layout)
    - Colorblind-safe palette
    - Min 7pt fonts
    - Uses constrained_layout for optimal spacing
    """
    setup_style()

    # Load data
    cell_df = load_cell_qc()
    allelic_df = load_allelic_counts()

    # Use constrained_layout for better handling of suptitle, legends, rotated labels
    fig, axes = plt.subplots(2, 2, figsize=(8, 7), layout="constrained")
    ax_a, ax_b = axes[0]
    ax_c, ax_d = axes[1]

    panel_a_cell_quality(ax_a, cell_df)
    panel_b_tss_enrichment(ax_b, cell_df)
    panel_c_allelic_ratio(ax_c, allelic_df)
    panel_d_coverage_vs_power(ax_d, allelic_df)

    # Add data source annotation at bottom left to avoid overlap with Panel D labels
    n_cells = len(cell_df) if cell_df is not None else 0
    n_variants = len(allelic_df) if allelic_df is not None else 0
    fig.text(0.01, 0.01, f'GM12878 scATAC-seq (pseudobulk) | {n_cells:,} cells, {n_variants:,} variants',
             ha='left', va='bottom', fontsize=6, color='gray')

    fig.suptitle('Figure 4: Single-Cell ATAC Analysis', fontsize=10, fontweight='bold')

    # Save
    out_dir = get_plot_path(4, 'figure4').parent
    out_dir.mkdir(parents=True, exist_ok=True)

    plt.savefig(get_plot_path(4, 'figure4', 'png'), bbox_inches='tight', facecolor='white')
    plt.savefig(get_plot_path(4, 'figure4', 'pdf'), bbox_inches='tight', facecolor='white')

    print(f"Saved: {get_plot_path(4, 'figure4', 'png')}")
    print(f"Saved: {get_plot_path(4, 'figure4', 'pdf')}")

    return fig


if __name__ == '__main__':
    generate_figure4()
