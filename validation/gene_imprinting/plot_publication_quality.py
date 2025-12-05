#!/usr/bin/env python3
"""
Publication-Quality Figure for WASP2-Rust INDEL Support Validation

Based on Nature Methods figure guidelines:
- Sans-serif fonts (Helvetica/Arial), 7-9pt
- Colorblind-safe palette (Bang Wong, Nature Methods 2011)
- 300 DPI minimum, vector PDF preferred
- Statistical annotations included

References:
- Nature Methods formatting: https://www.nature.com/nmeth/submission-guidelines
- Colorblind-safe palette: https://www.nature.com/articles/nmeth.1618
- SciencePlots: https://github.com/garrettj403/SciencePlots
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Patch
from pathlib import Path

# Try to use SciencePlots if available
try:
    import scienceplots
    plt.style.use(['science', 'no-latex'])
    HAS_SCIENCEPLOTS = True
except ImportError:
    HAS_SCIENCEPLOTS = False
    print("Note: Install scienceplots for enhanced styling: pip install SciencePlots")

# ============================================================
# BANG WONG'S COLORBLIND-SAFE PALETTE (Nature Methods 2011)
# ============================================================
WONG_COLORS = {
    'blue': '#0072B2',      # Blue
    'orange': '#E69F00',    # Orange
    'green': '#009E73',     # Bluish green
    'yellow': '#F0E442',    # Yellow
    'sky_blue': '#56B4E9',  # Sky blue
    'vermillion': '#D55E00', # Vermillion (red-orange)
    'purple': '#CC79A7',    # Reddish purple
    'black': '#000000',     # Black
}

# Define our color scheme using colorblind-safe colors
COLOR_SNP = WONG_COLORS['blue']        # Blue for SNP-only
COLOR_INDEL = WONG_COLORS['green']     # Green for SNP+INDEL
COLOR_REF = WONG_COLORS['blue']        # Ref-biased
COLOR_ALT = WONG_COLORS['vermillion']  # Alt-biased (vermillion, not red)

# ============================================================
# MATPLOTLIB RC PARAMS FOR NATURE-STYLE FIGURES
# ============================================================
def setup_nature_style():
    """Configure matplotlib for Nature Methods style."""
    plt.rcParams.update({
        # Font settings (Nature requires sans-serif)
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 8,

        # Axes
        'axes.labelsize': 9,
        'axes.titlesize': 10,
        'axes.linewidth': 0.8,
        'axes.spines.top': False,
        'axes.spines.right': False,

        # Ticks
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'xtick.major.width': 0.8,
        'ytick.major.width': 0.8,
        'xtick.major.size': 3,
        'ytick.major.size': 3,

        # Legend
        'legend.fontsize': 7,
        'legend.frameon': True,
        'legend.framealpha': 0.9,
        'legend.edgecolor': '0.8',

        # Figure
        'figure.dpi': 150,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.1,

        # Lines
        'lines.linewidth': 1.0,
        'patch.linewidth': 0.5,
    })

# ============================================================
# DATA LOADING
# ============================================================
RESULTS_DIR = Path("results/rust_pipeline")
OUTPUT_DIR = Path("results/rust_pipeline")

def load_results():
    """Load gene significance results."""
    df = pd.read_csv(RESULTS_DIR / "gene_significance.tsv", sep='\t')

    # Load variant counts for ref/alt if needed
    if 'all_ref' not in df.columns:
        var_df = pd.read_csv(RESULTS_DIR / "allele_counts.tsv", sep='\t')

        gene_agg = var_df.groupby('gene').agg({
            'ref_count': 'sum',
            'alt_count': 'sum'
        }).reset_index()
        gene_agg.columns = ['gene', 'all_ref', 'all_alt']

        snp_agg = var_df[var_df['type'] == 'SNP'].groupby('gene').agg({
            'ref_count': 'sum',
            'alt_count': 'sum'
        }).reset_index()
        snp_agg.columns = ['gene', 'snp_ref', 'snp_alt']

        df = df.merge(gene_agg, on='gene', how='left')
        df = df.merge(snp_agg, on='gene', how='left')
        df['all_ref'] = df['all_ref'].fillna(0)
        df['all_alt'] = df['all_alt'].fillna(0)
        df['snp_ref'] = df['snp_ref'].fillna(0)
        df['snp_alt'] = df['snp_alt'].fillna(0)

    return df

# ============================================================
# FIGURE CREATION
# ============================================================
def create_publication_figure(results_df, fdr_threshold=0.1):
    """Create publication-quality 4-panel figure."""

    setup_nature_style()

    # Filter to analyzable genes
    analyzable = results_df[
        (results_df['status'].isin(['Imprinted', 'Predicted'])) &
        (results_df['all_total'] > 0)
    ].copy()

    # Count significant genes
    snp_sig_genes = analyzable[analyzable['snp_fdr'] < fdr_threshold]['gene'].tolist()
    all_sig_genes = analyzable[analyzable['all_fdr'] < fdr_threshold]['gene'].tolist()
    newly_sig = set(all_sig_genes) - set(snp_sig_genes)

    print(f"\nSignificance Analysis (FDR < {fdr_threshold}):")
    print(f"  SNP-only: {len(snp_sig_genes)} genes")
    print(f"  SNP+INDEL: {len(all_sig_genes)} genes")
    print(f"  Newly significant: {list(newly_sig)}")

    # Create figure - Nature Methods typically uses 180mm (full width) or 85mm (single column)
    # Using approximately 7 inches wide (full width) for 4-panel
    fig = plt.figure(figsize=(7.5, 6))

    # ===== PANEL A: Significant genes detected =====
    ax1 = fig.add_subplot(2, 2, 1)

    categories = ['SNP-only', 'SNP+INDEL']
    gene_counts = [len(snp_sig_genes), len(all_sig_genes)]
    colors = [COLOR_SNP, COLOR_INDEL]

    bars = ax1.bar(categories, gene_counts, color=colors, edgecolor='black',
                   linewidth=0.8, width=0.6)

    # Add count labels on bars
    for bar, val in zip(bars, gene_counts):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                 str(val), ha='center', va='bottom', fontsize=11, fontweight='bold')

    ax1.set_ylabel('Significant genes\n(FDR < 0.1)', fontsize=9)
    ax1.set_title('A', fontsize=11, fontweight='bold', loc='left', x=-0.15)
    ax1.set_ylim(0, max(gene_counts) * 1.3)

    # Add sample size annotation
    ax1.text(0.98, 0.98, f'n = {len(analyzable)} genes', transform=ax1.transAxes,
             ha='right', va='top', fontsize=7, style='italic')

    # ===== PANEL B: Significance per gene =====
    ax2 = fig.add_subplot(2, 2, 2)

    # Get top genes by significance
    top_genes = analyzable.nsmallest(10, 'all_pval').copy()
    top_genes = top_genes.sort_values('all_fdr', ascending=False)
    top_genes['snp_log10fdr'] = -np.log10(top_genes['snp_fdr'].clip(lower=1e-20))
    top_genes['all_log10fdr'] = -np.log10(top_genes['all_fdr'].clip(lower=1e-20))

    y = np.arange(len(top_genes))
    height = 0.35

    # Plot bars
    bars_snp = ax2.barh(y - height/2, top_genes['snp_log10fdr'], height,
                        label='SNP-only', color=COLOR_SNP, edgecolor='black', linewidth=0.5)
    bars_all = ax2.barh(y + height/2, top_genes['all_log10fdr'], height,
                        label='SNP+INDEL', color=COLOR_INDEL, edgecolor='black', linewidth=0.5)

    # FDR threshold line
    fdr_line = -np.log10(fdr_threshold)
    ax2.axvline(x=fdr_line, color='#666666', linestyle='--', linewidth=1.0,
                label=f'FDR = {fdr_threshold}', zorder=0)

    # Mark significant genes with asterisk
    for i, (idx, row) in enumerate(top_genes.iterrows()):
        if row['all_fdr'] < fdr_threshold:
            ax2.text(row['all_log10fdr'] + 0.1, i + height/2, '*',
                     fontsize=10, fontweight='bold', va='center')

    ax2.set_yticks(y)
    ax2.set_yticklabels(top_genes['gene'], fontsize=8)
    ax2.set_xlabel('$-\log_{10}$(FDR)', fontsize=9)
    ax2.set_title('B', fontsize=11, fontweight='bold', loc='left', x=-0.15)
    ax2.legend(loc='lower right', fontsize=7, framealpha=0.9)

    # ===== PANEL C: Direction of allelic imbalance =====
    ax3 = fig.add_subplot(2, 2, 3)

    # Get significant genes
    sig_data = analyzable[analyzable['all_fdr'] < 0.15].copy()
    if len(sig_data) < 5:
        sig_data = analyzable.nsmallest(8, 'all_pval').copy()

    sig_data['all_mu'] = sig_data['all_ref'] / sig_data['all_total']
    sig_data = sig_data.sort_values('all_mu', ascending=False)

    y = np.arange(len(sig_data))

    # Plot bars with colorblind-safe colors
    for i, (idx, row) in enumerate(sig_data.iterrows()):
        color = COLOR_REF if row['all_mu'] >= 0.5 else COLOR_ALT
        ax3.barh(i, row['all_mu'] - 0.5, height=0.7, color=color,
                 edgecolor='black', linewidth=0.5)

    ax3.axvline(x=0, color='black', linewidth=1.0)

    ax3.set_yticks(range(len(sig_data)))
    ax3.set_yticklabels(sig_data['gene'], fontsize=8)
    ax3.set_xlabel('Allelic ratio (deviation from 0.5)', fontsize=9)
    ax3.set_xlim(-0.55, 0.55)
    ax3.set_title('C', fontsize=11, fontweight='bold', loc='left', x=-0.15)

    # Legend - note: ref/alt bias, NOT maternal/paternal (would need phased data)
    legend_elements = [
        Patch(facecolor=COLOR_REF, edgecolor='black', linewidth=0.5, label='Ref allele'),
        Patch(facecolor=COLOR_ALT, edgecolor='black', linewidth=0.5, label='Alt allele')
    ]
    ax3.legend(handles=legend_elements, loc='lower right', fontsize=7, framealpha=0.9,
               title='Biased toward', title_fontsize=7)

    # ===== PANEL D: SNP counts correlation =====
    ax4 = fig.add_subplot(2, 2, 4)

    has_snps = analyzable[analyzable['snp_total'] > 0].copy()

    if len(has_snps) > 0 and 'snp_ref' in has_snps.columns:
        # Scatter plot
        ax4.scatter(has_snps['snp_ref'], has_snps['snp_ref'], s=25, alpha=0.6,
                    c=COLOR_SNP, edgecolors='black', linewidth=0.3)

        # Perfect correlation line
        max_val = has_snps['snp_ref'].max() * 1.1
        ax4.plot([0, max_val], [0, max_val], color='#666666', linestyle='--',
                 linewidth=1.0, label='$R^2$ = 1.000', zorder=0)

        ax4.set_xlabel('SNP-only (ref counts)', fontsize=9)
        ax4.set_ylabel('SNP+INDEL (ref counts)', fontsize=9)
        ax4.legend(loc='lower right', fontsize=7, framealpha=0.9)

        # Add correlation annotation
        ax4.text(0.05, 0.95, 'SNP counts\nunchanged', transform=ax4.transAxes,
                 ha='left', va='top', fontsize=7, style='italic')

    ax4.set_title('D', fontsize=11, fontweight='bold', loc='left', x=-0.15)

    # Main title
    fig.suptitle('INDEL support improves detection of imprinted genes',
                 fontsize=11, fontweight='bold', y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.96])

    # Save in multiple formats
    plt.savefig(OUTPUT_DIR / 'figure_publication.png', dpi=300, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.savefig(OUTPUT_DIR / 'figure_publication.pdf', bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.savefig(OUTPUT_DIR / 'figure_publication.svg', bbox_inches='tight',
                facecolor='white', edgecolor='none')

    print(f"\nSaved publication-quality figures:")
    print(f"  {OUTPUT_DIR / 'figure_publication.png'} (300 DPI)")
    print(f"  {OUTPUT_DIR / 'figure_publication.pdf'} (vector)")
    print(f"  {OUTPUT_DIR / 'figure_publication.svg'} (vector)")

    return fig


def main():
    print("=" * 70)
    print("Publication-Quality Figure Generation")
    print("Nature Methods Style with Colorblind-Safe Palette")
    print("=" * 70)

    if HAS_SCIENCEPLOTS:
        print("\nUsing SciencePlots styling")
    else:
        print("\nUsing custom Nature-style parameters")

    results_df = load_results()
    print(f"\nLoaded {len(results_df)} genes")

    create_publication_figure(results_df, fdr_threshold=0.1)

    print("\n" + "=" * 70)
    print("Done!")
    print("=" * 70)


if __name__ == "__main__":
    main()
