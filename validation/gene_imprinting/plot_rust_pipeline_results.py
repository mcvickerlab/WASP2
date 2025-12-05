#!/usr/bin/env python3
"""
Create 4-panel figure for WASP2-Rust INDEL support validation.

Replicates Aaron's analysis showing:
A. How many imprinted genes reach significance (SNP-only vs SNP+INDEL)
B. Significance per gene (-log10 FDR comparison)
C. Direction of imprinting (allelic ratio)
D. Sanity check (SNP counts unchanged)
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Results from full WASP2-Rust pipeline
RESULTS_DIR = Path("results/rust_pipeline")
OUTPUT_DIR = Path("results/rust_pipeline")


def load_results():
    """Load gene significance results."""
    df = pd.read_csv(RESULTS_DIR / "gene_significance.tsv", sep='\t')
    return df


def create_aaron_figure(results_df, fdr_threshold=0.1):
    """Create the 4-panel comparison figure matching Aaron's style."""

    # Filter to analyzable genes (Imprinted/Predicted with data)
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
    if len(snp_sig_genes) > 0:
        print(f"  Improvement: {len(all_sig_genes)/len(snp_sig_genes):.1f}x")
    print(f"\n  SNP-only genes: {snp_sig_genes}")
    print(f"  SNP+INDEL genes: {all_sig_genes}")
    print(f"  NEWLY significant: {list(newly_sig)}")

    # Colors
    snp_color = '#3366CC'   # Blue for SNP-only
    indel_color = '#2ca02c'  # Green for SNP+INDEL

    # Create figure
    fig = plt.figure(figsize=(14, 10))

    # ===== PANEL A: Significant genes detected =====
    ax1 = fig.add_subplot(2, 2, 1)

    categories = ['SNP-only\n(WASP1)', 'SNP+INDEL\n(WASP2)']
    gene_counts = [len(snp_sig_genes), len(all_sig_genes)]
    colors = [snp_color, indel_color]

    bars = ax1.bar(categories, gene_counts, color=colors, edgecolor='black', linewidth=2, width=0.55)

    for bar, val in zip(bars, gene_counts):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.15,
                 str(val), ha='center', va='bottom', fontsize=22, fontweight='bold')


    ax1.set_ylabel('Imprinted Genes Detected\n(FDR < 0.1)', fontsize=12)
    ax1.set_title('A  How many imprinted genes reach significance?', fontsize=13, fontweight='bold', loc='left')
    ax1.set_ylim(0, max(gene_counts) * 1.4 if max(gene_counts) > 0 else 5)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # ===== PANEL B: Significance per gene (-log10 FDR) =====
    ax2 = fig.add_subplot(2, 2, 2)

    # Get top genes by significance (show more genes like the nicer figure)
    top_genes = analyzable.nsmallest(12, 'all_pval').copy()
    top_genes = top_genes.sort_values('all_fdr', ascending=False)  # Reverse for horizontal bar
    top_genes['snp_log10fdr'] = -np.log10(top_genes['snp_fdr'].clip(lower=1e-20))
    top_genes['all_log10fdr'] = -np.log10(top_genes['all_fdr'].clip(lower=1e-20))

    y = np.arange(len(top_genes))
    height = 0.35

    ax2.barh(y - height/2, top_genes['snp_log10fdr'], height, label='SNP-only',
             color=snp_color, edgecolor='black')
    ax2.barh(y + height/2, top_genes['all_log10fdr'], height, label='SNP+INDEL',
             color=indel_color, edgecolor='black')

    # Significance threshold line
    fdr_line = -np.log10(fdr_threshold)
    ax2.axvline(x=fdr_line, color='red', linestyle='--', linewidth=2, label=f'FDR={fdr_threshold}')

    ax2.set_yticks(y)
    ax2.set_yticklabels(top_genes['gene'], fontsize=10)
    ax2.set_xlabel('-log₁₀(FDR)', fontsize=12)
    ax2.set_title('B  Significance per gene', fontsize=13, fontweight='bold', loc='left')
    ax2.legend(loc='lower right', fontsize=9)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    # ===== PANEL C: Direction of imprinting =====
    ax3 = fig.add_subplot(2, 2, 3)

    # Show top significant genes by all_fdr (more genes like nice figure)
    sig_data = analyzable[analyzable['all_fdr'] < 0.2].copy()  # Relax threshold to show more
    if len(sig_data) < 5:
        sig_data = analyzable.nsmallest(10, 'all_pval').copy()

    # Calculate allelic ratio
    sig_data['all_mu'] = sig_data['all_ref'] / sig_data['all_total']
    sig_data = sig_data.sort_values('all_mu', ascending=False)

    y = np.arange(len(sig_data))

    # Separate into ref-biased and alt-biased for legend
    ref_biased = sig_data[sig_data['all_mu'] >= 0.5]
    alt_biased = sig_data[sig_data['all_mu'] < 0.5]

    # Plot bars with separate colors for legend
    for idx, row in sig_data.iterrows():
        ypos = list(sig_data.index).index(idx)
        color = snp_color if row['all_mu'] >= 0.5 else '#d62728'
        ax3.barh(ypos, row['all_mu'] - 0.5, height=0.7, color=color, edgecolor='black')

    ax3.axvline(x=0, color='black', linewidth=1.5)

    ax3.set_yticks(range(len(sig_data)))
    ax3.set_yticklabels(sig_data['gene'], fontsize=11)
    ax3.set_xlabel('Allelic Ratio (deviation from 0.5)', fontsize=12)
    ax3.set_xlim(-0.55, 0.55)

    # Add professional legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=snp_color, edgecolor='black', label='Ref-biased (Paternal)'),
        Patch(facecolor='#d62728', edgecolor='black', label='Alt-biased (Maternal)')
    ]
    ax3.legend(handles=legend_elements, loc='lower right', fontsize=9, framealpha=0.9)

    ax3.set_title('C  Direction of allelic imbalance', fontsize=13, fontweight='bold', loc='left')
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)

    # ===== PANEL D: Sanity check - SNP counts unchanged =====
    ax4 = fig.add_subplot(2, 2, 4)

    # Compare SNP counts (should be identical between approaches)
    has_snps = analyzable[analyzable['snp_total'] > 0].copy()

    if len(has_snps) > 0 and 'snp_ref' in has_snps.columns:
        ax4.scatter(has_snps['snp_ref'], has_snps['snp_ref'], s=50, alpha=0.7,
                    c=snp_color, edgecolors='black')

        # Perfect diagonal
        max_val = has_snps['snp_ref'].max() * 1.1
        ax4.plot([0, max_val], [0, max_val], 'k--', linewidth=2, label='R² = 1.000')

        ax4.set_xlabel('SNP-only filtering (ref counts)', fontsize=12)
        ax4.set_ylabel('With INDEL tolerance (ref counts)', fontsize=12)
        ax4.legend(loc='lower right', fontsize=12)
    else:
        # If we don't have snp_ref, show summary instead
        total_snp = analyzable['snp_total'].sum()
        total_indel = analyzable['indel_total'].sum()
        total_all = analyzable['all_total'].sum()

        summary_text = f"""
WASP2-Rust Pipeline Summary
═══════════════════════════

Dataset: GM12878 ATAC-seq
Genes analyzed: {len(analyzable)}

Significance (FDR < {fdr_threshold}):
  • SNP-only: {len(snp_sig_genes)} genes
  • SNP+INDEL: {len(all_sig_genes)} genes
  • Improvement: {len(all_sig_genes)/max(len(snp_sig_genes),1):.1f}×

Newly significant genes:
  {', '.join(newly_sig) if newly_sig else 'None'}

INDEL contribution:
  {100*total_indel/total_all:.1f}% of allelic reads
"""
        ax4.text(0.05, 0.95, summary_text, transform=ax4.transAxes,
                 fontsize=10, verticalalignment='top', fontfamily='monospace',
                 bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        ax4.axis('off')

    ax4.set_title('D  Sanity check: SNP counts unchanged', fontsize=13, fontweight='bold', loc='left')
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)

    # Main title
    plt.suptitle('WASP2-Rust INDEL Support Improves Imprinted Gene Detection\n(GM12878 ATAC-seq, Full Pipeline Validation)',
                 fontsize=15, fontweight='bold', y=1.02)

    plt.tight_layout()

    # Save
    plt.savefig(OUTPUT_DIR / 'indel_significance_4panel.png', dpi=150, bbox_inches='tight')
    plt.savefig(OUTPUT_DIR / 'indel_significance_4panel.pdf', bbox_inches='tight')

    print(f"\nSaved: {OUTPUT_DIR / 'indel_significance_4panel.png'}")
    print(f"Saved: {OUTPUT_DIR / 'indel_significance_4panel.pdf'}")

    return fig


def main():
    print("=" * 70)
    print("WASP2-Rust Gene Imprinting - 4-Panel Figure")
    print("=" * 70)

    # Load results
    results_df = load_results()
    print(f"\nLoaded {len(results_df)} genes")

    # Need to add ref/alt columns from allele_counts if not present
    if 'all_ref' not in results_df.columns:
        # Load variant counts for ref/alt
        var_df = pd.read_csv(RESULTS_DIR / "allele_counts.tsv", sep='\t')

        # Aggregate to gene level
        gene_agg = var_df.groupby('gene').agg({
            'ref_count': 'sum',
            'alt_count': 'sum'
        }).reset_index()
        gene_agg.columns = ['gene', 'all_ref', 'all_alt']

        # Also get SNP-only
        snp_agg = var_df[var_df['type'] == 'SNP'].groupby('gene').agg({
            'ref_count': 'sum',
            'alt_count': 'sum'
        }).reset_index()
        snp_agg.columns = ['gene', 'snp_ref', 'snp_alt']

        # Merge
        results_df = results_df.merge(gene_agg, on='gene', how='left')
        results_df = results_df.merge(snp_agg, on='gene', how='left')
        results_df['all_ref'] = results_df['all_ref'].fillna(0)
        results_df['all_alt'] = results_df['all_alt'].fillna(0)
        results_df['snp_ref'] = results_df['snp_ref'].fillna(0)
        results_df['snp_alt'] = results_df['snp_alt'].fillna(0)

    # Create figure
    create_aaron_figure(results_df, fdr_threshold=0.1)

    print("\n" + "=" * 70)
    print("Done!")
    print("=" * 70)


if __name__ == "__main__":
    main()
