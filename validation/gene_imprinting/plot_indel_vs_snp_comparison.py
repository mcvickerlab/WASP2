#!/usr/bin/env python3
"""
INDEL vs SNP-only Comparison Plot for Gene Imprinting Validation

This plot demonstrates the value of WASP2's INDEL support by comparing:
- SNP-only analysis (what WASP1 could do)
- SNP+INDEL analysis (what WASP2 can do)

Key metrics:
- Number of genes with coverage
- Number of variants with reads
- Total read counts
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


def main():
    # Load the Rust counts
    df = pd.read_csv('results/wasp2_rust_counts.tsv', sep='\t')

    # Classify variant type
    df['variant_type'] = df.apply(
        lambda x: 'INDEL' if len(str(x['ref'])) != len(str(x['alt'])) else 'SNP', axis=1
    )

    # Filter to variants with reads
    df_with_reads = df[df['total'] > 0].copy()

    print("=" * 70)
    print("WASP2 INDEL Support: Gene Imprinting Validation")
    print("=" * 70)

    # ========== Calculate metrics ==========
    # All variants
    total_variants = len(df)
    total_snps = len(df[df['variant_type'] == 'SNP'])
    total_indels = len(df[df['variant_type'] == 'INDEL'])

    # Variants with reads
    snp_with_reads = df_with_reads[df_with_reads['variant_type'] == 'SNP']
    indel_with_reads = df_with_reads[df_with_reads['variant_type'] == 'INDEL']

    n_snps_with_reads = len(snp_with_reads)
    n_indels_with_reads = len(indel_with_reads)

    # Read counts
    snp_reads = snp_with_reads['total'].sum()
    indel_reads = indel_with_reads['total'].sum()
    total_reads = snp_reads + indel_reads

    # Genes covered
    genes_snp_only = set(snp_with_reads['gene'].unique())
    genes_indel_only = set(indel_with_reads['gene'].unique())
    genes_both = genes_snp_only | genes_indel_only  # Union = all genes with any coverage

    # Genes ONLY covered by INDELs (not by SNPs)
    genes_indel_exclusive = genes_indel_only - genes_snp_only

    print(f"\n{'Category':<30} {'SNP-only':<15} {'+INDEL':<15} {'Improvement':<15}")
    print("-" * 70)
    print(f"{'Genes with coverage:':<30} {len(genes_snp_only):<15} {len(genes_both):<15} {len(genes_both)/max(len(genes_snp_only),1):.1f}x")
    print(f"{'Variants with reads:':<30} {n_snps_with_reads:<15} {n_snps_with_reads + n_indels_with_reads:<15} {(n_snps_with_reads + n_indels_with_reads)/max(n_snps_with_reads,1):.1f}x")
    print(f"{'Total reads:':<30} {snp_reads:<15} {total_reads:<15} {total_reads/max(snp_reads,1):.1f}x")
    print("-" * 70)

    if genes_indel_exclusive:
        print(f"\nGenes ONLY accessible via INDELs (not covered by SNPs):")
        for g in sorted(genes_indel_exclusive):
            gene_indel_df = indel_with_reads[indel_with_reads['gene'] == g]
            gene_reads = gene_indel_df['total'].sum()
            print(f"  - {g}: {len(gene_indel_df)} INDEL variants, {gene_reads} reads")

    # ========== Create the comparison figure ==========
    fig = plt.figure(figsize=(14, 10))

    # Color scheme
    snp_color = '#4393c3'  # Blue for SNPs
    indel_color = '#d6604d'  # Red/Orange for INDELs
    combined_color = '#2d7c2d'  # Green for combined benefit

    # ----- Panel A: Genes with Coverage -----
    ax1 = fig.add_subplot(2, 2, 1)

    categories = ['SNP-only\n(WASP1)', 'SNP+INDEL\n(WASP2)']
    gene_counts = [len(genes_snp_only), len(genes_both)]
    colors = [snp_color, combined_color]

    bars = ax1.bar(categories, gene_counts, color=colors, edgecolor='black', linewidth=1.5, width=0.6)

    # Add value labels
    for bar, val in zip(bars, gene_counts):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                 str(val), ha='center', va='bottom', fontsize=16, fontweight='bold')

    # Add improvement annotation
    if len(genes_snp_only) > 0:
        improvement = len(genes_both) / len(genes_snp_only)
        ax1.annotate(f'{improvement:.1f}× more genes',
                     xy=(1, gene_counts[1]),
                     xytext=(1.3, gene_counts[1] * 0.8),
                     fontsize=12, fontweight='bold', color=combined_color,
                     arrowprops=dict(arrowstyle='->', color=combined_color, lw=2))

    ax1.set_ylabel('Number of Genes', fontsize=14)
    ax1.set_title('A  Imprinted Genes with Read Coverage', fontsize=14, fontweight='bold', loc='left')
    ax1.set_ylim(0, max(gene_counts) * 1.3)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # ----- Panel B: Variants with Reads -----
    ax2 = fig.add_subplot(2, 2, 2)

    variant_counts = [n_snps_with_reads, n_snps_with_reads + n_indels_with_reads]
    bars = ax2.bar(categories, variant_counts, color=colors, edgecolor='black', linewidth=1.5, width=0.6)

    for bar, val in zip(bars, variant_counts):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                 str(val), ha='center', va='bottom', fontsize=16, fontweight='bold')

    if n_snps_with_reads > 0:
        improvement = (n_snps_with_reads + n_indels_with_reads) / n_snps_with_reads
        ax2.annotate(f'{improvement:.1f}× more variants',
                     xy=(1, variant_counts[1]),
                     xytext=(1.3, variant_counts[1] * 0.8),
                     fontsize=12, fontweight='bold', color=combined_color,
                     arrowprops=dict(arrowstyle='->', color=combined_color, lw=2))

    ax2.set_ylabel('Number of Variants', fontsize=14)
    ax2.set_title('B  Heterozygous Variants with Reads', fontsize=14, fontweight='bold', loc='left')
    ax2.set_ylim(0, max(variant_counts) * 1.3)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    # ----- Panel C: Read Counts (Stacked Bar) -----
    ax3 = fig.add_subplot(2, 2, 3)

    x = np.arange(2)
    width = 0.6

    # SNP-only shows only SNP reads, WASP2 shows both stacked
    snp_bars = ax3.bar(x, [snp_reads, snp_reads], width, label='SNP reads', color=snp_color, edgecolor='black')
    indel_bars = ax3.bar(x[1], indel_reads, width, bottom=snp_reads, label='INDEL reads', color=indel_color, edgecolor='black')

    # Value labels
    ax3.text(0, snp_reads + 2, str(snp_reads), ha='center', va='bottom', fontsize=14, fontweight='bold')
    ax3.text(1, total_reads + 2, str(total_reads), ha='center', va='bottom', fontsize=14, fontweight='bold')

    # Percentage annotation for INDEL contribution
    indel_pct = 100 * indel_reads / total_reads if total_reads > 0 else 0
    ax3.text(1, snp_reads + indel_reads/2, f'+{int(indel_reads)}\n({indel_pct:.0f}%)',
             ha='center', va='center', fontsize=11, fontweight='bold', color='white')

    ax3.set_xticks(x)
    ax3.set_xticklabels(categories)
    ax3.set_ylabel('Total Read Count', fontsize=14)
    ax3.set_title('C  Total Allelic Reads Captured', fontsize=14, fontweight='bold', loc='left')
    ax3.legend(loc='upper left', fontsize=11)
    ax3.set_ylim(0, total_reads * 1.25)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)

    # ----- Panel D: Per-Gene Breakdown -----
    ax4 = fig.add_subplot(2, 2, 4)

    # Get per-gene stats
    gene_stats = []
    for gene in df['gene'].unique():
        gene_df = df_with_reads[df_with_reads['gene'] == gene]
        snp_df = gene_df[gene_df['variant_type'] == 'SNP']
        indel_df = gene_df[gene_df['variant_type'] == 'INDEL']
        gene_stats.append({
            'gene': gene,
            'snp_reads': snp_df['total'].sum(),
            'indel_reads': indel_df['total'].sum(),
            'total_reads': snp_df['total'].sum() + indel_df['total'].sum()
        })

    gene_stats_df = pd.DataFrame(gene_stats)
    gene_stats_df = gene_stats_df.sort_values('total_reads', ascending=True)

    # Horizontal stacked bar
    y = np.arange(len(gene_stats_df))
    height = 0.7

    ax4.barh(y, gene_stats_df['snp_reads'], height, label='SNP reads', color=snp_color, edgecolor='black')
    ax4.barh(y, gene_stats_df['indel_reads'], height, left=gene_stats_df['snp_reads'],
             label='INDEL reads', color=indel_color, edgecolor='black')

    ax4.set_yticks(y)
    ax4.set_yticklabels(gene_stats_df['gene'], fontsize=11)
    ax4.set_xlabel('Read Count', fontsize=14)
    ax4.set_title('D  Per-Gene Read Contribution by Variant Type', fontsize=14, fontweight='bold', loc='left')
    ax4.legend(loc='lower right', fontsize=11)
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)

    # Add total labels
    for i, row in gene_stats_df.iterrows():
        total = row['total_reads']
        if total > 0:
            idx = list(gene_stats_df.index).index(i)
            ax4.text(total + 2, idx, str(int(total)), va='center', fontsize=10)

    plt.tight_layout()

    # Save
    plt.savefig('results/indel_vs_snp_comparison.png', dpi=150, bbox_inches='tight')
    plt.savefig('results/indel_vs_snp_comparison.pdf', bbox_inches='tight')

    print("\n" + "=" * 70)
    print("KEY FINDING: INDEL Support Value")
    print("=" * 70)
    print(f"  • SNP-only covers {len(genes_snp_only)} genes, SNP+INDEL covers {len(genes_both)} genes")
    if genes_indel_exclusive:
        print(f"  • {len(genes_indel_exclusive)} gene(s) ONLY accessible via INDELs: {', '.join(sorted(genes_indel_exclusive))}")
    print(f"  • INDELs contribute {indel_pct:.0f}% of total reads ({indel_reads}/{total_reads})")
    print(f"  • WASP2 captures {(n_snps_with_reads + n_indels_with_reads)/max(n_snps_with_reads,1):.1f}x more informative variants")
    print("=" * 70)

    print("\nSaved plots:")
    print("  results/indel_vs_snp_comparison.png")
    print("  results/indel_vs_snp_comparison.pdf")


if __name__ == "__main__":
    main()
