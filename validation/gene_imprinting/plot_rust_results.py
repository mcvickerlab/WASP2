#!/usr/bin/env python3
"""Plot gene imprinting results from WASP2-Rust BamCounter.

Replicates Aaron's analysis style using matplotlib.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats


def plot_imprinted_scatter(df, ax, title=None):
    """Create Aaron-style scatter plot: mu (ref proportion) vs -log10p"""
    if title is None:
        title = "Allele-Specific Expression - Imprinted Genes (WASP2-Rust)"

    # Color map matching Aaron's notebook
    cmap = {
        "Imprinted": "#DC3912",  # Red
        "Predicted": "#3366CC",  # Blue
        "Not Imprinted": "gray"
    }

    for status, color in cmap.items():
        subset = df[df['Status'] == status]
        if len(subset) > 0:
            ax.scatter(subset['mu'], subset['log10p'],
                       c=color, label=status, s=60, alpha=0.7, edgecolors='black', linewidth=0.5)

    ax.set_xlabel("Reference Proportion (μ)", fontsize=12)
    ax.set_ylabel("-log₁₀(p)", fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.set_xlim(-0.05, 1.05)
    ax.axhline(y=-np.log10(0.05), color='gray', linestyle='--', alpha=0.5, label='p=0.05')
    ax.legend()
    ax.grid(True, alpha=0.3)


def main():
    # Load the Rust counts
    df = pd.read_csv('results/wasp2_rust_counts.tsv', sep='\t')

    # Filter to variants with reads
    df_with_reads = df[df['total'] > 0].copy()

    # Compute mu (reference proportion)
    df_with_reads['mu'] = df_with_reads['ref_count'] / df_with_reads['total']

    # Classify variant type
    df_with_reads['variant_type'] = df_with_reads.apply(
        lambda x: 'INDEL' if len(x['ref']) != len(x['alt']) else 'SNP', axis=1
    )

    # Compute simple binomial p-value for allelic imbalance
    def compute_pval(row):
        n = row['ref_count'] + row['alt_count']
        k = row['ref_count']
        if n == 0:
            return 1.0
        # Two-sided binomial test against 0.5
        return stats.binom_test(k, n, 0.5)

    df_with_reads['pval'] = df_with_reads.apply(compute_pval, axis=1)
    df_with_reads['log10p'] = -np.log10(df_with_reads['pval'].clip(lower=1e-300))

    # Add imprinting status (all these genes are known imprinted)
    df_with_reads['Status'] = 'Imprinted'

    print("=" * 60)
    print("WASP2-Rust Gene Imprinting Analysis")
    print("=" * 60)

    # Summary per gene
    print("\n--- Per-Gene Summary ---")
    gene_summary = []
    for gene in df_with_reads['gene'].unique():
        gene_df = df_with_reads[df_with_reads['gene'] == gene]
        total_ref = gene_df['ref_count'].sum()
        total_alt = gene_df['alt_count'].sum()
        total_reads = total_ref + total_alt
        n_variants = len(gene_df)
        n_snps = len(gene_df[gene_df['variant_type'] == 'SNP'])
        n_indels = len(gene_df[gene_df['variant_type'] == 'INDEL'])
        mu = total_ref / total_reads if total_reads > 0 else 0.5

        # Gene-level p-value
        if total_reads > 0:
            pval = stats.binom_test(total_ref, total_reads, 0.5)
        else:
            pval = 1.0

        gene_summary.append({
            'gene': gene,
            'ref_count': total_ref,
            'alt_count': total_alt,
            'N': total_reads,
            'n_variants': n_variants,
            'n_snps': n_snps,
            'n_indels': n_indels,
            'mu': mu,
            'pval': pval,
            'log10p': -np.log10(max(pval, 1e-300)),
            'Status': 'Imprinted',
            'variant_type': 'mixed'
        })

        print(f"\n{gene}:")
        print(f"  Variants with reads: {n_variants} ({n_snps} SNPs, {n_indels} INDELs)")
        print(f"  Total reads: {total_reads} (ref={total_ref}, alt={total_alt})")
        print(f"  Gene-level μ: {mu:.3f}, p={pval:.2e}")

    gene_df = pd.DataFrame(gene_summary)
    gene_df = gene_df[gene_df['N'] > 0]

    # Create figure with multiple plots
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # 1. Per-variant scatter plot
    plot_imprinted_scatter(df_with_reads, axes[0, 0],
                           title="Per-Variant Allelic Imbalance (WASP2-Rust)")

    # Add gene labels for significant variants
    sig_variants = df_with_reads[df_with_reads['log10p'] > 2]
    for _, row in sig_variants.iterrows():
        axes[0, 0].annotate(row['gene'], (row['mu'], row['log10p']),
                            fontsize=8, alpha=0.7)

    # 2. Gene-level scatter plot
    plot_imprinted_scatter(gene_df, axes[0, 1],
                           title="Gene-Level Allelic Imbalance (WASP2-Rust)")

    # Add gene labels
    for _, row in gene_df.iterrows():
        axes[0, 1].annotate(row['gene'], (row['mu'], row['log10p']),
                            fontsize=10, fontweight='bold')

    # 3. Bar chart: reads per gene (stacked by allele)
    ax3 = axes[1, 0]
    genes = gene_df['gene'].tolist()
    x = np.arange(len(genes))
    width = 0.7

    ax3.bar(x, gene_df['ref_count'], width, label='Reference allele', color='#2166ac')
    ax3.bar(x, gene_df['alt_count'], width, bottom=gene_df['ref_count'],
            label='Alternative allele', color='#b2182b')
    ax3.set_xticks(x)
    ax3.set_xticklabels(genes, rotation=45, ha='right')
    ax3.set_ylabel('Read Count')
    ax3.set_title('Allelic Counts per Imprinted Gene')
    ax3.legend()

    # Add total labels
    for i, row in gene_df.iterrows():
        total = row['ref_count'] + row['alt_count']
        if total > 0:
            ax3.text(list(gene_df.index).index(i), total + 2,
                     str(int(total)), ha='center', va='bottom', fontsize=9)

    # 4. SNP vs INDEL contribution
    ax4 = axes[1, 1]
    snp_reads = df_with_reads[df_with_reads['variant_type'] == 'SNP']['total'].sum()
    indel_reads = df_with_reads[df_with_reads['variant_type'] == 'INDEL']['total'].sum()
    total = snp_reads + indel_reads

    snp_sites = len(df_with_reads[df_with_reads['variant_type'] == 'SNP'])
    indel_sites = len(df_with_reads[df_with_reads['variant_type'] == 'INDEL'])

    categories = ['Sites with reads', 'Total reads']
    snp_vals = [snp_sites, snp_reads]
    indel_vals = [indel_sites, indel_reads]

    x4 = np.arange(len(categories))
    width4 = 0.35

    ax4.bar(x4 - width4/2, snp_vals, width4, label='SNPs', color='#4393c3')
    ax4.bar(x4 + width4/2, indel_vals, width4, label='INDELs', color='#d6604d')
    ax4.set_xticks(x4)
    ax4.set_xticklabels(categories)
    ax4.set_ylabel('Count')
    ax4.set_title(f'SNP vs INDEL Contribution\n(INDELs: {100*indel_reads/total:.0f}% of reads)')
    ax4.legend()

    # Add labels
    for i, (s, d) in enumerate(zip(snp_vals, indel_vals)):
        ax4.text(i - width4/2, s + 1, str(s), ha='center', va='bottom', fontsize=9)
        ax4.text(i + width4/2, d + 1, str(d), ha='center', va='bottom', fontsize=9)

    plt.tight_layout()
    plt.savefig('results/gene_imprinting_rust_results.png', dpi=150, bbox_inches='tight')
    plt.savefig('results/gene_imprinting_rust_results.pdf', bbox_inches='tight')

    print("\n" + "=" * 60)
    print("Saved plots:")
    print("  results/gene_imprinting_rust_results.png")
    print("  results/gene_imprinting_rust_results.pdf")

    # Print gene summary table
    print("\n--- Gene Summary Table ---")
    print(gene_df[['gene', 'N', 'n_snps', 'n_indels', 'mu', 'pval']].to_string(index=False))

    # INDEL summary
    print("\n--- INDEL Contribution ---")
    print(f"SNP reads: {snp_reads} ({100*snp_reads/total:.1f}%)")
    print(f"INDEL reads: {indel_reads} ({100*indel_reads/total:.1f}%)")

    print("\n" + "=" * 60)
    print("Done!")


if __name__ == "__main__":
    main()
