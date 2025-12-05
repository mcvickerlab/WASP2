#!/usr/bin/env python3
"""
Replicate Aaron's INDEL Support Analysis for Gene Imprinting

This script properly replicates the "3× more genes" finding by:
1. Using Aaron's full imprinted gene list (187 genes)
2. Using the full NA12878 VCF
3. Analyzing SNP-only vs SNP+INDEL gene coverage

The key insight: Some genes only have INDEL heterozygous variants,
making them invisible to SNP-only analysis but visible with INDEL support.
"""
import sys
import gzip
import subprocess
from pathlib import Path
from collections import defaultdict

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# Paths to Aaron's data
AARON_DATA = Path("/iblm/netapp/home/aho/projects/wasp/testing/performance/data/GM12878_ATACseq_50k_merged")
IMPRINTED_GENES = AARON_DATA / "test_imprinted/geneimprint_coords.tsv"

# VCF and sample info
FULL_VCF = "/iblm/netapp/data1/aho/variants/NA12878.vcf.gz"
SAMPLE = "NA12878"

# BAM for counting (Aaron's WASP-filtered ATAC-seq)
ATAC_BAM = AARON_DATA / "remap_results" / "GM12878_ATACseq_50k_merged_wasp_filt.bam"
# Alternative: Use raw BAM for comparison
RAW_ATAC_BAM = "/iblm/netapp/data1/aho/atac/Buenrostro2013/merged/GM12878_ATACseq_50k/GM12878_ATACseq_50k_merged.sorted.bam"


def load_imprinted_genes():
    """Load Aaron's full imprinted gene list."""
    df = pd.read_csv(IMPRINTED_GENES, sep='\t')
    print(f"Loaded {len(df)} imprinted gene entries")
    print(f"Unique genes: {df['Gene'].nunique()}")
    return df


def get_het_variants_for_gene(vcf_path, chrom, start, end, sample="NA12878"):
    """Extract heterozygous variants in a genomic region."""
    region = f"{chrom}:{start}-{end}"
    cmd = f"bcftools view -r {region} -s {sample} {vcf_path} 2>/dev/null | bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%GT]\\n'"

    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=30)
        variants = []

        for line in result.stdout.strip().split('\n'):
            if not line:
                continue
            fields = line.split('\t')
            if len(fields) < 5:
                continue

            chrom, pos, ref, alt, gt = fields

            # Check if heterozygous (0/1, 0|1, 1/0, 1|0)
            if '/' in gt or '|' in gt:
                alleles = gt.replace('|', '/').split('/')
                if len(alleles) == 2 and alleles[0] != alleles[1]:
                    # Classify as SNP or INDEL
                    is_indel = len(ref) != len(alt.split(',')[0])
                    variants.append({
                        'chrom': chrom,
                        'pos': int(pos),
                        'ref': ref,
                        'alt': alt,
                        'gt': gt,
                        'is_indel': is_indel,
                        'type': 'INDEL' if is_indel else 'SNP'
                    })

        return variants
    except Exception as e:
        print(f"  Error extracting variants from {region}: {e}")
        return []


def analyze_gene_coverage():
    """Analyze which genes have SNP-only vs SNP+INDEL coverage."""
    genes_df = load_imprinted_genes()

    results = []

    print("\nAnalyzing heterozygous variants per gene...")
    print("-" * 70)

    # Group by unique gene regions
    unique_genes = genes_df.drop_duplicates(subset=['Gene', 'chrom', 'start', 'end'])

    for idx, row in unique_genes.iterrows():
        gene = row['Gene']
        chrom = row['chrom']
        start = row['start']
        end = row['end']
        status = row['Status']

        variants = get_het_variants_for_gene(FULL_VCF, chrom, start, end, SAMPLE)

        n_snps = sum(1 for v in variants if not v['is_indel'])
        n_indels = sum(1 for v in variants if v['is_indel'])
        n_total = len(variants)

        results.append({
            'gene': gene,
            'chrom': chrom,
            'start': start,
            'end': end,
            'status': status,
            'n_snps': n_snps,
            'n_indels': n_indels,
            'n_total': n_total,
            'has_snp': n_snps > 0,
            'has_indel': n_indels > 0,
            'has_any': n_total > 0,
            'indel_only': n_indels > 0 and n_snps == 0
        })

        if idx % 20 == 0:
            print(f"  Processed {idx+1}/{len(unique_genes)} genes...")

    results_df = pd.DataFrame(results)
    return results_df


def create_comparison_figure(results_df):
    """Create the 4-panel comparison figure like Aaron's analysis."""

    # Filter to genes with actual coverage (status = Imprinted or Predicted)
    analyzable = results_df[results_df['status'].isin(['Imprinted', 'Predicted'])]

    # Calculate key metrics
    genes_with_snp = analyzable[analyzable['has_snp']]['gene'].nunique()
    genes_with_any = analyzable[analyzable['has_any']]['gene'].nunique()
    genes_indel_only = analyzable[analyzable['indel_only']]['gene'].nunique()

    print("\n" + "=" * 70)
    print("INDEL SUPPORT VALUE ANALYSIS")
    print("=" * 70)
    print(f"\nTotal analyzable imprinted/predicted genes: {len(analyzable)}")
    print(f"\nGenes with heterozygous variants:")
    print(f"  SNP-only accessible:     {genes_with_snp}")
    print(f"  SNP+INDEL accessible:    {genes_with_any}")
    print(f"  INDEL-only genes:        {genes_indel_only}")

    if genes_with_snp > 0:
        improvement = genes_with_any / genes_with_snp
        print(f"\n  Improvement factor: {improvement:.1f}×")

    # Create figure
    fig = plt.figure(figsize=(14, 10))

    # Colors
    snp_color = '#3366CC'  # Blue for SNP-only
    indel_color = '#2ca02c'  # Green for SNP+INDEL

    # ----- Panel A: Genes Detected -----
    ax1 = fig.add_subplot(2, 2, 1)

    categories = ['SNP-only\n(WASP1)', 'SNP+INDEL\n(WASP2)']
    gene_counts = [genes_with_snp, genes_with_any]
    colors = [snp_color, indel_color]

    bars = ax1.bar(categories, gene_counts, color=colors, edgecolor='black', linewidth=1.5, width=0.6)

    for bar, val in zip(bars, gene_counts):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                 str(val), ha='center', va='bottom', fontsize=16, fontweight='bold')

    # Improvement annotation
    if genes_with_snp > 0:
        ax1.annotate(f'{improvement:.1f}× more genes',
                     xy=(1, gene_counts[1]),
                     xytext=(1.35, gene_counts[1] * 0.85),
                     fontsize=14, fontweight='bold', color=indel_color,
                     arrowprops=dict(arrowstyle='->', color=indel_color, lw=2))

    ax1.set_ylabel('Imprinted Genes Detected', fontsize=14)
    ax1.set_title('A  How many imprinted genes can we analyze?', fontsize=13, fontweight='bold', loc='left')
    ax1.set_ylim(0, max(gene_counts) * 1.35)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # ----- Panel B: Variant Counts per Gene -----
    ax2 = fig.add_subplot(2, 2, 2)

    # For genes with any variants, show SNP vs INDEL distribution
    genes_with_variants = analyzable[analyzable['has_any']].copy()
    genes_with_variants = genes_with_variants.sort_values('n_total', ascending=True).tail(20)  # Top 20

    y = np.arange(len(genes_with_variants))
    height = 0.7

    ax2.barh(y, genes_with_variants['n_snps'], height, label='SNPs', color=snp_color, edgecolor='black')
    ax2.barh(y, genes_with_variants['n_indels'], height, left=genes_with_variants['n_snps'],
             label='INDELs', color='#d62728', edgecolor='black')

    ax2.set_yticks(y)
    ax2.set_yticklabels(genes_with_variants['gene'], fontsize=9)
    ax2.set_xlabel('Heterozygous Variants in NA12878', fontsize=12)
    ax2.set_title('B  SNP vs INDEL variants per gene (top 20)', fontsize=13, fontweight='bold', loc='left')
    ax2.legend(loc='lower right', fontsize=10)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    # ----- Panel C: INDEL-only Genes -----
    ax3 = fig.add_subplot(2, 2, 3)

    indel_only_genes = analyzable[analyzable['indel_only']]

    if len(indel_only_genes) > 0:
        # Show these genes
        genes_list = indel_only_genes.sort_values('n_indels', ascending=False).head(10)
        y = np.arange(len(genes_list))

        ax3.barh(y, genes_list['n_indels'], height=0.7, color='#d62728', edgecolor='black')
        ax3.set_yticks(y)
        ax3.set_yticklabels(genes_list['gene'], fontsize=10)
        ax3.set_xlabel('INDEL Variants', fontsize=12)
        ax3.set_title(f'C  Genes ONLY accessible via INDELs ({len(indel_only_genes)} total)',
                      fontsize=13, fontweight='bold', loc='left')
    else:
        ax3.text(0.5, 0.5, 'No INDEL-only genes\nin this dataset',
                 ha='center', va='center', fontsize=14)
        ax3.set_title('C  Genes ONLY accessible via INDELs', fontsize=13, fontweight='bold', loc='left')

    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)

    # ----- Panel D: Summary Statistics -----
    ax4 = fig.add_subplot(2, 2, 4)

    # Summary text
    summary_text = f"""
WASP2 INDEL Support Summary
═══════════════════════════

Dataset: NA12878 imprinted genes
Total genes analyzed: {len(analyzable)}

Heterozygous Variant Coverage:
  • SNP-only:      {genes_with_snp} genes ({100*genes_with_snp/len(analyzable):.0f}%)
  • SNP+INDEL:     {genes_with_any} genes ({100*genes_with_any/len(analyzable):.0f}%)
  • INDEL-only:    {genes_indel_only} genes

Key Finding:
  INDEL support enables analysis of
  {improvement:.1f}× more imprinted genes

Biological Implication:
  Genes with only INDEL het variants
  are INVISIBLE to SNP-only pipelines
  but fully accessible with WASP2.
"""

    ax4.text(0.05, 0.95, summary_text, transform=ax4.transAxes,
             fontsize=11, verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    ax4.axis('off')
    ax4.set_title('D  Summary', fontsize=13, fontweight='bold', loc='left')

    plt.tight_layout()
    plt.savefig('results/indel_support_gene_coverage.png', dpi=150, bbox_inches='tight')
    plt.savefig('results/indel_support_gene_coverage.pdf', bbox_inches='tight')

    print("\nSaved: results/indel_support_gene_coverage.png")
    print("Saved: results/indel_support_gene_coverage.pdf")

    return results_df


def main():
    print("=" * 70)
    print("Replicating Aaron's INDEL Support Analysis")
    print("=" * 70)

    # Check if we have bcftools
    try:
        subprocess.run(['bcftools', '--version'], capture_output=True, check=True)
    except:
        print("ERROR: bcftools not found. Please load the module.")
        return

    # Analyze gene coverage
    results_df = analyze_gene_coverage()

    # Save results
    results_df.to_csv('results/gene_variant_coverage.tsv', sep='\t', index=False)
    print("\nSaved: results/gene_variant_coverage.tsv")

    # Create comparison figure
    create_comparison_figure(results_df)

    print("\n" + "=" * 70)
    print("Analysis complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
