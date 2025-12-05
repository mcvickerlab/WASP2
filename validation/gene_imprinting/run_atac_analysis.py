#!/usr/bin/env python3
"""
WASP2-Rust Gene Imprinting Analysis - CORRECT VERSION

Using Aaron's ATAC-seq data and full imprinted gene list to replicate
the "3× more genes" finding with INDEL support.
"""
import sys
import gzip
import subprocess
from pathlib import Path
from datetime import datetime
from collections import defaultdict

sys.path.insert(0, '/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp')

import pandas as pd
import numpy as np

# CORRECT: Aaron's ATAC-seq data
ATAC_BAM = "/iblm/netapp/data1/aho/atac/Buenrostro2013/merged/GM12878_ATACseq_50k/GM12878_ATACseq_50k_merged.sorted.bam"

# Full NA12878 VCF
FULL_VCF = "/iblm/netapp/data1/aho/variants/NA12878.vcf.gz"
SAMPLE = "NA12878"

# Aaron's full imprinted gene list (187 genes)
AARON_IMPRINTED = "/iblm/netapp/home/aho/projects/wasp/testing/performance/data/GM12878_ATACseq_50k_merged/test_imprinted/geneimprint_coords.tsv"

OUTPUT_DIR = Path("results/atac_analysis")


def extract_het_variants_for_genes(genes_df, vcf_path, sample):
    """Extract all heterozygous variants overlapping imprinted genes."""
    print("\nExtracting heterozygous variants from VCF...")

    all_variants = []
    gene_variant_map = defaultdict(list)

    # Get unique gene regions
    unique_regions = genes_df.drop_duplicates(subset=['chrom', 'start', 'end'])[['chrom', 'start', 'end', 'Gene', 'Status']].values

    for i, (chrom, start, end, gene, status) in enumerate(unique_regions):
        if i % 20 == 0:
            print(f"  Processing region {i+1}/{len(unique_regions)}...")

        region = f"{chrom}:{start}-{end}"
        cmd = f"bcftools view -r {region} -s {sample} {vcf_path} 2>/dev/null | bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%GT]\\n'"

        try:
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=60)

            for line in result.stdout.strip().split('\n'):
                if not line:
                    continue
                fields = line.split('\t')
                if len(fields) < 5:
                    continue

                v_chrom, pos, ref, alt, gt = fields

                # Check if heterozygous
                if '/' in gt or '|' in gt:
                    alleles = gt.replace('|', '/').split('/')
                    if len(alleles) == 2 and alleles[0] != alleles[1] and '.' not in alleles:
                        is_indel = len(ref) != len(alt.split(',')[0])

                        variant = {
                            'chrom': v_chrom,
                            'pos': int(pos),
                            'ref': ref,
                            'alt': alt.split(',')[0],  # Take first alt for multiallelic
                            'gt': gt,
                            'is_indel': is_indel,
                            'type': 'INDEL' if is_indel else 'SNP',
                            'gene': gene,
                            'status': status
                        }
                        all_variants.append(variant)
                        gene_variant_map[gene].append(variant)

        except Exception as e:
            print(f"  Warning: Error processing {region}: {e}")
            continue

    return all_variants, gene_variant_map


def run_rust_bamcounter(variants, bam_path):
    """Run WASP2-Rust BamCounter on the variants."""
    from wasp2_rust import BamCounter

    print(f"\nRunning WASP2-Rust BamCounter on {len(variants)} variants...")
    print(f"BAM: {bam_path}")

    # Prepare regions for BamCounter: (chrom, pos, ref, alt)
    regions = [(v['chrom'], v['pos'], v['ref'], v['alt']) for v in variants]

    counter = BamCounter(bam_path)
    counts = counter.count_alleles(regions, min_qual=20, threads=8)

    # Add counts to variants
    for i, (ref_cnt, alt_cnt, other_cnt) in enumerate(counts):
        variants[i]['ref_count'] = ref_cnt
        variants[i]['alt_count'] = alt_cnt
        variants[i]['other_count'] = other_cnt
        variants[i]['total'] = ref_cnt + alt_cnt

    return variants


def analyze_results(variants, gene_variant_map):
    """Analyze SNP-only vs SNP+INDEL coverage."""

    # Convert to DataFrame
    df = pd.DataFrame(variants)

    # Filter to variants with reads
    df_with_reads = df[df['total'] > 0].copy()

    print("\n" + "=" * 70)
    print("VARIANT COVERAGE SUMMARY")
    print("=" * 70)
    print(f"Total het variants: {len(df)}")
    print(f"  SNPs: {len(df[df['type'] == 'SNP'])}")
    print(f"  INDELs: {len(df[df['type'] == 'INDEL'])}")
    print(f"\nVariants with reads: {len(df_with_reads)}")
    print(f"  SNPs: {len(df_with_reads[df_with_reads['type'] == 'SNP'])}")
    print(f"  INDELs: {len(df_with_reads[df_with_reads['type'] == 'INDEL'])}")

    # Gene-level analysis
    print("\n" + "=" * 70)
    print("GENE-LEVEL ANALYSIS (SNP-only vs SNP+INDEL)")
    print("=" * 70)

    gene_stats = []
    for gene in df['gene'].unique():
        gene_df = df[df['gene'] == gene]
        gene_with_reads = gene_df[gene_df['total'] > 0]

        snp_variants = gene_df[gene_df['type'] == 'SNP']
        indel_variants = gene_df[gene_df['type'] == 'INDEL']

        snp_with_reads = gene_with_reads[gene_with_reads['type'] == 'SNP']
        indel_with_reads = gene_with_reads[gene_with_reads['type'] == 'INDEL']

        gene_stats.append({
            'gene': gene,
            'status': gene_df['status'].iloc[0],
            'n_snp_variants': len(snp_variants),
            'n_indel_variants': len(indel_variants),
            'n_snp_with_reads': len(snp_with_reads),
            'n_indel_with_reads': len(indel_with_reads),
            'snp_reads': snp_with_reads['total'].sum(),
            'indel_reads': indel_with_reads['total'].sum(),
            'total_reads': gene_with_reads['total'].sum(),
            'has_snp_coverage': len(snp_with_reads) > 0,
            'has_indel_coverage': len(indel_with_reads) > 0,
            'has_any_coverage': len(gene_with_reads) > 0,
            'indel_only': len(indel_with_reads) > 0 and len(snp_with_reads) == 0
        })

    gene_stats_df = pd.DataFrame(gene_stats)

    # Filter to Imprinted/Predicted genes
    analyzable = gene_stats_df[gene_stats_df['status'].isin(['Imprinted', 'Predicted'])]

    genes_snp_only = analyzable[analyzable['has_snp_coverage']]['gene'].nunique()
    genes_any = analyzable[analyzable['has_any_coverage']]['gene'].nunique()
    genes_indel_only = analyzable[analyzable['indel_only']]['gene'].nunique()

    print(f"\nImprinted/Predicted genes analyzed: {len(analyzable)}")
    print(f"\nGenes with read coverage:")
    print(f"  SNP-only accessible:  {genes_snp_only}")
    print(f"  SNP+INDEL accessible: {genes_any}")
    print(f"  INDEL-only genes:     {genes_indel_only}")

    if genes_snp_only > 0:
        improvement = genes_any / genes_snp_only
        print(f"\n  >>> {improvement:.1f}× MORE GENES with INDEL support! <<<")

    # List INDEL-only genes
    indel_only_genes = analyzable[analyzable['indel_only']]
    if len(indel_only_genes) > 0:
        print(f"\nGenes ONLY accessible via INDELs:")
        for _, row in indel_only_genes.iterrows():
            print(f"  - {row['gene']}: {row['n_indel_with_reads']} INDEL variants, {row['indel_reads']} reads")

    return df, gene_stats_df


def create_plot(df, gene_stats_df):
    """Create the 4-panel comparison figure."""
    import matplotlib.pyplot as plt

    # Filter to analyzable genes
    analyzable = gene_stats_df[gene_stats_df['status'].isin(['Imprinted', 'Predicted'])]

    genes_snp_only = analyzable[analyzable['has_snp_coverage']]['gene'].nunique()
    genes_any = analyzable[analyzable['has_any_coverage']]['gene'].nunique()
    genes_indel_only = analyzable[analyzable['indel_only']]['gene'].nunique()

    improvement = genes_any / genes_snp_only if genes_snp_only > 0 else 1.0

    fig = plt.figure(figsize=(14, 10))

    snp_color = '#3366CC'
    indel_color = '#2ca02c'

    # ----- Panel A: Genes Detected -----
    ax1 = fig.add_subplot(2, 2, 1)

    categories = ['SNP-only\n(WASP1)', 'SNP+INDEL\n(WASP2)']
    gene_counts = [genes_snp_only, genes_any]
    colors = [snp_color, indel_color]

    bars = ax1.bar(categories, gene_counts, color=colors, edgecolor='black', linewidth=1.5, width=0.6)

    for bar, val in zip(bars, gene_counts):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                 str(val), ha='center', va='bottom', fontsize=18, fontweight='bold')

    if improvement > 1:
        ax1.annotate(f'{improvement:.1f}× more genes',
                     xy=(1, gene_counts[1]),
                     xytext=(1.35, gene_counts[1] * 0.8),
                     fontsize=14, fontweight='bold', color=indel_color,
                     arrowprops=dict(arrowstyle='->', color=indel_color, lw=2))

    ax1.set_ylabel('Imprinted Genes with Coverage', fontsize=14)
    ax1.set_title('A  How many imprinted genes do we detect?', fontsize=13, fontweight='bold', loc='left')
    ax1.set_ylim(0, max(gene_counts) * 1.4)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # ----- Panel B: Significance per gene (reads as proxy) -----
    ax2 = fig.add_subplot(2, 2, 2)

    # Top genes by read count
    top_genes = analyzable[analyzable['has_any_coverage']].nlargest(15, 'total_reads')

    y = np.arange(len(top_genes))
    height = 0.7

    ax2.barh(y, top_genes['snp_reads'], height, label='SNP reads', color=snp_color, edgecolor='black')
    ax2.barh(y, top_genes['indel_reads'], height, left=top_genes['snp_reads'],
             label='INDEL reads', color='#d62728', edgecolor='black')

    ax2.set_yticks(y)
    ax2.set_yticklabels(top_genes['gene'], fontsize=10)
    ax2.set_xlabel('Total Reads', fontsize=12)
    ax2.set_title('B  Read coverage per gene (top 15)', fontsize=13, fontweight='bold', loc='left')
    ax2.legend(loc='lower right', fontsize=10)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    # ----- Panel C: INDEL-only Genes -----
    ax3 = fig.add_subplot(2, 2, 3)

    indel_only = analyzable[analyzable['indel_only']].nlargest(10, 'indel_reads')

    if len(indel_only) > 0:
        y = np.arange(len(indel_only))
        ax3.barh(y, indel_only['indel_reads'], height=0.7, color='#d62728', edgecolor='black')
        ax3.set_yticks(y)
        ax3.set_yticklabels(indel_only['gene'], fontsize=11)
        ax3.set_xlabel('INDEL Reads', fontsize=12)
        ax3.set_title(f'C  Genes ONLY accessible via INDELs ({genes_indel_only} genes)',
                      fontsize=13, fontweight='bold', loc='left')
    else:
        ax3.text(0.5, 0.5, 'No INDEL-only genes\nwith coverage', ha='center', va='center', fontsize=14)
        ax3.set_title('C  Genes ONLY accessible via INDELs', fontsize=13, fontweight='bold', loc='left')

    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)

    # ----- Panel D: Read contribution -----
    ax4 = fig.add_subplot(2, 2, 4)

    total_snp_reads = analyzable['snp_reads'].sum()
    total_indel_reads = analyzable['indel_reads'].sum()
    total_reads = total_snp_reads + total_indel_reads

    x = np.arange(2)
    ax4.bar(x, [total_snp_reads, total_snp_reads], width=0.6, label='SNP reads', color=snp_color, edgecolor='black')
    ax4.bar(x[1], total_indel_reads, width=0.6, bottom=total_snp_reads, label='INDEL reads', color='#d62728', edgecolor='black')

    ax4.text(0, total_snp_reads + 50, f'{total_snp_reads:,}', ha='center', va='bottom', fontsize=12, fontweight='bold')
    ax4.text(1, total_reads + 50, f'{total_reads:,}', ha='center', va='bottom', fontsize=12, fontweight='bold')

    if total_reads > 0:
        pct = 100 * total_indel_reads / total_reads
        ax4.text(1, total_snp_reads + total_indel_reads/2, f'+{pct:.0f}%\nINDEL',
                 ha='center', va='center', fontsize=10, fontweight='bold', color='white')

    ax4.set_xticks(x)
    ax4.set_xticklabels(['SNP-only\n(WASP1)', 'SNP+INDEL\n(WASP2)'])
    ax4.set_ylabel('Total Reads', fontsize=14)
    ax4.set_title('D  Total read contribution', fontsize=13, fontweight='bold', loc='left')
    ax4.legend(loc='upper left', fontsize=10)
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)

    plt.suptitle('WASP2 INDEL Support Improves Imprinted Gene Detection\n(GM12878 ATAC-seq, validated Aho dataset)',
                 fontsize=15, fontweight='bold', y=1.02)

    plt.tight_layout()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    plt.savefig(OUTPUT_DIR / 'indel_support_imprinting.png', dpi=150, bbox_inches='tight')
    plt.savefig(OUTPUT_DIR / 'indel_support_imprinting.pdf', bbox_inches='tight')

    print(f"\nSaved: {OUTPUT_DIR / 'indel_support_imprinting.png'}")
    print(f"Saved: {OUTPUT_DIR / 'indel_support_imprinting.pdf'}")


def main():
    print("=" * 70)
    print("WASP2-Rust Gene Imprinting Analysis - CORRECT VERSION")
    print(f"Time: {datetime.now()}")
    print("=" * 70)
    print(f"\nATAC-seq BAM: {ATAC_BAM}")
    print(f"VCF: {FULL_VCF}")
    print(f"Imprinted genes: {AARON_IMPRINTED}")

    # Load imprinted genes
    genes_df = pd.read_csv(AARON_IMPRINTED, sep='\t')
    print(f"\nLoaded {len(genes_df)} imprinted gene entries")
    print(f"Unique genes: {genes_df['Gene'].nunique()}")

    # Extract het variants
    variants, gene_map = extract_het_variants_for_genes(genes_df, FULL_VCF, SAMPLE)
    print(f"\nTotal het variants extracted: {len(variants)}")

    # Run BamCounter
    variants = run_rust_bamcounter(variants, ATAC_BAM)

    # Analyze results
    df, gene_stats = analyze_results(variants, gene_map)

    # Save data
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUTPUT_DIR / 'variant_counts.tsv', sep='\t', index=False)
    gene_stats.to_csv(OUTPUT_DIR / 'gene_stats.tsv', sep='\t', index=False)
    print(f"\nSaved: {OUTPUT_DIR / 'variant_counts.tsv'}")
    print(f"Saved: {OUTPUT_DIR / 'gene_stats.tsv'}")

    # Create plot
    create_plot(df, gene_stats)

    print("\n" + "=" * 70)
    print(f"COMPLETED: {datetime.now()}")
    print("=" * 70)


if __name__ == "__main__":
    main()
