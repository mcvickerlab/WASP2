#!/usr/bin/env python3
"""
Full WASP2-Rust Pipeline for Gene Imprinting Validation

Proper WASP pipeline workflow:
1. filter_bam_by_variants - Split BAM into to_remap + keep
2. unified_make_reads - Generate swapped-allele FASTQs from to_remap.bam
3. BWA remap - Remap swapped reads
4. filter_bam_wasp - Keep reads that map back to same position
5. Merge keep + filtered_remap → wasp_filtered.bam
6. Count alleles and analyze significance
"""
import sys
import os
import subprocess
from pathlib import Path
from datetime import datetime

sys.path.insert(0, '/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp')

import wasp2_rust
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

# ============ PATHS ============
RAW_BAM = "/iblm/netapp/data1/aho/atac/Buenrostro2013/merged/GM12878_ATACseq_50k/GM12878_ATACseq_50k_merged.sorted.bam"
FULL_VCF = "/iblm/netapp/data1/aho/variants/NA12878.vcf.gz"
SAMPLE = "NA12878"
REF_GENOME = "/iblm/netapp/data1/aho/ref_genomes/index/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa"
BWA = "/iblm/netapp/oldhome/j3gu/anaconda2/bin/bwa"
SAMTOOLS = "/iblm/netapp/home/jjaureguy/mambaforge/bin/samtools"
AARON_IMPRINTED = "/iblm/netapp/home/aho/projects/wasp/testing/performance/data/GM12878_ATACseq_50k_merged/test_imprinted/geneimprint_coords.tsv"

OUTPUT_DIR = Path("results/rust_pipeline")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def create_variant_bed():
    """Extract het variants from VCF and create BED."""
    print("\n" + "=" * 70)
    print("STEP 1: Creating variant BED from VCF")
    print("=" * 70)

    genes_df = pd.read_csv(AARON_IMPRINTED, sep='\t')
    print(f"Loaded {len(genes_df)} imprinted gene entries")

    bed_path = OUTPUT_DIR / "het_variants.bed"
    variants = []

    unique_regions = genes_df.drop_duplicates(subset=['chrom', 'start', 'end'])[['chrom', 'start', 'end', 'Gene', 'Status']].values
    print(f"Extracting het variants from {len(unique_regions)} gene regions...")

    for i, (chrom, start, end, gene, status) in enumerate(unique_regions):
        if i % 30 == 0:
            print(f"  Processing region {i+1}/{len(unique_regions)}...")

        region = f"{chrom}:{start}-{end}"
        cmd = f"bcftools view -r {region} -s {SAMPLE} {FULL_VCF} 2>/dev/null | bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%GT]\\n'"

        try:
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=60)
            for line in result.stdout.strip().split('\n'):
                if not line:
                    continue
                fields = line.split('\t')
                if len(fields) < 5:
                    continue
                v_chrom, pos, ref, alt, gt = fields
                pos = int(pos)
                if '/' in gt or '|' in gt:
                    alleles = gt.replace('|', '/').split('/')
                    if len(alleles) == 2 and alleles[0] != alleles[1] and '.' not in alleles:
                        alt = alt.split(',')[0]
                        is_indel = len(ref) != len(alt)
                        bed_start = pos - 1
                        bed_stop = bed_start + len(ref)
                        variants.append({
                            'chrom': v_chrom, 'start': bed_start, 'stop': bed_stop,
                            'ref': ref, 'alt': alt, 'gt': gt, 'pos': pos,
                            'gene': gene, 'status': status,
                            'type': 'INDEL' if is_indel else 'SNP'
                        })
        except Exception as e:
            print(f"  Warning: Error processing {region}: {e}")

    with open(bed_path, 'w') as f:
        for v in variants:
            f.write(f"{v['chrom']}\t{v['start']}\t{v['stop']}\t{v['ref']}\t{v['alt']}\t{v['gt']}\n")

    var_df = pd.DataFrame(variants)
    var_df.to_csv(OUTPUT_DIR / "variant_info.tsv", sep='\t', index=False)

    n_snps = sum(1 for v in variants if v['type'] == 'SNP')
    n_indels = sum(1 for v in variants if v['type'] == 'INDEL')
    print(f"\nCreated {bed_path}")
    print(f"Total het variants: {len(variants)} ({n_snps} SNPs, {n_indels} INDELs)")

    return bed_path, var_df


def split_bam_by_variants(bed_path):
    """Split BAM into reads overlapping variants (to_remap) and others (keep)."""
    print("\n" + "=" * 70)
    print("STEP 2: Splitting BAM by variant overlap")
    print("=" * 70)

    to_remap_bam = OUTPUT_DIR / "to_remap.bam"
    keep_bam = OUTPUT_DIR / "keep.bam"

    print(f"Input BAM: {RAW_BAM}")
    print(f"Variant BED: {bed_path}")

    remap_count, keep_count, unique_names = wasp2_rust.filter_bam_by_variants_py(
        str(RAW_BAM),
        str(bed_path),
        str(to_remap_bam),
        str(keep_bam),
        is_paired=True,
        threads=8
    )

    print(f"\nSplit results:")
    print(f"  Reads to remap: {remap_count}")
    print(f"  Reads to keep: {keep_count}")
    print(f"  Unique read names with variants: {unique_names}")

    # Index to_remap.bam
    subprocess.run(f"{SAMTOOLS} index {to_remap_bam}", shell=True)

    return to_remap_bam, keep_bam


def run_unified_pipeline(to_remap_bam, bed_path):
    """Generate swapped-allele reads using Rust unified pipeline."""
    print("\n" + "=" * 70)
    print("STEP 3: Running unified pipeline (generate swapped reads)")
    print("=" * 70)

    out_r1 = OUTPUT_DIR / "remap_r1.fq.gz"
    out_r2 = OUTPUT_DIR / "remap_r2.fq.gz"

    print(f"Input BAM: {to_remap_bam}")
    print(f"Output: {out_r1}, {out_r2}")

    stats = wasp2_rust.unified_make_reads_parallel_py(
        str(to_remap_bam),
        str(bed_path),
        str(out_r1),
        str(out_r2),
        max_seqs=64,
        threads=8
    )

    print(f"\nUnified pipeline stats:")
    for key, value in stats.items():
        print(f"  {key}: {value}")

    return out_r1, out_r2


def remap_with_bwa(r1_path, r2_path):
    """Remap swapped-allele reads with BWA."""
    print("\n" + "=" * 70)
    print("STEP 4: Remapping with BWA")
    print("=" * 70)

    sorted_bam = OUTPUT_DIR / "remapped.sorted.bam"
    print(f"Reference: {REF_GENOME}")

    bwa_cmd = f"{BWA} mem -t 8 {REF_GENOME} {r1_path} {r2_path} 2>/dev/null | {SAMTOOLS} sort -@ 4 -o {sorted_bam}"
    print(f"Running BWA...")

    result = subprocess.run(bwa_cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"BWA error: {result.stderr}")
        raise RuntimeError("BWA remapping failed")

    subprocess.run(f"{SAMTOOLS} index {sorted_bam}", shell=True)
    print(f"Created: {sorted_bam}")

    return sorted_bam


def wasp_filter(to_remap_bam, remapped_bam):
    """Filter remapped reads using WASP criteria."""
    print("\n" + "=" * 70)
    print("STEP 5: WASP filtering")
    print("=" * 70)

    remap_keep_bam = OUTPUT_DIR / "remap_keep.bam"

    print(f"Original reads: {to_remap_bam}")
    print(f"Remapped reads: {remapped_bam}")

    kept, removed_moved, removed_missing = wasp2_rust.filter_bam_wasp(
        str(to_remap_bam),
        str(remapped_bam),
        str(remap_keep_bam),
        threads=8
    )

    print(f"\nWASP filter results:")
    print(f"  Kept reads: {kept}")
    print(f"  Removed (moved): {removed_moved}")
    print(f"  Removed (missing): {removed_missing}")

    return remap_keep_bam


def merge_and_finalize(keep_bam, remap_keep_bam):
    """Merge keep + filtered_remap into final WASP-filtered BAM."""
    print("\n" + "=" * 70)
    print("STEP 6: Merging final WASP-filtered BAM")
    print("=" * 70)

    final_bam = OUTPUT_DIR / "wasp_filtered.bam"

    # Merge the two BAMs
    merge_cmd = f"{SAMTOOLS} merge -f {final_bam} {keep_bam} {remap_keep_bam}"
    subprocess.run(merge_cmd, shell=True)

    # Sort and index
    sorted_final = OUTPUT_DIR / "wasp_filtered.sorted.bam"
    subprocess.run(f"{SAMTOOLS} sort -@ 4 -o {sorted_final} {final_bam}", shell=True)
    subprocess.run(f"{SAMTOOLS} index {sorted_final}", shell=True)

    print(f"Created: {sorted_final}")

    return sorted_final


def count_alleles(filtered_bam, var_df):
    """Count alleles at each variant position."""
    print("\n" + "=" * 70)
    print("STEP 7: Counting alleles")
    print("=" * 70)

    regions = [(row['chrom'], row['pos'], row['ref'], row['alt']) for _, row in var_df.iterrows()]
    print(f"Counting at {len(regions)} variant positions...")

    counter = wasp2_rust.BamCounter(str(filtered_bam))
    counts = counter.count_alleles(regions, min_qual=20, threads=8)

    var_df = var_df.copy()
    var_df['ref_count'] = [c[0] for c in counts]
    var_df['alt_count'] = [c[1] for c in counts]
    var_df['total'] = var_df['ref_count'] + var_df['alt_count']

    var_df.to_csv(OUTPUT_DIR / "allele_counts.tsv", sep='\t', index=False)

    n_with_reads = len(var_df[var_df['total'] > 0])
    print(f"Variants with reads: {n_with_reads}/{len(var_df)}")
    print(f"Total allelic reads: {var_df['total'].sum()}")

    return var_df


def run_significance_analysis(var_df, fdr_threshold=0.1, min_reads=10):
    """Calculate gene-level significance."""
    print("\n" + "=" * 70)
    print("STEP 8: Significance Analysis")
    print("=" * 70)

    results = []
    for gene in var_df['gene'].unique():
        gene_df = var_df[var_df['gene'] == gene]
        gene_with_reads = gene_df[gene_df['total'] > 0]

        snp_df = gene_with_reads[gene_with_reads['type'] == 'SNP']
        snp_ref, snp_alt = snp_df['ref_count'].sum(), snp_df['alt_count'].sum()
        snp_total = snp_ref + snp_alt

        all_ref = gene_with_reads['ref_count'].sum()
        all_alt = gene_with_reads['alt_count'].sum()
        all_total = all_ref + all_alt

        indel_total = gene_with_reads[gene_with_reads['type'] == 'INDEL']['total'].sum()

        snp_pval = stats.binomtest(int(snp_ref), int(snp_total), 0.5).pvalue if snp_total >= min_reads else 1.0
        all_pval = stats.binomtest(int(all_ref), int(all_total), 0.5).pvalue if all_total >= min_reads else 1.0

        results.append({
            'gene': gene, 'status': gene_df['status'].iloc[0],
            'snp_total': snp_total, 'snp_pval': snp_pval,
            'all_total': all_total, 'all_pval': all_pval,
            'indel_total': indel_total
        })

    results_df = pd.DataFrame(results)

    # FDR correction
    for col in ['snp', 'all']:
        valid = results_df[f'{col}_pval'] < 1.0
        if valid.sum() > 0:
            _, fdr, _, _ = multipletests(results_df.loc[valid, f'{col}_pval'], method='fdr_bh')
            results_df.loc[valid, f'{col}_fdr'] = fdr
        results_df[f'{col}_fdr'] = results_df.get(f'{col}_fdr', 1.0).fillna(1.0)

    results_df.to_csv(OUTPUT_DIR / "gene_significance.tsv", sep='\t', index=False)

    analyzable = results_df[results_df['status'].isin(['Imprinted', 'Predicted'])]
    print(f"\nAnalyzable genes: {len(analyzable)}")

    for fdr in [0.05, 0.1]:
        snp_sig = analyzable[analyzable['snp_fdr'] < fdr]['gene'].tolist()
        all_sig = analyzable[analyzable['all_fdr'] < fdr]['gene'].tolist()
        ratio = len(all_sig) / len(snp_sig) if len(snp_sig) > 0 else float('inf')

        print(f"\nFDR < {fdr}:")
        print(f"  SNP-only: {len(snp_sig)} genes")
        print(f"  SNP+INDEL: {len(all_sig)} genes")
        print(f"  Ratio: {ratio:.1f}×")

        if fdr == 0.1:
            newly_sig = set(all_sig) - set(snp_sig)
            print(f"  SNP-only genes: {snp_sig}")
            print(f"  SNP+INDEL genes: {all_sig}")
            if newly_sig:
                print(f"  NEWLY significant: {list(newly_sig)}")

    return results_df


def main():
    print("=" * 70)
    print("FULL WASP2-RUST PIPELINE - Gene Imprinting")
    print(f"Started: {datetime.now()}")
    print("=" * 70)

    bed_path, var_df = create_variant_bed()
    to_remap_bam, keep_bam = split_bam_by_variants(bed_path)
    r1_path, r2_path = run_unified_pipeline(to_remap_bam, bed_path)
    remapped_bam = remap_with_bwa(r1_path, r2_path)
    remap_keep_bam = wasp_filter(to_remap_bam, remapped_bam)
    final_bam = merge_and_finalize(keep_bam, remap_keep_bam)
    var_df = count_alleles(final_bam, var_df)
    results_df = run_significance_analysis(var_df)

    print("\n" + "=" * 70)
    print(f"COMPLETED: {datetime.now()}")
    print(f"Results: {OUTPUT_DIR}")
    print("=" * 70)


if __name__ == "__main__":
    main()
