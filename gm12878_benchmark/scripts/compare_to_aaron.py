#!/usr/bin/env python3
"""
Compare our WASP2 indel results to Aaron's SNP-only baseline
"""

import pandas as pd
import sys
from pathlib import Path

def load_counts(filepath):
    """Load variant counts file"""
    print(f"Loading: {filepath}")
    df = pd.read_csv(filepath, sep="\t")
    print(f"  Loaded {len(df):,} variants")
    return df

def analyze_variant_types(df):
    """Determine SNPs vs indels from ref/alt alleles"""
    # Extract ref and alt allele columns
    # Format: chr1:12345 A>T
    df['variant_info'] = df['variant'].str.split(' ', n=1).str[1]
    df['ref_allele'] = df['variant_info'].str.split('>').str[0]
    df['alt_allele'] = df['variant_info'].str.split('>').str[1]

    df['ref_len'] = df['ref_allele'].str.len()
    df['alt_len'] = df['alt_allele'].str.len()

    # Classify variants
    df['variant_type'] = 'SNP'
    df.loc[df['ref_len'] != df['alt_len'], 'variant_type'] = 'indel'
    df.loc[(df['ref_len'] > 1) & (df['alt_len'] > 1) & (df['ref_len'] == df['alt_len']), 'variant_type'] = 'MNP'

    return df

def main():
    print("="*60)
    print("WASP2 Indel Validation: Comparison to Aaron's Baseline")
    print("="*60)
    print()

    # File paths
    aaron_file = Path("/iblm/netapp/data3/aho/project_data/wasp2/imprinted_rna/outputs/GM12878_geneimprint_transcript_gene_counts.tsv")

    our_counting_only = Path("gm12878_benchmark/results/counts_WITH_INDELS.tsv")
    our_wasp_filtered = Path("gm12878_benchmark/test_full_pipeline/counts_WITH_INDELS_WASP_FILTERED.tsv")

    # Load Aaron's baseline
    print("1. Aaron's SNP-only baseline (no WASP remapping):")
    print("-" * 60)
    if aaron_file.exists():
        aaron_df = load_counts(aaron_file)
        aaron_df = analyze_variant_types(aaron_df)

        print(f"\n  Variant types:")
        print(aaron_df['variant_type'].value_counts())
        print(f"\n  Total genes: {aaron_df['gene_name'].nunique()}")
        print(f"  Method: Counting only (no WASP remapping)")
    else:
        print(f"  ❌ File not found: {aaron_file}")
        aaron_df = None
    print()

    # Load our counting-only results
    print("2. Our SNP+indel counting (no WASP remapping):")
    print("-" * 60)
    if our_counting_only.exists():
        counting_df = load_counts(our_counting_only)
        counting_df = analyze_variant_types(counting_df)

        print(f"\n  Variant types:")
        print(counting_df['variant_type'].value_counts())
        print(f"\n  Total genes: {counting_df['gene_name'].nunique()}")
        print(f"  Method: Counting only (no WASP remapping)")

        if aaron_df is not None:
            snp_only_count = (counting_df['variant_type'] == 'SNP').sum()
            indel_count = (counting_df['variant_type'] == 'indel').sum()
            aaron_var_count = len(aaron_df)

            print(f"\n  Comparison to Aaron:")
            print(f"    Aaron SNPs:     {aaron_var_count:,}")
            print(f"    Our SNPs:       {snp_only_count:,} ({'same' if snp_only_count == aaron_var_count else 'DIFFERENT!'})")
            print(f"    Our indels:     {indel_count:,} (NEW!)")
            print(f"    Total increase: +{indel_count:,} variants ({100*indel_count/aaron_var_count:.1f}%)")
    else:
        print(f"  ⚠️  File not found: {our_counting_only}")
        print(f"     Run: bash gm12878_benchmark/scripts/run_benchmark.sh")
        counting_df = None
    print()

    # Load our WASP-filtered results
    print("3. Our SNP+indel with WASP remapping (bias removal):")
    print("-" * 60)
    if our_wasp_filtered.exists():
        wasp_df = load_counts(our_wasp_filtered)
        wasp_df = analyze_variant_types(wasp_df)

        print(f"\n  Variant types:")
        print(wasp_df['variant_type'].value_counts())
        print(f"\n  Total genes: {wasp_df['gene_name'].nunique()}")
        print(f"  Method: Full WASP2 pipeline with indel support")

        if counting_df is not None:
            wasp_snps = (wasp_df['variant_type'] == 'SNP').sum()
            wasp_indels = (wasp_df['variant_type'] == 'indel').sum()
            count_snps = (counting_df['variant_type'] == 'SNP').sum()
            count_indels = (counting_df['variant_type'] == 'indel').sum()

            print(f"\n  Effect of WASP filtering:")
            print(f"    SNPs before:    {count_snps:,}")
            print(f"    SNPs after:     {wasp_snps:,} (-{count_snps - wasp_snps:,}, {100*(count_snps - wasp_snps)/count_snps:.1f}% filtered)")
            print(f"    Indels before:  {count_indels:,}")
            print(f"    Indels after:   {wasp_indels:,} (-{count_indels - wasp_indels:,}, {100*(count_indels - wasp_indels)/count_indels:.1f}% filtered)")
    else:
        print(f"  ⚠️  File not found: {our_wasp_filtered}")
        print(f"     Run: bash gm12878_benchmark/scripts/run_full_pipeline_tier1.sh")
    print()

    # Check for classic imprinted genes
    print("4. Classic imprinted gene validation:")
    print("-" * 60)
    classic_genes = ['H19', 'IGF2', 'SNRPN', 'XIST', 'MEST', 'PEG3']

    for gene in classic_genes:
        print(f"\n  {gene}:")
        if aaron_df is not None:
            aaron_count = len(aaron_df[aaron_df['gene_name'] == gene])
            print(f"    Aaron:  {aaron_count:3d} variants")

        if counting_df is not None:
            our_count = len(counting_df[counting_df['gene_name'] == gene])
            our_indels = len(counting_df[(counting_df['gene_name'] == gene) & (counting_df['variant_type'] == 'indel')])
            print(f"    Ours:   {our_count:3d} variants ({our_indels} indels)")

    print()
    print("="*60)
    print("Analysis complete!")
    print("="*60)

if __name__ == "__main__":
    main()
