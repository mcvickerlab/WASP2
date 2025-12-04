#!/usr/bin/env python3
"""
Create ground truth CSV from VCF file for benchmark comparison.

For heterozygous variants (0|1), assumes true_ref_ratio = 0.5 (equal expression).
For homozygous ref (0|0), assumes true_ref_ratio = 1.0.
For homozygous alt (1|1), assumes true_ref_ratio = 0.0.

Usage:
    python simulation/create_ground_truth_from_vcf.py variants.vcf.gz output.csv
"""

import sys
import gzip
import pandas as pd
from pathlib import Path


def parse_vcf_to_ground_truth(vcf_file: str, output_csv: str):
    """
    Parse VCF and create ground truth CSV.

    Args:
        vcf_file: Path to VCF file (can be gzipped)
        output_csv: Path to output CSV
    """
    records = []

    # Open VCF (handle both gzipped and plain)
    if vcf_file.endswith('.gz'):
        f = gzip.open(vcf_file, 'rt')
    else:
        f = open(vcf_file, 'r')

    for line in f:
        # Skip header lines
        if line.startswith('#'):
            continue

        fields = line.strip().split('\t')
        chrom = fields[0]
        pos = int(fields[1])
        variant_id = fields[2]
        ref = fields[3]
        alt = fields[4]
        genotype = fields[9].split(':')[0]  # GT field

        # Determine true ref ratio based on genotype
        if genotype == '0|1' or genotype == '1|0' or genotype == '0/1':
            # Heterozygous: 50:50 ratio
            true_ref_ratio = 0.5
        elif genotype == '0|0' or genotype == '0/0':
            # Homozygous ref: all ref
            true_ref_ratio = 1.0
        elif genotype == '1|1' or genotype == '1/1':
            # Homozygous alt: all alt
            true_ref_ratio = 0.0
        else:
            print(f"Warning: Unknown genotype {genotype} at {chrom}:{pos}, skipping", file=sys.stderr)
            continue

        records.append({
            'chrom': chrom,
            'pos': pos,
            'variant_id': variant_id,
            'ref_allele': ref,
            'alt_allele': alt,
            'genotype': genotype,
            'true_ref_ratio': true_ref_ratio
        })

    f.close()

    # Create DataFrame and save
    df = pd.DataFrame(records)
    df.to_csv(output_csv, index=False)

    print(f"Created ground truth CSV with {len(df)} variants")
    print(f"  Heterozygous sites: {(df['true_ref_ratio'] == 0.5).sum()}")
    print(f"  Homozygous ref: {(df['true_ref_ratio'] == 1.0).sum()}")
    print(f"  Homozygous alt: {(df['true_ref_ratio'] == 0.0).sum()}")
    print(f"Saved to: {output_csv}")


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python create_ground_truth_from_vcf.py <vcf_file> <output_csv>")
        sys.exit(1)

    vcf_file = sys.argv[1]
    output_csv = sys.argv[2]

    if not Path(vcf_file).exists():
        print(f"ERROR: VCF file not found: {vcf_file}", file=sys.stderr)
        sys.exit(1)

    parse_vcf_to_ground_truth(vcf_file, output_csv)
