#!/usr/bin/env python3
"""
Generate count comparison data for Figure 2 Panel B

Compares allele counts from:
- WASP2-Rust
- GATK ASEReadCounter
- phASER

For both HG00731 RNA-seq and GM12878 ATAC-seq datasets.
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import argparse


def parse_wasp2_counts(file_path):
    """
    Parse WASP2-Rust count output.

    Expected format (TSV):
    chrom  pos  ref  alt  ref_count  alt_count
    """
    df = pd.read_csv(file_path, sep='\t')

    # Standardize column names
    df = df.rename(columns={
        'chrom': 'chrom',
        'pos': 'pos',
        'ref': 'ref_allele',
        'alt': 'alt_allele',
        'ref_count': 'ref',
        'alt_count': 'alt'
    })

    # Create unique variant ID
    df['variant_id'] = df['chrom'].astype(str) + ':' + df['pos'].astype(str)

    return df[['variant_id', 'chrom', 'pos', 'ref_allele', 'alt_allele', 'ref', 'alt']]


def parse_gatk_counts(file_path):
    """
    Parse GATK ASEReadCounter output.

    Expected format (.table):
    contig  position  variantID  refAllele  altAllele  refCount  altCount  totalCount  ...
    """
    # Read GATK table (skip comment lines)
    with open(file_path, 'r') as f:
        lines = [line for line in f if not line.startswith('#')]

    # Parse as TSV
    from io import StringIO
    df = pd.read_csv(StringIO(''.join(lines)), sep='\t')

    # Standardize column names
    df = df.rename(columns={
        'contig': 'chrom',
        'position': 'pos',
        'refAllele': 'ref_allele',
        'altAllele': 'alt_allele',
        'refCount': 'ref',
        'altCount': 'alt'
    })

    # Create unique variant ID
    df['variant_id'] = df['chrom'].astype(str) + ':' + df['pos'].astype(str)

    return df[['variant_id', 'chrom', 'pos', 'ref_allele', 'alt_allele', 'ref', 'alt']]


def parse_phaser_counts(file_path):
    """
    Parse phASER per-variant allelic counts.

    Prefer `*.allelic_counts.txt` (per-variant ref/alt counts). If a
    `*.haplotypic_counts.txt` path is provided, this will automatically try to
    locate the corresponding `*.allelic_counts.txt` next to it.
    """
    # Prefer allelic_counts (per-variant table).
    file_path = Path(file_path)
    if file_path.name.endswith(".haplotypic_counts.txt"):
        base_path = file_path.parent / file_path.stem.replace(".haplotypic_counts", "")
        candidate = Path(str(base_path) + ".allelic_counts.txt")
        if candidate.exists():
            file_path = candidate

    try:
        df = pd.read_csv(file_path, sep="\t")
    except Exception as e:
        print(f"Error parsing phASER output: {e}")
        return pd.DataFrame(columns=["variant_id", "chrom", "pos", "ref", "alt"])

    # phASER allelic_counts columns typically mirror GATK's.
    if "contig" not in df.columns:
        print(f"Warning: Could not parse phASER output (missing contig): {file_path}")
        return pd.DataFrame(columns=["variant_id", "chrom", "pos", "ref", "alt"])

    pos_col = "position" if "position" in df.columns else ("pos" if "pos" in df.columns else None)
    if pos_col is None:
        print(f"Warning: Could not parse phASER output (missing position): {file_path}")
        return pd.DataFrame(columns=["variant_id", "chrom", "pos", "ref", "alt"])

    if "refCount" in df.columns and "altCount" in df.columns:
        ref_col, alt_col = "refCount", "altCount"
    elif "ref_count" in df.columns and "alt_count" in df.columns:
        ref_col, alt_col = "ref_count", "alt_count"
    else:
        print(f"Warning: Could not parse phASER output (missing ref/alt counts): {file_path}")
        return pd.DataFrame(columns=["variant_id", "chrom", "pos", "ref", "alt"])

    out = pd.DataFrame(
        {
            "chrom": df["contig"].astype(str),
            "pos": df[pos_col].astype(int),
            "ref": df[ref_col].astype(int),
            "alt": df[alt_col].astype(int),
        }
    )
    out["variant_id"] = out["chrom"].astype(str) + ":" + out["pos"].astype(str)
    return out[["variant_id", "chrom", "pos", "ref", "alt"]]


def merge_counts(wasp2_df, gatk_df, phaser_df):
    """Merge counts from all three tools."""

    # Start with WASP2 as base
    merged = wasp2_df[['variant_id', 'chrom', 'pos']].copy()

    # Add WASP2 counts
    merged['wasp2_ref'] = wasp2_df['ref']
    merged['wasp2_alt'] = wasp2_df['alt']

    # Merge GATK counts
    gatk_subset = gatk_df[['variant_id', 'ref', 'alt']].rename(
        columns={'ref': 'gatk_ref', 'alt': 'gatk_alt'}
    )
    merged = merged.merge(gatk_subset, on='variant_id', how='left')

    # Merge phASER counts
    if len(phaser_df) > 0:
        phaser_subset = phaser_df[['variant_id', 'ref', 'alt']].rename(
            columns={'ref': 'phaser_ref', 'alt': 'phaser_alt'}
        )
        merged = merged.merge(phaser_subset, on='variant_id', how='left')
    else:
        merged['phaser_ref'] = np.nan
        merged['phaser_alt'] = np.nan

    # Fill NaN with 0
    merged = merged.fillna(0)

    return merged


def calculate_correlations(df):
    """Calculate Pearson correlations between tools."""

    results = {}

    # Filter to variants with counts in both tools
    mask_wasp2_gatk = (df['wasp2_ref'] + df['wasp2_alt'] > 0) & \
                      (df['gatk_ref'] + df['gatk_alt'] > 0)

    if mask_wasp2_gatk.sum() > 0:
        # WASP2 vs GATK - ref counts
        r_ref, p_ref = pearsonr(
            df.loc[mask_wasp2_gatk, 'wasp2_ref'],
            df.loc[mask_wasp2_gatk, 'gatk_ref']
        )

        # WASP2 vs GATK - alt counts
        r_alt, p_alt = pearsonr(
            df.loc[mask_wasp2_gatk, 'wasp2_alt'],
            df.loc[mask_wasp2_gatk, 'gatk_alt']
        )

        results['wasp2_vs_gatk'] = {
            'ref_r': r_ref,
            'ref_p': p_ref,
            'alt_r': r_alt,
            'alt_p': p_alt,
            'n_variants': mask_wasp2_gatk.sum()
        }

    # WASP2 vs phASER
    mask_wasp2_phaser = (df['wasp2_ref'] + df['wasp2_alt'] > 0) & \
                        (df['phaser_ref'] + df['phaser_alt'] > 0)

    if mask_wasp2_phaser.sum() > 0:
        r_ref, p_ref = pearsonr(
            df.loc[mask_wasp2_phaser, 'wasp2_ref'],
            df.loc[mask_wasp2_phaser, 'phaser_ref']
        )

        r_alt, p_alt = pearsonr(
            df.loc[mask_wasp2_phaser, 'wasp2_alt'],
            df.loc[mask_wasp2_phaser, 'phaser_alt']
        )

        results['wasp2_vs_phaser'] = {
            'ref_r': r_ref,
            'ref_p': p_ref,
            'alt_r': r_alt,
            'alt_p': p_alt,
            'n_variants': mask_wasp2_phaser.sum()
        }

    return results


def process_dataset(dataset_name, wasp2_file, gatk_file, phaser_file, output_file):
    """Process a single dataset."""

    print(f"\n{'='*60}")
    print(f"Processing: {dataset_name}")
    print(f"{'='*60}")

    # Parse counts
    print("Parsing WASP2 counts...")
    wasp2_df = parse_wasp2_counts(wasp2_file)
    print(f"  WASP2: {len(wasp2_df)} variants")

    print("Parsing GATK counts...")
    gatk_df = parse_gatk_counts(gatk_file)
    print(f"  GATK: {len(gatk_df)} variants")

    print("Parsing phASER counts...")
    phaser_df = parse_phaser_counts(phaser_file)
    print(f"  phASER: {len(phaser_df)} variants")

    # Merge counts
    print("\nMerging counts...")
    merged_df = merge_counts(wasp2_df, gatk_df, phaser_df)
    print(f"  Merged: {len(merged_df)} variants")

    # Calculate correlations
    print("\nCalculating correlations...")
    correlations = calculate_correlations(merged_df)

    for comparison, stats in correlations.items():
        print(f"\n  {comparison}:")
        print(f"    Ref counts: r={stats['ref_r']:.4f}, p={stats['ref_p']:.2e}")
        print(f"    Alt counts: r={stats['alt_r']:.4f}, p={stats['alt_p']:.2e}")
        print(f"    N variants: {stats['n_variants']}")

    # Save merged counts
    print(f"\nSaving to: {output_file}")
    merged_df.to_csv(output_file, sep='\t', index=False)

    # Print summary stats
    print("\nSummary statistics:")
    print(merged_df[['wasp2_ref', 'wasp2_alt', 'gatk_ref', 'gatk_alt']].describe())

    return merged_df, correlations


def main():
    parser = argparse.ArgumentParser(
        description='Generate count comparison data for Figure 2 Panel B'
    )
    parser.add_argument(
        '--dataset',
        choices=['hg00731', 'gm12878', 'both'],
        default='hg00731',
        help='Dataset to process'
    )
    parser.add_argument(
        '--data-dir',
        type=Path,
        default=Path(__file__).parent.parent / 'data',
        help='Base data directory'
    )

    args = parser.parse_args()

    print("="*60)
    print("Count Comparison Data Generation")
    print("="*60)

    # Process HG00731
    if args.dataset in ['hg00731', 'both']:
        hg00731_dir = args.data_dir / 'hg00731'

        if not hg00731_dir.exists():
            print(f"\nWarning: HG00731 directory not found: {hg00731_dir}")
        else:
            wasp2_file = hg00731_dir / 'wasp2_counts.tsv'
            gatk_file = hg00731_dir / 'gatk_counts.filtered.table'
            phaser_file = hg00731_dir / 'phaser.haplotypic_counts.txt'
            output_file = hg00731_dir / 'count_comparison.tsv'

            if all(f.exists() for f in [wasp2_file, gatk_file]):
                process_dataset(
                    'HG00731 RNA-seq',
                    wasp2_file,
                    gatk_file,
                    phaser_file,
                    output_file
                )
            else:
                print(f"\nWarning: Missing input files for HG00731")
                print(f"  WASP2: {wasp2_file.exists()}")
                print(f"  GATK: {gatk_file.exists()}")
                print(f"  phASER: {phaser_file.exists()}")

    # Process GM12878
    if args.dataset in ['gm12878', 'both']:
        gm12878_dir = args.data_dir / 'gm12878'

        if not gm12878_dir.exists():
            print(f"\nWarning: GM12878 directory not found: {gm12878_dir}")
        else:
            wasp2_file = gm12878_dir / 'wasp2_counts.tsv'
            gatk_file = gm12878_dir / 'gatk_counts.table'
            phaser_file = gm12878_dir / 'phaser.haplotypic_counts.txt'
            output_file = gm12878_dir / 'count_comparison.tsv'

            if all(f.exists() for f in [wasp2_file, gatk_file]):
                process_dataset(
                    'GM12878 ATAC-seq',
                    wasp2_file,
                    gatk_file,
                    phaser_file,
                    output_file
                )
            else:
                print(f"\nWarning: Missing input files for GM12878")
                print(f"  WASP2: {wasp2_file.exists()}")
                print(f"  GATK: {gatk_file.exists()}")
                print(f"  phASER: {phaser_file.exists()}")

    print("\n" + "="*60)
    print("Count comparison data generation complete!")
    print("="*60)


if __name__ == '__main__':
    main()
