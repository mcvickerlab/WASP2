#!/usr/bin/env python3
"""
Count reference vs alternate alleles at heterozygous variants.

Used for validation: measures allelic imbalance after WASP filtering.
"""

import pysam
import polars as pl
import argparse
from pathlib import Path
from scipy.stats import binom_test
from typing import Optional


def count_alleles_at_variant(bam, variant, min_baseq=20, min_mapq=10):
    """
    Count REF vs ALT alleles at a single variant position.

    Args:
        bam: pysam.AlignmentFile
        variant: pysam.VariantRecord
        min_baseq: Minimum base quality
        min_mapq: Minimum mapping quality

    Returns:
        (ref_count, alt_count, total_depth)
    """
    ref_count = 0
    alt_count = 0
    other_count = 0

    chrom = variant.contig
    pos = variant.pos - 1  # pysam uses 0-based
    ref_allele = variant.ref
    alt_allele = variant.alts[0] if variant.alts else None

    if alt_allele is None:
        return 0, 0, 0

    # Determine variant type
    is_snp = len(ref_allele) == 1 and len(alt_allele) == 1
    is_ins = len(alt_allele) > len(ref_allele)
    is_del = len(ref_allele) > len(alt_allele)

    for pileupcolumn in bam.pileup(chrom, pos, pos + 1,
                                    truncate=True,
                                    min_base_quality=min_baseq,
                                    min_mapping_quality=min_mapq,
                                    ignore_overlaps=False):

        if pileupcolumn.pos != pos:
            continue

        for pileupread in pileupcolumn.pileups:
            if pileupread.is_del or pileupread.is_refskip:
                continue

            read = pileupread.alignment

            # Get query position
            query_pos = pileupread.query_position

            if query_pos is None:
                continue

            if is_snp:
                # Simple SNP check
                base = read.query_sequence[query_pos]
                if base == ref_allele:
                    ref_count += 1
                elif base == alt_allele:
                    alt_count += 1
                else:
                    other_count += 1

            elif is_ins:
                # Insertion: check if read has extra bases
                # Get sequence starting at variant position
                seq_window = read.query_sequence[query_pos:query_pos + len(alt_allele)]

                if seq_window == alt_allele:
                    alt_count += 1
                elif read.query_sequence[query_pos] == ref_allele:
                    ref_count += 1
                else:
                    other_count += 1

            elif is_del:
                # Deletion: check if read is missing bases
                # This is tricky - need to parse CIGAR
                # Simplified: check sequence length
                seq_window = read.query_sequence[query_pos:query_pos + len(ref_allele)]

                if len(seq_window) == len(ref_allele) and seq_window == ref_allele:
                    ref_count += 1
                elif len(seq_window) < len(ref_allele):
                    # Read has deletion
                    alt_count += 1
                else:
                    other_count += 1

    total_depth = ref_count + alt_count + other_count

    return ref_count, alt_count, total_depth


def count_all_variants(bam_file: str, vcf_file: str,
                       variant_type: Optional[str] = None,
                       min_depth: int = 10) -> pl.DataFrame:
    """
    Count alleles at all heterozygous variants in VCF.

    Args:
        bam_file: Path to BAM file (WASP-filtered)
        vcf_file: Path to VCF with variants
        variant_type: Filter by type ("SNP", "INS", "DEL", or None for all)
        min_depth: Minimum depth to include variant

    Returns:
        Polars DataFrame with allelic counts
    """

    bam = pysam.AlignmentFile(bam_file)
    vcf = pysam.VariantFile(vcf_file)

    results = []

    for variant in vcf:
        # Check if heterozygous
        if not any(sample.get('GT', (None,))[0] != sample.get('GT', (None,))[1]
                   for sample in variant.samples.values()):
            continue  # Skip homozygous

        # Determine variant type
        ref = variant.ref
        alt = variant.alts[0] if variant.alts else None

        if alt is None:
            continue

        if len(ref) == 1 and len(alt) == 1:
            vtype = "SNP"
        elif len(alt) > len(ref):
            vtype = "INS"
        elif len(ref) > len(alt):
            vtype = "DEL"
        else:
            vtype = "COMPLEX"

        # Filter by type if requested
        if variant_type and vtype != variant_type:
            continue

        # Count alleles
        ref_count, alt_count, total_depth = count_alleles_at_variant(bam, variant)

        if total_depth < min_depth:
            continue  # Skip low coverage

        # Calculate metrics
        if alt_count > 0:
            ratio = ref_count / alt_count
        else:
            ratio = float('inf')

        # Binomial test for allelic imbalance
        p_value = binom_test(ref_count, ref_count + alt_count, p=0.5, alternative='two-sided')

        results.append({
            'chrom': variant.contig,
            'pos': variant.pos,
            'id': variant.id if variant.id else f"{variant.contig}:{variant.pos}",
            'ref': ref,
            'alt': alt,
            'type': vtype,
            'ref_count': ref_count,
            'alt_count': alt_count,
            'total_depth': total_depth,
            'ref_freq': ref_count / total_depth if total_depth > 0 else 0,
            'alt_freq': alt_count / total_depth if total_depth > 0 else 0,
            'ratio': ratio,
            'log2_ratio': float('inf') if ratio == float('inf') else (
                -float('inf') if ratio == 0 else pl.Series([ratio]).log(2).item()
            ),
            'p_value': p_value,
            'significant': p_value < 0.05,
        })

    return pl.DataFrame(results)


def main():
    parser = argparse.ArgumentParser(description='Count alleles at heterozygous variants')
    parser.add_argument('bam', help='BAM file (WASP-filtered)')
    parser.add_argument('vcf', help='VCF file with variants')
    parser.add_argument('-o', '--output', required=True, help='Output CSV file')
    parser.add_argument('--variant-type', choices=['SNP', 'INS', 'DEL'],
                        help='Filter by variant type')
    parser.add_argument('--min-depth', type=int, default=10,
                        help='Minimum depth (default: 10)')

    args = parser.parse_args()

    print(f"Counting alleles in {args.bam}")
    print(f"Using variants from {args.vcf}")

    if args.variant_type:
        print(f"Filtering for {args.variant_type} only")

    df = count_all_variants(args.bam, args.vcf,
                             variant_type=args.variant_type,
                             min_depth=args.min_depth)

    print(f"\nResults:")
    print(f"  Total variants: {len(df)}")
    print(f"  Significant ASE (p<0.05): {df.filter(pl.col('significant')).height}")
    print(f"  Mean depth: {df['total_depth'].mean():.1f}")
    print(f"  Mean REF frequency: {df['ref_freq'].mean():.3f}")

    if len(df) > 0:
        print(f"\nVariant type breakdown:")
        print(df.group_by('type').agg([
            pl.count().alias('count'),
            pl.col('ref_freq').mean().alias('mean_ref_freq')
        ]))

    df.write_csv(args.output)
    print(f"\n✅ Wrote results to {args.output}")


if __name__ == '__main__':
    main()
