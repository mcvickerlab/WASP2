#!/usr/bin/env python3
"""
Compare variant BED files between WASP1/WASP2-Python and WASP2-Rust.
Identifies which variants are missing in WASP2-Rust and analyzes patterns.
"""

import pandas as pd
from pathlib import Path
from collections import Counter

# Paths to BED files
RESULTS_DIR = Path(__file__).parent.parent / "results"

WASP1_BED = RESULTS_DIR / "wasp1_snp_FIXED_2025-12-07_12-27-45/variants.bed"
PYTHON_BED = RESULTS_DIR / "wasp2python_snp_DEV_MT_2025-12-07_15-52-30/variants.bed"
RUST_BED = RESULTS_DIR / "wasp2rust_snp_fixed_2025-12-06_22-38-30/variants.bed"


def load_bed(path: Path) -> set:
    """Load BED file and return set of (chrom, start, end) tuples."""
    variants = set()
    with open(path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                variants.add((chrom, start, end))
    return variants


def load_bed_with_alleles(path: Path) -> dict:
    """Load BED file and return dict of (chrom, start, end) -> allele info."""
    variants = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                key = (chrom, start, end)
                # Store full line for analysis
                variants[key] = parts[3:] if len(parts) > 3 else []
    return variants


def main():
    print("=" * 60)
    print("Variant BED File Comparison: WASP1/Python vs WASP2-Rust")
    print("=" * 60)

    # Load variant sets
    print("\nLoading BED files...")
    wasp1_variants = load_bed(WASP1_BED)
    python_variants = load_bed(PYTHON_BED)
    rust_variants = load_bed(RUST_BED)

    print(f"  WASP1:        {len(wasp1_variants):,} variants")
    print(f"  WASP2-Python: {len(python_variants):,} variants")
    print(f"  WASP2-Rust:   {len(rust_variants):,} variants")

    # Check WASP1 vs Python (should be identical)
    wasp1_python_diff = wasp1_variants.symmetric_difference(python_variants)
    print(f"\n  WASP1 vs Python difference: {len(wasp1_python_diff)} variants")

    # Find missing in Rust
    missing_in_rust = wasp1_variants - rust_variants
    extra_in_rust = rust_variants - wasp1_variants

    print(f"\n  Missing in WASP2-Rust: {len(missing_in_rust):,} variants")
    print(f"  Extra in WASP2-Rust:   {len(extra_in_rust):,} variants")

    # Analyze missing variants by chromosome
    print("\n" + "=" * 60)
    print("Missing Variants by Chromosome")
    print("=" * 60)

    chrom_counts = Counter(v[0] for v in missing_in_rust)
    for chrom, count in sorted(chrom_counts.items(), key=lambda x: -x[1])[:15]:
        print(f"  {chrom}: {count:,}")

    # Sample missing variants
    print("\n" + "=" * 60)
    print("Sample of Missing Variants (first 20)")
    print("=" * 60)

    # Load with allele info for detailed analysis
    wasp1_full = load_bed_with_alleles(WASP1_BED)

    for i, var in enumerate(sorted(missing_in_rust)[:20]):
        chrom, start, end = var
        allele_info = wasp1_full.get(var, [])
        allele_str = '\t'.join(allele_info[:4])
        print(f"  {chrom}:{start}-{end}\t{allele_str}")

    # Analyze patterns in missing variants
    print("\n" + "=" * 60)
    print("Pattern Analysis")
    print("=" * 60)

    # Check variant sizes (SNV vs longer)
    snv_count = 0
    non_snv_count = 0
    for chrom, start, end in missing_in_rust:
        if end - start == 1:
            snv_count += 1
        else:
            non_snv_count += 1

    print(f"  SNVs (size=1bp):     {snv_count:,}")
    print(f"  Non-SNVs (size>1bp): {non_snv_count:,}")

    # Check for position patterns (e.g., near chromosome boundaries)
    edge_count = 0
    for chrom, start, end in missing_in_rust:
        if start < 10000:
            edge_count += 1

    print(f"  Near chromosome start (<10kb): {edge_count:,}")

    # Save missing variants to file
    output_file = RESULTS_DIR / "missing_variants_in_rust.bed"
    with open(output_file, 'w') as f:
        for var in sorted(missing_in_rust):
            chrom, start, end = var
            allele_info = wasp1_full.get(var, [])
            allele_str = '\t'.join(allele_info) if allele_info else ''
            f.write(f"{chrom}\t{start}\t{end}\t{allele_str}\n")

    print(f"\n  Saved missing variants to: {output_file}")

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"""
    Total missing variants: {len(missing_in_rust):,}
    Impact: ~0.1% of variants, but may affect more reads due to:
      - High-coverage regions
      - Multi-allelic sites
      - Variants in repetitive regions

    Next steps:
    1. Examine vcf_to_bed.rs filtering logic
    2. Check for multi-allelic site handling
    3. Verify genotype parsing (0/1 vs 1/0 vs 0|1)
    """)


if __name__ == "__main__":
    main()
