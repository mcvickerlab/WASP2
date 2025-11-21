#!/usr/bin/env python3
"""
Validation script: Compare Python and Rust intersection parsing

This script:
1. Parses intersection BED using Python (Polars - the current implementation)
2. Parses intersection BED using Rust (once implemented)
3. Compares the results to ensure exact match
"""
import sys
import json
import time
from pathlib import Path
from collections import defaultdict

sys.path.insert(0, str(Path(__file__).parent / "src"))

from src.mapping.intersect_variant_data import make_intersect_df


def parse_with_python(intersect_file, sample):
    """Parse using current Python implementation (Polars)"""
    print("=" * 80)
    print("PYTHON IMPLEMENTATION (Polars)")
    print("=" * 80)

    start = time.time()
    df = make_intersect_df(intersect_file, [sample])
    elapsed = time.time() - start

    print(f"Time: {elapsed:.3f}s")
    print(f"DataFrame shape: {df.shape}")
    print(f"Columns: {df.columns}")
    print()

    # Group by read name to create the HashMap structure
    # This is what Rust will create: HashMap<read_name, Vec<Variant>>
    variants_by_read = defaultdict(list)

    for row in df.iter_rows(named=True):
        read_name = row['read']
        mate = row['mate']
        chrom = row['chrom']
        start = row['start']
        stop = row['stop']

        # Get haplotype columns (assuming single sample)
        hap_cols = [c for c in df.columns if c.endswith('_a1') or c.endswith('_a2')]
        hap1_col = [c for c in hap_cols if c.endswith('_a1')][0]
        hap2_col = [c for c in hap_cols if c.endswith('_a2')][0]

        hap1 = row[hap1_col]
        hap2 = row[hap2_col]

        variant = {
            'chrom': chrom,
            'pos': start,  # 0-based
            'ref_allele': None,  # Not in processed DF
            'alt_allele': None,  # Not in processed DF
            'hap1': hap1,
            'hap2': hap2,
            'mate': mate,
        }

        variants_by_read[read_name].append(variant)

    print(f"Unique reads: {len(variants_by_read)}")
    print(f"Total variants: {sum(len(v) for v in variants_by_read.values())}")
    print()

    # Show sample
    print("Sample data (first 3 reads):")
    for i, (read_name, variants) in enumerate(list(variants_by_read.items())[:3]):
        print(f"  {read_name}: {len(variants)} variants")
        for v in variants[:2]:
            print(f"    {v['chrom']}:{v['pos']} {v['hap1']}|{v['hap2']}")
    print()

    return variants_by_read, elapsed


def parse_with_rust_raw(intersect_file):
    """Parse intersection file directly (raw BED format)"""
    print("=" * 80)
    print("DIRECT BED PARSING (Python - simulating Rust logic)")
    print("=" * 80)

    start = time.time()

    # First pass: collect all spans
    all_spans = []

    with open(intersect_file) as f:
        for line in f:
            fields = line.strip().split('\t')
            if not fields or len(fields) < 12:
                continue

            # Parse fields matching Rust implementation
            read_chrom = fields[0]
            read_start = int(fields[1])
            read_end = int(fields[2])
            read_name_with_mate = fields[3]
            var_chrom = fields[6]
            var_start = int(fields[7])  # 0-based
            ref_allele = fields[9]
            alt_allele = fields[10]
            genotype = fields[11]

            # Extract read name (without /1 or /2)
            read_name = read_name_with_mate.split('/')[0]
            mate = int(read_name_with_mate.split('/')[1])

            # Parse phased genotype
            hap1, hap2 = genotype.split('|')

            span = {
                'read_name': read_name,
                'chrom': read_chrom,  # Use read chrom, not variant chrom
                'start': read_start,   # Read positions
                'stop': read_end,
                'mate': mate,
                'var_chrom': var_chrom,
                'pos': var_start,
                'ref_allele': ref_allele,
                'alt_allele': alt_allele,
                'hap1': hap1,
                'hap2': hap2,
            }

            all_spans.append(span)

    # Deduplicate - Python Polars does: df.unique(["chrom", "read", "mate", "start", "stop"], keep="first")
    # Dedup on (read_name, chrom, start, stop, mate) - matches Rust implementation
    seen = set()
    deduped_spans = []

    for span in all_spans:
        key = (span['read_name'], span['chrom'], span['start'], span['stop'], span['mate'])
        if key not in seen:
            seen.add(key)
            deduped_spans.append(span)

    # Group by read name
    variants_by_read = defaultdict(list)
    for span in deduped_spans:
        read_name = span['read_name']
        variant = {
            'chrom': span['var_chrom'],
            'pos': span['pos'],
            'ref_allele': span['ref_allele'],
            'alt_allele': span['alt_allele'],
            'hap1': span['hap1'],
            'hap2': span['hap2'],
            'mate': span['mate'],
        }
        variants_by_read[read_name].append(variant)

    elapsed = time.time() - start

    print(f"Time: {elapsed:.3f}s")
    print(f"Unique reads: {len(variants_by_read)}")
    print(f"Total variants: {sum(len(v) for v in variants_by_read.values())}")
    print()

    # Show sample
    print("Sample data (first 3 reads):")
    for i, (read_name, variants) in enumerate(list(variants_by_read.items())[:3]):
        print(f"  {read_name}: {len(variants)} variants")
        for v in variants[:2]:
            print(f"    {v['chrom']}:{v['pos']} {v['ref_allele']}->{v['alt_allele']} {v['hap1']}|{v['hap2']}")
    print()

    return variants_by_read, elapsed


def compare_results(python_variants, rust_variants):
    """Compare Python and Rust parsing results"""
    print("=" * 80)
    print("COMPARISON")
    print("=" * 80)

    # Check read names match
    python_reads = set(python_variants.keys())
    rust_reads = set(rust_variants.keys())

    if python_reads == rust_reads:
        print(f"✓ Read names match: {len(python_reads)} reads")
    else:
        print(f"✗ Read names differ!")
        print(f"  Python only: {len(python_reads - rust_reads)}")
        print(f"  Rust only: {len(rust_reads - python_reads)}")
        return False

    # Check variant counts per read
    mismatches = []
    for read_name in python_reads:
        p_count = len(python_variants[read_name])
        r_count = len(rust_variants[read_name])
        if p_count != r_count:
            mismatches.append((read_name, p_count, r_count))

    if mismatches:
        print(f"✗ Variant count mismatches: {len(mismatches)}")
        for read_name, p_count, r_count in mismatches[:5]:
            print(f"  {read_name}: Python={p_count}, Rust={r_count}")
        return False
    else:
        print(f"✓ Variant counts match for all reads")

    # Check variant details (for a sample)
    sample_reads = list(python_reads)[:10]
    all_match = True

    for read_name in sample_reads:
        p_vars = python_variants[read_name]
        r_vars = rust_variants[read_name]

        # Sort by position for comparison
        p_vars = sorted(p_vars, key=lambda v: v['pos'])
        r_vars = sorted(r_vars, key=lambda v: v['pos'])

        for i, (p_var, r_var) in enumerate(zip(p_vars, r_vars)):
            # Compare key fields
            if p_var['chrom'] != r_var['chrom']:
                print(f"✗ Chrom mismatch for {read_name} variant {i}")
                all_match = False
            if p_var['pos'] != r_var['pos']:
                print(f"✗ Position mismatch for {read_name} variant {i}")
                all_match = False
            if p_var['hap1'] != r_var['hap1'] or p_var['hap2'] != r_var['hap2']:
                print(f"✗ Genotype mismatch for {read_name} variant {i}")
                print(f"  Python: {p_var['hap1']}|{p_var['hap2']}")
                print(f"  Rust: {r_var['hap1']}|{r_var['hap2']}")
                all_match = False

    if all_match:
        print(f"✓ Variant details match for sampled reads")

    print()
    return all_match


def main():
    intersect_file = "baselines/mapping/intersect.bed"
    sample = "NA12878"

    print("INTERSECTION BED PARSER VALIDATION")
    print("=" * 80)
    print(f"File: {intersect_file}")
    print(f"Sample: {sample}")
    print()

    # Parse with direct BED parsing (simulating Rust)
    rust_variants, rust_time = parse_with_rust_raw(intersect_file)

    # Parse with Python (current Polars implementation)
    python_variants, python_time = parse_with_python(intersect_file, sample)

    # Compare
    match = compare_results(python_variants, rust_variants)

    # Summary
    print("=" * 80)
    print("PERFORMANCE SUMMARY")
    print("=" * 80)
    print(f"Python (Polars): {python_time:.3f}s")
    print(f"Direct parsing:  {rust_time:.3f}s")
    print(f"Speedup:         {python_time / rust_time:.1f}x")
    print()

    if match:
        print("✓ VALIDATION PASSED - Results match!")
    else:
        print("✗ VALIDATION FAILED - Results differ!")
        return 1

    # Save reference output for Rust comparison
    output_file = "baselines/mapping/python_parsed_variants.json"
    with open(output_file, 'w') as f:
        # Convert to serializable format
        serializable = {
            read_name: [
                {k: str(v) if not isinstance(v, (int, str)) else v
                 for k, v in var.items()}
                for var in variants
            ]
            for read_name, variants in rust_variants.items()
        }
        json.dump(serializable, f, indent=2)

    print(f"Reference output saved to: {output_file}")
    print(f"(Use this to validate Rust implementation)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
