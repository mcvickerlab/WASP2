#!/usr/bin/env python3
"""
Validate Rust counting implementation against Python baseline.
Tests equivalence and benchmarks performance.
"""
import timeit
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from counting.run_counting import run_count_variants
import polars as pl


def test_equivalence():
    """Test that Rust and Python produce identical results."""
    print("=" * 70)
    print("RUST COUNTING VALIDATION TEST")
    print("=" * 70)
    print()

    # Test files
    bam_file = "test_data/CD4_ATACseq_Day1_merged_filtered.sort.bam"
    vcf_file = "test_data/filter_chr10.vcf"
    region_file = "test_data/NA12878_snps_chr10.bed"
    sample = "NA12878"

    # Output files
    rust_output = "test_rust_counts.tsv"
    python_output = "test_python_counts.tsv"

    print(f"Test Data:")
    print(f"  BAM:     {bam_file}")
    print(f"  VCF:     {vcf_file}")
    print(f"  Regions: {region_file}")
    print(f"  Sample:  {sample}")
    print()

    # Test 1: Python baseline
    print("-" * 70)
    print("TEST 1: Python Baseline (--no-rust)")
    print("-" * 70)
    start = timeit.default_timer()

    run_count_variants(
        bam_file=bam_file,
        vcf_file=vcf_file,
        region_file=region_file,
        samples=sample,
        out_file=python_output,
        use_rust=False
    )

    python_time = timeit.default_timer() - start
    print(f"\n✓ Python completed in {python_time:.2f} seconds")
    print()

    # Test 2: Rust implementation
    print("-" * 70)
    print("TEST 2: Rust Implementation (--use-rust)")
    print("-" * 70)
    start = timeit.default_timer()

    run_count_variants(
        bam_file=bam_file,
        vcf_file=vcf_file,
        region_file=region_file,
        samples=sample,
        out_file=rust_output,
        use_rust=True
    )

    rust_time = timeit.default_timer() - start
    print(f"\n✓ Rust completed in {rust_time:.2f} seconds")
    print()

    # Load results
    print("-" * 70)
    print("TEST 3: Equivalence Check")
    print("-" * 70)

    python_df = pl.read_csv(python_output, separator="\t")
    rust_df = pl.read_csv(rust_output, separator="\t")

    print(f"Python output: {python_df.height} rows, {python_df.width} columns")
    print(f"Rust output:   {rust_df.height} rows, {rust_df.width} columns")
    print()

    # Check row counts match
    if python_df.height != rust_df.height:
        print(f"❌ FAILED: Row count mismatch!")
        print(f"   Python: {python_df.height} rows")
        print(f"   Rust:   {rust_df.height} rows")
        return False

    # Check columns match
    if set(python_df.columns) != set(rust_df.columns):
        print(f"❌ FAILED: Column mismatch!")
        print(f"   Python cols: {python_df.columns}")
        print(f"   Rust cols:   {rust_df.columns}")
        return False

    # Sort both for comparison
    sort_cols = ["chrom", "pos"]
    python_sorted = python_df.sort(sort_cols)
    rust_sorted = rust_df.sort(sort_cols)

    # Compare count columns
    count_cols = ["ref_count", "alt_count", "other_count"]

    mismatches = 0
    for col in count_cols:
        if col in python_sorted.columns:
            diff = (python_sorted[col] - rust_sorted[col]).abs().sum()
            if diff > 0:
                print(f"⚠️  {col}: {diff} total difference")
                mismatches += 1
                # Show first few mismatches
                mask = python_sorted[col] != rust_sorted[col]
                diff_rows = python_sorted.filter(mask).head(5)
                print(f"   First mismatches:")
                print(f"   {diff_rows[['chrom', 'pos', col]]}")
            else:
                print(f"✓ {col}: Perfect match!")

    print()

    # Performance summary
    print("-" * 70)
    print("PERFORMANCE SUMMARY")
    print("-" * 70)
    speedup = python_time / rust_time if rust_time > 0 else 0
    print(f"Python time: {python_time:.2f}s")
    print(f"Rust time:   {rust_time:.2f}s")
    print(f"Speedup:     {speedup:.2f}x")
    print()

    if speedup >= 4:
        print(f"🚀 EXCELLENT: {speedup:.1f}x speedup achieved!")
    elif speedup >= 2:
        print(f"✓ GOOD: {speedup:.1f}x speedup achieved")
    elif speedup >= 1.5:
        print(f"⚠️  MODEST: Only {speedup:.1f}x speedup")
    else:
        print(f"❌ SLOWER: Rust is slower than Python!")

    print()

    # Final verdict
    if mismatches == 0:
        print("=" * 70)
        print("✅ SUCCESS: Rust implementation is correct and faster!")
        print("=" * 70)
        return True
    else:
        print("=" * 70)
        print(f"❌ FAILED: {mismatches} count columns have mismatches")
        print("=" * 70)
        return False


if __name__ == "__main__":
    try:
        success = test_equivalence()
        sys.exit(0 if success else 1)
    except Exception as e:
        print(f"\n❌ ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
