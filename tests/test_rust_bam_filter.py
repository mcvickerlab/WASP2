"""Test Rust BAM filter against samtools ground truth.

Uses existing validation benchmark data from star_wasp_comparison to verify
that Rust filter_bam_by_variants produces identical read sets to samtools.
"""
import os
import sys
import tempfile
from pathlib import Path

import pysam

# Add src to path for wasp2_rust import
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

# Test data paths (existing validation benchmark)
BENCHMARK_DIR = Path(__file__).parent.parent / "benchmarking" / "star_wasp_comparison" / "results" / "wasp2_run"
INPUT_BAM = BENCHMARK_DIR / "A_sorted.bam"
VARIANT_BED = BENCHMARK_DIR / "HG00731_het_only_chr.bed"
GROUND_TRUTH_REMAP = BENCHMARK_DIR / "A_sorted_to_remap.bam"
GROUND_TRUTH_KEEP = BENCHMARK_DIR / "A_sorted_keep.bam"


def get_read_names_from_bam(bam_path: str) -> set:
    """Extract unique read names from a BAM file."""
    names = set()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            names.add(read.query_name)
    return names


def test_rust_filter_matches_samtools():
    """Verify Rust filter output matches samtools ground truth."""
    # Skip if test data doesn't exist
    if not INPUT_BAM.exists():
        print(f"SKIP: Test data not found at {INPUT_BAM}")
        return

    # Import Rust function
    try:
        from wasp2_rust import filter_bam_by_variants_py as rust_filter
    except ImportError as e:
        print(f"SKIP: wasp2_rust not available: {e}")
        return

    print("=== Rust BAM Filter vs Samtools Comparison ===")
    print(f"Input BAM: {INPUT_BAM}")
    print(f"Variant BED: {VARIANT_BED}")
    print(f"Ground truth remap: {GROUND_TRUTH_REMAP}")
    print(f"Ground truth keep: {GROUND_TRUTH_KEEP}")

    # Run Rust filter to temp files
    with tempfile.TemporaryDirectory() as tmpdir:
        rust_remap = os.path.join(tmpdir, "rust_remap.bam")
        rust_keep = os.path.join(tmpdir, "rust_keep.bam")

        print("\n--- Running Rust filter ---")
        import time
        start = time.time()

        remap_reads, keep_reads, unique_remap_names = rust_filter(
            str(INPUT_BAM),
            str(VARIANT_BED),
            rust_remap,
            rust_keep,
            is_paired=True,
            threads=8
        )

        elapsed = time.time() - start
        print(f"Rust filter completed in {elapsed:.2f}s")
        print(f"  Remap reads: {remap_reads}")
        print(f"  Keep reads: {keep_reads}")
        print(f"  Unique remap names: {unique_remap_names}")

        # Extract read names from outputs
        print("\n--- Extracting read names ---")

        print("Reading Rust remap BAM...")
        rust_remap_names = get_read_names_from_bam(rust_remap)
        print(f"  Rust remap: {len(rust_remap_names)} unique names")

        print("Reading ground truth remap BAM...")
        gt_remap_names = get_read_names_from_bam(str(GROUND_TRUTH_REMAP))
        print(f"  Ground truth remap: {len(gt_remap_names)} unique names")

        # Compare
        print("\n--- Comparison ---")

        # Check for exact match
        if rust_remap_names == gt_remap_names:
            print("✅ PASS: Rust remap names EXACTLY match samtools ground truth!")
        else:
            # Find differences
            only_rust = rust_remap_names - gt_remap_names
            only_gt = gt_remap_names - rust_remap_names

            print(f"❌ FAIL: Read name mismatch!")
            print(f"  In Rust but not ground truth: {len(only_rust)}")
            print(f"  In ground truth but not Rust: {len(only_gt)}")

            # Show some examples
            if only_rust:
                print(f"\n  Sample Rust-only names: {list(only_rust)[:5]}")
            if only_gt:
                print(f"\n  Sample ground-truth-only names: {list(only_gt)[:5]}")

            # Overlap percentage
            overlap = len(rust_remap_names & gt_remap_names)
            total = len(rust_remap_names | gt_remap_names)
            print(f"\n  Overlap: {overlap}/{total} = {100*overlap/total:.2f}%")

            return False

    return True


def main():
    """Run the test."""
    success = test_rust_filter_matches_samtools()
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
