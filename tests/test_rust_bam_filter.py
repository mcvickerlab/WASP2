"""Test Rust BAM filter against samtools ground truth.

Uses existing validation benchmark data from star_wasp_comparison to verify
that Rust filter_bam_by_variants produces identical read sets to samtools.
"""

import os
import sys
import tempfile
from pathlib import Path

import pysam
import pytest

# Add src to path for wasp2_rust import
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

# Test data paths (existing validation benchmark)
BENCHMARK_DIR = (
    Path(__file__).parent.parent / "benchmarking" / "star_wasp_comparison" / "results" / "wasp2_run"
)
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
    if not INPUT_BAM.exists():
        pytest.skip(f"Test data not found at {INPUT_BAM}")

    try:
        from wasp2_rust import filter_bam_by_variants_py as rust_filter
    except ImportError as e:
        pytest.skip(f"wasp2_rust not available: {e}")

    with tempfile.TemporaryDirectory() as tmpdir:
        rust_remap = os.path.join(tmpdir, "rust_remap.bam")
        rust_keep = os.path.join(tmpdir, "rust_keep.bam")

        remap_reads, keep_reads, unique_remap_names = rust_filter(
            str(INPUT_BAM), str(VARIANT_BED), rust_remap, rust_keep, is_paired=True, threads=8
        )

        rust_remap_names = get_read_names_from_bam(rust_remap)
        gt_remap_names = get_read_names_from_bam(str(GROUND_TRUTH_REMAP))

        only_rust = rust_remap_names - gt_remap_names
        only_gt = gt_remap_names - rust_remap_names

        assert rust_remap_names == gt_remap_names, (
            f"Read name mismatch!\n"
            f"  In Rust but not ground truth: {len(only_rust)}\n"
            f"  In ground truth but not Rust: {len(only_gt)}\n"
            f"  Sample Rust-only: {list(only_rust)[:5]}\n"
            f"  Sample GT-only: {list(only_gt)[:5]}"
        )
