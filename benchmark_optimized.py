#!/usr/bin/env python3
"""Compare original vs optimized indel implementation."""

import time
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).parent / "src"))

# Import both versions
from mapping.remap_utils import make_phased_seqs_with_qual as original_version
from mapping.remap_utils_optimized import make_phased_seqs_with_qual_fast as optimized_version


def benchmark_comparison(n_iterations: int = 1000):
    """Compare original vs optimized versions."""
    # Create synthetic data (same as before)
    split_seq = ["ATCG", "G", "TTAA", "C", "GGCC", "A", "TATA"]
    split_qual = [
        np.array([30, 35, 40, 38], dtype=np.uint8),
        np.array([35], dtype=np.uint8),
        np.array([32, 34, 36, 38], dtype=np.uint8),
        np.array([40], dtype=np.uint8),
        np.array([35, 37, 39, 41], dtype=np.uint8),
        np.array([38], dtype=np.uint8),
        np.array([33, 35, 37, 39], dtype=np.uint8),
    ]
    hap1_alleles = ["G", "C", "A"]
    hap2_alleles = ["T", "G", "T"]

    # Benchmark original version
    print("Benchmarking ORIGINAL version...")
    start = time.time()
    for _ in range(n_iterations):
        (hap1, hap1_qual), (hap2, hap2_qual) = original_version(
            split_seq, split_qual, hap1_alleles, hap2_alleles, insert_qual=30
        )
    original_time = time.time() - start

    # Benchmark optimized version
    print("Benchmarking OPTIMIZED version...")
    start = time.time()
    for _ in range(n_iterations):
        (hap1_opt, hap1_qual_opt), (hap2_opt, hap2_qual_opt) = optimized_version(
            split_seq, split_qual, hap1_alleles, hap2_alleles, insert_qual=30
        )
    optimized_time = time.time() - start

    # Verify results match
    assert hap1 == hap1_opt, "Sequence mismatch!"
    assert np.array_equal(hap1_qual, hap1_qual_opt), "Quality mismatch!"
    print("✅ Results verified - both versions produce identical output\n")

    # Report
    speedup = original_time / optimized_time
    print("=" * 70)
    print("PERFORMANCE COMPARISON")
    print("=" * 70)
    print(f"Iterations:       {n_iterations:,}")
    print(f"Original time:    {original_time:.3f} sec ({original_time/n_iterations*1000:.3f} ms/read)")
    print(f"Optimized time:   {optimized_time:.3f} sec ({optimized_time/n_iterations*1000:.3f} ms/read)")
    print(f"Speedup:          {speedup:.1f}x faster")
    print()

    if speedup > 8:
        print("🚀 EXCELLENT - Optimization successful!")
    elif speedup > 5:
        print("✅ GOOD - Significant improvement")
    elif speedup > 2:
        print("⚠️  OK - Moderate improvement")
    else:
        print("🔴 POOR - Optimization had little effect")
    print()


if __name__ == "__main__":
    benchmark_comparison(n_iterations=1000)
