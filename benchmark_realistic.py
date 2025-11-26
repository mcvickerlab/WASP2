#!/usr/bin/env python3
"""
Realistic benchmark using actual read lengths and variant counts.

Tests performance with:
- 150bp reads (realistic Illumina length)
- 5-20 heterozygous variants per read (realistic for high-diversity regions)
"""

import time
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).parent / "src"))

from mapping.remap_utils import make_phased_seqs_with_qual as original_version
from mapping.remap_utils_optimized import make_phased_seqs_with_qual_fast as optimized_version


def create_realistic_data(n_variants: int = 10, read_len: int = 150):
    """Create realistic read data for benchmarking."""
    # Simulate a 150bp read with n_variants heterozygous sites
    # Each variant splits the read into segments
    split_seq = []
    split_qual = []

    # Distribute variants evenly across read
    segment_size = read_len // (n_variants + 1)

    for i in range(n_variants + 1):
        if i < n_variants:
            # Non-variant segment
            seq = "A" * segment_size
            qual = np.random.randint(20, 40, size=segment_size, dtype=np.uint8)
            split_seq.append(seq)
            split_qual.append(qual)

            # Variant segment (1bp originally)
            split_seq.append("G")
            split_qual.append(np.array([35], dtype=np.uint8))
        else:
            # Last segment
            remaining = read_len - len("".join(split_seq))
            if remaining > 0:
                split_seq.append("A" * remaining)
                split_qual.append(np.random.randint(20, 40, size=remaining, dtype=np.uint8))

    # Create realistic alleles (mix of SNPs and small indels)
    hap1_alleles = []
    hap2_alleles = []
    for i in range(n_variants):
        allele_type = np.random.choice(["snp", "ins", "del"], p=[0.7, 0.15, 0.15])
        if allele_type == "snp":
            hap1_alleles.append("A")
            hap2_alleles.append("T")
        elif allele_type == "ins":
            # 1-3bp insertion
            ins_len = np.random.randint(1, 4)
            hap1_alleles.append("G" * ins_len)
            hap2_alleles.append("G")
        else:  # del
            # Deletion (reference has more bases)
            hap1_alleles.append("G")
            hap2_alleles.append("GGG")

    return split_seq, split_qual, hap1_alleles, hap2_alleles


def benchmark_with_variants(n_variants: int, n_iterations: int = 100):
    """Benchmark both versions with varying variant counts."""
    print(f"\nBenchmarking with {n_variants} variants per 150bp read:")

    split_seq, split_qual, hap1_alleles, hap2_alleles = create_realistic_data(n_variants)

    # Warmup
    for _ in range(10):
        original_version(split_seq, split_qual, hap1_alleles, hap2_alleles, insert_qual=30)

    # Benchmark original
    start = time.time()
    for _ in range(n_iterations):
        (h1, q1), (h2, q2) = original_version(
            split_seq, split_qual, hap1_alleles, hap2_alleles, insert_qual=30
        )
    original_time = time.time() - start

    # Warmup optimized
    for _ in range(10):
        optimized_version(split_seq, split_qual, hap1_alleles, hap2_alleles, insert_qual=30)

    # Benchmark optimized
    start = time.time()
    for _ in range(n_iterations):
        (h1_opt, q1_opt), (h2_opt, q2_opt) = optimized_version(
            split_seq, split_qual, hap1_alleles, hap2_alleles, insert_qual=30
        )
    optimized_time = time.time() - start

    # Verify correctness
    assert h1 == h1_opt and np.array_equal(q1, q1_opt), "Results don't match!"

    speedup = original_time / optimized_time
    print(f"  Original:  {original_time:.4f}s ({original_time/n_iterations*1000:.3f} ms/read)")
    print(f"  Optimized: {optimized_time:.4f}s ({optimized_time/n_iterations*1000:.3f} ms/read)")
    print(f"  Speedup:   {speedup:.2f}x")

    return original_time, optimized_time, speedup


def main():
    print("=" * 70)
    print("REALISTIC PERFORMANCE BENCHMARK")
    print("=" * 70)
    print("\nTesting with 150bp reads (Illumina-like)")
    print("Varying number of heterozygous variants per read")
    print()

    variant_counts = [5, 10, 15, 20, 30]
    results = []

    for n_variants in variant_counts:
        orig_time, opt_time, speedup = benchmark_with_variants(n_variants, n_iterations=100)
        results.append((n_variants, speedup))

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"{'Variants/Read':<15} {'Speedup':<15}")
    print("-" * 30)
    for n_variants, speedup in results:
        status = "✅" if speedup > 1.5 else "⚠️ " if speedup > 1.0 else "🔴"
        print(f"{n_variants:<15} {speedup:>6.2f}x {status}")

    print()
    avg_speedup = np.mean([s for _, s in results])
    if avg_speedup > 1.5:
        print(f"✅ Optimization SUCCESSFUL: Average {avg_speedup:.2f}x faster")
    elif avg_speedup > 1.0:
        print(f"⚠️  Optimization MARGINAL: Average {avg_speedup:.2f}x faster")
    else:
        print(f"🔴 Optimization INEFFECTIVE: Current code is already efficient")
        print("   Bottleneck is likely elsewhere (I/O, BAM parsing, etc.)")


if __name__ == "__main__":
    main()
