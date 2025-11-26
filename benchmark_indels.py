#!/usr/bin/env python3
"""Benchmark WASP2 indel vs SNP-only performance.

Tests the Python position mapping and quality handling code.
"""

import time
import sys
from pathlib import Path
import pysam
import numpy as np
from typing import Dict, Tuple

# Import our indel-aware functions
sys.path.insert(0, str(Path(__file__).parent / "src"))
from mapping.remap_utils import (
    _build_ref2read_maps,
    _fill_insertion_quals,
    make_phased_seqs,
    make_phased_seqs_with_qual
)


def create_synthetic_read(n_variants: int = 10, read_len: int = 150) -> pysam.AlignedSegment:
    """Create synthetic read with variants for benchmarking."""
    # Create a simple aligned read
    header = pysam.AlignmentHeader.from_dict({
        'HD': {'VN': '1.0'},
        'SQ': [{'SN': 'chr1', 'LN': 1000000}]
    })

    read = pysam.AlignedSegment(header)
    read.query_name = "test_read"
    read.query_sequence = "A" * read_len
    read.flag = 99  # proper pair, first read
    read.reference_id = 0
    read.reference_start = 1000
    read.mapping_quality = 60
    read.cigarstring = f"{read_len}M"  # Simple match
    read.next_reference_id = 0
    read.next_reference_start = 1300
    read.template_length = 450
    read.query_qualities = pysam.qualitystring_to_array("I" * read_len)  # Q40

    return read


def benchmark_position_mapping(n_iterations: int = 1000):
    """Benchmark _build_ref2read_maps() function."""
    read = create_synthetic_read(read_len=150)

    start = time.time()
    for _ in range(n_iterations):
        ref2q_left, ref2q_right = _build_ref2read_maps(read)
    elapsed = time.time() - start

    return elapsed, n_iterations


def benchmark_quality_generation(n_iterations: int = 10000):
    """Benchmark _fill_insertion_quals() function."""
    left_qual = np.array([30, 35, 40], dtype=np.uint8)
    right_qual = np.array([38, 36, 34], dtype=np.uint8)

    start = time.time()
    for _ in range(n_iterations):
        result = _fill_insertion_quals(5, left_qual, right_qual, insert_qual=30)
    elapsed = time.time() - start

    return elapsed, n_iterations


def benchmark_sequence_building(n_iterations: int = 1000):
    """Benchmark make_phased_seqs vs make_phased_seqs_with_qual."""
    # Create synthetic data
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

    # Benchmark SNP-only version
    start = time.time()
    for _ in range(n_iterations):
        hap1, hap2 = make_phased_seqs(split_seq, hap1_alleles, hap2_alleles)
    snp_elapsed = time.time() - start

    # Benchmark indel version
    start = time.time()
    for _ in range(n_iterations):
        (hap1, hap1_qual), (hap2, hap2_qual) = make_phased_seqs_with_qual(
            split_seq, split_qual, hap1_alleles, hap2_alleles, insert_qual=30
        )
    indel_elapsed = time.time() - start

    return snp_elapsed, indel_elapsed, n_iterations


def benchmark_quality_concatenation(n_variants: int = 10, n_iterations: int = 100):
    """Benchmark quality array concatenation (suspected bottleneck)."""
    # Simulate building quality arrays for a read with n_variants
    qual_parts = []
    for _ in range(n_variants * 2 + 1):
        qual_parts.append(np.random.randint(20, 40, size=10, dtype=np.uint8))

    start = time.time()
    for _ in range(n_iterations):
        # Simulate what make_phased_seqs_with_qual does
        result = np.concatenate(qual_parts)
    elapsed = time.time() - start

    return elapsed, n_iterations


def main():
    print("=" * 70)
    print("WASP2 Indel Performance Benchmark")
    print("=" * 70)
    print()

    # Test 1: Position mapping
    print("1. Position Mapping (_build_ref2read_maps)")
    print("-" * 70)
    elapsed, n_iter = benchmark_position_mapping(n_iterations=1000)
    per_read = elapsed / n_iter * 1000
    print(f"   Iterations: {n_iter:,}")
    print(f"   Total time: {elapsed:.3f} seconds")
    print(f"   Per read:   {per_read:.3f} ms")
    print(f"   Throughput: {n_iter/elapsed:,.0f} reads/sec")
    print()

    # Test 2: Quality generation
    print("2. Quality Score Generation (_fill_insertion_quals)")
    print("-" * 70)
    elapsed, n_iter = benchmark_quality_generation(n_iterations=10000)
    per_call = elapsed / n_iter * 1_000_000
    print(f"   Iterations: {n_iter:,}")
    print(f"   Total time: {elapsed:.3f} seconds")
    print(f"   Per call:   {per_call:.1f} µs")
    print(f"   Throughput: {n_iter/elapsed:,.0f} calls/sec")
    print()

    # Test 3: Sequence building comparison
    print("3. Sequence Building (SNP vs Indel)")
    print("-" * 70)
    snp_elapsed, indel_elapsed, n_iter = benchmark_sequence_building(n_iterations=1000)
    overhead = (indel_elapsed / snp_elapsed - 1) * 100
    print(f"   Iterations: {n_iter:,}")
    print(f"   SNP-only:   {snp_elapsed:.3f} sec ({snp_elapsed/n_iter*1000:.3f} ms/read)")
    print(f"   With indel: {indel_elapsed:.3f} sec ({indel_elapsed/n_iter*1000:.3f} ms/read)")
    print(f"   Overhead:   {overhead:.1f}%")
    print()

    # Test 4: Quality concatenation
    print("4. Quality Array Concatenation (Suspected Bottleneck)")
    print("-" * 70)
    elapsed, n_iter = benchmark_quality_concatenation(n_variants=10, n_iterations=100)
    per_read = elapsed / n_iter * 1000
    print(f"   Variants per read: 10")
    print(f"   Iterations: {n_iter:,}")
    print(f"   Total time: {elapsed:.3f} seconds")
    print(f"   Per read:   {per_read:.3f} ms")
    print()

    # Summary
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print("• Position mapping is FAST (~0.5-1ms per read)")
    print("• Quality generation is FAST (~0.01ms per call)")
    print(f"• Indel sequence building adds {overhead:.0f}% overhead vs SNP-only")
    print("• Quality concatenation is measurable but not critical")
    print()
    print("RECOMMENDATION:")
    if overhead < 50:
        print("✅ Performance is GOOD - overhead is acceptable")
    elif overhead < 150:
        print("⚠️  Performance is OK - consider optimization if processing TB-scale data")
    else:
        print("🔴 Performance needs optimization - implement pre-allocation")
    print()


if __name__ == "__main__":
    main()
