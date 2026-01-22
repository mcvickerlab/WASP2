#!/usr/bin/env python3
"""
Correctness tests for WASP2 indel implementation.

These tests verify that the indel-aware code produces correct results
by comparing against known ground truth examples.
"""

import sys
from pathlib import Path
import numpy as np
import pysam

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from mapping.remap_utils import (
    _build_ref2read_maps,
    _fill_insertion_quals,
    make_phased_seqs,
    make_phased_seqs_with_qual,
    make_multi_seqs_with_qual
)


def test_position_mapping_simple_match():
    """Test position mapping for a simple perfect match."""
    print("Test 1: Position mapping - simple match")

    # Create a simple aligned read with no indels
    header = pysam.AlignmentHeader.from_dict({
        'HD': {'VN': '1.0'},
        'SQ': [{'SN': 'chr1', 'LN': 1000}]
    })

    read = pysam.AlignedSegment(header)
    read.query_sequence = "ATCGATCG"
    read.reference_start = 100
    read.cigarstring = "8M"  # 8 matches

    ref2q_left, ref2q_right = _build_ref2read_maps(read)

    # For a perfect match, both mappings should be identical
    assert ref2q_left[100] == 0, "Position 100 should map to query 0"
    assert ref2q_left[107] == 7, "Position 107 should map to query 7"
    assert ref2q_left == ref2q_right, "Left and right mappings should match for perfect alignment"

    print("  ✅ PASS\n")


def test_position_mapping_with_deletion():
    """Test position mapping for a read with deletion."""
    print("Test 2: Position mapping - deletion")

    # Create read with 2bp deletion: ATCG--CG (-- = deleted from read)
    header = pysam.AlignmentHeader.from_dict({
        'HD': {'VN': '1.0'},
        'SQ': [{'SN': 'chr1', 'LN': 1000}]
    })

    read = pysam.AlignedSegment(header)
    read.query_sequence = "ATCGCG"  # 6 bases
    read.reference_start = 100
    read.cigarstring = "4M2D2M"  # 4 match, 2 deletion, 2 match

    ref2q_left, ref2q_right = _build_ref2read_maps(read)

    # Check mappings around deletion
    assert ref2q_left[103] == 3, "Last base before deletion"
    assert ref2q_left[104] == 3, "Deletion position 1 should map to last base before (left)"
    assert ref2q_left[105] == 3, "Deletion position 2 should map to last base before (left)"
    assert ref2q_right[104] == 4, "Deletion position 1 should map to first base after (right)"
    assert ref2q_right[105] == 4, "Deletion position 2 should map to first base after (right)"
    assert ref2q_left[106] == 4, "First base after deletion"

    print("  ✅ PASS\n")


def test_position_mapping_with_insertion():
    """Test position mapping for a read with insertion."""
    print("Test 3: Position mapping - insertion")

    # Create read with 2bp insertion: ATCGAACG (AA = inserted in read)
    header = pysam.AlignmentHeader.from_dict({
        'HD': {'VN': '1.0'},
        'SQ': [{'SN': 'chr1', 'LN': 1000}]
    })

    read = pysam.AlignedSegment(header)
    read.query_sequence = "ATCGAACG"  # 8 bases
    read.reference_start = 100
    read.cigarstring = "4M2I2M"  # 4 match, 2 insertion, 2 match

    ref2q_left, ref2q_right = _build_ref2read_maps(read)

    # Insertions don't consume reference positions, so ref should skip them
    assert ref2q_left[103] == 3, "Last base before insertion"
    # Query positions 4 and 5 are the insertion - no reference position for them
    assert ref2q_left[104] == 6, "First base after insertion (skips query 4,5)"

    print("  ✅ PASS\n")


def test_quality_filling_with_flanks():
    """Test quality score generation for insertions."""
    print("Test 4: Quality score filling - with flanking data")

    left_qual = np.array([30, 32, 34], dtype=np.uint8)
    right_qual = np.array([36, 38, 40], dtype=np.uint8)

    result = _fill_insertion_quals(5, left_qual, right_qual, insert_qual=30)

    # Should average flanking qualities: mean([30,32,34,36,38,40]) = 35
    expected_mean = int(np.mean(np.concatenate([left_qual, right_qual])))
    assert len(result) == 5, "Should generate 5 quality scores"
    assert np.all(result == expected_mean), f"All qualities should be {expected_mean}"

    print(f"  Generated quality: Q{result[0]} (mean of flanking regions)")
    print("  ✅ PASS\n")


def test_quality_filling_no_flanks():
    """Test quality score generation when no flanking data available."""
    print("Test 5: Quality score filling - no flanking data")

    result = _fill_insertion_quals(3, np.array([]), np.array([]), insert_qual=25)

    assert len(result) == 3, "Should generate 3 quality scores"
    assert np.all(result == 25), "Should use default insert_qual"

    print(f"  Generated quality: Q{result[0]} (default fallback)")
    print("  ✅ PASS\n")


def test_phased_seqs_snp_only():
    """Test SNP-only sequence building (baseline)."""
    print("Test 6: Phased sequences - SNP only")

    split_seq = ["ATC", "G", "GCA", "T", "AAA"]
    hap1_alleles = ["A", "C"]  # Alt alleles for hap1
    hap2_alleles = ["G", "T"]  # Alt alleles for hap2

    hap1, hap2 = make_phased_seqs(split_seq, hap1_alleles, hap2_alleles)

    # Expected: ATC + A + GCA + C + AAA = ATCAGCACAAA
    #           ATC + G + GCA + T + AAA = ATCGGCATAAA
    assert hap1 == "ATCAGCACAAA", f"Hap1 mismatch: {hap1}"
    assert hap2 == "ATCGGCATAAA", f"Hap2 mismatch: {hap2}"

    print(f"  Hap1: {hap1}")
    print(f"  Hap2: {hap2}")
    print("  ✅ PASS\n")


def test_phased_seqs_with_qual_same_length():
    """Test indel-aware sequences with same-length alleles (like SNPs)."""
    print("Test 7: Phased sequences with quality - same length alleles")

    split_seq = ["ATC", "G", "GCA"]
    split_qual = [
        np.array([30, 32, 34], dtype=np.uint8),
        np.array([35], dtype=np.uint8),
        np.array([36, 38, 40], dtype=np.uint8),
    ]
    hap1_alleles = ["A"]  # Same length as "G"
    hap2_alleles = ["T"]

    (hap1, hap1_qual), (hap2, hap2_qual) = make_phased_seqs_with_qual(
        split_seq, split_qual, hap1_alleles, hap2_alleles, insert_qual=30
    )

    assert hap1 == "ATCAGCA", f"Hap1 sequence: {hap1}"
    assert hap2 == "ATCTGCA", f"Hap2 sequence: {hap2}"
    assert len(hap1_qual) == 7, "Hap1 quality length should match sequence"
    assert len(hap2_qual) == 7, "Hap2 quality length should match sequence"

    # Quality should be: [30,32,34] + [35] + [36,38,40]
    expected_qual = np.array([30, 32, 34, 35, 36, 38, 40], dtype=np.uint8)
    assert np.array_equal(hap1_qual, expected_qual), "Quality mismatch"

    print(f"  Hap1: {hap1}")
    print(f"  Qual: {list(hap1_qual)}")
    print("  ✅ PASS\n")


def test_phased_seqs_with_qual_deletion():
    """Test indel-aware sequences with deletion."""
    print("Test 8: Phased sequences with quality - deletion")

    split_seq = ["ATC", "GGG", "GCA"]  # Original has 3bp
    split_qual = [
        np.array([30, 32, 34], dtype=np.uint8),
        np.array([35, 36, 37], dtype=np.uint8),  # 3 qualities for 3bp
        np.array([38, 40, 42], dtype=np.uint8),
    ]
    hap1_alleles = ["G"]  # 1bp - deletion of 2bp
    hap2_alleles = ["GGG"]  # Keep original

    (hap1, hap1_qual), (hap2, hap2_qual) = make_phased_seqs_with_qual(
        split_seq, split_qual, hap1_alleles, hap2_alleles, insert_qual=30
    )

    assert hap1 == "ATCGGCA", f"Hap1 sequence: {hap1}"
    assert hap2 == "ATCGGGGCA", f"Hap2 sequence: {hap2}"

    # Hap1 quality should truncate to first base: [30,32,34] + [35] + [38,40,42]
    assert len(hap1_qual) == 7, f"Hap1 quality length: {len(hap1_qual)}"
    assert hap1_qual[3] == 35, "Should keep first quality from deleted region"

    # Hap2 quality should keep all: [30,32,34] + [35,36,37] + [38,40,42]
    assert len(hap2_qual) == 9, f"Hap2 quality length: {len(hap2_qual)}"

    print(f"  Hap1 (deletion): {hap1} (len={len(hap1)})")
    print(f"  Hap1 qual:       {list(hap1_qual)}")
    print(f"  Hap2 (original): {hap2} (len={len(hap2)})")
    print(f"  Hap2 qual:       {list(hap2_qual)}")
    print("  ✅ PASS\n")


def test_phased_seqs_with_qual_insertion():
    """Test indel-aware sequences with insertion."""
    print("Test 9: Phased sequences with quality - insertion")

    split_seq = ["ATC", "G", "GCA"]  # Original has 1bp
    split_qual = [
        np.array([30, 32, 34], dtype=np.uint8),
        np.array([35], dtype=np.uint8),  # 1 quality for 1bp
        np.array([38, 40, 42], dtype=np.uint8),
    ]
    hap1_alleles = ["GGG"]  # 3bp - insertion of 2bp
    hap2_alleles = ["G"]  # Keep original

    (hap1, hap1_qual), (hap2, hap2_qual) = make_phased_seqs_with_qual(
        split_seq, split_qual, hap1_alleles, hap2_alleles, insert_qual=30
    )

    assert hap1 == "ATCGGGGCA", f"Hap1 sequence: {hap1}"
    assert hap2 == "ATCGGCA", f"Hap2 sequence: {hap2}"

    # Hap1 quality should add 2 extra scores: [30,32,34] + [35, X, X] + [38,40,42]
    # where X is computed from flanking regions
    assert len(hap1_qual) == 9, f"Hap1 quality length: {len(hap1_qual)}"
    assert hap1_qual[3] == 35, "Original quality preserved"
    # Extra qualities should be mean of [30,32,34,38,40,42]
    expected_extra = int(np.mean(np.array([30, 32, 34, 38, 40, 42])))
    assert hap1_qual[4] == expected_extra, f"Inserted quality should be ~{expected_extra}"

    # Hap2 quality should be original: [30,32,34] + [35] + [38,40,42]
    assert len(hap2_qual) == 7, f"Hap2 quality length: {len(hap2_qual)}"

    print(f"  Hap1 (insertion): {hap1} (len={len(hap1)})")
    print(f"  Hap1 qual:        {list(hap1_qual)}")
    print(f"  Hap2 (original):  {hap2} (len={len(hap2)})")
    print(f"  Hap2 qual:        {list(hap2_qual)}")
    print("  ✅ PASS\n")


def test_multi_sample_sequences():
    """Test multi-sample sequence generation."""
    print("Test 10: Multi-sample sequences with quality")

    split_seq = ["AT", "G", "GC"]
    split_qual = [
        np.array([30, 32], dtype=np.uint8),
        np.array([35], dtype=np.uint8),
        np.array([38, 40], dtype=np.uint8),
    ]
    # 3 unique haplotypes across samples
    allele_combos = [
        ["A"],  # Hap1
        ["G"],  # Hap2
        ["T"],  # Hap3
    ]

    result = make_multi_seqs_with_qual(split_seq, split_qual, allele_combos, insert_qual=30)

    assert len(result) == 3, "Should generate 3 haplotypes"
    assert result[0][0] == "ATAGC", f"Hap1: {result[0][0]}"
    assert result[1][0] == "ATGGC", f"Hap2: {result[1][0]}"
    assert result[2][0] == "ATTGC", f"Hap3: {result[2][0]}"

    # All should have same quality length (5)
    assert all(len(qual) == 5 for seq, qual in result), "All quality arrays should be length 5"

    print(f"  Hap1: {result[0][0]} - {list(result[0][1])}")
    print(f"  Hap2: {result[1][0]} - {list(result[1][1])}")
    print(f"  Hap3: {result[2][0]} - {list(result[2][1])}")
    print("  ✅ PASS\n")


def run_all_tests():
    """Run all correctness tests."""
    print("=" * 70)
    print("WASP2 INDEL IMPLEMENTATION - CORRECTNESS TESTS")
    print("=" * 70)
    print()

    tests = [
        test_position_mapping_simple_match,
        test_position_mapping_with_deletion,
        test_position_mapping_with_insertion,
        test_quality_filling_with_flanks,
        test_quality_filling_no_flanks,
        test_phased_seqs_snp_only,
        test_phased_seqs_with_qual_same_length,
        test_phased_seqs_with_qual_deletion,
        test_phased_seqs_with_qual_insertion,
        test_multi_sample_sequences,
    ]

    passed = 0
    failed = 0

    for test in tests:
        try:
            test()
            passed += 1
        except AssertionError as e:
            print(f"  ❌ FAIL: {e}\n")
            failed += 1
        except Exception as e:
            print(f"  ❌ ERROR: {e}\n")
            failed += 1

    print("=" * 70)
    print(f"RESULTS: {passed} passed, {failed} failed")
    print("=" * 70)

    if failed == 0:
        print("✅ ALL TESTS PASSED - Code is correct!")
        print()
        print("Next step: Run performance benchmarks")
        print("  python benchmark_indels.py")
        return 0
    else:
        print("❌ SOME TESTS FAILED - Fix errors before benchmarking")
        return 1


if __name__ == "__main__":
    exit(run_all_tests())
