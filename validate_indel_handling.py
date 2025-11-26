#!/usr/bin/env python3
"""
Validation script for WASP2 indel handling.

Tests Python implementation against ground truth expected behavior.
See INDEL_VALIDATION_PLAN.md for specifications.
"""

import sys
import numpy as np
from typing import Dict, Tuple, List
import pysam
from pysam import AlignedSegment, AlignmentFile
import polars as pl
import tempfile
import os

# Import functions to test
sys.path.insert(0, 'src')
from mapping.remap_utils import (
    _build_ref2read_maps,
    _fill_insertion_quals,
    make_phased_seqs_with_qual,
    get_read_het_data
)


class Colors:
    """ANSI color codes for terminal output"""
    GREEN = '\033[92m'
    RED = '\033[91m'
    BLUE = '\033[94m'
    YELLOW = '\033[93m'
    RESET = '\033[0m'
    BOLD = '\033[1m'


def create_mock_read(
    sequence: str,
    qualities: List[int],
    reference_start: int,
    cigarstring: str,
    name: str = "test_read"
) -> AlignedSegment:
    """Create a mock pysam AlignedSegment for testing.

    Args:
        sequence: Read sequence
        qualities: Phred quality scores
        reference_start: 0-based reference start position
        cigarstring: CIGAR string (e.g., "10M2I10M")
        name: Read name

    Returns:
        AlignedSegment object
    """
    # Create temporary BAM to get proper AlignedSegment
    header = {'HD': {'VN': '1.0'},
              'SQ': [{'SN': 'chr1', 'LN': 1000}]}

    with tempfile.NamedTemporaryFile(suffix='.bam', delete=False) as f:
        bam_path = f.name

    try:
        with pysam.AlignmentFile(bam_path, 'wb', header=header) as bam_out:
            read = pysam.AlignedSegment(bam_out.header)
            read.query_name = name
            read.query_sequence = sequence
            read.query_qualities = qualities
            read.reference_id = 0
            read.reference_start = reference_start
            read.cigarstring = cigarstring
            read.mapping_quality = 60
            read.flag = 99  # Paired, proper pair, first in pair
            bam_out.write(read)

        # Read it back to get proper object
        with pysam.AlignmentFile(bam_path, 'rb') as bam_in:
            read = next(bam_in)
            return read
    finally:
        if os.path.exists(bam_path):
            os.unlink(bam_path)


def test_position_mapping_snp():
    """Test 1: Position mapping for SNPs (baseline)"""
    print(f"\n{Colors.BLUE}Test 1: Position mapping for SNPs{Colors.RESET}")

    # Create read with simple match
    sequence = "ACGTACGTACGT"
    qualities = [30] * 12
    read = create_mock_read(sequence, qualities, reference_start=100, cigarstring="12M")

    # Build position maps
    ref2q_left, ref2q_right = _build_ref2read_maps(read)

    # For pure matches, left and right should be identical
    assert ref2q_left == ref2q_right, "SNP mapping: left and right maps should match"

    # Check specific positions
    assert ref2q_left[100] == 0, f"Ref pos 100 should map to read pos 0, got {ref2q_left[100]}"
    assert ref2q_left[105] == 5, f"Ref pos 105 should map to read pos 5, got {ref2q_left[105]}"
    assert ref2q_left[111] == 11, f"Ref pos 111 should map to read pos 11, got {ref2q_left[111]}"

    print(f"  {Colors.GREEN}✅ PASS{Colors.RESET}: SNP position mapping correct")
    return True


def test_position_mapping_insertion():
    """Test 2: Position mapping with insertion"""
    print(f"\n{Colors.BLUE}Test 2: Position mapping with 3bp insertion{Colors.RESET}")

    # Read with insertion: 5 matches, 3bp insert, 4 matches
    # Reference: positions 100-108 (9 bases)
    # Read: 12 bases (5M + 3I + 4M)
    sequence = "ACGTAXXXCGTA"  # XXX = insertion
    qualities = [30] * 12
    read = create_mock_read(sequence, qualities, reference_start=100, cigarstring="5M3I4M")

    ref2q_left, ref2q_right = _build_ref2read_maps(read)

    # Before insertion
    assert ref2q_left[100] == 0, "Ref 100 -> read 0"
    assert ref2q_left[104] == 4, "Ref 104 -> read 4"

    # After insertion (ref pos 105 starts after the 3bp insert in the read)
    # Ref: A C G T A | C G T A
    # Pos: 100 101 102 103 104 | 105 106 107 108
    # Read: A C G T A X X X | C G T A
    # Pos:  0 1 2 3 4 5 6 7 | 8 9 10 11
    assert ref2q_left[105] == 8, f"Ref 105 should map to read 8, got {ref2q_left[105]}"
    assert ref2q_left[108] == 11, f"Ref 108 should map to read 11, got {ref2q_left[108]}"

    print(f"  {Colors.GREEN}✅ PASS{Colors.RESET}: Insertion position mapping correct")
    return True


def test_position_mapping_deletion():
    """Test 3: Position mapping with deletion"""
    print(f"\n{Colors.BLUE}Test 3: Position mapping with 3bp deletion{Colors.RESET}")

    # Read with deletion: 5 matches, 3bp deletion, 4 matches
    # Reference: positions 100-111 (12 bases)
    # Read: 9 bases (5M + 3D + 4M)
    sequence = "ACGTACGTA"
    qualities = [30] * 9
    read = create_mock_read(sequence, qualities, reference_start=100, cigarstring="5M3D4M")

    ref2q_left, ref2q_right = _build_ref2read_maps(read)

    # Before deletion
    assert ref2q_left[100] == 0, "Ref 100 -> read 0"
    assert ref2q_left[104] == 4, "Ref 104 -> read 4"

    # Deletion region (ref positions 105, 106, 107 are deleted in read)
    # These should have flanking positions
    assert ref2q_left[105] == 4, f"Deleted ref 105 left flank should be 4, got {ref2q_left[105]}"
    assert ref2q_right[105] == 5, f"Deleted ref 105 right flank should be 5, got {ref2q_right[105]}"

    assert ref2q_left[106] == 4, "Deleted ref 106 left flank should be 4"
    assert ref2q_right[106] == 5, "Deleted ref 106 right flank should be 5"

    # After deletion
    assert ref2q_left[108] == 5, f"Ref 108 should map to read 5, got {ref2q_left[108]}"
    assert ref2q_left[111] == 8, f"Ref 111 should map to read 8, got {ref2q_left[111]}"

    print(f"  {Colors.GREEN}✅ PASS{Colors.RESET}: Deletion position mapping correct (left/right flanking)")
    return True


def test_quality_score_insertion():
    """Test 4: Quality score generation for insertions"""
    print(f"\n{Colors.BLUE}Test 4: Quality score generation for insertions{Colors.RESET}")

    # Test with flanking qualities
    left_qual = np.array([35, 30, 28], dtype=np.uint8)
    right_qual = np.array([32, 30], dtype=np.uint8)

    result = _fill_insertion_quals(3, left_qual, right_qual, insert_qual=30)

    # Should average: (35+30+28+32+30)/5 = 31
    expected_avg = int(np.mean(np.concatenate([left_qual, right_qual])))
    assert len(result) == 3, f"Should generate 3 qualities, got {len(result)}"
    assert all(q == expected_avg for q in result), f"All qualities should be {expected_avg}, got {result}"

    # Test with no flanking data
    result_no_flank = _fill_insertion_quals(5, np.array([]), np.array([]), insert_qual=30)
    assert len(result_no_flank) == 5, "Should generate 5 qualities"
    assert all(q == 30 for q in result_no_flank), "Should use default Q30"

    print(f"  {Colors.GREEN}✅ PASS{Colors.RESET}: Quality score generation correct")
    return True


def test_sequence_transformation_insertion():
    """Test 5: Full sequence transformation with 3bp insertion"""
    print(f"\n{Colors.BLUE}Test 5: Sequence transformation with 3bp insertion{Colors.RESET}")

    # Original read: ACGTACGTACGT (12 bases)
    # Insertion at position 5 (0-based): insert "GGG" after "ACGTA"
    # Expected:
    #   Hap1 (with insertion): ACGTA + AGGG + CGTACGT = ACGTAAGGGCGTACGT (16 bases)
    #   Hap2 (reference):      ACGTA + A    + CGTACGT = ACGTAACGTACGT (13 bases)

    sequence = "ACGTAACGTACGT"
    qualities = np.array([30] * 13, dtype=np.uint8)

    # Split sequence manually (simulating get_read_het_data output)
    split_seq = ["ACGTA", "A", "CGTACGT"]
    split_qual = [
        np.array([30, 30, 30, 30, 30], dtype=np.uint8),
        np.array([30], dtype=np.uint8),
        np.array([30, 30, 30, 30, 30, 30, 30], dtype=np.uint8)
    ]

    # Haplotypes: one with insertion, one without
    hap1_alleles = ["AGGG"]  # 3bp insertion
    hap2_alleles = ["A"]     # reference

    (hap1_seq, hap1_qual), (hap2_seq, hap2_qual) = make_phased_seqs_with_qual(
        split_seq, split_qual, hap1_alleles, hap2_alleles, insert_qual=30
    )

    # Check sequences
    assert hap1_seq == "ACGTAAGGGCGTACGT", f"Hap1 should be 'ACGTAAGGGCGTACGT', got '{hap1_seq}'"
    assert hap2_seq == "ACGTAACGTACGT", f"Hap2 should be 'ACGTAACGTACGT', got '{hap2_seq}'"

    # Check lengths
    assert len(hap1_seq) == 16, f"Hap1 should be 16 bases, got {len(hap1_seq)}"
    assert len(hap2_seq) == 13, f"Hap2 should be 13 bases, got {len(hap2_seq)}"

    # Check quality lengths match sequence lengths
    assert len(hap1_qual) == len(hap1_seq), f"Hap1 quality length {len(hap1_qual)} != seq length {len(hap1_seq)}"
    assert len(hap2_qual) == len(hap2_seq), f"Hap2 quality length {len(hap2_qual)} != seq length {len(hap2_seq)}"

    # Check inserted qualities are reasonable (>= 20)
    assert all(q >= 20 for q in hap1_qual), "All hap1 qualities should be >= Q20"

    print(f"  {Colors.GREEN}✅ PASS{Colors.RESET}: Insertion transformation correct")
    return True


def test_sequence_transformation_deletion():
    """Test 6: Full sequence transformation with 3bp deletion"""
    print(f"\n{Colors.BLUE}Test 6: Sequence transformation with 3bp deletion{Colors.RESET}")

    # Original read: ACGGAAACGTACGT (14 bases)
    # Deletion: delete "AAA" from "GAAA"
    # Split segments represent: [before_variant, variant_bases, after_variant]
    # Expected:
    #   Hap1 (with deletion):  ACGG + G    + CGTACGT = ACGGGCGTACGT (12 bases)
    #   Hap2 (reference):      ACGG + GAAA + CGTACGT = ACGGGAAACGTACGT (14 bases)

    split_seq = ["ACGG", "GAAA", "CGTACGT"]
    split_qual = [
        np.array([30, 30, 30, 30], dtype=np.uint8),
        np.array([30, 30, 30, 30], dtype=np.uint8),
        np.array([30, 30, 30, 30, 30, 30, 30], dtype=np.uint8)
    ]

    # Haplotypes: one with deletion, one without
    hap1_alleles = ["G"]      # 3bp deletion (only G left)
    hap2_alleles = ["GAAA"]   # reference (all 4 bases)

    (hap1_seq, hap1_qual), (hap2_seq, hap2_qual) = make_phased_seqs_with_qual(
        split_seq, split_qual, hap1_alleles, hap2_alleles, insert_qual=30
    )

    # Check sequences
    assert hap1_seq == "ACGGGCGTACGT", f"Hap1 should be 'ACGGGCGTACGT', got '{hap1_seq}'"
    assert hap2_seq == "ACGGGAAACGTACGT", f"Hap2 should be 'ACGGGAAACGTACGT', got '{hap2_seq}'"

    # Check lengths
    assert len(hap1_seq) == 12, f"Hap1 should be 12 bases, got {len(hap1_seq)}"
    assert len(hap2_seq) == 15, f"Hap2 should be 15 bases, got {len(hap2_seq)}"

    # Check quality lengths
    assert len(hap1_qual) == 12, f"Hap1 quality should be 12, got {len(hap1_qual)}"
    assert len(hap2_qual) == 15, f"Hap2 quality should be 15, got {len(hap2_qual)}"

    # Check all qualities are preserved
    assert all(q == 30 for q in hap1_qual), "All hap1 qualities should be 30"
    assert all(q == 30 for q in hap2_qual), "All hap2 qualities should be 30"

    print(f"  {Colors.GREEN}✅ PASS{Colors.RESET}: Deletion transformation correct")
    return True


def test_edge_case_multiple_variants():
    """Test 7: Multiple variants on same read"""
    print(f"\n{Colors.BLUE}Test 7: Multiple variants on same read{Colors.RESET}")

    # Read with 2 variants: 1 SNP and 1 insertion
    # Original: ACGTACGTACGT
    # Variant 1 at pos 2: G->T (SNP)
    # Variant 2 at pos 7: T->TGG (2bp insertion)

    split_seq = ["AC", "G", "TACG", "T", "ACGT"]
    split_qual = [
        np.array([30, 30], dtype=np.uint8),
        np.array([30], dtype=np.uint8),
        np.array([30, 30, 30, 30], dtype=np.uint8),
        np.array([30], dtype=np.uint8),
        np.array([30, 30, 30, 30], dtype=np.uint8)
    ]

    hap1_alleles = ["T", "TGG"]    # SNP + insertion
    hap2_alleles = ["G", "T"]      # reference

    (hap1_seq, hap1_qual), (hap2_seq, hap2_qual) = make_phased_seqs_with_qual(
        split_seq, split_qual, hap1_alleles, hap2_alleles, insert_qual=30
    )

    # Expected:
    # Hap1: AC + T + TACG + TGG + ACGT = ACTTACGTGGACGT (14 bases)
    # Hap2: AC + G + TACG + T   + ACGT = ACGTACGTACGT (12 bases)

    assert hap1_seq == "ACTTACGTGGACGT", f"Hap1 wrong: {hap1_seq}"
    assert hap2_seq == "ACGTACGTACGT", f"Hap2 wrong: {hap2_seq}"
    assert len(hap1_qual) == 14, "Hap1 quality length wrong"
    assert len(hap2_qual) == 12, "Hap2 quality length wrong"

    print(f"  {Colors.GREEN}✅ PASS{Colors.RESET}: Multiple variants handled correctly")
    return True


def main():
    """Run all validation tests"""
    print(f"\n{Colors.BOLD}{'='*60}{Colors.RESET}")
    print(f"{Colors.BOLD}WASP2 Indel Validation Suite{Colors.RESET}")
    print(f"{Colors.BOLD}{'='*60}{Colors.RESET}")

    tests = [
        test_position_mapping_snp,
        test_position_mapping_insertion,
        test_position_mapping_deletion,
        test_quality_score_insertion,
        test_sequence_transformation_insertion,
        test_sequence_transformation_deletion,
        test_edge_case_multiple_variants,
    ]

    passed = 0
    failed = 0

    for test_func in tests:
        try:
            if test_func():
                passed += 1
        except AssertionError as e:
            failed += 1
            print(f"  {Colors.RED}❌ FAIL{Colors.RESET}: {e}")
        except Exception as e:
            failed += 1
            print(f"  {Colors.RED}❌ ERROR{Colors.RESET}: {e}")

    # Summary
    print(f"\n{Colors.BOLD}{'='*60}{Colors.RESET}")
    print(f"{Colors.BOLD}Summary:{Colors.RESET}")
    print(f"  {Colors.GREEN}Passed: {passed}{Colors.RESET}")
    print(f"  {Colors.RED if failed > 0 else Colors.GREEN}Failed: {failed}{Colors.RESET}")
    print(f"  Total:  {passed + failed}")
    print(f"{Colors.BOLD}{'='*60}{Colors.RESET}\n")

    if failed == 0:
        print(f"{Colors.GREEN}{Colors.BOLD}✅ ALL TESTS PASSED{Colors.RESET}")
        print(f"{Colors.GREEN}Python implementation validated - ready for Rust implementation{Colors.RESET}\n")
        return 0
    else:
        print(f"{Colors.RED}{Colors.BOLD}❌ TESTS FAILED{Colors.RESET}")
        print(f"{Colors.RED}Fix Python implementation before proceeding to Rust{Colors.RESET}\n")
        return 1


if __name__ == "__main__":
    sys.exit(main())
