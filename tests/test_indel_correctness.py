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
import pytest

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from mapping.remap_utils import (
    _build_ref2read_maps,
    _fill_insertion_quals,
    make_phased_seqs,
    make_phased_seqs_with_qual,
    make_multi_seqs_with_qual
)


@pytest.fixture
def alignment_header():
    """Create a pysam alignment header for test reads."""
    return pysam.AlignmentHeader.from_dict({
        'HD': {'VN': '1.0'},
        'SQ': [{'SN': 'chr1', 'LN': 1000}]
    })


def test_position_mapping_simple_match(alignment_header):
    """Test position mapping for a simple perfect match."""
    # Create a simple aligned read with no indels
    read = pysam.AlignedSegment(alignment_header)
    read.query_sequence = "ATCGATCG"
    read.reference_start = 100
    read.cigarstring = "8M"  # 8 matches

    ref2q_left, ref2q_right = _build_ref2read_maps(read)

    # For a perfect match, both mappings should be identical
    assert ref2q_left[100] == 0, "Position 100 should map to query 0"
    assert ref2q_left[107] == 7, "Position 107 should map to query 7"
    assert ref2q_left == ref2q_right, "Left and right mappings should match for perfect alignment"


def test_position_mapping_with_deletion(alignment_header):
    """Test position mapping for a read with deletion."""
    # Create read with 2bp deletion: ATCG--CG (-- = deleted from read)
    read = pysam.AlignedSegment(alignment_header)
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


def test_position_mapping_with_insertion(alignment_header):
    """Test position mapping for a read with insertion."""
    # Create read with 2bp insertion: ATCGAACG (AA = inserted in read)
    read = pysam.AlignedSegment(alignment_header)
    read.query_sequence = "ATCGAACG"  # 8 bases
    read.reference_start = 100
    read.cigarstring = "4M2I2M"  # 4 match, 2 insertion, 2 match

    ref2q_left, ref2q_right = _build_ref2read_maps(read)

    # Insertions don't consume reference positions, so ref should skip them
    assert ref2q_left[103] == 3, "Last base before insertion"
    # Query positions 4 and 5 are the insertion - no reference position for them
    assert ref2q_left[104] == 6, "First base after insertion (skips query 4,5)"


def test_quality_filling_with_flanks():
    """Test quality score generation for insertions."""
    left_qual = np.array([30, 32, 34], dtype=np.uint8)
    right_qual = np.array([36, 38, 40], dtype=np.uint8)

    result = _fill_insertion_quals(5, left_qual, right_qual, insert_qual=30)

    # Should average flanking qualities: mean([30,32,34,36,38,40]) = 35
    expected_mean = int(np.mean(np.concatenate([left_qual, right_qual])))
    assert len(result) == 5, "Should generate 5 quality scores"
    assert np.all(result == expected_mean), f"All qualities should be {expected_mean}"


def test_quality_filling_no_flanks():
    """Test quality score generation when no flanking data available."""
    result = _fill_insertion_quals(3, np.array([]), np.array([]), insert_qual=25)

    assert len(result) == 3, "Should generate 3 quality scores"
    assert np.all(result == 25), "Should use default insert_qual"


def test_phased_seqs_snp_only():
    """Test SNP-only sequence building (baseline)."""
    split_seq = ["ATC", "G", "GCA", "T", "AAA"]
    hap1_alleles = ["A", "C"]  # Alt alleles for hap1
    hap2_alleles = ["G", "T"]  # Alt alleles for hap2

    hap1, hap2 = make_phased_seqs(split_seq, hap1_alleles, hap2_alleles)

    # Expected: ATC + A + GCA + C + AAA = ATCAGCACAAA
    #           ATC + G + GCA + T + AAA = ATCGGCATAAA
    assert hap1 == "ATCAGCACAAA", f"Hap1 mismatch: {hap1}"
    assert hap2 == "ATCGGCATAAA", f"Hap2 mismatch: {hap2}"


def test_phased_seqs_with_qual_same_length():
    """Test indel-aware sequences with same-length alleles (like SNPs)."""
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


def test_phased_seqs_with_qual_deletion():
    """Test indel-aware sequences with deletion."""
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


def test_phased_seqs_with_qual_insertion():
    """Test indel-aware sequences with insertion."""
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


def test_multi_sample_sequences():
    """Test multi-sample sequence generation."""
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
