#!/usr/bin/env python3
"""
Direct comparison: Verify Rust and Python INDEL algorithms match.
Uses the SAME test cases as Rust unit tests in multi_sample.rs
"""
import sys
from pathlib import Path

import numpy as np
import pysam
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from mapping.remap_utils import make_phased_seqs_with_qual, _build_ref2read_maps


@pytest.fixture
def pysam_header():
    """Create a pysam alignment header for tests."""
    return pysam.AlignmentHeader.from_dict({
        'HD': {'VN': '1.0'},
        'SQ': [{'SN': 'chr1', 'LN': 1000}]
    })


class TestRustPythonMatch:
    """
    These are the EXACT same test cases from Rust: multi_sample.rs lines 960-1097
    """

    def test_deletion_substitution_alt(self):
        """
        Rust test_cigar_aware_deletion_substitution:
        Sequence: AAACGAAAA (9 bases)
        Variant at pos 3: ACG -> A (delete CG)
        Expected output: AAAAAAA (7 bases)
        """
        # [before, variant_region, after]
        split_seq = ["AAA", "CGA", "AAA"]
        split_qual = [np.array([30, 30, 30]), np.array([30, 30, 30]), np.array([30, 30, 30])]
        hap1_alleles = ["A"]  # alt allele (deletion: CGA -> A)
        hap2_alleles = ["CGA"]  # keep original read content

        (seq1, qual1), (seq2, qual2) = make_phased_seqs_with_qual(
            split_seq, split_qual, hap1_alleles, hap2_alleles
        )

        assert seq1 == "AAAAAAA", "Replace 3bp region with A"

    def test_deletion_substitution_ref(self):
        """
        Rust test_cigar_aware_deletion_substitution:
        Keep original 3bp region (ref allele)
        """
        split_seq = ["AAA", "CGA", "AAA"]
        split_qual = [np.array([30, 30, 30]), np.array([30, 30, 30]), np.array([30, 30, 30])]
        hap1_alleles = ["A"]
        hap2_alleles = ["CGA"]

        (seq1, qual1), (seq2, qual2) = make_phased_seqs_with_qual(
            split_seq, split_qual, hap1_alleles, hap2_alleles
        )

        assert seq2 == "AAACGAAAA", "Keep original 3bp region"

    def test_insertion_substitution_alt(self):
        """
        Rust test_cigar_aware_insertion_substitution:
        Sequence: AAAAAAA (7 bases)
        Variant at pos 3: A -> ACGT (insert CGT)
        Expected output: AAAACGTAAA (10 bases)
        """
        split_seq = ["AAA", "A", "AAA"]
        split_qual = [np.array([30, 30, 30]), np.array([30]), np.array([30, 30, 30])]
        hap1_alleles = ["ACGT"]  # alt allele (insertion)
        hap2_alleles = ["A"]  # ref allele

        (seq1, qual1), (seq2, qual2) = make_phased_seqs_with_qual(
            split_seq, split_qual, hap1_alleles, hap2_alleles
        )

        assert seq1 == "AAAACGTAAA", "A->ACGT at pos 3"

    def test_insertion_substitution_ref(self):
        """
        Rust test_cigar_aware_insertion_substitution:
        Keep A at pos 3 (ref allele)
        """
        split_seq = ["AAA", "A", "AAA"]
        split_qual = [np.array([30, 30, 30]), np.array([30]), np.array([30, 30, 30])]
        hap1_alleles = ["ACGT"]
        hap2_alleles = ["A"]

        (seq1, qual1), (seq2, qual2) = make_phased_seqs_with_qual(
            split_seq, split_qual, hap1_alleles, hap2_alleles
        )

        assert seq2 == "AAAAAAA", "Keep A at pos 3"

    def test_multiple_snps_alt_alt(self):
        """
        Rust test_cigar_aware_multiple_variants:
        Sequence: AAAAAAAAA (9 bases)
        Variant at pos 2: A -> G
        Variant at pos 6: A -> T
        Expected output: AAGAAATAA
        """
        # Two variants: [before, v1, between, v2, after]
        split_seq = ["AA", "A", "AAA", "A", "AA"]
        split_qual = [
            np.array([30, 30]),
            np.array([30]),
            np.array([30, 30, 30]),
            np.array([30]),
            np.array([30, 30])
        ]
        hap1_alleles = ["G", "T"]  # both alt
        hap2_alleles = ["A", "A"]  # both ref

        (seq1, qual1), (seq2, qual2) = make_phased_seqs_with_qual(
            split_seq, split_qual, hap1_alleles, hap2_alleles
        )

        assert seq1 == "AAGAAATAA", "Both variants applied"

    def test_multiple_snps_ref_ref(self):
        """
        Rust test_cigar_aware_multiple_variants:
        No variants applied (ref/ref)
        """
        split_seq = ["AA", "A", "AAA", "A", "AA"]
        split_qual = [
            np.array([30, 30]),
            np.array([30]),
            np.array([30, 30, 30]),
            np.array([30]),
            np.array([30, 30])
        ]
        hap1_alleles = ["G", "T"]
        hap2_alleles = ["A", "A"]

        (seq1, qual1), (seq2, qual2) = make_phased_seqs_with_qual(
            split_seq, split_qual, hap1_alleles, hap2_alleles
        )

        assert seq2 == "AAAAAAAAA", "No variants applied"


class TestCigarAwareDeletionMapping:
    """
    Rust test_cigar_aware_with_deletion_in_cigar:
    Read: AAAAABBBBB (10 bp) with CIGAR 5M2D5M (deletion at ref 5-6)
    Variant at ref pos 7: B -> X
    Expected: AAAAAXBBBB (X at query pos 5, not 7!)

    This tests that CIGAR-aware position mapping correctly handles deletions.
    """

    @pytest.fixture
    def read_with_deletion(self, pysam_header):
        """Create a pysam read with deletion for CIGAR-aware tests."""
        read = pysam.AlignedSegment(pysam_header)
        read.query_sequence = "AAAAABBBBB"
        read.reference_start = 0
        read.cigarstring = "5M2D5M"  # 5 match, 2 deletion, 5 match
        read.query_qualities = pysam.qualitystring_to_array("?" * 10)
        return read

    @pytest.fixture
    def ref2q_left(self, read_with_deletion):
        """Build the left position map using Python's CIGAR-aware function."""
        ref2q_left, _ = _build_ref2read_maps(read_with_deletion)
        return ref2q_left

    def test_ref_pos_0_maps_to_query_pos_0(self, ref2q_left):
        """Reference position 0 should map to query position 0."""
        assert ref2q_left.get(0, -1) == 0

    def test_ref_pos_4_maps_to_query_pos_4(self, ref2q_left):
        """Reference position 4 should map to query position 4."""
        assert ref2q_left.get(4, -1) == 4

    def test_ref_pos_7_maps_to_query_pos_5(self, ref2q_left):
        """
        Key test: ref 7 should map to query 5 due to 2bp deletion.
        Positions 5-6 are deleted in ref, so ref 7 should map to query 5.
        """
        assert ref2q_left.get(7, -1) == 5

    def test_ref_pos_8_maps_to_query_pos_6(self, ref2q_left):
        """Reference position 8 should map to query position 6."""
        assert ref2q_left.get(8, -1) == 6
