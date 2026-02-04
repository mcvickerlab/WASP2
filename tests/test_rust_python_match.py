"""
Rust vs Python parity tests for INDEL algorithms.

Verifies that Python's make_phased_seqs_with_qual and _build_ref2read_maps
produce identical results to the Rust unit tests in multi_sample.rs.
"""

import sys
from pathlib import Path

import numpy as np
import pysam

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from mapping.remap_utils import _build_ref2read_maps, make_phased_seqs_with_qual


class TestRustPythonParity:
    """Verify Python produces the same outputs as Rust test cases in multi_sample.rs."""

    def test_deletion_substitution(self):
        """Match Rust test_cigar_aware_deletion_substitution.

        Sequence: AAACGAAAA (9 bases)
        Variant at pos 3: ACG -> A (delete CG)
        Expected alt output: AAAAAAA (7 bases)
        """
        split_seq = ["AAA", "CGA", "AAA"]
        split_qual = [np.array([30, 30, 30]), np.array([30, 30, 30]), np.array([30, 30, 30])]
        hap1_alleles = ["A"]  # alt allele (deletion: CGA -> A)
        hap2_alleles = ["CGA"]  # keep original read content

        (seq1, qual1), (seq2, qual2) = make_phased_seqs_with_qual(
            split_seq, split_qual, hap1_alleles, hap2_alleles
        )

        assert seq1 == "AAAAAAA", f"Deletion alt: expected AAAAAAA, got {seq1}"
        assert seq2 == "AAACGAAAA", f"Deletion ref: expected AAACGAAAA, got {seq2}"

    def test_insertion_substitution(self):
        """Match Rust test_cigar_aware_insertion_substitution.

        Sequence: AAAAAAA (7 bases)
        Variant at pos 3: A -> ACGT (insert CGT)
        Expected alt output: AAAACGTAAA (10 bases)
        """
        split_seq = ["AAA", "A", "AAA"]
        split_qual = [np.array([30, 30, 30]), np.array([30]), np.array([30, 30, 30])]
        hap1_alleles = ["ACGT"]  # alt allele (insertion)
        hap2_alleles = ["A"]  # ref allele

        (seq1, qual1), (seq2, qual2) = make_phased_seqs_with_qual(
            split_seq, split_qual, hap1_alleles, hap2_alleles
        )

        assert seq1 == "AAAACGTAAA", f"Insertion alt: expected AAAACGTAAA, got {seq1}"
        assert seq2 == "AAAAAAA", f"Insertion ref: expected AAAAAAA, got {seq2}"

    def test_multiple_snps(self):
        """Match Rust test_cigar_aware_multiple_variants.

        Sequence: AAAAAAAAA (9 bases)
        Variant at pos 2: A -> G
        Variant at pos 6: A -> T
        Expected alt output: AAGAAATAA
        """
        split_seq = ["AA", "A", "AAA", "A", "AA"]
        split_qual = [
            np.array([30, 30]),
            np.array([30]),
            np.array([30, 30, 30]),
            np.array([30]),
            np.array([30, 30]),
        ]
        hap1_alleles = ["G", "T"]  # both alt
        hap2_alleles = ["A", "A"]  # both ref

        (seq1, qual1), (seq2, qual2) = make_phased_seqs_with_qual(
            split_seq, split_qual, hap1_alleles, hap2_alleles
        )

        assert seq1 == "AAGAAATAA", f"Multi-SNP alt: expected AAGAAATAA, got {seq1}"
        assert seq2 == "AAAAAAAAA", f"Multi-SNP ref: expected AAAAAAAAA, got {seq2}"

    def test_cigar_aware_deletion_mapping(self):
        """Match Rust test_cigar_aware_with_deletion_in_cigar.

        Read: AAAAABBBBB (10 bp) with CIGAR 5M2D5M (deletion at ref 5-6)
        Ref pos 7 should map to query pos 5 (not 7) due to the 2bp deletion.
        """
        header = pysam.AlignmentHeader.from_dict(
            {"HD": {"VN": "1.0"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
        )
        read = pysam.AlignedSegment(header)
        read.query_sequence = "AAAAABBBBB"
        read.reference_start = 0
        read.cigarstring = "5M2D5M"
        read.query_qualities = pysam.qualitystring_to_array("?" * 10)

        ref2q_left, ref2q_right = _build_ref2read_maps(read)

        assert ref2q_left.get(0, -1) == 0, "ref pos 0 -> query pos 0"
        assert ref2q_left.get(4, -1) == 4, "ref pos 4 -> query pos 4"
        assert ref2q_left.get(7, -1) == 5, (
            "ref pos 7 -> query pos 5 (key test: accounts for 2bp deletion)"
        )
        assert ref2q_left.get(8, -1) == 6, "ref pos 8 -> query pos 6"
