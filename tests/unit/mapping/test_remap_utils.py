"""
Unit tests for the mapping.remap_utils module.

Tests cover:
- Paired read generators
- Alignment position generation
- Reference to read position mapping (indel support)
- Phased sequence creation
- Quality score handling for indels
"""

from unittest.mock import MagicMock, PropertyMock
from typing import List, Tuple

import numpy as np
import polars as pl
import pytest


class TestPairedReadGen:
    """Tests for paired read generators."""

    def test_yields_proper_pairs(self):
        """Test that proper pairs are yielded correctly."""
        from mapping.remap_utils import paired_read_gen

        # Create mock BAM
        mock_bam = MagicMock()

        # Create mock reads
        read1 = MagicMock()
        read1.query_name = "read1"
        read1.is_proper_pair = True
        read1.is_secondary = False
        read1.is_supplementary = False
        read1.is_read1 = True

        read2 = MagicMock()
        read2.query_name = "read1"
        read2.is_proper_pair = True
        read2.is_secondary = False
        read2.is_supplementary = False
        read2.is_read1 = False

        mock_bam.fetch.return_value = iter([read1, read2])

        pairs = list(paired_read_gen(mock_bam))

        assert len(pairs) == 1
        assert pairs[0] == (read1, read2)

    def test_filters_improper_pairs(self):
        """Test that improper pairs are filtered."""
        from mapping.remap_utils import paired_read_gen

        mock_bam = MagicMock()

        read1 = MagicMock()
        read1.query_name = "read1"
        read1.is_proper_pair = False  # Improper
        read1.is_secondary = False
        read1.is_supplementary = False

        mock_bam.fetch.return_value = iter([read1])

        pairs = list(paired_read_gen(mock_bam))

        assert len(pairs) == 0

    def test_filters_secondary_reads(self):
        """Test that secondary alignments are filtered."""
        from mapping.remap_utils import paired_read_gen

        mock_bam = MagicMock()

        read1 = MagicMock()
        read1.query_name = "read1"
        read1.is_proper_pair = True
        read1.is_secondary = True  # Secondary
        read1.is_supplementary = False

        mock_bam.fetch.return_value = iter([read1])

        pairs = list(paired_read_gen(mock_bam))

        assert len(pairs) == 0

    def test_filters_supplementary_reads(self):
        """Test that supplementary alignments are filtered."""
        from mapping.remap_utils import paired_read_gen

        mock_bam = MagicMock()

        read1 = MagicMock()
        read1.query_name = "read1"
        read1.is_proper_pair = True
        read1.is_secondary = False
        read1.is_supplementary = True  # Supplementary

        mock_bam.fetch.return_value = iter([read1])

        pairs = list(paired_read_gen(mock_bam))

        assert len(pairs) == 0

    def test_handles_read2_first(self):
        """Test correct pairing when read2 comes before read1."""
        from mapping.remap_utils import paired_read_gen

        mock_bam = MagicMock()

        # Read2 comes first in file
        read2 = MagicMock()
        read2.query_name = "read1"
        read2.is_proper_pair = True
        read2.is_secondary = False
        read2.is_supplementary = False
        read2.is_read1 = False

        read1 = MagicMock()
        read1.query_name = "read1"
        read1.is_proper_pair = True
        read1.is_secondary = False
        read1.is_supplementary = False
        read1.is_read1 = True

        mock_bam.fetch.return_value = iter([read2, read1])

        pairs = list(paired_read_gen(mock_bam))

        assert len(pairs) == 1
        # Should always be (read1, read2) order
        assert pairs[0] == (read1, read2)


class TestAlignPosGen:
    """Tests for alignment position generator."""

    def test_generates_positions(self):
        """Test basic position generation."""
        from mapping.remap_utils import align_pos_gen

        mock_read = MagicMock()
        mock_read.query_sequence = "ACGTACGTACGT"  # 12 bases

        align_dict = {100: 0, 101: 1, 102: 2, 103: 3, 104: 4, 105: 5,
                      106: 6, 107: 7, 108: 8, 109: 9, 110: 10, 111: 11}

        pos_list = [(102, 103), (106, 107)]  # Two SNP positions

        positions = list(align_pos_gen(mock_read, align_dict, pos_list))

        # Should yield: 0, start1, stop1, start2, stop2, len(seq)
        assert positions[0] == 0
        assert positions[-1] == 12  # len(query_sequence)

    def test_empty_pos_list(self):
        """Test with empty position list."""
        from mapping.remap_utils import align_pos_gen

        mock_read = MagicMock()
        mock_read.query_sequence = "ACGT"

        align_dict = {}
        pos_list = []

        positions = list(align_pos_gen(mock_read, align_dict, pos_list))

        assert positions == [0, 4]  # Just start and end


class TestBuildRef2ReadMaps:
    """Tests for reference-to-read position mapping (indel support)."""

    def test_simple_mapping(self):
        """Test mapping with no gaps."""
        from mapping.remap_utils import _build_ref2read_maps

        mock_read = MagicMock()
        # Simple alignment: all positions mapped
        mock_read.get_aligned_pairs.return_value = [
            (0, 100), (1, 101), (2, 102), (3, 103)
        ]

        left_map, right_map = _build_ref2read_maps(mock_read)

        # All positions should map directly
        assert left_map == {100: 0, 101: 1, 102: 2, 103: 3}
        assert right_map == {100: 0, 101: 1, 102: 2, 103: 3}

    def test_deletion_mapping(self):
        """Test mapping with deletion (ref pos without query pos)."""
        from mapping.remap_utils import _build_ref2read_maps

        mock_read = MagicMock()
        # Deletion at ref 101-102: query skips from 0 to 1
        mock_read.get_aligned_pairs.return_value = [
            (0, 100),
            (None, 101),  # Deleted base
            (None, 102),  # Deleted base
            (1, 103)
        ]

        left_map, right_map = _build_ref2read_maps(mock_read)

        # Left map uses last known query pos
        assert left_map[100] == 0
        assert left_map[101] == 0  # Uses last known (0)
        assert left_map[102] == 0  # Uses last known (0)
        assert left_map[103] == 1

        # Right map uses next known query pos
        assert right_map[100] == 0
        assert right_map[101] == 1  # Uses next known (1)
        assert right_map[102] == 1  # Uses next known (1)
        assert right_map[103] == 1

    def test_insertion_not_in_map(self):
        """Test that insertions (query pos without ref pos) don't appear in ref map."""
        from mapping.remap_utils import _build_ref2read_maps

        mock_read = MagicMock()
        # Insertion at query 1: ref skips from 100 to 101
        mock_read.get_aligned_pairs.return_value = [
            (0, 100),
            (1, None),  # Inserted base
            (2, 101)
        ]

        left_map, right_map = _build_ref2read_maps(mock_read)

        # Only ref positions appear in maps
        assert 100 in left_map
        assert 101 in left_map
        assert len(left_map) == 2


class TestFillInsertionQuals:
    """Tests for insertion quality score generation."""

    def test_uses_constant_without_flanks(self):
        """Test that constant quality is used without flanking data."""
        from mapping.remap_utils import _fill_insertion_quals

        result = _fill_insertion_quals(
            insert_len=3,
            left_qual=np.array([], dtype=np.uint8),
            right_qual=np.array([], dtype=np.uint8),
            insert_qual=25
        )

        assert len(result) == 3
        assert all(q == 25 for q in result)

    def test_averages_flanking_qualities(self):
        """Test that flanking qualities are averaged."""
        from mapping.remap_utils import _fill_insertion_quals

        result = _fill_insertion_quals(
            insert_len=2,
            left_qual=np.array([30, 30], dtype=np.uint8),
            right_qual=np.array([40, 40], dtype=np.uint8),
            insert_qual=20
        )

        assert len(result) == 2
        # Average of [30, 30, 40, 40] = 35
        assert all(q == 35 for q in result)

    def test_handles_one_sided_flanks(self):
        """Test with only one side of flanking data."""
        from mapping.remap_utils import _fill_insertion_quals

        result = _fill_insertion_quals(
            insert_len=2,
            left_qual=np.array([20, 20], dtype=np.uint8),
            right_qual=np.array([], dtype=np.uint8),
            insert_qual=10
        )

        assert len(result) == 2
        assert all(q == 20 for q in result)


class TestMakePhasedSeqs:
    """Tests for phased sequence creation (SNP-only version)."""

    def test_swaps_alleles(self):
        """Test basic allele swapping."""
        from mapping.remap_utils import make_phased_seqs

        # Split sequence: "ACG" + variant + "TT" + variant + "GA"
        split_seq = ["ACG", "C", "TT", "A", "GA"]

        hap1_alleles = ["T", "G"]  # Swap C->T, A->G
        hap2_alleles = ["C", "A"]  # Original

        hap1_seq, hap2_seq = make_phased_seqs(split_seq, hap1_alleles, hap2_alleles)

        assert hap1_seq == "ACGTTTAGGA"  # With swapped alleles
        assert hap2_seq == "ACGCTTAGA"   # Original alleles

    def test_preserves_flanking_sequence(self):
        """Test that non-variant sequence is preserved."""
        from mapping.remap_utils import make_phased_seqs

        split_seq = ["AAAA", "X", "TTTT"]
        hap1_alleles = ["G"]
        hap2_alleles = ["C"]

        hap1_seq, hap2_seq = make_phased_seqs(split_seq, hap1_alleles, hap2_alleles)

        assert hap1_seq.startswith("AAAA")
        assert hap1_seq.endswith("TTTT")
        assert hap2_seq.startswith("AAAA")
        assert hap2_seq.endswith("TTTT")


class TestMakePhasedSeqsWithQual:
    """Tests for phased sequence creation with quality scores (indel-aware)."""

    def test_handles_same_length_alleles(self):
        """Test SNPs where allele lengths match."""
        from mapping.remap_utils import make_phased_seqs_with_qual

        split_seq = ["AA", "C", "GG"]
        split_qual = [
            np.array([30, 30], dtype=np.uint8),
            np.array([35], dtype=np.uint8),
            np.array([40, 40], dtype=np.uint8)
        ]

        hap1_alleles = ["T"]  # Same length
        hap2_alleles = ["C"]

        (hap1_seq, hap1_qual), (hap2_seq, hap2_qual) = make_phased_seqs_with_qual(
            split_seq, split_qual, hap1_alleles, hap2_alleles
        )

        assert hap1_seq == "AATGG"
        assert hap2_seq == "AACGG"
        assert len(hap1_qual) == 5
        assert len(hap2_qual) == 5

    def test_handles_deletion(self):
        """Test deletion where alt allele is shorter."""
        from mapping.remap_utils import make_phased_seqs_with_qual

        split_seq = ["AA", "CCC", "GG"]  # Variant is 3 bases
        split_qual = [
            np.array([30, 30], dtype=np.uint8),
            np.array([35, 35, 35], dtype=np.uint8),
            np.array([40, 40], dtype=np.uint8)
        ]

        hap1_alleles = ["T"]  # Deletion: 3 -> 1 base
        hap2_alleles = ["CCC"]

        (hap1_seq, hap1_qual), (hap2_seq, hap2_qual) = make_phased_seqs_with_qual(
            split_seq, split_qual, hap1_alleles, hap2_alleles
        )

        assert hap1_seq == "AATGG"  # Shorter by 2
        assert hap2_seq == "AACCCGG"
        assert len(hap1_qual) == 5
        assert len(hap2_qual) == 7

    def test_handles_insertion(self):
        """Test insertion where alt allele is longer."""
        from mapping.remap_utils import make_phased_seqs_with_qual

        split_seq = ["AA", "C", "GG"]  # Variant is 1 base
        split_qual = [
            np.array([30, 30], dtype=np.uint8),
            np.array([35], dtype=np.uint8),
            np.array([40, 40], dtype=np.uint8)
        ]

        hap1_alleles = ["TTTT"]  # Insertion: 1 -> 4 bases
        hap2_alleles = ["C"]

        (hap1_seq, hap1_qual), (hap2_seq, hap2_qual) = make_phased_seqs_with_qual(
            split_seq, split_qual, hap1_alleles, hap2_alleles, insert_qual=25
        )

        assert hap1_seq == "AATTTTGG"  # Longer by 3
        assert hap2_seq == "AACGG"
        assert len(hap1_qual) == 8  # Extra qualities generated
        assert len(hap2_qual) == 5


class TestMakeMultiSeqs:
    """Tests for multi-sample sequence creation."""

    def test_creates_multiple_sequences(self):
        """Test creating sequences for multiple haplotypes."""
        from mapping.remap_utils import make_multi_seqs

        split_seq = ["AA", "X", "TT"]

        allele_combos = [
            ["G"],  # Haplotype 1
            ["C"],  # Haplotype 2
            ["A"],  # Haplotype 3
        ]

        seqs = make_multi_seqs(split_seq, allele_combos)

        assert len(seqs) == 3
        assert seqs[0] == "AAGTT"
        assert seqs[1] == "AACTT"
        assert seqs[2] == "AAATT"

    def test_handles_multiple_variants(self):
        """Test with multiple variants per read."""
        from mapping.remap_utils import make_multi_seqs

        split_seq = ["A", "X", "T", "Y", "G"]

        allele_combos = [
            ["C", "A"],  # Two variants
            ["G", "T"],
        ]

        seqs = make_multi_seqs(split_seq, allele_combos)

        assert seqs[0] == "ACTAG"
        assert seqs[1] == "AGTTG"


class TestWriteRead:
    """Tests for BAM write functionality."""

    def test_writes_snp_mode(self):
        """Test writing read in SNP mode (no quality array)."""
        from mapping.remap_utils import write_read

        mock_bam = MagicMock()
        mock_read = MagicMock()
        mock_read.query_qualities = np.array([30, 30, 30], dtype=np.uint8)

        write_read(mock_bam, mock_read, "ACG", "new_name", new_qual=None)

        mock_bam.write.assert_called_once_with(mock_read)
        assert mock_read.query_sequence == "ACG"
        assert mock_read.query_name == "new_name"

    def test_writes_indel_mode(self):
        """Test writing read in indel mode (with quality array)."""
        from mapping.remap_utils import write_read

        mock_bam = MagicMock()
        mock_read = MagicMock()
        mock_read.query_length = 3  # Original length
        mock_read.cigartuples = [(0, 3)]

        new_qual = np.array([25, 25, 25, 25], dtype=np.uint8)

        write_read(mock_bam, mock_read, "ACGT", "new_name", new_qual=new_qual)

        mock_bam.write.assert_called_once_with(mock_read)
        assert mock_read.query_sequence == "ACGT"
        # CIGAR should be updated for length change
        assert mock_read.cigartuples == [(0, 4)]
