"""
Unit tests for the counting.filter_variant_data module.

Tests cover:
- GTF to BED conversion
- BED intersection file parsing
- VCF to BED conversion (via wasp2.io)
"""

import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import polars as pl
import pytest


class TestGtfToBed:
    """Tests for GTF to BED conversion."""

    def test_converts_gtf_to_bed(self, tmp_path):
        """Test basic GTF to BED conversion."""
        from counting.filter_variant_data import gtf_to_bed

        # Create minimal GTF content
        gtf_content = """\
chr1\tensembl\texon\t100\t200\t.\t+\t.\tgene_id "GENE1";
chr1\tensembl\tgene\t100\t500\t.\t+\t.\tgene_id "GENE1";
chr1\tensembl\texon\t300\t400\t.\t+\t.\tgene_id "GENE1";
chr2\tensembl\texon\t1000\t2000\t.\t-\t.\tgene_id "GENE2";
"""
        gtf_file = tmp_path / "test.gtf"
        gtf_file.write_text(gtf_content)

        out_bed = tmp_path / "output.bed"

        result = gtf_to_bed(gtf_file, out_bed, feature="exon", attribute="gene_id")

        assert Path(result).exists()

        # Read and verify output
        df = pl.read_csv(out_bed, separator="\t", has_header=False)
        assert len(df) == 3  # 3 exon features

        # Check coordinates are 0-based (start - 1)
        assert df[0, 1] == 99  # First exon start

    def test_filters_by_feature(self, tmp_path):
        """Test that GTF conversion filters by feature type."""
        from counting.filter_variant_data import gtf_to_bed

        gtf_content = """\
chr1\tensembl\texon\t100\t200\t.\t+\t.\tgene_id "GENE1";
chr1\tensembl\tCDS\t150\t190\t.\t+\t.\tgene_id "GENE1";
chr1\tensembl\texon\t300\t400\t.\t+\t.\tgene_id "GENE1";
"""
        gtf_file = tmp_path / "test.gtf"
        gtf_file.write_text(gtf_content)

        out_bed = tmp_path / "output.bed"

        gtf_to_bed(gtf_file, out_bed, feature="CDS", attribute="gene_id")

        df = pl.read_csv(out_bed, separator="\t", has_header=False)
        assert len(df) == 1  # Only 1 CDS feature

    def test_extracts_different_attributes(self, tmp_path):
        """Test extraction of different attribute types."""
        from counting.filter_variant_data import gtf_to_bed

        gtf_content = """\
chr1\tensembl\texon\t100\t200\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1";
chr1\tensembl\texon\t300\t400\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX2";
"""
        gtf_file = tmp_path / "test.gtf"
        gtf_file.write_text(gtf_content)

        out_bed = tmp_path / "output.bed"

        gtf_to_bed(gtf_file, out_bed, feature="exon", attribute="transcript_id")

        df = pl.read_csv(out_bed, separator="\t", has_header=False)
        assert len(df) == 2

        # Fourth column should have transcript IDs
        assert df[0, 3] == "TX1"
        assert df[1, 3] == "TX2"


class TestParseIntersectRegion:
    """Tests for BED intersection file parsing."""

    def test_parses_basic_intersect(self, tmp_path):
        """Test parsing of basic intersection file."""
        from counting.filter_variant_data import parse_intersect_region

        # Minimal intersection output: chrom, pos0, pos, ref, alt
        intersect_content = """\
chr1\t99\t100\tA\tG
chr1\t199\t200\tC\tT
chr2\t999\t1000\tG\tA
"""
        intersect_file = tmp_path / "intersect.bed"
        intersect_file.write_text(intersect_content)

        df = parse_intersect_region(str(intersect_file))

        assert len(df) == 3
        assert list(df.columns) == ["chrom", "pos", "ref", "alt"]
        assert df["pos"].to_list() == [100, 200, 1000]

    def test_parses_intersect_with_region_coords(self, tmp_path):
        """Test parsing intersection with region coordinates."""
        from counting.filter_variant_data import parse_intersect_region

        # With region coords: chrom, pos0, pos, ref, alt, region_chrom, region_start, region_end
        intersect_content = """\
chr1\t99\t100\tA\tG\tchr1\t50\t150
chr1\t199\t200\tC\tT\tchr1\t150\t250
"""
        intersect_file = tmp_path / "intersect.bed"
        intersect_file.write_text(intersect_content)

        df = parse_intersect_region(str(intersect_file))

        assert len(df) == 2
        assert "region" in df.columns
        # Region should be formatted as chrom_start_end
        assert df["region"][0] == "chr1_50_150"

    def test_parses_intersect_with_region_names(self, tmp_path):
        """Test parsing intersection with named regions."""
        from counting.filter_variant_data import parse_intersect_region

        # With region name in col 9: chrom, pos0, pos, ref, alt, region_chrom, region_start, region_end, name
        intersect_content = """\
chr1\t99\t100\tA\tG\tchr1\t50\t150\tPEAK1
chr1\t199\t200\tC\tT\tchr1\t150\t250\tPEAK2
"""
        intersect_file = tmp_path / "intersect.bed"
        intersect_file.write_text(intersect_content)

        df = parse_intersect_region(str(intersect_file), use_region_names=True)

        assert len(df) == 2
        assert "region" in df.columns
        assert df["region"][0] == "PEAK1"
        assert df["region"][1] == "PEAK2"

    def test_removes_duplicates(self, tmp_path):
        """Test that duplicate entries are removed."""
        from counting.filter_variant_data import parse_intersect_region

        intersect_content = """\
chr1\t99\t100\tA\tG
chr1\t99\t100\tA\tG
chr1\t199\t200\tC\tT
"""
        intersect_file = tmp_path / "intersect.bed"
        intersect_file.write_text(intersect_content)

        df = parse_intersect_region(str(intersect_file))

        assert len(df) == 2  # Duplicates removed

    def test_custom_region_column_name(self, tmp_path):
        """Test using custom region column name."""
        from counting.filter_variant_data import parse_intersect_region

        intersect_content = """\
chr1\t99\t100\tA\tG\tchr1\t50\t150
"""
        intersect_file = tmp_path / "intersect.bed"
        intersect_file.write_text(intersect_content)

        df = parse_intersect_region(str(intersect_file), region_col="peak")

        assert "peak" in df.columns


class TestParseIntersectRegionNew:
    """Tests for the newer parse_intersect_region_new function."""

    def test_parses_with_samples(self, tmp_path):
        """Test parsing intersection with sample genotype columns."""
        from counting.filter_variant_data import parse_intersect_region_new

        # With sample GT columns: chrom, pos0, pos, ref, alt, sample1_gt, sample2_gt
        intersect_content = """\
chr1\t99\t100\tA\tG\t0/1\t0/0
chr1\t199\t200\tC\tT\t1/1\t0/1
"""
        intersect_file = tmp_path / "intersect.bed"
        intersect_file.write_text(intersect_content)

        df = parse_intersect_region_new(
            str(intersect_file),
            samples=["sample1", "sample2"]
        )

        assert len(df) == 2
        assert "sample1" in df.columns
        assert "sample2" in df.columns


class TestVcfToBed:
    """Tests for VCF to BED conversion wrapper."""

    @patch("counting.filter_variant_data._variants_to_bed")
    def test_calls_variants_to_bed(self, mock_variants_to_bed, tmp_path):
        """Test that vcf_to_bed correctly calls the underlying function."""
        from counting.filter_variant_data import vcf_to_bed

        mock_variants_to_bed.return_value = tmp_path / "output.bed"

        vcf_file = tmp_path / "test.vcf"
        out_bed = tmp_path / "output.bed"

        result = vcf_to_bed(vcf_file, out_bed)

        mock_variants_to_bed.assert_called_once()
        call_kwargs = mock_variants_to_bed.call_args[1]
        assert call_kwargs["variant_file"] == vcf_file
        assert call_kwargs["out_bed"] == out_bed
        assert call_kwargs["het_only"] is False  # No samples provided

    @patch("counting.filter_variant_data._variants_to_bed")
    def test_enables_het_only_with_samples(self, mock_variants_to_bed, tmp_path):
        """Test that het_only is enabled when samples provided."""
        from counting.filter_variant_data import vcf_to_bed

        mock_variants_to_bed.return_value = tmp_path / "output.bed"

        vcf_file = tmp_path / "test.vcf"
        out_bed = tmp_path / "output.bed"

        vcf_to_bed(vcf_file, out_bed, samples=["sample1"])

        call_kwargs = mock_variants_to_bed.call_args[1]
        assert call_kwargs["het_only"] is True
        assert call_kwargs["samples"] == ["sample1"]

    @patch("counting.filter_variant_data._variants_to_bed")
    def test_passes_indel_parameters(self, mock_variants_to_bed, tmp_path):
        """Test that indel parameters are passed correctly."""
        from counting.filter_variant_data import vcf_to_bed

        mock_variants_to_bed.return_value = tmp_path / "output.bed"

        vcf_file = tmp_path / "test.vcf"
        out_bed = tmp_path / "output.bed"

        vcf_to_bed(
            vcf_file, out_bed,
            include_indels=True,
            max_indel_len=20
        )

        call_kwargs = mock_variants_to_bed.call_args[1]
        assert call_kwargs["include_indels"] is True
        assert call_kwargs["max_indel_len"] == 20
