"""
Tests for VCFSource implementation.

These tests focus on VCF-specific functionality and don't require plink2.
Run with: pytest tests/io/test_vcf_source.py -v
"""

import pytest
from pathlib import Path

from wasp2.io.variant_source import VariantSource, Variant, Genotype, VariantGenotype
from wasp2.io.vcf_source import VCFSource


class TestVCFSourceBasics:
    """Basic VCFSource tests."""

    def test_open_vcf_file(self, sample_vcf):
        """Test opening a VCF file."""
        with VariantSource.open(sample_vcf) as source:
            assert isinstance(source, VCFSource)
            assert source.validate() is True

    def test_open_vcf_gz_file(self, sample_vcf_gz):
        """Test opening a compressed VCF file."""
        with VariantSource.open(sample_vcf_gz) as source:
            assert isinstance(source, VCFSource)
            assert source.validate() is True

    def test_samples_property(self, sample_vcf):
        """Test getting sample list."""
        with VariantSource.open(sample_vcf) as source:
            samples = source.samples
            assert samples == ["sample1", "sample2"]

    def test_sample_count(self, sample_vcf):
        """Test sample count."""
        with VariantSource.open(sample_vcf) as source:
            assert source.sample_count == 2

    def test_variant_count(self, sample_vcf):
        """Test variant count."""
        with VariantSource.open(sample_vcf) as source:
            assert source.variant_count == 6


class TestVCFSourceIteration:
    """Tests for iterating over VCF variants."""

    def test_iter_all_variants(self, sample_vcf, vcf_expected_variants):
        """Test iterating over all variants."""
        with VariantSource.open(sample_vcf) as source:
            variants = list(source.iter_variants())

            assert len(variants) == 6

            # Check first variant
            first = variants[0]
            assert first.variant.chrom == "chr1"
            assert first.variant.pos == 100
            assert first.variant.ref == "A"
            assert first.variant.alt == "G"
            assert first.variant.id == "rs1"

    def test_iter_variants_het_only(self, sample_vcf, vcf_expected_het_sites_sample1):
        """Test iterating over het sites for sample1."""
        with VariantSource.open(sample_vcf) as source:
            het_sites = list(source.iter_variants(samples=["sample1"], het_only=True))

            # sample1 has 3 het sites: rs1, rs4, rs5
            assert len(het_sites) == 3

            for vg in het_sites:
                assert vg.genotype == Genotype.HET

    def test_iter_variants_single_sample(self, sample_vcf):
        """Test iterating for a specific sample."""
        with VariantSource.open(sample_vcf) as source:
            variants = list(source.iter_variants(samples=["sample2"]))

            # Should get all 6 variants for sample2
            assert len(variants) == 6

            # Check genotypes for sample2 based on our test VCF:
            # rs1: 0/0 (HOM_REF), rs2: 0/1 (HET), rs3: 1/1 (HOM_ALT)
            # rs4: 0/1 (HET), rs5: 0/0 (HOM_REF), rs6: 0/1 (HET)
            genotypes = [v.genotype for v in variants]
            assert genotypes[0] == Genotype.HOM_REF  # rs1
            assert genotypes[1] == Genotype.HET      # rs2
            assert genotypes[2] == Genotype.HOM_ALT  # rs3
            assert genotypes[3] == Genotype.HET      # rs4
            assert genotypes[4] == Genotype.HOM_REF  # rs5
            assert genotypes[5] == Genotype.HET      # rs6

    def test_get_sample_idx(self, sample_vcf):
        """Test getting sample index."""
        with VariantSource.open(sample_vcf) as source:
            assert source.get_sample_idx("sample1") == 0
            assert source.get_sample_idx("sample2") == 1

    def test_get_sample_idx_invalid(self, sample_vcf):
        """Test invalid sample ID raises error."""
        with VariantSource.open(sample_vcf) as source:
            with pytest.raises(ValueError, match="not found"):
                source.get_sample_idx("nonexistent")


class TestVCFSourceToBed:
    """Tests for BED output functionality."""

    def test_to_bed_all_variants(self, sample_vcf, tmp_output_dir):
        """Test exporting all variants to BED."""
        output = tmp_output_dir / "all.bed"

        with VariantSource.open(sample_vcf) as source:
            result = source.to_bed(output, het_only=False, include_genotypes=False)

        assert result == output
        assert output.exists()

        lines = output.read_text().strip().split('\n')
        assert len(lines) == 6

        # Check format of first line
        fields = lines[0].split('\t')
        assert fields[0] == "chr1"
        assert fields[1] == "99"   # 0-based start
        assert fields[2] == "100"  # 1-based end
        assert fields[3] == "A"
        assert fields[4] == "G"

    def test_to_bed_het_only(self, sample_vcf, tmp_output_dir):
        """Test exporting het sites only."""
        output = tmp_output_dir / "het.bed"

        with VariantSource.open(sample_vcf) as source:
            source.to_bed(output, samples=["sample1"], het_only=True)

        lines = output.read_text().strip().split('\n')
        # sample1 has het at rs1, rs4, rs5
        assert len(lines) == 3

    def test_to_bed_with_genotypes(self, sample_vcf, tmp_output_dir):
        """Test BED with genotype columns."""
        output = tmp_output_dir / "with_gt.bed"

        with VariantSource.open(sample_vcf) as source:
            source.to_bed(
                output,
                samples=["sample1"],
                het_only=False,
                include_genotypes=True
            )

        lines = output.read_text().strip().split('\n')
        fields = lines[0].split('\t')

        # Should have at least 6 columns with genotype
        assert len(fields) >= 6


class TestVCFSourceQueryRegion:
    """Tests for region queries."""

    def test_query_region(self, sample_vcf_gz):
        """Test querying a region."""
        with VariantSource.open(sample_vcf_gz) as source:
            variants = list(source.query_region("chr1", 100, 300))

            positions = [v.variant.pos for v in variants]
            assert 100 in positions
            assert 200 in positions
            assert 300 in positions

    def test_query_region_empty(self, sample_vcf_gz):
        """Test querying empty region."""
        with VariantSource.open(sample_vcf_gz) as source:
            variants = list(source.query_region("chr1", 500, 600))
            assert len(variants) == 0

    def test_query_region_single_variant(self, sample_vcf_gz):
        """Test querying single position."""
        with VariantSource.open(sample_vcf_gz) as source:
            variants = list(source.query_region("chr1", 100, 100))
            assert len(variants) == 1
            assert variants[0].variant.pos == 100


class TestVCFSourceMissingData:
    """Tests for handling missing genotype data."""

    def test_missing_genotype(self, sample_vcf):
        """Test handling of missing genotype (./.)."""
        with VariantSource.open(sample_vcf) as source:
            # rs6 at chr2:200 has ./. for sample1
            variants = list(source.iter_variants(samples=["sample1"]))

            # Find rs6
            rs6 = next(v for v in variants if v.variant.id == "rs6")
            assert rs6.genotype == Genotype.MISSING

    def test_het_only_excludes_missing(self, sample_vcf):
        """Test that het_only filters out missing genotypes."""
        with VariantSource.open(sample_vcf) as source:
            het_sites = list(source.iter_variants(samples=["sample1"], het_only=True))

            # Should not include missing sites
            for vg in het_sites:
                assert vg.genotype != Genotype.MISSING
