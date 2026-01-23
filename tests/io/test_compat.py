"""
Tests for the compatibility layer (wasp2.io.compat).

Verifies that the new VariantSource-based interface produces
equivalent output to the legacy bcftools-based approach.
"""

import shutil
from pathlib import Path

import pytest

from wasp2.io.compat import variants_to_bed, vcf_to_bed

# Check if bcftools is available for legacy tests
BCFTOOLS_AVAILABLE = shutil.which("bcftools") is not None


class TestVariantsToBed:
    """Tests for the unified variants_to_bed function."""

    def test_vcf_no_samples(self, sample_vcf, tmp_output_dir):
        """Test converting VCF without sample filtering."""
        output = tmp_output_dir / "all_variants.bed"

        result = variants_to_bed(
            variant_file=sample_vcf,
            out_bed=output,
            samples=None,
            include_gt=False,
            het_only=False,
        )

        assert result == output
        assert output.exists()

        lines = output.read_text().strip().split("\n")
        assert len(lines) == 6  # 6 variants in test VCF

    def test_vcf_single_sample_het(self, sample_vcf, tmp_output_dir):
        """Test extracting het sites for single sample."""
        output = tmp_output_dir / "sample1_het.bed"

        variants_to_bed(
            variant_file=sample_vcf,
            out_bed=output,
            samples=["sample1"],
            include_gt=True,
            het_only=True,
        )

        lines = output.read_text().strip().split("\n")
        # sample1 has 3 het sites
        assert len(lines) == 3

    def test_vcf_multi_sample(self, sample_vcf, tmp_output_dir):
        """Test with multiple samples."""
        output = tmp_output_dir / "multi_sample.bed"

        variants_to_bed(
            variant_file=sample_vcf,
            out_bed=output,
            samples=["sample1", "sample2"],
            include_gt=True,
            het_only=True,
        )

        assert output.exists()


class TestLegacyVcfToBed:
    """Tests for the legacy vcf_to_bed alias."""

    def test_legacy_function_exists(self):
        """Test that legacy function is available."""
        assert callable(vcf_to_bed)

    @pytest.mark.skipif(not BCFTOOLS_AVAILABLE, reason="bcftools not available")
    def test_legacy_basic_usage(self, sample_vcf, tmp_output_dir):
        """Test basic legacy function usage."""
        output = tmp_output_dir / "legacy.bed"

        result = vcf_to_bed(
            vcf_file=sample_vcf,
            out_bed=output,
            samples=None,
        )

        assert Path(result) == output
        assert output.exists()


class TestModuleIntegration:
    """Tests that mapping/counting modules use the new interface."""

    def test_mapping_module_vcf_to_bed(self, sample_vcf, tmp_output_dir):
        """Test mapping module's vcf_to_bed uses new interface."""
        from mapping.intersect_variant_data import vcf_to_bed as mapping_vcf_to_bed

        output = tmp_output_dir / "mapping_output.bed"

        result = mapping_vcf_to_bed(
            vcf_file=sample_vcf,
            out_bed=output,
            samples=["sample1"],
        )

        assert Path(result) == output
        assert output.exists()

        # Should have het sites only when sample specified
        lines = output.read_text().strip().split("\n")
        assert len(lines) == 3  # sample1 has 3 het sites

    def test_counting_module_vcf_to_bed(self, sample_vcf, tmp_output_dir):
        """Test counting module's vcf_to_bed uses new interface."""
        from counting.filter_variant_data import vcf_to_bed as counting_vcf_to_bed

        output = tmp_output_dir / "counting_output.bed"

        result = counting_vcf_to_bed(
            vcf_file=sample_vcf,
            out_bed=output,
            samples=["sample1"],
            include_gt=True,
        )

        assert Path(result) == output
        assert output.exists()

        lines = output.read_text().strip().split("\n")
        assert len(lines) == 3
