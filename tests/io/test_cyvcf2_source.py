"""
Tests for CyVCF2Source implementation.

These tests verify the high-performance cyvcf2-based VCF reader.
Tests are skipped if cyvcf2 is not installed.

Run with: pytest tests/io/test_cyvcf2_source.py -v
"""

import pytest
from pathlib import Path

from wasp2.io.variant_source import VariantSource, Variant, Genotype, VariantGenotype

# Check if cyvcf2 is available
try:
    import cyvcf2
    from wasp2.io.cyvcf2_source import CyVCF2Source
    CYVCF2_AVAILABLE = True
except ImportError:
    CYVCF2_AVAILABLE = False

pytestmark = pytest.mark.skipif(
    not CYVCF2_AVAILABLE,
    reason="cyvcf2 not installed - install with: pip install wasp2[cyvcf2]"
)


@pytest.mark.skipif(not CYVCF2_AVAILABLE, reason="cyvcf2 not available")
class TestCyVCF2SourceBasics:
    """Basic CyVCF2Source tests."""

    def test_direct_instantiation(self, sample_vcf_gz):
        """Test direct instantiation of CyVCF2Source."""
        source = CyVCF2Source(sample_vcf_gz)
        assert source is not None
        assert len(source.samples) == 2
        source.close()

    def test_open_vcf_gz_file(self, sample_vcf_gz):
        """Test opening a compressed VCF file with CyVCF2Source."""
        # Note: Need to use special extension to force cyvcf2 usage
        # or test direct instantiation
        source = CyVCF2Source(sample_vcf_gz)
        try:
            assert source.samples == ["sample1", "sample2"]
        finally:
            source.close()

    def test_samples_property(self, sample_vcf_gz):
        """Test getting sample list."""
        with CyVCF2Source(sample_vcf_gz) as source:
            samples = source.samples
            assert samples == ["sample1", "sample2"]

    def test_sample_count(self, sample_vcf_gz):
        """Test sample count."""
        with CyVCF2Source(sample_vcf_gz) as source:
            assert source.sample_count == 2

    def test_variant_count(self, sample_vcf_gz):
        """Test variant count."""
        with CyVCF2Source(sample_vcf_gz) as source:
            assert source.variant_count == 6


@pytest.mark.skipif(not CYVCF2_AVAILABLE, reason="cyvcf2 not available")
class TestCyVCF2SourceIteration:
    """Tests for iterating over VCF variants with cyvcf2."""

    def test_iter_all_variants(self, sample_vcf_gz):
        """Test iterating over all variants."""
        with CyVCF2Source(sample_vcf_gz) as source:
            variants = list(source.iter_variants())

            assert len(variants) == 6

            # Check first variant
            first = variants[0]
            assert first.variant.chrom == "chr1"
            assert first.variant.pos == 100
            assert first.variant.ref == "A"
            assert first.variant.alt == "G"
            assert first.variant.id == "rs1"

    def test_iter_variants_het_only(self, sample_vcf_gz):
        """Test iterating over het sites for sample1."""
        with CyVCF2Source(sample_vcf_gz) as source:
            het_sites = list(source.iter_variants(samples=["sample1"], het_only=True))

            # sample1 has 3 het sites: rs1, rs4, rs5
            assert len(het_sites) == 3

            for vg in het_sites:
                assert vg.genotype == Genotype.HET

            # Verify it's the right variants
            ids = [vg.variant.id for vg in het_sites]
            assert "rs1" in ids
            assert "rs4" in ids
            assert "rs5" in ids

    def test_iter_variants_single_sample(self, sample_vcf_gz):
        """Test iterating for a specific sample."""
        with CyVCF2Source(sample_vcf_gz) as source:
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

    def test_allele_extraction(self, sample_vcf_gz):
        """Test that alleles are correctly extracted."""
        with CyVCF2Source(sample_vcf_gz) as source:
            variants = list(source.iter_variants(samples=["sample1"]))

            # rs1: 0/1 for sample1 (A/G)
            first = variants[0]
            assert first.allele1 == "A"
            assert first.allele2 == "G"

            # rs2: 1/1 for sample1 (T/T)
            second = variants[1]
            assert second.allele1 == "T"
            assert second.allele2 == "T"

    def test_missing_genotype(self, sample_vcf_gz):
        """Test handling of missing genotypes."""
        with CyVCF2Source(sample_vcf_gz) as source:
            # rs6 has missing genotype (./.) for sample1
            variants = list(source.iter_variants(samples=["sample1"]))
            rs6 = variants[5]  # Last variant

            assert rs6.variant.id == "rs6"
            assert rs6.genotype == Genotype.MISSING


@pytest.mark.skipif(not CYVCF2_AVAILABLE, reason="cyvcf2 not available")
class TestCyVCF2SourceQueries:
    """Tests for querying specific positions and regions."""

    def test_get_genotype(self, sample_vcf_gz):
        """Test getting genotype at a specific position."""
        with CyVCF2Source(sample_vcf_gz) as source:
            # rs1 at chr1:100 is 0/1 for sample1
            gt = source.get_genotype("sample1", "chr1", 100)
            assert gt == Genotype.HET

            # rs2 at chr1:200 is 1/1 for sample1
            gt = source.get_genotype("sample1", "chr1", 200)
            assert gt == Genotype.HOM_ALT

            # rs3 at chr1:300 is 0/0 for sample1
            gt = source.get_genotype("sample1", "chr1", 300)
            assert gt == Genotype.HOM_REF

    def test_query_region(self, sample_vcf_gz):
        """Test querying a genomic region."""
        with CyVCF2Source(sample_vcf_gz) as source:
            # Query chr1:100-300 (should get rs1, rs2, rs3)
            variants = list(source.query_region("chr1", 100, 300, samples=["sample1"]))

            assert len(variants) == 3
            ids = [v.variant.id for v in variants]
            assert ids == ["rs1", "rs2", "rs3"]

    def test_query_region_single_variant(self, sample_vcf_gz):
        """Test querying a region with a single variant."""
        with CyVCF2Source(sample_vcf_gz) as source:
            # Query chr1:100-100 (should get only rs1)
            variants = list(source.query_region("chr1", 100, 100, samples=["sample1"]))

            assert len(variants) == 1
            assert variants[0].variant.id == "rs1"

    def test_query_region_chromosome(self, sample_vcf_gz):
        """Test querying different chromosomes."""
        with CyVCF2Source(sample_vcf_gz) as source:
            # chr2 has 2 variants: rs5, rs6
            variants = list(source.query_region("chr2", 1, 1000, samples=["sample1"]))

            assert len(variants) == 2
            ids = [v.variant.id for v in variants]
            assert "rs5" in ids
            assert "rs6" in ids


@pytest.mark.skipif(not CYVCF2_AVAILABLE, reason="cyvcf2 not available")
class TestCyVCF2SourceBED:
    """Tests for BED export functionality."""

    def test_to_bed_basic(self, sample_vcf_gz, tmp_path):
        """Test basic BED export."""
        with CyVCF2Source(sample_vcf_gz) as source:
            bed_path = tmp_path / "test.bed"
            result = source.to_bed(
                bed_path,
                samples=["sample1"],
                het_only=False,
                include_genotypes=False
            )

            assert result.exists()
            assert result == bed_path

            # Read and check content
            lines = bed_path.read_text().strip().split("\n")
            assert len(lines) > 0

    def test_to_bed_het_only(self, sample_vcf_gz, tmp_path):
        """Test BED export with het_only filter."""
        with CyVCF2Source(sample_vcf_gz) as source:
            bed_path = tmp_path / "test_het.bed"
            source.to_bed(
                bed_path,
                samples=["sample1"],
                het_only=True,
                include_genotypes=True
            )

            assert bed_path.exists()

            # Should have het sites for sample1: rs1, rs4, rs5
            lines = bed_path.read_text().strip().split("\n")
            # Note: bcftools filters, so exact count depends on filtering
            assert len(lines) > 0


@pytest.mark.skipif(not CYVCF2_AVAILABLE, reason="cyvcf2 not available")
class TestCyVCF2SourceComparison:
    """Tests comparing CyVCF2Source with VCFSource for correctness."""

    def test_same_variants_as_vcfsource(self, sample_vcf_gz):
        """Verify CyVCF2Source returns same variants as VCFSource."""
        from wasp2.io.vcf_source import VCFSource

        # Get variants from pysam VCFSource
        with VCFSource(sample_vcf_gz) as pysam_source:
            pysam_variants = list(pysam_source.iter_variants())

        # Get variants from cyvcf2 CyVCF2Source
        with CyVCF2Source(sample_vcf_gz) as cyvcf2_source:
            cyvcf2_variants = list(cyvcf2_source.iter_variants())

        # Should have same number of variants
        assert len(pysam_variants) == len(cyvcf2_variants)

        # Check each variant matches
        for pv, cv in zip(pysam_variants, cyvcf2_variants):
            assert pv.variant.chrom == cv.variant.chrom
            assert pv.variant.pos == cv.variant.pos
            assert pv.variant.ref == cv.variant.ref
            assert pv.variant.alt == cv.variant.alt
            assert pv.variant.id == cv.variant.id
            assert pv.genotype == cv.genotype

    def test_same_het_sites_as_vcfsource(self, sample_vcf_gz):
        """Verify CyVCF2Source returns same het sites as VCFSource."""
        from wasp2.io.vcf_source import VCFSource

        # Get het sites from pysam VCFSource
        with VCFSource(sample_vcf_gz) as pysam_source:
            pysam_hets = list(pysam_source.iter_variants(samples=["sample1"], het_only=True))

        # Get het sites from cyvcf2 CyVCF2Source
        with CyVCF2Source(sample_vcf_gz) as cyvcf2_source:
            cyvcf2_hets = list(cyvcf2_source.iter_variants(samples=["sample1"], het_only=True))

        # Should have same het sites
        assert len(pysam_hets) == len(cyvcf2_hets)

        # Check positions match
        pysam_positions = [(v.variant.chrom, v.variant.pos) for v in pysam_hets]
        cyvcf2_positions = [(v.variant.chrom, v.variant.pos) for v in cyvcf2_hets]
        assert pysam_positions == cyvcf2_positions


@pytest.mark.skipif(not CYVCF2_AVAILABLE, reason="cyvcf2 not available")
class TestCyVCF2SourceErrors:
    """Tests for error handling."""

    def test_invalid_sample(self, sample_vcf_gz):
        """Test error when requesting invalid sample."""
        with CyVCF2Source(sample_vcf_gz) as source:
            with pytest.raises(ValueError, match="not found"):
                list(source.iter_variants(samples=["nonexistent"]))

    def test_nonexistent_file(self):
        """Test error when file doesn't exist."""
        with pytest.raises(ValueError):
            CyVCF2Source("/nonexistent/file.vcf.gz")

    def test_invalid_position(self, sample_vcf_gz):
        """Test error when querying invalid position."""
        with CyVCF2Source(sample_vcf_gz) as source:
            with pytest.raises(ValueError):
                source.get_genotype("sample1", "chrNONE", 999999)
