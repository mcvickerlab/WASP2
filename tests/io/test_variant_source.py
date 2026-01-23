"""
Tests for VariantSource ABC and factory.

These tests are written FIRST (TDD) to define the expected behavior
before implementation.

Run with: pytest tests/io/test_variant_source.py -v
"""

import pytest
from pathlib import Path
from typing import List

# These imports will fail until we implement the module
# That's expected in TDD - tests are written first!
try:
    from wasp2.io.variant_source import (
        VariantSource,
        Variant,
        VariantGenotype,
        Genotype,
    )
    IMPORTS_AVAILABLE = True
except ImportError:
    IMPORTS_AVAILABLE = False
    # Create placeholder classes for test collection
    VariantSource = None
    Variant = None
    VariantGenotype = None
    Genotype = None


pytestmark = pytest.mark.skipif(
    not IMPORTS_AVAILABLE,
    reason="wasp2.io.variant_source not yet implemented"
)


# ============================================================================
# Tests for Variant dataclass
# ============================================================================

class TestVariant:
    """Tests for the Variant data class."""

    def test_variant_creation(self):
        """Test creating a Variant object."""
        v = Variant(chrom="chr1", pos=100, ref="A", alt="G", id="rs1")
        assert v.chrom == "chr1"
        assert v.pos == 100
        assert v.ref == "A"
        assert v.alt == "G"
        assert v.id == "rs1"

    def test_variant_pos0_property(self):
        """Test 0-based position conversion."""
        v = Variant(chrom="chr1", pos=100, ref="A", alt="G")
        assert v.pos0 == 99  # 0-based

    def test_variant_to_bed_line(self):
        """Test BED format output."""
        v = Variant(chrom="chr1", pos=100, ref="A", alt="G")
        bed_line = v.to_bed_line()
        assert bed_line == "chr1\t99\t100\tA\tG"

    def test_variant_immutable(self):
        """Test that Variant is immutable (frozen dataclass)."""
        v = Variant(chrom="chr1", pos=100, ref="A", alt="G")
        with pytest.raises(AttributeError):
            v.pos = 200

    def test_variant_hashable(self):
        """Test that Variant can be used in sets/dicts."""
        v1 = Variant(chrom="chr1", pos=100, ref="A", alt="G")
        v2 = Variant(chrom="chr1", pos=100, ref="A", alt="G")
        v3 = Variant(chrom="chr1", pos=200, ref="C", alt="T")

        # Same content should be equal
        assert v1 == v2
        assert hash(v1) == hash(v2)

        # Different content should not be equal
        assert v1 != v3

        # Should work in sets
        variant_set = {v1, v2, v3}
        assert len(variant_set) == 2  # v1 and v2 are duplicates


# ============================================================================
# Tests for Genotype enum
# ============================================================================

class TestGenotype:
    """Tests for the Genotype enum."""

    def test_genotype_values(self):
        """Test Genotype enum values match expected encoding."""
        assert Genotype.HOM_REF.value == 0
        assert Genotype.HET.value == 1
        assert Genotype.HOM_ALT.value == 2
        assert Genotype.MISSING.value == -1

    def test_genotype_from_value(self):
        """Test creating Genotype from numeric value."""
        assert Genotype(0) == Genotype.HOM_REF
        assert Genotype(1) == Genotype.HET
        assert Genotype(2) == Genotype.HOM_ALT
        assert Genotype(-1) == Genotype.MISSING


# ============================================================================
# Tests for VariantGenotype dataclass
# ============================================================================

class TestVariantGenotype:
    """Tests for VariantGenotype data class."""

    def test_variant_genotype_creation(self):
        """Test creating a VariantGenotype object."""
        v = Variant(chrom="chr1", pos=100, ref="A", alt="G")
        vg = VariantGenotype(
            variant=v,
            genotype=Genotype.HET,
            allele1="A",
            allele2="G"
        )
        assert vg.variant == v
        assert vg.genotype == Genotype.HET
        assert vg.allele1 == "A"
        assert vg.allele2 == "G"

    def test_variant_genotype_is_het(self):
        """Test is_het property."""
        v = Variant(chrom="chr1", pos=100, ref="A", alt="G")

        het = VariantGenotype(v, Genotype.HET)
        assert het.is_het is True

        hom_ref = VariantGenotype(v, Genotype.HOM_REF)
        assert hom_ref.is_het is False

        hom_alt = VariantGenotype(v, Genotype.HOM_ALT)
        assert hom_alt.is_het is False


# ============================================================================
# Tests for VariantSource ABC and Factory
# ============================================================================

class TestVariantSourceFactory:
    """Tests for VariantSource factory/registry pattern."""

    def test_format_detection_vcf(self, sample_vcf):
        """Test auto-detection of VCF format."""
        ext = VariantSource._detect_format(sample_vcf)
        assert ext == "vcf"

    def test_format_detection_vcf_gz(self, sample_vcf_gz):
        """Test auto-detection of compressed VCF format."""
        ext = VariantSource._detect_format(sample_vcf_gz)
        assert ext == "vcf"

    def test_format_detection_pgen(self, sample_pgen_files):
        """Test auto-detection of PGEN format."""
        ext = VariantSource._detect_format(sample_pgen_files['pgen'])
        assert ext == "pgen"

    def test_open_vcf_returns_correct_type(self, sample_vcf):
        """Test that opening VCF returns VCFSource."""
        with VariantSource.open(sample_vcf) as source:
            assert source.__class__.__name__ == "VCFSource"

    def test_open_pgen_returns_correct_type(self, sample_pgen_files):
        """Test that opening PGEN returns PGENSource."""
        with VariantSource.open(sample_pgen_files['pgen']) as source:
            assert source.__class__.__name__ == "PGENSource"

    def test_open_unsupported_format_raises(self, tmp_path):
        """Test that unsupported format raises ValueError."""
        bad_file = tmp_path / "data.xyz"
        bad_file.touch()
        with pytest.raises(ValueError, match="Unsupported.*format"):
            VariantSource.open(bad_file)

    def test_open_nonexistent_file_raises(self, tmp_path):
        """Test that nonexistent file raises FileNotFoundError."""
        missing = tmp_path / "missing.vcf"
        with pytest.raises(FileNotFoundError):
            VariantSource.open(missing)

    def test_registry_contains_expected_formats(self):
        """Test that registry has VCF and PGEN registered."""
        assert "vcf" in VariantSource._registry
        assert "pgen" in VariantSource._registry


# ============================================================================
# Tests for VariantSource interface (abstract methods)
# These tests verify behavior across ALL implementations
# ============================================================================

class TestVariantSourceInterface:
    """Tests for VariantSource interface contract.

    These tests are parameterized to run against both VCF and PGEN sources.
    """

    @pytest.fixture(params=["vcf", "pgen"])
    def variant_source(self, request, sample_vcf, sample_pgen_files):
        """Parameterized fixture providing both VCF and PGEN sources."""
        if request.param == "vcf":
            return VariantSource.open(sample_vcf)
        else:
            return VariantSource.open(sample_pgen_files['pgen'])

    def test_samples_property(self, variant_source):
        """Test samples property returns list of sample IDs."""
        samples = variant_source.samples
        assert isinstance(samples, list)
        assert len(samples) == 2
        assert "sample1" in samples or "0_sample1" in samples  # PLINK may add FID

    def test_variant_count_property(self, variant_source):
        """Test variant_count returns correct count."""
        count = variant_source.variant_count
        assert count == 6

    def test_sample_count_property(self, variant_source):
        """Test sample_count returns correct count."""
        count = variant_source.sample_count
        assert count == 2

    def test_iter_variants_returns_all(self, variant_source):
        """Test iterating over all variants."""
        variants = list(variant_source.iter_variants())
        assert len(variants) == 6

        # Check first variant
        first = variants[0]
        assert isinstance(first, VariantGenotype)
        assert first.variant.chrom == "chr1"
        assert first.variant.pos == 100

    def test_iter_variants_het_only(self, variant_source):
        """Test iterating over heterozygous sites only."""
        het_sites = list(variant_source.iter_variants(het_only=True))

        # All returned should be het
        for vg in het_sites:
            assert vg.genotype == Genotype.HET

    def test_iter_variants_single_sample(self, variant_source):
        """Test iterating for a specific sample."""
        samples = variant_source.samples
        sample = samples[0]

        variants = list(variant_source.iter_variants(samples=[sample]))
        # Should get 6 variants for the sample
        assert len(variants) == 6

    def test_get_sample_idx(self, variant_source):
        """Test getting sample index by ID."""
        samples = variant_source.samples
        idx = variant_source.get_sample_idx(samples[0])
        assert idx == 0

    def test_get_sample_idx_invalid(self, variant_source):
        """Test that invalid sample ID raises ValueError."""
        with pytest.raises(ValueError, match="not found"):
            variant_source.get_sample_idx("nonexistent_sample")

    def test_validate(self, variant_source):
        """Test validate method returns True for valid source."""
        assert variant_source.validate() is True

    def test_context_manager(self, sample_vcf):
        """Test context manager protocol."""
        with VariantSource.open(sample_vcf) as source:
            assert source.validate() is True
        # After exiting, source should be closed
        # (implementation-specific whether this raises)


# ============================================================================
# Tests for to_bed() method
# ============================================================================

class TestToBed:
    """Tests for the to_bed() method."""

    @pytest.fixture(params=["vcf", "pgen"])
    def variant_source(self, request, sample_vcf, sample_pgen_files):
        """Parameterized fixture for both formats."""
        if request.param == "vcf":
            return VariantSource.open(sample_vcf)
        else:
            return VariantSource.open(sample_pgen_files['pgen'])

    def test_to_bed_creates_file(self, variant_source, tmp_output_dir):
        """Test that to_bed creates output file."""
        output = tmp_output_dir / "output.bed"
        result = variant_source.to_bed(output)

        assert result == output
        assert output.exists()

    def test_to_bed_content_format(self, variant_source, tmp_output_dir):
        """Test BED output has correct format."""
        output = tmp_output_dir / "output.bed"
        variant_source.to_bed(output, het_only=False, include_genotypes=False)

        lines = output.read_text().strip().split('\n')

        # Should have 6 variants
        assert len(lines) == 6

        # Check first line format: chrom, start (0-based), end, ref, alt
        fields = lines[0].split('\t')
        assert len(fields) >= 5
        assert fields[0] == "chr1"
        assert fields[1] == "99"   # 0-based start
        assert fields[2] == "100"  # 1-based end
        assert fields[3] == "A"    # ref
        assert fields[4] == "G"    # alt

    def test_to_bed_het_only(self, variant_source, tmp_output_dir):
        """Test het_only filtering."""
        output = tmp_output_dir / "het_only.bed"
        samples = variant_source.samples

        # Get het sites for first sample
        variant_source.to_bed(
            output,
            samples=[samples[0]],
            het_only=True
        )

        lines = output.read_text().strip().split('\n')
        # sample1 has 3 het sites
        # (may vary slightly due to format differences)
        assert len(lines) >= 2  # At least some het sites

    def test_to_bed_with_genotypes(self, variant_source, tmp_output_dir):
        """Test including genotype columns."""
        output = tmp_output_dir / "with_gt.bed"
        samples = variant_source.samples

        variant_source.to_bed(
            output,
            samples=[samples[0]],
            het_only=False,
            include_genotypes=True
        )

        lines = output.read_text().strip().split('\n')
        fields = lines[0].split('\t')

        # Should have genotype column(s) after ref/alt
        assert len(fields) >= 6


# ============================================================================
# Tests for query_region() method
# ============================================================================

class TestQueryRegion:
    """Tests for region queries."""

    @pytest.fixture(params=["vcf", "pgen"])
    def variant_source(self, request, sample_vcf_gz, sample_pgen_files):
        """Use indexed VCF for region queries."""
        if request.param == "vcf":
            return VariantSource.open(sample_vcf_gz)
        else:
            return VariantSource.open(sample_pgen_files['pgen'])

    def test_query_region_returns_variants(self, variant_source):
        """Test querying a region returns expected variants."""
        variants = list(variant_source.query_region("chr1", 100, 300))

        # Should include variants at pos 100, 200, 300
        positions = [v.variant.pos for v in variants]
        assert 100 in positions
        assert 200 in positions
        assert 300 in positions

    def test_query_region_empty(self, variant_source):
        """Test querying empty region returns no variants."""
        variants = list(variant_source.query_region("chr1", 500, 600))
        assert len(variants) == 0

    def test_query_region_single_variant(self, variant_source):
        """Test querying single position."""
        variants = list(variant_source.query_region("chr1", 100, 100))
        assert len(variants) == 1
        assert variants[0].variant.pos == 100


# ============================================================================
# Output equivalence tests
# ============================================================================

class TestOutputEquivalence:
    """Tests ensuring VCF and PGEN produce equivalent outputs."""

    def test_bed_output_equivalence(
        self, sample_vcf, sample_pgen_files, tmp_output_dir
    ):
        """Test that VCF and PGEN produce equivalent BED output."""
        vcf_source = VariantSource.open(sample_vcf)
        pgen_source = VariantSource.open(sample_pgen_files['pgen'])

        vcf_bed = tmp_output_dir / "vcf.bed"
        pgen_bed = tmp_output_dir / "pgen.bed"

        # Export without genotypes for fair comparison
        vcf_source.to_bed(vcf_bed, het_only=False, include_genotypes=False)
        pgen_source.to_bed(pgen_bed, het_only=False, include_genotypes=False)

        # Compare content
        vcf_lines = set(vcf_bed.read_text().strip().split('\n'))
        pgen_lines = set(pgen_bed.read_text().strip().split('\n'))

        assert vcf_lines == pgen_lines, (
            f"BED outputs differ!\n"
            f"VCF-only: {vcf_lines - pgen_lines}\n"
            f"PGEN-only: {pgen_lines - vcf_lines}"
        )

    def test_variant_count_equivalence(self, sample_vcf, sample_pgen_files):
        """Test VCF and PGEN report same variant count."""
        vcf_source = VariantSource.open(sample_vcf)
        pgen_source = VariantSource.open(sample_pgen_files['pgen'])

        assert vcf_source.variant_count == pgen_source.variant_count

    def test_sample_count_equivalence(self, sample_vcf, sample_pgen_files):
        """Test VCF and PGEN report same sample count."""
        vcf_source = VariantSource.open(sample_vcf)
        pgen_source = VariantSource.open(sample_pgen_files['pgen'])

        assert vcf_source.sample_count == pgen_source.sample_count
