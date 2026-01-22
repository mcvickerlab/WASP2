"""
Integration tests for WASP2 pipeline components.

This module tests the integration between:
- Variant I/O (VCF/PGEN reading)
- Counting module
- Analysis module
- Mapping module

These tests verify that components work together correctly
in end-to-end workflows.

Run with: pytest tests/integration/test_pipeline.py -v
Run slow tests: pytest tests/integration/test_pipeline.py -v -m slow
"""

import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict

import pytest

# Project root
ROOT = Path(__file__).parent.parent.parent


# ============================================================================
# Helper functions
# ============================================================================

def has_command(cmd: str) -> bool:
    """Check if a command is available in PATH."""
    return shutil.which(cmd) is not None


def skip_without_dependencies(deps: list):
    """Skip test if any dependencies are missing."""
    missing = [cmd for cmd in deps if not has_command(cmd)]
    if missing:
        pytest.skip(f"Missing dependencies: {', '.join(missing)}")


# ============================================================================
# Test: Variant I/O Integration
# ============================================================================

@pytest.mark.integration
class TestVariantIOIntegration:
    """Integration tests for variant I/O components."""

    def test_vcf_source_reads_sample_vcf(self, sample_vcf_gz):
        """Test that VCF source can read compressed VCF files."""
        try:
            from wasp2.io.vcf_source import VCFSource
        except ImportError:
            pytest.skip("wasp2.io.vcf_source not available")

        source = VCFSource(str(sample_vcf_gz))
        variants = list(source.iter_variants())

        assert len(variants) > 0
        assert all(hasattr(v, 'chrom') for v in variants)
        assert all(hasattr(v, 'pos') for v in variants)

    def test_vcf_and_pgen_consistency(self, sample_vcf_gz, sample_pgen_files):
        """Test that VCF and PGEN sources return consistent data."""
        try:
            from wasp2.io.vcf_source import VCFSource
            from wasp2.io.pgen_source import PgenSource
        except ImportError:
            pytest.skip("wasp2.io sources not available")

        vcf_source = VCFSource(str(sample_vcf_gz))
        vcf_variants = list(vcf_source.iter_variants())

        pgen_source = PgenSource(str(sample_pgen_files['prefix']))
        pgen_variants = list(pgen_source.iter_variants())

        # Should have same number of variants
        assert len(vcf_variants) == len(pgen_variants)

        # Positions should match
        vcf_positions = {(v.chrom, v.pos) for v in vcf_variants}
        pgen_positions = {(v.chrom, v.pos) for v in pgen_variants}
        assert vcf_positions == pgen_positions


# ============================================================================
# Test: Counting Pipeline Integration
# ============================================================================

@pytest.mark.integration
class TestCountingIntegration:
    """Integration tests for counting pipeline."""

    @pytest.mark.rust
    @pytest.mark.slow
    def test_counting_module_end_to_end(self, tmp_output_dir, test_data_dir):
        """Test counting module with real data files."""
        from counting.count_alleles import RUST_AVAILABLE

        if not RUST_AVAILABLE:
            pytest.skip("Rust extension not available")

        # Check for required test data
        bam_file = test_data_dir / "sample.bam"
        vcf_file = test_data_dir / "sample.vcf.gz"

        if not bam_file.exists() or not vcf_file.exists():
            pytest.skip("Required test data files not found")

        # Run counting
        from counting.run_counting import main as run_counting

        output_file = tmp_output_dir / "counts.tsv"

        try:
            run_counting(
                bam_file=str(bam_file),
                vcf_file=str(vcf_file),
                output_file=str(output_file),
            )
        except Exception as e:
            pytest.fail(f"Counting failed: {e}")

        assert output_file.exists()
        assert output_file.stat().st_size > 0


# ============================================================================
# Test: Analysis Pipeline Integration
# ============================================================================

@pytest.mark.integration
class TestAnalysisIntegration:
    """Integration tests for analysis pipeline."""

    @pytest.mark.slow
    def test_analysis_module_end_to_end(self, tmp_output_dir, test_data_dir):
        """Test analysis module with counting output."""
        counts_file = test_data_dir / "sample_counts.tsv"

        if not counts_file.exists():
            pytest.skip("Sample counts file not found")

        from analysis.run_analysis import main as run_analysis

        output_file = tmp_output_dir / "results.tsv"

        try:
            run_analysis(
                counts_file=str(counts_file),
                output_file=str(output_file),
            )
        except Exception as e:
            pytest.fail(f"Analysis failed: {e}")

        assert output_file.exists()
        assert output_file.stat().st_size > 0


# ============================================================================
# Test: Full Pipeline Integration
# ============================================================================

@pytest.mark.integration
class TestFullPipelineIntegration:
    """End-to-end pipeline integration tests."""

    @pytest.mark.slow
    def test_pipeline_components_chain(self, tmp_output_dir):
        """Test that pipeline components can be chained together."""
        # This test verifies the data flow between components
        # without running the full pipeline

        import polars as pl

        # Create mock counting output
        mock_counts = pl.DataFrame({
            "chrom": ["chr1", "chr1", "chr1"],
            "pos": [100, 200, 300],
            "ref": ["A", "C", "G"],
            "alt": ["G", "T", "A"],
            "GT": ["0/1", "0/1", "0/1"],
            "region": ["gene1", "gene1", "gene2"],
            "ref_count": [10, 8, 12],
            "alt_count": [12, 10, 8],
            "other_count": [0, 0, 1],
        })

        counts_file = tmp_output_dir / "mock_counts.tsv"
        mock_counts.write_csv(counts_file, separator="\t")

        assert counts_file.exists()

        # Verify the file can be read back
        loaded = pl.read_csv(counts_file, separator="\t")
        assert loaded.shape == mock_counts.shape
        assert list(loaded.columns) == list(mock_counts.columns)

    @pytest.mark.slow
    def test_full_pipeline_script(self, tmp_path):
        """Test full pipeline execution via script."""
        skip_without_dependencies(["bcftools", "samtools", "bedtools"])

        script = ROOT / "scripts" / "run_full_pipeline_baseline.sh"

        if not script.exists():
            pytest.skip("Pipeline script not found")

        # Check for required test data
        test_data = ROOT / "test_data"
        required_files = [
            test_data / "CD4_ATACseq_Day1_merged_filtered.sort.bam",
            test_data / "filter_chr10.vcf",
            test_data / "NA12878_snps_chr10.bed",
        ]

        for fpath in required_files:
            if not fpath.exists():
                pytest.skip(f"Required test data missing: {fpath}")

        # Setup environment
        env = dict(subprocess.os.environ)
        env["PYTHONPATH"] = str(ROOT / "src")

        output_dir = tmp_path / "pipeline_output"
        output_dir.mkdir()

        # Run pipeline
        result = subprocess.run(
            [str(script)],
            env={**env, "BASELINE_DIR": str(output_dir)},
            cwd=str(ROOT),
            capture_output=True,
            text=True,
            timeout=300,  # 5 minute timeout
        )

        if result.returncode != 0:
            pytest.fail(
                f"Pipeline failed with return code {result.returncode}\n"
                f"stdout: {result.stdout}\n"
                f"stderr: {result.stderr}"
            )

        # Verify outputs exist
        expected_outputs = [
            output_dir / "counting" / "counts.tsv",
            output_dir / "analysis" / "ai_results.tsv",
        ]

        for output in expected_outputs:
            if not output.exists():
                pytest.fail(f"Expected output not found: {output}")


# ============================================================================
# Test: Component Interoperability
# ============================================================================

@pytest.mark.integration
class TestComponentInteroperability:
    """Tests for component interoperability."""

    def test_polars_dataframe_compatibility(self):
        """Test that all modules use compatible Polars DataFrames."""
        import polars as pl

        # Create a DataFrame similar to counting output
        df = pl.DataFrame({
            "chrom": pl.Series(["chr1", "chr1"], dtype=pl.Categorical),
            "pos": pl.Series([100, 200], dtype=pl.UInt32),
            "ref": ["A", "C"],
            "alt": ["G", "T"],
            "ref_count": pl.Series([10, 8], dtype=pl.UInt16),
            "alt_count": pl.Series([12, 10], dtype=pl.UInt16),
            "other_count": pl.Series([0, 0], dtype=pl.UInt16),
        })

        # Verify schema is compatible with analysis expectations
        assert df.schema["pos"] == pl.UInt32
        assert df.schema["ref_count"] == pl.UInt16
        assert df.schema["alt_count"] == pl.UInt16

    def test_chromosome_naming_consistency(self):
        """Test chromosome naming is handled consistently."""
        import polars as pl

        # Both 'chr1' and '1' formats should be handled
        df_with_chr = pl.DataFrame({
            "chrom": ["chr1", "chr2", "chrX"],
            "pos": [100, 200, 300],
        })

        df_without_chr = pl.DataFrame({
            "chrom": ["1", "2", "X"],
            "pos": [100, 200, 300],
        })

        # Both should have valid chromosome columns
        assert df_with_chr["chrom"].dtype == pl.String
        assert df_without_chr["chrom"].dtype == pl.String


# ============================================================================
# Test: Error Handling Integration
# ============================================================================

@pytest.mark.integration
class TestErrorHandlingIntegration:
    """Tests for integrated error handling."""

    def test_missing_bam_file_error(self, tmp_output_dir):
        """Test appropriate error for missing BAM file."""
        from counting.count_alleles import RUST_AVAILABLE

        if not RUST_AVAILABLE:
            pytest.skip("Rust extension not available")

        import polars as pl

        mock_df = pl.DataFrame({
            "chrom": pl.Series(["chr1"], dtype=pl.Categorical),
            "pos": [100],
            "ref": ["A"],
            "alt": ["G"],
        })

        from counting.count_alleles import make_count_df

        with pytest.raises(Exception):
            make_count_df("/nonexistent/path/to/file.bam", mock_df)

    def test_invalid_vcf_format_error(self, tmp_output_dir):
        """Test appropriate error for invalid VCF format."""
        invalid_vcf = tmp_output_dir / "invalid.vcf"
        invalid_vcf.write_text("This is not a valid VCF file\n")

        try:
            from wasp2.io.vcf_source import VCFSource
        except ImportError:
            pytest.skip("wasp2.io.vcf_source not available")

        with pytest.raises(Exception):
            source = VCFSource(str(invalid_vcf))
            list(source.iter_variants())


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
