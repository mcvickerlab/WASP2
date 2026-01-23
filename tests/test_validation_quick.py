#!/usr/bin/env python3
"""
Quick validation tests for WASP2 pipeline.

These tests validate:
1. Unit tests pass (Rust vs Python parity)
2. INDEL correctness tests pass
3. Module imports work correctly

Run with: pytest tests/test_validation_quick.py -v
"""
import pytest
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"

if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


class TestQuickValidation:
    """Quick validation tests that don't require large test data."""

    def test_rust_module_imports(self):
        """Test that Rust module can be imported."""
        try:
            import wasp2_rust
            assert hasattr(wasp2_rust, 'remap_all_chromosomes')
            assert hasattr(wasp2_rust, 'filter_bam_rust')
        except ImportError as e:
            pytest.skip(f"Rust module not available: {e}")

    def test_python_module_imports(self):
        """Test that Python modules can be imported."""
        from mapping import run_mapping
        from counting import run_counting
        from wasp2.io import vcf_source
        assert callable(run_mapping.make_reads_pipeline)

    def test_rust_python_parity(self):
        """Run the Rust vs Python parity tests."""
        test_file = ROOT / "tests" / "test_rust_python_match.py"
        if not test_file.exists():
            pytest.skip("test_rust_python_match.py not found")

        result = subprocess.run(
            [sys.executable, "-m", "pytest", str(test_file), "-v", "--tb=short"],
            capture_output=True,
            text=True,
            cwd=ROOT
        )

        if result.returncode != 0:
            print(result.stdout)
            print(result.stderr)

        assert result.returncode == 0, f"Rust/Python parity tests failed:\n{result.stdout}\n{result.stderr}"

    def test_indel_correctness(self):
        """Run the INDEL correctness tests."""
        test_file = ROOT / "tests" / "test_indel_correctness.py"
        if not test_file.exists():
            pytest.skip("test_indel_correctness.py not found")

        result = subprocess.run(
            [sys.executable, "-m", "pytest", str(test_file), "-v", "--tb=short"],
            capture_output=True,
            text=True,
            cwd=ROOT
        )

        if result.returncode != 0:
            print(result.stdout)
            print(result.stderr)

        assert result.returncode == 0, f"INDEL correctness tests failed:\n{result.stdout}\n{result.stderr}"


class TestExpectedCounts:
    """Tests that validate expected pipeline output counts."""

    EXPECTED_COUNTS_FILE = ROOT / "baselines" / "mapping" / "expected_counts.json"

    def test_expected_counts_file_exists(self):
        """Verify expected counts baseline file exists."""
        assert self.EXPECTED_COUNTS_FILE.exists(), \
            f"Expected counts file not found: {self.EXPECTED_COUNTS_FILE}"

    def test_expected_counts_structure(self):
        """Verify expected counts file has correct structure."""
        import json

        if not self.EXPECTED_COUNTS_FILE.exists():
            pytest.skip("Expected counts file not found")

        with open(self.EXPECTED_COUNTS_FILE) as f:
            data = json.load(f)

        # Check required fields
        assert "expected_counts" in data
        counts = data["expected_counts"]

        required_fields = [
            "vcf_variants",
            "r1_fastq_reads",
            "r2_fastq_reads",
            "total_haplotypes"
        ]

        for field in required_fields:
            assert field in counts, f"Missing required field: {field}"
            assert isinstance(counts[field], int), f"{field} should be an integer"
            assert counts[field] > 0, f"{field} should be > 0"

    def test_fastq_count_consistency(self):
        """Verify R1 and R2 FASTQ counts match."""
        import json

        if not self.EXPECTED_COUNTS_FILE.exists():
            pytest.skip("Expected counts file not found")

        with open(self.EXPECTED_COUNTS_FILE) as f:
            data = json.load(f)

        counts = data["expected_counts"]
        assert counts["r1_fastq_reads"] == counts["r2_fastq_reads"], \
            "R1 and R2 FASTQ read counts should match for paired-end data"

    def test_haplotype_count_consistency(self):
        """Verify total haplotypes = 2 * FASTQ reads."""
        import json

        if not self.EXPECTED_COUNTS_FILE.exists():
            pytest.skip("Expected counts file not found")

        with open(self.EXPECTED_COUNTS_FILE) as f:
            data = json.load(f)

        counts = data["expected_counts"]
        expected_haps = counts["r1_fastq_reads"] * 2
        assert counts["total_haplotypes"] == expected_haps, \
            f"Total haplotypes ({counts['total_haplotypes']}) should be 2 * R1 reads ({expected_haps})"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
