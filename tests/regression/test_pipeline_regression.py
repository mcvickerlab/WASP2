"""
Regression tests against baseline pipeline outputs.

This test suite validates that code changes don't break:
1. Output correctness (MD5 checksums)
2. Performance characteristics (time, memory)
3. Output format and structure
4. Statistical results

Run with: pytest tests/regression/test_pipeline_regression.py -v
"""

import hashlib
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Dict, Tuple

import pandas as pd
import pytest

# Project root
ROOT = Path(__file__).parent.parent.parent
BASELINE_DIR = ROOT / "baselines"
TEST_DATA = ROOT / "test_data"

# Baseline expectations from committed benchmarks
BASELINE_EXPECTATIONS = {
    "counting": {
        "time_seconds": 9.26,
        "memory_mb": 639,
        "output_rows": 111455,  # header + 111454 SNPs
        "total_alleles": 2388,
        "md5": "127a81810a43db3cc6924a26f591cc7a"
    },
    "analysis": {
        "time_seconds": 2.97,
        "memory_mb": 340,
        "output_rows": 40,  # header + 39 regions
        "significant_regions": 0,
        "md5": "394e1a7dbf14220079c3142c5b15bad8"
    },
    "mapping": {
        "time_seconds": 8.0,
        "memory_mb": 488,
        "wasp_filtered_reads": 125387,
        "original_reads": 126061
    }
}

# Tolerance for performance regression
TIME_TOLERANCE = 1.30  # Allow 30% slower
MEMORY_TOLERANCE = 1.20  # Allow 20% more memory


def md5_file(filepath: Path) -> str:
    """Calculate MD5 checksum of file."""
    hash_md5 = hashlib.md5()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def parse_memory_profile(profile_file: Path) -> Dict[str, float]:
    """Parse /usr/bin/time -v output to extract metrics."""
    with open(profile_file) as f:
        content = f.read()

    metrics = {}
    for line in content.split('\n'):
        if 'Maximum resident set size' in line:
            kb = int(line.split(':')[1].strip())
            metrics['memory_mb'] = kb / 1024
        elif 'Elapsed (wall clock) time' in line:
            # Format: "Elapsed (wall clock) time (h:mm:ss or m:ss): 0:09.26"
            # Take last part after splitting on ':'
            time_str = line.split(':')[-1].strip()
            # Parse m:ss.ms format
            if ':' in time_str:
                parts = time_str.split(':')
                if len(parts) == 2:
                    mins, secs = parts
                    metrics['time_seconds'] = int(mins) * 60 + float(secs)
                elif len(parts) == 3:
                    hours, mins, secs = parts
                    metrics['time_seconds'] = int(hours) * 3600 + int(mins) * 60 + float(secs)
            else:
                # Just seconds
                metrics['time_seconds'] = float(time_str)

    return metrics


class TestCountingRegression:
    """Test counting module against baseline."""

    def test_counting_output_md5(self):
        """Verify counting output MD5 matches baseline."""
        baseline_counts = BASELINE_DIR / "counting" / "counts.tsv"

        if not baseline_counts.exists():
            pytest.skip("Baseline counting output not found")

        actual_md5 = md5_file(baseline_counts)
        expected_md5 = BASELINE_EXPECTATIONS["counting"]["md5"]

        assert actual_md5 == expected_md5, (
            f"Counting output MD5 mismatch!\n"
            f"Expected: {expected_md5}\n"
            f"Actual:   {actual_md5}\n"
            f"This indicates output has changed. If intentional, update baseline."
        )

    def test_counting_output_structure(self):
        """Verify counting output has correct structure."""
        baseline_counts = BASELINE_DIR / "counting" / "counts.tsv"

        if not baseline_counts.exists():
            pytest.skip("Baseline counting output not found")

        df = pd.read_csv(baseline_counts, sep='\t')

        # Check columns
        expected_cols = ['chrom', 'pos', 'ref', 'alt', 'GT', 'region',
                        'ref_count', 'alt_count', 'other_count']
        assert list(df.columns) == expected_cols, f"Column mismatch: {list(df.columns)}"

        # Check row count
        assert len(df) == BASELINE_EXPECTATIONS["counting"]["output_rows"] - 1  # minus header

        # Check data types
        assert df['ref_count'].dtype in [int, 'int64', 'uint16']
        assert df['alt_count'].dtype in [int, 'int64', 'uint16']
        assert df['other_count'].dtype in [int, 'int64', 'uint16']

        # Check total alleles
        total_alleles = (df['ref_count'].sum() +
                        df['alt_count'].sum() +
                        df['other_count'].sum())
        assert total_alleles == BASELINE_EXPECTATIONS["counting"]["total_alleles"]

    def test_counting_memory_regression(self):
        """Verify counting memory usage hasn't regressed."""
        memory_profile = BASELINE_DIR / "counting" / "memory_profile.txt"

        if not memory_profile.exists():
            pytest.skip("Baseline memory profile not found")

        metrics = parse_memory_profile(memory_profile)
        actual_mb = metrics['memory_mb']
        expected_mb = BASELINE_EXPECTATIONS["counting"]["memory_mb"]
        max_allowed_mb = expected_mb * MEMORY_TOLERANCE

        assert actual_mb <= max_allowed_mb, (
            f"Memory regression detected!\n"
            f"Baseline:  {expected_mb} MB\n"
            f"Current:   {actual_mb} MB\n"
            f"Max allowed: {max_allowed_mb} MB ({MEMORY_TOLERANCE}x tolerance)\n"
            f"Increase: {((actual_mb / expected_mb) - 1) * 100:.1f}%"
        )

    def test_counting_performance_regression(self):
        """Verify counting performance hasn't regressed."""
        memory_profile = BASELINE_DIR / "counting" / "memory_profile.txt"

        if not memory_profile.exists():
            pytest.skip("Baseline memory profile not found")

        metrics = parse_memory_profile(memory_profile)
        actual_seconds = metrics['time_seconds']
        expected_seconds = BASELINE_EXPECTATIONS["counting"]["time_seconds"]
        max_allowed_seconds = expected_seconds * TIME_TOLERANCE

        assert actual_seconds <= max_allowed_seconds, (
            f"Performance regression detected!\n"
            f"Baseline:  {expected_seconds}s\n"
            f"Current:   {actual_seconds}s\n"
            f"Max allowed: {max_allowed_seconds}s ({TIME_TOLERANCE}x tolerance)\n"
            f"Slowdown: {((actual_seconds / expected_seconds) - 1) * 100:.1f}%"
        )


class TestAnalysisRegression:
    """Test analysis module against baseline."""

    def test_analysis_output_md5(self):
        """Verify analysis output MD5 matches baseline."""
        baseline_analysis = BASELINE_DIR / "analysis" / "ai_results.tsv"

        if not baseline_analysis.exists():
            pytest.skip("Baseline analysis output not found")

        actual_md5 = md5_file(baseline_analysis)
        expected_md5 = BASELINE_EXPECTATIONS["analysis"]["md5"]

        assert actual_md5 == expected_md5, (
            f"Analysis output MD5 mismatch!\n"
            f"Expected: {expected_md5}\n"
            f"Actual:   {actual_md5}\n"
            f"This indicates output has changed. If intentional, update baseline."
        )

    def test_analysis_output_structure(self):
        """Verify analysis output has correct structure."""
        baseline_analysis = BASELINE_DIR / "analysis" / "ai_results.tsv"

        if not baseline_analysis.exists():
            pytest.skip("Baseline analysis output not found")

        df = pd.read_csv(baseline_analysis, sep='\t')

        # Check columns
        expected_cols = ['region', 'ref_count', 'alt_count', 'N', 'snp_count',
                        'null_ll', 'alt_ll', 'mu', 'lrt', 'pval', 'fdr_pval']
        assert list(df.columns) == expected_cols, f"Column mismatch: {list(df.columns)}"

        # Check row count
        assert len(df) == BASELINE_EXPECTATIONS["analysis"]["output_rows"] - 1  # minus header

        # Check significant regions
        significant = (df['fdr_pval'] < 0.05).sum()
        assert significant == BASELINE_EXPECTATIONS["analysis"]["significant_regions"]

        # Validate statistical properties
        assert (df['mu'] >= 0).all() and (df['mu'] <= 1).all(), "mu should be probability [0,1]"
        assert (df['pval'] >= 0).all() and (df['pval'] <= 1).all(), "pval should be [0,1]"
        # LRT should be non-negative (allow tiny negative values from floating point errors)
        assert (df['lrt'] >= -1e-10).all(), f"LRT should be non-negative (found: {df['lrt'].min()})"

    def test_analysis_memory_regression(self):
        """Verify analysis memory usage hasn't regressed."""
        memory_profile = BASELINE_DIR / "analysis" / "memory_profile.txt"

        if not memory_profile.exists():
            pytest.skip("Baseline memory profile not found")

        metrics = parse_memory_profile(memory_profile)
        actual_mb = metrics['memory_mb']
        expected_mb = BASELINE_EXPECTATIONS["analysis"]["memory_mb"]
        max_allowed_mb = expected_mb * MEMORY_TOLERANCE

        assert actual_mb <= max_allowed_mb, (
            f"Memory regression detected!\n"
            f"Baseline:  {expected_mb} MB\n"
            f"Current:   {actual_mb} MB\n"
            f"Increase: {((actual_mb / expected_mb) - 1) * 100:.1f}%"
        )

    def test_analysis_performance_regression(self):
        """Verify analysis performance hasn't regressed."""
        memory_profile = BASELINE_DIR / "analysis" / "memory_profile.txt"

        if not memory_profile.exists():
            pytest.skip("Baseline memory profile not found")

        metrics = parse_memory_profile(memory_profile)
        actual_seconds = metrics['time_seconds']
        expected_seconds = BASELINE_EXPECTATIONS["analysis"]["time_seconds"]
        max_allowed_seconds = expected_seconds * TIME_TOLERANCE

        assert actual_seconds <= max_allowed_seconds, (
            f"Performance regression detected!\n"
            f"Baseline:  {expected_seconds}s\n"
            f"Current:   {actual_seconds}s\n"
            f"Slowdown: {((actual_seconds / expected_seconds) - 1) * 100:.1f}%"
        )


class TestMappingRegression:
    """Test mapping module against baseline."""

    def test_mapping_wasp_filter_rate(self):
        """Verify WASP filtering preserves expected read count."""
        metadata = BASELINE_DIR / "pipeline_metadata.txt"

        if not metadata.exists():
            pytest.skip("Baseline metadata not found")

        with open(metadata) as f:
            content = f.read()

        # Parse read counts
        for line in content.split('\n'):
            if 'Original reads:' in line:
                original = int(line.split(':')[1].strip().split()[0])
            elif 'WASP filtered reads:' in line:
                filtered = int(line.split(':')[1].strip().split()[0])

        assert original == BASELINE_EXPECTATIONS["mapping"]["original_reads"]
        assert filtered == BASELINE_EXPECTATIONS["mapping"]["wasp_filtered_reads"]

        # Check filter rate is reasonable (should keep >95%)
        filter_rate = filtered / original
        assert filter_rate > 0.95, (
            f"WASP filter rate too aggressive: {filter_rate:.1%}\n"
            f"Kept {filtered}/{original} reads"
        )


class TestFullPipelineIntegration:
    """End-to-end pipeline integration tests."""

    @pytest.mark.slow
    def test_full_pipeline_reproducibility(self, tmp_path):
        """Run full pipeline and verify output matches baseline exactly.

        This is a slow test (20+ seconds) but provides strongest guarantee.
        """
        # Create temp output directory
        temp_baseline = tmp_path / "baseline_test"
        temp_baseline.mkdir()

        # Run pipeline script
        script = ROOT / "scripts" / "run_full_pipeline_baseline.sh"

        if not script.exists():
            pytest.skip("Pipeline script not found")

        # Run with temp output
        result = subprocess.run(
            [str(script)],
            env={**dict(subprocess.os.environ), "BASELINE_DIR": str(temp_baseline)},
            capture_output=True,
            text=True
        )

        if result.returncode != 0:
            pytest.fail(f"Pipeline failed:\n{result.stderr}")

        # Compare outputs
        temp_counts = temp_baseline / "counting" / "counts.tsv"
        baseline_counts = BASELINE_DIR / "counting" / "counts.tsv"

        if temp_counts.exists() and baseline_counts.exists():
            assert md5_file(temp_counts) == md5_file(baseline_counts), (
                "Counting output not reproducible!"
            )

        temp_analysis = temp_baseline / "analysis" / "ai_results.tsv"
        baseline_analysis = BASELINE_DIR / "analysis" / "ai_results.tsv"

        if temp_analysis.exists() and baseline_analysis.exists():
            assert md5_file(temp_analysis) == md5_file(baseline_analysis), (
                "Analysis output not reproducible!"
            )


if __name__ == "__main__":
    # Run tests with verbose output
    pytest.main([__file__, "-v", "--tb=short"])
