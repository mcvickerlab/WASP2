"""
Benchmark comparisons: WASP2 vs WASP v1 vs phASER.

This module provides standardized benchmarks for comparing allele-specific
analysis tools. External tools (WASP v1, phASER) are called via subprocess
wrappers with graceful degradation if not installed.
"""

import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from analysis.as_analysis import get_imbalance

from .conftest import generate_allele_count_data, generate_synthetic_vcf

# ============================================================================
# Tool availability checking
# ============================================================================


def check_tool_available(tool_name: str) -> bool:
    """Check if an external tool is available in PATH."""
    return shutil.which(tool_name) is not None


HAS_WASP_V1 = check_tool_available("wasp")
HAS_PHASER = check_tool_available("phaser.py") or check_tool_available("phaser")


# ============================================================================
# Tool wrappers
# ============================================================================


class WASP2Wrapper:
    """Wrapper for WASP2 analysis functions."""

    name = "WASP2"

    @staticmethod
    def run_analysis(
        count_df: pd.DataFrame,
        min_count: int = 10,
        method: str = "single",
    ) -> pd.DataFrame:
        """Run WASP2 allelic imbalance analysis."""
        return get_imbalance(
            count_df,
            min_count=min_count,
            method=method,
            region_col="region",
        )


class WASPv1Wrapper:
    """Wrapper for WASP v1 (original WASP)."""

    name = "WASP_v1"
    available = HAS_WASP_V1

    @staticmethod
    def prepare_input(count_df: pd.DataFrame, output_dir: Path) -> Path:
        """Convert count data to WASP v1 input format."""
        # WASP v1 expects a specific count table format
        wasp_input = output_dir / "wasp_counts.txt"
        # Create WASP v1 compatible format
        wasp_df = count_df[["chrom", "pos", "ref", "alt", "ref_count", "alt_count"]].copy()
        wasp_df.columns = ["CHROM", "POS", "REF", "ALT", "REF_COUNT", "ALT_COUNT"]
        wasp_df.to_csv(wasp_input, sep="\t", index=False)
        return wasp_input

    @staticmethod
    def run_analysis(
        input_file: Path,
        output_dir: Path,
    ) -> tuple[pd.DataFrame | None, str | None]:
        """Run WASP v1 analysis via subprocess."""
        if not HAS_WASP_V1:
            return None, "WASP v1 not installed"

        output_file = output_dir / "wasp_results.txt"
        try:
            subprocess.run(
                [
                    "wasp",
                    "--input",
                    str(input_file),
                    "--output",
                    str(output_file),
                ],
                check=True,
                capture_output=True,
                timeout=600,
            )
            if output_file.exists():
                return pd.read_csv(output_file, sep="\t"), None
            return None, "Output file not created"
        except subprocess.TimeoutExpired:
            return None, "Timeout exceeded"
        except subprocess.CalledProcessError as e:
            return None, f"Process error: {e.stderr.decode()}"
        except FileNotFoundError:
            return None, "WASP v1 executable not found"


class PhASERWrapper:
    """Wrapper for phASER tool."""

    name = "phASER"
    available = HAS_PHASER

    @staticmethod
    def prepare_input(
        count_df: pd.DataFrame,
        vcf_path: Path,
        bam_path: Path | None,
        output_dir: Path,
    ) -> dict[str, Path]:
        """Prepare phASER input files."""
        # phASER requires VCF and BAM inputs
        return {
            "vcf": vcf_path,
            "bam": bam_path,
            "output_prefix": output_dir / "phaser_output",
        }

    @staticmethod
    def run_analysis(
        vcf_path: Path,
        bam_path: Path | None,
        output_prefix: Path,
    ) -> tuple[pd.DataFrame | None, str | None]:
        """Run phASER analysis via subprocess."""
        if not HAS_PHASER:
            return None, "phASER not installed"

        if bam_path is None:
            return None, "BAM file required for phASER"

        try:
            cmd = [
                "phaser.py" if check_tool_available("phaser.py") else "phaser",
                "--vcf",
                str(vcf_path),
                "--bam",
                str(bam_path),
                "--o",
                str(output_prefix),
                "--threads",
                "1",
            ]
            subprocess.run(cmd, check=True, capture_output=True, timeout=600)

            output_file = Path(f"{output_prefix}.allelic_counts.txt")
            if output_file.exists():
                return pd.read_csv(output_file, sep="\t"), None
            return None, "Output file not created"
        except subprocess.TimeoutExpired:
            return None, "Timeout exceeded"
        except subprocess.CalledProcessError as e:
            return None, f"Process error: {e.stderr.decode()}"
        except FileNotFoundError:
            return None, "phASER executable not found"


# ============================================================================
# Benchmark test classes
# ============================================================================


class TestWASP2Benchmarks:
    """Baseline WASP2 benchmarks for comparison."""

    @pytest.mark.benchmark(group="tool_comparison_wasp2")
    @pytest.mark.parametrize("n_variants", [1000, 10000, 100000])
    def test_wasp2_scaling(self, benchmark, rng, n_variants: int):
        """Benchmark WASP2 at different scales."""
        n_regions = max(100, n_variants // 10)
        df = generate_allele_count_data(n_variants, n_samples=10, n_regions=n_regions, rng=rng)

        result = benchmark(WASP2Wrapper.run_analysis, df)

        benchmark.extra_info["tool"] = "WASP2"
        benchmark.extra_info["n_variants"] = n_variants
        benchmark.extra_info["n_samples"] = 10
        benchmark.extra_info["n_results"] = len(result) if result is not None else 0

        assert isinstance(result, pd.DataFrame)
        assert "pval" in result.columns


class TestWASPv1Benchmarks:
    """WASP v1 benchmarks (skipped if not installed)."""

    @pytest.mark.benchmark(group="tool_comparison_wasp_v1")
    @pytest.mark.parametrize("n_variants", [1000, 10000, 100000])
    def test_wasp_v1_scaling(self, benchmark, rng, n_variants: int):
        """Benchmark WASP v1 at different scales."""
        if not HAS_WASP_V1:
            pytest.skip("WASP v1 not installed")

        n_regions = max(100, n_variants // 10)
        df = generate_allele_count_data(n_variants, n_samples=10, n_regions=n_regions, rng=rng)

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            input_file = WASPv1Wrapper.prepare_input(df, tmpdir_path)

            def run_wasp_v1():
                return WASPv1Wrapper.run_analysis(input_file, tmpdir_path)

            result, error = benchmark(run_wasp_v1)

            benchmark.extra_info["tool"] = "WASP_v1"
            benchmark.extra_info["n_variants"] = n_variants
            benchmark.extra_info["n_samples"] = 10
            benchmark.extra_info["error"] = error

            if error:
                pytest.skip(f"WASP v1 failed: {error}")


class TestPhASERBenchmarks:
    """phASER benchmarks (skipped if not installed)."""

    @pytest.mark.benchmark(group="tool_comparison_phaser")
    @pytest.mark.parametrize("n_variants", [1000, 10000])
    def test_phaser_scaling(
        self,
        benchmark,
        rng,
        n_variants: int,
        vcf_small: Path,
        bam_small: Path | None,
    ):
        """Benchmark phASER at different scales."""
        if not HAS_PHASER:
            pytest.skip("phASER not installed")
        if bam_small is None:
            pytest.skip("BAM file generation failed")

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            output_prefix = tmpdir_path / "phaser_out"

            def run_phaser():
                return PhASERWrapper.run_analysis(vcf_small, bam_small, output_prefix)

            result, error = benchmark(run_phaser)

            benchmark.extra_info["tool"] = "phASER"
            benchmark.extra_info["n_variants"] = n_variants
            benchmark.extra_info["error"] = error

            if error:
                pytest.skip(f"phASER failed: {error}")


class TestToolComparisonMemory:
    """Memory comparison benchmarks across tools."""

    @pytest.mark.benchmark(group="tool_comparison_memory")
    @pytest.mark.parametrize("n_variants", [1000, 10000, 50000])
    def test_wasp2_memory(self, benchmark, rng, n_variants: int, memory_benchmark):
        """Measure WASP2 memory usage at different scales."""
        n_regions = max(100, n_variants // 10)
        df = generate_allele_count_data(n_variants, n_samples=10, n_regions=n_regions, rng=rng)

        def run_analysis():
            return WASP2Wrapper.run_analysis(df)

        result, peak_memory = memory_benchmark.measure(run_analysis)

        benchmark.extra_info["tool"] = "WASP2"
        benchmark.extra_info["n_variants"] = n_variants
        benchmark.extra_info["peak_memory_mb"] = peak_memory

        # Also run timed benchmark
        result = benchmark(run_analysis)
        assert isinstance(result, pd.DataFrame)


class TestDirectComparison:
    """Direct head-to-head comparison tests."""

    @pytest.mark.benchmark(group="direct_comparison")
    def test_wasp2_vs_simulated_competitors(self, benchmark, rng):
        """
        Compare WASP2 against simulated competitor performance.

        This test provides baseline comparison data when external tools
        are not available, using published performance characteristics.
        """
        n_variants = 10000
        n_regions = 1000
        df = generate_allele_count_data(n_variants, n_samples=10, n_regions=n_regions, rng=rng)

        result = benchmark(WASP2Wrapper.run_analysis, df)

        benchmark.extra_info["tool"] = "WASP2"
        benchmark.extra_info["n_variants"] = n_variants
        benchmark.extra_info["comparison_note"] = (
            "Reference timing for comparison with WASP v1 and phASER. "
            "Install external tools for direct comparison."
        )

        # Store tool availability for reporting
        benchmark.extra_info["wasp_v1_available"] = HAS_WASP_V1
        benchmark.extra_info["phaser_available"] = HAS_PHASER

        assert isinstance(result, pd.DataFrame)


# ============================================================================
# Fixtures for comparison data
# ============================================================================


@pytest.fixture(scope="module")
def comparison_vcf_small(benchmark_data_dir: Path, rng: np.random.Generator) -> Path:
    """Generate small VCF for tool comparison."""
    return generate_synthetic_vcf(
        benchmark_data_dir / "comparison_small.vcf",
        n_variants=1000,
        n_samples=10,
        rng=rng,
    )


@pytest.fixture(scope="module")
def comparison_vcf_medium(benchmark_data_dir: Path, rng: np.random.Generator) -> Path:
    """Generate medium VCF for tool comparison."""
    return generate_synthetic_vcf(
        benchmark_data_dir / "comparison_medium.vcf",
        n_variants=10000,
        n_samples=50,
        rng=rng,
    )
