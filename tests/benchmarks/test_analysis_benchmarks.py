"""
Performance benchmarks for WASP2 analysis module.

Tests statistical analysis functions including allelic imbalance calculation.
"""

import numpy as np
import pandas as pd
import pytest

from analysis.as_analysis import get_imbalance, linear_model, opt_linear, opt_prob, single_model


def generate_allele_count_df(
    n_variants: int,
    n_regions: int = 100,
    rng: np.random.Generator | None = None,
    include_phasing: bool = False,
) -> pd.DataFrame:
    """Generate synthetic allele count data for benchmarking."""
    if rng is None:
        rng = np.random.default_rng(42)

    chroms = rng.choice([f"chr{i}" for i in range(1, 23)], size=n_variants)
    positions = rng.integers(1, 250_000_000, size=n_variants)
    bases = ["A", "C", "G", "T"]
    refs = rng.choice(bases, size=n_variants)
    alts = np.array([rng.choice([b for b in bases if b != r]) for r in refs])

    total_counts = rng.exponential(scale=30, size=n_variants).astype(int) + 10
    ratios = rng.beta(10, 10, size=n_variants)
    ref_counts = (total_counts * ratios).astype(int)
    alt_counts = total_counts - ref_counts

    region_names = [f"region_{i:05d}" for i in range(n_regions)]
    regions = rng.choice(region_names, size=n_variants)

    df = pd.DataFrame(
        {
            "chrom": pd.Categorical(chroms),
            "pos": positions.astype(np.uint32),
            "ref": pd.Categorical(refs),
            "alt": pd.Categorical(alts),
            "ref_count": ref_counts.astype(np.uint16),
            "alt_count": alt_counts.astype(np.uint16),
            "other_count": np.zeros(n_variants, dtype=np.uint16),
            "region": regions,
        }
    )

    if include_phasing:
        df["GT"] = rng.choice(["0|1", "1|0"], size=n_variants)

    return df


class TestOptimizationBenchmarks:
    """Benchmark tests for optimization functions."""

    @pytest.mark.benchmark(group="optimization")
    def test_opt_linear_small(self, benchmark, rng):
        """Benchmark opt_linear with small dataset."""
        n = 100
        ref_counts = rng.integers(5, 50, size=n)
        n_array = rng.integers(20, 100, size=n)
        disp_params = np.array([0.0, 0.0])

        result = benchmark(opt_linear, disp_params, ref_counts, n_array)
        assert isinstance(result, float)

    @pytest.mark.benchmark(group="optimization")
    def test_opt_linear_large(self, benchmark, rng):
        """Benchmark opt_linear with large dataset."""
        n = 10000
        ref_counts = rng.integers(5, 50, size=n)
        n_array = rng.integers(20, 100, size=n)
        disp_params = np.array([0.0, 0.0])

        result = benchmark(opt_linear, disp_params, ref_counts, n_array)
        assert isinstance(result, float)

    @pytest.mark.benchmark(group="optimization")
    def test_opt_prob_large(self, benchmark, rng):
        """Benchmark opt_prob with large dataset."""
        n = 10000
        ref_counts = rng.integers(5, 50, size=n)
        n_array = rng.integers(20, 100, size=n)

        result = benchmark(opt_prob, 0.1, 0.05, ref_counts, n_array)
        assert isinstance(result, float | np.floating | np.ndarray)


class TestModelBenchmarks:
    """Benchmark tests for statistical models."""

    @pytest.mark.benchmark(group="models")
    def test_single_model_1k_variants(self, benchmark, rng):
        """Benchmark single_model with 1K variants."""
        df = generate_allele_count_df(1000, n_regions=50, rng=rng)
        df["N"] = df["ref_count"] + df["alt_count"]

        result = benchmark(single_model, df, "region", False)
        assert isinstance(result, pd.DataFrame)

    @pytest.mark.benchmark(group="models")
    def test_linear_model_1k_variants(self, benchmark, rng):
        """Benchmark linear_model with 1K variants."""
        df = generate_allele_count_df(1000, n_regions=50, rng=rng)
        df["N"] = df["ref_count"] + df["alt_count"]

        result = benchmark(linear_model, df, "region", False)
        assert isinstance(result, pd.DataFrame)


class TestGetImbalanceBenchmarks:
    """Benchmark tests for the main get_imbalance function."""

    @pytest.mark.benchmark(group="get_imbalance")
    def test_get_imbalance_single_1k(self, benchmark, rng):
        """Benchmark get_imbalance (single model) with 1K variants."""
        df = generate_allele_count_df(1000, n_regions=100, rng=rng)

        result = benchmark(
            get_imbalance,
            df,
            min_count=10,
            pseudocount=1,
            method="single",
            region_col="region",
        )
        assert isinstance(result, pd.DataFrame)
        assert "pval" in result.columns

    @pytest.mark.benchmark(group="get_imbalance")
    def test_get_imbalance_linear_1k(self, benchmark, rng):
        """Benchmark get_imbalance (linear model) with 1K variants."""
        df = generate_allele_count_df(1000, n_regions=100, rng=rng)

        result = benchmark(
            get_imbalance,
            df,
            min_count=10,
            pseudocount=1,
            method="linear",
            region_col="region",
        )
        assert isinstance(result, pd.DataFrame)

    @pytest.mark.benchmark(group="get_imbalance")
    @pytest.mark.slow
    def test_get_imbalance_single_10k(self, benchmark, rng):
        """Benchmark get_imbalance (single model) with 10K variants."""
        df = generate_allele_count_df(10000, n_regions=1000, rng=rng)

        result = benchmark(
            get_imbalance,
            df,
            min_count=10,
            pseudocount=1,
            method="single",
            region_col="region",
        )
        assert isinstance(result, pd.DataFrame)


class TestAnalysisMemoryBenchmarks:
    """Memory usage benchmarks for analysis functions."""

    @pytest.mark.benchmark(group="memory")
    def test_get_imbalance_memory_1k(self, benchmark, rng, memory_benchmark):
        """Measure memory usage for 1K variant analysis."""
        df = generate_allele_count_df(1000, n_regions=100, rng=rng)

        def run_analysis():
            return get_imbalance(df, min_count=10, method="single", region_col="region")

        result, peak_memory = memory_benchmark.measure(run_analysis)
        benchmark.extra_info["peak_memory_mb"] = peak_memory
        result = benchmark(run_analysis)
        assert isinstance(result, pd.DataFrame)
