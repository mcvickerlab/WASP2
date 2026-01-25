"""
Sample size scaling benchmarks for WASP2.

Tests how WASP2 performance scales with the number of samples,
which is critical for large cohort studies.
"""

import numpy as np
import pandas as pd
import pytest

from analysis.as_analysis import get_imbalance, single_model

from .conftest import generate_allele_count_data


class TestSampleScaling:
    """Benchmark how performance scales with sample count."""

    @pytest.mark.benchmark(group="sample_scaling")
    @pytest.mark.parametrize("n_samples", [1, 10, 50, 100, 500])
    def test_get_imbalance_sample_scaling(self, benchmark, rng, n_samples: int):
        """Benchmark get_imbalance scaling with sample count."""
        n_variants = 5000
        n_regions = 500
        df = generate_allele_count_data(n_variants, n_samples, n_regions, rng=rng)

        result = benchmark(
            get_imbalance,
            df,
            min_count=10,
            method="single",
            region_col="region",
        )

        benchmark.extra_info["n_samples"] = n_samples
        benchmark.extra_info["n_variants"] = n_variants
        benchmark.extra_info["n_regions"] = n_regions
        benchmark.extra_info["total_coverage"] = int(df["ref_count"].sum() + df["alt_count"].sum())

        assert isinstance(result, pd.DataFrame)
        assert "pval" in result.columns

    @pytest.mark.benchmark(group="sample_scaling")
    @pytest.mark.parametrize("n_samples", [1, 10, 50, 100])
    def test_single_model_sample_scaling(self, benchmark, rng, n_samples: int):
        """Benchmark single_model scaling with sample count."""
        n_variants = 5000
        n_regions = 500
        df = generate_allele_count_data(n_variants, n_samples, n_regions, rng=rng)
        df["N"] = df["ref_count"] + df["alt_count"]

        result = benchmark(single_model, df, "region", False)

        benchmark.extra_info["n_samples"] = n_samples
        benchmark.extra_info["n_variants"] = n_variants
        benchmark.extra_info["mean_coverage"] = float(df["N"].mean())

        assert isinstance(result, pd.DataFrame)


class TestSampleVariantMatrix:
    """Combined sample × variant scaling matrix."""

    @pytest.mark.benchmark(group="sample_variant_matrix")
    @pytest.mark.parametrize(
        "n_samples,n_variants",
        [
            (1, 1000),
            (1, 10000),
            (10, 1000),
            (10, 10000),
            (50, 1000),
            (50, 10000),
            (100, 1000),
            (100, 10000),
        ],
    )
    def test_scaling_matrix(self, benchmark, rng, n_samples: int, n_variants: int):
        """Benchmark across sample × variant combinations."""
        n_regions = max(100, n_variants // 10)
        df = generate_allele_count_data(n_variants, n_samples, n_regions, rng=rng)

        result = benchmark(
            get_imbalance,
            df,
            min_count=10,
            method="single",
            region_col="region",
        )

        benchmark.extra_info["n_samples"] = n_samples
        benchmark.extra_info["n_variants"] = n_variants
        benchmark.extra_info["n_regions"] = n_regions
        benchmark.extra_info["matrix_cell"] = f"{n_samples}x{n_variants}"

        assert isinstance(result, pd.DataFrame)


class TestSampleMemoryScaling:
    """Memory scaling with sample count."""

    @pytest.mark.benchmark(group="sample_memory_scaling")
    @pytest.mark.parametrize("n_samples", [1, 10, 50, 100])
    def test_memory_sample_scaling(self, benchmark, rng, n_samples: int, memory_benchmark):
        """Measure memory scaling with sample count."""
        n_variants = 10000
        n_regions = 1000
        df = generate_allele_count_data(n_variants, n_samples, n_regions, rng=rng)

        def run_analysis():
            return get_imbalance(df, min_count=10, method="single", region_col="region")

        result, peak_memory = memory_benchmark.measure(run_analysis)

        benchmark.extra_info["n_samples"] = n_samples
        benchmark.extra_info["n_variants"] = n_variants
        benchmark.extra_info["peak_memory_mb"] = peak_memory
        benchmark.extra_info["data_size_mb"] = df.memory_usage(deep=True).sum() / (1024 * 1024)

        # Also run timed benchmark
        result = benchmark(run_analysis)
        assert isinstance(result, pd.DataFrame)


class TestCohortSimulation:
    """Simulate realistic cohort study scaling."""

    @pytest.mark.benchmark(group="cohort_simulation")
    @pytest.mark.slow
    @pytest.mark.parametrize(
        "cohort_size,variants_per_sample",
        [
            ("small", 10),  # Small study: 10 samples, 10K variants
            ("medium", 50),  # Medium study: 50 samples, 50K variants
            ("large", 100),  # Large study: 100 samples, 100K variants
        ],
    )
    def test_cohort_study_simulation(
        self,
        benchmark,
        rng,
        cohort_size: str,
        variants_per_sample: int,
    ):
        """Simulate realistic cohort study workloads."""
        cohort_config = {
            "small": {"n_samples": 10, "n_variants": 10000},
            "medium": {"n_samples": 50, "n_variants": 50000},
            "large": {"n_samples": 100, "n_variants": 100000},
        }

        config = cohort_config[cohort_size]
        n_samples = config["n_samples"]
        n_variants = config["n_variants"]
        n_regions = n_variants // 10

        df = generate_allele_count_data(n_variants, n_samples, n_regions, rng=rng)

        result = benchmark(
            get_imbalance,
            df,
            min_count=10,
            method="single",
            region_col="region",
        )

        benchmark.extra_info["cohort_size"] = cohort_size
        benchmark.extra_info["n_samples"] = n_samples
        benchmark.extra_info["n_variants"] = n_variants
        benchmark.extra_info["n_regions"] = n_regions
        benchmark.extra_info["total_coverage"] = int(df["ref_count"].sum() + df["alt_count"].sum())

        assert isinstance(result, pd.DataFrame)


class TestHighThroughputScaling:
    """High-throughput sequencing depth scaling."""

    @pytest.mark.benchmark(group="coverage_scaling")
    @pytest.mark.parametrize("coverage_multiplier", [1, 5, 10, 20])
    def test_coverage_depth_scaling(self, benchmark, rng, coverage_multiplier: int):
        """Test performance with varying sequencing depth."""
        n_variants = 5000
        n_samples = 10
        n_regions = 500

        # Generate base data
        df = generate_allele_count_data(n_variants, n_samples, n_regions, rng=rng)

        # Scale coverage
        df["ref_count"] = (df["ref_count"] * coverage_multiplier).astype(np.uint32)
        df["alt_count"] = (df["alt_count"] * coverage_multiplier).astype(np.uint32)

        result = benchmark(
            get_imbalance,
            df,
            min_count=10,
            method="single",
            region_col="region",
        )

        benchmark.extra_info["coverage_multiplier"] = coverage_multiplier
        benchmark.extra_info["mean_coverage"] = float(
            (df["ref_count"] + df["alt_count"]).mean()
        )
        benchmark.extra_info["n_variants"] = n_variants

        assert isinstance(result, pd.DataFrame)
