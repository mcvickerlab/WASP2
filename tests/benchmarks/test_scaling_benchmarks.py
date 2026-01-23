"""
Scaling benchmarks for WASP2.

Tests how performance scales with variant count, region count, and sample count.
"""

import numpy as np
import pandas as pd
import pytest

from analysis.as_analysis import get_imbalance, single_model


def generate_scaled_data(
    n_variants: int,
    n_samples: int,
    n_regions: int,
    rng: np.random.Generator,
    include_phasing: bool = False,
) -> pd.DataFrame:
    """Generate synthetic allele count data at specified scale."""
    chroms = rng.choice([f"chr{i}" for i in range(1, 23)], size=n_variants)
    positions = rng.integers(1, 250_000_000, size=n_variants)
    bases = ["A", "C", "G", "T"]
    refs = rng.choice(bases, size=n_variants)
    alts = np.array([rng.choice([b for b in bases if b != r]) for r in refs])

    mean_coverage = 30 * n_samples / 10
    total_counts = rng.exponential(scale=mean_coverage, size=n_variants).astype(int) + 10
    ratios = rng.beta(10, 10, size=n_variants)
    ref_counts = (total_counts * ratios).astype(int)
    alt_counts = total_counts - ref_counts

    region_names = [f"region_{i:06d}" for i in range(n_regions)]
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


class TestVariantScaling:
    """Benchmark how performance scales with variant count."""

    @pytest.mark.benchmark(group="variant_scaling")
    @pytest.mark.parametrize("n_variants", [100, 1000, 10000])
    def test_single_model_variant_scaling(self, benchmark, rng, n_variants: int):
        """Benchmark single_model scaling with variant count."""
        df = generate_scaled_data(n_variants, n_samples=10, n_regions=100, rng=rng)
        df["N"] = df["ref_count"] + df["alt_count"]

        result = benchmark(single_model, df, "region", False)
        benchmark.extra_info["n_variants"] = n_variants
        assert isinstance(result, pd.DataFrame)

    @pytest.mark.benchmark(group="variant_scaling")
    @pytest.mark.parametrize("n_variants", [100, 1000, 10000])
    def test_get_imbalance_variant_scaling(self, benchmark, rng, n_variants: int):
        """Benchmark get_imbalance scaling with variant count."""
        df = generate_scaled_data(n_variants, n_samples=10, n_regions=100, rng=rng)

        result = benchmark(get_imbalance, df, min_count=10, method="single", region_col="region")
        benchmark.extra_info["n_variants"] = n_variants
        assert isinstance(result, pd.DataFrame)


class TestRegionScaling:
    """Benchmark how performance scales with region count."""

    @pytest.mark.benchmark(group="region_scaling")
    @pytest.mark.parametrize("n_regions", [10, 100, 1000])
    def test_single_model_region_scaling(self, benchmark, rng, n_regions: int):
        """Benchmark single_model scaling with region count."""
        n_variants = n_regions * 10
        df = generate_scaled_data(n_variants, n_samples=10, n_regions=n_regions, rng=rng)
        df["N"] = df["ref_count"] + df["alt_count"]

        result = benchmark(single_model, df, "region", False)
        benchmark.extra_info["n_regions"] = n_regions
        assert isinstance(result, pd.DataFrame)


class TestMethodComparison:
    """Compare single vs linear model performance at various scales."""

    @pytest.mark.benchmark(group="method_comparison")
    @pytest.mark.parametrize(
        "method,n_variants",
        [
            ("single", 100),
            ("single", 1000),
            ("linear", 100),
            ("linear", 1000),
        ],
    )
    def test_method_scaling_comparison(self, benchmark, rng, method: str, n_variants: int):
        """Compare single vs linear model at different scales."""
        n_regions = max(10, n_variants // 10)
        df = generate_scaled_data(n_variants, n_samples=10, n_regions=n_regions, rng=rng)

        result = benchmark(get_imbalance, df, min_count=10, method=method, region_col="region")
        benchmark.extra_info["method"] = method
        benchmark.extra_info["n_variants"] = n_variants
        assert isinstance(result, pd.DataFrame)


class TestMemoryScaling:
    """Benchmark memory usage at different scales."""

    @pytest.mark.benchmark(group="memory_scaling")
    @pytest.mark.parametrize("n_variants", [1000, 10000])
    def test_memory_scaling(self, benchmark, rng, n_variants: int, memory_benchmark):
        """Measure memory scaling with variant count."""
        n_regions = max(100, n_variants // 10)
        df = generate_scaled_data(n_variants, n_samples=10, n_regions=n_regions, rng=rng)

        def run_analysis():
            return get_imbalance(df, min_count=10, method="single", region_col="region")

        result, peak_memory = memory_benchmark.measure(run_analysis)
        benchmark.extra_info["n_variants"] = n_variants
        benchmark.extra_info["peak_memory_mb"] = peak_memory

        result = benchmark(run_analysis)
        assert isinstance(result, pd.DataFrame)
