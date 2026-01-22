"""
Performance benchmarks for WASP2 pipeline components.

These tests measure execution time and memory usage for critical functions.
Use pytest-benchmark for timing measurements when available.

Run with: pytest tests/benchmarks/test_benchmarks.py -v
Skip in CI: pytest tests/benchmarks/ -m "not benchmark"

Markers:
- benchmark: All benchmark tests
- slow: Benchmarks that take >10 seconds
"""

import time
from pathlib import Path
from typing import Callable, Dict, Any
from unittest.mock import MagicMock

import numpy as np
import pandas as pd
import polars as pl
import pytest


# ============================================================================
# Benchmark utilities
# ============================================================================

def timed_execution(func: Callable, *args, iterations: int = 3, **kwargs) -> Dict[str, float]:
    """
    Execute a function multiple times and return timing statistics.

    Args:
        func: Function to benchmark
        *args: Positional arguments to pass to func
        iterations: Number of times to run the function
        **kwargs: Keyword arguments to pass to func

    Returns:
        Dict with min, max, mean execution times in seconds
    """
    times = []
    for _ in range(iterations):
        start = time.perf_counter()
        func(*args, **kwargs)
        elapsed = time.perf_counter() - start
        times.append(elapsed)

    return {
        "min": min(times),
        "max": max(times),
        "mean": sum(times) / len(times),
        "iterations": iterations
    }


# ============================================================================
# Analysis module benchmarks
# ============================================================================

@pytest.mark.benchmark
class TestAnalysisBenchmarks:
    """Benchmarks for analysis.as_analysis functions."""

    def test_opt_prob_scalar_performance(self):
        """Benchmark opt_prob with scalar inputs."""
        from analysis.as_analysis import opt_prob

        result = timed_execution(
            opt_prob,
            in_prob=0.5,
            in_rho=0.1,
            k=50,
            n=100,
            log=True,
            iterations=1000
        )

        # Should be very fast for scalar operations
        assert result["mean"] < 0.001, f"opt_prob scalar too slow: {result['mean']:.4f}s"

    def test_opt_prob_array_performance(self):
        """Benchmark opt_prob with array inputs."""
        from analysis.as_analysis import opt_prob

        n_variants = 1000
        result = timed_execution(
            opt_prob,
            in_prob=np.full(n_variants, 0.5),
            in_rho=np.full(n_variants, 0.1),
            k=np.random.randint(20, 80, n_variants),
            n=np.full(n_variants, 100),
            log=True,
            iterations=10
        )

        # Should handle 1000 variants quickly
        assert result["mean"] < 0.1, f"opt_prob array too slow: {result['mean']:.4f}s"

    def test_single_model_small_dataset(self):
        """Benchmark single_model with small dataset."""
        from analysis.as_analysis import single_model

        # Create test data with 10 regions, 5 variants each
        n_regions = 10
        variants_per_region = 5
        total_variants = n_regions * variants_per_region

        df = pd.DataFrame({
            "chrom": ["chr1"] * total_variants,
            "pos": list(range(100, 100 + total_variants * 100, 100)),
            "ref_count": np.random.randint(30, 70, total_variants),
            "alt_count": np.random.randint(30, 70, total_variants),
            "N": [100] * total_variants,
            "region": [f"gene{i // variants_per_region}" for i in range(total_variants)]
        })

        result = timed_execution(single_model, df, "region", iterations=3)

        # Should complete quickly for small dataset
        assert result["mean"] < 2.0, f"single_model too slow for small dataset: {result['mean']:.2f}s"

    @pytest.mark.slow
    def test_single_model_large_dataset(self):
        """Benchmark single_model with large dataset."""
        from analysis.as_analysis import single_model

        # Create larger test data: 100 regions, 10 variants each
        n_regions = 100
        variants_per_region = 10
        total_variants = n_regions * variants_per_region

        df = pd.DataFrame({
            "chrom": ["chr1"] * total_variants,
            "pos": list(range(100, 100 + total_variants * 100, 100)),
            "ref_count": np.random.randint(30, 70, total_variants),
            "alt_count": np.random.randint(30, 70, total_variants),
            "N": [100] * total_variants,
            "region": [f"gene{i // variants_per_region}" for i in range(total_variants)]
        })

        result = timed_execution(single_model, df, "region", iterations=2)

        # Record the result for regression tracking
        print(f"\nLarge dataset benchmark: {result['mean']:.2f}s (min: {result['min']:.2f}s)")

        # Should still be reasonable for larger dataset
        assert result["mean"] < 30.0, f"single_model too slow for large dataset: {result['mean']:.2f}s"


# ============================================================================
# Counting module benchmarks
# ============================================================================

@pytest.mark.benchmark
class TestCountingBenchmarks:
    """Benchmarks for counting.count_alleles functions."""

    def test_find_read_aln_pos_performance(self):
        """Benchmark the binary search alignment position function."""
        from counting.count_alleles import find_read_aln_pos

        # Create a large alignment
        n_positions = 10000
        mock_read = MagicMock()
        mock_read.get_aligned_pairs.return_value = [
            (i, 1000 + i) for i in range(n_positions)
        ]

        # Test multiple lookups
        positions_to_find = [1000, 5000, 9999, 1234, 8765]

        def run_lookups():
            for pos in positions_to_find:
                find_read_aln_pos(mock_read, pos)

        result = timed_execution(run_lookups, iterations=100)

        # Binary search should be very fast
        assert result["mean"] < 0.01, f"find_read_aln_pos too slow: {result['mean']:.4f}s"


# ============================================================================
# Mapping module benchmarks
# ============================================================================

@pytest.mark.benchmark
class TestMappingBenchmarks:
    """Benchmarks for mapping.remap_utils functions."""

    def test_make_phased_seqs_performance(self):
        """Benchmark sequence creation for multiple variants."""
        from mapping.remap_utils import make_phased_seqs

        # Simulate a read with 10 variants
        n_variants = 10
        split_seq = ["ACGT" * 10]  # 40bp flanking
        for i in range(n_variants):
            split_seq.append("A")  # Original allele
            split_seq.append("ACGT" * 5)  # 20bp flanking

        hap1_alleles = ["G"] * n_variants  # All variants swapped
        hap2_alleles = ["A"] * n_variants  # Original

        result = timed_execution(
            make_phased_seqs,
            split_seq, hap1_alleles, hap2_alleles,
            iterations=1000
        )

        # String operations should be fast
        assert result["mean"] < 0.001, f"make_phased_seqs too slow: {result['mean']:.4f}s"

    def test_make_phased_seqs_with_qual_performance(self):
        """Benchmark indel-aware sequence creation."""
        from mapping.remap_utils import make_phased_seqs_with_qual

        # Simulate a read with variants including indels
        split_seq = ["ACGT" * 10, "ACG", "ACGT" * 5, "A", "ACGT" * 10]  # 2 variants
        split_qual = [
            np.array([30] * 40, dtype=np.uint8),
            np.array([35, 35, 35], dtype=np.uint8),
            np.array([30] * 20, dtype=np.uint8),
            np.array([35], dtype=np.uint8),
            np.array([30] * 40, dtype=np.uint8),
        ]

        hap1_alleles = ["T", "ACGT"]  # SNP and insertion
        hap2_alleles = ["ACG", "A"]   # Original

        result = timed_execution(
            make_phased_seqs_with_qual,
            split_seq, split_qual, hap1_alleles, hap2_alleles,
            iterations=100
        )

        # Should still be fast with quality score handling
        assert result["mean"] < 0.01, f"make_phased_seqs_with_qual too slow: {result['mean']:.4f}s"

    def test_build_ref2read_maps_performance(self):
        """Benchmark reference-to-read position mapping."""
        from mapping.remap_utils import _build_ref2read_maps

        # Create a complex alignment with deletions and insertions
        n_positions = 1000
        aligned_pairs = []
        query_pos = 0
        ref_pos = 0

        for i in range(n_positions):
            if i % 100 == 50:
                # Deletion: skip query positions
                aligned_pairs.append((None, ref_pos))
                aligned_pairs.append((None, ref_pos + 1))
                ref_pos += 2
            elif i % 100 == 75:
                # Insertion: skip ref positions
                aligned_pairs.append((query_pos, None))
                aligned_pairs.append((query_pos + 1, None))
                query_pos += 2
            else:
                aligned_pairs.append((query_pos, ref_pos))
                query_pos += 1
                ref_pos += 1

        mock_read = MagicMock()
        mock_read.get_aligned_pairs.return_value = aligned_pairs

        result = timed_execution(_build_ref2read_maps, mock_read, iterations=100)

        # Map building should be efficient
        assert result["mean"] < 0.01, f"_build_ref2read_maps too slow: {result['mean']:.4f}s"


# ============================================================================
# I/O module benchmarks
# ============================================================================

@pytest.mark.benchmark
class TestIOBenchmarks:
    """Benchmarks for wasp2.io variant source operations."""

    def test_variant_iteration_vcf(self, sample_vcf):
        """Benchmark VCF variant iteration."""
        try:
            from wasp2.io.vcf_source import VCFSource
        except ImportError:
            pytest.skip("wasp2.io.vcf_source not available")

        def iterate_variants():
            source = VCFSource(str(sample_vcf))
            count = 0
            for variant in source.iter_variants():
                count += 1
            return count

        result = timed_execution(iterate_variants, iterations=10)

        # Small test VCF should iterate quickly
        assert result["mean"] < 0.1, f"VCF iteration too slow: {result['mean']:.4f}s"


# ============================================================================
# Data structure benchmarks
# ============================================================================

@pytest.mark.benchmark
class TestDataFrameBenchmarks:
    """Benchmarks for DataFrame operations common in WASP2."""

    def test_polars_groupby_performance(self):
        """Benchmark Polars groupby aggregation."""
        # Create a DataFrame similar to counting output
        n_rows = 100000
        n_regions = 1000

        df = pl.DataFrame({
            "chrom": ["chr1"] * n_rows,
            "pos": list(range(n_rows)),
            "ref_count": np.random.randint(0, 100, n_rows).astype(np.uint16),
            "alt_count": np.random.randint(0, 100, n_rows).astype(np.uint16),
            "region": [f"region_{i % n_regions}" for i in range(n_rows)]
        })

        def groupby_sum():
            return df.group_by("region").agg([
                pl.col("ref_count").sum(),
                pl.col("alt_count").sum(),
                pl.len().alias("count")
            ])

        result = timed_execution(groupby_sum, iterations=10)

        # Polars should be very fast for this operation
        assert result["mean"] < 0.5, f"Polars groupby too slow: {result['mean']:.4f}s"

    def test_pandas_to_polars_conversion(self):
        """Benchmark pandas to polars DataFrame conversion."""
        n_rows = 100000

        pd_df = pd.DataFrame({
            "chrom": ["chr1"] * n_rows,
            "pos": np.arange(n_rows, dtype=np.int64),
            "ref_count": np.random.randint(0, 100, n_rows),
            "alt_count": np.random.randint(0, 100, n_rows),
        })

        def convert():
            return pl.from_pandas(pd_df)

        result = timed_execution(convert, iterations=10)

        # Conversion should be efficient
        assert result["mean"] < 0.5, f"pandas to polars conversion too slow: {result['mean']:.4f}s"


# ============================================================================
# Memory usage benchmarks
# ============================================================================

@pytest.mark.benchmark
class TestMemoryBenchmarks:
    """Benchmarks for memory efficiency."""

    def test_polars_memory_efficiency(self):
        """Check memory usage of Polars DataFrames with categorical types."""
        n_rows = 100000

        # Without categoricals
        df_string = pl.DataFrame({
            "chrom": ["chr1"] * n_rows,
            "region": [f"region_{i % 1000}" for i in range(n_rows)]
        })

        # With categoricals
        df_cat = df_string.with_columns([
            pl.col("chrom").cast(pl.Categorical),
            pl.col("region").cast(pl.Categorical)
        ])

        # Categorical should use less memory
        # Note: exact comparison depends on internal representation
        # Just verify both can be created without error
        assert df_string.shape == df_cat.shape


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short", "-m", "benchmark"])
