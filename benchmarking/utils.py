"""
Benchmark utilities for WASP2 performance validation.

Provides:
- Timer context manager with statistical analysis
- Tool availability checking
- Result formatting and reporting
- Data generation helpers
"""

import gc
import json
import shutil
import statistics
import subprocess
import time
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np


@dataclass
class BenchmarkResult:
    """Container for benchmark timing results with statistics."""

    name: str
    iterations: list[float] = field(default_factory=list)
    tool: str = ""
    parameters: dict[str, Any] = field(default_factory=dict)
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def mean(self) -> float:
        return statistics.mean(self.iterations) if self.iterations else 0.0

    @property
    def std(self) -> float:
        return statistics.stdev(self.iterations) if len(self.iterations) > 1 else 0.0

    @property
    def min(self) -> float:
        return min(self.iterations) if self.iterations else 0.0

    @property
    def max(self) -> float:
        return max(self.iterations) if self.iterations else 0.0

    @property
    def median(self) -> float:
        return statistics.median(self.iterations) if self.iterations else 0.0

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "tool": self.tool,
            "mean": self.mean,
            "std": self.std,
            "min": self.min,
            "max": self.max,
            "median": self.median,
            "iterations": len(self.iterations),
            "raw_times": self.iterations,
            "parameters": self.parameters,
            "metadata": self.metadata,
        }


class BenchmarkTimer:
    """
    Context manager for timing benchmark operations.

    Usage:
        timer = BenchmarkTimer("my_operation", warmup=2, iterations=5)
        for t in timer:
            with t:
                run_operation()
        print(timer.result)
    """

    def __init__(
        self,
        name: str,
        warmup: int = 2,
        iterations: int = 5,
        gc_collect: bool = True,
    ):
        self.name = name
        self.warmup = warmup
        self.iterations = iterations
        self.gc_collect = gc_collect
        self.result = BenchmarkResult(name=name)
        self._current_iteration = 0
        self._is_warmup = True
        self._start_time: float = 0.0

    def __iter__(self):
        self._current_iteration = 0
        self._is_warmup = True
        total = self.warmup + self.iterations
        for i in range(total):
            self._is_warmup = i < self.warmup
            self._current_iteration = i
            if self.gc_collect:
                gc.collect()
            yield self

    def __enter__(self):
        self._start_time = time.perf_counter()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        elapsed = time.perf_counter() - self._start_time
        if not self._is_warmup:
            self.result.iterations.append(elapsed)
        return False


def check_tool(tool_name: str) -> bool:
    """Check if an external tool is available in PATH."""
    return shutil.which(tool_name) is not None


def get_tool_version(tool_name: str, version_flag: str = "--version") -> str | None:
    """Get version string for an external tool."""
    try:
        result = subprocess.run(
            [tool_name, version_flag],
            capture_output=True,
            text=True,
            timeout=10,
        )
        return result.stdout.strip() or result.stderr.strip()
    except (subprocess.SubprocessError, FileNotFoundError):
        return None


def format_time(seconds: float) -> str:
    """Format time in human-readable units."""
    if seconds < 0.001:
        return f"{seconds * 1_000_000:.2f} μs"
    elif seconds < 1:
        return f"{seconds * 1000:.2f} ms"
    elif seconds < 60:
        return f"{seconds:.3f} s"
    else:
        minutes = int(seconds // 60)
        secs = seconds % 60
        return f"{minutes}m {secs:.1f}s"


def print_comparison_table(
    results: list[BenchmarkResult],
    baseline_tool: str = "WASP2",
) -> None:
    """Print a formatted comparison table of benchmark results."""
    print("\n" + "=" * 70)
    print("BENCHMARK COMPARISON")
    print("=" * 70)

    baseline = next((r for r in results if r.tool == baseline_tool), None)

    header = f"{'Tool':<15} {'Mean':>12} {'Std':>10} {'Speedup':>10}"
    print(header)
    print("-" * 50)

    for r in sorted(results, key=lambda x: x.mean):
        speedup = ""
        if baseline and r.tool != baseline_tool and r.mean > 0:
            speedup = (
                f"{baseline.mean / r.mean:.2f}x"
                if baseline.mean < r.mean
                else f"{r.mean / baseline.mean:.2f}x slower"
            )
        print(f"{r.tool:<15} {format_time(r.mean):>12} {format_time(r.std):>10} {speedup:>10}")

    print("=" * 70)


def save_results(
    results: list[BenchmarkResult],
    output_path: Path,
    include_raw: bool = False,
) -> None:
    """Save benchmark results to JSON file."""
    benchmarks = [r.to_dict() for r in results]
    if not include_raw:
        for b in benchmarks:
            b.pop("raw_times", None)
    output_data: dict[str, Any] = {
        "timestamp": datetime.now().isoformat(),
        "benchmarks": benchmarks,
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(output_data, f, indent=2)
    print(f"Results saved to {output_path}")


def generate_synthetic_counts(
    n_variants: int,
    n_regions: int,
    seed: int = 42,
):
    """Generate synthetic allele count data for benchmarking."""
    import pandas as pd

    rng = np.random.default_rng(seed)
    chroms = rng.choice([f"chr{i}" for i in range(1, 23)], size=n_variants)
    positions = rng.integers(1, 250_000_000, size=n_variants)
    bases = ["A", "C", "G", "T"]
    refs = rng.choice(bases, size=n_variants)
    alts = np.array([rng.choice([b for b in bases if b != r]) for r in refs])

    total_counts = rng.exponential(scale=30, size=n_variants).astype(int) + 10
    ratios = rng.beta(10, 10, size=n_variants)
    ref_counts = (total_counts * ratios).astype(int)
    alt_counts = total_counts - ref_counts

    region_names = [f"region_{i:06d}" for i in range(n_regions)]
    regions = rng.choice(region_names, size=n_variants)

    return pd.DataFrame(
        {
            "chrom": pd.Categorical(chroms),
            "pos": positions.astype(np.uint32),
            "ref": pd.Categorical(refs),
            "alt": pd.Categorical(alts),
            "ref_count": ref_counts.astype(np.uint32),
            "alt_count": alt_counts.astype(np.uint32),
            "other_count": np.zeros(n_variants, dtype=np.uint16),
            "region": regions,
        }
    )
