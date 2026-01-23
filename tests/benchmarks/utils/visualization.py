"""
Benchmark visualization utilities for generating publication-quality figures.
"""

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

try:
    import matplotlib.pyplot as plt

    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

try:
    import seaborn as sns

    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False


@dataclass
class BenchmarkResult:
    """Container for a single benchmark result."""

    name: str
    group: str
    mean: float
    stddev: float
    min: float
    max: float
    rounds: int
    extra_info: dict[str, Any]

    @classmethod
    def from_dict(cls, data: dict) -> "BenchmarkResult":
        """Create from pytest-benchmark JSON format."""
        stats = data.get("stats", {})
        return cls(
            name=data.get("name", ""),
            group=data.get("group", ""),
            mean=stats.get("mean", 0),
            stddev=stats.get("stddev", 0),
            min=stats.get("min", 0),
            max=stats.get("max", 0),
            rounds=stats.get("rounds", 0),
            extra_info=data.get("extra_info", {}),
        )


def load_benchmark_results(json_path: Path) -> list[BenchmarkResult]:
    """Load benchmark results from pytest-benchmark JSON file."""
    with open(json_path) as f:
        data = json.load(f)

    return [BenchmarkResult.from_dict(b) for b in data.get("benchmarks", [])]


def filter_by_group(results: list[BenchmarkResult], group: str) -> list[BenchmarkResult]:
    """Filter results by benchmark group."""
    return [r for r in results if r.group == group]


def setup_publication_style():
    """Configure matplotlib for publication-quality figures."""
    if not HAS_MATPLOTLIB:
        raise ImportError("matplotlib required for visualization")

    if HAS_SEABORN:
        sns.set_theme(style="whitegrid", context="paper", font_scale=1.2)
    else:
        plt.style.use("seaborn-v0_8-whitegrid")

    plt.rcParams.update(
        {
            "figure.figsize": (8, 6),
            "figure.dpi": 150,
            "savefig.dpi": 300,
            "font.size": 12,
            "axes.labelsize": 14,
            "axes.titlesize": 16,
        }
    )


def plot_scaling(
    results: list[BenchmarkResult],
    x_param: str,
    title: str = "Performance Scaling",
    xlabel: str = "Problem Size",
    ylabel: str = "Time (seconds)",
    log_scale: bool = True,
    output_path: Path | None = None,
):
    """Plot performance scaling with problem size."""
    if not HAS_MATPLOTLIB:
        raise ImportError("matplotlib required for visualization")

    setup_publication_style()

    x_values, y_values, y_errors = [], [], []
    for r in results:
        if x_param in r.extra_info:
            x_values.append(r.extra_info[x_param])
            y_values.append(r.mean)
            y_errors.append(r.stddev)

    if not x_values:
        raise ValueError(f"No results found with {x_param} in extra_info")

    sorted_indices = np.argsort(x_values)
    x_values = np.array(x_values)[sorted_indices]
    y_values = np.array(y_values)[sorted_indices]
    y_errors = np.array(y_errors)[sorted_indices]

    fig, ax = plt.subplots()
    ax.errorbar(x_values, y_values, yerr=y_errors, marker="o", capsize=5, color="steelblue")

    if log_scale:
        ax.set_xscale("log")
        ax.set_yscale("log")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, bbox_inches="tight")
        print(f"Saved figure to {output_path}")

    return fig, ax


def generate_all_figures(
    results_path: Path,
    output_dir: Path,
    formats: list[str] | None = None,
) -> None:
    """Generate all standard benchmark figures."""
    if formats is None:
        formats = ["png", "pdf"]

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    results = load_benchmark_results(results_path)

    variant_scaling = filter_by_group(results, "variant_scaling")
    if variant_scaling:
        for fmt in formats:
            plot_scaling(
                variant_scaling,
                x_param="n_variants",
                title="Performance vs. Variant Count",
                xlabel="Number of Variants",
                output_path=output_dir / f"variant_scaling.{fmt}",
            )
            plt.close()

    print(f"Generated figures in {output_dir}")
