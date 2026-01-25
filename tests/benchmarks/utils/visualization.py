"""
Benchmark visualization utilities for generating publication-quality figures.

Supports:
- Multi-tool comparison plots (bar charts, grouped comparisons)
- Scaling plots with error bars (log-scale capable)
- Heatmaps for sample × variant matrices
- Memory profiling visualizations
- Nature/Cell-style publication formatting
"""

import json
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

try:
    import matplotlib.pyplot as plt

    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    plt = None

try:
    import seaborn as sns

    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False
    sns = None


# ============================================================================
# Color palettes for publication figures
# ============================================================================

# Nature-style color palette
NATURE_COLORS = {
    "wasp2": "#E64B35",  # Red
    "wasp_v1": "#4DBBD5",  # Cyan
    "phaser": "#00A087",  # Teal
    "default": "#3C5488",  # Blue
    "secondary": "#F39B7F",  # Salmon
    "tertiary": "#8491B4",  # Gray-blue
}

# Colorblind-friendly palette
COLORBLIND_SAFE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#F0E442", "#56B4E9"]


# ============================================================================
# Data classes
# ============================================================================


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
    extra_info: dict[str, Any] = field(default_factory=dict)

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


@dataclass
class FigureConfig:
    """Configuration for figure generation."""

    figsize: tuple[float, float] = (8, 6)
    dpi: int = 150
    save_dpi: int = 300
    font_size: int = 12
    label_size: int = 14
    title_size: int = 16
    line_width: float = 2.0
    marker_size: int = 8
    colorblind_safe: bool = True


# ============================================================================
# Result loading and filtering
# ============================================================================


def load_benchmark_results(json_path: Path) -> list[BenchmarkResult]:
    """Load benchmark results from pytest-benchmark JSON file."""
    with open(json_path) as f:
        data = json.load(f)

    return [BenchmarkResult.from_dict(b) for b in data.get("benchmarks", [])]


def load_multiple_results(json_paths: list[Path]) -> list[BenchmarkResult]:
    """Load and combine results from multiple JSON files."""
    all_results = []
    for path in json_paths:
        all_results.extend(load_benchmark_results(path))
    return all_results


def filter_by_group(results: list[BenchmarkResult], group: str) -> list[BenchmarkResult]:
    """Filter results by benchmark group."""
    return [r for r in results if r.group == group]


def filter_by_extra_info(
    results: list[BenchmarkResult],
    key: str,
    value: Any,
) -> list[BenchmarkResult]:
    """Filter results by extra_info key-value pair."""
    return [r for r in results if r.extra_info.get(key) == value]


def group_by_extra_info(
    results: list[BenchmarkResult],
    key: str,
) -> dict[Any, list[BenchmarkResult]]:
    """Group results by an extra_info key."""
    grouped = defaultdict(list)
    for r in results:
        if key in r.extra_info:
            grouped[r.extra_info[key]].append(r)
    return dict(grouped)


# ============================================================================
# Style configuration
# ============================================================================


def setup_publication_style(config: FigureConfig | None = None):
    """Configure matplotlib for publication-quality figures (Nature/Cell style)."""
    if not HAS_MATPLOTLIB:
        raise ImportError("matplotlib required for visualization")

    if config is None:
        config = FigureConfig()

    if HAS_SEABORN:
        sns.set_theme(style="whitegrid", context="paper", font_scale=1.2)
    else:
        try:
            plt.style.use("seaborn-v0_8-whitegrid")
        except OSError:
            plt.style.use("ggplot")

    plt.rcParams.update(
        {
            "figure.figsize": config.figsize,
            "figure.dpi": config.dpi,
            "savefig.dpi": config.save_dpi,
            "font.size": config.font_size,
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "axes.labelsize": config.label_size,
            "axes.titlesize": config.title_size,
            "axes.linewidth": 1.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.major.width": 1.5,
            "ytick.major.width": 1.5,
            "legend.frameon": False,
            "legend.fontsize": config.font_size - 1,
            "lines.linewidth": config.line_width,
            "lines.markersize": config.marker_size,
        }
    )


def get_color_palette(n_colors: int, colorblind_safe: bool = True) -> list[str]:
    """Get a color palette for plotting."""
    if colorblind_safe:
        return COLORBLIND_SAFE[:n_colors]
    return list(NATURE_COLORS.values())[:n_colors]


def get_tool_color(tool_name: str) -> str:
    """Get consistent color for a tool."""
    tool_lower = tool_name.lower().replace(" ", "_")
    return NATURE_COLORS.get(tool_lower, NATURE_COLORS["default"])


# ============================================================================
# Basic scaling plots
# ============================================================================


def plot_scaling(
    results: list[BenchmarkResult],
    x_param: str,
    title: str = "Performance Scaling",
    xlabel: str = "Problem Size",
    ylabel: str = "Time (seconds)",
    log_scale: bool = True,
    output_path: Path | None = None,
    color: str | None = None,
    label: str | None = None,
    ax: Any | None = None,
) -> tuple[Any, Any]:
    """Plot performance scaling with problem size."""
    if not HAS_MATPLOTLIB:
        raise ImportError("matplotlib required for visualization")

    if ax is None:
        setup_publication_style()
        fig, ax = plt.subplots()
    else:
        fig = ax.figure

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

    plot_color = color or NATURE_COLORS["default"]
    ax.errorbar(
        x_values,
        y_values,
        yerr=y_errors,
        marker="o",
        capsize=5,
        color=plot_color,
        label=label,
    )

    if log_scale:
        ax.set_xscale("log")
        ax.set_yscale("log")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3, linestyle="--")

    if label:
        ax.legend()

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, bbox_inches="tight", facecolor="white")
        print(f"Saved figure to {output_path}")

    return fig, ax


# ============================================================================
# Tool comparison plots
# ============================================================================


def plot_tool_comparison(
    results: list[BenchmarkResult],
    x_param: str = "n_variants",
    title: str = "Tool Performance Comparison",
    xlabel: str = "Number of Variants",
    ylabel: str = "Time (seconds)",
    log_scale: bool = True,
    output_path: Path | None = None,
) -> tuple[Any, Any]:
    """Plot multi-tool comparison with different lines per tool."""
    if not HAS_MATPLOTLIB:
        raise ImportError("matplotlib required for visualization")

    setup_publication_style()
    fig, ax = plt.subplots()

    # Group results by tool
    tool_groups = group_by_extra_info(results, "tool")

    for tool_name, tool_results in tool_groups.items():
        x_values, y_values, y_errors = [], [], []
        for r in tool_results:
            if x_param in r.extra_info:
                x_values.append(r.extra_info[x_param])
                y_values.append(r.mean)
                y_errors.append(r.stddev)

        if not x_values:
            continue

        sorted_indices = np.argsort(x_values)
        x_values = np.array(x_values)[sorted_indices]
        y_values = np.array(y_values)[sorted_indices]
        y_errors = np.array(y_errors)[sorted_indices]

        color = get_tool_color(tool_name)
        ax.errorbar(
            x_values,
            y_values,
            yerr=y_errors,
            marker="o",
            capsize=5,
            color=color,
            label=tool_name,
        )

    if log_scale:
        ax.set_xscale("log")
        ax.set_yscale("log")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3, linestyle="--")
    ax.legend(loc="upper left")
    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, bbox_inches="tight", facecolor="white")
        print(f"Saved figure to {output_path}")

    return fig, ax


def plot_tool_comparison_bars(
    results: list[BenchmarkResult],
    group_param: str = "n_variants",
    title: str = "Tool Performance Comparison",
    ylabel: str = "Time (seconds)",
    log_scale: bool = True,
    output_path: Path | None = None,
) -> tuple[Any, Any]:
    """Plot grouped bar chart comparing tools at different scales."""
    if not HAS_MATPLOTLIB:
        raise ImportError("matplotlib required for visualization")

    setup_publication_style()
    fig, ax = plt.subplots(figsize=(10, 6))

    # Organize data: {scale: {tool: (mean, std)}}
    data = defaultdict(dict)
    for r in results:
        if group_param in r.extra_info and "tool" in r.extra_info:
            scale = r.extra_info[group_param]
            tool = r.extra_info["tool"]
            data[scale][tool] = (r.mean, r.stddev)

    scales = sorted(data.keys())
    tools = sorted({tool for scale_data in data.values() for tool in scale_data})

    x = np.arange(len(scales))
    width = 0.8 / len(tools)

    for i, tool in enumerate(tools):
        means = [data[s].get(tool, (0, 0))[0] for s in scales]
        stds = [data[s].get(tool, (0, 0))[1] for s in scales]
        offset = (i - len(tools) / 2 + 0.5) * width
        color = get_tool_color(tool)
        ax.bar(x + offset, means, width, yerr=stds, label=tool, color=color, capsize=3)

    ax.set_xlabel(f"Number of {group_param.replace('n_', '').title()}")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{s:,}" for s in scales])
    ax.legend()

    if log_scale:
        ax.set_yscale("log")

    ax.grid(True, alpha=0.3, axis="y", linestyle="--")
    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, bbox_inches="tight", facecolor="white")
        print(f"Saved figure to {output_path}")

    return fig, ax


# ============================================================================
# Heatmap plots
# ============================================================================


def plot_scaling_heatmap(
    results: list[BenchmarkResult],
    x_param: str = "n_samples",
    y_param: str = "n_variants",
    value_key: str = "mean",
    title: str = "Performance Scaling Matrix",
    xlabel: str = "Number of Samples",
    ylabel: str = "Number of Variants",
    cbar_label: str = "Time (seconds)",
    log_scale: bool = True,
    output_path: Path | None = None,
) -> tuple[Any, Any]:
    """Plot heatmap for 2D parameter scaling (e.g., sample × variant matrix)."""
    if not HAS_MATPLOTLIB:
        raise ImportError("matplotlib required for visualization")

    setup_publication_style()

    # Extract unique parameter values
    x_values = sorted({r.extra_info.get(x_param) for r in results if x_param in r.extra_info})
    y_values = sorted({r.extra_info.get(y_param) for r in results if y_param in r.extra_info})

    if not x_values or not y_values:
        raise ValueError(f"No results found with {x_param} and {y_param} in extra_info")

    # Create matrix
    matrix = np.full((len(y_values), len(x_values)), np.nan)
    for r in results:
        if x_param in r.extra_info and y_param in r.extra_info:
            x_idx = x_values.index(r.extra_info[x_param])
            y_idx = y_values.index(r.extra_info[y_param])
            if value_key == "mean":
                matrix[y_idx, x_idx] = r.mean
            elif value_key in r.extra_info:
                matrix[y_idx, x_idx] = r.extra_info[value_key]

    fig, ax = plt.subplots(figsize=(10, 8))

    # Apply log transform if requested
    plot_data = np.log10(matrix) if log_scale else matrix
    cbar_label_final = f"log10({cbar_label})" if log_scale else cbar_label

    # Create heatmap
    if HAS_SEABORN:
        heatmap = sns.heatmap(
            plot_data,
            ax=ax,
            cmap="viridis",
            annot=True,
            fmt=".2f",
            xticklabels=[f"{v:,}" for v in x_values],
            yticklabels=[f"{v:,}" for v in y_values],
            cbar_kws={"label": cbar_label_final},
        )
    else:
        im = ax.imshow(plot_data, cmap="viridis", aspect="auto")
        ax.set_xticks(range(len(x_values)))
        ax.set_xticklabels([f"{v:,}" for v in x_values])
        ax.set_yticks(range(len(y_values)))
        ax.set_yticklabels([f"{v:,}" for v in y_values])
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label(cbar_label_final)

        # Add annotations
        for i in range(len(y_values)):
            for j in range(len(x_values)):
                if not np.isnan(plot_data[i, j]):
                    ax.text(j, i, f"{plot_data[i, j]:.2f}", ha="center", va="center", color="white")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, bbox_inches="tight", facecolor="white")
        print(f"Saved figure to {output_path}")

    return fig, ax


# ============================================================================
# Memory profiling plots
# ============================================================================


def plot_memory_scaling(
    results: list[BenchmarkResult],
    x_param: str = "n_variants",
    title: str = "Memory Usage Scaling",
    xlabel: str = "Problem Size",
    ylabel: str = "Peak Memory (MB)",
    log_scale: bool = True,
    output_path: Path | None = None,
) -> tuple[Any, Any]:
    """Plot memory usage scaling."""
    if not HAS_MATPLOTLIB:
        raise ImportError("matplotlib required for visualization")

    setup_publication_style()
    fig, ax = plt.subplots()

    x_values, y_values = [], []
    for r in results:
        if x_param in r.extra_info and "peak_memory_mb" in r.extra_info:
            x_values.append(r.extra_info[x_param])
            y_values.append(r.extra_info["peak_memory_mb"])

    if not x_values:
        raise ValueError(f"No results with {x_param} and peak_memory_mb")

    sorted_indices = np.argsort(x_values)
    x_values = np.array(x_values)[sorted_indices]
    y_values = np.array(y_values)[sorted_indices]

    ax.plot(x_values, y_values, marker="s", color=NATURE_COLORS["secondary"], linewidth=2)
    ax.fill_between(x_values, 0, y_values, alpha=0.3, color=NATURE_COLORS["secondary"])

    if log_scale:
        ax.set_xscale("log")
        ax.set_yscale("log")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3, linestyle="--")
    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, bbox_inches="tight", facecolor="white")
        print(f"Saved figure to {output_path}")

    return fig, ax


def plot_time_memory_comparison(
    results: list[BenchmarkResult],
    x_param: str = "n_variants",
    title: str = "Time vs Memory Trade-off",
    xlabel: str = "Problem Size",
    output_path: Path | None = None,
) -> tuple[Any, Any]:
    """Plot time and memory on dual y-axes."""
    if not HAS_MATPLOTLIB:
        raise ImportError("matplotlib required for visualization")

    setup_publication_style()
    fig, ax1 = plt.subplots()

    # Extract data
    x_values, time_values, memory_values = [], [], []
    for r in results:
        if x_param in r.extra_info and "peak_memory_mb" in r.extra_info:
            x_values.append(r.extra_info[x_param])
            time_values.append(r.mean)
            memory_values.append(r.extra_info["peak_memory_mb"])

    if not x_values:
        raise ValueError("No results with required data")

    sorted_indices = np.argsort(x_values)
    x_values = np.array(x_values)[sorted_indices]
    time_values = np.array(time_values)[sorted_indices]
    memory_values = np.array(memory_values)[sorted_indices]

    # Plot time on left axis
    color1 = NATURE_COLORS["wasp2"]
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel("Time (seconds)", color=color1)
    line1 = ax1.plot(x_values, time_values, marker="o", color=color1, label="Time")
    ax1.tick_params(axis="y", labelcolor=color1)
    ax1.set_xscale("log")
    ax1.set_yscale("log")

    # Plot memory on right axis
    ax2 = ax1.twinx()
    color2 = NATURE_COLORS["phaser"]
    ax2.set_ylabel("Peak Memory (MB)", color=color2)
    line2 = ax2.plot(x_values, memory_values, marker="s", color=color2, label="Memory")
    ax2.tick_params(axis="y", labelcolor=color2)
    ax2.set_yscale("log")

    # Combined legend
    lines = line1 + line2
    labels = [l.get_label() for l in lines]
    ax1.legend(lines, labels, loc="upper left")

    ax1.set_title(title)
    ax1.grid(True, alpha=0.3, linestyle="--")
    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, bbox_inches="tight", facecolor="white")
        print(f"Saved figure to {output_path}")

    return fig, (ax1, ax2)


# ============================================================================
# Multi-panel figures for papers
# ============================================================================


def generate_paper_figure(
    results: list[BenchmarkResult],
    output_path: Path,
    title: str = "WASP2 Performance Benchmarks",
) -> Any:
    """Generate a multi-panel figure suitable for publication."""
    if not HAS_MATPLOTLIB:
        raise ImportError("matplotlib required for visualization")

    setup_publication_style()
    fig = plt.figure(figsize=(14, 10))

    # Panel A: Variant scaling
    ax1 = fig.add_subplot(2, 2, 1)
    variant_results = filter_by_group(results, "variant_scaling")
    if variant_results:
        plot_scaling(
            variant_results,
            x_param="n_variants",
            title="A. Variant Scaling",
            xlabel="Number of Variants",
            ax=ax1,
        )

    # Panel B: Sample scaling
    ax2 = fig.add_subplot(2, 2, 2)
    sample_results = filter_by_group(results, "sample_scaling")
    if sample_results:
        plot_scaling(
            sample_results,
            x_param="n_samples",
            title="B. Sample Scaling",
            xlabel="Number of Samples",
            ax=ax2,
        )

    # Panel C: Tool comparison
    ax3 = fig.add_subplot(2, 2, 3)
    comparison_results = [
        r for r in results if r.group.startswith("tool_comparison") and "tool" in r.extra_info
    ]
    if comparison_results:
        # Plot on existing axes
        tool_groups = group_by_extra_info(comparison_results, "tool")
        for tool_name, tool_results in tool_groups.items():
            x_values, y_values = [], []
            for r in tool_results:
                if "n_variants" in r.extra_info:
                    x_values.append(r.extra_info["n_variants"])
                    y_values.append(r.mean)
            if x_values:
                sorted_indices = np.argsort(x_values)
                x_values = np.array(x_values)[sorted_indices]
                y_values = np.array(y_values)[sorted_indices]
                ax3.plot(x_values, y_values, marker="o", label=tool_name, color=get_tool_color(tool_name))
        ax3.set_xscale("log")
        ax3.set_yscale("log")
        ax3.set_xlabel("Number of Variants")
        ax3.set_ylabel("Time (seconds)")
        ax3.set_title("C. Tool Comparison")
        ax3.legend()
        ax3.grid(True, alpha=0.3, linestyle="--")

    # Panel D: Memory scaling
    ax4 = fig.add_subplot(2, 2, 4)
    memory_results = [r for r in results if "peak_memory_mb" in r.extra_info]
    if memory_results:
        x_values, y_values = [], []
        for r in memory_results:
            if "n_variants" in r.extra_info:
                x_values.append(r.extra_info["n_variants"])
                y_values.append(r.extra_info["peak_memory_mb"])
        if x_values:
            sorted_indices = np.argsort(x_values)
            x_values = np.array(x_values)[sorted_indices]
            y_values = np.array(y_values)[sorted_indices]
            ax4.plot(x_values, y_values, marker="s", color=NATURE_COLORS["secondary"])
            ax4.fill_between(x_values, 0, y_values, alpha=0.3, color=NATURE_COLORS["secondary"])
            ax4.set_xscale("log")
            ax4.set_yscale("log")
    ax4.set_xlabel("Number of Variants")
    ax4.set_ylabel("Peak Memory (MB)")
    ax4.set_title("D. Memory Usage")
    ax4.grid(True, alpha=0.3, linestyle="--")

    fig.suptitle(title, fontsize=18, y=1.02)
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches="tight", facecolor="white", dpi=300)
    print(f"Saved paper figure to {output_path}")

    return fig


# ============================================================================
# Main figure generation
# ============================================================================


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

    # Variant scaling
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

    # Sample scaling
    sample_scaling = filter_by_group(results, "sample_scaling")
    if sample_scaling:
        for fmt in formats:
            plot_scaling(
                sample_scaling,
                x_param="n_samples",
                title="Performance vs. Sample Count",
                xlabel="Number of Samples",
                output_path=output_dir / f"sample_scaling.{fmt}",
            )
            plt.close()

    # Sample × variant matrix heatmap
    matrix_results = filter_by_group(results, "sample_variant_matrix")
    if matrix_results:
        for fmt in formats:
            plot_scaling_heatmap(
                matrix_results,
                x_param="n_samples",
                y_param="n_variants",
                title="Sample × Variant Scaling Matrix",
                xlabel="Number of Samples",
                ylabel="Number of Variants",
                output_path=output_dir / f"scaling_heatmap.{fmt}",
            )
            plt.close()

    # Tool comparison
    tool_results = [
        r for r in results if r.group.startswith("tool_comparison") and "tool" in r.extra_info
    ]
    if tool_results:
        for fmt in formats:
            plot_tool_comparison(
                tool_results,
                title="WASP2 vs Competitors",
                output_path=output_dir / f"tool_comparison.{fmt}",
            )
            plt.close()

            plot_tool_comparison_bars(
                tool_results,
                title="Tool Performance by Scale",
                output_path=output_dir / f"tool_comparison_bars.{fmt}",
            )
            plt.close()

    # Memory scaling
    memory_results = [r for r in results if "peak_memory_mb" in r.extra_info]
    if memory_results:
        for fmt in formats:
            plot_memory_scaling(
                memory_results,
                title="Memory Usage Scaling",
                output_path=output_dir / f"memory_scaling.{fmt}",
            )
            plt.close()

    # Time vs memory comparison
    time_memory_results = [
        r for r in results
        if "peak_memory_mb" in r.extra_info and "n_variants" in r.extra_info
    ]
    if time_memory_results:
        for fmt in formats:
            plot_time_memory_comparison(
                time_memory_results,
                title="Time vs Memory Trade-off",
                output_path=output_dir / f"time_memory_comparison.{fmt}",
            )
            plt.close()

    # Combined paper figure
    generate_paper_figure(
        results,
        output_path=output_dir / "paper_figure.pdf",
        title="WASP2 Performance Benchmarks",
    )
    plt.close()

    print(f"Generated all figures in {output_dir}")
