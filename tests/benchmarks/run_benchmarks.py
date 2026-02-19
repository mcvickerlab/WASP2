#!/usr/bin/env python3
"""
WASP2 Benchmark Runner

A CLI tool for running the complete WASP2 benchmark suite and generating
publication-quality figures.

Usage:
    python run_benchmarks.py                    # Run all benchmarks
    python run_benchmarks.py --quick            # Run quick subset
    python run_benchmarks.py --groups scaling   # Run specific groups
    python run_benchmarks.py --figures-only     # Generate figures from existing results
"""

import argparse
import importlib.util
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

# ============================================================================
# Configuration
# ============================================================================

BENCHMARK_GROUPS = {
    "variant_scaling": "tests/benchmarks/test_scaling_benchmarks.py::TestVariantScaling",
    "region_scaling": "tests/benchmarks/test_scaling_benchmarks.py::TestRegionScaling",
    "method_comparison": "tests/benchmarks/test_scaling_benchmarks.py::TestMethodComparison",
    "memory_scaling": "tests/benchmarks/test_scaling_benchmarks.py::TestMemoryScaling",
    "sample_scaling": "tests/benchmarks/test_sample_scaling.py::TestSampleScaling",
    "sample_memory": "tests/benchmarks/test_sample_scaling.py::TestSampleMemoryScaling",
    "sample_variant_matrix": "tests/benchmarks/test_sample_scaling.py::TestSampleVariantMatrix",
    "cohort_simulation": "tests/benchmarks/test_sample_scaling.py::TestCohortSimulation",
    "coverage_scaling": "tests/benchmarks/test_sample_scaling.py::TestHighThroughputScaling",
    "tool_comparison": "tests/benchmarks/test_tool_comparison.py",
    "analysis": "tests/benchmarks/test_analysis_benchmarks.py",
}

QUICK_GROUPS = ["variant_scaling", "sample_scaling", "analysis"]

DEFAULT_OUTPUT_DIR = Path(".benchmarks")
FIGURES_OUTPUT_DIR = Path("benchmark_figures")


# ============================================================================
# Runner functions
# ============================================================================


def get_project_root() -> Path:
    """Get the project root directory."""
    current = Path(__file__).resolve()
    for parent in current.parents:
        if (parent / "pyproject.toml").exists():
            return parent
    return current.parent.parent.parent


def check_dependencies() -> dict[str, bool]:
    """Check for required and optional dependencies."""
    deps = {
        "pytest": importlib.util.find_spec("pytest") is not None,
        "pytest-benchmark": importlib.util.find_spec("pytest_benchmark") is not None,
        "memory-profiler": importlib.util.find_spec("memory_profiler") is not None,
        "matplotlib": importlib.util.find_spec("matplotlib") is not None,
        "seaborn": importlib.util.find_spec("seaborn") is not None,
    }
    return deps


def run_pytest_benchmarks(
    groups: list[str] | None = None,
    output_dir: Path | None = None,
    extra_args: list[str] | None = None,
    verbose: bool = True,
    skip_slow: bool = False,
) -> tuple[bool, Path | None]:
    """
    Run pytest benchmarks for specified groups.

    Returns:
        Tuple of (success, results_path)
    """
    project_root = get_project_root()
    output_dir = output_dir or DEFAULT_OUTPUT_DIR
    output_dir = project_root / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    # Build pytest command
    cmd = ["python", "-m", "pytest"]

    # Add benchmark-specific options
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    json_output = output_dir / f"benchmark_{timestamp}.json"
    cmd.extend(
        [
            "--benchmark-only",
            "--benchmark-json",
            str(json_output),
            "--benchmark-columns=mean,stddev,min,max,rounds",
            "--benchmark-sort=mean",
        ]
    )

    # Add test paths for specified groups
    if groups:
        for group in groups:
            if group in BENCHMARK_GROUPS:
                cmd.append(str(project_root / BENCHMARK_GROUPS[group]))
            else:
                print(f"Warning: Unknown benchmark group '{group}'")
    else:
        # Run all benchmark tests
        cmd.append(str(project_root / "tests/benchmarks/"))

    # Skip slow tests if requested
    if skip_slow:
        cmd.extend(["-m", "not slow"])

    # Add extra args
    if extra_args:
        cmd.extend(extra_args)

    # Add verbosity
    if verbose:
        cmd.append("-v")

    print(f"Running: {' '.join(cmd)}")
    print("-" * 60)

    try:
        result = subprocess.run(
            cmd,
            cwd=str(project_root),
            check=False,
        )

        if json_output.exists():
            print(f"\nBenchmark results saved to: {json_output}")
            return result.returncode == 0, json_output
        else:
            print("\nWarning: No benchmark results file generated")
            return False, None

    except Exception as e:
        print(f"Error running benchmarks: {e}")
        return False, None


def find_latest_results(output_dir: Path | None = None) -> Path | None:
    """Find the most recent benchmark results file."""
    output_dir = output_dir or DEFAULT_OUTPUT_DIR
    project_root = get_project_root()
    results_dir = project_root / output_dir

    if not results_dir.exists():
        return None

    json_files = list(results_dir.glob("benchmark_*.json"))
    if not json_files:
        # Try pytest-benchmark default location
        json_files = list(results_dir.glob("*/*.json"))

    if not json_files:
        return None

    return max(json_files, key=lambda p: p.stat().st_mtime)


def generate_figures(
    results_path: Path | None = None,
    output_dir: Path | None = None,
    formats: list[str] | None = None,
) -> bool:
    """Generate benchmark figures from results."""
    project_root = get_project_root()

    if results_path is None:
        results_path = find_latest_results()
        if results_path is None:
            print("Error: No benchmark results found. Run benchmarks first.")
            return False

    output_dir = project_root / (output_dir or FIGURES_OUTPUT_DIR)
    formats = formats or ["png", "pdf"]

    print(f"Generating figures from: {results_path}")
    print(f"Output directory: {output_dir}")

    try:
        from .utils.visualization import generate_all_figures

        generate_all_figures(results_path, output_dir, formats)
        return True
    except ImportError as e:
        print(f"Error: Missing visualization dependencies: {e}")
        print("Install with: pip install matplotlib seaborn")
        return False
    except Exception as e:
        print(f"Error generating figures: {e}")
        return False


def print_results_summary(results_path: Path) -> None:
    """Print a summary of benchmark results."""
    with open(results_path) as f:
        data = json.load(f)

    benchmarks = data.get("benchmarks", [])
    if not benchmarks:
        print("No benchmark results found")
        return

    print("\n" + "=" * 70)
    print("BENCHMARK RESULTS SUMMARY")
    print("=" * 70)

    # Group by benchmark group
    from collections import defaultdict

    groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for b in benchmarks:
        groups[b.get("group", "ungrouped")].append(b)

    for group_name, group_benchmarks in sorted(groups.items()):
        print(f"\n{group_name.upper()}")
        print("-" * 40)

        for b in sorted(group_benchmarks, key=lambda x: x["stats"]["mean"]):
            name = b["name"].split("::")[-1]
            mean = b["stats"]["mean"]
            stddev = b["stats"]["stddev"]
            extra = b.get("extra_info", {})

            # Format time appropriately
            if mean < 0.001:
                time_str = f"{mean * 1_000_000:.2f} μs"
            elif mean < 1:
                time_str = f"{mean * 1000:.2f} ms"
            else:
                time_str = f"{mean:.3f} s"

            # Build extra info string
            extra_parts = []
            if "n_variants" in extra:
                extra_parts.append(f"variants={extra['n_variants']:,}")
            if "n_samples" in extra:
                extra_parts.append(f"samples={extra['n_samples']}")
            if "peak_memory_mb" in extra:
                extra_parts.append(f"mem={extra['peak_memory_mb']:.1f}MB")
            if "tool" in extra:
                extra_parts.append(f"tool={extra['tool']}")

            extra_str = f" ({', '.join(extra_parts)})" if extra_parts else ""
            print(f"  {name}: {time_str} ± {stddev * 1000:.2f}ms{extra_str}")

    print("\n" + "=" * 70)

    # Print machine info
    machine = data.get("machine_info", {})
    if machine:
        print(f"Machine: {machine.get('node', 'unknown')}")
        print(f"CPU: {machine.get('processor', 'unknown')}")
        print(f"Python: {machine.get('python_version', 'unknown')}")


def list_groups() -> None:
    """List available benchmark groups."""
    print("Available benchmark groups:")
    print("-" * 40)
    for name, path in sorted(BENCHMARK_GROUPS.items()):
        print(f"  {name:25s} -> {path}")
    print("\nQuick groups (--quick):", ", ".join(QUICK_GROUPS))


# ============================================================================
# Main CLI
# ============================================================================


def main() -> int:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="WASP2 Benchmark Runner",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    parser.add_argument(
        "--groups",
        "-g",
        nargs="+",
        help="Specific benchmark groups to run",
    )
    parser.add_argument(
        "--quick",
        "-q",
        action="store_true",
        help="Run quick subset of benchmarks",
    )
    parser.add_argument(
        "--skip-slow",
        action="store_true",
        help="Skip benchmarks marked as slow",
    )
    parser.add_argument(
        "--figures-only",
        "-f",
        action="store_true",
        help="Only generate figures from existing results",
    )
    parser.add_argument(
        "--results",
        "-r",
        type=Path,
        help="Path to benchmark results JSON file",
    )
    parser.add_argument(
        "--output-dir",
        "-o",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for benchmark results",
    )
    parser.add_argument(
        "--figures-dir",
        type=Path,
        default=FIGURES_OUTPUT_DIR,
        help="Directory for generated figures",
    )
    parser.add_argument(
        "--formats",
        nargs="+",
        default=["png", "pdf"],
        help="Figure output formats (default: png pdf)",
    )
    parser.add_argument(
        "--list-groups",
        "-l",
        action="store_true",
        help="List available benchmark groups",
    )
    parser.add_argument(
        "--check-deps",
        action="store_true",
        help="Check required dependencies",
    )
    parser.add_argument(
        "--no-figures",
        action="store_true",
        help="Skip figure generation after benchmarks",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        default=True,
        help="Verbose output",
    )
    parser.add_argument(
        "extra_args",
        nargs="*",
        help="Additional arguments to pass to pytest",
    )

    args = parser.parse_args()

    # Handle special commands
    if args.list_groups:
        list_groups()
        return 0

    if args.check_deps:
        deps = check_dependencies()
        print("Dependency status:")
        for dep, available in deps.items():
            status = "✓" if available else "✗"
            print(f"  {status} {dep}")
        missing = [d for d, a in deps.items() if not a]
        if missing:
            print(f"\nMissing dependencies: {', '.join(missing)}")
            print("Install with: pip install wasp2[benchmark]")
            return 1
        return 0

    # Check dependencies
    deps = check_dependencies()
    if not deps["pytest"] or not deps["pytest-benchmark"]:
        print("Error: pytest and pytest-benchmark are required")
        print("Install with: pip install pytest pytest-benchmark")
        return 1

    # Figures only mode
    if args.figures_only:
        success = generate_figures(
            results_path=args.results,
            output_dir=args.figures_dir,
            formats=args.formats,
        )
        return 0 if success else 1

    # Determine groups to run
    groups = args.groups
    if args.quick:
        groups = QUICK_GROUPS

    # Run benchmarks
    success, results_path = run_pytest_benchmarks(
        groups=groups,
        output_dir=args.output_dir,
        extra_args=args.extra_args,
        verbose=args.verbose,
        skip_slow=args.skip_slow,
    )

    # Print summary
    if results_path and results_path.exists():
        print_results_summary(results_path)

    # Generate figures
    if success and not args.no_figures and results_path:
        print("\nGenerating figures...")
        generate_figures(
            results_path=results_path,
            output_dir=args.figures_dir,
            formats=args.formats,
        )

    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
