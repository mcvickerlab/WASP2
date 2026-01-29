#!/usr/bin/env python3
"""
WASP2 Benchmark Runner

Standalone script for running reproducible benchmarks to validate
WASP2's performance claims.

Performance Claims to Validate:
- 61x faster WASP filtering (vs WASP v1)
- 6.4x faster counting (vs phASER)
- r² > 0.99 concordance with GATK ASEReadCounter

Usage:
    python benchmarking/run_benchmarks.py --all
    python benchmarking/run_benchmarks.py --counting
    python benchmarking/run_benchmarks.py --mapping
    python benchmarking/run_benchmarks.py --concordance
"""

import argparse
import sys
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from benchmarking.utils import (
    check_tool,
    print_comparison_table,
    save_results,
)


def run_counting_benchmarks(
    n_variants: int = 10000,
    n_regions: int = 1000,
    iterations: int = 5,
    warmup: int = 2,
) -> list:
    """Run counting speed benchmarks: WASP2 vs phASER vs GATK."""
    from benchmarking.scripts.benchmark_counting import (
        benchmark_counting_speed,
    )

    print("\n" + "=" * 70)
    print("COUNTING SPEED BENCHMARKS")
    print(f"Variants: {n_variants:,}, Regions: {n_regions:,}")
    print("=" * 70)

    results = benchmark_counting_speed(
        n_variants=n_variants,
        n_regions=n_regions,
        iterations=iterations,
        warmup=warmup,
    )

    print_comparison_table(results)
    return results


def run_mapping_benchmarks(
    n_reads: int = 10000,
    iterations: int = 5,
    warmup: int = 2,
) -> list:
    """Run mapping filter benchmarks: WASP2 vs WASP v1."""
    from benchmarking.scripts.benchmark_mapping import (
        benchmark_mapping_filter,
    )

    print("\n" + "=" * 70)
    print("MAPPING FILTER BENCHMARKS")
    print(f"Reads: {n_reads:,}")
    print("=" * 70)

    results = benchmark_mapping_filter(
        n_reads=n_reads,
        iterations=iterations,
        warmup=warmup,
    )

    print_comparison_table(results)
    return results


def run_concordance_benchmarks(
    n_variants: int = 1000,
) -> dict:
    """Run concordance validation: WASP2 vs GATK ASEReadCounter."""
    from benchmarking.scripts.benchmark_concordance import (
        validate_concordance,
    )

    print("\n" + "=" * 70)
    print("CONCORDANCE VALIDATION")
    print(f"Variants: {n_variants:,}")
    print("=" * 70)

    results = validate_concordance(n_variants=n_variants)
    return results


def print_tool_availability() -> None:
    """Print availability status of external tools."""
    tools = {
        "phASER": ["phaser.py", "phaser"],
        "WASP v1": ["wasp", "rmdup_pe.py"],
        "GATK": ["gatk"],
        "samtools": ["samtools"],
        "bcftools": ["bcftools"],
    }

    print("\nExternal Tool Availability:")
    print("-" * 40)

    for name, executables in tools.items():
        available = any(check_tool(exe) for exe in executables)
        status = "✓" if available else "✗"
        print(f"  {status} {name}")

    print()


def main() -> int:
    parser = argparse.ArgumentParser(
        description="WASP2 Benchmark Runner",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    parser.add_argument(
        "--all",
        "-a",
        action="store_true",
        help="Run all benchmarks",
    )
    parser.add_argument(
        "--counting",
        "-c",
        action="store_true",
        help="Run counting speed benchmarks",
    )
    parser.add_argument(
        "--mapping",
        "-m",
        action="store_true",
        help="Run mapping filter benchmarks",
    )
    parser.add_argument(
        "--concordance",
        action="store_true",
        help="Run concordance validation",
    )
    parser.add_argument(
        "--quick",
        "-q",
        action="store_true",
        help="Run quick benchmarks with reduced iterations",
    )
    parser.add_argument(
        "--n-variants",
        type=int,
        default=10000,
        help="Number of variants for benchmarks (default: 10000)",
    )
    parser.add_argument(
        "--n-reads",
        type=int,
        default=10000,
        help="Number of reads for mapping benchmarks (default: 10000)",
    )
    parser.add_argument(
        "--iterations",
        type=int,
        default=5,
        help="Number of benchmark iterations (default: 5)",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        default=Path("benchmarking/results/benchmark_results.json"),
        help="Output file for results",
    )
    parser.add_argument(
        "--check-tools",
        action="store_true",
        help="Check external tool availability and exit",
    )

    args = parser.parse_args()

    if args.check_tools:
        print_tool_availability()
        return 0

    if args.quick:
        args.iterations = 2
        args.n_variants = 1000
        args.n_reads = 1000

    # Default to all if no specific benchmark selected
    if not any([args.all, args.counting, args.mapping, args.concordance]):
        args.all = True

    print_tool_availability()

    all_results = []
    warmup = 1 if args.quick else 2

    try:
        if args.all or args.counting:
            results = run_counting_benchmarks(
                n_variants=args.n_variants,
                n_regions=max(100, args.n_variants // 10),
                iterations=args.iterations,
                warmup=warmup,
            )
            all_results.extend(results)

        if args.all or args.mapping:
            results = run_mapping_benchmarks(
                n_reads=args.n_reads,
                iterations=args.iterations,
                warmup=warmup,
            )
            all_results.extend(results)

        if args.all or args.concordance:
            run_concordance_benchmarks(n_variants=min(1000, args.n_variants))

    except ImportError as e:
        print(f"Error: Missing dependency - {e}")
        print("Install with: pip install wasp2[benchmark]")
        return 1
    except Exception as e:
        import traceback

        print(f"Error running benchmarks: {e}")
        traceback.print_exc()
        return 1

    if all_results:
        save_results(all_results, args.output)

    print("\n✓ Benchmarks completed successfully")
    return 0


if __name__ == "__main__":
    sys.exit(main())
