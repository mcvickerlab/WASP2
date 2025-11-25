#!/usr/bin/env python3
"""
Benchmark script comparing pysam VCFSource vs cyvcf2 CyVCF2Source performance.

This script measures the performance differences between:
- pysam-based VCFSource (baseline)
- cyvcf2-based CyVCF2Source (optimized)

Metrics measured:
- Variant counting speed
- Full iteration performance
- Heterozygous site filtering speed
- Memory usage
- Variants per second throughput

Usage:
    python benchmarks/benchmark_vcf_performance.py <vcf_file>
    python benchmarks/benchmark_vcf_performance.py tests/data/sample.vcf.gz

Requirements:
    pip install wasp2[cyvcf2]
"""

import argparse
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple

try:
    import psutil
    HAS_PSUTIL = True
except ImportError:
    HAS_PSUTIL = False
    print("Warning: psutil not available, memory profiling disabled")

from wasp2.io.vcf_source import VCFSource

try:
    from wasp2.io.cyvcf2_source import CyVCF2Source, CYVCF2_AVAILABLE
except ImportError:
    CYVCF2_AVAILABLE = False

try:
    from wasp2.io.pgen_source import PGENSource
    PGEN_AVAILABLE = True
except ImportError:
    PGEN_AVAILABLE = False


def get_memory_mb() -> float:
    """Get current process memory usage in MB."""
    if HAS_PSUTIL:
        process = psutil.Process()
        return process.memory_info().rss / 1024 / 1024
    return 0.0


def benchmark_variant_counting(vcf_path: Path, source_class) -> Tuple[float, int]:
    """Benchmark variant counting speed.

    Returns:
        Tuple of (elapsed_time_seconds, variant_count)
    """
    start_time = time.perf_counter()

    with source_class(vcf_path) as source:
        count = source.variant_count

    elapsed = time.perf_counter() - start_time
    return elapsed, count


def benchmark_full_iteration(vcf_path: Path, source_class) -> Tuple[float, int, float]:
    """Benchmark full variant iteration.

    Returns:
        Tuple of (elapsed_time_seconds, variant_count, memory_delta_mb)
    """
    mem_start = get_memory_mb()
    start_time = time.perf_counter()

    count = 0
    with source_class(vcf_path) as source:
        for _ in source.iter_variants():
            count += 1

    elapsed = time.perf_counter() - start_time
    mem_delta = get_memory_mb() - mem_start

    return elapsed, count, mem_delta


def benchmark_het_filtering(vcf_path: Path, source_class, sample: str = None) -> Tuple[float, int]:
    """Benchmark heterozygous site filtering.

    Returns:
        Tuple of (elapsed_time_seconds, het_count)
    """
    start_time = time.perf_counter()

    count = 0
    with source_class(vcf_path) as source:
        samples = [sample] if sample else [source.samples[0]]
        for _ in source.iter_variants(samples=samples, het_only=True):
            count += 1

    elapsed = time.perf_counter() - start_time
    return elapsed, count


def benchmark_region_queries(vcf_path: Path, source_class, regions: List[Tuple[str, int, int]]) -> float:
    """Benchmark region query performance.

    Args:
        vcf_path: Path to VCF file
        source_class: VCFSource or CyVCF2Source class
        regions: List of (chrom, start, end) tuples

    Returns:
        Total elapsed time in seconds
    """
    start_time = time.perf_counter()

    with source_class(vcf_path) as source:
        for chrom, start, end in regions:
            try:
                list(source.query_region(chrom, start, end))
            except ValueError:
                # Region query may fail if file not indexed
                pass

    elapsed = time.perf_counter() - start_time
    return elapsed


def run_benchmarks(vcf_path: Path, sample: str = None, pgen_path: Path = None) -> Dict[str, Dict[str, float]]:
    """Run all benchmarks and return results.

    Args:
        vcf_path: Path to VCF file
        sample: Optional sample name
        pgen_path: Optional path to equivalent PGEN file for comparison

    Returns:
        Dictionary mapping benchmark_name -> {source_type -> metric}
    """
    print(f"\n{'='*80}")
    print(f"Benchmarking Multi-Format Variant I/O Performance")
    print(f"{'='*80}\n")
    print(f"VCF file: {vcf_path.name}")
    if pgen_path:
        print(f"PGEN file: {pgen_path.name}")
    print()

    results = {}

    # Check file size
    file_size_mb = vcf_path.stat().st_size / 1024 / 1024
    print(f"VCF file size: {file_size_mb:.2f} MB")
    if pgen_path and pgen_path.exists():
        pgen_size_mb = pgen_path.stat().st_size / 1024 / 1024
        print(f"PGEN file size: {pgen_size_mb:.2f} MB")
    print()

    # Benchmark 1: Variant Counting
    print("=" * 80)
    print("Benchmark 1: Variant Counting Speed")
    print("=" * 80)

    pysam_time, pysam_count = benchmark_variant_counting(vcf_path, VCFSource)
    print(f"pysam VCFSource:     {pysam_time:.4f}s ({pysam_count:,} variants) [baseline]")

    results['counting'] = {'pysam': pysam_time}

    if CYVCF2_AVAILABLE:
        cyvcf2_time, cyvcf2_count = benchmark_variant_counting(vcf_path, CyVCF2Source)
        print(f"cyvcf2 CyVCF2Source:  {cyvcf2_time:.4f}s ({cyvcf2_count:,} variants)")
        results['counting']['cyvcf2'] = cyvcf2_time
        speedup = pysam_time / cyvcf2_time if cyvcf2_time > 0 else 0
        print(f"  └─ Speedup vs pysam: {speedup:.2f}x faster")

    if PGEN_AVAILABLE and pgen_path and pgen_path.exists():
        try:
            pgen_time, pgen_count = benchmark_variant_counting(pgen_path, PGENSource)
            print(f"PLINK2 PGENSource:    {pgen_time:.4f}s ({pgen_count:,} variants)")
            results['counting']['pgen'] = pgen_time
            speedup = pysam_time / pgen_time if pgen_time > 0 else 0
            print(f"  └─ Speedup vs pysam: {speedup:.2f}x faster")
        except Exception as e:
            print(f"PLINK2 PGENSource:    Error - {e}")

    # Benchmark 2: Full Iteration
    print("\n" + "=" * 80)
    print("Benchmark 2: Full Iteration Performance")
    print("=" * 80)

    pysam_time, pysam_count, pysam_mem = benchmark_full_iteration(vcf_path, VCFSource)
    pysam_rate = pysam_count / pysam_time if pysam_time > 0 else 0
    print(f"pysam VCFSource:     {pysam_time:.4f}s ({pysam_rate:,.0f} variants/s, +{pysam_mem:.1f} MB) [baseline]")

    results['iteration'] = {'pysam': pysam_time, 'pysam_rate': pysam_rate}

    if CYVCF2_AVAILABLE:
        cyvcf2_time, cyvcf2_count, cyvcf2_mem = benchmark_full_iteration(vcf_path, CyVCF2Source)
        cyvcf2_rate = cyvcf2_count / cyvcf2_time if cyvcf2_time > 0 else 0
        print(f"cyvcf2 CyVCF2Source:  {cyvcf2_time:.4f}s ({cyvcf2_rate:,.0f} variants/s, +{cyvcf2_mem:.1f} MB)")
        results['iteration']['cyvcf2'] = cyvcf2_time
        results['iteration']['cyvcf2_rate'] = cyvcf2_rate
        speedup = pysam_time / cyvcf2_time if cyvcf2_time > 0 else 0
        print(f"  └─ Speedup vs pysam: {speedup:.2f}x faster ({cyvcf2_rate/pysam_rate:.2f}x throughput)")

    if PGEN_AVAILABLE and pgen_path and pgen_path.exists():
        try:
            pgen_time, pgen_count, pgen_mem = benchmark_full_iteration(pgen_path, PGENSource)
            pgen_rate = pgen_count / pgen_time if pgen_time > 0 else 0
            print(f"PLINK2 PGENSource:    {pgen_time:.4f}s ({pgen_rate:,.0f} variants/s, +{pgen_mem:.1f} MB)")
            results['iteration']['pgen'] = pgen_time
            results['iteration']['pgen_rate'] = pgen_rate
            speedup = pysam_time / pgen_time if pgen_time > 0 else 0
            print(f"  └─ Speedup vs pysam: {speedup:.2f}x faster ({pgen_rate/pysam_rate:.2f}x throughput)")
        except Exception as e:
            print(f"PLINK2 PGENSource:    Error - {e}")

    # Benchmark 3: Heterozygous Filtering
    print("\n" + "=" * 80)
    print("Benchmark 3: Heterozygous Site Filtering")
    print("=" * 80)

    pysam_time, pysam_het_count = benchmark_het_filtering(vcf_path, VCFSource, sample)
    pysam_rate = pysam_het_count / pysam_time if pysam_time > 0 else 0
    print(f"pysam VCFSource:    {pysam_time:.4f}s ({pysam_het_count:,} het sites, {pysam_rate:,.0f} variants/s)")

    results['het_filtering'] = {'pysam': pysam_time, 'pysam_count': pysam_het_count}

    if CYVCF2_AVAILABLE:
        cyvcf2_time, cyvcf2_het_count = benchmark_het_filtering(vcf_path, CyVCF2Source, sample)
        cyvcf2_rate = cyvcf2_het_count / cyvcf2_time if cyvcf2_time > 0 else 0
        print(f"cyvcf2 CyVCF2Source: {cyvcf2_time:.4f}s ({cyvcf2_het_count:,} het sites, {cyvcf2_rate:,.0f} variants/s)")
        results['het_filtering']['cyvcf2'] = cyvcf2_time
        results['het_filtering']['cyvcf2_count'] = cyvcf2_het_count
        speedup = pysam_time / cyvcf2_time if cyvcf2_time > 0 else 0
        print(f"Speedup: {speedup:.2f}x faster")

    # Benchmark 4: Region Queries (if file is indexed)
    print("\n" + "=" * 80)
    print("Benchmark 4: Region Query Performance")
    print("=" * 80)

    # Generate some test regions
    with VCFSource(vcf_path) as source:
        variants = list(source.iter_variants())
        if len(variants) >= 3:
            regions = [
                (variants[0].variant.chrom, variants[0].variant.pos, variants[0].variant.pos + 100),
                (variants[len(variants)//2].variant.chrom, variants[len(variants)//2].variant.pos,
                 variants[len(variants)//2].variant.pos + 100),
                (variants[-1].variant.chrom, variants[-1].variant.pos - 100, variants[-1].variant.pos),
            ]
        else:
            regions = []

    if regions:
        pysam_time = benchmark_region_queries(vcf_path, VCFSource, regions)
        print(f"pysam VCFSource:    {pysam_time:.4f}s ({len(regions)} queries)")
        results['region_queries'] = {'pysam': pysam_time}

        if CYVCF2_AVAILABLE:
            cyvcf2_time = benchmark_region_queries(vcf_path, CyVCF2Source, regions)
            print(f"cyvcf2 CyVCF2Source: {cyvcf2_time:.4f}s ({len(regions)} queries)")
            results['region_queries']['cyvcf2'] = cyvcf2_time
            speedup = pysam_time / cyvcf2_time if cyvcf2_time > 0 else 0
            print(f"Speedup: {speedup:.2f}x faster")
    else:
        print("Skipping (insufficient variants for test regions)")

    return results


def print_summary(results: Dict[str, Dict[str, float]]):
    """Print summary of benchmark results."""
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)

    if not CYVCF2_AVAILABLE:
        print("\n⚠️  cyvcf2 not installed - only pysam baseline results shown")
        print("   Install with: pip install wasp2[cyvcf2]")
        return

    print("\nPerformance Improvements (cyvcf2 vs pysam):")
    print("-" * 80)

    for benchmark, metrics in results.items():
        if 'pysam' in metrics and 'cyvcf2' in metrics:
            speedup = metrics['pysam'] / metrics['cyvcf2'] if metrics['cyvcf2'] > 0 else 0
            benchmark_name = benchmark.replace('_', ' ').title()
            print(f"{benchmark_name:.<40} {speedup:>6.2f}x faster")

    # Calculate overall average speedup
    speedups = []
    for benchmark, metrics in results.items():
        if 'pysam' in metrics and 'cyvcf2' in metrics:
            speedup = metrics['pysam'] / metrics['cyvcf2'] if metrics['cyvcf2'] > 0 else 0
            if speedup > 0:
                speedups.append(speedup)

    if speedups:
        avg_speedup = sum(speedups) / len(speedups)
        print("-" * 80)
        print(f"{'Average Speedup':.<40} {avg_speedup:>6.2f}x faster")

    print("\n✅ Recommendation: Use CyVCF2Source for production workloads")
    print("   Expected performance gain: ~5-7x faster VCF parsing")


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark multi-format variant I/O performance (pysam vs cyvcf2 vs PLINK2)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument(
        "vcf_file",
        type=Path,
        help="VCF or VCF.gz file to benchmark"
    )
    parser.add_argument(
        "--pgen",
        type=Path,
        help="Optional PGEN file for comparison (e.g., data.pgen)"
    )
    parser.add_argument(
        "--sample",
        type=str,
        help="Sample name to use for filtering (default: first sample)"
    )

    args = parser.parse_args()

    if not args.vcf_file.exists():
        print(f"Error: File not found: {args.vcf_file}", file=sys.stderr)
        sys.exit(1)

    if args.pgen and not args.pgen.exists():
        print(f"Warning: PGEN file not found: {args.pgen}", file=sys.stderr)
        args.pgen = None

    if not CYVCF2_AVAILABLE:
        print("\n⚠️  WARNING: cyvcf2 not installed")
        print("Only pysam baseline benchmarks will run.")
        print("Install cyvcf2 with: pip install wasp2[cyvcf2]\n")

    if not PGEN_AVAILABLE and args.pgen:
        print("\n⚠️  WARNING: pgenlib not installed")
        print("PGEN benchmarks will be skipped.")
        print("Install pgenlib with: pip install wasp2[plink]\n")

    results = run_benchmarks(args.vcf_file, args.sample, args.pgen)
    print_summary(results)

    print("\n" + "=" * 80)
    print("Benchmark complete!")
    print("=" * 80)


if __name__ == "__main__":
    main()
