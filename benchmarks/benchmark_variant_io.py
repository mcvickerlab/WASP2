#!/usr/bin/env python
"""
Benchmark: VCF vs PGEN Variant I/O Performance

This script benchmarks the performance of reading variant files and
extracting heterozygous sites - a critical operation in both the
WASP2 mapping and counting pipelines.

Usage:
    python benchmarks/benchmark_variant_io.py --vcf data.vcf.gz --pgen data.pgen --sample SAMPLE_ID

Output:
    - Console summary with timing results
    - CSV file with detailed metrics
    - Publication-quality figure (PDF/PNG) following Nature guidelines
"""

import argparse
import gc
import os
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from wasp2.io import VariantSource, Genotype


@dataclass
class BenchmarkResult:
    """Results from a single benchmark run."""
    format: str
    operation: str
    file_size_mb: float
    variant_count: int
    sample_count: int
    het_count: int
    time_seconds: float
    memory_mb: float
    variants_per_second: float


def get_file_size_mb(path: Path) -> float:
    """Get file size in MB."""
    if path.exists():
        return path.stat().st_size / (1024 * 1024)
    return 0.0


def get_memory_usage_mb() -> float:
    """Get current memory usage in MB."""
    try:
        import resource
        return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    except:
        return 0.0


def benchmark_variant_count(source: VariantSource, format_name: str) -> BenchmarkResult:
    """Benchmark counting all variants."""
    gc.collect()
    start_mem = get_memory_usage_mb()

    start_time = time.perf_counter()
    count = source.variant_count
    elapsed = time.perf_counter() - start_time

    end_mem = get_memory_usage_mb()

    return BenchmarkResult(
        format=format_name,
        operation="variant_count",
        file_size_mb=0,  # Set later
        variant_count=count,
        sample_count=source.sample_count,
        het_count=0,
        time_seconds=elapsed,
        memory_mb=end_mem - start_mem,
        variants_per_second=count / elapsed if elapsed > 0 else 0
    )


def benchmark_iter_all(source: VariantSource, format_name: str, sample: str) -> BenchmarkResult:
    """Benchmark iterating over all variants for a sample."""
    gc.collect()
    start_mem = get_memory_usage_mb()

    start_time = time.perf_counter()
    count = 0
    for vg in source.iter_variants(samples=[sample], het_only=False):
        count += 1
    elapsed = time.perf_counter() - start_time

    end_mem = get_memory_usage_mb()

    return BenchmarkResult(
        format=format_name,
        operation="iter_all_variants",
        file_size_mb=0,
        variant_count=count,
        sample_count=1,
        het_count=0,
        time_seconds=elapsed,
        memory_mb=end_mem - start_mem,
        variants_per_second=count / elapsed if elapsed > 0 else 0
    )


def benchmark_iter_het(source: VariantSource, format_name: str, sample: str) -> BenchmarkResult:
    """Benchmark iterating over heterozygous variants only."""
    gc.collect()
    start_mem = get_memory_usage_mb()

    start_time = time.perf_counter()
    het_count = 0
    for vg in source.iter_variants(samples=[sample], het_only=True):
        het_count += 1
    elapsed = time.perf_counter() - start_time

    end_mem = get_memory_usage_mb()

    return BenchmarkResult(
        format=format_name,
        operation="iter_het_only",
        file_size_mb=0,
        variant_count=source.variant_count,
        sample_count=1,
        het_count=het_count,
        time_seconds=elapsed,
        memory_mb=end_mem - start_mem,
        variants_per_second=source.variant_count / elapsed if elapsed > 0 else 0
    )


def benchmark_to_bed(source: VariantSource, format_name: str, sample: str, output_dir: Path) -> BenchmarkResult:
    """Benchmark exporting to BED file (the actual pipeline operation)."""
    gc.collect()
    start_mem = get_memory_usage_mb()

    output_bed = output_dir / f"{format_name}_het.bed"

    start_time = time.perf_counter()
    source.to_bed(output_bed, samples=[sample], het_only=True, include_genotypes=True)
    elapsed = time.perf_counter() - start_time

    end_mem = get_memory_usage_mb()

    # Count output lines
    het_count = sum(1 for _ in open(output_bed))

    return BenchmarkResult(
        format=format_name,
        operation="to_bed_het",
        file_size_mb=0,
        variant_count=source.variant_count,
        sample_count=1,
        het_count=het_count,
        time_seconds=elapsed,
        memory_mb=end_mem - start_mem,
        variants_per_second=source.variant_count / elapsed if elapsed > 0 else 0
    )


def run_benchmarks(
    vcf_path: Optional[Path],
    pgen_path: Optional[Path],
    sample: str,
    n_runs: int = 3,
    output_dir: Path = None
) -> pd.DataFrame:
    """Run all benchmarks and return results DataFrame."""

    if output_dir is None:
        output_dir = Path(tempfile.mkdtemp())

    results = []

    # Benchmark VCF if provided
    if vcf_path and vcf_path.exists():
        print(f"\n{'='*60}")
        print(f"Benchmarking VCF: {vcf_path}")
        print(f"{'='*60}")

        vcf_size = get_file_size_mb(vcf_path)

        for run in range(n_runs):
            print(f"\n  Run {run + 1}/{n_runs}...")

            with VariantSource.open(vcf_path) as source:
                # Check sample exists
                if sample not in source.samples:
                    print(f"  Warning: Sample '{sample}' not in VCF. Using first sample: {source.samples[0]}")
                    sample = source.samples[0]

                # Run benchmarks
                r1 = benchmark_variant_count(source, "VCF")
                r1.file_size_mb = vcf_size
                results.append(r1)
                print(f"    variant_count: {r1.time_seconds:.4f}s")

                r2 = benchmark_iter_all(source, "VCF", sample)
                r2.file_size_mb = vcf_size
                results.append(r2)
                print(f"    iter_all: {r2.time_seconds:.4f}s ({r2.variant_count} variants)")

                r3 = benchmark_iter_het(source, "VCF", sample)
                r3.file_size_mb = vcf_size
                results.append(r3)
                print(f"    iter_het: {r3.time_seconds:.4f}s ({r3.het_count} het sites)")

                r4 = benchmark_to_bed(source, "VCF", sample, output_dir)
                r4.file_size_mb = vcf_size
                results.append(r4)
                print(f"    to_bed: {r4.time_seconds:.4f}s ({r4.het_count} het sites)")

    # Benchmark PGEN if provided
    if pgen_path and pgen_path.exists():
        print(f"\n{'='*60}")
        print(f"Benchmarking PGEN: {pgen_path}")
        print(f"{'='*60}")

        # PGEN file size includes .pgen, .pvar, .psam
        pgen_size = get_file_size_mb(pgen_path)
        pvar_path = pgen_path.with_suffix('.pvar')
        psam_path = pgen_path.with_suffix('.psam')
        if pvar_path.exists():
            pgen_size += get_file_size_mb(pvar_path)
        if psam_path.exists():
            pgen_size += get_file_size_mb(psam_path)

        for run in range(n_runs):
            print(f"\n  Run {run + 1}/{n_runs}...")

            with VariantSource.open(pgen_path) as source:
                # Check sample exists
                if sample not in source.samples:
                    print(f"  Warning: Sample '{sample}' not in PGEN. Using first sample: {source.samples[0]}")
                    sample = source.samples[0]

                # Run benchmarks
                r1 = benchmark_variant_count(source, "PGEN")
                r1.file_size_mb = pgen_size
                results.append(r1)
                print(f"    variant_count: {r1.time_seconds:.4f}s")

                r2 = benchmark_iter_all(source, "PGEN", sample)
                r2.file_size_mb = pgen_size
                results.append(r2)
                print(f"    iter_all: {r2.time_seconds:.4f}s ({r2.variant_count} variants)")

                r3 = benchmark_iter_het(source, "PGEN", sample)
                r3.file_size_mb = pgen_size
                results.append(r3)
                print(f"    iter_het: {r3.time_seconds:.4f}s ({r3.het_count} het sites)")

                r4 = benchmark_to_bed(source, "PGEN", sample, output_dir)
                r4.file_size_mb = pgen_size
                results.append(r4)
                print(f"    to_bed: {r4.time_seconds:.4f}s ({r4.het_count} het sites)")

    # Convert to DataFrame
    df = pd.DataFrame([vars(r) for r in results])

    return df


def create_benchmark_figure(df: pd.DataFrame, output_path: Path):
    """Create publication-quality benchmark figure following Nature guidelines.

    Style specifications based on Nature Genetics requirements:
    - Font: Helvetica/Arial, 5-7pt for labels, 8pt bold for panel labels
    - Colors: Okabe-Ito colorblind-friendly palette
    - Resolution: 300+ dpi
    - Width: max 180mm (single column ~89mm, double ~183mm)
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib as mpl
        from matplotlib.patches import Patch

    except ImportError:
        print("Warning: matplotlib not available, skipping figure generation")
        return

    # ==========================================================================
    # Nature-style configuration
    # ==========================================================================

    # Okabe-Ito colorblind-friendly palette
    COLORS = {
        'VCF': '#E69F00',    # Orange
        'PGEN': '#56B4E9',   # Sky blue
        'accent': '#009E73', # Bluish green
        'gray': '#999999',   # Gray
    }

    # Figure dimensions (Nature: max 180mm width, we use ~170mm for margins)
    # Single column = 89mm, 1.5 column = 120mm, double = 183mm
    FIGURE_WIDTH_MM = 170
    FIGURE_HEIGHT_MM = 70
    MM_TO_INCH = 1 / 25.4

    fig_width = FIGURE_WIDTH_MM * MM_TO_INCH
    fig_height = FIGURE_HEIGHT_MM * MM_TO_INCH

    # Font settings (Nature: Helvetica/Arial, 5-7pt text, 8pt bold panel labels)
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
        'font.size': 7,
        'axes.labelsize': 7,
        'axes.titlesize': 8,
        'axes.titleweight': 'bold',
        'xtick.labelsize': 6,
        'ytick.labelsize': 6,
        'legend.fontsize': 6,
        'axes.linewidth': 0.5,
        'xtick.major.width': 0.5,
        'ytick.major.width': 0.5,
        'xtick.major.size': 2,
        'ytick.major.size': 2,
        'lines.linewidth': 0.75,
        'patch.linewidth': 0.5,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.02,
    })

    # ==========================================================================
    # Data preparation
    # ==========================================================================

    # Aggregate results by format and operation
    summary = df.groupby(['format', 'operation']).agg({
        'time_seconds': ['mean', 'std'],
        'variants_per_second': ['mean', 'std'],
        'variant_count': 'first',
        'het_count': 'first',
        'file_size_mb': 'first'
    }).reset_index()

    summary.columns = ['format', 'operation', 'time_mean', 'time_std',
                       'vps_mean', 'vps_std', 'variant_count', 'het_count', 'file_size_mb']

    # ==========================================================================
    # Create figure
    # ==========================================================================

    fig, axes = plt.subplots(1, 2, figsize=(fig_width, fig_height))
    plt.subplots_adjust(wspace=0.35)

    # Operations to show (pipeline-relevant)
    operations = ['iter_het_only', 'to_bed_het']
    op_labels = ['Iterate\nhet sites', 'Export\nto BED']

    x = np.arange(len(operations))
    width = 0.35

    # --------------------------------------------------------------------------
    # Panel a: Time comparison
    # --------------------------------------------------------------------------
    ax1 = axes[0]

    for i, fmt in enumerate(['VCF', 'PGEN']):
        fmt_data = summary[summary['format'] == fmt]
        times = []
        stds = []
        for op in operations:
            op_data = fmt_data[fmt_data['operation'] == op]
            if len(op_data) > 0:
                times.append(op_data['time_mean'].values[0] * 1000)  # Convert to ms
                stds.append(op_data['time_std'].values[0] * 1000)
            else:
                times.append(0)
                stds.append(0)

        bars = ax1.bar(x + (i - 0.5) * width, times, width,
                       label=fmt, color=COLORS[fmt],
                       edgecolor='black', linewidth=0.5,
                       yerr=stds, capsize=2, error_kw={'linewidth': 0.5})

    ax1.set_ylabel('Time (ms)')
    ax1.set_xticks(x)
    ax1.set_xticklabels(op_labels)
    ax1.set_ylim(bottom=0)
    ax1.legend(frameon=False, loc='upper left')

    # Panel label (Nature style: 8pt bold lowercase)
    ax1.text(-0.15, 1.05, 'a', transform=ax1.transAxes,
             fontsize=8, fontweight='bold', va='bottom')

    # --------------------------------------------------------------------------
    # Panel b: Speedup factor
    # --------------------------------------------------------------------------
    ax2 = axes[1]

    speedups = []
    for op in operations:
        vcf_time = summary[(summary['format'] == 'VCF') & (summary['operation'] == op)]['time_mean'].values
        pgen_time = summary[(summary['format'] == 'PGEN') & (summary['operation'] == op)]['time_mean'].values

        if len(vcf_time) > 0 and len(pgen_time) > 0 and pgen_time[0] > 0:
            speedups.append(vcf_time[0] / pgen_time[0])
        else:
            speedups.append(1.0)

    bars = ax2.bar(x, speedups, width * 1.5,
                   color=COLORS['PGEN'], edgecolor='black', linewidth=0.5)

    # Reference line at 1x
    ax2.axhline(y=1.0, color=COLORS['gray'], linestyle='--', linewidth=0.5, zorder=0)

    ax2.set_ylabel('Fold change (VCF / PGEN)')
    ax2.set_xticks(x)
    ax2.set_xticklabels(op_labels)
    ax2.set_ylim(bottom=0)

    # Add fold-change labels on bars
    for bar, speedup in zip(bars, speedups):
        height = bar.get_height()
        label = f'{speedup:.1f}×' if speedup >= 1 else f'{speedup:.2f}×'
        ax2.annotate(label,
                     xy=(bar.get_x() + bar.get_width() / 2, height),
                     xytext=(0, 2), textcoords='offset points',
                     ha='center', va='bottom', fontsize=6, fontweight='bold')

    # Panel label
    ax2.text(-0.15, 1.05, 'b', transform=ax2.transAxes,
             fontsize=8, fontweight='bold', va='bottom')

    # ==========================================================================
    # Add figure title
    # ==========================================================================

    fig.suptitle('Variant I/O performance: VCF vs PGEN format',
                 fontsize=8, fontweight='bold', y=1.02)

    # ==========================================================================
    # Save figure
    # ==========================================================================

    # Save as PNG (for preview) and PDF (for publication)
    plt.savefig(output_path, format='png', dpi=300, bbox_inches='tight')
    plt.savefig(output_path.with_suffix('.pdf'), format='pdf', bbox_inches='tight')

    # Also save as EPS (Nature preferred vector format)
    try:
        plt.savefig(output_path.with_suffix('.eps'), format='eps', bbox_inches='tight')
        print(f"EPS saved to: {output_path.with_suffix('.eps')}")
    except Exception as e:
        print(f"Note: EPS export skipped ({e})")

    plt.close()

    print(f"\nFigure saved to: {output_path}")
    print(f"PDF saved to: {output_path.with_suffix('.pdf')}")


def print_summary(df: pd.DataFrame):
    """Print a summary of benchmark results."""
    print("\n" + "=" * 70)
    print("BENCHMARK SUMMARY")
    print("=" * 70)

    # Aggregate
    summary = df.groupby(['format', 'operation']).agg({
        'time_seconds': ['mean', 'std'],
        'variant_count': 'first',
        'het_count': 'first',
        'file_size_mb': 'first'
    }).reset_index()

    summary.columns = ['Format', 'Operation', 'Time (s)', 'Std',
                       'Variants', 'Het Sites', 'File Size (MB)']

    print(f"\n{summary.to_string(index=False)}")

    # Calculate speedups
    print("\n" + "-" * 70)
    print("SPEEDUP ANALYSIS (PGEN vs VCF)")
    print("-" * 70)

    for op in df['operation'].unique():
        vcf_time = df[(df['format'] == 'VCF') & (df['operation'] == op)]['time_seconds'].mean()
        pgen_time = df[(df['format'] == 'PGEN') & (df['operation'] == op)]['time_seconds'].mean()

        if pgen_time > 0:
            speedup = vcf_time / pgen_time
            print(f"  {op:20s}: {speedup:6.2f}× faster with PGEN")


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark VCF vs PGEN variant I/O performance",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Benchmark with existing files
  python benchmark_variant_io.py --vcf variants.vcf.gz --pgen variants.pgen --sample NA12878

  # Run more iterations for stable results
  python benchmark_variant_io.py --vcf variants.vcf.gz --pgen variants.pgen --sample NA12878 --runs 5
        """
    )

    parser.add_argument('--vcf', type=Path, help='Path to VCF/VCF.GZ file')
    parser.add_argument('--pgen', type=Path, help='Path to PGEN file')
    parser.add_argument('--sample', type=str, required=True, help='Sample ID to benchmark')
    parser.add_argument('--runs', type=int, default=3, help='Number of benchmark runs (default: 3)')
    parser.add_argument('--output', type=Path, default=Path('benchmark_results'),
                        help='Output directory for results')

    args = parser.parse_args()

    if not args.vcf and not args.pgen:
        parser.error("At least one of --vcf or --pgen must be provided")

    # Create output directory
    args.output.mkdir(parents=True, exist_ok=True)

    # Run benchmarks
    df = run_benchmarks(
        vcf_path=args.vcf,
        pgen_path=args.pgen,
        sample=args.sample,
        n_runs=args.runs,
        output_dir=args.output
    )

    # Save raw results
    csv_path = args.output / 'benchmark_results.csv'
    df.to_csv(csv_path, index=False)
    print(f"\nResults saved to: {csv_path}")

    # Print summary
    print_summary(df)

    # Create figure
    fig_path = args.output / 'benchmark_comparison.png'
    create_benchmark_figure(df, fig_path)


if __name__ == '__main__':
    main()
