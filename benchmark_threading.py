#!/usr/bin/env python3
"""
Benchmark script to demonstrate threading optimization in BWA alignment.

This script runs BWA alignment with different thread counts to show
the performance improvement from dynamic threading.
"""

import subprocess
import tempfile
import time
from pathlib import Path
import multiprocessing

def create_test_data(workdir: Path, num_reads: int = 5000):
    """Create minimal test data for benchmarking."""
    ref_fasta = workdir / "ref.fa"
    fastq_file = workdir / "reads.fq"

    # Create simple reference
    with open(ref_fasta, 'w') as f:
        f.write(">chr1\n")
        f.write("ATCG" * 250 + "\n")  # 1kb reference

    # Index reference
    subprocess.run(['samtools', 'faidx', str(ref_fasta)],
                   check=True, capture_output=True)
    subprocess.run(['bwa', 'index', str(ref_fasta)],
                   check=True, capture_output=True)

    # Create synthetic reads
    with open(fastq_file, 'w') as f:
        for i in range(num_reads):
            f.write(f"@read_{i}\n")
            f.write("ATCGATCGATCGATCGATCGATCGATCGATCG\n")  # 32bp read
            f.write("+\n")
            f.write("I" * 32 + "\n")  # Quality scores

    return ref_fasta, fastq_file

def benchmark_bwa(ref_fasta: Path, fastq_file: Path, output_bam: Path, threads: int):
    """Run BWA alignment and measure time."""
    start = time.time()

    # Run BWA
    sam_file = output_bam.with_suffix('.sam')
    with open(sam_file, 'w') as sam:
        subprocess.run([
            'bwa', 'mem',
            '-t', str(threads),
            str(ref_fasta),
            str(fastq_file)
        ], stdout=sam, stderr=subprocess.PIPE, check=True)

    elapsed = time.time() - start

    # Cleanup
    sam_file.unlink()

    return elapsed

def main():
    print("=" * 80)
    print("BWA THREADING BENCHMARK")
    print("=" * 80)
    print()

    # Get system info
    cpu_count = multiprocessing.cpu_count()
    optimal_threads = min(cpu_count, 16)

    print(f"System CPU count: {cpu_count}")
    print(f"Optimal thread count (capped at 16): {optimal_threads}")
    print()

    # Create test data
    print("Creating test data...")
    with tempfile.TemporaryDirectory(prefix='bwa_bench_') as tmpdir:
        workdir = Path(tmpdir)
        ref_fasta, fastq_file = create_test_data(workdir, num_reads=10000)
        print(f"  Created test files in {workdir}")
        print()

        # Benchmark different thread counts
        thread_counts = [1, 4, 8, optimal_threads]
        results = []

        print("Running benchmarks...")
        print("-" * 80)
        print(f"{'Threads':<10} {'Time (s)':<15} {'Speedup':<15}")
        print("-" * 80)

        baseline_time = None
        for threads in thread_counts:
            if threads > cpu_count:
                continue

            output_bam = workdir / f"test_{threads}t.bam"
            elapsed = benchmark_bwa(ref_fasta, fastq_file, output_bam, threads)

            if baseline_time is None:
                baseline_time = elapsed
                speedup = 1.0
            else:
                speedup = baseline_time / elapsed

            results.append((threads, elapsed, speedup))
            print(f"{threads:<10} {elapsed:<15.3f} {speedup:<15.2f}x")

        print("-" * 80)
        print()

        # Summary
        optimal_result = [r for r in results if r[0] == optimal_threads][0]
        old_result = [r for r in results if r[0] == 4][0]

        print("SUMMARY:")
        print(f"  Old default (4 threads):     {old_result[1]:.3f}s")
        print(f"  New optimized ({optimal_threads} threads): {optimal_result[1]:.3f}s")
        print(f"  Speedup:                     {old_result[1]/optimal_result[1]:.2f}x faster")
        print()

        print("=" * 80)
        if optimal_result[2] >= 2.0:
            print("✅ OPTIMIZATION SUCCESSFUL: Achieved 2-3x speedup!")
        else:
            print(f"⚠️  Speedup: {optimal_result[2]:.2f}x (goal: 2-3x)")
        print("=" * 80)

if __name__ == '__main__':
    main()
