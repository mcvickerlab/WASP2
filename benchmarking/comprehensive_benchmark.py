#!/usr/bin/env python3
"""
Comprehensive WASP2 Benchmark: Accuracy + Speed

Benchmarks WASP2-Rust against GATK, phASER, and biastools.

Part 1: ACCURACY (Simulation)
- Ground truth allele counts vs measured
- R², RMSE, MAE, bias metrics
- SNP vs INDEL performance

Part 2: SPEED (Real Data - HG00731)
- Wall clock time
- Reads per second
- Memory usage
"""

import subprocess
import time
import os
import sys
import json
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
from scipy import stats

# Paths
REPO_ROOT = Path("/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp")
BENCHMARK_DIR = REPO_ROOT / "benchmarking"
SIM_DIR = REPO_ROOT / "simulation_results" / "benchmark_v3"
REAL_DATA_DIR = BENCHMARK_DIR / "star_wasp_comparison" / "data"
RESULTS_DIR = BENCHMARK_DIR / "comprehensive_results"

# Real data paths (HG00731 from STAR-WASP benchmark)
HG00731_R1 = REAL_DATA_DIR / "ERR1050079_1.fastq.gz"
HG00731_R2 = REAL_DATA_DIR / "ERR1050079_2.fastq.gz"
HG00731_VCF = REAL_DATA_DIR / "HG00731_het_only_chr.vcf.gz"


def run_command(cmd, description, timeout=3600):
    """Run command and return timing."""
    print(f"\n{'='*60}")
    print(f"Running: {description}")
    print(f"Command: {' '.join(cmd) if isinstance(cmd, list) else cmd}")
    print(f"{'='*60}")

    start = time.time()
    try:
        result = subprocess.run(
            cmd, shell=isinstance(cmd, str),
            capture_output=True, text=True, timeout=timeout
        )
        elapsed = time.time() - start

        if result.returncode != 0:
            print(f"ERROR: {result.stderr[:500]}")
            return None, elapsed

        print(f"Completed in {elapsed:.2f}s")
        return result.stdout, elapsed
    except subprocess.TimeoutExpired:
        print(f"TIMEOUT after {timeout}s")
        return None, timeout


class AccuracyBenchmark:
    """Benchmark accuracy using simulation with ground truth."""

    def __init__(self, sim_dir: Path, output_dir: Path):
        self.sim_dir = sim_dir
        self.output_dir = output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Load ground truth
        self.ground_truth = pd.read_csv(sim_dir / "ground_truth.csv")
        print(f"Loaded {len(self.ground_truth)} ground truth records")

    def run_wasp2_rust(self):
        """Run WASP2 Rust implementation."""
        print("\n" + "="*60)
        print("Running WASP2-Rust")
        print("="*60)

        bam = self.sim_dir / "aligned.sorted.bam"
        vcf = self.sim_dir / "variants.vcf.gz"
        output = self.output_dir / "wasp2_output.bam"

        # Use the Rust bam_intersect module
        cmd = f"""
        cd {REPO_ROOT} && source ~/.bashrc && conda activate WASP2_dev2 && \
        python -c "
import sys
sys.path.insert(0, 'src')
from wasp2_rust import bam_intersect_variants

result = bam_intersect_variants(
    '{bam}',
    '{vcf}',
    '{output}',
    sample_name='SIMULATED',
    is_paired=True,
    threads=8
)
print(f'Processed {{result}} reads')
"
        """

        start = time.time()
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        elapsed = time.time() - start

        if result.returncode == 0:
            print(f"WASP2-Rust completed in {elapsed:.2f}s")
            print(result.stdout)
        else:
            print(f"WASP2-Rust error: {result.stderr[:500]}")

        return elapsed

    def count_alleles_from_reads(self, bam_path: Path) -> pd.DataFrame:
        """
        Count alleles by parsing read names (ground truth encoded).

        Read names format: @chr1_50000_0001_HAP1/1
        HAP1 = REF allele, HAP2 = ALT allele
        """
        import pysam

        counts = {}
        bam = pysam.AlignmentFile(str(bam_path), 'rb')

        for read in bam:
            if read.is_secondary or read.is_supplementary:
                continue

            # Parse read name: chr1_50000_0001_HAP1
            name = read.query_name
            parts = name.split('_')
            if len(parts) >= 4:
                chrom = parts[0]
                pos = int(parts[1])
                haplotype = parts[3]  # HAP1 or HAP2

                key = (chrom, pos)
                if key not in counts:
                    counts[key] = {'ref': 0, 'alt': 0}

                if haplotype == 'HAP1':
                    counts[key]['ref'] += 1
                elif haplotype == 'HAP2':
                    counts[key]['alt'] += 1

        bam.close()

        # Convert to DataFrame
        rows = []
        for (chrom, pos), c in counts.items():
            total = c['ref'] + c['alt']
            ratio = c['ref'] / total if total > 0 else 0.5
            rows.append({
                'chrom': chrom,
                'pos': pos,
                'wasp2_ref': c['ref'],
                'wasp2_alt': c['alt'],
                'wasp2_total': total,
                'wasp2_ratio': ratio
            })

        return pd.DataFrame(rows)

    def run_gatk(self):
        """Run GATK ASEReadCounter."""
        print("\n" + "="*60)
        print("Running GATK ASEReadCounter")
        print("="*60)

        cmd = f"""
        cd {REPO_ROOT} && source ~/.bashrc && conda activate WASP2_dev2 && \
        python simulation/competitors/run_gatk.py \
            --bam {self.sim_dir}/aligned.sorted.bam \
            --vcf {self.sim_dir}/variants.vcf.gz \
            --ref {self.sim_dir}/reference.fa \
            --output {self.output_dir}/gatk/
        """

        start = time.time()
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        elapsed = time.time() - start

        print(result.stdout)
        if result.returncode != 0:
            print(f"GATK error: {result.stderr[:500]}")

        return elapsed

    def run_phaser(self):
        """Run phASER."""
        print("\n" + "="*60)
        print("Running phASER")
        print("="*60)

        cmd = f"""
        cd {REPO_ROOT} && source ~/.bashrc && conda activate WASP2_dev2 && \
        python simulation/competitors/run_phaser.py \
            --bam {self.sim_dir}/aligned.sorted.bam \
            --vcf {self.sim_dir}/variants.vcf.gz \
            --output {self.output_dir}/phaser/
        """

        start = time.time()
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        elapsed = time.time() - start

        print(result.stdout)
        if result.returncode != 0:
            print(f"phASER error: {result.stderr[:500]}")

        return elapsed

    def run_biastools(self):
        """Run biastools."""
        print("\n" + "="*60)
        print("Running biastools")
        print("="*60)

        cmd = f"""
        cd {REPO_ROOT} && source ~/.bashrc && conda activate WASP2_dev2 && \
        python simulation/competitors/run_biastools.py \
            --bam {self.sim_dir}/aligned.sorted.bam \
            --vcf {self.sim_dir}/variants.vcf.gz \
            --output {self.output_dir}/biastools/
        """

        start = time.time()
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        elapsed = time.time() - start

        print(result.stdout)
        if result.returncode != 0:
            print(f"biastools error: {result.stderr[:500]}")

        return elapsed

    def calculate_metrics(self, true_ratios, measured_ratios) -> dict:
        """Calculate accuracy metrics."""
        # Remove NaN
        mask = ~(np.isnan(true_ratios) | np.isnan(measured_ratios))
        true = true_ratios[mask]
        measured = measured_ratios[mask]

        if len(true) < 2:
            return {'n': 0, 'r': np.nan, 'r2': np.nan, 'rmse': np.nan, 'mae': np.nan}

        r, p = stats.pearsonr(true, measured)
        rmse = np.sqrt(np.mean((true - measured) ** 2))
        mae = np.mean(np.abs(true - measured))

        return {
            'n': len(true),
            'r': r,
            'r2': r ** 2,
            'rmse': rmse,
            'mae': mae,
            'mean_bias': np.mean(measured) - 0.5
        }

    def run_all(self):
        """Run complete accuracy benchmark."""
        results = {'timing': {}, 'accuracy': {}}

        # Run all tools
        results['timing']['wasp2'] = self.run_wasp2_rust()
        results['timing']['gatk'] = self.run_gatk()
        results['timing']['phaser'] = self.run_phaser()
        results['timing']['biastools'] = self.run_biastools()

        # TODO: Load results and calculate accuracy metrics
        # This requires the tools to complete successfully

        return results


class SpeedBenchmark:
    """Benchmark speed using real HG00731 data."""

    def __init__(self, output_dir: Path):
        self.output_dir = output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def get_read_count(self, fastq_path: Path) -> int:
        """Count reads in FASTQ file."""
        cmd = f"zcat {fastq_path} | wc -l"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = int(result.stdout.strip())
        return lines // 4  # 4 lines per read

    def run_wasp2_rust(self, bam_path: Path, vcf_path: Path, threads: int = 8):
        """Run WASP2-Rust on real data."""
        print("\n" + "="*60)
        print(f"Running WASP2-Rust ({threads} threads)")
        print("="*60)

        output = self.output_dir / "wasp2_hg00731.bam"

        cmd = f"""
        cd {REPO_ROOT} && source ~/.bashrc && conda activate WASP2_dev2 && \
        python -c "
import sys
sys.path.insert(0, 'src')
from wasp2_rust import bam_intersect_variants

result = bam_intersect_variants(
    '{bam_path}',
    '{vcf_path}',
    '{output}',
    sample_name='HG00731',
    is_paired=True,
    threads={threads}
)
print(f'Processed {{result}} reads')
"
        """

        start = time.time()
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        elapsed = time.time() - start

        if result.returncode == 0:
            print(f"WASP2-Rust completed in {elapsed:.2f}s")
            print(result.stdout)
        else:
            print(f"Error: {result.stderr[:500]}")

        return elapsed

    def run_scaling_test(self, bam_path: Path, vcf_path: Path):
        """Test scaling across thread counts."""
        thread_counts = [1, 2, 4, 8, 16]
        results = []

        for threads in thread_counts:
            elapsed = self.run_wasp2_rust(bam_path, vcf_path, threads)
            results.append({
                'threads': threads,
                'time_seconds': elapsed,
                'speedup': results[0]['time_seconds'] / elapsed if results else 1.0
            })

        return pd.DataFrame(results)


def main():
    """Run comprehensive benchmark."""
    print("="*60)
    print("WASP2 Comprehensive Benchmark")
    print(f"Started: {datetime.now().isoformat()}")
    print("="*60)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # Check simulation data exists
    if not (SIM_DIR / "aligned.sorted.bam").exists():
        print(f"ERROR: Simulation data not found at {SIM_DIR}")
        print("Run the simulator first:")
        print("  python simulation/ase_simulator_v3.py --output simulation_results/benchmark_v3 ...")
        sys.exit(1)

    # Part 1: Accuracy Benchmark
    print("\n" + "#"*60)
    print("# PART 1: ACCURACY BENCHMARK (Simulation)")
    print("#"*60)

    accuracy_bench = AccuracyBenchmark(SIM_DIR, RESULTS_DIR / "accuracy")
    accuracy_results = accuracy_bench.run_all()

    # Part 2: Speed Benchmark (if real data available)
    if HG00731_R1.exists():
        print("\n" + "#"*60)
        print("# PART 2: SPEED BENCHMARK (HG00731)")
        print("#"*60)

        speed_bench = SpeedBenchmark(RESULTS_DIR / "speed")

        # Need aligned BAM for speed test
        aligned_bam = BENCHMARK_DIR / "star_wasp_comparison" / "results" / "aligned.sorted.bam"
        if aligned_bam.exists():
            speed_results = speed_bench.run_scaling_test(aligned_bam, HG00731_VCF)
            speed_results.to_csv(RESULTS_DIR / "speed" / "scaling_results.csv", index=False)
            print("\nScaling Results:")
            print(speed_results.to_string())
        else:
            print(f"Aligned BAM not found at {aligned_bam}")
            print("Run STAR alignment first")
    else:
        print(f"\nReal data not found at {HG00731_R1}")
        print("Skipping speed benchmark")

    # Save results
    with open(RESULTS_DIR / "benchmark_results.json", 'w') as f:
        json.dump({
            'timestamp': datetime.now().isoformat(),
            'accuracy': accuracy_results
        }, f, indent=2, default=str)

    print("\n" + "="*60)
    print(f"Benchmark complete. Results saved to {RESULTS_DIR}")
    print("="*60)


if __name__ == '__main__':
    main()
