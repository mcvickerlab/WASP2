#!/bin/bash
#$ -N wasp2_hg00731
#$ -pe iblm 8
#$ -l h_vmem=40G
#$ -l h_rt=12:00:00
#$ -cwd
#$ -o logs/wasp2_benchmark_$JOB_ID.out
#$ -e logs/wasp2_benchmark_$JOB_ID.err
#$ -j n
#$ -V

# WASP2-Rust RNA-seq benchmark on HG00731
# Matches STAR+WASP paper methodology exactly:
# - 8 threads
# - Same STAR parameters
# - /usr/bin/time -v memory profiling

echo "============================================"
echo "WASP2-Rust Benchmark Job"
echo "Job ID: ${JOB_ID}"
echo "Host: $(hostname)"
echo "Start: $(date)"
echo "============================================"

# Run the benchmark
bash scripts/run_wasp2_benchmark.sh

echo ""
echo "============================================"
echo "Job complete: $(date)"
echo "============================================"
