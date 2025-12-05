#!/bin/bash
#$ -N nature_methods_bench
#$ -V
#$ -pe iblm 8
#$ -l h_vmem=16G
#$ -j y
#$ -o /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/logs/
#$ -cwd

# =============================================================================
# SGE submission for Nature Methods Benchmark
# WASP2 vs GATK vs phASER vs biastools
# 8 threads, 16GB memory
# =============================================================================

set -e

source /iblm/netapp/home/jjaureguy/mambaforge/etc/profile.d/conda.sh
conda activate WASP2_dev2

cd /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp

./benchmarking/nature_methods_benchmark.sh
