#!/bin/bash
#$ -N wasp2_starwasp_unified
#$ -V
#$ -pe iblm 8
#$ -l h_vmem=16G
#$ -j y
#$ -o /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/logs/
#$ -cwd

# =============================================================================
# WASP2 UNIFIED PIPELINE - Star-wasp comparison benchmark (56M reads)
# Full Rust pipeline: vcf_to_bed + unified_make_reads + filter_bam_wasp
# =============================================================================

set -e

source /iblm/netapp/home/jjaureguy/mambaforge/etc/profile.d/conda.sh
conda activate WASP2_dev2

# Run the benchmark script
cd /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp
./benchmarking/star_wasp_comparison/scripts/run_wasp2_unified_benchmark.sh
