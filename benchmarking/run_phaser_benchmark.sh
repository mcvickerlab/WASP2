#!/bin/bash
#$ -N phaser_bench
#$ -cwd
#$ -q iblm.q
#$ -pe iblm 8
#$ -l h_rt=04:00:00
#$ -l h_vmem=4G
#$ -o logs/phaser_benchmark.o$JOB_ID
#$ -e logs/phaser_benchmark.e$JOB_ID

# phASER Benchmark with 8 threads
# Fair comparison against WASP2-Rust

set -e

# Paths
REPO_ROOT="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp"
PHASER_DIR="${REPO_ROOT}/benchmarking/phaser_tool/phaser"
DATA_DIR="${REPO_ROOT}/benchmarking/star_wasp_comparison/data"
RESULTS_DIR="${REPO_ROOT}/benchmarking/star_wasp_comparison/results/unified_2025-12-04_00-29-39"
OUTPUT_DIR="${REPO_ROOT}/benchmarking/results/phaser_8threads"

BAM="${RESULTS_DIR}/A_sorted.bam"
VCF="${DATA_DIR}/HG00731_het_only_chr.vcf.gz"
THREADS=8

# Activate environment
source ~/.bashrc
conda activate WASP2_dev2

echo "=============================================="
echo "phASER Benchmark (${THREADS} threads)"
echo "=============================================="
echo "BAM: ${BAM}"
echo "VCF: ${VCF}"
echo "Threads: ${THREADS}"
echo "Start: $(date)"
echo ""

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Check inputs
if [[ ! -f "${BAM}" ]]; then
    echo "ERROR: BAM not found: ${BAM}"
    exit 1
fi
if [[ ! -f "${VCF}" ]]; then
    echo "ERROR: VCF not found: ${VCF}"
    exit 1
fi

# Run phASER
echo "Running phASER..."
START=$(date +%s.%N)

cd ${PHASER_DIR}
python3 phaser.py \
    --vcf ${VCF} \
    --bam ${BAM} \
    --sample HG00731 \
    --mapq 10 \
    --baseq 20 \
    --paired_end 1 \
    --threads ${THREADS} \
    --o ${OUTPUT_DIR}/phaser_result \
    2>&1 | tee ${OUTPUT_DIR}/phaser.log

END=$(date +%s.%N)
ELAPSED=$(echo "${END} - ${START}" | bc)

echo ""
echo "=============================================="
echo "RESULTS"
echo "=============================================="
echo "phASER time: ${ELAPSED}s"
echo "Threads used: ${THREADS}"
echo "End: $(date)"

# Save timing
echo "phaser_8threads,${ELAPSED}" > ${OUTPUT_DIR}/timing.csv

echo ""
echo "Output files:"
ls -la ${OUTPUT_DIR}/
