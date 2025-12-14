#!/bin/bash
# =============================================================================
# WASP2 UNIFIED PIPELINE BENCHMARK - 15M reads (direct run, no SGE)
# Uses the fast single-pass unified pipeline (Rust throughout, no fallbacks)
# =============================================================================

set -e

source /iblm/netapp/home/jjaureguy/mambaforge/etc/profile.d/conda.sh
conda activate WASP2_dev2

# Timestamp for this run
TIMESTAMP=$(date +%Y-%m-%d_%H-%M-%S)
echo "=============================================="
echo "WASP2 UNIFIED - 15M reads benchmark"
echo "Timestamp: ${TIMESTAMP}"
echo "=============================================="

WASP2_DIR="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp"

# Timestamped output directory
OUTPUT_DIR="${WASP2_DIR}/benchmarking/results/unified_15m_${TIMESTAMP}"
mkdir -p "${OUTPUT_DIR}"
echo "Output directory: ${OUTPUT_DIR}"

# Inputs
genome_index="/iblm/netapp/data1/aho/ref_genomes/index/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/genome.fa"
input_bam="/iblm/netapp/data1/aho/atac/Buenrostro2013/merged/GM12878_ATACseq_50k/GM12878_ATACseq_50k_merged.sorted.bam"
input_vcf="/iblm/netapp/data1/aho/variants/NA12878.vcf.gz"
sample="NA12878"
THREADS="${THREADS:-8}"
COMPRESSION_THREADS="${COMPRESSION_THREADS:-1}"

# 15M reads
n_subset=15000000
bam_reads=159133307
bam_prop=$(echo "scale=7; ${n_subset}/${bam_reads}" | bc)
seed=$RANDOM

# Verify Rust module
echo ""
echo "Verifying Rust module..."
python -c "import wasp2_rust; print('wasp2_rust module: AVAILABLE')" || {
    echo "ERROR: wasp2_rust module not available"
    exit 1
}
python -c "from wasp2_rust import unified_make_reads_parallel_py; print('Unified pipeline: AVAILABLE')" || {
    echo "ERROR: Unified pipeline not available"
    exit 1
}
echo ""

# Working directory
cd "${OUTPUT_DIR}"

# Subset BAM
echo "Subsetting BAM to ${n_subset} reads..."
subset_bam="${OUTPUT_DIR}/subset_15m.bam"
samtools view -h --bam --subsample ${bam_prop} --subsample-seed ${seed} -o ${subset_bam} ${input_bam}
samtools index ${subset_bam}

# Total time tracking
export PYTHONPATH="${WASP2_DIR}/src:${PYTHONPATH}"
TOTAL_START=$(date +%s.%N)

# -----------------------------------------------------------------------------
# STEP 1: WASP2 unified make-reads
# -----------------------------------------------------------------------------
echo ""
echo "STEP 1: WASP2 unified make-reads..."
STEP1_START=$(date +%s.%N)

python -c "
from mapping.run_mapping import run_make_remap_reads_unified
import json

stats = run_make_remap_reads_unified(
    bam_file='${subset_bam}',
    variant_file='${input_vcf}',
    samples='${sample}',
    out_dir='${OUTPUT_DIR}',
    threads=${THREADS},
    compression_threads=${COMPRESSION_THREADS},
    use_parallel=True
)

with open('${OUTPUT_DIR}/unified_stats.json', 'w') as f:
    json.dump(stats, f, indent=2)

print(f'Pairs processed: {stats[\"pairs_processed\"]:,}')
print(f'Haplotypes written: {stats[\"haplotypes_written\"]:,}')
"

STEP1_END=$(date +%s.%N)
STEP1_TIME=$(echo "${STEP1_END} - ${STEP1_START}" | bc)
echo "STEP 1 completed in ${STEP1_TIME} seconds"

# Find FASTQ files
r1_reads=$(ls ${OUTPUT_DIR}/*_remap_r1.fq.gz 2>/dev/null | head -1)
r2_reads=$(ls ${OUTPUT_DIR}/*_remap_r2.fq.gz 2>/dev/null | head -1)

# -----------------------------------------------------------------------------
# STEP 2: BWA remap
# -----------------------------------------------------------------------------
echo ""
echo "STEP 2: BWA remap..."
STEP2_START=$(date +%s.%N)

zcat ${r1_reads} > ${OUTPUT_DIR}/remap_r1.fq
zcat ${r2_reads} > ${OUTPUT_DIR}/remap_r2.fq

bwa mem -t 16 -M ${genome_index} ${OUTPUT_DIR}/remap_r1.fq ${OUTPUT_DIR}/remap_r2.fq | \
    samtools view -S -b -h -F 4 - > ${OUTPUT_DIR}/remapped.bam
samtools sort -o ${OUTPUT_DIR}/remapped_sorted.bam ${OUTPUT_DIR}/remapped.bam
mv ${OUTPUT_DIR}/remapped_sorted.bam ${OUTPUT_DIR}/remapped.bam
samtools index ${OUTPUT_DIR}/remapped.bam

rm -f ${OUTPUT_DIR}/remap_r1.fq ${OUTPUT_DIR}/remap_r2.fq

STEP2_END=$(date +%s.%N)
STEP2_TIME=$(echo "${STEP2_END} - ${STEP2_START}" | bc)
echo "STEP 2 completed in ${STEP2_TIME} seconds"

# -----------------------------------------------------------------------------
# STEP 3: WASP2 filter-remapped
# -----------------------------------------------------------------------------
echo ""
echo "STEP 3: WASP2 filter-remapped..."
STEP3_START=$(date +%s.%N)

python -c "
from wasp2_rust import filter_bam_wasp

kept, removed_moved, removed_missing = filter_bam_wasp(
    '${subset_bam}',
    '${OUTPUT_DIR}/remapped.bam',
    '${OUTPUT_DIR}/remap_keep.bam',
    threads=8
)

print(f'Kept: {kept}, Removed: {removed_moved}, Missing: {removed_missing}')
"

samtools index ${OUTPUT_DIR}/remap_keep.bam

STEP3_END=$(date +%s.%N)
STEP3_TIME=$(echo "${STEP3_END} - ${STEP3_START}" | bc)
echo "STEP 3 completed in ${STEP3_TIME} seconds"

# -----------------------------------------------------------------------------
# SUMMARY
# -----------------------------------------------------------------------------
TOTAL_END=$(date +%s.%N)
TOTAL_TIME=$(echo "${TOTAL_END} - ${TOTAL_START}" | bc)
WASP_ONLY=$(echo "${STEP1_TIME} + ${STEP3_TIME}" | bc)

echo ""
echo "=============================================="
echo "BENCHMARK RESULTS - 15M reads (UNIFIED)"
echo "=============================================="
echo "  Unified make-reads: ${STEP1_TIME} s"
echo "  BWA remap:          ${STEP2_TIME} s"
echo "  Filter:             ${STEP3_TIME} s"
echo "  WASP-only (1+3):    ${WASP_ONLY} s"
echo "  TOTAL:              ${TOTAL_TIME} s"
echo "=============================================="
echo "Output: ${OUTPUT_DIR}"

# Save results
cat > ${OUTPUT_DIR}/benchmark_results.json << EOF
{
    "timestamp": "${TIMESTAMP}",
    "pipeline": "unified",
    "n_reads": ${n_subset},
    "seed": ${seed},
    "step1_unified_s": ${STEP1_TIME},
    "step2_bwa_remap_s": ${STEP2_TIME},
    "step3_filter_s": ${STEP3_TIME},
    "wasp_only_s": ${WASP_ONLY},
    "total_s": ${TOTAL_TIME}
}
EOF

echo "Results saved to: ${OUTPUT_DIR}/benchmark_results.json"
