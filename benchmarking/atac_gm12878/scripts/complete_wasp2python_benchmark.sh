#!/bin/bash
# Complete WASP2-Python benchmark from Step 2 onwards
# Uses existing Step 1 output (Step 1 time: 2824.846381888s)
#
#$ -N wasp2py_complete
#$ -V
#$ -pe iblm 8
#$ -l h_vmem=48G
#$ -j y
#$ -o /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/atac_gm12878/logs/
#$ -cwd

set -e

# Use existing output directory
OUTPUT_DIR="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/atac_gm12878/results/wasp2python_snp_fixed_2025-12-07_13-12-00"
WASP2_PYTHON_DIR="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-python-bench"
REF_GENOME="/iblm/netapp/data1/aho/ref_genomes/index/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/genome.fa"
THREADS=8
SAMPLE="NA12878"
TIMESTAMP="2025-12-07_13-12-00"

# Step 1 time from original run
STEP1_TIME=2824.846381888

# Conda
source /iblm/netapp/home/jjaureguy/mambaforge/etc/profile.d/conda.sh
conda activate WASP2_dev2

# Tools
BWA=$(which bwa)
SAMTOOLS=$(which samtools)
PYTHON=$(which python)
TIME="/usr/bin/time"

# Set PYTHONPATH
export PYTHONPATH="${WASP2_PYTHON_DIR}/src:${PYTHONPATH}"

cd "${OUTPUT_DIR}"

# Log file
PROFILE_LOG="${OUTPUT_DIR}/profile.log"

echo "========================================"
echo "Completing WASP2-Python benchmark from Step 2"
echo "Using existing Step 1 output"
echo "========================================"
echo "Start time: $(date)"

# Find existing files
REMAP_FQ1="${OUTPUT_DIR}/original_swapped_alleles_r1.fq"
REMAP_FQ2="${OUTPUT_DIR}/original_swapped_alleles_r2.fq"
KEEP_BAM="${OUTPUT_DIR}/original_keep.bam"
TO_REMAP_BAM="${OUTPUT_DIR}/original_to_remap.bam"
WASP_JSON="${OUTPUT_DIR}/original_wasp_data_files.json"

echo "Input files:"
echo "  REMAP_FQ1: ${REMAP_FQ1}"
echo "  REMAP_FQ2: ${REMAP_FQ2}"
echo "  KEEP_BAM: ${KEEP_BAM}"
echo "  WASP_JSON: ${WASP_JSON}"

# Verify files exist
[ -f "${REMAP_FQ1}" ] || { echo "ERROR: ${REMAP_FQ1} not found"; exit 1; }
[ -f "${REMAP_FQ2}" ] || { echo "ERROR: ${REMAP_FQ2} not found"; exit 1; }
[ -f "${KEEP_BAM}" ] || { echo "ERROR: ${KEEP_BAM} not found"; exit 1; }
[ -f "${WASP_JSON}" ] || { echo "ERROR: ${WASP_JSON} not found"; exit 1; }

ORIGINAL_BAM="${OUTPUT_DIR}/original.bam"
ORIGINAL_READS=$(${SAMTOOLS} view -c ${ORIGINAL_BAM})
echo "Original reads: ${ORIGINAL_READS}"

KEEP_READS=$(${SAMTOOLS} view -c ${KEEP_BAM})
echo "Keep reads: ${KEEP_READS}"

# Count reads to remap
REMAP_LINES=$(wc -l < ${REMAP_FQ1})
REMAP_PAIRS=$((REMAP_LINES / 4))
echo "Reads to remap: ${REMAP_PAIRS} pairs"

# -----------------------------------------------------------------------------
# STEP 2: BWA-MEM remap flipped reads
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 2: BWA-MEM remap flipped reads"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 2: BWA-MEM Remap (completion run)" >> ${PROFILE_LOG}

STEP2_START=$(date +%s.%N)

${TIME} -v ${BWA} mem -t ${THREADS} -M ${REF_GENOME} \
    ${REMAP_FQ1} ${REMAP_FQ2} 2>> ${PROFILE_LOG} | \
    ${SAMTOOLS} view -S -b -h -F 4 - > ${OUTPUT_DIR}/remapped_unsorted.bam

${SAMTOOLS} sort -@ ${THREADS} -o ${OUTPUT_DIR}/remapped.bam ${OUTPUT_DIR}/remapped_unsorted.bam
${SAMTOOLS} index ${OUTPUT_DIR}/remapped.bam
rm -f ${OUTPUT_DIR}/remapped_unsorted.bam

STEP2_END=$(date +%s.%N)
STEP2_TIME=$(echo "${STEP2_END} - ${STEP2_START}" | bc)
echo "STEP 2 completed in ${STEP2_TIME} seconds"

# -----------------------------------------------------------------------------
# STEP 3: WASP2-Python filter-remapped
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 3: WASP2-Python filter-remapped"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 3: WASP2-Python filter-remapped" >> ${PROFILE_LOG}

STEP3_START=$(date +%s.%N)

cd "${WASP2_PYTHON_DIR}/src/mapping"

${TIME} -v ${PYTHON} __main__.py filter-remapped \
    "${OUTPUT_DIR}/remapped.bam" \
    --wasp_data_json "${WASP_JSON}" \
    --out_bam "${OUTPUT_DIR}/wasp_filtered.bam" \
    --remap_keep_bam "${OUTPUT_DIR}/remap_keep.bam" 2>> ${PROFILE_LOG}

cd "${OUTPUT_DIR}"

${SAMTOOLS} index ${OUTPUT_DIR}/wasp_filtered.bam 2>/dev/null || true
${SAMTOOLS} index ${OUTPUT_DIR}/remap_keep.bam 2>/dev/null || true

STEP3_END=$(date +%s.%N)
STEP3_TIME=$(echo "${STEP3_END} - ${STEP3_START}" | bc)
echo "STEP 3 completed in ${STEP3_TIME} seconds"

# -----------------------------------------------------------------------------
# SUMMARY
# -----------------------------------------------------------------------------
TOTAL_TIME=$(echo "${STEP1_TIME} + ${STEP2_TIME} + ${STEP3_TIME}" | bc)

# Count final reads
FINAL_READS=$(${SAMTOOLS} view -c ${OUTPUT_DIR}/wasp_filtered.bam 2>/dev/null || echo "0")
REMAP_KEEP_READS=$(${SAMTOOLS} view -c ${OUTPUT_DIR}/remap_keep.bam 2>/dev/null || echo "0")
if [ "${ORIGINAL_READS}" -gt 0 ]; then
    PASS_RATE=$(echo "scale=2; ${FINAL_READS} * 100 / ${ORIGINAL_READS}" | bc)
else
    PASS_RATE="0"
fi

# WASP-only time (steps 1+3, excluding BWA alignment)
WASP_ONLY=$(echo "${STEP1_TIME} + ${STEP3_TIME}" | bc)

# Variant count
VARIANT_COUNT=$(wc -l < ${OUTPUT_DIR}/variants.bed)

echo "========================================"
echo "BENCHMARK RESULTS (WASP2-Python FIXED)"
echo "========================================"
echo ""
echo "Individual step times:"
echo "  STEP 1 (make-reads):        ${STEP1_TIME} s"
echo "  STEP 2 (BWA-MEM remap):     ${STEP2_TIME} s"
echo "  STEP 3 (filter-remapped):   ${STEP3_TIME} s"
echo ""
echo "WASP-only time (steps 1+3):   ${WASP_ONLY} s"
echo "TOTAL time:                   ${TOTAL_TIME} s"
echo ""
echo "Read counts:"
echo "  Original reads:       ${ORIGINAL_READS}"
echo "  Keep reads:           ${KEEP_READS}"
echo "  Remap keep reads:     ${REMAP_KEEP_READS}"
echo "  FINAL (filtered):     ${FINAL_READS}"
echo "  Pass rate:            ${PASS_RATE}%"
echo ""
echo "Variant counts (het-only):"
echo "  Het SNPs:             ${VARIANT_COUNT}"
echo ""
echo "Output directory: ${OUTPUT_DIR}"
echo "End time: $(date)"

# Save results to JSON
cat > ${OUTPUT_DIR}/benchmark_results.json << EOF
{
    "timestamp": "${TIMESTAMP}",
    "pipeline": "wasp2python_snp_FIXED",
    "sample": "${SAMPLE}",
    "data_type": "ATAC-seq",
    "threads": ${THREADS},
    "het_only": true,
    "step1_make_reads_s": ${STEP1_TIME},
    "step2_bwa_remap_s": ${STEP2_TIME},
    "step3_filter_remapped_s": ${STEP3_TIME},
    "wasp_only_s": ${WASP_ONLY},
    "total_s": ${TOTAL_TIME},
    "original_reads": ${ORIGINAL_READS},
    "keep_reads": ${KEEP_READS},
    "remap_keep_reads": ${REMAP_KEEP_READS},
    "final_reads": ${FINAL_READS},
    "pass_rate_percent": ${PASS_RATE},
    "snv_count": ${VARIANT_COUNT}
}
EOF

echo ""
echo "Results saved to: ${OUTPUT_DIR}/benchmark_results.json"
