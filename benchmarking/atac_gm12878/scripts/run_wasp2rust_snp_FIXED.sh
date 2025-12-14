#!/bin/bash
# WASP2-Rust FIXED v2 benchmark on GM12878 ATAC-seq data - SNP only mode
# CORRECT pipeline: filter → unified → remap → filter → extract_keep_no_flip → MERGE
#
# This produces EQUIVALENT output to WASP1/Python:
# - keep.bam: reads that don't overlap any variants
# - keep_no_flip.bam: reads that overlap variants but already have ref allele (no flip needed)
# - remap_keep.bam: reads that needed allele flipping and passed remapping
#
# Final output = keep + keep_no_flip + remap_keep
#
#$ -N wasp2rust_gm12878_snp_fixed
#$ -V
#$ -pe iblm 8
#$ -l h_vmem=32G
#$ -j y
#$ -o /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/atac_gm12878/logs/
#$ -cwd

set -e

# Timestamp for this run
TIMESTAMP=$(date +%Y-%m-%d_%H-%M-%S)
echo "Benchmark timestamp: ${TIMESTAMP}"

# Paths
WASP2_DIR="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp"
BENCHMARK_DIR="${WASP2_DIR}/benchmarking/atac_gm12878"
FINAL_OUTPUT_DIR="${BENCHMARK_DIR}/results/wasp2rust_snp_fixed_${TIMESTAMP}"
LOCAL_SCRATCH="${TMPDIR:-/tmp}"
WORK_OUTPUT_DIR="${LOCAL_SCRATCH}/wasp2rust_gm12878_snp_${JOB_ID:-local}_${TIMESTAMP}"
OUTPUT_DIR="${WORK_OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}"

# Data files
INPUT_BAM="/iblm/netapp/data1/aho/atac/Buenrostro2013/merged/GM12878_ATACseq_50k/GM12878_ATACseq_50k_merged.sorted.bam"
VCF="/iblm/netapp/data1/aho/variants/NA12878.vcf.gz"
REF_GENOME="/iblm/netapp/data1/aho/ref_genomes/index/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/genome.fa"

# Sample
SAMPLE="NA12878"
THREADS=8

# Conda
source /iblm/netapp/home/jjaureguy/mambaforge/etc/profile.d/conda.sh
conda activate WASP2_dev2

# Tools
BWA=$(which bwa)
SAMTOOLS=$(which samtools)
PYTHON=$(which python)
TIME="/usr/bin/time"

# Set PYTHONPATH
export PYTHONPATH="${WASP2_DIR}/src:${PYTHONPATH}"

cd "${OUTPUT_DIR}"

# Log file
PROFILE_LOG="${OUTPUT_DIR}/profile.log"
echo "WASP2-Rust FIXED Pipeline Benchmark - GM12878 ATAC-seq" > ${PROFILE_LOG}
echo "========================================================" >> ${PROFILE_LOG}
echo "Timestamp: ${TIMESTAMP}" >> ${PROFILE_LOG}
echo "Date: $(date)" >> ${PROFILE_LOG}
echo "Sample: ${SAMPLE}" >> ${PROFILE_LOG}
echo "Threads: ${THREADS}" >> ${PROFILE_LOG}
echo "" >> ${PROFILE_LOG}
echo "FIXED PIPELINE: filter → unified → remap → filter → MERGE" >> ${PROFILE_LOG}
echo "" >> ${PROFILE_LOG}

echo "========================================"
echo "WASP2-Rust FIXED Benchmark on GM12878 ATAC-seq"
echo "Timestamp: ${TIMESTAMP}"
echo "========================================"
echo "CORRECT pipeline: filter → unified → remap → filter → MERGE"
echo "Start time: $(date)"
echo "Input BAM: ${INPUT_BAM}"
echo "VCF: ${VCF}"
echo "Threads: ${THREADS}"
echo ""

# Verify Rust module
echo "Verifying Rust module..."
${PYTHON} -c "import wasp2_rust; print('wasp2_rust module: AVAILABLE')" || {
    echo "ERROR: wasp2_rust module not available"
    exit 1
}
${PYTHON} -c "from wasp2_rust import filter_bam_by_variants_py, unified_make_reads_parallel_py; print('All functions: AVAILABLE')" || {
    echo "ERROR: Required functions not available"
    exit 1
}
echo ""

# -----------------------------------------------------------------------------
# PRE-BENCHMARK: Copy input BAM to working directory
# -----------------------------------------------------------------------------
echo "Copying input BAM to output directory..."
cp "${INPUT_BAM}" "${OUTPUT_DIR}/original.bam"
${SAMTOOLS} index "${OUTPUT_DIR}/original.bam"
ORIGINAL_BAM="${OUTPUT_DIR}/original.bam"

# Count original reads
ORIGINAL_READS=$(${SAMTOOLS} view -c ${ORIGINAL_BAM})
echo "Original reads: ${ORIGINAL_READS}"

# Create variant BED file using Python vcf_to_bed (produces correct 6-column format)
echo "Creating variant BED file using vcf_to_bed..."
${PYTHON} -c "
from mapping.intersect_variant_data import vcf_to_bed

vcf_to_bed(
    vcf_file='${VCF}',
    out_bed='${OUTPUT_DIR}/variants.bed',
    samples=['${SAMPLE}'],
    include_indels=False  # SNPs only
)
print('BED file created successfully')
"
VARIANT_COUNT=$(wc -l < ${OUTPUT_DIR}/variants.bed)
echo "SNP variants in BED: ${VARIANT_COUNT}"

# -----------------------------------------------------------------------------
# BENCHMARK STARTS HERE
# -----------------------------------------------------------------------------
echo "========================================"
echo "BENCHMARK TIMING STARTS NOW"
echo "========================================"

TOTAL_START=$(date +%s.%N)

# -----------------------------------------------------------------------------
# STEP 1: Filter BAM by variants → get keep.bam + to_remap.bam
# THIS IS THE MISSING STEP IN THE OLD BENCHMARK!
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 1: Filter BAM by variants (split into keep + to_remap)"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 1: filter_bam_by_variants (split BAM)" >> ${PROFILE_LOG}

STEP1_START=$(date +%s.%N)

${TIME} -v ${PYTHON} -c "
import wasp2_rust

remap_count, keep_count, unique_names = wasp2_rust.filter_bam_by_variants_py(
    '${ORIGINAL_BAM}',
    '${OUTPUT_DIR}/variants.bed',
    '${OUTPUT_DIR}/to_remap.bam',
    '${OUTPUT_DIR}/keep.bam',
    is_paired=True,
    threads=${THREADS}
)

print(f'Reads to remap: {remap_count:,}')
print(f'Reads to keep (no variants): {keep_count:,}')
print(f'Unique read names: {unique_names:,}')
" 2>> ${PROFILE_LOG}

${SAMTOOLS} index ${OUTPUT_DIR}/to_remap.bam
${SAMTOOLS} index ${OUTPUT_DIR}/keep.bam

STEP1_END=$(date +%s.%N)
STEP1_TIME=$(echo "${STEP1_END} - ${STEP1_START}" | bc)
echo "STEP 1 completed in ${STEP1_TIME} seconds"

# Count keep reads
KEEP_READS=$(${SAMTOOLS} view -c ${OUTPUT_DIR}/keep.bam)
echo "Keep reads (no variants): ${KEEP_READS}"

# -----------------------------------------------------------------------------
# STEP 2: WASP2-Rust unified make-reads (on to_remap.bam ONLY!)
# NOW WITH keep_no_flip_names_path to capture reads that overlap variants but don't need flipping
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 2: WASP2-Rust unified make-reads (on to_remap.bam only)"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 2: unified make-reads (from to_remap.bam)" >> ${PROFILE_LOG}

STEP2_START=$(date +%s.%N)

${TIME} -v ${PYTHON} -c "
from wasp2_rust import unified_make_reads_parallel_py
import json

stats = unified_make_reads_parallel_py(
    '${OUTPUT_DIR}/to_remap.bam',
    '${OUTPUT_DIR}/variants.bed',
    '${OUTPUT_DIR}/remap_r1.fq',
    '${OUTPUT_DIR}/remap_r2.fq',
    max_seqs=64,
    threads=${THREADS},
    keep_no_flip_names_path='${OUTPUT_DIR}/keep_no_flip_names.txt',
    remap_names_path='${OUTPUT_DIR}/remap_names.txt'
)

with open('${OUTPUT_DIR}/unified_stats.json', 'w') as f:
    json.dump(stats, f, indent=2)

print(f'Pairs processed: {stats[\"pairs_processed\"]:,}')
print(f'Pairs keep-no-flip: {stats[\"pairs_keep_no_flip\"]:,}')
print(f'Haplotypes written: {stats[\"haplotypes_written\"]:,}')
" 2>> ${PROFILE_LOG}

STEP2_END=$(date +%s.%N)
STEP2_TIME=$(echo "${STEP2_END} - ${STEP2_START}" | bc)
echo "STEP 2 completed in ${STEP2_TIME} seconds"

# Count reads to remap
if [ -f "${OUTPUT_DIR}/remap_r1.fq" ]; then
    REMAP_COUNT=$(wc -l < ${OUTPUT_DIR}/remap_r1.fq)
    REMAP_PAIRS=$((REMAP_COUNT / 4))
    echo "Reads to remap: ${REMAP_PAIRS} pairs"
fi

# Count keep-no-flip reads
if [ -f "${OUTPUT_DIR}/keep_no_flip_names.txt" ]; then
    KEEP_NO_FLIP_NAMES=$(wc -l < ${OUTPUT_DIR}/keep_no_flip_names.txt)
    echo "Keep-no-flip read names: ${KEEP_NO_FLIP_NAMES}"
fi

# Count remap reads
if [ -f "${OUTPUT_DIR}/remap_names.txt" ]; then
    REMAP_NAMES=$(wc -l < ${OUTPUT_DIR}/remap_names.txt)
    echo "Remap read names: ${REMAP_NAMES}"
fi

# -----------------------------------------------------------------------------
# STEP 2b: Create to_remap_actual.bam with only reads that were sent for remapping
# This is CRITICAL for correct filtering - the filter compares remapped BAM to
# to_remap BAM, so they must have the same read names!
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 2b: Extract to_remap_actual.bam (only reads sent for remapping)"
echo "========================================"

STEP2B_START=$(date +%s.%N)

if [ -s "${OUTPUT_DIR}/remap_names.txt" ]; then
    ${SAMTOOLS} view -@ ${THREADS} -N ${OUTPUT_DIR}/remap_names.txt \
        -o ${OUTPUT_DIR}/to_remap_actual.bam ${OUTPUT_DIR}/to_remap.bam
    ${SAMTOOLS} index ${OUTPUT_DIR}/to_remap_actual.bam
    TO_REMAP_ACTUAL_READS=$(${SAMTOOLS} view -c ${OUTPUT_DIR}/to_remap_actual.bam)
    echo "to_remap_actual.bam reads: ${TO_REMAP_ACTUAL_READS}"
else
    echo "No remap names found - using full to_remap.bam"
    cp ${OUTPUT_DIR}/to_remap.bam ${OUTPUT_DIR}/to_remap_actual.bam
    ${SAMTOOLS} index ${OUTPUT_DIR}/to_remap_actual.bam
fi

STEP2B_END=$(date +%s.%N)
STEP2B_TIME=$(echo "${STEP2B_END} - ${STEP2B_START}" | bc)
echo "STEP 2b completed in ${STEP2B_TIME} seconds"

# -----------------------------------------------------------------------------
# STEP 3: BWA-MEM remap flipped reads
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 3: BWA-MEM remap flipped reads"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 3: BWA-MEM Remap" >> ${PROFILE_LOG}

STEP3_START=$(date +%s.%N)

${TIME} -v ${BWA} mem -t ${THREADS} -M ${REF_GENOME} \
    ${OUTPUT_DIR}/remap_r1.fq ${OUTPUT_DIR}/remap_r2.fq 2>> ${PROFILE_LOG} | \
    ${SAMTOOLS} view -S -b -h -F 4 - > ${OUTPUT_DIR}/remapped_unsorted.bam

${SAMTOOLS} sort -@ ${THREADS} -o ${OUTPUT_DIR}/remapped.bam ${OUTPUT_DIR}/remapped_unsorted.bam
${SAMTOOLS} index ${OUTPUT_DIR}/remapped.bam
rm -f ${OUTPUT_DIR}/remapped_unsorted.bam

STEP3_END=$(date +%s.%N)
STEP3_TIME=$(echo "${STEP3_END} - ${STEP3_START}" | bc)
echo "STEP 3 completed in ${STEP3_TIME} seconds"

# -----------------------------------------------------------------------------
# STEP 4: WASP2-Rust filter-remapped
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 4: WASP2-Rust filter-remapped"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 4: filter_bam_wasp" >> ${PROFILE_LOG}

STEP4_START=$(date +%s.%N)

${TIME} -v ${PYTHON} -c "
from wasp2_rust import filter_bam_wasp

# Use to_remap_actual.bam which contains only reads that were sent for remapping
# This matches Python pipeline behavior where to_remap.bam only has reads needing flipping
kept, removed_moved, removed_missing = filter_bam_wasp(
    '${OUTPUT_DIR}/to_remap_actual.bam',
    '${OUTPUT_DIR}/remapped.bam',
    '${OUTPUT_DIR}/remap_keep.bam',
    threads=${THREADS}
)

print(f'Kept: {kept}, Removed (moved): {removed_moved}, Removed (missing): {removed_missing}')
" 2>> ${PROFILE_LOG}

${SAMTOOLS} index ${OUTPUT_DIR}/remap_keep.bam

STEP4_END=$(date +%s.%N)
STEP4_TIME=$(echo "${STEP4_END} - ${STEP4_START}" | bc)
echo "STEP 4 completed in ${STEP4_TIME} seconds"

# -----------------------------------------------------------------------------
# STEP 5: Extract keep-no-flip reads from to_remap.bam
# These are reads that overlap variants but don't need allele flipping
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 5: Extract keep-no-flip reads from to_remap.bam"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 5: Extract keep-no-flip reads" >> ${PROFILE_LOG}

STEP5_START=$(date +%s.%N)

if [ -s "${OUTPUT_DIR}/keep_no_flip_names.txt" ]; then
    # Extract reads by name from to_remap.bam
    ${TIME} -v ${SAMTOOLS} view -@ ${THREADS} -N ${OUTPUT_DIR}/keep_no_flip_names.txt \
        -o ${OUTPUT_DIR}/keep_no_flip.bam ${OUTPUT_DIR}/to_remap.bam 2>> ${PROFILE_LOG}
    ${SAMTOOLS} index ${OUTPUT_DIR}/keep_no_flip.bam
    KEEP_NO_FLIP_READS=$(${SAMTOOLS} view -c ${OUTPUT_DIR}/keep_no_flip.bam)
    echo "Keep-no-flip reads extracted: ${KEEP_NO_FLIP_READS}"
else
    echo "No keep-no-flip reads to extract"
    KEEP_NO_FLIP_READS=0
    touch ${OUTPUT_DIR}/keep_no_flip.bam  # Empty placeholder
fi

STEP5_END=$(date +%s.%N)
STEP5_TIME=$(echo "${STEP5_END} - ${STEP5_START}" | bc)
echo "STEP 5 completed in ${STEP5_TIME} seconds"

# -----------------------------------------------------------------------------
# STEP 6: MERGE keep.bam + keep_no_flip.bam + remap_keep.bam → wasp_filtered.bam
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 6: Merge keep + keep_no_flip + remap_keep → wasp_filtered.bam"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 6: Merge BAMs" >> ${PROFILE_LOG}

STEP6_START=$(date +%s.%N)

# Build merge command - only include non-empty BAMs
MERGE_INPUTS="${OUTPUT_DIR}/keep.bam"
if [ -s "${OUTPUT_DIR}/keep_no_flip.bam" ]; then
    MERGE_INPUTS="${MERGE_INPUTS} ${OUTPUT_DIR}/keep_no_flip.bam"
fi
MERGE_INPUTS="${MERGE_INPUTS} ${OUTPUT_DIR}/remap_keep.bam"

${TIME} -v ${SAMTOOLS} merge -@ ${THREADS} -o ${OUTPUT_DIR}/wasp_filtered_unsorted.bam \
    ${MERGE_INPUTS} 2>> ${PROFILE_LOG}

${SAMTOOLS} sort -@ ${THREADS} -o ${OUTPUT_DIR}/wasp_filtered.bam ${OUTPUT_DIR}/wasp_filtered_unsorted.bam
${SAMTOOLS} index ${OUTPUT_DIR}/wasp_filtered.bam
rm -f ${OUTPUT_DIR}/wasp_filtered_unsorted.bam

STEP6_END=$(date +%s.%N)
STEP6_TIME=$(echo "${STEP6_END} - ${STEP6_START}" | bc)
echo "STEP 6 completed in ${STEP6_TIME} seconds"

# -----------------------------------------------------------------------------
# SUMMARY
# -----------------------------------------------------------------------------
TOTAL_END=$(date +%s.%N)
TOTAL_TIME=$(echo "${TOTAL_END} - ${TOTAL_START}" | bc)

# Count final reads
FINAL_READS=$(${SAMTOOLS} view -c ${OUTPUT_DIR}/wasp_filtered.bam)
REMAP_KEEP_READS=$(${SAMTOOLS} view -c ${OUTPUT_DIR}/remap_keep.bam)
PASS_RATE=$(echo "scale=2; ${FINAL_READS} * 100 / ${ORIGINAL_READS}" | bc)

# WASP-only time (steps 1+2+2b+4+5+6, excluding BWA alignment)
WASP_ONLY=$(echo "${STEP1_TIME} + ${STEP2_TIME} + ${STEP2B_TIME} + ${STEP4_TIME} + ${STEP5_TIME} + ${STEP6_TIME}" | bc)

echo "========================================"
echo "BENCHMARK RESULTS (FIXED PIPELINE v3 - with correct to_remap reference)"
echo "========================================"
echo ""
echo "Individual step times:"
echo "  STEP 1 (filter BAM by variants):   ${STEP1_TIME} s"
echo "  STEP 2 (unified make-reads):       ${STEP2_TIME} s"
echo "  STEP 2b (extract to_remap_actual): ${STEP2B_TIME} s"
echo "  STEP 3 (BWA-MEM remap):            ${STEP3_TIME} s"
echo "  STEP 4 (filter remapped):          ${STEP4_TIME} s"
echo "  STEP 5 (extract keep-no-flip):     ${STEP5_TIME} s"
echo "  STEP 6 (merge BAMs):               ${STEP6_TIME} s"
echo ""
echo "WASP-only time (steps 1+2+2b+4+5+6): ${WASP_ONLY} s"
echo "TOTAL time:                          ${TOTAL_TIME} s"
echo ""
echo "Read counts:"
echo "  Original reads:       ${ORIGINAL_READS}"
echo "  Keep reads (no var):  ${KEEP_READS}"
echo "  Keep-no-flip reads:   ${KEEP_NO_FLIP_READS}"
echo "  Remap keep reads:     ${REMAP_KEEP_READS}"
echo "  FINAL (merged):       ${FINAL_READS}"
echo "  Pass rate:            ${PASS_RATE}%"
echo ""
echo "Work directory: ${OUTPUT_DIR}"
echo "Final output directory: ${FINAL_OUTPUT_DIR}"
echo "End time: $(date)"

# Save results to JSON
cat > ${OUTPUT_DIR}/benchmark_results.json << EOF
{
    "timestamp": "${TIMESTAMP}",
    "pipeline": "wasp2rust_snp_FIXED_v3",
    "sample": "${SAMPLE}",
    "data_type": "ATAC-seq",
    "threads": ${THREADS},
    "step1_filter_bam_s": ${STEP1_TIME},
    "step2_unified_make_reads_s": ${STEP2_TIME},
    "step2b_extract_to_remap_actual_s": ${STEP2B_TIME},
    "step3_bwa_remap_s": ${STEP3_TIME},
    "step4_filter_remapped_s": ${STEP4_TIME},
    "step5_extract_keep_no_flip_s": ${STEP5_TIME},
    "step6_merge_bams_s": ${STEP6_TIME},
    "wasp_only_s": ${WASP_ONLY},
    "total_s": ${TOTAL_TIME},
    "original_reads": ${ORIGINAL_READS},
    "keep_reads": ${KEEP_READS},
    "keep_no_flip_reads": ${KEEP_NO_FLIP_READS},
    "remap_keep_reads": ${REMAP_KEEP_READS},
    "final_reads": ${FINAL_READS},
    "pass_rate_percent": ${PASS_RATE},
    "snv_count": ${VARIANT_COUNT}
}
EOF

echo ""
echo "Results saved to: ${OUTPUT_DIR}/benchmark_results.json"
echo ""
echo "COMPARISON TO WASP1:"
echo "  WASP1 total:    5523.97 s (147M reads output)"
echo "  WASP2 total:    ${TOTAL_TIME} s (${FINAL_READS} reads output)"
echo "  WASP2 speedup:  $(echo "scale=2; 5523.97 / ${TOTAL_TIME}" | bc)x"

# -----------------------------------------------------------------------------
# Copy artifacts back to final output directory (NFS)
# -----------------------------------------------------------------------------
echo ""
echo "Copying benchmark artifacts to final output directory..."
mkdir -p "${FINAL_OUTPUT_DIR}"
for f in \
    "${OUTPUT_DIR}/benchmark_results.json" \
    "${PROFILE_LOG}" \
    "${OUTPUT_DIR}/unified_stats.json" \
    "${OUTPUT_DIR}/remap_keep.bam" \
    "${OUTPUT_DIR}/remap_keep.bam.bai" \
    "${OUTPUT_DIR}/wasp_filtered.bam" \
    "${OUTPUT_DIR}/wasp_filtered.bam.bai"
do
    if [ -f "${f}" ]; then
        rsync -a "${f}" "${FINAL_OUTPUT_DIR}/"
    fi
done
echo "Final output directory: ${FINAL_OUTPUT_DIR}"
cd /
rm -rf "${WORK_OUTPUT_DIR}"
