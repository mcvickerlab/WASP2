#!/bin/bash
# WASP2-Rust FAIR benchmark on HG00731 RNA-seq data
#
# FAIR COMPARISON vs STAR+WASP:
# 1. Uses --genomeLoad LoadAndKeep (load genome once, reuse for remap)
# 2. Excludes samtools index from STEP 1 timing (STAR+WASP doesn't do this)
# 3. Same STAR parameters as STAR+WASP benchmark
# 4. Correct reverse complement handling (FIXED code)
#
#$ -N wasp2rust_rnaseq_fair
#$ -V
#$ -pe iblm 8
#$ -l h_vmem=32G
#$ -j y
#$ -o /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/logs/
#$ -cwd

set -e

# Increase file descriptor limit for STAR BAM sorting
ulimit -n 10000

# Timestamp for this run
TIMESTAMP=$(date +%Y-%m-%d_%H-%M-%S)
echo "Benchmark timestamp: ${TIMESTAMP}"

# Paths
WASP2_DIR="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp"
BENCHMARK_DIR="${WASP2_DIR}/benchmarking/star_wasp_comparison"
DATA_DIR="${BENCHMARK_DIR}/data"
FINAL_OUTPUT_DIR="${BENCHMARK_DIR}/results/wasp2rust_fair_${TIMESTAMP}"

# Data files
FASTQ_R1="${DATA_DIR}/ERR1050079_1.fastq.gz"
FASTQ_R2="${DATA_DIR}/ERR1050079_2.fastq.gz"
VCF="${DATA_DIR}/HG00731_het_only_chr.vcf.gz"

# Reference - STAR index
STAR_INDEX_REMOTE="/iblm/netapp/data1/external/GRC38/combined/google_cloud/star_index"
LOCAL_SCRATCH="${TMPDIR:-/tmp}"
WORK_OUTPUT_DIR="${LOCAL_SCRATCH}/wasp2rust_fair_${JOB_ID:-local}_${TIMESTAMP}"
OUTPUT_DIR="${WORK_OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}"

# Sample
SAMPLE="HG00731"
THREADS=8

# Conda
source /iblm/netapp/home/jjaureguy/mambaforge/etc/profile.d/conda.sh
conda activate WASP2_dev2

# Tools
STAR=$(which STAR)
SAMTOOLS=$(which samtools)
PYTHON=$(which python)
TIME="/usr/bin/time"

# Set PYTHONPATH
export PYTHONPATH="${WASP2_DIR}/src:${PYTHONPATH}"

cd "${OUTPUT_DIR}"

# Log file
PROFILE_LOG="${OUTPUT_DIR}/profile.log"
echo "WASP2-Rust FAIR Pipeline Benchmark - HG00731 RNA-seq" > ${PROFILE_LOG}
echo "======================================================" >> ${PROFILE_LOG}
echo "Timestamp: ${TIMESTAMP}" >> ${PROFILE_LOG}
echo "Date: $(date)" >> ${PROFILE_LOG}
echo "Sample: ${SAMPLE}" >> ${PROFILE_LOG}
echo "Threads: ${THREADS}" >> ${PROFILE_LOG}
echo "" >> ${PROFILE_LOG}
echo "FAIR COMPARISON METHODOLOGY:" >> ${PROFILE_LOG}
echo "  - genomeLoad: LoadAndKeep (load once, reuse for remap)" >> ${PROFILE_LOG}
echo "  - samtools index NOT included in STAR timing" >> ${PROFILE_LOG}
echo "  - Same STAR params as STAR+WASP benchmark" >> ${PROFILE_LOG}
echo "" >> ${PROFILE_LOG}

echo "========================================"
echo "WASP2-Rust FAIR Benchmark on HG00731 RNA-seq"
echo "Timestamp: ${TIMESTAMP}"
echo "========================================"
echo "Start time: $(date)"
echo "FASTQs: ${FASTQ_R1}"
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
# PRE-BENCHMARK: Stage files to local scratch (NOT timed)
# -----------------------------------------------------------------------------
echo "========================================"
echo "PRE-BENCHMARK: Staging files (NOT timed)"
echo "========================================"
STAGE_START=$(date +%s.%N)

LOCAL_STAR_INDEX="${LOCAL_SCRATCH}/wasp2_fair_star_index_${JOB_ID:-local}"
if [ -d "${LOCAL_STAR_INDEX}" ]; then
    echo "Local STAR index already exists, reusing..."
else
    mkdir -p "${LOCAL_STAR_INDEX}"
    rsync -a "${STAR_INDEX_REMOTE}/" "${LOCAL_STAR_INDEX}/"
fi
STAR_INDEX="${LOCAL_STAR_INDEX}"

# Copy FASTQs to local scratch
LOCAL_FASTQ_DIR="${LOCAL_SCRATCH}/wasp2_fair_fastqs_${TIMESTAMP}"
mkdir -p "${LOCAL_FASTQ_DIR}"
LOCAL_FASTQ_R1="${LOCAL_FASTQ_DIR}/R1.fastq.gz"
LOCAL_FASTQ_R2="${LOCAL_FASTQ_DIR}/R2.fastq.gz"
cp "${FASTQ_R1}" "${LOCAL_FASTQ_R1}"
cp "${FASTQ_R2}" "${LOCAL_FASTQ_R2}"
echo "FASTQs staged to local scratch"

# Create variant BED file using Python vcf_to_bed
# Include indels to match STAR+WASP (which uses full VCF)
echo "Creating variant BED file using vcf_to_bed..."
${PYTHON} -c "
from mapping.intersect_variant_data import vcf_to_bed

vcf_to_bed(
    vcf_file='${VCF}',
    out_bed='${OUTPUT_DIR}/variants.bed',
    samples=['${SAMPLE}'],
    include_indels=True  # Include ALL variants for fair comparison
)
print('BED file created successfully')
"
VARIANT_COUNT=$(wc -l < ${OUTPUT_DIR}/variants.bed)
echo "Variants in BED (SNPs + indels): ${VARIANT_COUNT}"

STAGE_END=$(date +%s.%N)
STAGE_TIME=$(echo "${STAGE_END} - ${STAGE_START}" | bc)
echo "Staging time: ${STAGE_TIME}s (NOT included in benchmark)"
echo "" >> ${PROFILE_LOG}
echo "Staging time: ${STAGE_TIME}s (NOT included in benchmark)" >> ${PROFILE_LOG}

# STAR parameters - SAME as STAR+WASP benchmark (minus WASP-specific options)
# Key: --genomeLoad LoadAndKeep to cache genome for remap pass
# Note: limitBAMsortRAM required when using LoadAndKeep with SortedByCoordinate output
STARpar_initial="--runThreadN ${THREADS} \
         --genomeDir ${STAR_INDEX} \
         --genomeLoad LoadAndKeep \
         --limitBAMsortRAM 15000000000 \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMattributes NH HI AS nM NM MD jM jI \
         --alignEndsType EndToEnd \
         --outSAMunmapped Within \
         --outFilterMultimapNmax 1"

STARpar_remap="--runThreadN ${THREADS} \
         --genomeDir ${STAR_INDEX} \
         --genomeLoad LoadAndKeep \
         --outSAMtype BAM Unsorted \
         --outSAMattributes NH HI AS nM NM MD jM jI \
         --alignEndsType EndToEnd \
         --outSAMunmapped Within \
         --outFilterMultimapNmax 1"

# -----------------------------------------------------------------------------
# BENCHMARK STARTS HERE
# -----------------------------------------------------------------------------
echo "========================================"
echo "BENCHMARK TIMING STARTS NOW"
echo "========================================"

TOTAL_START=$(date +%s.%N)

# -----------------------------------------------------------------------------
# STEP 1: Initial STAR alignment (timing excludes samtools index)
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 1: Initial STAR alignment"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 1: STAR initial alignment" >> ${PROFILE_LOG}

STEP1_START=$(date +%s.%N)

${TIME} -v ${STAR} ${STARpar_initial} \
    --readFilesCommand gunzip -c \
    --readFilesIn ${LOCAL_FASTQ_R1} ${LOCAL_FASTQ_R2} 2>> ${PROFILE_LOG}

mv Aligned.sortedByCoord.out.bam original.bam

STEP1_END=$(date +%s.%N)
STEP1_TIME=$(echo "${STEP1_END} - ${STEP1_START}" | bc)
echo "STEP 1 (STAR only) completed in ${STEP1_TIME} seconds"

# Index BAM (timed separately, not part of STAR comparison)
INDEX_START=$(date +%s.%N)
${SAMTOOLS} index original.bam
INDEX_END=$(date +%s.%N)
INDEX_TIME=$(echo "${INDEX_END} - ${INDEX_START}" | bc)
echo "BAM indexing took ${INDEX_TIME} seconds (separate from STAR timing)"

ORIGINAL_READS=$(${SAMTOOLS} view -c original.bam)
echo "Original aligned reads: ${ORIGINAL_READS}"

# -----------------------------------------------------------------------------
# STEP 2: Filter BAM by variants → get keep.bam + to_remap.bam
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 2: Filter BAM by variants (split into keep + to_remap)"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 2: filter_bam_by_variants (split BAM)" >> ${PROFILE_LOG}

STEP2_START=$(date +%s.%N)

${TIME} -v ${PYTHON} -c "
import wasp2_rust

remap_count, keep_count, unique_names = wasp2_rust.filter_bam_by_variants_py(
    '${OUTPUT_DIR}/original.bam',
    '${OUTPUT_DIR}/variants.bed',
    '${OUTPUT_DIR}/to_remap.bam',
    '${OUTPUT_DIR}/keep.bam',
    is_paired=True,
    threads=${THREADS}
)

print(f'Reads to remap: {remap_count:,}')
print(f'Reads to keep (no variants): {keep_count:,}')
print(f'Unique read names: {unique_names}')
" 2>> ${PROFILE_LOG}

${SAMTOOLS} index ${OUTPUT_DIR}/to_remap.bam
${SAMTOOLS} index ${OUTPUT_DIR}/keep.bam

STEP2_END=$(date +%s.%N)
STEP2_TIME=$(echo "${STEP2_END} - ${STEP2_START}" | bc)
echo "STEP 2 completed in ${STEP2_TIME} seconds"

KEEP_READS=$(${SAMTOOLS} view -c ${OUTPUT_DIR}/keep.bam)
echo "Keep reads (no variants): ${KEEP_READS}"

# -----------------------------------------------------------------------------
# STEP 3: WASP2-Rust unified make-reads (on to_remap.bam ONLY!)
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 3: WASP2-Rust unified make-reads (on to_remap.bam only)"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 3: unified make-reads (from to_remap.bam)" >> ${PROFILE_LOG}

STEP3_START=$(date +%s.%N)

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
    compress_output=False
)

with open('${OUTPUT_DIR}/unified_stats.json', 'w') as f:
    json.dump(dict(stats), f, indent=2)

print(f'Pairs processed: {stats[\"pairs_processed\"]:,}')
print(f'Haplotypes written: {stats[\"haplotypes_written\"]:,}')
" 2>> ${PROFILE_LOG}

STEP3_END=$(date +%s.%N)
STEP3_TIME=$(echo "${STEP3_END} - ${STEP3_START}" | bc)
echo "STEP 3 completed in ${STEP3_TIME} seconds"

# Count reads to remap
if [ -f "${OUTPUT_DIR}/remap_r1.fq" ]; then
    REMAP_COUNT=$(wc -l < ${OUTPUT_DIR}/remap_r1.fq)
    REMAP_PAIRS=$((REMAP_COUNT / 4))
    echo "Reads to remap: ${REMAP_PAIRS} pairs"
fi

# -----------------------------------------------------------------------------
# STEP 4: STAR remap flipped reads (genome already loaded via LoadAndKeep)
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 4: STAR remap flipped reads (genome cached)"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 4: STAR Remap (genome cached via LoadAndKeep)" >> ${PROFILE_LOG}

STEP4_START=$(date +%s.%N)

${TIME} -v ${STAR} ${STARpar_remap} \
    --readFilesIn ${OUTPUT_DIR}/remap_r1.fq ${OUTPUT_DIR}/remap_r2.fq 2>> ${PROFILE_LOG}

mv Aligned.out.bam ${OUTPUT_DIR}/remapped.bam

# Clean up uncompressed FASTQs
rm -f ${OUTPUT_DIR}/remap_r1.fq ${OUTPUT_DIR}/remap_r2.fq

STEP4_END=$(date +%s.%N)
STEP4_TIME=$(echo "${STEP4_END} - ${STEP4_START}" | bc)
echo "STEP 4 completed in ${STEP4_TIME} seconds"

# -----------------------------------------------------------------------------
# STEP 5: WASP2-Rust filter-remapped
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 5: WASP2-Rust filter-remapped"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 5: filter_bam_wasp" >> ${PROFILE_LOG}

STEP5_START=$(date +%s.%N)

${TIME} -v ${PYTHON} -c "
from wasp2_rust import filter_bam_wasp

kept, removed_moved, removed_missing = filter_bam_wasp(
    '${OUTPUT_DIR}/to_remap.bam',
    '${OUTPUT_DIR}/remapped.bam',
    '${OUTPUT_DIR}/remap_keep.bam',
    threads=${THREADS}
)

print(f'Kept: {kept}, Removed (moved): {removed_moved}, Removed (missing): {removed_missing}')
" 2>> ${PROFILE_LOG}

${SAMTOOLS} index ${OUTPUT_DIR}/remap_keep.bam

STEP5_END=$(date +%s.%N)
STEP5_TIME=$(echo "${STEP5_END} - ${STEP5_START}" | bc)
echo "STEP 5 completed in ${STEP5_TIME} seconds"

# -----------------------------------------------------------------------------
# STEP 6: MERGE keep.bam + remap_keep.bam → wasp_filtered.bam
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 6: Merge keep + remap_keep → wasp_filtered.bam"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 6: Merge BAMs" >> ${PROFILE_LOG}

STEP6_START=$(date +%s.%N)

${TIME} -v ${SAMTOOLS} merge -@ ${THREADS} -o ${OUTPUT_DIR}/wasp_filtered_unsorted.bam \
    ${OUTPUT_DIR}/keep.bam ${OUTPUT_DIR}/remap_keep.bam 2>> ${PROFILE_LOG}

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

# WASP-only time (steps 2+3+5+6, excluding STAR alignments)
WASP_ONLY=$(echo "${STEP2_TIME} + ${STEP3_TIME} + ${STEP5_TIME} + ${STEP6_TIME}" | bc)

# STAR-only time (steps 1+4)
STAR_ONLY=$(echo "${STEP1_TIME} + ${STEP4_TIME}" | bc)

# Remove genome from shared memory
echo "Removing genome from shared memory..."
${STAR} --genomeDir ${STAR_INDEX} --genomeLoad Remove 2>/dev/null || true

# Cleanup local scratch
rm -rf "${LOCAL_FASTQ_DIR}"

echo "========================================"
echo "BENCHMARK RESULTS (FAIR COMPARISON)"
echo "========================================"
echo ""
echo "Individual step times:"
echo "  STEP 1 (STAR initial):            ${STEP1_TIME} s"
echo "  STEP 2 (filter BAM by variants):  ${STEP2_TIME} s"
echo "  STEP 3 (unified make-reads):      ${STEP3_TIME} s"
echo "  STEP 4 (STAR remap, cached):      ${STEP4_TIME} s"
echo "  STEP 5 (filter remapped):         ${STEP5_TIME} s"
echo "  STEP 6 (merge BAMs):              ${STEP6_TIME} s"
echo ""
echo "Aggregate times:"
echo "  STAR-only (steps 1+4):            ${STAR_ONLY} s"
echo "  WASP-only (steps 2+3+5+6):        ${WASP_ONLY} s"
echo "  TOTAL time:                       ${TOTAL_TIME} s"
echo ""
echo "Read counts:"
echo "  Original reads:       ${ORIGINAL_READS}"
echo "  Keep reads:           ${KEEP_READS}"
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
    "pipeline": "wasp2rust_FAIR",
    "methodology": "LoadAndKeep genome caching, samtools index excluded from STAR timing",
    "sample": "${SAMPLE}",
    "data_type": "RNA-seq",
    "threads": ${THREADS},
    "step1_star_initial_s": ${STEP1_TIME},
    "step1_bam_index_s": ${INDEX_TIME},
    "step2_filter_bam_s": ${STEP2_TIME},
    "step3_unified_make_reads_s": ${STEP3_TIME},
    "step4_star_remap_s": ${STEP4_TIME},
    "step5_filter_remapped_s": ${STEP5_TIME},
    "step6_merge_bams_s": ${STEP6_TIME},
    "star_only_s": ${STAR_ONLY},
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
echo ""
echo "FAIR COMPARISON NOTES:"
echo "  - STAR+WASP baseline: 532.6s (8.9m)"
echo "  - This run STAR-only: ${STAR_ONLY}s ($(echo "scale=1; ${STAR_ONLY}/60" | bc)m)"
echo "  - This run WASP-only: ${WASP_ONLY}s ($(echo "scale=1; ${WASP_ONLY}/60" | bc)m)"
echo "  - This run TOTAL:     ${TOTAL_TIME}s ($(echo "scale=1; ${TOTAL_TIME}/60" | bc)m)"

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
    "${OUTPUT_DIR}/Log.final.out"
do
    if [ -f "${f}" ]; then
        rsync -a "${f}" "${FINAL_OUTPUT_DIR}/"
    fi
done

# Optional: keep large BAM outputs (merged/kept reads)
if [ "${KEEP_BAMS:-0}" = "1" ]; then
    for f in "${OUTPUT_DIR}/wasp_filtered.bam" "${OUTPUT_DIR}/wasp_filtered.bam.bai"; do
        if [ -f "${f}" ]; then
            rsync -a "${f}" "${FINAL_OUTPUT_DIR}/"
        fi
    done
fi

echo "Final output directory: ${FINAL_OUTPUT_DIR}"
cd /
rm -rf "${WORK_OUTPUT_DIR}"
