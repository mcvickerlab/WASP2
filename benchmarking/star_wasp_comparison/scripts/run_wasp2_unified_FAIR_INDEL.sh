#!/bin/bash
# WASP2-Rust FAIR benchmark (SNV+INDEL) on HG00731 RNA-seq data
#
# Goal: apples-to-apples with `run_wasp2_unified_FAIR.sh` except INDEL mode is enabled.
# This avoids comparing the older "FAIR split+merge" pipeline to the newer v2 indel script
# (which has different steps and therefore different total_s).
#
# FAIR COMPARISON vs STAR+WASP:
# 1. Uses --genomeLoad LoadAndKeep (load genome once, reuse for remap)
# 2. Excludes samtools index from STEP 1 timing (STAR+WASP doesn't do this)
# 3. Same STAR parameters as STAR+WASP benchmark
#
#$ -N wasp2rust_rnaseq_fair_indel
#$ -V
#$ -pe iblm 8
#$ -l h_vmem=32G
#$ -j y
#$ -o /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/logs/
#$ -cwd

set -e

# Increase file descriptor limit for STAR BAM sorting
ulimit -n 10000

TIMESTAMP=$(date +%Y-%m-%d_%H-%M-%S)
echo "Benchmark timestamp: ${TIMESTAMP}"

# Paths
WASP2_DIR="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp"
BENCHMARK_DIR="${WASP2_DIR}/benchmarking/star_wasp_comparison"
DATA_DIR="${BENCHMARK_DIR}/data"
FINAL_OUTPUT_DIR="${BENCHMARK_DIR}/results/wasp2rust_fair_indel_${TIMESTAMP}"

# Provenance
GIT_SHA=$(git -C "${WASP2_DIR}" rev-parse HEAD 2>/dev/null || echo "unknown")

# Data files
FASTQ_R1="${DATA_DIR}/ERR1050079_1.fastq.gz"
FASTQ_R2="${DATA_DIR}/ERR1050079_2.fastq.gz"
VCF="${DATA_DIR}/HG00731_het_only_chr.vcf.gz"

# Reference - STAR index
STAR_INDEX_REMOTE="/iblm/netapp/data1/external/GRC38/combined/google_cloud/star_index"
LOCAL_SCRATCH="${TMPDIR:-/tmp}"
WORK_OUTPUT_DIR="${LOCAL_SCRATCH}/wasp2rust_fair_indel_${JOB_ID:-local}_${TIMESTAMP}"
OUTPUT_DIR="${WORK_OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}"

SAMPLE="HG00731"
THREADS=8
MAX_INDEL_LEN="${MAX_INDEL_LEN:-10}"
SAME_LOCUS_SLOP="${SAME_LOCUS_SLOP:-10}"

# Prevent hidden oversubscription (numpy/BLAS) inside helper Python steps.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
# For WASP2 parallel unified pipeline: disable per-worker htslib BAM threads unless explicitly overridden.
export WASP2_BAM_THREADS="${WASP2_BAM_THREADS:-0}"

ENABLE_WASP2_TIMING="${ENABLE_WASP2_TIMING:-0}"
if [ -n "${WASP2_TIMING+x}" ]; then
    ENABLE_WASP2_TIMING=1
fi
if [ "${ENABLE_WASP2_TIMING}" = "1" ]; then
    export WASP2_TIMING=1
fi

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
echo "WASP2-Rust FAIR Pipeline Benchmark (SNV+INDEL) - HG00731 RNA-seq" > ${PROFILE_LOG}
echo "==============================================================" >> ${PROFILE_LOG}
echo "Timestamp: ${TIMESTAMP}" >> ${PROFILE_LOG}
echo "Date: $(date)" >> ${PROFILE_LOG}
echo "Sample: ${SAMPLE}" >> ${PROFILE_LOG}
echo "Threads: ${THREADS}" >> ${PROFILE_LOG}
echo "WASP2_TIMING enabled: ${ENABLE_WASP2_TIMING}" >> ${PROFILE_LOG}
echo "WASP2_BAM_THREADS: ${WASP2_BAM_THREADS}" >> ${PROFILE_LOG}
echo "git_sha: ${GIT_SHA}" >> ${PROFILE_LOG}
echo "include_indels: true (max_indel_len=${MAX_INDEL_LEN})" >> ${PROFILE_LOG}
echo "same_locus_slop: ${SAME_LOCUS_SLOP}" >> ${PROFILE_LOG}
echo "" >> ${PROFILE_LOG}

echo "========================================"
echo "WASP2-Rust FAIR Benchmark (SNV+INDEL) on HG00731 RNA-seq"
echo "Timestamp: ${TIMESTAMP}"
echo "========================================"
echo "Start time: $(date)"
echo "FASTQs: ${FASTQ_R1}"
echo "VCF: ${VCF}"
echo "Threads: ${THREADS}"
echo "MAX_INDEL_LEN: ${MAX_INDEL_LEN}"
echo "SAME_LOCUS_SLOP: ${SAME_LOCUS_SLOP}"
echo "WASP2_TIMING enabled: ${ENABLE_WASP2_TIMING}"
echo ""

# Verify Rust module
echo "Verifying Rust module..."
${PYTHON} -c "import wasp2_rust; print('wasp2_rust module: AVAILABLE')" || exit 1
${PYTHON} -c "from wasp2_rust import filter_bam_by_variants_py, unified_make_reads_parallel_py, filter_bam_wasp_with_sidecar; print('All functions: AVAILABLE')" || exit 1
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

echo "Creating variant BED file using vcf_to_bed (SNVs + indels)..."
${PYTHON} -c "
from mapping.intersect_variant_data import vcf_to_bed

vcf_to_bed(
    vcf_file='${VCF}',
    out_bed='${OUTPUT_DIR}/variants.bed',
    samples=['${SAMPLE}'],
    include_indels=True,
    max_indel_len=${MAX_INDEL_LEN},
)
print('BED file created successfully')
"
VARIANT_COUNT=$(wc -l < ${OUTPUT_DIR}/variants.bed)
echo "Variants in BED (SNVs + indels): ${VARIANT_COUNT}"

STAGE_END=$(date +%s.%N)
STAGE_TIME=$(echo "${STAGE_END} - ${STAGE_START}" | bc)
echo "Staging time: ${STAGE_TIME}s (NOT included in benchmark)"
echo "" >> ${PROFILE_LOG}
echo "Staging time: ${STAGE_TIME}s (NOT included in benchmark)" >> ${PROFILE_LOG}

# STAR parameters - SAME as STAR+WASP benchmark (minus WASP-specific options)
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
# STEP 3: WASP2-Rust unified make-reads (INDEL MODE) on to_remap.bam ONLY
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 3: WASP2-Rust unified make-reads (INDEL MODE)"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 3: unified make-reads (INDEL MODE) from to_remap.bam" >> ${PROFILE_LOG}

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
    compression_threads=1,
    compress_output=False,
    indel_mode=True,
    max_indel_size=${MAX_INDEL_LEN},
)

with open('${OUTPUT_DIR}/unified_stats.json', 'w') as f:
    json.dump(dict(stats), f, indent=2)

print(f'Pairs processed: {stats[\"pairs_processed\"]:,}')
print(f'Haplotypes written: {stats[\"haplotypes_written\"]:,}')
" 2>> ${PROFILE_LOG}

STEP3_END=$(date +%s.%N)
STEP3_TIME=$(echo "${STEP3_END} - ${STEP3_START}" | bc)
echo "STEP 3 completed in ${STEP3_TIME} seconds"

${PYTHON} "${WASP2_DIR}/benchmarking/tools/log_unified_stats.py" "${OUTPUT_DIR}/unified_stats.json" \
    --label "rnaseq_indel_fair_step3" >> ${PROFILE_LOG} 2>&1 || echo "WARN: failed to log unified_stats" >> ${PROFILE_LOG}

# Derive pre-remap overlap type counts (pairs) from unified_stats.json
TYPE_PRE=$(${PYTHON} - <<PY
import json
d=json.load(open("${OUTPUT_DIR}/unified_stats.json"))
snv_only = int(d.get("pairs_with_snvs_only", 0))
indel_only = int(d.get("pairs_with_indels_only", 0))
both = int(d.get("pairs_with_snvs_and_indels", 0))
print(f"{snv_only} {indel_only} {both}")
PY
)
read SNV_ONLY_PAIRS_PRE INDEL_ONLY_PAIRS_PRE BOTH_PAIRS_PRE <<< "${TYPE_PRE}"
echo "Overlap type pairs PRE-filter:"
echo "  SNV-only pairs:   ${SNV_ONLY_PAIRS_PRE}"
echo "  INDEL-only pairs: ${INDEL_ONLY_PAIRS_PRE}"
echo "  BOTH pairs:       ${BOTH_PAIRS_PRE}"

# Sidecar expected positions (written alongside remap_r1.fq)
SIDECAR="${OUTPUT_DIR}/remap_r1.fq.expected_positions.tsv"
if [ ! -f "${SIDECAR}" ]; then
    echo "ERROR: Expected sidecar not found at ${SIDECAR}" >&2
    ls -1 "${OUTPUT_DIR}"/*expected_positions.tsv 2>/dev/null || true
    exit 1
fi

# Count reads to remap (pairs)
REMAP_COUNT=$(wc -l < ${OUTPUT_DIR}/remap_r1.fq)
REMAP_PAIRS=$((REMAP_COUNT / 4))
echo "Reads to remap: ${REMAP_PAIRS} pairs"

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

# Clean up uncompressed FASTQs (keep sidecar)
rm -f ${OUTPUT_DIR}/remap_r1.fq ${OUTPUT_DIR}/remap_r2.fq

STEP4_END=$(date +%s.%N)
STEP4_TIME=$(echo "${STEP4_END} - ${STEP4_START}" | bc)
echo "STEP 4 completed in ${STEP4_TIME} seconds"

# -----------------------------------------------------------------------------
# STEP 5: WASP2-Rust filter-remapped (with expected positions sidecar)
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 5: WASP2-Rust filter-remapped (CIGAR-aware expected positions)"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 5: filter_bam_wasp_with_sidecar" >> ${PROFILE_LOG}

STEP5_START=$(date +%s.%N)

${TIME} -v ${PYTHON} -c "
from wasp2_rust import filter_bam_wasp_with_sidecar

kept, removed_moved, removed_missing = filter_bam_wasp_with_sidecar(
    '${OUTPUT_DIR}/to_remap.bam',
    '${OUTPUT_DIR}/remapped.bam',
    '${OUTPUT_DIR}/remap_keep.bam',
    threads=${THREADS},
    same_locus_slop=${SAME_LOCUS_SLOP},
    expected_sidecar='${SIDECAR}',
)

print(f'Kept: {kept}, Removed (moved): {removed_moved}, Removed (missing): {removed_missing}')
" 2>> ${PROFILE_LOG}

${SAMTOOLS} index ${OUTPUT_DIR}/remap_keep.bam

STEP5_END=$(date +%s.%N)
STEP5_TIME=$(echo "${STEP5_END} - ${STEP5_START}" | bc)
echo "STEP 5 completed in ${STEP5_TIME} seconds"

REMAP_KEEP_READS=$(${SAMTOOLS} view -c ${OUTPUT_DIR}/remap_keep.bam)
echo "Remap keep reads: ${REMAP_KEEP_READS}"

# Derive post-filter overlap type counts (pairs) using sidecar mask + remap_keep.bam
TYPE_POST=$(${PYTHON} - <<PY
import pysam
from collections import defaultdict

sidecar_path = "${SIDECAR}"
mask_by_name = {}
with open(sidecar_path) as f:
    for line in f:
        parts = line.rstrip("\\n").split("\\t")
        if len(parts) < 6:
            continue
        q = parts[0]
        if "_WASP_" in q:
            base = q.split("_WASP_")[0]
        else:
            base = q.split("/")[0]
        try:
            mask = int(parts[5])
        except Exception:
            mask = 0
        if mask:
            mask_by_name[base] = mask

counts = defaultdict(int)
bam = pysam.AlignmentFile("${OUTPUT_DIR}/remap_keep.bam", "rb")
for rec in bam.fetch(until_eof=True):
    base = rec.query_name
    m = mask_by_name.get(base, 0)
    if m == 1:
        counts["snv"] += 1
    elif m == 2:
        counts["indel"] += 1
    elif m == 3:
        counts["both"] += 1

snv_pairs = counts["snv"] // 2
indel_pairs = counts["indel"] // 2
both_pairs = counts["both"] // 2
print(f"{snv_pairs} {indel_pairs} {both_pairs}")
PY
)
read SNV_ONLY_PAIRS_POST INDEL_ONLY_PAIRS_POST BOTH_PAIRS_POST <<< "${TYPE_POST}"
echo "Overlap type pairs POST-filter:"
echo "  SNV-only pairs:   ${SNV_ONLY_PAIRS_POST}"
echo "  INDEL-only pairs: ${INDEL_ONLY_PAIRS_POST}"
echo "  BOTH pairs:       ${BOTH_PAIRS_POST}"

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

FINAL_READS=$(${SAMTOOLS} view -c ${OUTPUT_DIR}/wasp_filtered.bam)
PASS_RATE=$(echo "scale=2; ${FINAL_READS} * 100 / ${ORIGINAL_READS}" | bc)

WASP_ONLY=$(echo "${STEP2_TIME} + ${STEP3_TIME} + ${STEP5_TIME} + ${STEP6_TIME}" | bc)
STAR_ONLY=$(echo "${STEP1_TIME} + ${STEP4_TIME}" | bc)

VARIANT_OVERLAP_PRE_TOTAL=$((SNV_ONLY_PAIRS_PRE + INDEL_ONLY_PAIRS_PRE + BOTH_PAIRS_PRE))
VARIANT_OVERLAP_POST_TOTAL=$((SNV_ONLY_PAIRS_POST + INDEL_ONLY_PAIRS_POST + BOTH_PAIRS_POST))
if [ "${VARIANT_OVERLAP_PRE_TOTAL}" -gt 0 ]; then
    RETENTION_OVERLAP=$(echo "scale=2; ${VARIANT_OVERLAP_POST_TOTAL} * 100 / ${VARIANT_OVERLAP_PRE_TOTAL}" | bc)
else
    RETENTION_OVERLAP="0"
fi

echo "========================================"
echo "BENCHMARK RESULTS (FAIR SNV+INDEL)"
echo "========================================"
echo ""
echo "Aggregate times:"
echo "  STAR-only (steps 1+4):            ${STAR_ONLY} s"
echo "  WASP-only (steps 2+3+5+6):        ${WASP_ONLY} s"
echo "  TOTAL time:                       ${TOTAL_TIME} s"
echo ""
echo "Variant-overlap retention:"
echo "  Overlap pairs pre:  ${VARIANT_OVERLAP_PRE_TOTAL}"
echo "  Overlap pairs post: ${VARIANT_OVERLAP_POST_TOTAL}"
echo "  Retention:          ${RETENTION_OVERLAP}%"
echo ""
echo "Read counts:"
echo "  Original reads:       ${ORIGINAL_READS}"
echo "  Keep reads:           ${KEEP_READS}"
echo "  Remap keep reads:     ${REMAP_KEEP_READS}"
echo "  FINAL (merged):       ${FINAL_READS}"
echo "  Pass rate (overall):  ${PASS_RATE}%"
echo ""
echo "Final output directory: ${FINAL_OUTPUT_DIR}"

# Remove genome from shared memory
echo "Removing genome from shared memory..."
${STAR} --genomeDir ${STAR_INDEX} --genomeLoad Remove 2>/dev/null || true

# Cleanup local scratch
rm -rf "${LOCAL_FASTQ_DIR}"

# Save results to JSON (keys align with SNV FAIR JSON where possible)
cat > ${OUTPUT_DIR}/benchmark_results.json << EOF
{
    "timestamp": "${TIMESTAMP}",
    "pipeline": "wasp2rust_fair_indel",
    "git_sha": "${GIT_SHA}",
    "methodology": "LoadAndKeep genome caching, samtools index excluded from STAR timing",
    "sample": "${SAMPLE}",
    "data_type": "RNA-seq",
    "threads": ${THREADS},
    "wasp2_bam_threads": ${WASP2_BAM_THREADS},
    "include_indels": true,
    "indel_mode": true,
    "max_indel_len": ${MAX_INDEL_LEN},
    "same_locus_slop": ${SAME_LOCUS_SLOP},
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
    "variant_count": ${VARIANT_COUNT},
    "snv_only_overlap_pairs_pre": ${SNV_ONLY_PAIRS_PRE},
    "indel_only_overlap_pairs_pre": ${INDEL_ONLY_PAIRS_PRE},
    "snv_indel_overlap_pairs_pre": ${BOTH_PAIRS_PRE},
    "snv_only_overlap_pairs_post": ${SNV_ONLY_PAIRS_POST},
    "indel_only_overlap_pairs_post": ${INDEL_ONLY_PAIRS_POST},
    "snv_indel_overlap_pairs_post": ${BOTH_PAIRS_POST},
    "variant_overlap_retention_percent": ${RETENTION_OVERLAP}
}
EOF

echo ""
echo "Results saved to: ${OUTPUT_DIR}/benchmark_results.json"

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

echo "Final output directory: ${FINAL_OUTPUT_DIR}"
cd /
rm -rf "${WORK_OUTPUT_DIR}"

