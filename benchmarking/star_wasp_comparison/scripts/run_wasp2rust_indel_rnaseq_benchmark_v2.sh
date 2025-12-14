#!/bin/bash
# WASP2-Rust INDEL Benchmark for RNA-seq (HG00731)
# Uses new INDEL mode with coordinated R1/R2 trim combinations
# Sample: HG00731 (1000 Genomes)
# Aligner: STAR (RNA-seq appropriate)
#
#$ -N wasp2rust_indel_rnaseq_v2
#$ -V
#$ -pe iblm 8
#$ -l h_vmem=32G
#$ -j y
#$ -o /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/star_wasp_comparison/logs/
#$ -cwd

set -e

# Increase file descriptor limit for STAR BAM sorting
ulimit -n 10000

# Timestamp
TIMESTAMP=$(date +%Y-%m-%d_%H-%M-%S)
echo "Benchmark timestamp: ${TIMESTAMP}"

# Paths
WASP2_DIR="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp"
BENCHMARK_DIR="${WASP2_DIR}/benchmarking/star_wasp_comparison"
DATA_DIR="${BENCHMARK_DIR}/data"
FINAL_OUTPUT_DIR="${BENCHMARK_DIR}/results/wasp2rust_indel_rnaseq_v2_${TIMESTAMP}"

# Input data
FASTQ_R1="${DATA_DIR}/ERR1050079_1.fastq.gz"
FASTQ_R2="${DATA_DIR}/ERR1050079_2.fastq.gz"
VCF="${DATA_DIR}/HG00731_het_only_chr.vcf.gz"

# Reference - STAR index
STAR_INDEX_REMOTE="/iblm/netapp/data1/external/GRC38/combined/google_cloud/star_index"
LOCAL_SCRATCH="${TMPDIR:-/tmp}"
WORK_OUTPUT_DIR="${LOCAL_SCRATCH}/wasp2rust_indel_rnaseq_v2_${JOB_ID:-local}_${TIMESTAMP}"
OUTPUT_DIR="${WORK_OUTPUT_DIR}"

SAMPLE="HG00731"
THREADS=8
COMPRESSION_THREADS="${COMPRESSION_THREADS:-1}"
MAX_INDEL_LEN=10  # Maximum INDEL length to include
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

# Set PYTHONPATH for WASP2 src AND Rust module
export PYTHONPATH="/iblm/netapp/home/jjaureguy/mambaforge/lib/python3.10/site-packages:${WASP2_DIR}/src:${WASP2_DIR}/rust:${PYTHONPATH}"
# Ensure conda OpenSSL is used (pysam/sidecar)
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH}"

mkdir -p "${OUTPUT_DIR}"
cd "${OUTPUT_DIR}"

# Provenance
GIT_SHA=$(git -C "${WASP2_DIR}" rev-parse HEAD 2>/dev/null || echo "unknown")

# Log file
PROFILE_LOG="${OUTPUT_DIR}/profile.log"
echo "WASP2-Rust INDEL Benchmark Profile (v2 - coordinated trims)" > ${PROFILE_LOG}
echo "Sample: HG00731 RNA-seq" >> ${PROFILE_LOG}
echo "Timestamp: ${TIMESTAMP}" >> ${PROFILE_LOG}
echo "Threads: ${THREADS}" >> ${PROFILE_LOG}
echo "Compression threads: ${COMPRESSION_THREADS}" >> ${PROFILE_LOG}
echo "WASP2_TIMING enabled: ${ENABLE_WASP2_TIMING}" >> ${PROFILE_LOG}
echo "INDEL Support: ENABLED (max ${MAX_INDEL_LEN}bp)" >> ${PROFILE_LOG}
echo "Pipeline: run_make_remap_reads_unified with indel_mode=True" >> ${PROFILE_LOG}

echo "========================================"
echo "WASP2-Rust INDEL RNA-seq Benchmark v2"
echo "Using coordinated R1/R2 trim combinations"
echo "Sample: HG00731"
echo "Timestamp: ${TIMESTAMP}"
echo "INDEL Support: ENABLED (max ${MAX_INDEL_LEN}bp)"
echo "========================================"
echo "Start time: $(date)"
echo "FASTQs: ${FASTQ_R1}"
echo "VCF: ${VCF}"
echo "Threads: ${THREADS}"
echo "Compression threads: ${COMPRESSION_THREADS}"
echo "WASP2_TIMING enabled: ${ENABLE_WASP2_TIMING}"

# Verify modules
echo "Verifying modules..."
${PYTHON} -c "from mapping.run_mapping import run_make_remap_reads_unified; print('run_make_remap_reads_unified: AVAILABLE')" || {
    echo "ERROR: run_make_remap_reads_unified not available"
    exit 1
}
${PYTHON} -c "from wasp2_rust import filter_bam_wasp; print('filter_bam_wasp: AVAILABLE')" || {
    echo "ERROR: filter_bam_wasp not available"
    exit 1
}
${PYTHON} -c "import wasp2_rust; import inspect; sig = inspect.signature(wasp2_rust.unified_make_reads_parallel_py); assert 'indel_mode' in sig.parameters; print('indel_mode parameter: AVAILABLE')" || {
    echo "ERROR: indel_mode parameter not available in wasp2_rust"
    exit 1
}
echo ""

# -----------------------------------------------------------------------------
# PRE-BENCHMARK: Stage STAR index to local scratch
# -----------------------------------------------------------------------------
echo "========================================"
echo "PRE-BENCHMARK: Staging STAR index"
echo "========================================"
LOCAL_STAR_INDEX="${LOCAL_SCRATCH}/wasp2_star_index_${JOB_ID:-local}"

if [ ! -d "${LOCAL_STAR_INDEX}" ]; then
    echo "Copying STAR index to local scratch..."
    mkdir -p "${LOCAL_STAR_INDEX}"
    rsync -a --progress "${STAR_INDEX_REMOTE}/" "${LOCAL_STAR_INDEX}/"
else
    echo "Using existing STAR index at ${LOCAL_STAR_INDEX}"
fi

# -----------------------------------------------------------------------------
# STEP 0: Initial STAR alignment
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 0: Initial STAR alignment"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 0: Initial STAR alignment" >> ${PROFILE_LOG}

STEP0_START=$(date +%s.%N)

STAR_OUT="${OUTPUT_DIR}/star_initial"
mkdir -p ${STAR_OUT}

# STAR parameters matching rasiimwe et al. STAR-WASP benchmark paper
# https://github.com/rasiimwe/STAR-WASP-reduces-reference-bias-in-the-allele-specific-mapping-of-RNA-seq-reads
${TIME} -v ${STAR} \
    --runThreadN ${THREADS} \
    --genomeDir ${LOCAL_STAR_INDEX} \
    --genomeLoad NoSharedMemory \
    --readFilesIn ${FASTQ_R1} ${FASTQ_R2} \
    --readFilesCommand zcat \
    --outFileNamePrefix ${STAR_OUT}/${SAMPLE}_ \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI AS nM NM MD jM jI \
    --alignEndsType EndToEnd \
    --outSAMunmapped Within \
    --outFilterMultimapNmax 1 2>> ${PROFILE_LOG}

STEP0_END=$(date +%s.%N)
STEP0_TIME=$(echo "${STEP0_END} - ${STEP0_START}" | bc)
echo "STEP 0 completed in ${STEP0_TIME} seconds"

ORIGINAL_BAM="${STAR_OUT}/${SAMPLE}_Aligned.sortedByCoord.out.bam"
${SAMTOOLS} index ${ORIGINAL_BAM}
ORIGINAL_READS=$(${SAMTOOLS} view -c ${ORIGINAL_BAM})
echo "Original reads: ${ORIGINAL_READS}"

echo "========================================"
echo "BENCHMARK TIMING STARTS NOW"
echo "========================================"

# -----------------------------------------------------------------------------
# STEP 1: WASP2-Rust unified make-reads with INDEL support
# Uses run_make_remap_reads_unified with include_indels=True
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 1: WASP2-Rust unified make-reads (with INDELs)"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 1: WASP2-Rust unified make-reads (INDEL support)" >> ${PROFILE_LOG}

STEP1_START=$(date +%s.%N)

${TIME} -v ${PYTHON} -c "
from mapping.run_mapping import run_make_remap_reads_unified
import json

stats = run_make_remap_reads_unified(
    bam_file='${ORIGINAL_BAM}',
    variant_file='${VCF}',
    samples='${SAMPLE}',
    out_dir='${OUTPUT_DIR}',
	    include_indels=True,    # INDEL SUPPORT ENABLED
	    max_indel_len=${MAX_INDEL_LEN},
	    threads=${THREADS},
	    compression_threads=${COMPRESSION_THREADS},
	    use_parallel=True
	)

with open('${OUTPUT_DIR}/unified_stats.json', 'w') as f:
    json.dump(stats, f, indent=2)

print(f'Pairs processed: {stats[\"pairs_processed\"]:,}')
print(f'Haplotypes written: {stats[\"haplotypes_written\"]:,}')
" 2>> ${PROFILE_LOG}

STEP1_END=$(date +%s.%N)
STEP1_TIME=$(echo "${STEP1_END} - ${STEP1_START}" | bc)
echo "STEP 1 completed in ${STEP1_TIME} seconds"

# Log unified stats (including timing breakdown)
${PYTHON} "${WASP2_DIR}/benchmarking/tools/log_unified_stats.py" "${OUTPUT_DIR}/unified_stats.json" \
    --label "rnaseq_indel_step1" >> ${PROFILE_LOG} 2>&1 || echo "WARN: failed to log unified_stats" >> ${PROFILE_LOG}

# Count SNV-only vs INDEL-only vs BOTH overlaps (PRE-filter, from unified stats; pairs)
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

# Find FASTQ files (named *_remap_r1.fq.gz by run_make_remap_reads_unified)
R1_READS=$(ls ${OUTPUT_DIR}/*_remap_r1.fq.gz 2>/dev/null | head -1)
R2_READS=$(ls ${OUTPUT_DIR}/*_remap_r2.fq.gz 2>/dev/null | head -1)

if [ -n "${R1_READS}" ]; then
    REMAP_COUNT=$(zcat ${R1_READS} | wc -l)
    REMAP_PAIRS=$((REMAP_COUNT / 4))
    echo "Reads to remap: ${REMAP_PAIRS} pairs"
    echo "SNV_READS_PRE: ${REMAP_PAIRS}"
else
    echo "ERROR: No remap FASTQ files found"
    exit 1
fi

# Sidecar expected positions produced by unified make-reads lives alongside R1_READS
SIDECAR="${R1_READS}.expected_positions.tsv"
if [ ! -f "${SIDECAR}" ]; then
    echo "ERROR: Expected sidecar not found at ${SIDECAR}"
    echo "Available sidecars in ${OUTPUT_DIR}:"
    ls -1 "${OUTPUT_DIR}"/*expected_positions.tsv 2>/dev/null || true
    exit 1
fi

# Create variant count from unified stats
VARIANT_COUNT=$(${PYTHON} -c "import json; d=json.load(open('${OUTPUT_DIR}/unified_stats.json')); print(d.get('variants_used', 'N/A'))" 2>/dev/null || echo "N/A")
echo "Variant count: ${VARIANT_COUNT}"

# -----------------------------------------------------------------------------
# STEP 2: STAR remap flipped reads
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 2: STAR remap flipped reads"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 2: STAR Remap" >> ${PROFILE_LOG}

STEP2_START=$(date +%s.%N)

STAR_REMAP="${OUTPUT_DIR}/star_remap"
mkdir -p ${STAR_REMAP}

${TIME} -v ${STAR} \
    --runThreadN ${THREADS} \
    --genomeDir ${LOCAL_STAR_INDEX} \
    --readFilesIn ${R1_READS} ${R2_READS} \
    --readFilesCommand zcat \
    --outFileNamePrefix ${STAR_REMAP}/${SAMPLE}_remap_ \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI AS nM NM MD \
    --outFilterMultimapNmax 1 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverLmax 0.04 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 2>> ${PROFILE_LOG}

REMAPPED_BAM="${STAR_REMAP}/${SAMPLE}_remap_Aligned.sortedByCoord.out.bam"
${SAMTOOLS} index ${REMAPPED_BAM}

STEP2_END=$(date +%s.%N)
STEP2_TIME=$(echo "${STEP2_END} - ${STEP2_START}" | bc)
echo "STEP 2 completed in ${STEP2_TIME} seconds"

# Check mapping stats from STAR log
if [ -f "${STAR_REMAP}/${SAMPLE}_remap_Log.final.out" ]; then
    echo "STAR remap stats:"
    grep -E "Uniquely mapped|Number of input reads|mapped reads %" ${STAR_REMAP}/${SAMPLE}_remap_Log.final.out || true
fi

# -----------------------------------------------------------------------------
# STEP 3: WASP2-Rust filter_bam_wasp
# filter_bam_wasp strips _WASP_* suffix from remapped reads to match originals
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 3: WASP2-Rust filter_bam_wasp"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 3: WASP2-Rust filter_bam_wasp" >> ${PROFILE_LOG}

STEP3_START=$(date +%s.%N)

${TIME} -v ${PYTHON} -c "
from wasp2_rust import filter_bam_wasp_with_sidecar

# filter_bam_wasp strips _WASP_* suffix from remapped reads
# to match against original read names in ORIGINAL_BAM
kept, removed_moved, removed_missing = filter_bam_wasp_with_sidecar(
    '${ORIGINAL_BAM}',
    '${REMAPPED_BAM}',
    '${OUTPUT_DIR}/remap_keep.bam',
    threads=${THREADS},
    same_locus_slop=10,
    expected_sidecar='${SIDECAR}'
)

print(f'Kept: {kept}, Removed (moved): {removed_moved}, Removed (missing): {removed_missing}')
" 2>> ${PROFILE_LOG}

${SAMTOOLS} index ${OUTPUT_DIR}/remap_keep.bam

STEP3_END=$(date +%s.%N)
STEP3_TIME=$(echo "${STEP3_END} - ${STEP3_START}" | bc)
echo "STEP 3 completed in ${STEP3_TIME} seconds"

REMAP_KEEP_READS=$(${SAMTOOLS} view -c ${OUTPUT_DIR}/remap_keep.bam 2>/dev/null || echo "0")
echo "Remap keep reads: ${REMAP_KEEP_READS}"

# Count SNV-only vs INDEL-only vs BOTH overlaps (POST-filter, from sidecar + remap_keep.bam)
TYPE_POST=$(${PYTHON} - <<PY
import pysam
from collections import defaultdict
sidecar_path = "${SIDECAR}"
mask_by_name = {}
with open(sidecar_path) as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
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
# Count SNV-overlapping reads (POST-filter)
# -----------------------------------------------------------------------------
echo "========================================"
echo "Count SNV-overlapping reads"
echo "========================================"

# SNV-overlapping read pairs POST-filter (remap_keep has the kept reads from remapping)
SNV_READS_POST=$((REMAP_KEEP_READS / 2))
echo "SNV_READS_POST: ${SNV_READS_POST}"

# Calculate pass rate for SNV-overlapping reads
if [ "${REMAP_PAIRS}" -gt 0 ]; then
    SNV_PASS_RATE=$(echo "scale=2; ${SNV_READS_POST} * 100 / ${REMAP_PAIRS}" | bc)
else
    SNV_PASS_RATE="0"
fi
echo "SNV_PASS_RATE: ${SNV_PASS_RATE}%"

# -----------------------------------------------------------------------------
# SUMMARY
# -----------------------------------------------------------------------------
TOTAL_TIME=$(echo "${STEP0_TIME} + ${STEP1_TIME} + ${STEP2_TIME} + ${STEP3_TIME}" | bc)
WASP_ONLY=$(echo "${STEP1_TIME} + ${STEP3_TIME}" | bc)

if [ "${ORIGINAL_READS}" -gt 0 ]; then
    OVERALL_PASS_RATE=$(echo "scale=2; ${REMAP_KEEP_READS} * 100 / ${ORIGINAL_READS}" | bc)
else
    OVERALL_PASS_RATE="0"
fi

echo "========================================"
echo "BENCHMARK RESULTS (WASP2-Rust INDEL v2)"
echo "========================================"
echo ""
echo "Individual step times:"
echo "  STEP 0 (STAR alignment):              ${STEP0_TIME} s"
echo "  STEP 1 (unified-make-reads w/INDELs): ${STEP1_TIME} s"
echo "  STEP 2 (STAR remap):                  ${STEP2_TIME} s"
echo "  STEP 3 (filter_bam_wasp):             ${STEP3_TIME} s"
echo ""
echo "WASP-only time (steps 1+3):             ${WASP_ONLY} s"
echo "TOTAL time:                             ${TOTAL_TIME} s"
echo ""
echo "SNV-Overlapping Read Pair Counts:"
echo "  SNV_READS_PRE:  ${REMAP_PAIRS}"
echo "  SNV_READS_POST: ${SNV_READS_POST}"
echo "  SNV_PASS_RATE:  ${SNV_PASS_RATE}%"
echo ""
echo "Overall Read counts:"
echo "  Original reads:       ${ORIGINAL_READS}"
echo "  Remap keep reads:     ${REMAP_KEEP_READS}"
echo "  Overall pass rate:    ${OVERALL_PASS_RATE}%"
echo ""
echo "Variant counts:"
echo "  Total variants:       ${VARIANT_COUNT}"
echo ""
echo "Work directory: ${OUTPUT_DIR}"
echo "Final output directory: ${FINAL_OUTPUT_DIR}"
echo "End time: $(date)"

# Save results to JSON
cat > ${OUTPUT_DIR}/benchmark_results.json << EOF
{
    "timestamp": "${TIMESTAMP}",
    "pipeline": "wasp2rust_indel_rnaseq_v2",
    "git_sha": "${GIT_SHA}",
    "method": "run_make_remap_reads_unified with indel_mode=True",
    "sample": "${SAMPLE}",
    "data_type": "RNA-seq",
    "aligner": "STAR",
    "threads": ${THREADS},
    "include_indels": true,
    "indel_mode": true,
    "cigar_aware_expected_pos": true,
    "sidecar_format_version": 2,
    "max_indel_len": ${MAX_INDEL_LEN},
    "step0_star_align_s": ${STEP0_TIME},
    "step1_unified_make_reads_s": ${STEP1_TIME},
    "step2_star_remap_s": ${STEP2_TIME},
    "step3_filter_wasp_s": ${STEP3_TIME},
    "wasp_only_s": ${WASP_ONLY},
    "total_s": ${TOTAL_TIME},
    "snv_reads_pre": ${REMAP_PAIRS},
    "snv_reads_post": ${SNV_READS_POST},
    "snv_pass_rate_percent": ${SNV_PASS_RATE},
    "snv_only_overlap_pairs_pre": ${SNV_ONLY_PAIRS_PRE},
    "indel_only_overlap_pairs_pre": ${INDEL_ONLY_PAIRS_PRE},
    "snv_indel_overlap_pairs_pre": ${BOTH_PAIRS_PRE},
    "snv_only_overlap_pairs_post": ${SNV_ONLY_PAIRS_POST},
    "indel_only_overlap_pairs_post": ${INDEL_ONLY_PAIRS_POST},
    "snv_indel_overlap_pairs_post": ${BOTH_PAIRS_POST},
    "original_reads": ${ORIGINAL_READS},
    "remap_keep_reads": ${REMAP_KEEP_READS},
    "overall_pass_rate_percent": ${OVERALL_PASS_RATE},
    "variant_count": "${VARIANT_COUNT}"
}
EOF

echo ""
echo "Results saved to: ${OUTPUT_DIR}/benchmark_results.json (work dir)"

# -----------------------------------------------------------------------------
# Copy artifacts back to final output directory (NFS)
# -----------------------------------------------------------------------------
echo "Copying benchmark artifacts to final output directory..."
mkdir -p "${FINAL_OUTPUT_DIR}"
for f in \
    "${OUTPUT_DIR}/benchmark_results.json" \
    "${OUTPUT_DIR}/profile.log" \
    "${OUTPUT_DIR}/unified_stats.json" \
    "${OUTPUT_DIR}/remap_keep.bam" \
    "${OUTPUT_DIR}/remap_keep.bam.bai" \
    "${OUTPUT_DIR}/star_initial/${SAMPLE}_Log.final.out" \
    "${OUTPUT_DIR}/star_remap/${SAMPLE}_remap_Log.final.out"
do
    if [ -f "${f}" ]; then
        rsync -a "${f}" "${FINAL_OUTPUT_DIR}/"
    fi
done
echo "Final output directory: ${FINAL_OUTPUT_DIR}"

# Cleanup STAR index from scratch
echo "Cleaning up STAR index from local scratch..."
rm -rf "${LOCAL_STAR_INDEX}" 2>/dev/null || true

echo "Done!"
