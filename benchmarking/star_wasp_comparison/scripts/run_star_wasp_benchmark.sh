#!/bin/bash
# STAR+WASP Benchmark on HG00731 RNA-seq data
# This runs the STAR+WASP method (Asiimwe & Dobin, bioRxiv 2024) on our hardware
# for a fair comparison with WASP2-Rust.
#
# Methodology matches their GitHub:
# https://github.com/rasiimwe/STAR-WASP-reduces-reference-bias-in-the-allele-specific-mapping-of-RNA-seq-reads
#
#$ -o /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/logs/run_star_wasp_benchmark.sh.o$JOB_ID
#$ -e /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/logs/run_star_wasp_benchmark.sh.e$JOB_ID

set -e

# Increase file descriptor limit for STAR BAM sorting
ulimit -n 10000

# Timestamp for this run
TIMESTAMP=$(date +%Y-%m-%d_%H-%M-%S)
echo "STAR+WASP Benchmark timestamp: ${TIMESTAMP}"

# =============================================================================
# PATHS
# =============================================================================
BENCHMARK_DIR="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/star_wasp_comparison"
DATA_DIR="${BENCHMARK_DIR}/data"
OUTPUT_BASE="${BENCHMARK_DIR}/results"
FINAL_OUTPUT_DIR="${OUTPUT_BASE}/star_wasp_${TIMESTAMP}"

# Input files (remote/NFS)
FASTQ_R1="${DATA_DIR}/ERR1050079_1.fastq.gz"
FASTQ_R2="${DATA_DIR}/ERR1050079_2.fastq.gz"
VCF_GZ="${DATA_DIR}/HG00731_het_only_chr.vcf.gz"
STAR_INDEX_REMOTE="/iblm/netapp/data1/external/GRC38/combined/google_cloud/star_index"

# Local scratch
LOCAL_SCRATCH="${TMPDIR:-/tmp}"
WORK_OUTPUT_DIR="${LOCAL_SCRATCH}/star_wasp_${JOB_ID:-local}_${TIMESTAMP}"
OUTPUT_DIR="${WORK_OUTPUT_DIR}"

mkdir -p "${OUTPUT_DIR}"
echo "Work directory: ${OUTPUT_DIR}"
echo "Final output directory: ${FINAL_OUTPUT_DIR}"

# =============================================================================
# CONFIGURATION (matches paper)
# =============================================================================
THREADS=8
SAMPLE="HG00731"

# Conda environment
source /iblm/netapp/home/jjaureguy/mambaforge/etc/profile.d/conda.sh
conda activate WASP2_dev2

# STAR executable
STAR=$(which STAR)

# GNU time for profiling
TIME="/usr/bin/time"

# Log file
PROFILE_LOG="${OUTPUT_DIR}/profile.log"
echo "STAR+WASP Benchmark" > ${PROFILE_LOG}
echo "===================" >> ${PROFILE_LOG}
echo "Timestamp: ${TIMESTAMP}" >> ${PROFILE_LOG}
echo "Date: $(date)" >> ${PROFILE_LOG}
echo "Sample: ${SAMPLE}" >> ${PROFILE_LOG}
echo "Threads: ${THREADS}" >> ${PROFILE_LOG}
echo "STAR: ${STAR}" >> ${PROFILE_LOG}
echo "STAR version: $(${STAR} --version 2>&1 | head -1)" >> ${PROFILE_LOG}
echo "" >> ${PROFILE_LOG}

# =============================================================================
# PRE-BENCHMARK: Stage files to local scratch (NOT TIMED)
# This matches their setup where data is already on local disk
# =============================================================================
echo "========================================"
echo "PRE-BENCHMARK: Staging files to local scratch"
echo "========================================"

# Stage STAR index
LOCAL_STAR_INDEX="${LOCAL_SCRATCH}/star_wasp_index_${JOB_ID:-local}"
echo "Staging STAR index to: ${LOCAL_STAR_INDEX}"
STAGE_START=$(date +%s.%N)
if [ -d "${LOCAL_STAR_INDEX}" ]; then
    echo "  Local index already exists, reusing..."
else
    mkdir -p "${LOCAL_STAR_INDEX}"
    rsync -a "${STAR_INDEX_REMOTE}/" "${LOCAL_STAR_INDEX}/"
fi
STAGE_INDEX_END=$(date +%s.%N)
STAGE_INDEX_TIME=$(echo "${STAGE_INDEX_END} - ${STAGE_START}" | bc)
echo "  Index staging time: ${STAGE_INDEX_TIME}s"

# Copy FASTQs to local scratch
LOCAL_FASTQ_DIR="${LOCAL_SCRATCH}/star_wasp_fastqs_${TIMESTAMP}"
mkdir -p "${LOCAL_FASTQ_DIR}"
LOCAL_FASTQ_R1="${LOCAL_FASTQ_DIR}/R1.fastq.gz"
LOCAL_FASTQ_R2="${LOCAL_FASTQ_DIR}/R2.fastq.gz"
echo "Copying FASTQs to: ${LOCAL_FASTQ_DIR}"
cp "${FASTQ_R1}" "${LOCAL_FASTQ_R1}"
cp "${FASTQ_R2}" "${LOCAL_FASTQ_R2}"
STAGE_FASTQ_END=$(date +%s.%N)
STAGE_FASTQ_TIME=$(echo "${STAGE_FASTQ_END} - ${STAGE_INDEX_END}" | bc)
echo "  FASTQ copy time: ${STAGE_FASTQ_TIME}s"

# Decompress VCF (STAR --varVCFfile requires uncompressed VCF)
LOCAL_VCF="${LOCAL_SCRATCH}/star_wasp_vcf_${TIMESTAMP}/HG00731_het.vcf"
mkdir -p "$(dirname ${LOCAL_VCF})"
echo "Decompressing VCF to: ${LOCAL_VCF}"
zcat "${VCF_GZ}" > "${LOCAL_VCF}"
STAGE_VCF_END=$(date +%s.%N)
STAGE_VCF_TIME=$(echo "${STAGE_VCF_END} - ${STAGE_FASTQ_END}" | bc)
echo "  VCF decompress time: ${STAGE_VCF_TIME}s"

TOTAL_STAGE_TIME=$(echo "${STAGE_VCF_END} - ${STAGE_START}" | bc)
echo ""
echo "Total staging time: ${TOTAL_STAGE_TIME}s (NOT counted in benchmark)"
echo "" >> ${PROFILE_LOG}
echo "PRE-BENCHMARK STAGING (not timed):" >> ${PROFILE_LOG}
echo "  Index staging: ${STAGE_INDEX_TIME}s" >> ${PROFILE_LOG}
echo "  FASTQ copy: ${STAGE_FASTQ_TIME}s" >> ${PROFILE_LOG}
echo "  VCF decompress: ${STAGE_VCF_TIME}s" >> ${PROFILE_LOG}
echo "  TOTAL: ${TOTAL_STAGE_TIME}s" >> ${PROFILE_LOG}
echo "" >> ${PROFILE_LOG}

# =============================================================================
# BENCHMARK: STAR+WASP (TIMED)
# Single STAR pass with WASP integrated - this is what the paper reports
# =============================================================================
echo ""
echo "========================================"
echo "BENCHMARK: STAR+WASP alignment"
echo "========================================"
echo "Start time: $(date)"
echo ""

cd "${OUTPUT_DIR}"

# STAR+WASP parameters - EXACT match to paper's GitHub
# Key additions vs plain STAR:
#   --waspOutputMode SAMtag  : Enable WASP filtering, add vW tags
#   --varVCFfile             : VCF with het variants for WASP
#   --outSAMattributes ... vA vG vW : Include WASP-specific tags
STARpar="--runThreadN ${THREADS} \
         --genomeDir ${LOCAL_STAR_INDEX} \
         --genomeLoad NoSharedMemory \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMattributes NH HI AS nM NM MD jM jI vA vG vW \
         --waspOutputMode SAMtag \
         --alignEndsType EndToEnd \
         --outSAMunmapped Within \
         --outFilterMultimapNmax 1"

readFiles="--readFilesCommand gunzip -c \
           --readFilesIn ${LOCAL_FASTQ_R1} ${LOCAL_FASTQ_R2}"

varVcf="--varVCFfile ${LOCAL_VCF}"

echo "STAR+WASP Command:" >> ${PROFILE_LOG}
echo "${STAR} ${STARpar} ${readFiles} ${varVcf}" >> ${PROFILE_LOG}
echo "" >> ${PROFILE_LOG}

# Run STAR+WASP with timing
BENCH_START=$(date +%s.%N)
echo "STAR+WASP Timing:" >> ${PROFILE_LOG}
${TIME} -v ${STAR} ${STARpar} ${readFiles} ${varVcf} 2>> ${PROFILE_LOG}
BENCH_END=$(date +%s.%N)

BENCH_TIME=$(echo "${BENCH_END} - ${BENCH_START}" | bc)

# Rename output
mv Aligned.sortedByCoord.out.bam star_wasp_output.bam

echo ""
echo "STAR+WASP completed in ${BENCH_TIME} seconds"
echo "End time: $(date)"

# =============================================================================
# POST-BENCHMARK: Index BAM and extract stats
# =============================================================================
echo ""
echo "========================================"
echo "POST-BENCHMARK: Indexing and stats"
echo "========================================"

samtools index star_wasp_output.bam

# Count WASP-tagged reads
echo "Counting WASP-filtered reads..."
TOTAL_READS=$(samtools view -c star_wasp_output.bam)
WASP_PASS=$(samtools view star_wasp_output.bam | awk '$0 ~ /vW:i:1/' | wc -l)
WASP_FAIL=$(samtools view star_wasp_output.bam | awk '$0 ~ /vW:i:[2-7]/' | wc -l)

echo "  Total aligned reads: ${TOTAL_READS}"
echo "  WASP pass (vW:i:1): ${WASP_PASS}"
echo "  WASP fail (vW:i:2-7): ${WASP_FAIL}"

# =============================================================================
# RESULTS
# =============================================================================
echo ""
echo "========================================"
echo "BENCHMARK RESULTS"
echo "========================================"
echo ""

# Format time
BENCH_MIN=$(echo "${BENCH_TIME} / 60" | bc)
BENCH_SEC=$(printf "%.2f" $(echo "${BENCH_TIME} - ${BENCH_MIN} * 60" | bc))

echo "STAR+WASP wall clock time: ${BENCH_TIME}s (${BENCH_MIN}:${BENCH_SEC})"
echo ""
echo "For comparison (same hardware):"
echo "  Run WASP2-Rust benchmark to get comparable time"
echo ""

# Save results JSON
cat > benchmark_results.json << EOF
{
    "timestamp": "${TIMESTAMP}",
    "method": "STAR+WASP",
    "sample": "${SAMPLE}",
    "threads": ${THREADS},
    "star_version": "$(${STAR} --version 2>&1 | head -1)",
    "wall_clock_s": ${BENCH_TIME},
    "wall_clock_mmss": "${BENCH_MIN}:${BENCH_SEC}",
    "total_aligned_reads": ${TOTAL_READS},
    "wasp_pass_reads": ${WASP_PASS},
    "wasp_fail_reads": ${WASP_FAIL},
    "staging_time_s": ${TOTAL_STAGE_TIME},
    "note": "Staging time NOT included in benchmark (matches paper methodology)"
}
EOF

echo "Results saved to: ${OUTPUT_DIR}/benchmark_results.json"

# Append summary to profile log
echo "" >> ${PROFILE_LOG}
echo "=======================================" >> ${PROFILE_LOG}
echo "SUMMARY" >> ${PROFILE_LOG}
echo "=======================================" >> ${PROFILE_LOG}
echo "STAR+WASP wall clock: ${BENCH_TIME}s (${BENCH_MIN}:${BENCH_SEC})" >> ${PROFILE_LOG}
echo "Total aligned reads: ${TOTAL_READS}" >> ${PROFILE_LOG}
echo "WASP pass (vW:i:1): ${WASP_PASS}" >> ${PROFILE_LOG}
echo "WASP fail (vW:i:2-7): ${WASP_FAIL}" >> ${PROFILE_LOG}

# =============================================================================
# COPY BACK (NFS)
# =============================================================================
echo ""
echo "Copying benchmark artifacts to final output directory..."
mkdir -p "${FINAL_OUTPUT_DIR}"
for f in \
    "${OUTPUT_DIR}/benchmark_results.json" \
    "${PROFILE_LOG}" \
    "${OUTPUT_DIR}/Log.final.out"
do
    if [ -f "${f}" ]; then
        rsync -a "${f}" "${FINAL_OUTPUT_DIR}/"
    fi
done

# Optional: keep BAM outputs (large)
if [ "${KEEP_BAMS:-0}" = "1" ]; then
    for f in "${OUTPUT_DIR}/star_wasp_output.bam" "${OUTPUT_DIR}/star_wasp_output.bam.bai"; do
        if [ -f "${f}" ]; then
            rsync -a "${f}" "${FINAL_OUTPUT_DIR}/"
        fi
    done
fi

echo "Final output directory: ${FINAL_OUTPUT_DIR}"

# =============================================================================
# CLEANUP
# =============================================================================
echo ""
echo "Cleaning up local scratch files..."
cd /
rm -rf "${LOCAL_FASTQ_DIR}"
rm -rf "$(dirname ${LOCAL_VCF})"
rm -rf "${LOCAL_STAR_INDEX}"
rm -rf "${WORK_OUTPUT_DIR}"
echo "Done."

echo ""
echo "========================================"
echo "STAR+WASP Benchmark Complete"
echo "Final output: ${FINAL_OUTPUT_DIR}"
echo "========================================"
