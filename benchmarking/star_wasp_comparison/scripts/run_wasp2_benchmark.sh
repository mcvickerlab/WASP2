#!/bin/bash
# WASP2-Rust benchmark on HG00731 RNA-seq data
# Matching exact STAR parameters from STAR+WASP paper (Asiimwe & Dobin 2024)
#
# This benchmarks WASP2-Rust against their published WASP1 time (25:41 = 1541s)
#
# Their pipeline steps:
#   1. STAR initial alignment
#   2. find_intersecting_snps.py (Python) -> WASP2 make-reads (Rust)
#   3. STAR remap
#   4. filter_remapped_reads.py (Python) -> WASP2 filter-remapped (Rust)
#
# Memory profiling: Using /usr/bin/time -v to match their methodology

set -e

# Increase file descriptor limit for STAR BAM sorting (same as original paper)
ulimit -n 10000

# Paths
BENCHMARK_DIR="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/star_wasp_comparison"
DATA_DIR="${BENCHMARK_DIR}/data"
OUTPUT_DIR="${BENCHMARK_DIR}/results"
WASP2_DIR="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp"

# Data files
FASTQ_R1="${DATA_DIR}/ERR1050079_1.fastq.gz"
FASTQ_R2="${DATA_DIR}/ERR1050079_2.fastq.gz"
VCF="${DATA_DIR}/HG00731_het_only_chr.vcf.gz"

# Reference - using google_cloud index (built with STAR 2.7.4a, compatible with 2.7.11b)
# GRCh38 + Gencode v26
STAR_INDEX="/iblm/netapp/data1/external/GRC38/combined/google_cloud/star_index"

# Conda environment
source /iblm/netapp/home/jjaureguy/mambaforge/etc/profile.d/conda.sh
conda activate WASP2_dev2

# STAR executable (from conda)
STAR=$(which STAR)

# Threads (match their benchmark: 8 threads)
THREADS=8

# Sample name
SAMPLE="HG00731"

# GNU time for memory profiling (same as paper)
TIME="/usr/bin/time"

mkdir -p ${OUTPUT_DIR}
cd ${OUTPUT_DIR}

# Create working directory
WORK_DIR="${OUTPUT_DIR}/wasp2_run"
rm -rf ${WORK_DIR}
mkdir -p ${WORK_DIR}
cd ${WORK_DIR}

# Log file for detailed profiling
PROFILE_LOG="${OUTPUT_DIR}/HG00731_wasp2_profile.log"
echo "WASP2-Rust Benchmark Profile Log" > ${PROFILE_LOG}
echo "================================" >> ${PROFILE_LOG}
echo "Date: $(date)" >> ${PROFILE_LOG}
echo "Sample: HG00731" >> ${PROFILE_LOG}
echo "Threads: ${THREADS}" >> ${PROFILE_LOG}
echo "STAR: ${STAR}" >> ${PROFILE_LOG}
echo "" >> ${PROFILE_LOG}

echo "========================================"
echo "WASP2-Rust Benchmark on HG00731"
echo "========================================"
echo "Start time: $(date)"
echo "STAR: ${STAR}"
echo "STAR version: $(${STAR} --version 2>&1 | head -1)"
echo "STAR index: ${STAR_INDEX}"
echo "Threads: ${THREADS}"
echo "WASP2 features: --phased (SNPs only)"
echo ""

# Record total wall clock time
TOTAL_START=$(date +%s.%N)

# STAR parameters - EXACT match to paper
STARpar="--runThreadN ${THREADS} \
         --genomeDir ${STAR_INDEX} \
         --genomeLoad NoSharedMemory \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMattributes NH HI AS nM NM MD jM jI \
         --alignEndsType EndToEnd \
         --outSAMunmapped Within \
         --outFilterMultimapNmax 1"

readFiles="--readFilesCommand gunzip -c \
           --readFilesIn ${FASTQ_R1} ${FASTQ_R2}"

# -----------------------------------------------------------------------------
# STEP 1: Initial STAR alignment (SAME as their pipeline)
# Using EXACT same parameters from STAR+WASP paper
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 1: Initial STAR alignment"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 1: STAR Initial Alignment" >> ${PROFILE_LOG}
echo "------------------------------" >> ${PROFILE_LOG}

STEP1_START=$(date +%s.%N)

# Run with /usr/bin/time -v for memory profiling
${TIME} -v ${STAR} ${STARpar} ${readFiles} 2>> ${PROFILE_LOG}

mv Aligned.sortedByCoord.out.bam A_sorted.bam
samtools index A_sorted.bam

STEP1_END=$(date +%s.%N)
STEP1_TIME=$(echo "${STEP1_END} - ${STEP1_START}" | bc)
echo "STEP 1 completed in ${STEP1_TIME} seconds"
echo ""

# -----------------------------------------------------------------------------
# STEP 2: WASP2 make-reads (replaces find_intersecting_snps.py)
# This is the Rust-optimized version with INDEL support
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 2: WASP2 make-reads (Rust + INDELs)"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 2: WASP2 make-reads (Rust)" >> ${PROFILE_LOG}
echo "-------------------------------" >> ${PROFILE_LOG}

STEP2_START=$(date +%s.%N)

# Run with /usr/bin/time -v for memory profiling
# Use PYTHONPATH to ensure proper module import
${TIME} -v env PYTHONPATH=${WASP2_DIR}/src python -m mapping make-reads \
    A_sorted.bam \
    ${VCF} \
    -t ${THREADS} \
    -s ${SAMPLE} \
    --out ./ \
    --paired \
    --phased \
    --temp_loc ./ 2>> ${PROFILE_LOG}

STEP2_END=$(date +%s.%N)
STEP2_TIME=$(echo "${STEP2_END} - ${STEP2_START}" | bc)
echo "STEP 2 completed in ${STEP2_TIME} seconds"
echo ""

# -----------------------------------------------------------------------------
# STEP 3: STAR remap flipped reads (SAME as their pipeline)
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 3: STAR remap flipped reads"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 3: STAR Remap" >> ${PROFILE_LOG}
echo "------------------" >> ${PROFILE_LOG}

STEP3_START=$(date +%s.%N)

# Find the swapped allele FASTQ files
REMAP_R1=$(ls *_swapped_alleles_r1.fq 2>/dev/null | head -1)
REMAP_R2=$(ls *_swapped_alleles_r2.fq 2>/dev/null | head -1)

if [ -z "${REMAP_R1}" ] || [ -z "${REMAP_R2}" ]; then
    echo "ERROR: Swapped allele FASTQ files not found"
    ls -la
    exit 1
fi

echo "Remapping files: ${REMAP_R1}, ${REMAP_R2}"
REMAP_READS=$(wc -l < ${REMAP_R1})
echo "Reads to remap: $((REMAP_READS / 4)) read pairs"

# Run with /usr/bin/time -v for memory profiling
${TIME} -v ${STAR} ${STARpar} \
    --readFilesIn ${REMAP_R1} ${REMAP_R2} 2>> ${PROFILE_LOG}

mv Aligned.sortedByCoord.out.bam remapped.bam
samtools index remapped.bam

STEP3_END=$(date +%s.%N)
STEP3_TIME=$(echo "${STEP3_END} - ${STEP3_START}" | bc)
echo "STEP 3 completed in ${STEP3_TIME} seconds"
echo ""

# -----------------------------------------------------------------------------
# STEP 4: WASP2 filter-remapped (replaces filter_remapped_reads.py)
# This is the Rust-optimized version
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 4: WASP2 filter-remapped (Rust)"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 4: WASP2 filter-remapped (Rust)" >> ${PROFILE_LOG}
echo "------------------------------------" >> ${PROFILE_LOG}

STEP4_START=$(date +%s.%N)

# Find the JSON file
JSON_FILE=$(ls *_wasp_data_files.json 2>/dev/null | head -1)

if [ -z "${JSON_FILE}" ]; then
    echo "ERROR: WASP JSON file not found"
    ls -la
    exit 1
fi

# Run with /usr/bin/time -v for memory profiling
# Use PYTHONPATH to ensure proper module import
${TIME} -v env PYTHONPATH=${WASP2_DIR}/src python -m mapping filter-remapped \
    remapped.bam \
    --json ${JSON_FILE} 2>> ${PROFILE_LOG}

STEP4_END=$(date +%s.%N)
STEP4_TIME=$(echo "${STEP4_END} - ${STEP4_START}" | bc)
echo "STEP 4 completed in ${STEP4_TIME} seconds"
echo ""

# -----------------------------------------------------------------------------
# SUMMARY
# -----------------------------------------------------------------------------
TOTAL_END=$(date +%s.%N)
TOTAL_TIME=$(echo "${TOTAL_END} - ${TOTAL_START}" | bc)

# Convert to mm:ss format
TOTAL_MIN=$(echo "${TOTAL_TIME} / 60" | bc)
TOTAL_SEC=$(echo "${TOTAL_TIME} - (${TOTAL_MIN} * 60)" | bc)

echo "========================================"
echo "BENCHMARK RESULTS"
echo "========================================"
echo ""
echo "Individual step times:"
echo "  STEP 1 (STAR initial):        ${STEP1_TIME} s"
echo "  STEP 2 (WASP2 make-reads):    ${STEP2_TIME} s"
echo "  STEP 3 (STAR remap):          ${STEP3_TIME} s"
echo "  STEP 4 (WASP2 filter):        ${STEP4_TIME} s"
echo ""
echo "TOTAL WASP2-Rust time: ${TOTAL_TIME} s (${TOTAL_MIN}:${TOTAL_SEC})"
echo ""
echo "Published comparison times (HG00731, 56M reads, 8 threads):"
echo "  WASP1 (Python):    1541.13 s (25:41)"
echo "  STAR+WASP (C++):    331.67 s ( 5:32)"
echo ""

# Calculate speedups
WASP1_TIME=1541.13
STAR_WASP_TIME=331.67
SPEEDUP_VS_WASP1=$(echo "scale=2; ${WASP1_TIME} / ${TOTAL_TIME}" | bc)
SPEEDUP_VS_STAR_WASP=$(echo "scale=2; ${STAR_WASP_TIME} / ${TOTAL_TIME}" | bc)

echo "WASP2-Rust speedup:"
echo "  vs WASP1:      ${SPEEDUP_VS_WASP1}x"
echo "  vs STAR+WASP:  ${SPEEDUP_VS_STAR_WASP}x"
echo ""
echo "End time: $(date)"
echo ""
echo "Detailed profile log: ${PROFILE_LOG}"

# Save results to TSV (matching paper format)
RESULTS_FILE="${OUTPUT_DIR}/HG00731_wasp2_timing.tsv"
echo -e "Sample\tMethod\tThreads\tStep\tTime_seconds\tTime_mmss" > ${RESULTS_FILE}
echo -e "HG00731\tWASP2-Rust\t${THREADS}\tSTAR_initial\t${STEP1_TIME}\t" >> ${RESULTS_FILE}
echo -e "HG00731\tWASP2-Rust\t${THREADS}\tmake_reads\t${STEP2_TIME}\t" >> ${RESULTS_FILE}
echo -e "HG00731\tWASP2-Rust\t${THREADS}\tSTAR_remap\t${STEP3_TIME}\t" >> ${RESULTS_FILE}
echo -e "HG00731\tWASP2-Rust\t${THREADS}\tfilter_remapped\t${STEP4_TIME}\t" >> ${RESULTS_FILE}
echo -e "HG00731\tWASP2-Rust\t${THREADS}\tTOTAL\t${TOTAL_TIME}\t${TOTAL_MIN}:${TOTAL_SEC}" >> ${RESULTS_FILE}
echo -e "HG00731\tWASP1\t8\tTOTAL\t1541.13\t25:41" >> ${RESULTS_FILE}
echo -e "HG00731\tSTAR+WASP\t8\tTOTAL\t331.67\t5:32" >> ${RESULTS_FILE}

echo "Results saved to: ${RESULTS_FILE}"

# Append summary to profile log
echo "" >> ${PROFILE_LOG}
echo "=======================================" >> ${PROFILE_LOG}
echo "SUMMARY" >> ${PROFILE_LOG}
echo "=======================================" >> ${PROFILE_LOG}
echo "STEP 1 (STAR initial):     ${STEP1_TIME} s" >> ${PROFILE_LOG}
echo "STEP 2 (WASP2 make-reads): ${STEP2_TIME} s" >> ${PROFILE_LOG}
echo "STEP 3 (STAR remap):       ${STEP3_TIME} s" >> ${PROFILE_LOG}
echo "STEP 4 (WASP2 filter):     ${STEP4_TIME} s" >> ${PROFILE_LOG}
echo "TOTAL:                     ${TOTAL_TIME} s" >> ${PROFILE_LOG}
echo "" >> ${PROFILE_LOG}
echo "Speedup vs WASP1:     ${SPEEDUP_VS_WASP1}x" >> ${PROFILE_LOG}
echo "Speedup vs STAR+WASP: ${SPEEDUP_VS_STAR_WASP}x" >> ${PROFILE_LOG}
