#!/bin/bash
# WASP1 (Python) benchmark on HG00731 RNA-seq data
# Full pipeline: STAR + WASP1 (find_intersecting_snps + filter_remapped_reads)
#
# This provides a fair comparison point for WASP2-Rust unified pipeline.
#
#$ -N wasp1_benchmark
#$ -V
#$ -pe iblm 8
#$ -l h_vmem=32G
#$ -j y
#$ -o /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/logs/
#$ -cwd

# Pipeline:
#   0. (Pre-benchmark) snp2h5 VCF -> HDF5 conversion
#   1. STAR initial alignment
#   2. WASP1 find_intersecting_snps.py
#   3. STAR remap
#   4. WASP1 filter_remapped_reads.py
#   5. Merge BAMs
#
# This is the original WASP1 Python pipeline for fair comparison.

set -e

# Increase file descriptor limit for STAR BAM sorting
ulimit -n 10000

# Timestamp for this run
TIMESTAMP=$(date +%Y-%m-%d_%H-%M-%S)
echo "Benchmark timestamp: ${TIMESTAMP}"

# Paths
BENCHMARK_DIR="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/star_wasp_comparison"
DATA_DIR="${BENCHMARK_DIR}/data"
OUTPUT_BASE="${BENCHMARK_DIR}/results"
OUTPUT_DIR="${OUTPUT_BASE}/wasp1_${TIMESTAMP}"

# WASP1 source (local copy in our repo)
WASP1_DIR="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/wasp1_source"
SNP2H5="${WASP1_DIR}/snp2h5/snp2h5"
FIND_INTERSECTING="${WASP1_DIR}/mapping/find_intersecting_snps.py"
FILTER_REMAPPED="${WASP1_DIR}/mapping/filter_remapped_reads.py"

# Create timestamped output directory
mkdir -p "${OUTPUT_DIR}"
echo "Output directory: ${OUTPUT_DIR}"

# Data files
FASTQ_R1="${DATA_DIR}/ERR1050079_1.fastq.gz"
FASTQ_R2="${DATA_DIR}/ERR1050079_2.fastq.gz"
VCF="${DATA_DIR}/HG00731_het_only_chr.vcf.gz"

# Reference - STAR index (GRCh38 + Gencode v26)
STAR_INDEX_REMOTE="/iblm/netapp/data1/external/GRC38/combined/google_cloud/star_index"
LOCAL_SCRATCH="${TMPDIR:-/tmp}"

# Chromsizes for snp2h5
CHROM_SIZES="${WASP1_DIR}/hg38_chromsizes.txt"

# Sample name
SAMPLE="HG00731"

# Threads
THREADS=8

# Conda environment
source /iblm/netapp/home/jjaureguy/mambaforge/etc/profile.d/conda.sh
conda activate WASP2_dev2

# STAR executable
STAR=$(which STAR)

# GNU time for memory profiling
TIME="/usr/bin/time"

cd "${OUTPUT_DIR}"

# Log file for detailed profiling
PROFILE_LOG="${OUTPUT_DIR}/profile.log"
echo "WASP1 (Python) Pipeline Benchmark" > ${PROFILE_LOG}
echo "==================================" >> ${PROFILE_LOG}
echo "Timestamp: ${TIMESTAMP}" >> ${PROFILE_LOG}
echo "Date: $(date)" >> ${PROFILE_LOG}
echo "Sample: ${SAMPLE}" >> ${PROFILE_LOG}
echo "Threads: ${THREADS}" >> ${PROFILE_LOG}
echo "STAR: ${STAR}" >> ${PROFILE_LOG}
echo "WASP1: ${WASP1_DIR}" >> ${PROFILE_LOG}
echo "" >> ${PROFILE_LOG}

echo "========================================"
echo "WASP1 (Python) Benchmark on HG00731"
echo "Timestamp: ${TIMESTAMP}"
echo "========================================"
echo "Start time: $(date)"
echo "STAR: ${STAR}"
echo "STAR version: $(${STAR} --version 2>&1 | head -1)"
echo "WASP1 directory: ${WASP1_DIR}"
echo "Threads: ${THREADS}"
echo ""

# -----------------------------------------------------------------------------
# PRE-BENCHMARK: Stage data to local scratch (NOT timed)
# -----------------------------------------------------------------------------
echo "========================================"
echo "PRE-BENCHMARK: Staging data to local scratch"
echo "========================================"

# Stage STAR index
LOCAL_STAR_INDEX="${LOCAL_SCRATCH}/wasp1_star_index_${JOB_ID:-local}"
echo "Staging STAR index to: ${LOCAL_STAR_INDEX}"
STAGE_START=$(date +%s.%N)
if [ -d "${LOCAL_STAR_INDEX}" ]; then
    echo "Local STAR index already exists, reusing..."
else
    mkdir -p "${LOCAL_STAR_INDEX}"
    rsync -a "${STAR_INDEX_REMOTE}/" "${LOCAL_STAR_INDEX}/"
fi
STAGE_END=$(date +%s.%N)
STAGE_TIME=$(echo "${STAGE_END} - ${STAGE_START}" | bc)
echo "Index staging time: ${STAGE_TIME}s (NOT counted in benchmark)"

# Copy FASTQs to local scratch
LOCAL_FASTQ_DIR="${LOCAL_SCRATCH}/wasp1_fastqs_${TIMESTAMP}"
mkdir -p "${LOCAL_FASTQ_DIR}"
LOCAL_FASTQ_R1="${LOCAL_FASTQ_DIR}/R1.fastq.gz"
LOCAL_FASTQ_R2="${LOCAL_FASTQ_DIR}/R2.fastq.gz"
echo "Copying FASTQs to local scratch: ${LOCAL_FASTQ_DIR}"
FASTQ_STAGE_START=$(date +%s.%N)
cp "${FASTQ_R1}" "${LOCAL_FASTQ_R1}"
cp "${FASTQ_R2}" "${LOCAL_FASTQ_R2}"
FASTQ_STAGE_END=$(date +%s.%N)
FASTQ_STAGE_TIME=$(echo "${FASTQ_STAGE_END} - ${FASTQ_STAGE_START}" | bc)
echo "FASTQ copy time: ${FASTQ_STAGE_TIME}s (NOT counted in benchmark)"

STAR_INDEX="${LOCAL_STAR_INDEX}"

# -----------------------------------------------------------------------------
# PRE-BENCHMARK: snp2h5 VCF -> HDF5 conversion (NOT timed in main benchmark)
# This is a one-time preprocessing step
# -----------------------------------------------------------------------------
echo "========================================"
echo "PRE-BENCHMARK: snp2h5 VCF -> HDF5 conversion"
echo "========================================"

SNP2H5_DIR="${OUTPUT_DIR}/snp2h5_data"
mkdir -p "${SNP2H5_DIR}"

# First split VCF by chromosome (required by snp2h5)
echo "Splitting VCF by chromosome..."
SNP2H5_PREP_START=$(date +%s.%N)

for i in {1..22}; do
    echo "  Processing chr${i}..."
    bcftools view -r "chr${i}" ${VCF} -Oz -o ${SNP2H5_DIR}/chr${i}.vcf.gz
done

echo "Running snp2h5..."
${SNP2H5} --chrom ${CHROM_SIZES} \
          --format "vcf" \
          --haplotype ${SNP2H5_DIR}/haplotypes.h5 \
          --snp_index ${SNP2H5_DIR}/snp_index.h5 \
          --snp_tab ${SNP2H5_DIR}/snp_tab.h5 \
          ${SNP2H5_DIR}/*.vcf.gz

SNP2H5_PREP_END=$(date +%s.%N)
SNP2H5_PREP_TIME=$(echo "${SNP2H5_PREP_END} - ${SNP2H5_PREP_START}" | bc)
echo "snp2h5 preprocessing time: ${SNP2H5_PREP_TIME}s (NOT counted in main benchmark)"
echo ""

# -----------------------------------------------------------------------------
# BENCHMARK STARTS HERE
# -----------------------------------------------------------------------------
echo "========================================"
echo "BENCHMARK TIMING STARTS NOW"
echo "========================================"

TOTAL_START=$(date +%s.%N)

# STAR parameters - use LoadAndKeep to cache genome for 2nd STAR pass
STARpar="--runThreadN ${THREADS} \
         --genomeDir ${STAR_INDEX} \
         --genomeLoad LoadAndKeep \
         --limitBAMsortRAM 20000000000 \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMattributes NH HI AS nM NM MD jM jI \
         --alignEndsType EndToEnd \
         --outSAMunmapped Within \
         --outFilterMultimapNmax 1"

readFiles="--readFilesCommand gunzip -c \
           --readFilesIn ${LOCAL_FASTQ_R1} ${LOCAL_FASTQ_R2}"

# -----------------------------------------------------------------------------
# STEP 1: Initial STAR alignment
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 1: Initial STAR alignment"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 1: STAR Initial Alignment" >> ${PROFILE_LOG}

STEP1_START=$(date +%s.%N)

${TIME} -v ${STAR} ${STARpar} ${readFiles} 2>> ${PROFILE_LOG}

mv Aligned.sortedByCoord.out.bam initial.bam
samtools index initial.bam

STEP1_END=$(date +%s.%N)
STEP1_TIME=$(echo "${STEP1_END} - ${STEP1_START}" | bc)
echo "STEP 1 completed in ${STEP1_TIME} seconds"
echo ""

# -----------------------------------------------------------------------------
# STEP 2: WASP1 find_intersecting_snps.py
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 2: WASP1 find_intersecting_snps.py"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 2: WASP1 find_intersecting_snps.py" >> ${PROFILE_LOG}

STEP2_START=$(date +%s.%N)

${TIME} -v python ${FIND_INTERSECTING} \
    --is_paired_end \
    --is_sorted \
    --output_dir ${OUTPUT_DIR} \
    --snp_index ${SNP2H5_DIR}/snp_index.h5 \
    --snp_tab ${SNP2H5_DIR}/snp_tab.h5 \
    --haplotype ${SNP2H5_DIR}/haplotypes.h5 \
    --samples ${SAMPLE} \
    initial.bam 2>> ${PROFILE_LOG}

STEP2_END=$(date +%s.%N)
STEP2_TIME=$(echo "${STEP2_END} - ${STEP2_START}" | bc)
echo "STEP 2 completed in ${STEP2_TIME} seconds"
echo ""

# Output files from find_intersecting_snps:
# - initial.keep.bam (reads that don't overlap SNPs - keep as-is)
# - initial.to.remap.bam (reads that need to be remapped)
# - initial.remap.fq1.gz / initial.remap.fq2.gz (FASTQs with flipped alleles)

# Report read counts
REMAP_R1=$(ls initial.remap.fq1.gz 2>/dev/null || echo "")
if [ -n "${REMAP_R1}" ]; then
    R1_READS=$(zcat initial.remap.fq1.gz | wc -l)
    R1_PAIRS=$((R1_READS / 4))
    echo "Reads to remap: ${R1_PAIRS} pairs"
fi

# -----------------------------------------------------------------------------
# STEP 3: STAR remap flipped reads
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 3: STAR remap flipped reads"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 3: STAR Remap" >> ${PROFILE_LOG}

STEP3_START=$(date +%s.%N)

# Decompress for STAR (STAR reads uncompressed faster)
echo "Decompressing FASTQs for STAR..."
zcat initial.remap.fq1.gz > remap_r1.fq
zcat initial.remap.fq2.gz > remap_r2.fq

${TIME} -v ${STAR} ${STARpar} \
    --readFilesIn remap_r1.fq remap_r2.fq 2>> ${PROFILE_LOG}

mv Aligned.sortedByCoord.out.bam remapped.bam
samtools index remapped.bam

# Clean up uncompressed FASTQs
rm -f remap_r1.fq remap_r2.fq

STEP3_END=$(date +%s.%N)
STEP3_TIME=$(echo "${STEP3_END} - ${STEP3_START}" | bc)
echo "STEP 3 completed in ${STEP3_TIME} seconds"
echo ""

# -----------------------------------------------------------------------------
# STEP 4: WASP1 filter_remapped_reads.py
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 4: WASP1 filter_remapped_reads.py"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 4: WASP1 filter_remapped_reads.py" >> ${PROFILE_LOG}

STEP4_START=$(date +%s.%N)

${TIME} -v python ${FILTER_REMAPPED} \
    initial.to.remap.bam \
    remapped.bam \
    remap_filtered.bam 2>> ${PROFILE_LOG}

STEP4_END=$(date +%s.%N)
STEP4_TIME=$(echo "${STEP4_END} - ${STEP4_START}" | bc)
echo "STEP 4 completed in ${STEP4_TIME} seconds"
echo ""

# -----------------------------------------------------------------------------
# STEP 5: Merge BAMs
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 5: Merge BAMs (keep + filtered)"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 5: Merge BAMs" >> ${PROFILE_LOG}

STEP5_START=$(date +%s.%N)

# Merge keep BAM with filtered remapped BAM
samtools merge -@ ${THREADS} -o wasp1_final.bam initial.keep.bam remap_filtered.bam
samtools sort -@ ${THREADS} -o wasp1_final_sorted.bam wasp1_final.bam
samtools index wasp1_final_sorted.bam
rm wasp1_final.bam

STEP5_END=$(date +%s.%N)
STEP5_TIME=$(echo "${STEP5_END} - ${STEP5_START}" | bc)
echo "STEP 5 completed in ${STEP5_TIME} seconds"
echo ""

# -----------------------------------------------------------------------------
# SUMMARY
# -----------------------------------------------------------------------------
TOTAL_END=$(date +%s.%N)
TOTAL_TIME=$(echo "${TOTAL_END} - ${TOTAL_START}" | bc)

# Convert to mm:ss format
TOTAL_MIN=$(echo "${TOTAL_TIME} / 60" | bc)
TOTAL_SEC=$(echo "${TOTAL_TIME} - (${TOTAL_MIN} * 60)" | bc | xargs printf "%.0f")

# Calculate WASP-only time (excluding STAR alignment)
WASP_ONLY=$(echo "${STEP2_TIME} + ${STEP4_TIME} + ${STEP5_TIME}" | bc)

echo "========================================"
echo "BENCHMARK RESULTS"
echo "========================================"
echo ""
echo "Individual step times:"
echo "  STEP 1 (STAR initial):             ${STEP1_TIME} s"
echo "  STEP 2 (WASP1 find_intersecting):  ${STEP2_TIME} s"
echo "  STEP 3 (STAR remap):               ${STEP3_TIME} s"
echo "  STEP 4 (WASP1 filter_remapped):    ${STEP4_TIME} s"
echo "  STEP 5 (Merge BAMs):               ${STEP5_TIME} s"
echo ""
echo "WASP-only time (steps 2+4+5):        ${WASP_ONLY} s"
echo "TOTAL time:                          ${TOTAL_TIME} s (${TOTAL_MIN}:${TOTAL_SEC})"
echo ""
echo "Pre-benchmark (NOT counted):"
echo "  snp2h5 preprocessing:              ${SNP2H5_PREP_TIME} s"
echo ""
echo "Published comparison times (HG00731, 56M reads, 8 threads):"
echo "  STAR+WASP (C++):    331.67 s ( 5:32)"
echo ""

# Calculate speedups
STAR_WASP_TIME=331.67
SPEEDUP_VS_STAR_WASP=$(echo "scale=2; ${TOTAL_TIME} / ${STAR_WASP_TIME}" | bc)

echo "WASP1 vs STAR+WASP: ${SPEEDUP_VS_STAR_WASP}x slower"
echo ""
echo "Output directory: ${OUTPUT_DIR}"
echo "End time: $(date)"

# Save results to JSON
cat > benchmark_results.json << EOF
{
    "timestamp": "${TIMESTAMP}",
    "pipeline": "wasp1_python",
    "sample": "${SAMPLE}",
    "threads": ${THREADS},
    "step1_star_initial_s": ${STEP1_TIME},
    "step2_wasp1_find_intersecting_s": ${STEP2_TIME},
    "step3_star_remap_s": ${STEP3_TIME},
    "step4_wasp1_filter_remapped_s": ${STEP4_TIME},
    "step5_merge_bams_s": ${STEP5_TIME},
    "wasp_only_s": ${WASP_ONLY},
    "total_s": ${TOTAL_TIME},
    "total_mmss": "${TOTAL_MIN}:${TOTAL_SEC}",
    "snp2h5_preprocessing_s": ${SNP2H5_PREP_TIME},
    "comparison": {
        "star_wasp_s": 331.67
    }
}
EOF

echo ""
echo "Results saved to: ${OUTPUT_DIR}/benchmark_results.json"

# Remove genome from shared memory before cleanup
echo "Removing genome from shared memory..."
${STAR} --genomeDir ${STAR_INDEX} --genomeLoad Remove 2>/dev/null || true

# Cleanup local scratch
rm -rf "${LOCAL_FASTQ_DIR}"
rm -rf "${LOCAL_STAR_INDEX}"
echo "Cleaned up local scratch files"

# Append summary to profile log
echo "" >> ${PROFILE_LOG}
echo "=======================================" >> ${PROFILE_LOG}
echo "SUMMARY" >> ${PROFILE_LOG}
echo "=======================================" >> ${PROFILE_LOG}
echo "STEP 1 (STAR initial):             ${STEP1_TIME} s" >> ${PROFILE_LOG}
echo "STEP 2 (WASP1 find_intersecting):  ${STEP2_TIME} s" >> ${PROFILE_LOG}
echo "STEP 3 (STAR remap):               ${STEP3_TIME} s" >> ${PROFILE_LOG}
echo "STEP 4 (WASP1 filter_remapped):    ${STEP4_TIME} s" >> ${PROFILE_LOG}
echo "STEP 5 (Merge BAMs):               ${STEP5_TIME} s" >> ${PROFILE_LOG}
echo "WASP-only (2+4+5):                 ${WASP_ONLY} s" >> ${PROFILE_LOG}
echo "TOTAL:                             ${TOTAL_TIME} s" >> ${PROFILE_LOG}
echo "" >> ${PROFILE_LOG}
echo "vs STAR+WASP: ${SPEEDUP_VS_STAR_WASP}x slower" >> ${PROFILE_LOG}
