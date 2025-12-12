#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -N wasp_compare
#$ -l h_vmem=64G
#$ -l h_rt=4:00:00

# Comparison benchmark: Run STAR+WASP and WASP2-Rust back-to-back
# This ensures both use the same cluster conditions for fair comparison

set -e

TIMESTAMP=$(date +%Y-%m-%d_%H-%M-%S)
WASP2_DIR="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp"
SCRIPT_DIR="${WASP2_DIR}/benchmarking/star_wasp_comparison/scripts"
LOG_DIR="${WASP2_DIR}/benchmarking/logs"
RESULTS_DIR="${WASP2_DIR}/benchmarking/star_wasp_comparison/results/comparison_${TIMESTAMP}"

mkdir -p "${RESULTS_DIR}"
mkdir -p "${LOG_DIR}"

echo "========================================"
echo "WASP2 vs STAR+WASP Comparison Benchmark"
echo "Timestamp: ${TIMESTAMP}"
echo "========================================"
echo ""

# Activate conda environment
source /iblm/netapp/home/jjaureguy/mambaforge/etc/profile.d/conda.sh
conda activate WASP2_dev2

# Common paths
STAR_INDEX="/iblm/netapp/data1/external/GRC38/combined/google_cloud/star_index"
FASTQ_DIR="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/wasp2_star_wasp_evaluation/HG00731"
VCF="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/wasp2_star_wasp_evaluation/HG00731_GRCh38_phased.vcf.gz"
THREADS=8

# Stage files to local scratch once
LOCAL_SCRATCH="/tmp/${JOB_ID}.1.iblm.q"
LOCAL_STAR_INDEX="${LOCAL_SCRATCH}/star_index"
LOCAL_FASTQS="${LOCAL_SCRATCH}/fastqs"

echo "Staging files to local scratch..."
mkdir -p "${LOCAL_STAR_INDEX}"
mkdir -p "${LOCAL_FASTQS}"

STAGE_START=$(date +%s.%N)
cp -r ${STAR_INDEX}/* ${LOCAL_STAR_INDEX}/
cp ${FASTQ_DIR}/R1.fastq.gz ${LOCAL_FASTQS}/
cp ${FASTQ_DIR}/R2.fastq.gz ${LOCAL_FASTQS}/
STAGE_END=$(date +%s.%N)
STAGE_TIME=$(echo "${STAGE_END} - ${STAGE_START}" | bc)
echo "Staging time: ${STAGE_TIME}s (not counted)"
echo ""

# Results file
RESULTS_FILE="${RESULTS_DIR}/comparison_results.txt"
echo "Comparison Benchmark Results - ${TIMESTAMP}" > ${RESULTS_FILE}
echo "==========================================" >> ${RESULTS_FILE}
echo "" >> ${RESULTS_FILE}

# Run STAR+WASP benchmark
echo "========================================"
echo "RUN 1: STAR+WASP (C++ integrated)"
echo "========================================"

STAR_WASP_START=$(date +%s.%N)

# STAR with waspOutputMode
STAR --runThreadN ${THREADS} \
     --genomeDir ${LOCAL_STAR_INDEX} \
     --genomeLoad LoadAndKeep \
     --limitBAMsortRAM 20000000000 \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMattributes NH HI AS nM NM MD jM jI vW \
     --alignEndsType EndToEnd \
     --outSAMunmapped Within \
     --outFilterMultimapNmax 1 \
     --waspOutputMode SAMtag \
     --varVCFfile ${VCF} \
     --readFilesCommand gunzip -c \
     --readFilesIn ${LOCAL_FASTQS}/R1.fastq.gz ${LOCAL_FASTQS}/R2.fastq.gz \
     --outFileNamePrefix ${RESULTS_DIR}/star_wasp_

# Filter WASP-passing reads
samtools view -h -b -q 255 ${RESULTS_DIR}/star_wasp_Aligned.sortedByCoord.out.bam | \
    samtools view -h -b -e '[vW]==1' - > ${RESULTS_DIR}/star_wasp_filtered.bam
samtools index ${RESULTS_DIR}/star_wasp_filtered.bam

STAR_WASP_END=$(date +%s.%N)
STAR_WASP_TIME=$(echo "${STAR_WASP_END} - ${STAR_WASP_START}" | bc)

echo "STAR+WASP completed in ${STAR_WASP_TIME}s"
echo "STAR+WASP: ${STAR_WASP_TIME}s" >> ${RESULTS_FILE}

# Unload genome
STAR --genomeDir ${LOCAL_STAR_INDEX} --genomeLoad Remove 2>/dev/null || true

echo ""
echo "========================================"
echo "RUN 2: WASP2-Rust OPTIMIZED"
echo "========================================"

cd ${RESULTS_DIR}

WASP2_START=$(date +%s.%N)

# STEP 1: Initial STAR alignment
STEP1_START=$(date +%s.%N)
STAR --runThreadN ${THREADS} \
     --genomeDir ${LOCAL_STAR_INDEX} \
     --genomeLoad LoadAndKeep \
     --limitBAMsortRAM 20000000000 \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMattributes NH HI AS nM NM MD jM jI \
     --alignEndsType EndToEnd \
     --outSAMunmapped Within \
     --outFilterMultimapNmax 1 \
     --readFilesCommand gunzip -c \
     --readFilesIn ${LOCAL_FASTQS}/R1.fastq.gz ${LOCAL_FASTQS}/R2.fastq.gz

mv Aligned.sortedByCoord.out.bam A_sorted.bam
samtools index A_sorted.bam
STEP1_END=$(date +%s.%N)
STEP1_TIME=$(echo "${STEP1_END} - ${STEP1_START}" | bc)
echo "STEP 1 (STAR initial): ${STEP1_TIME}s"

# STEP 2: WASP2 unified pipeline (uncompressed output)
STEP2_START=$(date +%s.%N)
PYTHONPATH=${WASP2_DIR}/src python -c "
from wasp2_rust import unified_make_reads_parallel_py
from mapping.intersect_variant_data import vcf_to_bed

bed_file = 'variants.bed'
vcf_to_bed('${VCF}', bed_file, samples=['HG00731'], include_indels=False)

stats = unified_make_reads_parallel_py(
    'A_sorted.bam',
    bed_file,
    'remap_r1.fq',
    'remap_r2.fq',
    max_seqs=64,
    threads=${THREADS},
    compress_output=False
)
print(f'Haplotypes written: {stats[\"haplotypes_written\"]:,}')
"
STEP2_END=$(date +%s.%N)
STEP2_TIME=$(echo "${STEP2_END} - ${STEP2_START}" | bc)
echo "STEP 2 (WASP2 unified): ${STEP2_TIME}s"

# STEP 3: STAR remap (unsorted BAM)
STEP3_START=$(date +%s.%N)
STAR --runThreadN ${THREADS} \
     --genomeDir ${LOCAL_STAR_INDEX} \
     --genomeLoad LoadAndKeep \
     --outSAMtype BAM Unsorted \
     --outSAMattributes NH HI AS nM NM MD jM jI \
     --alignEndsType EndToEnd \
     --outSAMunmapped Within \
     --outFilterMultimapNmax 1 \
     --readFilesIn remap_r1.fq remap_r2.fq

mv Aligned.out.bam remapped.bam
rm -f remap_r1.fq remap_r2.fq
STEP3_END=$(date +%s.%N)
STEP3_TIME=$(echo "${STEP3_END} - ${STEP3_START}" | bc)
echo "STEP 3 (STAR remap): ${STEP3_TIME}s"

# STEP 4: Filter remapped reads
STEP4_START=$(date +%s.%N)
PYTHONPATH=${WASP2_DIR}/src python -c "
from wasp2_rust import filter_bam_wasp
kept, removed_moved, removed_missing = filter_bam_wasp(
    'A_sorted.bam',
    'remapped.bam',
    'remap_keep.bam',
    threads=${THREADS}
)
print(f'Kept: {kept}, Removed: {removed_moved + removed_missing}')
"
samtools index remap_keep.bam
STEP4_END=$(date +%s.%N)
STEP4_TIME=$(echo "${STEP4_END} - ${STEP4_START}" | bc)
echo "STEP 4 (WASP2 filter): ${STEP4_TIME}s"

WASP2_END=$(date +%s.%N)
WASP2_TIME=$(echo "${WASP2_END} - ${WASP2_START}" | bc)

echo ""
echo "WASP2-Rust completed in ${WASP2_TIME}s"
echo "WASP2-Rust: ${WASP2_TIME}s" >> ${RESULTS_FILE}
echo "  STEP 1: ${STEP1_TIME}s" >> ${RESULTS_FILE}
echo "  STEP 2: ${STEP2_TIME}s" >> ${RESULTS_FILE}
echo "  STEP 3: ${STEP3_TIME}s" >> ${RESULTS_FILE}
echo "  STEP 4: ${STEP4_TIME}s" >> ${RESULTS_FILE}

# Cleanup
STAR --genomeDir ${LOCAL_STAR_INDEX} --genomeLoad Remove 2>/dev/null || true
rm -rf ${LOCAL_SCRATCH}

# Summary
echo ""
echo "========================================"
echo "COMPARISON SUMMARY"
echo "========================================"
echo ""
DIFF=$(echo "${WASP2_TIME} - ${STAR_WASP_TIME}" | bc)
echo "STAR+WASP:   ${STAR_WASP_TIME}s"
echo "WASP2-Rust:  ${WASP2_TIME}s"
echo "Difference:  ${DIFF}s (positive = WASP2 slower)"
echo ""
echo "Results saved to: ${RESULTS_FILE}"

echo "" >> ${RESULTS_FILE}
echo "Difference: ${DIFF}s" >> ${RESULTS_FILE}
