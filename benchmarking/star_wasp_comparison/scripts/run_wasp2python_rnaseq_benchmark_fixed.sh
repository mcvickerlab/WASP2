#!/bin/bash
# WASP2-Python SNP Benchmark for STAR-WASP RNA-seq Comparison (FIXED)
# Uses upstream/dev branch for true multithreaded Python processing
# FIX: Pre-sorts VCF to match BAM chromosome ordering before processing
# Sample: HG00731 (1000 Genomes)
# Aligner: STAR (RNA-seq appropriate)
#
#$ -N wasp2py_star_rnaseq_fixed
#$ -V
#$ -pe iblm 8
#$ -l h_vmem=48G
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
FINAL_OUTPUT_DIR="${BENCHMARK_DIR}/results/wasp2python_rnaseq_fixed_${TIMESTAMP}"

# Input data
FASTQ_R1="${DATA_DIR}/ERR1050079_1.fastq.gz"
FASTQ_R2="${DATA_DIR}/ERR1050079_2.fastq.gz"
VCF_ORIGINAL="${DATA_DIR}/HG00731_het_only_chr.vcf.gz"

# Reference - STAR index
STAR_INDEX_REMOTE="/iblm/netapp/data1/external/GRC38/combined/google_cloud/star_index"
LOCAL_SCRATCH="${TMPDIR:-/tmp}"
WORK_OUTPUT_DIR="${LOCAL_SCRATCH}/wasp2python_rnaseq_fixed_${JOB_ID:-local}_${TIMESTAMP}"
OUTPUT_DIR="${WORK_OUTPUT_DIR}"

# WASP2-Python dev branch location
WASP2_DEV_DIR="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-python-dev"

SAMPLE="HG00731"
THREADS=8

# Conda
source /iblm/netapp/home/jjaureguy/mambaforge/etc/profile.d/conda.sh
conda activate WASP2_dev2

# Tools
STAR=$(which STAR)
SAMTOOLS=$(which samtools)
BCFTOOLS=$(which bcftools)
PYTHON=$(which python)
TIME="/usr/bin/time"

mkdir -p "${OUTPUT_DIR}"
cd "${OUTPUT_DIR}"

# Log file
PROFILE_LOG="${OUTPUT_DIR}/profile.log"
echo "WASP2-Python DEV BRANCH (Multithreaded) Benchmark Profile (FIXED)" > ${PROFILE_LOG}
echo "Sample: HG00731 RNA-seq" >> ${PROFILE_LOG}
echo "Timestamp: ${TIMESTAMP}" >> ${PROFILE_LOG}
echo "Threads: ${THREADS}" >> ${PROFILE_LOG}
echo "FIX: Pre-sorted VCF for chromosome ordering compatibility" >> ${PROFILE_LOG}

echo "========================================"
echo "WASP2-Python RNA-seq Benchmark (FIXED)"
echo "Using MULTITHREADING (--threads ${THREADS})"
echo "Sample: HG00731"
echo "Timestamp: ${TIMESTAMP}"
echo "FIX: VCF pre-sorted to match BAM chromosome order"
echo "========================================"
echo "Start time: $(date)"
echo "FASTQs: ${FASTQ_R1}"
echo "VCF: ${VCF_ORIGINAL}"
echo "Threads: ${THREADS}"

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
# PRE-BENCHMARK: Sort VCF to match STAR BAM lexicographic chromosome order
# STAR produces BAM with lexicographic order: chr1, chr10, chr11... chr2, chr20...
# VCF is in natural order: chr1, chr2, chr3... chr10, chr11...
# bedtools requires matching sort order for intersection
# -----------------------------------------------------------------------------
echo "========================================"
echo "PRE-BENCHMARK: VCF Setup (Chromosome Re-sorting)"
echo "========================================"

# Create sorted VCF matching STAR's lexicographic chromosome order
VCF="${OUTPUT_DIR}/HG00731_het_only_chr_lexsorted.vcf.gz"
echo "Sorting VCF to match STAR BAM chromosome order (lexicographic)..."

# Extract chromosome order from STAR BAM header (will be created in Step 0)
# For now, create a pre-defined lexicographic order file
cat > "${OUTPUT_DIR}/chrom_lex_order.txt" << 'CHROMEOF'
chr1
chr10
chr11
chr12
chr13
chr14
chr15
chr16
chr17
chr18
chr19
chr2
chr20
chr21
chr22
chr3
chr4
chr5
chr6
chr7
chr8
chr9
chrM
chrX
chrY
CHROMEOF

# Sort VCF to match BAM lexicographic order
${BCFTOOLS} view -h "${VCF_ORIGINAL}" > "${OUTPUT_DIR}/vcf_header.txt"
${BCFTOOLS} view -H "${VCF_ORIGINAL}" | \
    awk -F'\t' 'BEGIN{
        # Read chromosome order
        while((getline line < "'"${OUTPUT_DIR}/chrom_lex_order.txt"'") > 0) {
            order[line] = ++n
        }
    }
    { print order[$1], $0 }' | \
    sort -k1,1n -k3,3n | \
    cut -d' ' -f2- > "${OUTPUT_DIR}/vcf_sorted_body.txt"

cat "${OUTPUT_DIR}/vcf_header.txt" "${OUTPUT_DIR}/vcf_sorted_body.txt" | ${BCFTOOLS} view -Oz -o "${VCF}"
${BCFTOOLS} index -t "${VCF}"

# Cleanup temp files
rm -f "${OUTPUT_DIR}/vcf_header.txt" "${OUTPUT_DIR}/vcf_sorted_body.txt" "${OUTPUT_DIR}/chrom_lex_order.txt"

echo "Created lexicographically-sorted VCF: ${VCF}"

# Verify VCF is valid
echo "Verifying VCF file..."
${BCFTOOLS} view -h "${VCF}" | head -5
echo ""
echo "First 3 variants (should be chr1, then chr10 etc):"
${BCFTOOLS} view -H "${VCF}" | head -3 | cut -f1-5
echo ""
echo "Sample of chr10 (should come right after chr1):"
${BCFTOOLS} view -H "${VCF}" | grep "^chr10" | head -1 | cut -f1-5
echo "VCF verified with lexicographic chromosome order"

# -----------------------------------------------------------------------------
# Setup: Verify dev branch
# -----------------------------------------------------------------------------
echo "========================================"
echo "SETUP: Verify upstream/dev branch"
echo "========================================"

if [ ! -d "${WASP2_DEV_DIR}" ]; then
    echo "Cloning WASP2 dev branch..."
    git clone https://github.com/mcvickerlab/WASP2.git "${WASP2_DEV_DIR}"
    cd "${WASP2_DEV_DIR}"
    git checkout dev
else
    echo "Using existing WASP2 dev checkout at ${WASP2_DEV_DIR}"
fi

cd "${OUTPUT_DIR}"

# Set PYTHONPATH for dev branch
export PYTHONPATH="${WASP2_DEV_DIR}/src:${PYTHONPATH}"

# Verify we have the multithreaded version
echo "Verifying dev branch has --threads option..."
${PYTHON} "${WASP2_DEV_DIR}/src/mapping/__main__.py" make-reads --help | grep -q "threads" || { echo "ERROR: Dev branch missing --threads"; exit 1; }
echo "Dev branch verified with --threads support"

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

# STAR params MATCH STAR+WASP benchmark exactly (except waspOutputMode)
# Key: EndToEnd alignment, no mismatch filter, output unmapped reads
${TIME} -v ${STAR} \
    --runThreadN ${THREADS} \
    --genomeDir ${LOCAL_STAR_INDEX} \
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

# Create variants.bed for SNV counting
echo "Creating variants.bed..."
${BCFTOOLS} query -f '%CHROM\t%POS0\t%END\n' ${VCF} > ${OUTPUT_DIR}/variants.bed
SNP_COUNT=$(wc -l < ${OUTPUT_DIR}/variants.bed)
echo "Het SNP count: ${SNP_COUNT}"

echo "========================================"
echo "BENCHMARK TIMING STARTS NOW"
echo "========================================"

# -----------------------------------------------------------------------------
# STEP 1: WASP2-Python make-reads (MULTITHREADED)
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 1: WASP2-Python make-reads (${THREADS} threads)"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 1: WASP2-Python make-reads (multithreaded)" >> ${PROFILE_LOG}

STEP1_START=$(date +%s.%N)

${TIME} -v ${PYTHON} "${WASP2_DEV_DIR}/src/mapping/__main__.py" make-reads \
    "${ORIGINAL_BAM}" \
    "${VCF}" \
    --samples ${SAMPLE} \
    --out_dir "${OUTPUT_DIR}" \
    --threads ${THREADS} 2>> ${PROFILE_LOG}

STEP1_END=$(date +%s.%N)
STEP1_TIME=$(echo "${STEP1_END} - ${STEP1_START}" | bc)
echo "STEP 1 completed in ${STEP1_TIME} seconds"

# Find output files
WASP_JSON=$(find ${OUTPUT_DIR} -name "*wasp_data_files.json" | head -1)
KEEP_BAM=$(find ${OUTPUT_DIR} -name "*keep.bam" | head -1)
REMAP_FQ1=$(find ${OUTPUT_DIR} -name "*swapped_alleles_r1*" -o -name "*_r1.fq*" 2>/dev/null | head -1)
REMAP_FQ2=$(find ${OUTPUT_DIR} -name "*swapped_alleles_r2*" -o -name "*_r2.fq*" 2>/dev/null | head -1)

echo "  WASP JSON: ${WASP_JSON}"
echo "  KEEP BAM: ${KEEP_BAM}"
echo "  REMAP FQ1: ${REMAP_FQ1}"
echo "  REMAP FQ2: ${REMAP_FQ2}"

KEEP_READS=$(${SAMTOOLS} view -c ${KEEP_BAM})
echo "Keep reads: ${KEEP_READS}"

# Count SNV-overlapping read pairs (PRE-filter) - proper FASTQ counting
if [[ "${REMAP_FQ1}" == *.gz ]]; then
    REMAP_LINES=$(zcat ${REMAP_FQ1} | wc -l)
else
    REMAP_LINES=$(wc -l < ${REMAP_FQ1})
fi
REMAP_PAIRS=$((REMAP_LINES / 4))
echo "SNV_READS_PRE: ${REMAP_PAIRS}"
echo "Reads to remap (SNV-overlapping pairs): ${REMAP_PAIRS}"

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

# Determine read command based on compression
if [[ "${REMAP_FQ1}" == *.gz ]]; then
    READ_CMD="zcat"
else
    READ_CMD="cat"
fi

# STAR remap params MATCH STAR+WASP benchmark exactly
${TIME} -v ${STAR} \
    --runThreadN ${THREADS} \
    --genomeDir ${LOCAL_STAR_INDEX} \
    --readFilesIn ${REMAP_FQ1} ${REMAP_FQ2} \
    --readFilesCommand ${READ_CMD} \
    --outFileNamePrefix ${STAR_REMAP}/${SAMPLE}_remap_ \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI AS nM NM MD jM jI \
    --alignEndsType EndToEnd \
    --outSAMunmapped Within \
    --outFilterMultimapNmax 1 2>> ${PROFILE_LOG}

REMAPPED_BAM="${STAR_REMAP}/${SAMPLE}_remap_Aligned.sortedByCoord.out.bam"
${SAMTOOLS} index ${REMAPPED_BAM}

STEP2_END=$(date +%s.%N)
STEP2_TIME=$(echo "${STEP2_END} - ${STEP2_START}" | bc)
echo "STEP 2 completed in ${STEP2_TIME} seconds"

# -----------------------------------------------------------------------------
# STEP 3: WASP2-Python filter-remapped (MULTITHREADED)
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 3: WASP2-Python filter-remapped (${THREADS} threads)"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 3: WASP2-Python filter-remapped (multithreaded)" >> ${PROFILE_LOG}

STEP3_START=$(date +%s.%N)

${TIME} -v ${PYTHON} "${WASP2_DEV_DIR}/src/mapping/__main__.py" filter-remapped \
    "${REMAPPED_BAM}" \
    --wasp_data_json "${WASP_JSON}" \
    --out_bam "${OUTPUT_DIR}/wasp_filtered.bam" \
    --remap_keep_bam "${OUTPUT_DIR}/remap_keep.bam" \
    --threads ${THREADS} 2>> ${PROFILE_LOG}

STEP3_END=$(date +%s.%N)
STEP3_TIME=$(echo "${STEP3_END} - ${STEP3_START}" | bc)
echo "STEP 3 completed in ${STEP3_TIME} seconds"

# Index output BAMs
${SAMTOOLS} index ${OUTPUT_DIR}/wasp_filtered.bam 2>/dev/null || true
${SAMTOOLS} index ${OUTPUT_DIR}/remap_keep.bam 2>/dev/null || true

# -----------------------------------------------------------------------------
# STEP 4: Count SNV-overlapping reads (POST-filter)
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 4: Count SNV-overlapping reads"
echo "========================================"

FINAL_READS=$(${SAMTOOLS} view -c ${OUTPUT_DIR}/wasp_filtered.bam 2>/dev/null || echo "0")
REMAP_KEEP_READS=$(${SAMTOOLS} view -c ${OUTPUT_DIR}/remap_keep.bam 2>/dev/null || echo "0")

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
    OVERALL_PASS_RATE=$(echo "scale=2; ${FINAL_READS} * 100 / ${ORIGINAL_READS}" | bc)
else
    OVERALL_PASS_RATE="0"
fi

echo "========================================"
echo "BENCHMARK RESULTS (WASP2-Python RNA-seq FIXED)"
echo "========================================"
echo ""
echo "Individual step times:"
echo "  STEP 0 (STAR alignment):    ${STEP0_TIME} s"
echo "  STEP 1 (make-reads MT):     ${STEP1_TIME} s"
echo "  STEP 2 (STAR remap):        ${STEP2_TIME} s"
echo "  STEP 3 (filter-remapped MT): ${STEP3_TIME} s"
echo ""
echo "WASP-only time (steps 1+3):   ${WASP_ONLY} s"
echo "TOTAL time:                   ${TOTAL_TIME} s"
echo ""
echo "SNV-Overlapping Read Pair Counts:"
echo "  SNV_READS_PRE:  ${REMAP_PAIRS}"
echo "  SNV_READS_POST: ${SNV_READS_POST}"
echo "  SNV_PASS_RATE:  ${SNV_PASS_RATE}%"
echo ""
echo "Overall Read counts:"
echo "  Original reads:       ${ORIGINAL_READS}"
echo "  Keep reads:           ${KEEP_READS}"
echo "  Remap keep reads:     ${REMAP_KEEP_READS}"
echo "  FINAL (filtered):     ${FINAL_READS}"
echo "  Overall pass rate:    ${OVERALL_PASS_RATE}%"
echo ""
echo "Variant counts:"
echo "  Het SNPs:             ${SNP_COUNT}"
echo ""
echo "Work directory: ${OUTPUT_DIR}"
echo "Final output directory: ${FINAL_OUTPUT_DIR}"
echo "End time: $(date)"

# Save results to JSON
cat > ${OUTPUT_DIR}/benchmark_results.json << EOF
{
    "timestamp": "${TIMESTAMP}",
    "pipeline": "wasp2python_rnaseq_fixed",
    "branch": "upstream/dev",
    "sample": "${SAMPLE}",
    "data_type": "RNA-seq",
    "aligner": "STAR",
    "threads": ${THREADS},
    "het_only": true,
    "step0_star_align_s": ${STEP0_TIME},
    "step1_make_reads_s": ${STEP1_TIME},
    "step2_star_remap_s": ${STEP2_TIME},
    "step3_filter_remapped_s": ${STEP3_TIME},
    "wasp_only_s": ${WASP_ONLY},
    "total_s": ${TOTAL_TIME},
    "snv_reads_pre": ${REMAP_PAIRS},
    "snv_reads_post": ${SNV_READS_POST},
    "snv_pass_rate_percent": ${SNV_PASS_RATE},
    "original_reads": ${ORIGINAL_READS},
    "keep_reads": ${KEEP_READS},
    "remap_keep_reads": ${REMAP_KEEP_READS},
    "final_reads": ${FINAL_READS},
    "overall_pass_rate_percent": ${OVERALL_PASS_RATE},
    "snv_count": ${SNP_COUNT}
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
    "${PROFILE_LOG}"
do
    if [ -f "${f}" ]; then
        rsync -a "${f}" "${FINAL_OUTPUT_DIR}/"
    fi
done

# Optional: keep large outputs
if [ "${KEEP_BAMS:-0}" = "1" ]; then
    for f in "${OUTPUT_DIR}/remap_keep.bam" "${OUTPUT_DIR}/remap_keep.bam.bai" "${OUTPUT_DIR}/wasp_filtered.bam" "${OUTPUT_DIR}/wasp_filtered.bam.bai"; do
        if [ -f "${f}" ]; then
            rsync -a "${f}" "${FINAL_OUTPUT_DIR}/"
        fi
    done
fi

echo "Final output directory: ${FINAL_OUTPUT_DIR}"

# Cleanup STAR index from scratch
echo "Cleaning up STAR index from local scratch..."
rm -rf "${LOCAL_STAR_INDEX}" 2>/dev/null || true

cd /
rm -rf "${WORK_OUTPUT_DIR}" 2>/dev/null || true

echo "Done!"
