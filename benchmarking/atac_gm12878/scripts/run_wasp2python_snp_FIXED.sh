#!/bin/bash
# WASP2-Python FIXED benchmark on GM12878 ATAC-seq data - SNP only mode
# CORRECT pipeline: make-reads → BWA remap → filter-remapped
# FIXED: Uses het-only SNPs (samples parameter triggers het filtering)
#
#$ -N wasp2py_gm12878_snp_fixed
#$ -V
#$ -pe iblm 8
#$ -l h_vmem=48G
#$ -j y
#$ -o /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/atac_gm12878/logs/
#$ -cwd

set -e

# Timestamp for this run
TIMESTAMP=$(date +%Y-%m-%d_%H-%M-%S)
echo "Benchmark timestamp: ${TIMESTAMP}"

# Paths
WASP2_DIR="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp"
WASP2_PYTHON_DIR="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-python-bench"
BENCHMARK_DIR="${WASP2_DIR}/benchmarking/atac_gm12878"
OUTPUT_DIR="${BENCHMARK_DIR}/results/wasp2python_snp_fixed_${TIMESTAMP}"
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
BCFTOOLS=$(which bcftools)
PYTHON=$(which python)
TIME="/usr/bin/time"

# Set PYTHONPATH for WASP2-Python
export PYTHONPATH="${WASP2_PYTHON_DIR}/src:${PYTHONPATH}"

cd "${OUTPUT_DIR}"

# Log file
PROFILE_LOG="${OUTPUT_DIR}/profile.log"
echo "WASP2-Python FIXED Pipeline Benchmark - GM12878 ATAC-seq" > ${PROFILE_LOG}
echo "=========================================================" >> ${PROFILE_LOG}
echo "Timestamp: ${TIMESTAMP}" >> ${PROFILE_LOG}
echo "Date: $(date)" >> ${PROFILE_LOG}
echo "Sample: ${SAMPLE}" >> ${PROFILE_LOG}
echo "Threads: ${THREADS}" >> ${PROFILE_LOG}
echo "" >> ${PROFILE_LOG}
echo "WASP2-Python source: ${WASP2_PYTHON_DIR}" >> ${PROFILE_LOG}
echo "" >> ${PROFILE_LOG}

echo "========================================"
echo "WASP2-Python FIXED Benchmark on GM12878 ATAC-seq"
echo "Timestamp: ${TIMESTAMP}"
echo "========================================"
echo "Pipeline: make-reads → BWA remap → filter-remapped"
echo "FIXED: Uses het-only SNPs (samples parameter)"
echo "Start time: $(date)"
echo "Input BAM: ${INPUT_BAM}"
echo "VCF: ${VCF}"
echo "Threads: ${THREADS}"
echo ""

# Verify WASP2-Python module (must run from mapping dir due to relative imports)
echo "Verifying WASP2-Python module..."
cd "${WASP2_PYTHON_DIR}/src/mapping"
${PYTHON} __main__.py --help > /dev/null 2>&1 && echo "WASP2-Python mapping module: AVAILABLE" || {
    echo "ERROR: WASP2-Python mapping module not available"
    exit 1
}
cd "${OUTPUT_DIR}"
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

# -----------------------------------------------------------------------------
# PRE-BENCHMARK: Filter VCF to het-only SNPs
# -----------------------------------------------------------------------------
echo "========================================"
echo "PRE-BENCHMARK: Filter VCF to het-only SNPs"
echo "========================================"

HET_VCF="${OUTPUT_DIR}/het_only.vcf.gz"

echo "Filtering VCF to heterozygous SNPs only..."
${BCFTOOLS} view -i 'GT="het"' -v snps ${VCF} -Oz -o ${HET_VCF}
${BCFTOOLS} index -t ${HET_VCF}

ALL_SNPS=$(${BCFTOOLS} view -v snps ${VCF} 2>/dev/null | grep -v "^#" | wc -l)
HET_SNPS=$(${BCFTOOLS} view ${HET_VCF} 2>/dev/null | grep -v "^#" | wc -l)
echo "  Total SNPs in VCF: ${ALL_SNPS}"
echo "  Heterozygous SNPs: ${HET_SNPS}"

# Create BED for post-analysis
${BCFTOOLS} query -f '%CHROM\t%POS0\t%POS\t%REF\t%ALT\n' ${HET_VCF} > ${OUTPUT_DIR}/variants.bed
VARIANT_COUNT=$(wc -l < ${OUTPUT_DIR}/variants.bed)
echo "Het-only SNP variants in BED: ${VARIANT_COUNT}"

# -----------------------------------------------------------------------------
# BENCHMARK STARTS HERE
# -----------------------------------------------------------------------------
echo "========================================"
echo "BENCHMARK TIMING STARTS NOW"
echo "========================================"

TOTAL_START=$(date +%s.%N)

# -----------------------------------------------------------------------------
# STEP 1: WASP2-Python make-reads
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 1: WASP2-Python make-reads"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 1: WASP2-Python make-reads" >> ${PROFILE_LOG}

STEP1_START=$(date +%s.%N)

cd "${WASP2_PYTHON_DIR}/src/mapping"

${TIME} -v ${PYTHON} __main__.py make-reads \
    "${ORIGINAL_BAM}" \
    "${HET_VCF}" \
    --samples "${SAMPLE}" \
    --out_dir "${OUTPUT_DIR}" \
    --paired 2>> ${PROFILE_LOG}

cd "${OUTPUT_DIR}"

STEP1_END=$(date +%s.%N)
STEP1_TIME=$(echo "${STEP1_END} - ${STEP1_START}" | bc)
echo "STEP 1 completed in ${STEP1_TIME} seconds"

# Find output files from make-reads
# WASP2-Python outputs: to_remap.bam, keep.bam, remap FASTQs
TO_REMAP_BAM=$(find ${OUTPUT_DIR} -name "*to_remap*.bam" | head -1)
KEEP_BAM=$(find ${OUTPUT_DIR} -name "*keep*.bam" | head -1)
# WASP2-Python outputs files like: original_swapped_alleles_r1.fq
REMAP_FQ1=$(find ${OUTPUT_DIR} -name "*swapped_alleles_r1*" -o -name "*remap*r1*" -o -name "*_r1.fq*" 2>/dev/null | head -1)
REMAP_FQ2=$(find ${OUTPUT_DIR} -name "*swapped_alleles_r2*" -o -name "*remap*r2*" -o -name "*_r2.fq*" 2>/dev/null | head -1)

echo "Found files:"
echo "  to_remap BAM: ${TO_REMAP_BAM}"
echo "  keep BAM: ${KEEP_BAM}"
echo "  remap FQ1: ${REMAP_FQ1}"
echo "  remap FQ2: ${REMAP_FQ2}"

# Index BAMs
${SAMTOOLS} index ${TO_REMAP_BAM} 2>/dev/null || true
${SAMTOOLS} index ${KEEP_BAM} 2>/dev/null || true

# Count keep reads
KEEP_READS=$(${SAMTOOLS} view -c ${KEEP_BAM} 2>/dev/null || echo "0")
echo "Keep reads (no variants): ${KEEP_READS}"

# Count reads to remap
if [ -n "${REMAP_FQ1}" ] && [ -f "${REMAP_FQ1}" ]; then
    if [[ "${REMAP_FQ1}" == *.gz ]]; then
        REMAP_LINES=$(zcat ${REMAP_FQ1} | wc -l)
    else
        REMAP_LINES=$(wc -l < ${REMAP_FQ1})
    fi
    REMAP_PAIRS=$((REMAP_LINES / 4))
    echo "Reads to remap: ${REMAP_PAIRS} pairs"
fi

# -----------------------------------------------------------------------------
# STEP 2: BWA-MEM remap flipped reads
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 2: BWA-MEM remap flipped reads"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 2: BWA-MEM Remap" >> ${PROFILE_LOG}

STEP2_START=$(date +%s.%N)

# Decompress if needed
if [[ "${REMAP_FQ1}" == *.gz ]]; then
    zcat ${REMAP_FQ1} > ${OUTPUT_DIR}/remap_r1.fq
    zcat ${REMAP_FQ2} > ${OUTPUT_DIR}/remap_r2.fq
    REMAP_R1="${OUTPUT_DIR}/remap_r1.fq"
    REMAP_R2="${OUTPUT_DIR}/remap_r2.fq"
else
    REMAP_R1="${REMAP_FQ1}"
    REMAP_R2="${REMAP_FQ2}"
fi

${TIME} -v ${BWA} mem -t ${THREADS} -M ${REF_GENOME} \
    ${REMAP_R1} ${REMAP_R2} 2>> ${PROFILE_LOG} | \
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

# Find the wasp_data_files.json if it exists
WASP_JSON=$(find ${OUTPUT_DIR} -name "*wasp_data_files*.json" | head -1)

if [ -n "${WASP_JSON}" ] && [ -f "${WASP_JSON}" ]; then
    ${TIME} -v ${PYTHON} __main__.py filter-remapped \
        "${OUTPUT_DIR}/remapped.bam" \
        --wasp_data_json "${WASP_JSON}" \
        --out_bam "${OUTPUT_DIR}/wasp_filtered.bam" \
        --remap_keep_bam "${OUTPUT_DIR}/remap_keep.bam" 2>> ${PROFILE_LOG}
else
    ${TIME} -v ${PYTHON} __main__.py filter-remapped \
        "${OUTPUT_DIR}/remapped.bam" \
        "${TO_REMAP_BAM}" \
        "${KEEP_BAM}" \
        --out_bam "${OUTPUT_DIR}/wasp_filtered.bam" \
        --remap_keep_bam "${OUTPUT_DIR}/remap_keep.bam" 2>> ${PROFILE_LOG}
fi

cd "${OUTPUT_DIR}"

${SAMTOOLS} index ${OUTPUT_DIR}/wasp_filtered.bam 2>/dev/null || true
${SAMTOOLS} index ${OUTPUT_DIR}/remap_keep.bam 2>/dev/null || true

STEP3_END=$(date +%s.%N)
STEP3_TIME=$(echo "${STEP3_END} - ${STEP3_START}" | bc)
echo "STEP 3 completed in ${STEP3_TIME} seconds"

# -----------------------------------------------------------------------------
# SUMMARY
# -----------------------------------------------------------------------------
TOTAL_END=$(date +%s.%N)
TOTAL_TIME=$(echo "${TOTAL_END} - ${TOTAL_START}" | bc)

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
