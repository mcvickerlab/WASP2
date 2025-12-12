#!/bin/bash
# WASP2-Python SNP Benchmark - DEV BRANCH with MULTITHREADING
# Uses upstream/dev branch for true multithreaded Python processing
#
#$ -N wasp2py_dev_mt
#$ -V
#$ -pe iblm 8
#$ -l h_vmem=48G
#$ -j y
#$ -o /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/atac_gm12878/logs/
#$ -cwd

set -e

# Timestamp
TIMESTAMP=$(date +%Y-%m-%d_%H-%M-%S)
echo "Benchmark timestamp: ${TIMESTAMP}"

# Paths
INPUT_BAM="/iblm/netapp/data1/aho/atac/Buenrostro2013/merged/GM12878_ATACseq_50k/GM12878_ATACseq_50k_merged.sorted.bam"
VCF="/iblm/netapp/data1/aho/variants/NA12878.vcf.gz"
REF_GENOME="/iblm/netapp/data1/aho/ref_genomes/index/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/genome.fa"
OUTPUT_BASE="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/atac_gm12878/results"
OUTPUT_DIR="${OUTPUT_BASE}/wasp2python_snp_DEV_MT_${TIMESTAMP}"

# WASP2-Python dev branch location - we'll create a fresh checkout
WASP2_DEV_DIR="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-python-dev"

THREADS=8
SAMPLE="NA12878"

# Conda
source /iblm/netapp/home/jjaureguy/mambaforge/etc/profile.d/conda.sh
conda activate WASP2_dev2

# Tools
BWA=$(which bwa)
SAMTOOLS=$(which samtools)
BCFTOOLS=$(which bcftools)
PYTHON=$(which python)
TIME="/usr/bin/time"

mkdir -p "${OUTPUT_DIR}"
cd "${OUTPUT_DIR}"

# Profile log
PROFILE_LOG="${OUTPUT_DIR}/profile.log"
echo "WASP2-Python DEV BRANCH (Multithreaded) Benchmark Profile" > ${PROFILE_LOG}
echo "Timestamp: ${TIMESTAMP}" >> ${PROFILE_LOG}
echo "Threads: ${THREADS}" >> ${PROFILE_LOG}

echo "========================================"
echo "WASP2-Python DEV BRANCH SNP Benchmark"
echo "Using MULTITHREADING (--threads ${THREADS})"
echo "Timestamp: ${TIMESTAMP}"
echo "========================================"
echo "Start time: $(date)"
echo "Input BAM: ${INPUT_BAM}"
echo "VCF: ${VCF}"
echo "Threads: ${THREADS}"

# -----------------------------------------------------------------------------
# Setup: Clone/checkout dev branch
# -----------------------------------------------------------------------------
echo "========================================"
echo "SETUP: Checkout upstream/dev branch"
echo "========================================"

if [ ! -d "${WASP2_DEV_DIR}" ]; then
    echo "Cloning WASP2 dev branch..."
    git clone https://github.com/mcvickerlab/WASP2.git "${WASP2_DEV_DIR}"
    cd "${WASP2_DEV_DIR}"
    git checkout dev
else
    echo "Using existing WASP2 dev checkout at ${WASP2_DEV_DIR}"
    cd "${WASP2_DEV_DIR}"
    git fetch origin
    git checkout dev
    git pull origin dev
fi

cd "${OUTPUT_DIR}"

# Set PYTHONPATH for dev branch
export PYTHONPATH="${WASP2_DEV_DIR}/src:${PYTHONPATH}"

# Verify we have the multithreaded version
echo "Verifying dev branch has --threads option..."
${PYTHON} "${WASP2_DEV_DIR}/src/mapping/__main__.py" make-reads --help | grep -q "threads" || { echo "ERROR: Dev branch missing --threads"; exit 1; }
echo "Dev branch verified with --threads support"

# -----------------------------------------------------------------------------
# PRE-BENCHMARK: Filter VCF to het-only SNPs
# -----------------------------------------------------------------------------
echo "========================================"
echo "PRE-BENCHMARK: Filter VCF to het-only SNPs"
echo "========================================"

HET_VCF="${OUTPUT_DIR}/het_only_snps.vcf.gz"
echo "Filtering VCF to heterozygous SNPs only..."
${BCFTOOLS} view -g het -v snps ${VCF} -Oz -o ${HET_VCF}
${BCFTOOLS} index -t ${HET_VCF}

SNP_COUNT=$(${BCFTOOLS} view -H ${HET_VCF} | wc -l)
echo "Het-only SNP count: ${SNP_COUNT}"

# Create variants.bed for comparison
${BCFTOOLS} query -f '%CHROM\t%POS0\t%END\n' ${HET_VCF} > ${OUTPUT_DIR}/variants.bed

# Copy input BAM
echo "Copying input BAM to output directory..."
cp ${INPUT_BAM} ${OUTPUT_DIR}/original.bam
${SAMTOOLS} index ${OUTPUT_DIR}/original.bam
ORIGINAL_READS=$(${SAMTOOLS} view -c ${OUTPUT_DIR}/original.bam)
echo "Original reads: ${ORIGINAL_READS}"

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
    "${OUTPUT_DIR}/original.bam" \
    "${HET_VCF}" \
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

REMAP_LINES=$(wc -l < ${REMAP_FQ1})
REMAP_PAIRS=$((REMAP_LINES / 4))
echo "Reads to remap: ${REMAP_PAIRS} pairs"

# -----------------------------------------------------------------------------
# STEP 2: BWA-MEM remap flipped reads (not counted in WASP-only)
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 2: BWA-MEM remap flipped reads"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 2: BWA-MEM Remap" >> ${PROFILE_LOG}

STEP2_START=$(date +%s.%N)

${TIME} -v ${BWA} mem -t ${THREADS} -M ${REF_GENOME} \
    ${REMAP_FQ1} ${REMAP_FQ2} 2>> ${PROFILE_LOG} | \
    ${SAMTOOLS} view -S -b -h -F 4 - > ${OUTPUT_DIR}/remapped_unsorted.bam

${SAMTOOLS} sort -@ ${THREADS} -o ${OUTPUT_DIR}/remapped.bam ${OUTPUT_DIR}/remapped_unsorted.bam
${SAMTOOLS} index ${OUTPUT_DIR}/remapped.bam
rm -f ${OUTPUT_DIR}/remapped_unsorted.bam

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
    "${OUTPUT_DIR}/remapped.bam" \
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
# SUMMARY
# -----------------------------------------------------------------------------
TOTAL_TIME=$(echo "${STEP1_TIME} + ${STEP2_TIME} + ${STEP3_TIME}" | bc)
WASP_ONLY=$(echo "${STEP1_TIME} + ${STEP3_TIME}" | bc)

FINAL_READS=$(${SAMTOOLS} view -c ${OUTPUT_DIR}/wasp_filtered.bam 2>/dev/null || echo "0")
REMAP_KEEP_READS=$(${SAMTOOLS} view -c ${OUTPUT_DIR}/remap_keep.bam 2>/dev/null || echo "0")

if [ "${ORIGINAL_READS}" -gt 0 ]; then
    PASS_RATE=$(echo "scale=2; ${FINAL_READS} * 100 / ${ORIGINAL_READS}" | bc)
else
    PASS_RATE="0"
fi

echo "========================================"
echo "BENCHMARK RESULTS (WASP2-Python DEV MT)"
echo "========================================"
echo ""
echo "Individual step times:"
echo "  STEP 1 (make-reads MT):     ${STEP1_TIME} s"
echo "  STEP 2 (BWA-MEM remap):     ${STEP2_TIME} s"
echo "  STEP 3 (filter-remapped MT): ${STEP3_TIME} s"
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
echo "  Het SNPs:             ${SNP_COUNT}"
echo ""
echo "Output directory: ${OUTPUT_DIR}"
echo "End time: $(date)"

# Save results to JSON
cat > ${OUTPUT_DIR}/benchmark_results.json << EOF
{
    "timestamp": "${TIMESTAMP}",
    "pipeline": "wasp2python_snp_DEV_MULTITHREADED",
    "branch": "upstream/dev",
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
    "snv_count": ${SNP_COUNT}
}
EOF

echo ""
echo "Results saved to: ${OUTPUT_DIR}/benchmark_results.json"
