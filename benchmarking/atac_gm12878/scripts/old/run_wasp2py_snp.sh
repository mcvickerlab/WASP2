#!/bin/bash
# WASP2-Python benchmark on GM12878 ATAC-seq data - SNP only mode
# Uses Python mapping pipeline (with Rust acceleration for filter step)
#
# For Figure 1B (timing) and Figure 1C (pre/post remap read counts)
#
#$ -N wasp2py_gm12878_snp
#$ -V
#$ -pe iblm 8
#$ -l h_vmem=32G
#$ -j y
#$ -o /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/atac_gm12878/logs/
#$ -cwd

set -e

# Timestamp for this run
TIMESTAMP=$(date +%Y-%m-%d_%H-%M-%S)
echo "Benchmark timestamp: ${TIMESTAMP}"

# Paths
WASP2_DIR="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp"
BENCHMARK_DIR="${WASP2_DIR}/benchmarking/atac_gm12878"
OUTPUT_DIR="${BENCHMARK_DIR}/results/wasp2py_snp_${TIMESTAMP}"
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
PYTHON=$(which python)
TIME="/usr/bin/time"

# Set PYTHONPATH to include src directory for mapping module
export PYTHONPATH="${WASP2_DIR}/src:${PYTHONPATH}"

cd "${WASP2_DIR}"

# Log file
PROFILE_LOG="${OUTPUT_DIR}/profile.log"
echo "WASP2-Python SNP-only Pipeline Benchmark - GM12878 ATAC-seq" > ${PROFILE_LOG}
echo "============================================================" >> ${PROFILE_LOG}
echo "Timestamp: ${TIMESTAMP}" >> ${PROFILE_LOG}
echo "Date: $(date)" >> ${PROFILE_LOG}
echo "Sample: ${SAMPLE}" >> ${PROFILE_LOG}
echo "Threads: ${THREADS}" >> ${PROFILE_LOG}
echo "" >> ${PROFILE_LOG}

echo "========================================"
echo "WASP2-Python SNP Benchmark on GM12878 ATAC-seq"
echo "Timestamp: ${TIMESTAMP}"
echo "========================================"
echo "Start time: $(date)"
echo "Input BAM: ${INPUT_BAM}"
echo "VCF: ${VCF}"
echo "Threads: ${THREADS}"
echo ""

# -----------------------------------------------------------------------------
# PRE-BENCHMARK: Copy input BAM to working directory (for output consistency)
# -----------------------------------------------------------------------------
echo "Copying input BAM to output directory..."
cp "${INPUT_BAM}" "${OUTPUT_DIR}/original.bam"
${SAMTOOLS} index "${OUTPUT_DIR}/original.bam"
ORIGINAL_BAM="${OUTPUT_DIR}/original.bam"

# Count original reads
ORIGINAL_READS=$(${SAMTOOLS} view -c ${ORIGINAL_BAM})
echo "Original reads: ${ORIGINAL_READS}"

# Extract variant BED for Figure 1C analysis
echo "Creating variant BED file for read counting..."
bcftools query -f '%CHROM\t%POS0\t%POS\t%REF\t%ALT\n' ${VCF} | \
    awk '$4 !~ /,/ && $5 !~ /,/ && length($4)==1 && length($5)==1' > ${OUTPUT_DIR}/variants.bed
VARIANT_COUNT=$(wc -l < ${OUTPUT_DIR}/variants.bed)
echo "SNP variants in BED: ${VARIANT_COUNT}"

# -----------------------------------------------------------------------------
# BENCHMARK STARTS HERE
# -----------------------------------------------------------------------------
echo "========================================"
echo "BENCHMARK TIMING STARTS NOW"
echo "========================================"

TOTAL_START=$(date +%s.%N)

# -----------------------------------------------------------------------------
# STEP 1: WASP2-Python make-reads (find intersecting variants)
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 1: WASP2-Python make-reads"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 1: WASP2-Python make-reads" >> ${PROFILE_LOG}

STEP1_START=$(date +%s.%N)

${TIME} -v ${PYTHON} -m mapping make-reads \
    ${ORIGINAL_BAM} \
    ${VCF} \
    -s ${SAMPLE} \
    -o ${OUTPUT_DIR} \
    --snps-only \
    --paired \
    --threads ${THREADS} \
    -j ${OUTPUT_DIR}/wasp_data.json 2>> ${PROFILE_LOG}

STEP1_END=$(date +%s.%N)
STEP1_TIME=$(echo "${STEP1_END} - ${STEP1_START}" | bc)
echo "STEP 1 completed in ${STEP1_TIME} seconds"

# Find the output files from make-reads using wasp_data.json
TO_REMAP_BAM=$(${PYTHON} -c "import json; d=json.load(open('${OUTPUT_DIR}/wasp_data.json')); print(d['to_remap_bam'])")
KEEP_BAM_ORIG=$(${PYTHON} -c "import json; d=json.load(open('${OUTPUT_DIR}/wasp_data.json')); print(d['keep_bam'])")
REMAP_FQ1=$(${PYTHON} -c "import json; d=json.load(open('${OUTPUT_DIR}/wasp_data.json')); print(d['remap_fq1'])")
REMAP_FQ2=$(${PYTHON} -c "import json; d=json.load(open('${OUTPUT_DIR}/wasp_data.json')); print(d['remap_fq2'])")

echo "TO_REMAP_BAM: ${TO_REMAP_BAM}"
echo "KEEP_BAM: ${KEEP_BAM_ORIG}"
echo "REMAP_FQ1: ${REMAP_FQ1}"
echo "REMAP_FQ2: ${REMAP_FQ2}"

if [ -n "${REMAP_FQ1}" ] && [ -f "${REMAP_FQ1}" ]; then
    R1_READS=$(wc -l < ${REMAP_FQ1})
    R1_PAIRS=$((R1_READS / 4))
    echo "Reads to remap: ${R1_PAIRS} pairs"
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

# Copy FASTQs (they are uncompressed .fq files from WASP2-Python)
cp ${REMAP_FQ1} ${OUTPUT_DIR}/remap_r1.fq
cp ${REMAP_FQ2} ${OUTPUT_DIR}/remap_r2.fq

${TIME} -v ${BWA} mem -t ${THREADS} ${REF_GENOME} \
    ${OUTPUT_DIR}/remap_r1.fq ${OUTPUT_DIR}/remap_r2.fq 2>> ${PROFILE_LOG} | \
    ${SAMTOOLS} view -bS - | \
    ${SAMTOOLS} sort -@ ${THREADS} -o ${OUTPUT_DIR}/remapped.bam -

${SAMTOOLS} index ${OUTPUT_DIR}/remapped.bam
rm -f ${OUTPUT_DIR}/remap_r1.fq ${OUTPUT_DIR}/remap_r2.fq

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

# Note: filter-remapped uses Rust acceleration (no pure Python fallback)
${TIME} -v ${PYTHON} -m mapping filter-remapped \
    ${OUTPUT_DIR}/remapped.bam \
    ${TO_REMAP_BAM} \
    ${KEEP_BAM_ORIG} \
    --threads ${THREADS} \
    -o ${OUTPUT_DIR}/remap_filtered.bam 2>> ${PROFILE_LOG}

STEP3_END=$(date +%s.%N)
STEP3_TIME=$(echo "${STEP3_END} - ${STEP3_START}" | bc)
echo "STEP 3 completed in ${STEP3_TIME} seconds"

# -----------------------------------------------------------------------------
# STEP 4: Merge BAMs to create final keep BAM
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 4: Merge BAMs (keep + filtered)"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 4: Merge BAMs" >> ${PROFILE_LOG}

STEP4_START=$(date +%s.%N)

# Merge keep.bam with filtered remapped reads
${SAMTOOLS} merge -@ ${THREADS} -o ${OUTPUT_DIR}/wasp2py_merged.bam \
    ${KEEP_BAM_ORIG} ${OUTPUT_DIR}/remap_filtered.bam
${SAMTOOLS} sort -@ ${THREADS} -o ${OUTPUT_DIR}/remap_keep.bam ${OUTPUT_DIR}/wasp2py_merged.bam
${SAMTOOLS} index ${OUTPUT_DIR}/remap_keep.bam
rm -f ${OUTPUT_DIR}/wasp2py_merged.bam

STEP4_END=$(date +%s.%N)
STEP4_TIME=$(echo "${STEP4_END} - ${STEP4_START}" | bc)
echo "STEP 4 completed in ${STEP4_TIME} seconds"

# -----------------------------------------------------------------------------
# SUMMARY
# -----------------------------------------------------------------------------
TOTAL_END=$(date +%s.%N)
TOTAL_TIME=$(echo "${TOTAL_END} - ${TOTAL_START}" | bc)

# Count kept reads
KEPT_READS=$(${SAMTOOLS} view -c ${OUTPUT_DIR}/remap_keep.bam)
PASS_RATE=$(echo "scale=2; ${KEPT_READS} * 100 / ${ORIGINAL_READS}" | bc)

# WASP-only time (excluding alignment)
WASP_ONLY=$(echo "${STEP1_TIME} + ${STEP3_TIME} + ${STEP4_TIME}" | bc)

echo "========================================"
echo "BENCHMARK RESULTS"
echo "========================================"
echo ""
echo "Individual step times:"
echo "  STEP 1 (WASP2-Python make-reads):   ${STEP1_TIME} s"
echo "  STEP 2 (BWA-MEM remap):             ${STEP2_TIME} s"
echo "  STEP 3 (WASP2-Python filter):       ${STEP3_TIME} s"
echo "  STEP 4 (Merge BAMs):                ${STEP4_TIME} s"
echo ""
echo "WASP-only time (steps 1+3+4):         ${WASP_ONLY} s"
echo "TOTAL time:                           ${TOTAL_TIME} s"
echo ""
echo "Read counts:"
echo "  Original reads:  ${ORIGINAL_READS}"
echo "  Kept reads:      ${KEPT_READS}"
echo "  Pass rate:       ${PASS_RATE}%"
echo ""
echo "Output directory: ${OUTPUT_DIR}"
echo "End time: $(date)"

# Save results to JSON
cat > ${OUTPUT_DIR}/benchmark_results.json << EOF
{
    "timestamp": "${TIMESTAMP}",
    "pipeline": "wasp2py_snp",
    "sample": "${SAMPLE}",
    "data_type": "ATAC-seq",
    "threads": ${THREADS},
    "step1_wasp2py_make_reads_s": ${STEP1_TIME},
    "step2_bwa_remap_s": ${STEP2_TIME},
    "step3_wasp2py_filter_s": ${STEP3_TIME},
    "step4_merge_bams_s": ${STEP4_TIME},
    "wasp_only_s": ${WASP_ONLY},
    "total_s": ${TOTAL_TIME},
    "original_reads": ${ORIGINAL_READS},
    "kept_reads": ${KEPT_READS},
    "pass_rate_percent": ${PASS_RATE},
    "variant_count": ${VARIANT_COUNT}
}
EOF

echo ""
echo "Results saved to: ${OUTPUT_DIR}/benchmark_results.json"
echo "Key output files for Figure 1C:"
echo "  - original.bam (pre-remap)"
echo "  - remap_keep.bam (post-remap)"
echo "  - variants.bed (SNP positions)"
