#!/bin/bash
# WASP2-Rust (Unified Pipeline) benchmark on GM12878 ATAC-seq data - SNP only mode
# Uses Rust unified_make_reads_parallel for maximum performance
#
# For Figure 1B (timing) and Figure 1C (pre/post remap read counts)
#
#$ -N wasp2rust_gm12878_snp
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
OUTPUT_DIR="${BENCHMARK_DIR}/results/wasp2rust_snp_${TIMESTAMP}"
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

# Set PYTHONPATH
export PYTHONPATH="${WASP2_DIR}/src:${PYTHONPATH}"

cd "${OUTPUT_DIR}"

# Log file
PROFILE_LOG="${OUTPUT_DIR}/profile.log"
echo "WASP2-Rust SNP-only Unified Pipeline Benchmark - GM12878 ATAC-seq" > ${PROFILE_LOG}
echo "===================================================================" >> ${PROFILE_LOG}
echo "Timestamp: ${TIMESTAMP}" >> ${PROFILE_LOG}
echo "Date: $(date)" >> ${PROFILE_LOG}
echo "Sample: ${SAMPLE}" >> ${PROFILE_LOG}
echo "Threads: ${THREADS}" >> ${PROFILE_LOG}
echo "" >> ${PROFILE_LOG}

echo "========================================"
echo "WASP2-Rust SNP Benchmark on GM12878 ATAC-seq"
echo "Timestamp: ${TIMESTAMP}"
echo "========================================"
echo "Start time: $(date)"
echo "Input BAM: ${INPUT_BAM}"
echo "VCF: ${VCF}"
echo "Threads: ${THREADS}"
echo ""

# Verify Rust module
echo "Verifying Rust module..."
${PYTHON} -c "import wasp2_rust; print('wasp2_rust module: AVAILABLE')" || {
    echo "ERROR: wasp2_rust module not available"
    exit 1
}
${PYTHON} -c "from wasp2_rust import unified_make_reads_parallel_py; print('Unified pipeline: AVAILABLE')" || {
    echo "ERROR: Unified pipeline not available"
    exit 1
}
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

# Extract variant BED for Figure 1C analysis (SNPs only)
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
# STEP 1: WASP2-Rust unified make-reads (SNP only)
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 1: WASP2-Rust unified make-reads (SNP only)"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 1: WASP2-Rust unified make-reads" >> ${PROFILE_LOG}

STEP1_START=$(date +%s.%N)

${TIME} -v ${PYTHON} -c "
from mapping.run_mapping import run_make_remap_reads_unified
import json

stats = run_make_remap_reads_unified(
    bam_file='${ORIGINAL_BAM}',
    variant_file='${VCF}',
    samples='${SAMPLE}',
    out_dir='${OUTPUT_DIR}',
    include_indels=False,   # SNPs ONLY
    threads=${THREADS},
    compression_threads=4,
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

# Find FASTQ files
R1_READS=$(ls ${OUTPUT_DIR}/*_remap_r1.fq.gz 2>/dev/null | head -1)
R2_READS=$(ls ${OUTPUT_DIR}/*_remap_r2.fq.gz 2>/dev/null | head -1)

if [ -n "${R1_READS}" ]; then
    REMAP_COUNT=$(zcat ${R1_READS} | wc -l)
    REMAP_PAIRS=$((REMAP_COUNT / 4))
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

# Decompress FASTQs
zcat ${R1_READS} > ${OUTPUT_DIR}/remap_r1.fq
zcat ${R2_READS} > ${OUTPUT_DIR}/remap_r2.fq

${TIME} -v ${BWA} mem -t ${THREADS} -M ${REF_GENOME} \
    ${OUTPUT_DIR}/remap_r1.fq ${OUTPUT_DIR}/remap_r2.fq 2>> ${PROFILE_LOG} | \
    ${SAMTOOLS} view -S -b -h -F 4 - > ${OUTPUT_DIR}/remapped_unsorted.bam

${SAMTOOLS} sort -@ ${THREADS} -o ${OUTPUT_DIR}/remapped.bam ${OUTPUT_DIR}/remapped_unsorted.bam
${SAMTOOLS} index ${OUTPUT_DIR}/remapped.bam
rm -f ${OUTPUT_DIR}/remap_r1.fq ${OUTPUT_DIR}/remap_r2.fq ${OUTPUT_DIR}/remapped_unsorted.bam

STEP2_END=$(date +%s.%N)
STEP2_TIME=$(echo "${STEP2_END} - ${STEP2_START}" | bc)
echo "STEP 2 completed in ${STEP2_TIME} seconds"

# -----------------------------------------------------------------------------
# STEP 3: WASP2-Rust filter-remapped
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 3: WASP2-Rust filter-remapped"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 3: WASP2-Rust filter-remapped" >> ${PROFILE_LOG}

STEP3_START=$(date +%s.%N)

${TIME} -v ${PYTHON} -c "
from wasp2_rust import filter_bam_wasp

kept, removed_moved, removed_missing = filter_bam_wasp(
    '${ORIGINAL_BAM}',
    '${OUTPUT_DIR}/remapped.bam',
    '${OUTPUT_DIR}/remap_keep.bam',
    threads=${THREADS}
)

print(f'Kept: {kept}, Removed (moved): {removed_moved}, Removed (missing): {removed_missing}')
" 2>> ${PROFILE_LOG}

${SAMTOOLS} index ${OUTPUT_DIR}/remap_keep.bam

STEP3_END=$(date +%s.%N)
STEP3_TIME=$(echo "${STEP3_END} - ${STEP3_START}" | bc)
echo "STEP 3 completed in ${STEP3_TIME} seconds"

# -----------------------------------------------------------------------------
# FIGURE 1C: SNV-specific read counting
# -----------------------------------------------------------------------------
echo "========================================"
echo "FIGURE 1C: Counting SNV-overlapping reads"
echo "========================================"

BEDTOOLS=$(which bedtools)

# Count unique read pairs overlapping SNVs in original BAM (pre-WASP)
echo "Counting SNV-overlapping reads in original BAM..."
SNV_READS_PRE=$(${BEDTOOLS} intersect -a ${ORIGINAL_BAM} -b ${OUTPUT_DIR}/variants.bed -u | \
    ${SAMTOOLS} view | cut -f1 | sort -u | wc -l)
echo "SNV-overlapping read pairs (pre-WASP): ${SNV_READS_PRE}"

# Count unique read pairs overlapping SNVs in filtered BAM (post-WASP)
echo "Counting SNV-overlapping reads in filtered BAM..."
SNV_READS_POST=$(${BEDTOOLS} intersect -a ${OUTPUT_DIR}/remap_keep.bam -b ${OUTPUT_DIR}/variants.bed -u | \
    ${SAMTOOLS} view | cut -f1 | sort -u | wc -l)
echo "SNV-overlapping read pairs (post-WASP): ${SNV_READS_POST}"

# Calculate SNV pass rate
if [ ${SNV_READS_PRE} -gt 0 ]; then
    SNV_PASS_RATE=$(echo "scale=2; ${SNV_READS_POST} * 100 / ${SNV_READS_PRE}" | bc)
else
    SNV_PASS_RATE=0
fi
echo "SNV read pass rate: ${SNV_PASS_RATE}%"

# -----------------------------------------------------------------------------
# SUMMARY
# -----------------------------------------------------------------------------
TOTAL_END=$(date +%s.%N)
TOTAL_TIME=$(echo "${TOTAL_END} - ${TOTAL_START}" | bc)

# Count kept reads
KEPT_READS=$(${SAMTOOLS} view -c ${OUTPUT_DIR}/remap_keep.bam)
PASS_RATE=$(echo "scale=2; ${KEPT_READS} * 100 / ${ORIGINAL_READS}" | bc)

# WASP-only time (excluding alignment)
WASP_ONLY=$(echo "${STEP1_TIME} + ${STEP3_TIME}" | bc)

echo "========================================"
echo "BENCHMARK RESULTS"
echo "========================================"
echo ""
echo "Individual step times:"
echo "  STEP 1 (WASP2-Rust make-reads):  ${STEP1_TIME} s"
echo "  STEP 2 (BWA-MEM remap):          ${STEP2_TIME} s"
echo "  STEP 3 (WASP2-Rust filter):      ${STEP3_TIME} s"
echo ""
echo "WASP-only time (steps 1+3):        ${WASP_ONLY} s"
echo "TOTAL time:                        ${TOTAL_TIME} s"
echo ""
echo "Read counts (total):"
echo "  Original reads:  ${ORIGINAL_READS}"
echo "  Kept reads:      ${KEPT_READS}"
echo "  Pass rate:       ${PASS_RATE}%"
echo ""
echo "SNV-overlapping reads (for Figure 1C):"
echo "  Pre-WASP:   ${SNV_READS_PRE} read pairs"
echo "  Post-WASP:  ${SNV_READS_POST} read pairs"
echo "  Pass rate:  ${SNV_PASS_RATE}%"
echo ""
echo "Output directory: ${OUTPUT_DIR}"
echo "End time: $(date)"

# Save results to JSON
cat > ${OUTPUT_DIR}/benchmark_results.json << EOF
{
    "timestamp": "${TIMESTAMP}",
    "pipeline": "wasp2rust_snp",
    "sample": "${SAMPLE}",
    "data_type": "ATAC-seq",
    "threads": ${THREADS},
    "step1_wasp2rust_make_reads_s": ${STEP1_TIME},
    "step2_bwa_remap_s": ${STEP2_TIME},
    "step3_wasp2rust_filter_s": ${STEP3_TIME},
    "wasp_only_s": ${WASP_ONLY},
    "total_s": ${TOTAL_TIME},
    "original_reads": ${ORIGINAL_READS},
    "kept_reads": ${KEPT_READS},
    "pass_rate_percent": ${PASS_RATE},
    "snv_count": ${VARIANT_COUNT},
    "snv_reads_pre": ${SNV_READS_PRE},
    "snv_reads_post": ${SNV_READS_POST},
    "snv_pass_rate_percent": ${SNV_PASS_RATE}
}
EOF

echo ""
echo "Results saved to: ${OUTPUT_DIR}/benchmark_results.json"
echo "Key output files for Figure 1C:"
echo "  - original.bam (pre-remap)"
echo "  - remap_keep.bam (post-remap)"
echo "  - variants.bed (SNP positions)"
