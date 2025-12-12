#!/bin/bash
# WASP2-Rust INDEL Benchmark for ATAC-seq (GM12878)
# Uses new INDEL mode with coordinated R1/R2 trim combinations
# Sample: GM12878 (50k ATAC-seq dataset, ~150M reads)
# Uses run_make_remap_reads_unified with include_indels=True
#
#$ -N wasp2rust_indel_atac
#$ -V
#$ -pe iblm 8
#$ -l h_vmem=32G
#$ -j y
#$ -o /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/atac_gm12878/logs/
#$ -cwd

set -e

# Timestamp
TIMESTAMP=$(date +%Y-%m-%d_%H-%M-%S)
echo "Benchmark timestamp: ${TIMESTAMP}"

# Paths
WASP2_DIR="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp"
BENCHMARK_DIR="${WASP2_DIR}/benchmarking/atac_gm12878"
OUTPUT_DIR="${BENCHMARK_DIR}/results/wasp2rust_indel_atac_${TIMESTAMP}"

# Data files - use pre-aligned BAM (same as existing ATAC-seq benchmarks)
INPUT_BAM="/iblm/netapp/data1/aho/atac/Buenrostro2013/merged/GM12878_ATACseq_50k/GM12878_ATACseq_50k_merged.sorted.bam"
VCF="/iblm/netapp/data1/aho/variants/NA12878.vcf.gz"
REF_GENOME="/iblm/netapp/data1/aho/ref_genomes/index/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/genome.fa"

SAMPLE="NA12878"
THREADS=8
MAX_INDEL_LEN=10  # Maximum INDEL length to include

# Conda
source /iblm/netapp/home/jjaureguy/mambaforge/etc/profile.d/conda.sh
conda activate WASP2_dev2

# Tools
BWA=$(which bwa)
SAMTOOLS=$(which samtools)
PYTHON=$(which python)
TIME="/usr/bin/time"

# Set PYTHONPATH for WASP2 src AND Rust module
export PYTHONPATH="${WASP2_DIR}/src:${WASP2_DIR}/rust:${PYTHONPATH}"

mkdir -p "${OUTPUT_DIR}"
cd "${OUTPUT_DIR}"

# Log file
PROFILE_LOG="${OUTPUT_DIR}/profile.log"
echo "WASP2-Rust INDEL Benchmark Profile (ATAC-seq)" > ${PROFILE_LOG}
echo "Sample: GM12878 ATAC-seq (50k, ~150M reads)" >> ${PROFILE_LOG}
echo "Timestamp: ${TIMESTAMP}" >> ${PROFILE_LOG}
echo "Threads: ${THREADS}" >> ${PROFILE_LOG}
echo "INDEL Support: ENABLED (max ${MAX_INDEL_LEN}bp)" >> ${PROFILE_LOG}
echo "Pipeline: run_make_remap_reads_unified with include_indels=True" >> ${PROFILE_LOG}

echo "========================================"
echo "WASP2-Rust INDEL ATAC-seq Benchmark"
echo "Using coordinated R1/R2 trim combinations"
echo "Sample: GM12878"
echo "Timestamp: ${TIMESTAMP}"
echo "INDEL Support: ENABLED (max ${MAX_INDEL_LEN}bp)"
echo "========================================"
echo "Start time: $(date)"
echo "Input BAM: ${INPUT_BAM}"
echo "VCF: ${VCF}"
echo "Threads: ${THREADS}"
echo "BWA: ${BWA}"

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

# Copy original BAM for reference
echo "Copying original BAM for reference..."
cp ${INPUT_BAM} ${OUTPUT_DIR}/original.bam
${SAMTOOLS} index ${OUTPUT_DIR}/original.bam
ORIGINAL_BAM="${OUTPUT_DIR}/original.bam"
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

# Create variant count from unified stats
VARIANT_COUNT=$(${PYTHON} -c "import json; d=json.load(open('${OUTPUT_DIR}/unified_stats.json')); print(d.get('variants_used', 'N/A'))" 2>/dev/null || echo "N/A")
echo "Variant count: ${VARIANT_COUNT}"

# -----------------------------------------------------------------------------
# STEP 2: BWA remap flipped reads
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 2: BWA remap flipped reads"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 2: BWA Remap" >> ${PROFILE_LOG}

STEP2_START=$(date +%s.%N)

# Decompress for BWA (BWA doesn't read .gz well in all versions)
R1_UNZIPPED="${OUTPUT_DIR}/remap_r1.fq"
R2_UNZIPPED="${OUTPUT_DIR}/remap_r2.fq"
zcat ${R1_READS} > ${R1_UNZIPPED}
zcat ${R2_READS} > ${R2_UNZIPPED}

${TIME} -v ${BWA} mem -t ${THREADS} -M ${REF_GENOME} ${R1_UNZIPPED} ${R2_UNZIPPED} 2>> ${PROFILE_LOG} | \
    ${SAMTOOLS} sort -@ ${THREADS} -o ${OUTPUT_DIR}/remapped.bam - 2>> ${PROFILE_LOG}

REMAPPED_BAM="${OUTPUT_DIR}/remapped.bam"
${SAMTOOLS} index ${REMAPPED_BAM}

# Clean up unzipped FASTQs
rm -f ${R1_UNZIPPED} ${R2_UNZIPPED}

STEP2_END=$(date +%s.%N)
STEP2_TIME=$(echo "${STEP2_END} - ${STEP2_START}" | bc)
echo "STEP 2 completed in ${STEP2_TIME} seconds"

# Check mapping stats
REMAPPED_READS=$(${SAMTOOLS} view -c ${REMAPPED_BAM})
REMAPPED_PROPER=$(${SAMTOOLS} view -c -f 2 ${REMAPPED_BAM})
echo "Remapped reads: ${REMAPPED_READS}"
echo "Remapped proper pairs: ${REMAPPED_PROPER}"

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
from wasp2_rust import filter_bam_wasp

# filter_bam_wasp strips _WASP_* suffix from remapped reads
# to match against original read names in ORIGINAL_BAM
kept, removed_moved, removed_missing = filter_bam_wasp(
    '${ORIGINAL_BAM}',
    '${REMAPPED_BAM}',
    '${OUTPUT_DIR}/remap_keep.bam',
    threads=${THREADS}
)

print(f'Kept: {kept}, Removed (moved): {removed_moved}, Removed (missing): {removed_missing}')
" 2>> ${PROFILE_LOG}

${SAMTOOLS} index ${OUTPUT_DIR}/remap_keep.bam

STEP3_END=$(date +%s.%N)
STEP3_TIME=$(echo "${STEP3_END} - ${STEP3_START}" | bc)
echo "STEP 3 completed in ${STEP3_TIME} seconds"

REMAP_KEEP_READS=$(${SAMTOOLS} view -c ${OUTPUT_DIR}/remap_keep.bam 2>/dev/null || echo "0")
echo "Remap keep reads: ${REMAP_KEEP_READS}"

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
TOTAL_TIME=$(echo "${STEP1_TIME} + ${STEP2_TIME} + ${STEP3_TIME}" | bc)
WASP_ONLY=$(echo "${STEP1_TIME} + ${STEP3_TIME}" | bc)

if [ "${ORIGINAL_READS}" -gt 0 ]; then
    OVERALL_PASS_RATE=$(echo "scale=2; ${REMAP_KEEP_READS} * 100 / ${ORIGINAL_READS}" | bc)
else
    OVERALL_PASS_RATE="0"
fi

echo "========================================"
echo "BENCHMARK RESULTS (WASP2-Rust INDEL ATAC)"
echo "========================================"
echo ""
echo "Individual step times:"
echo "  STEP 1 (unified-make-reads w/INDELs): ${STEP1_TIME} s"
echo "  STEP 2 (BWA remap):                   ${STEP2_TIME} s"
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
echo "  Remapped reads:       ${REMAPPED_READS}"
echo "  Remap keep reads:     ${REMAP_KEEP_READS}"
echo "  Overall pass rate:    ${OVERALL_PASS_RATE}%"
echo ""
echo "Variant counts:"
echo "  Total variants:       ${VARIANT_COUNT}"
echo ""
echo "Output directory: ${OUTPUT_DIR}"
echo "End time: $(date)"

# Save results to JSON
cat > ${OUTPUT_DIR}/benchmark_results.json << EOF
{
    "timestamp": "${TIMESTAMP}",
    "pipeline": "wasp2rust_indel_atac",
    "method": "run_make_remap_reads_unified with include_indels=True",
    "sample": "${SAMPLE}",
    "data_type": "ATAC-seq",
    "aligner": "BWA",
    "threads": ${THREADS},
    "include_indels": true,
    "max_indel_len": ${MAX_INDEL_LEN},
    "step1_unified_make_reads_s": ${STEP1_TIME},
    "step2_bwa_remap_s": ${STEP2_TIME},
    "step3_filter_wasp_s": ${STEP3_TIME},
    "wasp_only_s": ${WASP_ONLY},
    "total_s": ${TOTAL_TIME},
    "snv_reads_pre": ${REMAP_PAIRS},
    "snv_reads_post": ${SNV_READS_POST},
    "snv_pass_rate_percent": ${SNV_PASS_RATE},
    "original_reads": ${ORIGINAL_READS},
    "remapped_reads": ${REMAPPED_READS},
    "remap_keep_reads": ${REMAP_KEEP_READS},
    "overall_pass_rate_percent": ${OVERALL_PASS_RATE},
    "variant_count": "${VARIANT_COUNT}"
}
EOF

echo ""
echo "Results saved to: ${OUTPUT_DIR}/benchmark_results.json"

echo "Done!"
