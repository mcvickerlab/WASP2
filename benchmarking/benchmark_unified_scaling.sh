#!/bin/bash
#$ -N wasp2_unified_scaling
#$ -t 1-150
#$ -tc 10
#$ -V
#$ -pe iblm 8
#$ -l h_vmem=8G
#$ -j y
#$ -o /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/logs/
#$ -cwd

# =============================================================================
# WASP2 UNIFIED PIPELINE SCALING BENCHMARK - 1M to 150M reads
# Uses the fast single-pass unified pipeline (Rust throughout, no fallbacks)
# =============================================================================
# Task ID determines read count: task 1 = 1M, task 15 = 15M, task 150 = 150M
#
# Threading:
#   - SGE: 8 cores (-pe iblm 8)
#   - Unified pipeline: 8 threads + configurable compression threads
#   - BWA mem: -t 16
#   - filter: 8 threads
#
# Pipeline:
#   1. Subset BAM to N reads
#   2. WASP2 unified make-reads (single-pass Rust)
#   3. BWA remap
#   4. WASP2 filter-remapped (Rust)
# =============================================================================

set -e

source /iblm/netapp/home/jjaureguy/mambaforge/etc/profile.d/conda.sh
conda activate WASP2_dev2

# Threading controls
# - THREADS defaults to NSLOTS (SGE), else 8
# - COMPRESSION_THREADS is per FASTQ file (R1 and R2); defaults to 1 to avoid oversubscription
THREADS="${THREADS:-${NSLOTS:-8}}"
COMPRESSION_THREADS="${COMPRESSION_THREADS:-1}"

# Use a stable run ID shared across all array tasks.
# Prefer RUN_TIMESTAMP if passed at qsub time, else fall back to JOB_ID (shared
# across tasks), else a wall-clock timestamp (interactive use).
RUN_TIMESTAMP="${RUN_TIMESTAMP:-}"
if [ -z "${RUN_TIMESTAMP}" ]; then
    RUN_TIMESTAMP="${JOB_ID:-$(date +%Y-%m-%d_%H-%M-%S)}"
fi
TIMESTAMP="${RUN_TIMESTAMP}"

WASP2_DIR="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp"

# Timestamped results file
RESULTS_DIR="${WASP2_DIR}/benchmarking/results/unified_scaling_${TIMESTAMP}"
mkdir -p "${RESULTS_DIR}"
log_file="${RESULTS_DIR}/wasp2_unified_scaling.tsv"

# Write header exactly once across all tasks (atomic mkdir)
if mkdir "${RESULTS_DIR}/.init" 2>/dev/null; then
    echo -e "timestamp\tn_reads\tseed\ttotal_s\tunified_s\tremap_s\tfilter_s" > "${log_file}"
fi

# Same inputs as aho/SNP benchmark
genome_index="/iblm/netapp/data1/aho/ref_genomes/index/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/genome.fa"
input_bam="/iblm/netapp/data1/aho/atac/Buenrostro2013/merged/GM12878_ATACseq_50k/GM12878_ATACseq_50k_merged.sorted.bam"
input_vcf="/iblm/netapp/data1/aho/variants/NA12878.vcf.gz"
sample="NA12878"

# Read count scales with task ID
n_subset=$(($SGE_TASK_ID * 1000000))
bam_reads=159133307
bam_prop=$(echo "scale=7; ${n_subset}/${bam_reads}" | bc)
seed=$RANDOM

dir=$(mktemp -d)
trap 'rm -rf "$dir"' EXIT

echo "=============================================="
echo "WASP2 UNIFIED SCALING - Task ${SGE_TASK_ID}"
echo "Timestamp: ${TIMESTAMP}"
echo "Reads: ${n_subset} (${SGE_TASK_ID}M)"
echo "Pipeline: UNIFIED (single-pass Rust)"
echo "Date: $(date)"
echo "=============================================="

# Verify Rust module
echo "Verifying Rust module..."
python -c "import wasp2_rust; print('wasp2_rust module: AVAILABLE')" || {
    echo "ERROR: wasp2_rust module not available"
    exit 1
}

python -c "from wasp2_rust import unified_make_reads_parallel_py; print('Unified pipeline: AVAILABLE')" || {
    echo "ERROR: Unified pipeline not available"
    exit 1
}

# Subset BAM
prefix="bam_${n_subset}_${seed}"
subset_bam="${dir}/${prefix}.bam"
echo "Subsetting BAM to ${n_subset} reads..."
samtools view -h --bam --subsample ${bam_prop} --subsample-seed ${seed} -o ${subset_bam} ${input_bam}
samtools index ${subset_bam}

# Total time tracking
export PYTHONPATH="${WASP2_DIR}/src:${PYTHONPATH}"
start=$(date +%s)

# -----------------------------------------------------------------------------
# STEP 1: WASP2 unified make-reads (single-pass Rust)
# -----------------------------------------------------------------------------
echo ""
echo "STEP 1: WASP2 unified make-reads..."
step1_start=$(date +%s)

python -c "
from mapping.run_mapping import run_make_remap_reads_unified
import json

stats = run_make_remap_reads_unified(
    bam_file='${subset_bam}',
    variant_file='${input_vcf}',
    samples='${sample}',
    out_dir='${dir}',
    threads=${THREADS},
    compression_threads=${COMPRESSION_THREADS},
    use_parallel=True
)

# Save stats
with open('${dir}/unified_stats.json', 'w') as f:
    json.dump(stats, f, indent=2)

print(f'Pairs processed: {stats[\"pairs_processed\"]:,}')
print(f'Haplotypes written: {stats[\"haplotypes_written\"]:,}')
"

unified_end=$(date +%s)
unified_runtime=$((unified_end - step1_start))
echo "Unified make-reads completed in ${unified_runtime}s"

# Find output FASTQ files
r1_reads=$(ls ${dir}/*_remap_r1.fq.gz 2>/dev/null | head -1)
r2_reads=$(ls ${dir}/*_remap_r2.fq.gz 2>/dev/null | head -1)

if [ -z "${r1_reads}" ] || [ -z "${r2_reads}" ]; then
    echo "ERROR: Remap FASTQ files not found"
    ls -la ${dir}
    exit 1
fi

# -----------------------------------------------------------------------------
# STEP 2: BWA remap
# -----------------------------------------------------------------------------
echo ""
echo "STEP 2: BWA remap..."
step2_start=$(date +%s)

# Decompress for BWA
zcat ${r1_reads} > ${dir}/remap_r1.fq
zcat ${r2_reads} > ${dir}/remap_r2.fq

remapped_bam="${dir}/${prefix}_remapped.bam"
bwa mem -t 16 -M ${genome_index} ${dir}/remap_r1.fq ${dir}/remap_r2.fq | \
    samtools view -S -b -h -F 4 - > ${remapped_bam}
samtools sort -o ${remapped_bam} ${remapped_bam}
samtools index ${remapped_bam}

# Clean up uncompressed FASTQs
rm -f ${dir}/remap_r1.fq ${dir}/remap_r2.fq

remap_end=$(date +%s)
remap_runtime=$((remap_end - step2_start))
echo "BWA remap completed in ${remap_runtime}s"

# -----------------------------------------------------------------------------
# STEP 3: WASP2 filter-remapped (Rust)
# -----------------------------------------------------------------------------
echo ""
echo "STEP 3: WASP2 filter-remapped..."
step3_start=$(date +%s)

python -c "
from wasp2_rust import filter_bam_wasp

kept, removed_moved, removed_missing = filter_bam_wasp(
    '${subset_bam}',
    '${remapped_bam}',
    '${dir}/remap_keep.bam',
    threads=8
)

print(f'Kept: {kept}, Removed: {removed_moved}, Missing: {removed_missing}')
"

end=$(date +%s)
filt_runtime=$((end - step3_start))
total_runtime=$((end - start))

echo ""
echo "=============================================="
echo "RESULTS - ${SGE_TASK_ID}M reads (UNIFIED)"
echo "=============================================="
echo "  Unified make-reads: ${unified_runtime}s"
echo "  BWA remap:          ${remap_runtime}s"
echo "  Filter:             ${filt_runtime}s"
echo "  TOTAL:              ${total_runtime}s"
echo "=============================================="

# Write results
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${TIMESTAMP}" "${n_subset}" "${seed}" \
    "${total_runtime}" "${unified_runtime}" "${remap_runtime}" "${filt_runtime}" >> ${log_file}

echo "Task ${SGE_TASK_ID} complete"
echo "Results: ${log_file}"
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"
