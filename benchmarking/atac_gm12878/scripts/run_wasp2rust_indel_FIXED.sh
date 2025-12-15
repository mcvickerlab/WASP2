#!/bin/bash
# WASP2-Rust FIXED benchmark on GM12878 ATAC-seq data - SNP + INDEL mode
# CORRECT pipeline: filter → unified → remap → filter → MERGE
# FIXED: Uses het-only variants (samples parameter triggers het filtering)
#
# This produces EQUIVALENT output to WASP1 (keep + remap_keep merged)
# Plus indel support that WASP1 doesn't have
#
#$ -N wasp2rust_gm12878_indel_fixed
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
FINAL_OUTPUT_DIR="${BENCHMARK_DIR}/results/wasp2rust_indel_fixed_${TIMESTAMP}"
LOCAL_SCRATCH="${TMPDIR:-/tmp}"
WORK_OUTPUT_DIR="${LOCAL_SCRATCH}/wasp2rust_gm12878_indel_${JOB_ID:-local}_${TIMESTAMP}"
OUTPUT_DIR="${WORK_OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}"

# Provenance
GIT_SHA=$(git -C "${WASP2_DIR}" rev-parse HEAD 2>/dev/null || echo "unknown")

# Data files
INPUT_BAM="/iblm/netapp/data1/aho/atac/Buenrostro2013/merged/GM12878_ATACseq_50k/GM12878_ATACseq_50k_merged.sorted.bam"
VCF="/iblm/netapp/data1/aho/variants/NA12878.vcf.gz"
REF_GENOME="/iblm/netapp/data1/aho/ref_genomes/index/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/genome.fa"

# Sample
SAMPLE="NA12878"
THREADS=8
SAME_LOCUS_SLOP=10   # allow positional slop when comparing remap loci (needed for indels)
# Prevent hidden oversubscription (numpy/BLAS) inside helper Python steps.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
# For WASP2 parallel unified pipeline: disable per-worker htslib BAM threads unless explicitly overridden.
export WASP2_BAM_THREADS="${WASP2_BAM_THREADS:-0}"
ENABLE_WASP2_TIMING="${ENABLE_WASP2_TIMING:-0}"
if [ -n "${WASP2_TIMING+x}" ]; then
    ENABLE_WASP2_TIMING=1
fi
if [ "${ENABLE_WASP2_TIMING}" = "1" ]; then
    export WASP2_TIMING=1
fi

# Conda
source /iblm/netapp/home/jjaureguy/mambaforge/etc/profile.d/conda.sh
conda activate WASP2_dev2

# Ensure we pick up the installed wasp2_rust with sidecar binding
export PYTHONPATH="/iblm/netapp/home/jjaureguy/mambaforge/lib/python3.10/site-packages:${PYTHONPATH}"
# Ensure conda-provided OpenSSL is used (pysam depends on it)
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH}"

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
echo "WASP2-Rust INDEL FIXED Pipeline Benchmark - GM12878 ATAC-seq" > ${PROFILE_LOG}
echo "=============================================================" >> ${PROFILE_LOG}
echo "Timestamp: ${TIMESTAMP}" >> ${PROFILE_LOG}
echo "Date: $(date)" >> ${PROFILE_LOG}
echo "Sample: ${SAMPLE}" >> ${PROFILE_LOG}
echo "Threads: ${THREADS}" >> ${PROFILE_LOG}
echo "WASP2_TIMING enabled: ${ENABLE_WASP2_TIMING}" >> ${PROFILE_LOG}
echo "WASP2_BAM_THREADS: ${WASP2_BAM_THREADS}" >> ${PROFILE_LOG}
echo "" >> ${PROFILE_LOG}
echo "FIXED PIPELINE: filter → unified → remap → filter → MERGE" >> ${PROFILE_LOG}
echo "VARIANT MODE: SNPs + Indels (het-only)" >> ${PROFILE_LOG}
echo "" >> ${PROFILE_LOG}

echo "========================================"
echo "WASP2-Rust INDEL FIXED Benchmark on GM12878 ATAC-seq"
echo "Timestamp: ${TIMESTAMP}"
echo "========================================"
echo "CORRECT pipeline: filter → unified → remap → filter → MERGE"
echo "VARIANT MODE: SNPs + Indels (het-only filtering via samples param)"
echo "Start time: $(date)"
echo "Input BAM: ${INPUT_BAM}"
echo "VCF: ${VCF}"
echo "Threads: ${THREADS}"
echo "same_locus_slop (indels): ${SAME_LOCUS_SLOP}"
echo "WASP2_TIMING enabled: ${ENABLE_WASP2_TIMING}"
echo ""

# Verify Rust module
echo "Verifying Rust module..."
${PYTHON} -c "import wasp2_rust; print('wasp2_rust module: AVAILABLE')" || {
    echo "ERROR: wasp2_rust module not available"
    exit 1
}
${PYTHON} -c "from wasp2_rust import filter_bam_by_variants_py, unified_make_reads_parallel_py; print('All functions: AVAILABLE')" || {
    echo "ERROR: Required functions not available"
    exit 1
}
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

# Create variant BED file using Python vcf_to_bed
# KEY FIX: samples=['NA12878'] triggers het-only filtering
# include_indels=True adds indel support
echo "Creating variant BED file using vcf_to_bed (het-only SNPs + indels)..."
${PYTHON} -c "
from mapping.intersect_variant_data import vcf_to_bed

vcf_to_bed(
    vcf_file='${VCF}',
    out_bed='${OUTPUT_DIR}/variants.bed',
    samples=['${SAMPLE}'],      # This triggers het-only filtering!
    include_indels=True,        # Include indels (WASP2 feature)
    max_indel_len=10            # Max indel length
)
print('BED file created successfully (het-only SNPs + indels)')
"

# Count variants
VARIANT_COUNT=$(wc -l < ${OUTPUT_DIR}/variants.bed)
echo "Total het variants in BED: ${VARIANT_COUNT}"

# Count SNPs vs indels separately for reporting
SNP_COUNT=$(awk 'length($4)==1 && length($5)==1' ${OUTPUT_DIR}/variants.bed | wc -l)
INDEL_COUNT=$((VARIANT_COUNT - SNP_COUNT))
echo "  Het SNPs: ${SNP_COUNT}"
echo "  Het Indels: ${INDEL_COUNT}"

# Create separate BED files for post-analysis
awk 'length($4)==1 && length($5)==1' ${OUTPUT_DIR}/variants.bed > ${OUTPUT_DIR}/snps.bed
awk 'length($4)!=1 || length($5)!=1' ${OUTPUT_DIR}/variants.bed > ${OUTPUT_DIR}/indels.bed

# -----------------------------------------------------------------------------
# BENCHMARK STARTS HERE
# -----------------------------------------------------------------------------
echo "========================================"
echo "BENCHMARK TIMING STARTS NOW"
echo "========================================"

TOTAL_START=$(date +%s.%N)

# -----------------------------------------------------------------------------
# STEP 1: Filter BAM by variants → get keep.bam + to_remap.bam
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 1: Filter BAM by variants (split into keep + to_remap)"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 1: filter_bam_by_variants (split BAM)" >> ${PROFILE_LOG}

STEP1_START=$(date +%s.%N)

${TIME} -v ${PYTHON} -c "
import wasp2_rust

remap_count, keep_count, unique_names = wasp2_rust.filter_bam_by_variants_py(
    '${ORIGINAL_BAM}',
    '${OUTPUT_DIR}/variants.bed',
    '${OUTPUT_DIR}/to_remap.bam',
    '${OUTPUT_DIR}/keep.bam',
    is_paired=True,
    threads=${THREADS}
)

print(f'Reads to remap: {remap_count:,}')
print(f'Reads to keep (no variants): {keep_count:,}')
print(f'Unique read names: {unique_names:,}')
" 2>> ${PROFILE_LOG}

${SAMTOOLS} index ${OUTPUT_DIR}/to_remap.bam
${SAMTOOLS} index ${OUTPUT_DIR}/keep.bam

STEP1_END=$(date +%s.%N)
STEP1_TIME=$(echo "${STEP1_END} - ${STEP1_START}" | bc)
echo "STEP 1 completed in ${STEP1_TIME} seconds"

# Count keep reads
KEEP_READS=$(${SAMTOOLS} view -c ${OUTPUT_DIR}/keep.bam)
echo "Keep reads (no variants): ${KEEP_READS}"

# -----------------------------------------------------------------------------
# STEP 2: WASP2-Rust unified make-reads (on to_remap.bam ONLY!)
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 2: WASP2-Rust unified make-reads (on to_remap.bam only)"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 2: unified make-reads (from to_remap.bam)" >> ${PROFILE_LOG}

STEP2_START=$(date +%s.%N)

${TIME} -v ${PYTHON} -c "
from wasp2_rust import unified_make_reads_parallel_py
import json

stats = unified_make_reads_parallel_py(
    '${OUTPUT_DIR}/to_remap.bam',
    '${OUTPUT_DIR}/variants.bed',
    '${OUTPUT_DIR}/remap_r1.fq',
    '${OUTPUT_DIR}/remap_r2.fq',
    max_seqs=64,
    threads=${THREADS},
    compress_output=False,
    indel_mode=True
)

with open('${OUTPUT_DIR}/unified_stats.json', 'w') as f:
    json.dump(stats, f, indent=2)

print(f'Pairs processed: {stats[\"pairs_processed\"]:,}')
print(f'Haplotypes written: {stats[\"haplotypes_written\"]:,}')
" 2>> ${PROFILE_LOG}

STEP2_END=$(date +%s.%N)
STEP2_TIME=$(echo "${STEP2_END} - ${STEP2_START}" | bc)
echo "STEP 2 completed in ${STEP2_TIME} seconds"

${PYTHON} "${WASP2_DIR}/benchmarking/tools/log_unified_stats.py" "${OUTPUT_DIR}/unified_stats.json" \
    --label "atac_indel_step2" >> ${PROFILE_LOG} 2>&1 || echo "WARN: failed to log unified_stats" >> ${PROFILE_LOG}

# Count SNV-only vs INDEL-only vs BOTH overlaps (PRE-filter, from unified stats; convert pairs→reads)
TYPE_PRE=$(${PYTHON} - <<PY
import json
d=json.load(open("${OUTPUT_DIR}/unified_stats.json"))
snv_only = int(d.get("pairs_with_snvs_only", 0)) * 2
indel_only = int(d.get("pairs_with_indels_only", 0)) * 2
both = int(d.get("pairs_with_snvs_and_indels", 0)) * 2
print(f"{snv_only} {indel_only} {both}")
PY
)
read SNV_ONLY_READS_PRE INDEL_ONLY_READS_PRE BOTH_READS_PRE <<< "${TYPE_PRE}"
echo "Overlap type reads PRE-filter:"
echo "  SNV-only reads:   ${SNV_ONLY_READS_PRE}"
echo "  INDEL-only reads: ${INDEL_ONLY_READS_PRE}"
echo "  BOTH reads:       ${BOTH_READS_PRE}"

# Count reads to remap
if [ -f "${OUTPUT_DIR}/remap_r1.fq" ]; then
    REMAP_COUNT=$(wc -l < ${OUTPUT_DIR}/remap_r1.fq)
    REMAP_PAIRS=$((REMAP_COUNT / 4))
    echo "Reads to remap: ${REMAP_PAIRS} pairs"
fi

# -----------------------------------------------------------------------------
# STEP 3: BWA-MEM remap flipped reads
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 3: BWA-MEM remap flipped reads"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 3: BWA-MEM Remap" >> ${PROFILE_LOG}

STEP3_START=$(date +%s.%N)

${TIME} -v ${BWA} mem -t ${THREADS} -M ${REF_GENOME} \
    ${OUTPUT_DIR}/remap_r1.fq ${OUTPUT_DIR}/remap_r2.fq 2>> ${PROFILE_LOG} | \
    ${SAMTOOLS} view -S -b -h -F 4 - > ${OUTPUT_DIR}/remapped_unsorted.bam

${SAMTOOLS} sort -@ ${THREADS} -o ${OUTPUT_DIR}/remapped.bam ${OUTPUT_DIR}/remapped_unsorted.bam
${SAMTOOLS} index ${OUTPUT_DIR}/remapped.bam
rm -f ${OUTPUT_DIR}/remapped_unsorted.bam

STEP3_END=$(date +%s.%N)
STEP3_TIME=$(echo "${STEP3_END} - ${STEP3_START}" | bc)
echo "STEP 3 completed in ${STEP3_TIME} seconds"

# -----------------------------------------------------------------------------
# STEP 4: WASP2-Rust filter-remapped
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 4: WASP2-Rust filter-remapped"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 4: filter_bam_wasp" >> ${PROFILE_LOG}

STEP4_START=$(date +%s.%N)

echo "Sanity check: wasp2_rust signature" >> ${PROFILE_LOG}
${PYTHON} - <<'PY' >> ${PROFILE_LOG} 2>&1
import wasp2_rust, inspect
print("filter_bam_wasp signature:", inspect.signature(wasp2_rust.filter_bam_wasp))
print("has filter_bam_wasp_with_sidecar:", hasattr(wasp2_rust, "filter_bam_wasp_with_sidecar"))
PY

${TIME} -v ${PYTHON} -c "
from wasp2_rust import filter_bam_wasp
from wasp2_rust import filter_bam_wasp_with_sidecar

kept, removed_moved, removed_missing = filter_bam_wasp_with_sidecar(
    '${OUTPUT_DIR}/to_remap.bam',
    '${OUTPUT_DIR}/remapped.bam',
    '${OUTPUT_DIR}/remap_keep.bam',
    threads=${THREADS},
    same_locus_slop=${SAME_LOCUS_SLOP},
    expected_sidecar='${OUTPUT_DIR}/remap_r1.fq.expected_positions.tsv'
)

print(f'Kept: {kept}, Removed (moved): {removed_moved}, Removed (missing): {removed_missing}')
" 2>> ${PROFILE_LOG}

${SAMTOOLS} index ${OUTPUT_DIR}/remap_keep.bam

STEP4_END=$(date +%s.%N)
STEP4_TIME=$(echo "${STEP4_END} - ${STEP4_START}" | bc)
echo "STEP 4 completed in ${STEP4_TIME} seconds"

# Count SNV-only vs INDEL-only vs BOTH overlaps (POST-filter, from sidecar + remap_keep.bam)
TYPE_POST=$(${PYTHON} - <<PY
import pysam
from collections import defaultdict
sidecar_path = "${OUTPUT_DIR}/remap_r1.fq.expected_positions.tsv"
mask_by_name = {}
with open(sidecar_path) as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 6:
            continue
        q = parts[0]
        if "_WASP_" in q:
            base = q.split("_WASP_")[0]
        else:
            base = q.split("/")[0]
        try:
            mask = int(parts[5])
        except Exception:
            mask = 0
        if mask:
            mask_by_name[base] = mask

counts = defaultdict(int)
bam = pysam.AlignmentFile("${OUTPUT_DIR}/remap_keep.bam", "rb")
for rec in bam.fetch(until_eof=True):
    base = rec.query_name
    m = mask_by_name.get(base, 0)
    if m == 1:
        counts["snv"] += 1
    elif m == 2:
        counts["indel"] += 1
    elif m == 3:
        counts["both"] += 1

print(f"{counts['snv']} {counts['indel']} {counts['both']}")
PY
)
read SNV_ONLY_READS_POST INDEL_ONLY_READS_POST BOTH_READS_POST <<< "${TYPE_POST}"
echo "Overlap type reads POST-filter:"
echo "  SNV-only reads:   ${SNV_ONLY_READS_POST}"
echo "  INDEL-only reads: ${INDEL_ONLY_READS_POST}"
echo "  BOTH reads:       ${BOTH_READS_POST}"

# -----------------------------------------------------------------------------
# STEP 5: MERGE keep.bam + remap_keep.bam → wasp_filtered.bam
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 5: Merge keep + remap_keep → wasp_filtered.bam"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 5: Merge BAMs" >> ${PROFILE_LOG}

STEP5_START=$(date +%s.%N)

${TIME} -v ${SAMTOOLS} merge -@ ${THREADS} -o ${OUTPUT_DIR}/wasp_filtered_unsorted.bam \
    ${OUTPUT_DIR}/keep.bam ${OUTPUT_DIR}/remap_keep.bam 2>> ${PROFILE_LOG}

${SAMTOOLS} sort -@ ${THREADS} -o ${OUTPUT_DIR}/wasp_filtered.bam ${OUTPUT_DIR}/wasp_filtered_unsorted.bam
${SAMTOOLS} index ${OUTPUT_DIR}/wasp_filtered.bam
rm -f ${OUTPUT_DIR}/wasp_filtered_unsorted.bam

STEP5_END=$(date +%s.%N)
STEP5_TIME=$(echo "${STEP5_END} - ${STEP5_START}" | bc)
echo "STEP 5 completed in ${STEP5_TIME} seconds"

# -----------------------------------------------------------------------------
# SUMMARY
# -----------------------------------------------------------------------------
TOTAL_END=$(date +%s.%N)
TOTAL_TIME=$(echo "${TOTAL_END} - ${TOTAL_START}" | bc)

# Count final reads
FINAL_READS=$(${SAMTOOLS} view -c ${OUTPUT_DIR}/wasp_filtered.bam)
REMAP_KEEP_READS=$(${SAMTOOLS} view -c ${OUTPUT_DIR}/remap_keep.bam)
PASS_RATE=$(echo "scale=2; ${FINAL_READS} * 100 / ${ORIGINAL_READS}" | bc)

# WASP-only time (steps 1+2+4+5, excluding BWA alignment)
WASP_ONLY=$(echo "${STEP1_TIME} + ${STEP2_TIME} + ${STEP4_TIME} + ${STEP5_TIME}" | bc)

echo "========================================"
echo "BENCHMARK RESULTS (INDEL FIXED PIPELINE)"
echo "========================================"
echo ""
echo "Individual step times:"
echo "  STEP 1 (filter BAM by variants):  ${STEP1_TIME} s"
echo "  STEP 2 (unified make-reads):      ${STEP2_TIME} s"
echo "  STEP 3 (BWA-MEM remap):           ${STEP3_TIME} s"
echo "  STEP 4 (filter remapped):         ${STEP4_TIME} s"
echo "  STEP 5 (merge BAMs):              ${STEP5_TIME} s"
echo ""
echo "WASP-only time (steps 1+2+4+5):     ${WASP_ONLY} s"
echo "TOTAL time:                         ${TOTAL_TIME} s"
echo ""
echo "Read counts:"
echo "  Original reads:       ${ORIGINAL_READS}"
echo "  Keep reads:           ${KEEP_READS}"
echo "  Remap keep reads:     ${REMAP_KEEP_READS}"
echo "  FINAL (merged):       ${FINAL_READS}"
echo "  Pass rate:            ${PASS_RATE}%"
echo ""
echo "Variant counts (het-only):"
echo "  Total variants:       ${VARIANT_COUNT}"
echo "  SNPs:                 ${SNP_COUNT}"
echo "  Indels:               ${INDEL_COUNT}"
echo "same_locus_slop used:   ${SAME_LOCUS_SLOP}"
echo ""
echo "Work directory: ${OUTPUT_DIR}"
echo "Final output directory: ${FINAL_OUTPUT_DIR}"
echo "End time: $(date)"

# Save results to JSON
cat > ${OUTPUT_DIR}/benchmark_results.json << EOF
{
    "timestamp": "${TIMESTAMP}",
    "pipeline": "wasp2rust_indel_FIXED",
    "git_sha": "${GIT_SHA}",
    "sample": "${SAMPLE}",
    "data_type": "ATAC-seq",
    "threads": ${THREADS},
    "include_indels": true,
    "indel_mode": true,
    "cigar_aware_expected_pos": true,
    "sidecar_format_version": 2,
    "het_only": true,
    "step1_filter_bam_s": ${STEP1_TIME},
    "step2_unified_make_reads_s": ${STEP2_TIME},
    "step3_bwa_remap_s": ${STEP3_TIME},
    "step4_filter_remapped_s": ${STEP4_TIME},
    "step5_merge_bams_s": ${STEP5_TIME},
    "wasp_only_s": ${WASP_ONLY},
    "total_s": ${TOTAL_TIME},
    "original_reads": ${ORIGINAL_READS},
    "keep_reads": ${KEEP_READS},
    "remap_keep_reads": ${REMAP_KEEP_READS},
    "final_reads": ${FINAL_READS},
    "pass_rate_percent": ${PASS_RATE},
    "snv_only_overlap_reads_pre": ${SNV_ONLY_READS_PRE},
    "indel_only_overlap_reads_pre": ${INDEL_ONLY_READS_PRE},
    "snv_indel_overlap_reads_pre": ${BOTH_READS_PRE},
    "snv_only_overlap_reads_post": ${SNV_ONLY_READS_POST},
    "indel_only_overlap_reads_post": ${INDEL_ONLY_READS_POST},
    "snv_indel_overlap_reads_post": ${BOTH_READS_POST},
    "total_variants": ${VARIANT_COUNT},
    "snp_count": ${SNP_COUNT},
    "indel_count": ${INDEL_COUNT},
    "same_locus_slop": ${SAME_LOCUS_SLOP}
}
EOF

echo ""
echo "Results saved to: ${OUTPUT_DIR}/benchmark_results.json"
echo ""
echo "COMPARISON TO SNP-ONLY:"
echo "  SNP-only variants:    2,173,342"
echo "  SNP+Indel variants:   ${VARIANT_COUNT}"
echo "  Indel overhead:       ${INDEL_COUNT} additional variants"

# -----------------------------------------------------------------------------
# Copy artifacts back to final output directory (NFS)
# -----------------------------------------------------------------------------
echo ""
echo "Copying benchmark artifacts to final output directory..."
mkdir -p "${FINAL_OUTPUT_DIR}"
for f in \
    "${OUTPUT_DIR}/benchmark_results.json" \
    "${PROFILE_LOG}" \
    "${OUTPUT_DIR}/unified_stats.json" \
    "${OUTPUT_DIR}/remap_keep.bam" \
    "${OUTPUT_DIR}/remap_keep.bam.bai" \
    "${OUTPUT_DIR}/wasp_filtered.bam" \
    "${OUTPUT_DIR}/wasp_filtered.bam.bai"
do
    if [ -f "${f}" ]; then
        rsync -a "${f}" "${FINAL_OUTPUT_DIR}/"
    fi
done
echo "Final output directory: ${FINAL_OUTPUT_DIR}"
cd /
rm -rf "${WORK_OUTPUT_DIR}"
