#!/bin/bash
# WASP2-Rust Benchmark with Multi-Variant Overlap Fix
#
#$ -N wasp2rust_multivar_fix
#$ -V
#$ -pe iblm 8
#$ -l h_vmem=4G
#$ -q iblm.q
#$ -j y
#$ -o /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/atac_gm12878/logs/
#$ -cwd

set -e

# Configuration - use existing input BAM from previous run
INPUT_BAM="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/atac_gm12878/results/wasp2rust_snp_fixed_2025-12-08_08-31-23/original.bam"
INPUT_VCF="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/atac_gm12878/GM12878_chr_formatted.vcf.gz"
REFERENCE="/iblm/netapp/data1/aho/ref_genomes/index/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/genome.fa"
BWA_INDEX="${REFERENCE}"

# Output directory with timestamp
TIMESTAMP=$(date +%Y-%m-%d_%H-%M-%S)
OUTDIR="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/atac_gm12878/results/wasp2rust_multivar_fix_${TIMESTAMP}"
mkdir -p "${OUTDIR}"

THREADS=8
COMPRESSION_THREADS="${COMPRESSION_THREADS:-1}"
SAMPLE="GM12878"

echo "=========================================="
echo "WASP2-Rust Benchmark (Multi-Variant Fix)"
echo "Timestamp: $(date)"
echo "Output: ${OUTDIR}"
echo "=========================================="

cd /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/rust

# Step 1: Use existing variants.bed (VCF to BED was run previously)
echo ""
echo "=========================================="
echo "STEP 1: Using existing variants.bed"
echo "=========================================="
EXISTING_BED="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/atac_gm12878/results/wasp2rust_snp_fixed_2025-12-08_08-31-23/variants.bed"
cp "${EXISTING_BED}" "${OUTDIR}/variants.bed"
echo "  Variants: $(wc -l < ${OUTDIR}/variants.bed)"

# Step 2+3: Unified Make Reads (filter + intersect + remap in single pass)
echo ""
echo "=========================================="
echo "STEP 2+3: Unified Make Reads"
echo "=========================================="
start_time=$(date +%s)

python3 << EOF
import wasp2_rust
stats = wasp2_rust.unified_make_reads_parallel_py(
    "${INPUT_BAM}",
    "${OUTDIR}/variants.bed",
    "${OUTDIR}/remap_r1.fq.gz",
    "${OUTDIR}/remap_r2.fq.gz",
    max_seqs=64,
    threads=${THREADS},
    channel_buffer=50000,
    compression_threads=${COMPRESSION_THREADS},
    compress_output=True
)
print(f"Stats: {stats}")
EOF

end_time=$(date +%s)

# Count FASTQ entries (this is the KEY metric for the fix)
echo ""
echo "  Counting FASTQ entries..."
FASTQ_R1_COUNT=$(zcat "${OUTDIR}/remap_r1.fq.gz" | grep -c "^@" || echo "0")
echo "  R1 FASTQ entries: ${FASTQ_R1_COUNT}"
echo "  Time: $((end_time - start_time))s"

# Step 4: BWA-MEM2 Alignment
echo ""
echo "=========================================="
echo "STEP 4: BWA-MEM Alignment"
echo "=========================================="
start_time=$(date +%s)

bwa mem -t ${THREADS} -M \
    "${BWA_INDEX}" \
    "${OUTDIR}/remap_r1.fq.gz" \
    "${OUTDIR}/remap_r2.fq.gz" | \
    samtools sort -@ ${THREADS} -o "${OUTDIR}/remapped.bam" -

samtools index -@ ${THREADS} "${OUTDIR}/remapped.bam"

end_time=$(date +%s)
REMAPPED_READS=$(samtools view -c "${OUTDIR}/remapped.bam")
echo "  Remapped reads: ${REMAPPED_READS}"
echo "  Time: $((end_time - start_time))s"

# Step 4.5: Create pseudo to_remap.bam from original BAM (extract reads used for remapping)
echo ""
echo "=========================================="
echo "STEP 4.5: Create to_remap.bam from original BAM"
echo "=========================================="
# Extract read names from FASTQ (strip _WASP suffix and /1 or /2)
zcat "${OUTDIR}/remap_r1.fq.gz" | grep "^@" | sed 's/^@//' | sed 's/_WASP_.*//' | sort -u > "${OUTDIR}/remap_readnames.txt"
echo "  Unique read names for remapping: $(wc -l < "${OUTDIR}/remap_readnames.txt")"

# Extract these reads from original BAM to create to_remap.bam
samtools view -N "${OUTDIR}/remap_readnames.txt" -@ ${THREADS} -b -o "${OUTDIR}/to_remap.bam" "${INPUT_BAM}"
samtools index -@ ${THREADS} "${OUTDIR}/to_remap.bam"
echo "  to_remap.bam reads: $(samtools view -c "${OUTDIR}/to_remap.bam")"

# Step 5: Filter Remapped using filter_bam_wasp
echo ""
echo "=========================================="
echo "STEP 5: Filter Remapped (filter_bam_wasp)"
echo "=========================================="
start_time=$(date +%s)

python3 << EOF
from wasp2_rust import filter_bam_wasp
stats = filter_bam_wasp(
    "${OUTDIR}/to_remap.bam",
    "${OUTDIR}/remapped.bam",
    "${OUTDIR}/wasp_filtered.bam",
    threads=${THREADS}
)
print(f"Stats: {stats}")
EOF

samtools index -@ ${THREADS} "${OUTDIR}/wasp_filtered.bam"

end_time=$(date +%s)
FILTERED_READS=$(samtools view -c "${OUTDIR}/wasp_filtered.bam")
echo "  WASP filtered reads: ${FILTERED_READS}"
echo "  Time: $((end_time - start_time))s"

# Summary
echo ""
echo "=========================================="
echo "SUMMARY"
echo "=========================================="
echo "  VCF variants:      $(wc -l < ${OUTDIR}/variants.bed)"
echo "  FASTQ R1 entries:  ${FASTQ_R1_COUNT}"
echo "  Remapped reads:    ${REMAPPED_READS}"
echo "  WASP filtered:     ${FILTERED_READS}"
echo ""
echo "Output directory: ${OUTDIR}"
echo "Completed: $(date)"

# Save summary to file
cat << SUMMARY > "${OUTDIR}/summary.txt"
WASP2-Rust Multi-Variant Fix Benchmark
======================================
Date: $(date)
Input BAM: ${INPUT_BAM}
Input VCF: ${INPUT_VCF}

Results:
  VCF variants:      $(wc -l < ${OUTDIR}/variants.bed)
  FASTQ R1 entries:  ${FASTQ_R1_COUNT}
  Remapped reads:    ${REMAPPED_READS}
  WASP filtered:     ${FILTERED_READS}

Expected (Python DEV baseline):
  VCF variants:      2,175,509
  FASTQ R1 entries:  5,269,562
  WASP filtered:     158,340,301
SUMMARY

echo "Summary saved to ${OUTDIR}/summary.txt"
