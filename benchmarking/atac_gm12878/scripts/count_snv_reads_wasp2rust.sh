#!/bin/bash
# Count SNV-overlapping reads for WASP2-Rust benchmarks (Figure 1C)
#
#$ -N count_snv
#$ -V
#$ -pe iblm 4
#$ -l h_vmem=24G
#$ -j y
#$ -o /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/atac_gm12878/logs/
#$ -cwd

set -e

# Conda
source /iblm/netapp/home/jjaureguy/mambaforge/etc/profile.d/conda.sh
conda activate WASP2_dev2

BEDTOOLS=$(which bedtools)
SAMTOOLS=$(which samtools)

RESULTS_BASE="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/atac_gm12878/results"

echo "=========================================="
echo "SNV-overlapping Read Counting for Figure 1C"
echo "Timestamp: $(date)"
echo "=========================================="

# Function to count SNV-overlapping reads
count_snv_reads() {
    local BAM=$1
    local BED=$2
    local LABEL=$3

    echo "Counting SNV-overlapping reads in ${LABEL}..."
    local COUNT=$(${BEDTOOLS} intersect -a ${BAM} -b ${BED} -u | \
        ${SAMTOOLS} view | cut -f1 | sort -u | wc -l)
    echo "  ${LABEL}: ${COUNT} read pairs"
    echo ${COUNT}
}

# ==========================================
# WASP2-Rust SNP
# ==========================================
echo ""
echo "=========================================="
echo "WASP2-Rust SNP Benchmark"
echo "=========================================="

RUST_SNP_DIR="${RESULTS_BASE}/wasp2rust_snp_fixed_2025-12-07_19-02-45"

if [ -f "${RUST_SNP_DIR}/original.bam" ] && [ -f "${RUST_SNP_DIR}/variants.bed" ]; then
    echo "Found WASP2-Rust SNP files"

    # Count pre-WASP
    RUST_SNP_PRE=$(${BEDTOOLS} intersect -a ${RUST_SNP_DIR}/original.bam -b ${RUST_SNP_DIR}/variants.bed -u | \
        ${SAMTOOLS} view | cut -f1 | sort -u | wc -l)
    echo "  Pre-WASP SNV reads: ${RUST_SNP_PRE}"

    # Count post-WASP
    RUST_SNP_POST=$(${BEDTOOLS} intersect -a ${RUST_SNP_DIR}/wasp_filtered.bam -b ${RUST_SNP_DIR}/variants.bed -u | \
        ${SAMTOOLS} view | cut -f1 | sort -u | wc -l)
    echo "  Post-WASP SNV reads: ${RUST_SNP_POST}"

    # Calculate pass rate
    RUST_SNP_RATE=$(echo "scale=2; ${RUST_SNP_POST} * 100 / ${RUST_SNP_PRE}" | bc)
    echo "  SNV read pass rate: ${RUST_SNP_RATE}%"

    # Update benchmark_results.json
    RUST_SNP_JSON="${RUST_SNP_DIR}/benchmark_results.json"
    if [ -f "${RUST_SNP_JSON}" ]; then
        # Create updated JSON with SNV counts
        python3 << EOF
import json

with open("${RUST_SNP_JSON}", 'r') as f:
    data = json.load(f)

data['snv_reads_pre'] = ${RUST_SNP_PRE}
data['snv_reads_post'] = ${RUST_SNP_POST}
data['snv_read_pass_rate'] = ${RUST_SNP_RATE}

with open("${RUST_SNP_JSON}", 'w') as f:
    json.dump(data, f, indent=4)
    f.write('\n')

print("Updated ${RUST_SNP_JSON}")
EOF
    fi
else
    echo "ERROR: WASP2-Rust SNP files not found!"
fi

# ==========================================
# WASP2-Rust Indel
# ==========================================
echo ""
echo "=========================================="
echo "WASP2-Rust Indel Benchmark"
echo "=========================================="

RUST_INDEL_DIR="${RESULTS_BASE}/wasp2rust_indel_fixed_2025-12-07_19-02-46"

if [ -f "${RUST_INDEL_DIR}/original.bam" ] && [ -f "${RUST_INDEL_DIR}/snps.bed" ]; then
    echo "Found WASP2-Rust Indel files"

    # For indel benchmark, we count reads overlapping SNPs (snps.bed)
    # This is consistent with Figure 1C which shows SNV retention

    # Count pre-WASP (using snps.bed for SNV counting)
    RUST_INDEL_PRE=$(${BEDTOOLS} intersect -a ${RUST_INDEL_DIR}/original.bam -b ${RUST_INDEL_DIR}/snps.bed -u | \
        ${SAMTOOLS} view | cut -f1 | sort -u | wc -l)
    echo "  Pre-WASP SNV reads: ${RUST_INDEL_PRE}"

    # Count post-WASP
    RUST_INDEL_POST=$(${BEDTOOLS} intersect -a ${RUST_INDEL_DIR}/wasp_filtered.bam -b ${RUST_INDEL_DIR}/snps.bed -u | \
        ${SAMTOOLS} view | cut -f1 | sort -u | wc -l)
    echo "  Post-WASP SNV reads: ${RUST_INDEL_POST}"

    # Calculate pass rate
    RUST_INDEL_RATE=$(echo "scale=2; ${RUST_INDEL_POST} * 100 / ${RUST_INDEL_PRE}" | bc)
    echo "  SNV read pass rate: ${RUST_INDEL_RATE}%"

    # Update benchmark_results.json
    RUST_INDEL_JSON="${RUST_INDEL_DIR}/benchmark_results.json"
    if [ -f "${RUST_INDEL_JSON}" ]; then
        python3 << EOF
import json

with open("${RUST_INDEL_JSON}", 'r') as f:
    data = json.load(f)

data['snv_reads_pre'] = ${RUST_INDEL_PRE}
data['snv_reads_post'] = ${RUST_INDEL_POST}
data['snv_read_pass_rate'] = ${RUST_INDEL_RATE}

with open("${RUST_INDEL_JSON}", 'w') as f:
    json.dump(data, f, indent=4)
    f.write('\n')

print("Updated ${RUST_INDEL_JSON}")
EOF
    fi
else
    echo "ERROR: WASP2-Rust Indel files not found!"
fi

# ==========================================
# Summary
# ==========================================
echo ""
echo "=========================================="
echo "SUMMARY - Figure 1C SNV Read Counts"
echo "=========================================="
echo ""
echo "WASP2-Rust SNP:"
echo "  Pre-WASP:  ${RUST_SNP_PRE:-N/A}"
echo "  Post-WASP: ${RUST_SNP_POST:-N/A}"
echo "  Pass rate: ${RUST_SNP_RATE:-N/A}%"
echo ""
echo "WASP2-Rust Indel:"
echo "  Pre-WASP:  ${RUST_INDEL_PRE:-N/A}"
echo "  Post-WASP: ${RUST_INDEL_POST:-N/A}"
echo "  Pass rate: ${RUST_INDEL_RATE:-N/A}%"
echo ""
echo "Completed: $(date)"
