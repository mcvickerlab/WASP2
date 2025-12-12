#!/bin/bash
# Count SNV-overlapping reads for WASP2-Python benchmarks (Figure 1C)
#
#$ -N count_snv_py
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
echo "SNV-overlapping Read Counting for WASP2-Python"
echo "Timestamp: $(date)"
echo "=========================================="

# ==========================================
# WASP2-Python Master (single-threaded)
# ==========================================
echo ""
echo "=========================================="
echo "WASP2-Python Master Branch"
echo "=========================================="

PY_MASTER_DIR="${RESULTS_BASE}/wasp2python_snp_fixed_2025-12-07_13-12-00"

if [ -f "${PY_MASTER_DIR}/original.bam" ] && [ -f "${PY_MASTER_DIR}/het_only_snps.vcf.gz" ]; then
    echo "Found WASP2-Python master files"

    # Create BED from VCF if not exists
    BCFTOOLS=$(which bcftools)
    if [ ! -f "${PY_MASTER_DIR}/variants.bed" ]; then
        echo "Creating variants.bed from VCF..."
        ${BCFTOOLS} query -f '%CHROM\t%POS0\t%END\n' ${PY_MASTER_DIR}/het_only_snps.vcf.gz > ${PY_MASTER_DIR}/variants.bed
    fi

    # Count pre-WASP
    echo "Counting pre-WASP SNV reads..."
    PY_MASTER_PRE=$(${BEDTOOLS} intersect -a ${PY_MASTER_DIR}/original.bam -b ${PY_MASTER_DIR}/variants.bed -u | \
        ${SAMTOOLS} view | cut -f1 | sort -u | wc -l)
    echo "  Pre-WASP SNV reads: ${PY_MASTER_PRE}"

    # Count post-WASP
    echo "Counting post-WASP SNV reads..."
    PY_MASTER_POST=$(${BEDTOOLS} intersect -a ${PY_MASTER_DIR}/wasp_filtered.bam -b ${PY_MASTER_DIR}/variants.bed -u | \
        ${SAMTOOLS} view | cut -f1 | sort -u | wc -l)
    echo "  Post-WASP SNV reads: ${PY_MASTER_POST}"

    # Calculate pass rate
    PY_MASTER_RATE=$(echo "scale=2; ${PY_MASTER_POST} * 100 / ${PY_MASTER_PRE}" | bc)
    echo "  SNV read pass rate: ${PY_MASTER_RATE}%"

    # Update benchmark_results.json
    PY_MASTER_JSON="${PY_MASTER_DIR}/benchmark_results.json"
    if [ -f "${PY_MASTER_JSON}" ]; then
        python3 << EOF
import json

with open("${PY_MASTER_JSON}", 'r') as f:
    data = json.load(f)

data['snv_reads_pre'] = ${PY_MASTER_PRE}
data['snv_reads_post'] = ${PY_MASTER_POST}
data['snv_read_pass_rate'] = ${PY_MASTER_RATE}

with open("${PY_MASTER_JSON}", 'w') as f:
    json.dump(data, f, indent=4)
    f.write('\n')

print("Updated ${PY_MASTER_JSON}")
EOF
    fi
else
    echo "ERROR: WASP2-Python master files not found!"
fi

# ==========================================
# WASP2-Python Dev MT (multithreaded)
# ==========================================
echo ""
echo "=========================================="
echo "WASP2-Python Dev Branch (Multithreaded)"
echo "=========================================="

PY_DEV_DIR="${RESULTS_BASE}/wasp2python_snp_DEV_MT_2025-12-07_15-52-30"

if [ -f "${PY_DEV_DIR}/original.bam" ] && [ -f "${PY_DEV_DIR}/het_only_snps.vcf.gz" ]; then
    echo "Found WASP2-Python dev MT files"

    # Create BED from VCF if not exists
    BCFTOOLS=$(which bcftools)
    if [ ! -f "${PY_DEV_DIR}/variants.bed" ]; then
        echo "Creating variants.bed from VCF..."
        ${BCFTOOLS} query -f '%CHROM\t%POS0\t%END\n' ${PY_DEV_DIR}/het_only_snps.vcf.gz > ${PY_DEV_DIR}/variants.bed
    fi

    # Count pre-WASP
    echo "Counting pre-WASP SNV reads..."
    PY_DEV_PRE=$(${BEDTOOLS} intersect -a ${PY_DEV_DIR}/original.bam -b ${PY_DEV_DIR}/variants.bed -u | \
        ${SAMTOOLS} view | cut -f1 | sort -u | wc -l)
    echo "  Pre-WASP SNV reads: ${PY_DEV_PRE}"

    # Count post-WASP
    echo "Counting post-WASP SNV reads..."
    PY_DEV_POST=$(${BEDTOOLS} intersect -a ${PY_DEV_DIR}/wasp_filtered.bam -b ${PY_DEV_DIR}/variants.bed -u | \
        ${SAMTOOLS} view | cut -f1 | sort -u | wc -l)
    echo "  Post-WASP SNV reads: ${PY_DEV_POST}"

    # Calculate pass rate
    PY_DEV_RATE=$(echo "scale=2; ${PY_DEV_POST} * 100 / ${PY_DEV_PRE}" | bc)
    echo "  SNV read pass rate: ${PY_DEV_RATE}%"

    # Update benchmark_results.json
    PY_DEV_JSON="${PY_DEV_DIR}/benchmark_results.json"
    if [ -f "${PY_DEV_JSON}" ]; then
        python3 << EOF
import json

with open("${PY_DEV_JSON}", 'r') as f:
    data = json.load(f)

data['snv_reads_pre'] = ${PY_DEV_PRE}
data['snv_reads_post'] = ${PY_DEV_POST}
data['snv_read_pass_rate'] = ${PY_DEV_RATE}

with open("${PY_DEV_JSON}", 'w') as f:
    json.dump(data, f, indent=4)
    f.write('\n')

print("Updated ${PY_DEV_JSON}")
EOF
    fi
else
    echo "ERROR: WASP2-Python dev MT files not found!"
fi

# ==========================================
# Summary
# ==========================================
echo ""
echo "=========================================="
echo "SUMMARY - Figure 1C SNV Read Counts"
echo "=========================================="
echo ""
echo "WASP2-Python Master:"
echo "  Pre-WASP:  ${PY_MASTER_PRE:-N/A}"
echo "  Post-WASP: ${PY_MASTER_POST:-N/A}"
echo "  Pass rate: ${PY_MASTER_RATE:-N/A}%"
echo ""
echo "WASP2-Python Dev MT:"
echo "  Pre-WASP:  ${PY_DEV_PRE:-N/A}"
echo "  Post-WASP: ${PY_DEV_POST:-N/A}"
echo "  Pass rate: ${PY_DEV_RATE:-N/A}%"
echo ""
echo "Completed: $(date)"
