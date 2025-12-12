#!/bin/bash
# Compare WASP2-Rust vs WASP2-Python DEV outputs at each step
#
#$ -N compare_chr22
#$ -V
#$ -l h_vmem=8G
#$ -j y
#$ -o /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/atac_gm12878/logs/
#$ -cwd

set -e

WORKDIR="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/atac_gm12878/chr22_comparison"
RUST_DIR="${WORKDIR}/rust"
PYTHON_DIR="${WORKDIR}/python"
ANALYSIS_DIR="${WORKDIR}/analysis"

mkdir -p "${ANALYSIS_DIR}"
cd "${WORKDIR}"

echo "========================================"
echo "WASP2 Pipeline Comparison"
echo "Timestamp: $(date)"
echo "========================================"

# Create comparison report
REPORT="${ANALYSIS_DIR}/comparison_report.txt"
echo "WASP2-Rust vs WASP2-Python DEV Comparison" > "${REPORT}"
echo "Generated: $(date)" >> "${REPORT}"
echo "========================================" >> "${REPORT}"

echo ""
echo "========================================"
echo "STEP 1: VCF to BED Comparison"
echo "========================================"
echo "" >> "${REPORT}"
echo "STEP 1: VCF to BED" >> "${REPORT}"
echo "----------------------------------------" >> "${REPORT}"

RUST_BED="${RUST_DIR}/step1_vcf_to_bed/variants.bed"
PYTHON_BED="${PYTHON_DIR}/step1_vcf_to_bed/variants.bed"

if [[ -f "${RUST_BED}" && -f "${PYTHON_BED}" ]]; then
    rust_count=$(wc -l < "${RUST_BED}")
    python_count=$(wc -l < "${PYTHON_BED}")

    echo "  Rust variants:   ${rust_count}"
    echo "  Python variants: ${python_count}"
    echo "Rust variants:   ${rust_count}" >> "${REPORT}"
    echo "Python variants: ${python_count}" >> "${REPORT}"

    if diff -q "${RUST_BED}" "${PYTHON_BED}" > /dev/null 2>&1; then
        echo "  STATUS: IDENTICAL"
        echo "STATUS: IDENTICAL" >> "${REPORT}"
    else
        echo "  STATUS: DIFFERENT"
        echo "STATUS: DIFFERENT" >> "${REPORT}"
        diff "${RUST_BED}" "${PYTHON_BED}" > "${ANALYSIS_DIR}/step1_bed_diff.txt" 2>&1 || true
        head -20 "${ANALYSIS_DIR}/step1_bed_diff.txt" >> "${REPORT}" 2>/dev/null || true
    fi
else
    echo "  WARNING: One or both BED files missing"
    echo "WARNING: Missing files" >> "${REPORT}"
fi

echo ""
echo "========================================"
echo "STEP 2: Filter BAM Comparison"
echo "========================================"
echo "" >> "${REPORT}"
echo "STEP 2: Filter BAM" >> "${REPORT}"
echo "----------------------------------------" >> "${REPORT}"

RUST_NAMES="${RUST_DIR}/step2_filter_bam/read_names.txt"
PYTHON_NAMES="${PYTHON_DIR}/step2_filter_bam/read_names.txt"

if [[ -f "${RUST_NAMES}" && -f "${PYTHON_NAMES}" ]]; then
    rust_reads=$(wc -l < "${RUST_NAMES}")
    python_reads=$(wc -l < "${PYTHON_NAMES}")

    echo "  Rust unique reads to remap:   ${rust_reads}"
    echo "  Python unique reads to remap: ${python_reads}"
    echo "Rust unique reads:   ${rust_reads}" >> "${REPORT}"
    echo "Python unique reads: ${python_reads}" >> "${REPORT}"

    sort "${RUST_NAMES}" > "${ANALYSIS_DIR}/rust_names_sorted.txt"
    sort "${PYTHON_NAMES}" > "${ANALYSIS_DIR}/python_names_sorted.txt"

    comm -23 "${ANALYSIS_DIR}/rust_names_sorted.txt" "${ANALYSIS_DIR}/python_names_sorted.txt" > "${ANALYSIS_DIR}/step2_rust_only.txt"
    rust_only=$(wc -l < "${ANALYSIS_DIR}/step2_rust_only.txt")

    comm -13 "${ANALYSIS_DIR}/rust_names_sorted.txt" "${ANALYSIS_DIR}/python_names_sorted.txt" > "${ANALYSIS_DIR}/step2_python_only.txt"
    python_only=$(wc -l < "${ANALYSIS_DIR}/step2_python_only.txt")

    comm -12 "${ANALYSIS_DIR}/rust_names_sorted.txt" "${ANALYSIS_DIR}/python_names_sorted.txt" > "${ANALYSIS_DIR}/step2_common.txt"
    common=$(wc -l < "${ANALYSIS_DIR}/step2_common.txt")

    echo "  Common reads:       ${common}"
    echo "  Rust-only reads:    ${rust_only}"
    echo "  Python-only reads:  ${python_only}"
    echo "Common reads:       ${common}" >> "${REPORT}"
    echo "Rust-only reads:    ${rust_only}" >> "${REPORT}"
    echo "Python-only reads:  ${python_only}" >> "${REPORT}"
else
    echo "  WARNING: One or both read name files missing"
    echo "WARNING: Missing files" >> "${REPORT}"
fi

echo ""
echo "========================================"
echo "STEP 3: Make Reads Comparison"
echo "========================================"
echo "" >> "${REPORT}"
echo "STEP 3: Make Reads (CRITICAL)" >> "${REPORT}"
echo "----------------------------------------" >> "${REPORT}"

RUST_WASP="${RUST_DIR}/step3_make_reads/wasp_names.txt"
PYTHON_WASP="${PYTHON_DIR}/step3_make_reads/wasp_names.txt"

if [[ -f "${RUST_WASP}" && -f "${PYTHON_WASP}" ]]; then
    rust_wasp=$(wc -l < "${RUST_WASP}")
    python_wasp=$(wc -l < "${PYTHON_WASP}")

    echo "  Rust WASP reads:   ${rust_wasp}"
    echo "  Python WASP reads: ${python_wasp}"
    echo "Rust WASP reads:   ${rust_wasp}" >> "${REPORT}"
    echo "Python WASP reads: ${python_wasp}" >> "${REPORT}"

    # Analyze haplotype totals
    awk -F'_' '{print $NF}' "${RUST_WASP}" | sort | uniq -c | sort -rn > "${ANALYSIS_DIR}/rust_haplotype_totals.txt"
    awk -F'_' '{print $NF}' "${PYTHON_WASP}" | sort | uniq -c | sort -rn > "${ANALYSIS_DIR}/python_haplotype_totals.txt"

    echo "  Rust haplotype total distribution:"
    head -5 "${ANALYSIS_DIR}/rust_haplotype_totals.txt"
    echo "  Python haplotype total distribution:"
    head -5 "${ANALYSIS_DIR}/python_haplotype_totals.txt"
else
    echo "  WARNING: One or both WASP name files missing"
    echo "WARNING: Missing files" >> "${REPORT}"
fi

echo ""
echo "========================================"
echo "STEP 5: Filter Remapped Comparison"
echo "========================================"
echo "" >> "${REPORT}"
echo "STEP 5: Filter Remapped (CRITICAL)" >> "${REPORT}"
echo "----------------------------------------" >> "${REPORT}"

RUST_KEPT="${RUST_DIR}/step5_filter_remapped/kept_names.txt"
PYTHON_KEPT="${PYTHON_DIR}/step5_filter_remapped/kept_read_names.txt"

if [[ -f "${RUST_KEPT}" && -f "${PYTHON_KEPT}" ]]; then
    rust_kept=$(wc -l < "${RUST_KEPT}")
    python_kept=$(wc -l < "${PYTHON_KEPT}")

    echo "  Rust kept unique reads:   ${rust_kept}"
    echo "  Python kept unique reads: ${python_kept}"
    echo "Rust kept reads:   ${rust_kept}" >> "${REPORT}"
    echo "Python kept reads: ${python_kept}" >> "${REPORT}"

    sort "${RUST_KEPT}" > "${ANALYSIS_DIR}/rust_kept_sorted.txt"
    sort "${PYTHON_KEPT}" > "${ANALYSIS_DIR}/python_kept_sorted.txt"

    comm -23 "${ANALYSIS_DIR}/rust_kept_sorted.txt" "${ANALYSIS_DIR}/python_kept_sorted.txt" > "${ANALYSIS_DIR}/step5_rust_keeps_python_removes.txt"
    rust_only_kept=$(wc -l < "${ANALYSIS_DIR}/step5_rust_keeps_python_removes.txt")

    comm -13 "${ANALYSIS_DIR}/rust_kept_sorted.txt" "${ANALYSIS_DIR}/python_kept_sorted.txt" > "${ANALYSIS_DIR}/step5_python_keeps_rust_removes.txt"
    python_only_kept=$(wc -l < "${ANALYSIS_DIR}/step5_python_keeps_rust_removes.txt")

    echo "  Rust keeps, Python removes: ${rust_only_kept}"
    echo "  Python keeps, Rust removes: ${python_only_kept}"
    echo "Rust-only kept: ${rust_only_kept}" >> "${REPORT}"
    echo "Python-only kept: ${python_only_kept}" >> "${REPORT}"

    gap=$((python_kept - rust_kept))
    if [[ ${python_kept} -gt 0 ]]; then
        gap_pct=$(echo "scale=2; ${gap} * 100 / ${python_kept}" | bc)
        echo ""
        echo "  GAP: Python retains ${gap} more reads (${gap_pct}%)"
        echo "GAP: ${gap} reads (${gap_pct}%)" >> "${REPORT}"
    fi
else
    echo "  WARNING: One or both kept name files missing"
    echo "WARNING: Missing files" >> "${REPORT}"
fi

echo ""
echo "========================================"
echo "Summary Report"
echo "========================================"
cat "${REPORT}"

echo ""
echo "Detailed analysis files in: ${ANALYSIS_DIR}"
ls -la "${ANALYSIS_DIR}"

echo ""
echo "Completed: $(date)"
