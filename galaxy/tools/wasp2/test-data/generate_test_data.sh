#!/bin/bash
# =============================================================================
# WASP2 Galaxy Tool Test Data Generator
# =============================================================================
# Creates test fixtures for Galaxy tool XML tests. Galaxy's planemo test runner
# looks for files in test-data/ relative to the XML files.
#
# Required test files (from Galaxy XML test sections):
#   test.bam        - Aligned reads (from shared sample1.bam)
#   test.vcf        - Uncompressed phased VCF (from shared variants.vcf)
#   remapped.bam    - Simulated remapped BAM (copy of test.bam)
#   counts.tsv      - WASP2 counting output
#
# Prerequisites: samtools (WASP2_dev2 conda env)
#
# Usage:
#   cd galaxy/tools/wasp2/test-data
#   bash generate_test_data.sh
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

SHARED_DATA="../../../../tests/shared_data"

echo "==================================================================="
echo " WASP2 Galaxy Tool Test Data Generator"
echo "==================================================================="

# Validate shared core data exists
if [[ ! -f "$SHARED_DATA/sample1.bam" ]]; then
    echo "ERROR: Shared core data not found at $SHARED_DATA"
    echo "  Run: cd tests/shared_data && bash generate_core_data.sh"
    exit 1
fi

# -----------------------------------------------------------------------------
# Copy test files (Galaxy needs actual files, not symlinks)
# -----------------------------------------------------------------------------
echo "[1/2] Copying test data from shared core..."

# BAM for count-variants and make-reads tests
if [[ ! -f "test.bam" ]]; then
    cp "$SHARED_DATA/sample1.bam" test.bam
    cp "$SHARED_DATA/sample1.bam.bai" test.bam.bai
    echo "  ✓ Copied test.bam + .bai"
else
    echo "  - test.bam already exists"
fi

# Uncompressed VCF for Galaxy (Galaxy handles format internally)
if [[ ! -f "test.vcf" ]]; then
    cp "$SHARED_DATA/variants.vcf" test.vcf
    sed -i 's#0/1#0|1#g' test.vcf
    echo "  ✓ Copied and phased test.vcf"
else
    echo "  - test.vcf already exists"
fi

# Remapped BAM for filter-remapped test
if [[ ! -f "remapped.bam" ]]; then
    cp "$SHARED_DATA/sample1.bam" remapped.bam
    cp "$SHARED_DATA/sample1.bam.bai" remapped.bam.bai
    echo "  ✓ Created remapped.bam (copy of sample1.bam)"
else
    echo "  - remapped.bam already exists"
fi

echo ""

# -----------------------------------------------------------------------------
# Create counts TSV (minimal valid structure for find-imbalance test)
# -----------------------------------------------------------------------------
echo "[2/2] Creating counts TSV..."

if [[ ! -f "counts.tsv" ]]; then
    cat > counts.tsv << 'EOTSV'
chrom	pos0	pos	ref	alt	GT	ref_count	alt_count	other_count	sample
chr_test	749	750	C	T	0|1	15	12	0	SAMPLE1
chr_test	1199	1200	T	G	0|1	20	18	0	SAMPLE1
chr_test	2799	2800	A	C	0|1	8	14	0	SAMPLE1
chr_test	3199	3200	G	A	0|1	25	22	0	SAMPLE1
chr_test	4999	5000	G	T	0|1	11	13	0	SAMPLE1
chr_test	10799	10800	T	C	1|0	18	15	0	SAMPLE1
chr_test	11199	11200	A	G	1|0	22	19	0	SAMPLE1
chr_test	12799	12800	C	A	1|0	16	20	0	SAMPLE1
chr_test	13199	13200	G	T	1|0	14	12	0	SAMPLE1
chr_test	14999	15000	A	C	1|0	19	17	0	SAMPLE1
EOTSV
    echo "  ✓ Created counts.tsv (10 variants)"
else
    echo "  - counts.tsv already exists"
fi

echo ""
echo "==================================================================="
echo " SUCCESS! Galaxy test data generated."
echo "==================================================================="
echo "Total: $(du -sh . | cut -f1)"
echo ""
echo "To test Galaxy tools:"
echo "  planemo test galaxy/tools/wasp2/"
echo ""
