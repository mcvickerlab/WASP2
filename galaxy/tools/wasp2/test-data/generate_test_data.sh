#!/bin/bash
# =============================================================================
# WASP2 Galaxy Tool Test Data Generator
# =============================================================================
# Creates test fixtures for Galaxy tool XML tests. Galaxy's planemo test runner
# looks for files in test-data/ relative to the XML files.
#
# Required test files (from Galaxy XML test sections):
#   test.bam        - Aligned reads (from shared sample1.bam)
#   test.vcf        - Uncompressed VCF (from shared variants.vcf)
#   remapped.bam    - Simulated remapped BAM (copy of test.bam)
#   wasp_data.json  - WASP mapping metadata
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
echo "[1/3] Copying test data from shared core..."

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
    echo "  ✓ Copied test.vcf"
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
# Create WASP data JSON (minimal valid structure)
# -----------------------------------------------------------------------------
echo "[2/3] Creating WASP data JSON..."

if [[ ! -f "wasp_data.json" ]]; then
    cat > wasp_data.json << 'EOJSON'
{
    "keep_bam": "test_keep.bam",
    "to_remap_bam": "test_to_remap.bam",
    "to_remap_fq1": "test_to_remap_R1.fq.gz",
    "to_remap_fq2": "test_to_remap_R2.fq.gz",
    "remap_num_reads": 50,
    "keep_num_reads": 450,
    "variant_source": "test.vcf",
    "samples": ["SAMPLE1"],
    "is_paired": true,
    "is_phased": true
}
EOJSON
    echo "  ✓ Created wasp_data.json"
else
    echo "  - wasp_data.json already exists"
fi

echo ""

# -----------------------------------------------------------------------------
# Create counts TSV (minimal valid structure for find-imbalance test)
# -----------------------------------------------------------------------------
echo "[3/3] Creating counts TSV..."

if [[ ! -f "counts.tsv" ]]; then
    cat > counts.tsv << 'EOTSV'
chrom	pos	ref	alt	genotype	ref_count	alt_count	total_count	sample
chr_test	750	C	T	0/1	15	12	27	SAMPLE1
chr_test	1200	T	G	0/1	20	18	38	SAMPLE1
chr_test	2800	A	C	0/1	8	14	22	SAMPLE1
chr_test	3200	G	A	0/1	25	22	47	SAMPLE1
chr_test	5000	G	T	0/1	11	13	24	SAMPLE1
chr_test	10800	T	C	0/1	18	15	33	SAMPLE1
chr_test	11200	A	G	0/1	22	19	41	SAMPLE1
chr_test	12800	C	A	0/1	16	20	36	SAMPLE1
chr_test	13200	G	T	0/1	14	12	26	SAMPLE1
chr_test	15000	A	C	0/1	19	17	36	SAMPLE1
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
