#!/bin/bash
# =============================================================================
# WASP2 nf-outrider Test Data Generator
# =============================================================================
# Creates OUTRIDER pipeline test data by symlinking 3 BAMs from shared core
# (OUTRIDER requires >= 3 samples) plus annotation and variant data.
#
# Prerequisites: Shared core data must exist
#
# Usage:
#   cd pipelines/nf-outrider/tests/data
#   bash generate_test_data.sh
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

SHARED_DATA="../../../../tests/shared_data"

echo "==================================================================="
echo " WASP2 nf-outrider Test Data Generator"
echo "==================================================================="

# Validate shared core data exists
if [[ ! -f "$SHARED_DATA/sample1.bam" ]]; then
    echo "ERROR: Shared core data not found at $SHARED_DATA"
    echo "  Run: cd tests/shared_data && bash generate_core_data.sh"
    exit 1
fi

# -----------------------------------------------------------------------------
# Symlink shared BAMs, variants, and annotation
# -----------------------------------------------------------------------------
echo "[1/2] Symlinking shared data..."

# BAM files (3 samples for OUTRIDER minimum)
for i in 1 2 3; do
    for ext in bam bam.bai; do
        src="$SHARED_DATA/sample${i}.${ext}"
        dst="sample${i}.${ext}"
        if [[ ! -e "$dst" ]]; then
            ln -sf "$src" "$dst"
            echo "  ✓ Linked $dst"
        else
            echo "  - $dst already exists"
        fi
    done
done

# Variants and annotation
for f in variants.vcf.gz variants.vcf.gz.tbi annotation.gtf; do
    if [[ ! -e "$f" ]]; then
        ln -sf "$SHARED_DATA/$f" "$f"
        echo "  ✓ Linked $f"
    else
        echo "  - $f already exists"
    fi
done

echo ""

# -----------------------------------------------------------------------------
# Create test samplesheet
# -----------------------------------------------------------------------------
echo "[2/2] Creating test samplesheet..."

SAMPLESHEET="samplesheet_test.csv"
if [[ -f "$SAMPLESHEET" ]]; then
    echo "  $SAMPLESHEET already exists, skipping"
else
    cat > "$SAMPLESHEET" << EOF
sample,bam,bai
sample1,${SCRIPT_DIR}/sample1.bam,${SCRIPT_DIR}/sample1.bam.bai
sample2,${SCRIPT_DIR}/sample2.bam,${SCRIPT_DIR}/sample2.bam.bai
sample3,${SCRIPT_DIR}/sample3.bam,${SCRIPT_DIR}/sample3.bam.bai
EOF
    echo "  ✓ Created $SAMPLESHEET"
fi

echo ""
echo "==================================================================="
echo " SUCCESS! nf-outrider test data generated."
echo "==================================================================="
echo "Total: $(du -sh . | cut -f1)"
echo ""
