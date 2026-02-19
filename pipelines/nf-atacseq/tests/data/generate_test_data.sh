#!/bin/bash
# =============================================================================
# WASP2 nf-atacseq Test Data Generator
# =============================================================================
# Creates ATAC-seq-like test data by symlinking shared core data and generating
# pipeline-specific files (shorter fragment FASTQs, BWA index, samplesheet).
#
# Prerequisites: samtools, bgzip, tabix, wgsim, bwa (WASP2_dev2 conda env)
#
# Usage:
#   cd pipelines/nf-atacseq/tests/data
#   bash generate_test_data.sh
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

SHARED_DATA="../../../../tests/shared_data"

echo "==================================================================="
echo " WASP2 nf-atacseq Test Data Generator"
echo "==================================================================="

# Validate shared core data exists
if [[ ! -f "$SHARED_DATA/chr_test.fa" ]]; then
    echo "ERROR: Shared core data not found at $SHARED_DATA"
    echo "  Run: cd tests/shared_data && bash generate_core_data.sh"
    exit 1
fi

# -----------------------------------------------------------------------------
# Symlink shared reference and variants
# -----------------------------------------------------------------------------
echo "[1/4] Symlinking shared reference data..."

for f in chr_test.fa chr_test.fa.fai variants.vcf.gz variants.vcf.gz.tbi annotation.gtf regions.bed; do
    if [[ ! -e "$f" ]]; then
        ln -sf "$SHARED_DATA/$f" "$f"
        echo "  ✓ Linked $f"
    else
        echo "  - $f already exists"
    fi
done

echo ""

# -----------------------------------------------------------------------------
# Simulate ATAC-seq-like reads (shorter fragments, 150-250bp)
# -----------------------------------------------------------------------------
echo "[2/4] Simulating ATAC-seq reads..."

NUM_READS=500
READ_LEN=75
FRAG_SIZE=180
FRAG_STD=30
ERROR_RATE=0.001
SEED=100

if [[ -f "sample1_R1.fq.gz" && -f "sample1_R2.fq.gz" ]]; then
    echo "  FASTQs already exist, skipping"
else
    wgsim -N $NUM_READS \
          -1 $READ_LEN \
          -2 $READ_LEN \
          -r 0 -R 0 -X 0 \
          -e $ERROR_RATE \
          -S $SEED \
          -d $FRAG_SIZE \
          -s $FRAG_STD \
          "$SHARED_DATA/chr_test.fa" \
          sample1_R1.fq \
          sample1_R2.fq \
          > /dev/null 2>&1

    gzip -f sample1_R1.fq
    gzip -f sample1_R2.fq
    echo "  ✓ Created sample1_R{1,2}.fq.gz (${NUM_READS} pairs, ${READ_LEN}bp, ${FRAG_SIZE}bp frags)"
fi

echo ""

# -----------------------------------------------------------------------------
# Build BWA index (for local testing)
# -----------------------------------------------------------------------------
echo "[3/4] Building BWA index..."

BWA_INDEX_DIR="bwa_index"
if [[ -f "${BWA_INDEX_DIR}/chr_test.fa.bwt" ]]; then
    echo "  BWA index already exists, skipping"
else
    mkdir -p "$BWA_INDEX_DIR"
    cp "$SHARED_DATA/chr_test.fa" "$BWA_INDEX_DIR/"
    bwa index "$BWA_INDEX_DIR/chr_test.fa" 2>&1 | tail -2
    echo "  ✓ Created BWA index ($(du -sh $BWA_INDEX_DIR | cut -f1))"
fi

echo ""

# -----------------------------------------------------------------------------
# Create test samplesheet
# -----------------------------------------------------------------------------
echo "[4/4] Creating test samplesheet..."

SAMPLESHEET="samplesheet_test.csv"
if [[ -f "$SAMPLESHEET" ]]; then
    echo "  $SAMPLESHEET already exists, skipping"
else
    cat > "$SAMPLESHEET" << EOF
sample,fastq_1,fastq_2,sample_name
test_sample1,${SCRIPT_DIR}/sample1_R1.fq.gz,${SCRIPT_DIR}/sample1_R2.fq.gz,SAMPLE1
EOF
    echo "  ✓ Created $SAMPLESHEET"
fi

echo ""
echo "==================================================================="
echo " SUCCESS! nf-atacseq test data generated."
echo "==================================================================="
echo "Total: $(du -sh . | cut -f1)"
echo ""
