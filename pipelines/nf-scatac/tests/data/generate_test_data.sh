#!/bin/bash
# =============================================================================
# WASP2 nf-scatac Test Data Generator
# =============================================================================
# Creates single-cell ATAC-seq test data: synthetic fragment file with cell
# barcodes, peaks BED, and local samplesheet.
#
# Prerequisites: samtools, bgzip, tabix (WASP2_dev2 conda env)
#
# Usage:
#   cd pipelines/nf-scatac/tests/data
#   bash generate_test_data.sh
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

SHARED_DATA="../../../../tests/shared_data"

echo "==================================================================="
echo " WASP2 nf-scatac Test Data Generator"
echo "==================================================================="

# Validate shared core data exists
if [[ ! -f "$SHARED_DATA/chr_test.fa" ]]; then
    echo "ERROR: Shared core data not found at $SHARED_DATA"
    echo "  Run: cd tests/shared_data && bash generate_core_data.sh"
    exit 1
fi

# -----------------------------------------------------------------------------
# Symlink shared variant data
# -----------------------------------------------------------------------------
echo "[1/5] Symlinking shared data..."

for f in variants.vcf.gz variants.vcf.gz.tbi chr_test.fa chr_test.fa.fai; do
    if [[ ! -e "$f" ]]; then
        ln -sf "$SHARED_DATA/$f" "$f"
        echo "  ✓ Linked $f"
    else
        echo "  - $f already exists"
    fi
done

echo ""

# -----------------------------------------------------------------------------
# Create cell barcodes
# -----------------------------------------------------------------------------
echo "[2/5] Creating cell barcodes..."

if [[ -f "barcodes.txt" ]]; then
    echo "  barcodes.txt already exists, skipping"
else
    cat > barcodes.txt << 'EOBC'
AAACCTGAGAAACCAT
AAACCTGAGAAACCTA
AAACCTGAGAAACGAG
AAACCTGAGAAACTGT
AAACCTGAGAAAGACA
AAACCTGCAGTCAGCC
AAACCTGCATACCATG
AAACCTGCATGCTAGT
AAACCTGTCAGTCCCT
AAACCTGTCATGTAGC
EOBC
    echo "  ✓ Created barcodes.txt (10 cell barcodes)"
fi

echo ""

# -----------------------------------------------------------------------------
# Create peaks BED from gene regions
# -----------------------------------------------------------------------------
echo "[3/5] Creating peaks BED..."

if [[ -f "peaks.bed" ]]; then
    echo "  peaks.bed already exists, skipping"
else
    # ATAC-seq peaks at promoter/gene body regions
    cat > peaks.bed << 'EOPEAKS'
chr_test	400	1600	peak1	500	.
chr_test	2400	3600	peak2	400	.
chr_test	4400	5600	peak3	350	.
chr_test	10400	11600	peak4	450	.
chr_test	12400	13600	peak5	380	.
chr_test	14400	15600	peak6	420	.
EOPEAKS
    echo "  ✓ Created peaks.bed (6 peaks)"
fi

echo ""

# -----------------------------------------------------------------------------
# Generate synthetic fragment file
# -----------------------------------------------------------------------------
echo "[4/5] Generating synthetic fragments..."

if [[ -f "fragments.tsv.gz" ]]; then
    echo "  fragments.tsv.gz already exists, skipping"
else
    # Create synthetic fragment file mimicking 10x scATAC-seq format:
    # chr  start  end  barcode  duplicate_count
    # Place fragments at variant positions to ensure they're captured
    BARCODES=($(cat barcodes.txt))
    NUM_BARCODES=${#BARCODES[@]}

    # Generate ~100 fragments spread across peaks
    {
        # Fragments covering Gene 1 SNPs
        for i in $(seq 0 $((NUM_BARCODES - 1))); do
            bc=${BARCODES[$i]}
            # Fragments near SNP positions with realistic sizes (150-500bp)
            echo -e "chr_test\t650\t950\t${bc}\t1"
            echo -e "chr_test\t1100\t1400\t${bc}\t1"
            echo -e "chr_test\t2700\t3000\t${bc}\t1"
            echo -e "chr_test\t4900\t5200\t${bc}\t1"
        done

        # Fragments covering Gene 2 SNPs
        for i in $(seq 0 $((NUM_BARCODES - 1))); do
            bc=${BARCODES[$i]}
            echo -e "chr_test\t10700\t11000\t${bc}\t1"
            echo -e "chr_test\t12700\t13000\t${bc}\t1"
            echo -e "chr_test\t14900\t15200\t${bc}\t1"
        done
    } | sort -k1,1 -k2,2n > fragments.tsv

    bgzip -c fragments.tsv > fragments.tsv.gz
    tabix -p bed fragments.tsv.gz
    rm -f fragments.tsv
    echo "  ✓ Created fragments.tsv.gz + .tbi ($(du -h fragments.tsv.gz | cut -f1))"
fi

echo ""

# -----------------------------------------------------------------------------
# Create test samplesheet
# -----------------------------------------------------------------------------
echo "[5/5] Creating test samplesheet..."

SAMPLESHEET="samplesheet_local.csv"
if [[ -f "$SAMPLESHEET" ]]; then
    echo "  $SAMPLESHEET already exists, skipping"
else
    cat > "$SAMPLESHEET" << EOF
sample,fragments,cellranger_dir,bam,barcode_tag,chemistry,barcodes,peaks
test_sample,${SCRIPT_DIR}/fragments.tsv.gz,,,CB,10x-atac-v2,${SCRIPT_DIR}/barcodes.txt,${SCRIPT_DIR}/peaks.bed
EOF
    echo "  ✓ Created $SAMPLESHEET"
fi

echo ""
echo "==================================================================="
echo " SUCCESS! nf-scatac test data generated."
echo "==================================================================="
echo "Total: $(du -sh . | cut -f1)"
echo ""
