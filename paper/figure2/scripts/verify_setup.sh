#!/bin/bash
###############################################################################
# Verify Figure 2 Setup
#
# Checks that all required files, tools, and dependencies are available
# before running benchmarks.
###############################################################################

# NOTE: Avoid `set -u` here because conda's shell hooks are not nounset-safe.
set -eo pipefail

REPO_ROOT="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp"
cd "$REPO_ROOT"

# Ensure the in-repo Python package (`src/wasp2`) is importable inside conda.
export PYTHONPATH="${REPO_ROOT}/src:${PYTHONPATH:-}"

echo "=================================================="
echo "Figure 2 Setup Verification"
echo "=================================================="
echo ""

# Track errors
ERRORS=0

# Which dataset(s) to verify. Main-text Figure 2 is RNA-seq (HG00731).
# ATAC-seq support is optional/supplementary.
DATASET="${DATASET:-hg00731}"  # hg00731 | gm12878 | both

# -----------------------------------------------------------------------------
# Check Input Files
# -----------------------------------------------------------------------------

echo "Checking input files..."
echo ""

# Reference genome
REF_GENOME="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/paper/figure2/tools/hg38.fa"
if [ -f "$REF_GENOME" ]; then
    echo "✓ Reference genome: $REF_GENOME"
else
    echo "✗ Missing reference genome: $REF_GENOME"
    ERRORS=$((ERRORS + 1))
fi

# HG00731 files
if [[ "$DATASET" == "hg00731" || "$DATASET" == "both" ]]; then
    echo ""
    echo "HG00731 RNA-seq files:"
    HG00731_BAM="$REPO_ROOT/benchmarking/star_wasp_comparison/results/wasp2rust_fair_2025-12-15_07-06-54/remap_keep.bam"
    HG00731_VCF="$REPO_ROOT/benchmarking/star_wasp_comparison/data/HG00731_het_only_chr.vcf.gz"

    if [ -f "$HG00731_BAM" ]; then
        echo "  ✓ BAM: $HG00731_BAM"
        SIZE=$(du -h "$HG00731_BAM" | cut -f1)
        echo "    Size: $SIZE"
    else
        echo "  ✗ Missing BAM: $HG00731_BAM"
        ERRORS=$((ERRORS + 1))
    fi

    if [ -f "${HG00731_BAM}.bai" ]; then
        echo "  ✓ BAM index: ${HG00731_BAM}.bai"
    else
        echo "  ✗ Missing BAM index: ${HG00731_BAM}.bai"
        ERRORS=$((ERRORS + 1))
    fi

    if [ -f "$HG00731_VCF" ]; then
        echo "  ✓ VCF: $HG00731_VCF"
        SIZE=$(du -h "$HG00731_VCF" | cut -f1)
        VARS=$(zcat "$HG00731_VCF" | grep -v "^#" | wc -l)
        echo "    Size: $SIZE"
        echo "    Variants: $VARS"
    else
        echo "  ✗ Missing VCF: $HG00731_VCF"
        ERRORS=$((ERRORS + 1))
    fi

    if [ -f "${HG00731_VCF}.tbi" ]; then
        echo "  ✓ VCF index: ${HG00731_VCF}.tbi"
    else
        echo "  ✗ Missing VCF index: ${HG00731_VCF}.tbi"
        ERRORS=$((ERRORS + 1))
    fi
fi

# GM12878 files
if [[ "$DATASET" == "gm12878" || "$DATASET" == "both" ]]; then
    echo ""
    echo "GM12878 ATAC-seq files:"
    GM12878_BAM_ORIG="/iblm/netapp/data1/aho/atac/Buenrostro2013/merged/GM12878_ATACseq_50k/GM12878_ATACseq_50k_merged.sorted.bam"
    GM12878_BAM_REMAP="$REPO_ROOT/benchmarking/atac_gm12878/results/wasp2rust_snp_fixed_2025-12-15_03-43-24/remap_keep.bam"
    GM12878_VCF="/iblm/netapp/data1/aho/variants/NA12878.vcf.gz"

    if [ -f "$GM12878_BAM_ORIG" ]; then
        echo "  ✓ Original BAM: $GM12878_BAM_ORIG"
        SIZE=$(du -h "$GM12878_BAM_ORIG" | cut -f1)
        echo "    Size: $SIZE"
    else
        echo "  ✗ Missing original BAM: $GM12878_BAM_ORIG"
        ERRORS=$((ERRORS + 1))
    fi

    if [ -f "$GM12878_BAM_REMAP" ]; then
        echo "  ✓ Remapped BAM: $GM12878_BAM_REMAP"
        SIZE=$(du -h "$GM12878_BAM_REMAP" | cut -f1)
        echo "    Size: $SIZE"
    else
        echo "  ✗ Missing remapped BAM: $GM12878_BAM_REMAP"
        ERRORS=$((ERRORS + 1))
    fi

    if [ -f "$GM12878_VCF" ]; then
        echo "  ✓ VCF: $GM12878_VCF"
        SIZE=$(du -h "$GM12878_VCF" | cut -f1)
        echo "    Size: $SIZE"
    else
        echo "  ✗ Missing VCF: $GM12878_VCF"
        ERRORS=$((ERRORS + 1))
    fi
fi

# -----------------------------------------------------------------------------
# Check Tools
# -----------------------------------------------------------------------------

echo ""
echo "Checking tools and environment..."
echo ""

# Activate conda environment
source /iblm/netapp/home/jjaureguy/mambaforge/etc/profile.d/conda.sh
conda activate WASP2_dev2 2>/dev/null || true

# Python
if command -v python &> /dev/null; then
    PYTHON_VERSION=$(python --version 2>&1)
    echo "✓ Python: $PYTHON_VERSION"
else
    echo "✗ Python not found"
    ERRORS=$((ERRORS + 1))
fi

# GATK
if command -v gatk &> /dev/null; then
    GATK_BIN="$(command -v gatk)"
    # The GATK wrapper in this environment does not support `--version` and can emit a lot of output.
    # Prefer a lightweight capability check instead.
    if gatk ASEReadCounter --help >/dev/null 2>&1; then
        echo "✓ GATK: $GATK_BIN (ASEReadCounter available)"
    else
        echo "✓ GATK: $GATK_BIN"
    fi
else
    echo "✗ GATK not found (required for benchmarks)"
    ERRORS=$((ERRORS + 1))
fi

# phASER
PHASER_PATH="$REPO_ROOT/benchmarking/phaser_tool/phaser/phaser.py"
if [ -f "$PHASER_PATH" ]; then
    echo "✓ phASER: $PHASER_PATH"
else
    echo "✗ phASER not found: $PHASER_PATH"
    ERRORS=$((ERRORS + 1))
fi

# WASP2 counting module
if python -c "from src.counting import run_counting" 2>/dev/null; then
    echo "✓ WASP2 counting module"
else
    echo "✗ WASP2 counting module not importable"
    ERRORS=$((ERRORS + 1))
fi

# -----------------------------------------------------------------------------
# Check Python Dependencies
# -----------------------------------------------------------------------------

echo ""
echo "Checking Python dependencies..."
echo ""

for pkg in numpy pandas scipy matplotlib pysam; do
    if python -c "import $pkg" 2>/dev/null; then
        VERSION=$(python -c "import $pkg; print($pkg.__version__)" 2>/dev/null || echo "unknown")
        echo "  ✓ $pkg ($VERSION)"
    else
        echo "  ✗ $pkg not found"
        ERRORS=$((ERRORS + 1))
    fi
done

# -----------------------------------------------------------------------------
# Check Directory Structure
# -----------------------------------------------------------------------------

echo ""
echo "Checking directory structure..."
echo ""

DIRS=(
    "paper/figure2/scripts"
    "paper/figure2/data"
    "paper/figure2/data/hg00731"
    "paper/figure2/data/gm12878"
    "paper/figure2/plots"
    "paper/figure2/logs"
)

for dir in "${DIRS[@]}"; do
    FULL_PATH="$REPO_ROOT/$dir"
    if [ -d "$FULL_PATH" ]; then
        echo "  ✓ $dir"
    else
        echo "  ✗ Missing directory: $dir"
        echo "    Creating: $FULL_PATH"
        mkdir -p "$FULL_PATH"
    fi
done

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------

echo ""
echo "=================================================="
if [ $ERRORS -eq 0 ]; then
    echo "✓ All checks passed! Ready to run benchmarks."
    echo ""
    echo "Next steps:"
    echo "  1. Submit benchmark job: qsub paper/figure2/scripts/run_figure2_benchmarks.sh"
    echo "  2. After completion, run: python paper/figure2/scripts/generate_count_comparison.py"
    echo "  3. Then run: python paper/figure2/scripts/generate_before_after_counts.py --dataset hg00731 --min-total 10"
    echo "  4. Finally: python paper/figure2/scripts/generate_figure2.py --dataset combined"
    exit 0
else
    echo "✗ Found $ERRORS error(s). Please fix before running benchmarks."
    exit 1
fi
echo "=================================================="
