#!/bin/bash
###############################################################################
# Figure 2 Quick Start
#
# Complete pipeline to generate Figure 2 from scratch.
#
# Usage:
#   ./quickstart.sh                    # Run everything
#   ./quickstart.sh --skip-benchmarks  # Skip benchmarks (use existing data)
#   ./quickstart.sh --dataset hg00731  # Process only HG00731
###############################################################################

# NOTE: Avoid `set -u` here because conda's shell hooks are not nounset-safe.
set -eo pipefail

REPO_ROOT="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp"
cd "$REPO_ROOT/paper/figure2"

# Ensure the in-repo Python package (`src/wasp2`) is importable inside conda.
export PYTHONPATH="${REPO_ROOT}/src:${PYTHONPATH:-}"

# Parse arguments
SKIP_BENCHMARKS=false
# Main-text Figure 2 is RNA-seq (HG00731) by default.
DATASET="hg00731"

while [[ $# -gt 0 ]]; do
    case $1 in
        --skip-benchmarks)
            SKIP_BENCHMARKS=true
            shift
            ;;
        --dataset)
            DATASET="$2"
            shift 2
            ;;
        --help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --skip-benchmarks      Skip running benchmarks (use existing data)"
            echo "  --dataset DATASET      Process specific dataset: hg00731, gm12878, or both (default: hg00731)"
            echo "  --help                 Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

echo "=================================================="
echo "Figure 2 Quick Start"
echo "=================================================="
echo "Dataset: $DATASET"
echo "Skip benchmarks: $SKIP_BENCHMARKS"
echo ""

# Step 0: Verify setup
echo "Step 0: Verifying setup..."
if ! ./scripts/verify_setup.sh; then
    echo "ERROR: Setup verification failed. Please fix errors above."
    exit 1
fi
echo ""

# Step 1: Run benchmarks (if not skipped)
if [ "$SKIP_BENCHMARKS" = false ]; then
    echo "=================================================="
    echo "Step 1: Submitting benchmark job to SGE..."
    echo "=================================================="

    mkdir -p logs

    JOB_ID=$(qsub -v DATASET="$DATASET" scripts/run_figure2_benchmarks.sh | grep -oP '\d+')
    echo "Job submitted: $JOB_ID"
    LOG_PATH="logs/figure2_benchmarks.${JOB_ID}.log"
    echo "Log: ${LOG_PATH}"
    echo ""
    echo "Waiting for job to complete..."

    # Wait for job to complete
    while qstat -j $JOB_ID &> /dev/null; do
        sleep 30
        echo -n "."
    done
    echo ""
    echo "Job completed!"
    echo ""

    # Check if job succeeded (SGE stdout is directed to the per-job log file)
    if [ ! -f "${LOG_PATH}" ] || ! grep -q "Figure 2 Benchmarks Complete" "${LOG_PATH}"; then
        echo "ERROR: Benchmark job may have failed. Check log file:"
        echo "  tail -50 ${LOG_PATH}"
        exit 1
    fi
else
    echo "=================================================="
    echo "Step 1: Skipping benchmarks (using existing data)"
    echo "=================================================="
    echo ""
fi

# Step 2: Generate count comparisons
echo "=================================================="
echo "Step 2: Generating count comparisons..."
echo "=================================================="

source /iblm/netapp/home/jjaureguy/mambaforge/etc/profile.d/conda.sh
conda activate WASP2_dev2

python scripts/generate_count_comparison.py --dataset "$DATASET"

if [ $? -ne 0 ]; then
    echo "ERROR: Count comparison failed"
    exit 1
fi
echo ""

# Step 3: Generate before/after table (Panel C)
echo "=================================================="
echo "Step 3: Generating before/after table..."
echo "=================================================="

if [ "$DATASET" = "hg00731" ]; then
    python scripts/generate_before_after_counts.py --dataset "$DATASET" --min-total 10
    if [ $? -ne 0 ]; then
        echo "WARNING: Before/after table generation failed (continuing anyway)"
    fi
    python scripts/generate_delta_counts.py --dataset "$DATASET" --min-total 10
    if [ $? -ne 0 ]; then
        echo "WARNING: Delta-counts table generation failed (continuing anyway)"
    fi
else
    echo "Skipping before/after table for dataset '$DATASET' (currently supported for hg00731 only)."
fi
echo ""

# Step 4: Generate plots
echo "=================================================="
echo "Step 4: Generating plots..."
echo "=================================================="

python scripts/generate_figure2.py --dataset "$DATASET"

if [ $? -ne 0 ]; then
    echo "ERROR: Plot generation failed"
    exit 1
fi
echo ""

# Summary
echo "=================================================="
echo "✓ Figure 2 generation complete!"
echo "=================================================="
echo ""
echo "Output files:"
echo "  Data:  paper/figure2/data/"
echo "  Plots: paper/figure2/plots/"
echo ""

# List generated plots
echo "Generated plots:"
ls -lh plots/*.png 2>/dev/null || echo "  (none found)"
echo ""

# Show sample plot
if command -v display &> /dev/null; then
    echo "Open plot with: display plots/figure2_combined.png"
else
    echo "View plots at: $REPO_ROOT/paper/figure2/plots/"
fi

echo "=================================================="
