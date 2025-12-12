#!/bin/bash
# Master script to run WASP2-Rust vs Python DEV comparison
# Submits jobs in sequence with dependencies

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

echo "========================================"
echo "WASP2 Pipeline Comparison Runner"
echo "========================================"
echo "Script directory: ${SCRIPT_DIR}"
echo ""

# Step 1: Prepare chr22 subset
echo "Submitting Step 1: Prepare chr22 subset..."
JOB1=$(qsub -terse prepare_chr22_subset.sh)
echo "  Job ID: ${JOB1}"

# Step 2 & 3: Run both pipelines (after step 1 completes)
echo "Submitting Step 2: Rust pipeline (depends on ${JOB1})..."
JOB2=$(qsub -terse -hold_jid "${JOB1}" run_rust_chr22.sh)
echo "  Job ID: ${JOB2}"

echo "Submitting Step 3: Python pipeline (depends on ${JOB1})..."
JOB3=$(qsub -terse -hold_jid "${JOB1}" run_python_chr22.sh)
echo "  Job ID: ${JOB3}"

# Step 4: Compare outputs (after steps 2 & 3 complete)
echo "Submitting Step 4: Compare outputs (depends on ${JOB2},${JOB3})..."
JOB4=$(qsub -terse -hold_jid "${JOB2},${JOB3}" compare_chr22_outputs.sh)
echo "  Job ID: ${JOB4}"

echo ""
echo "========================================"
echo "All jobs submitted!"
echo "========================================"
echo ""
echo "Job chain:"
echo "  1. Prepare chr22: ${JOB1}"
echo "  2. Rust pipeline: ${JOB2} (waits for ${JOB1})"
echo "  3. Python pipeline: ${JOB3} (waits for ${JOB1})"
echo "  4. Compare outputs: ${JOB4} (waits for ${JOB2}, ${JOB3})"
echo ""
echo "Monitor with: qstat | grep -E '${JOB1}|${JOB2}|${JOB3}|${JOB4}'"
echo ""
echo "After completion, run detailed analysis:"
echo "  cd ../chr22_comparison"
echo "  python ../scripts/05_detailed_analysis.py"
echo ""
echo "Or use the Jupyter notebook:"
echo "  jupyter notebook analysis/chr22_pipeline_comparison.ipynb"
