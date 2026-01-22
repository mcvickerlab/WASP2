#!/bin/bash
###############################################################################
# Test Figure 2 Scripts
#
# Tests all Python scripts for syntax errors and import issues
# without running the full benchmarks.
###############################################################################

# NOTE: Avoid `set -u` here because conda's shell hooks are not nounset-safe.
set -eo pipefail

REPO_ROOT="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp"
cd "$REPO_ROOT"

echo "=================================================="
echo "Testing Figure 2 Scripts"
echo "=================================================="
echo ""

# Activate conda environment
source /iblm/netapp/home/jjaureguy/mambaforge/etc/profile.d/conda.sh
conda activate WASP2_dev2

ERRORS=0

# Test each Python script
SCRIPTS=(
    "paper/figure2/scripts/generate_count_comparison.py"
    "paper/figure2/scripts/generate_before_after_counts.py"
    "paper/figure2/scripts/generate_delta_counts.py"
    "paper/figure2/scripts/generate_figure2.py"
)

for script in "${SCRIPTS[@]}"; do
    echo "Testing: $script"

    # Check syntax
    if python -m py_compile "$script" 2>/dev/null; then
        echo "  ✓ Syntax OK"
    else
        echo "  ✗ Syntax error"
        python -m py_compile "$script"
        ERRORS=$((ERRORS + 1))
        continue
    fi

    # Check imports
    if python -c "import sys; sys.path.insert(0, 'paper'); exec(open('$script').read().split('if __name__')[0])" 2>/dev/null; then
        echo "  ✓ Imports OK"
    else
        echo "  ✗ Import error"
        python -c "import sys; sys.path.insert(0, 'paper'); exec(open('$script').read().split('if __name__')[0])"
        ERRORS=$((ERRORS + 1))
    fi

    # Check help (if argparse is used)
    if grep -q "argparse" "$script"; then
        if python "$script" --help > /dev/null 2>&1; then
            echo "  ✓ CLI OK"
        else
            echo "  ✗ CLI error"
            ERRORS=$((ERRORS + 1))
        fi
    fi

    echo ""
done

# Test figure generation with dummy data
echo "Testing figure generation with placeholders..."
if python paper/figure2/scripts/generate_figure2.py --dataset hg00731 2>&1 | grep -q "Saved:"; then
    echo "  ✓ Figure generation works (with placeholders)"
else
    echo "  ✗ Figure generation failed"
    ERRORS=$((ERRORS + 1))
fi
echo ""

# Summary
echo "=================================================="
if [ $ERRORS -eq 0 ]; then
    echo "✓ All tests passed!"
    echo ""
    echo "Scripts are ready to use. Next steps:"
    echo "  1. Run verify_setup.sh to check data files"
    echo "  2. Submit benchmarks: qsub paper/figure2/scripts/run_figure2_benchmarks.sh"
    exit 0
else
    echo "✗ Found $ERRORS error(s)"
    exit 1
fi
echo "=================================================="
