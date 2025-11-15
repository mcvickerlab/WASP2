#!/bin/bash
# WASP2 Regression Test - Validate Against Baseline
# Purpose: Run pipeline and compare outputs against saved baseline

set -e
set -o pipefail

echo "======================================"
echo " WASP2 Regression Test"
echo "======================================"
echo "Started: $(date)"
echo ""

# Check if baseline exists
if [ ! -f "baselines/counting_baseline.md5" ]; then
    echo "✗ No baseline found!"
    echo "Please run: ./scripts/run_baseline.sh first"
    exit 1
fi

echo "Baseline found:"
cat baselines/baseline_metadata.txt | grep "Date:"
echo ""

# Create test output directory
TEST_DIR="baselines/test_run_$(date +%Y%m%d_%H%M%S)"
mkdir -p $TEST_DIR/{counting,analysis}

echo "Test outputs will be saved to: $TEST_DIR"
echo ""

# STEP 1: Run Counting
echo "======================================"
echo " Running Counting Module"
echo "======================================"

python -m src.counting count-variants \
    test_data/CD4_ATACseq_Day1_merged_filtered.sort.bam \
    test_data/filter_chr10.vcf \
    --samples NA12878 \
    --region test_data/NA12878_snps_chr10.bed \
    --out $TEST_DIR/counting/counts.tsv \
    --temp $TEST_DIR/counting/temp

# STEP 2: Run Analysis
echo ""
echo "======================================"
echo " Running Analysis Module"
echo "======================================"

python -m src.analysis find-imbalance \
    $TEST_DIR/counting/counts.tsv \
    --out $TEST_DIR/analysis/ai_results.tsv \
    --min 10 \
    --pseudocount 1

# STEP 3: Compare Outputs
echo ""
echo "======================================"
echo " Comparing Against Baseline"
echo "======================================"

# Counting comparison
BASELINE_COUNT_MD5=$(cat baselines/counting_baseline.md5 | awk '{print $1}')
TEST_COUNT_MD5=$(md5sum $TEST_DIR/counting/counts.tsv | awk '{print $1}')

echo "Counting Module:"
echo "  Baseline MD5: $BASELINE_COUNT_MD5"
echo "  Test MD5:     $TEST_COUNT_MD5"

if [ "$BASELINE_COUNT_MD5" = "$TEST_COUNT_MD5" ]; then
    echo "  ✓ PASS - Outputs are identical"
    COUNT_PASS=true
else
    echo "  ✗ FAIL - Outputs differ!"
    COUNT_PASS=false

    # Show differences
    echo ""
    echo "Detailed diff (first 20 lines):"
    diff baselines/counting/counts.tsv $TEST_DIR/counting/counts.tsv | head -20 || true
fi

echo ""

# Analysis comparison
BASELINE_ANALYSIS_MD5=$(cat baselines/analysis_baseline.md5 | awk '{print $1}')
TEST_ANALYSIS_MD5=$(md5sum $TEST_DIR/analysis/ai_results.tsv | awk '{print $1}')

echo "Analysis Module:"
echo "  Baseline MD5: $BASELINE_ANALYSIS_MD5"
echo "  Test MD5:     $TEST_ANALYSIS_MD5"

if [ "$BASELINE_ANALYSIS_MD5" = "$TEST_ANALYSIS_MD5" ]; then
    echo "  ✓ PASS - Outputs are identical"
    ANALYSIS_PASS=true
else
    echo "  ✗ FAIL - Outputs differ!"
    ANALYSIS_PASS=false

    # Show differences
    echo ""
    echo "Detailed diff (first 20 lines):"
    diff baselines/analysis/ai_results.tsv $TEST_DIR/analysis/ai_results.tsv | head -20 || true
fi

# Row count comparison
echo ""
echo "Row Count Comparison:"
BASELINE_COUNT_ROWS=$(wc -l < baselines/counting/counts.tsv)
TEST_COUNT_ROWS=$(wc -l < $TEST_DIR/counting/counts.tsv)
BASELINE_ANALYSIS_ROWS=$(wc -l < baselines/analysis/ai_results.tsv)
TEST_ANALYSIS_ROWS=$(wc -l < $TEST_DIR/analysis/ai_results.tsv)

echo "  Counting:  Baseline=$BASELINE_COUNT_ROWS, Test=$TEST_COUNT_ROWS"
echo "  Analysis:  Baseline=$BASELINE_ANALYSIS_ROWS, Test=$TEST_ANALYSIS_ROWS"

# Final result
echo ""
echo "======================================"
echo " Regression Test Results"
echo "======================================"

if [ "$COUNT_PASS" = true ] && [ "$ANALYSIS_PASS" = true ]; then
    echo "✓ ALL TESTS PASSED"
    echo ""
    echo "No regressions detected - all outputs match baseline!"
    exit 0
else
    echo "✗ TESTS FAILED"
    echo ""
    if [ "$COUNT_PASS" = false ]; then
        echo "  - Counting module output differs from baseline"
    fi
    if [ "$ANALYSIS_PASS" = false ]; then
        echo "  - Analysis module output differs from baseline"
    fi
    echo ""
    echo "This indicates either:"
    echo "  1. A regression bug was introduced"
    echo "  2. Intentional changes that require baseline update"
    echo ""
    echo "Review diff output above and determine if changes are expected."
    echo "If changes are intentional, update baseline: ./scripts/run_baseline.sh"
    exit 1
fi
