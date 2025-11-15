#!/bin/bash
# WASP2 Baseline Pipeline Execution
# Purpose: Run full counting→analysis pipeline and save baseline outputs

set -e  # Exit on error
set -o pipefail

echo "======================================"
echo " WASP2 Baseline Pipeline Execution"
echo "======================================"
echo "Started: $(date)"
echo ""

# Check dependencies
echo "Checking dependencies..."
DEPS_OK=true

for cmd in bcftools bedtools samtools python; do
    if ! command -v $cmd >/dev/null 2>&1; then
        echo "✗ $cmd not found!"
        DEPS_OK=false
    else
        echo "✓ $cmd found: $(command -v $cmd)"
    fi
done

if [ "$DEPS_OK" = false ]; then
    echo ""
    echo "ERROR: Missing required dependencies!"
    echo "Please install via: conda env create -f environment.yml"
    exit 1
fi

# Check Python imports
echo ""
echo "Checking Python dependencies..."
python -c "
import sys
try:
    import pysam, polars, pandas, scipy, typer, anndata
    print('✓ All Python packages imported successfully')
except ImportError as e:
    print(f'✗ Import error: {e}')
    sys.exit(1)
"

# Create output directories
echo ""
echo "Creating output directories..."
mkdir -p baselines/{mapping,counting,analysis}
mkdir -p baselines/counting/temp
mkdir -p baselines/analysis/temp

# Check test data exists
echo ""
echo "Verifying test data..."
if [ ! -f "test_data/CD4_ATACseq_Day1_merged_filtered.sort.bam" ]; then
    echo "✗ Test BAM file not found!"
    exit 1
fi
if [ ! -f "test_data/filter_chr10.vcf" ]; then
    echo "✗ Test VCF file not found!"
    exit 1
fi
echo "✓ Test data found"

# STEP 1: Counting
echo ""
echo "======================================"
echo " Step 1: Counting Alleles"
echo "======================================"
START_COUNT=$(date +%s)

python -m src.counting count-variants \
    test_data/CD4_ATACseq_Day1_merged_filtered.sort.bam \
    test_data/filter_chr10.vcf \
    --samples NA12878 \
    --region test_data/NA12878_snps_chr10.bed \
    --out baselines/counting/counts.tsv \
    --temp baselines/counting/temp

END_COUNT=$(date +%s)
COUNT_TIME=$((END_COUNT - START_COUNT))

if [ -f "baselines/counting/counts.tsv" ]; then
    COUNT_ROWS=$(wc -l < baselines/counting/counts.tsv)
    COUNT_MD5=$(md5sum baselines/counting/counts.tsv | awk '{print $1}')
    echo "✓ Counting complete: $COUNT_ROWS rows, MD5: $COUNT_MD5"
    echo "  Time: ${COUNT_TIME}s"
else
    echo "✗ Counting failed - output file not created!"
    exit 1
fi

# STEP 2: Analysis
echo ""
echo "======================================"
echo " Step 2: Analyzing Allelic Imbalance"
echo "======================================"
START_ANALYSIS=$(date +%s)

python -m src.analysis find-imbalance \
    baselines/counting/counts.tsv \
    --out baselines/analysis/ai_results.tsv \
    --min 10 \
    --pseudocount 1

END_ANALYSIS=$(date +%s)
ANALYSIS_TIME=$((END_ANALYSIS - START_ANALYSIS))

if [ -f "baselines/analysis/ai_results.tsv" ]; then
    ANALYSIS_ROWS=$(wc -l < baselines/analysis/ai_results.tsv)
    ANALYSIS_MD5=$(md5sum baselines/analysis/ai_results.tsv | awk '{print $1}')
    echo "✓ Analysis complete: $ANALYSIS_ROWS rows, MD5: $ANALYSIS_MD5"
    echo "  Time: ${ANALYSIS_TIME}s"
else
    echo "✗ Analysis failed - output file not created!"
    exit 1
fi

# Save baseline metadata
echo ""
echo "======================================"
echo " Saving Baseline Metadata"
echo "======================================"

cat > baselines/baseline_metadata.txt <<EOF
WASP2 Baseline Execution
========================
Date: $(date)
Host: $(hostname)
User: $(whoami)

Dependencies:
-------------
bcftools: $(bcftools --version | head -1)
bedtools: $(bedtools --version)
samtools: $(samtools --version | head -1)
python: $(python --version)

Test Data:
----------
BAM: test_data/CD4_ATACseq_Day1_merged_filtered.sort.bam ($(stat -c%s test_data/CD4_ATACseq_Day1_merged_filtered.sort.bam) bytes)
VCF: test_data/filter_chr10.vcf ($(stat -c%s test_data/filter_chr10.vcf) bytes)
BED: test_data/NA12878_snps_chr10.bed ($(stat -c%s test_data/NA12878_snps_chr10.bed) bytes)

Outputs:
--------
Counting: baselines/counting/counts.tsv
  Rows: $COUNT_ROWS
  MD5: $COUNT_MD5
  Time: ${COUNT_TIME}s

Analysis: baselines/analysis/ai_results.tsv
  Rows: $ANALYSIS_ROWS
  MD5: $ANALYSIS_MD5
  Time: ${ANALYSIS_TIME}s

Total Time: $((COUNT_TIME + ANALYSIS_TIME))s
EOF

echo "✓ Metadata saved to baselines/baseline_metadata.txt"

# Save checksums for regression testing
md5sum baselines/counting/counts.tsv > baselines/counting_baseline.md5
md5sum baselines/analysis/ai_results.tsv > baselines/analysis_baseline.md5

# Save first 10 lines of each output for manual inspection
head -10 baselines/counting/counts.tsv > baselines/counting/counts_head10.txt
head -10 baselines/analysis/ai_results.tsv > baselines/analysis/ai_results_head10.txt

echo ""
echo "======================================"
echo " Baseline Execution Complete!"
echo "======================================"
echo "Finished: $(date)"
echo ""
echo "Results Summary:"
echo "  Counting: $COUNT_ROWS rows (${COUNT_TIME}s)"
echo "  Analysis: $ANALYSIS_ROWS rows (${ANALYSIS_TIME}s)"
echo "  Total time: $((COUNT_TIME + ANALYSIS_TIME))s"
echo ""
echo "Baseline files saved in: baselines/"
echo "To validate future runs: ./scripts/validate_against_baseline.sh"
