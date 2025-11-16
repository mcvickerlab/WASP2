#!/bin/bash
# WASP2 Full Pipeline Baseline Execution
# Purpose: Run complete mapping→counting→analysis pipeline with memory profiling
# This creates comprehensive baseline outputs for regression testing

set -e  # Exit on error
set -o pipefail

echo "=========================================="
echo " WASP2 Full Pipeline Baseline"
echo "=========================================="
echo "Started: $(date)"
echo ""

# Configuration
SAMPLE="NA12878"
TEST_BAM="test_data/CD4_ATACseq_Day1_merged_filtered.sort.bam"
TEST_VCF="test_data/filter_chr10.vcf"
TEST_BED="test_data/NA12878_snps_chr10.bed"
GENOME_REF="test_data/genome.fa"  # Update if different location

# Output directories
BASELINE_DIR="baselines"
MAP_DIR="$BASELINE_DIR/mapping"
COUNT_DIR="$BASELINE_DIR/counting"
ANALYSIS_DIR="$BASELINE_DIR/analysis"

# Create output directories
echo "Creating output directories..."
mkdir -p $MAP_DIR $COUNT_DIR $ANALYSIS_DIR
mkdir -p $COUNT_DIR/temp $ANALYSIS_DIR/temp
echo "✓ Directories created"

# Check dependencies
echo ""
echo "Checking dependencies..."
DEPS_OK=true

for cmd in bcftools bedtools samtools python bwa; do
    if ! command -v $cmd >/dev/null 2>&1; then
        if [ "$cmd" = "bwa" ]; then
            echo "⚠ bwa not found - will skip mapping step"
            SKIP_MAPPING=true
        else
            echo "✗ $cmd not found!"
            DEPS_OK=false
        fi
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
    import pysam, polars, pandas, scipy, typer, anndata, pybedtools
    print('✓ All Python packages imported successfully')
except ImportError as e:
    print(f'✗ Import error: {e}')
    sys.exit(1)
"

# Check test data exists
echo ""
echo "Verifying test data..."
if [ ! -f "$TEST_BAM" ]; then
    echo "✗ Test BAM file not found: $TEST_BAM"
    exit 1
fi
if [ ! -f "$TEST_VCF" ]; then
    echo "✗ Test VCF file not found: $TEST_VCF"
    exit 1
fi
if [ ! -f "$TEST_BED" ]; then
    echo "✗ Test BED file not found: $TEST_BED"
    exit 1
fi
echo "✓ Test data found"

# Check for genome reference (optional for mapping)
if [ ! -f "$GENOME_REF" ]; then
    echo "⚠ Genome reference not found: $GENOME_REF"
    echo "  Mapping step will be skipped"
    SKIP_MAPPING=true
fi

# ===========================================
# STEP 1: WASP Mapping (if genome available)
# ===========================================

if [ "$SKIP_MAPPING" = true ]; then
    echo ""
    echo "=========================================="
    echo " Step 1: WASP Mapping (SKIPPED)"
    echo "=========================================="
    echo "⚠ Using original BAM for counting (no WASP filtering)"
    echo ""
    WASP_BAM="$TEST_BAM"  # Use original
else
    echo ""
    echo "=========================================="
    echo " Step 1: WASP Mapping"
    echo "=========================================="
    START_MAP=$(date +%s)

    # Step 1a: Generate swapped allele reads
    echo "1a. Generating swapped allele reads..."
    /usr/bin/time -v -o $MAP_DIR/make_reads_memory.txt \
        python src/mapping/__main__.py make-reads \
        $TEST_BAM \
        $TEST_VCF \
        --samples $SAMPLE \
        --out_dir $MAP_DIR \
        --out_json $MAP_DIR/wasp_data.json

    MAKE_READS_MEM=$(grep "Maximum resident set size" $MAP_DIR/make_reads_memory.txt | awk '{print $6}')
    MAKE_READS_MEM_MB=$((MAKE_READS_MEM / 1024))
    echo "✓ Swapped reads generated (Peak Memory: ${MAKE_READS_MEM_MB} MB)"

    # Check output files
    if [ ! -f "$MAP_DIR/CD4_ATACseq_Day1_merged_filtered.sort_swapped_alleles_r1.fq" ]; then
        echo "✗ Swapped reads not generated!"
        exit 1
    fi

    # Step 1b: Remap with BWA
    echo "1b. Remapping with BWA..."
    START_REMAP=$(date +%s)

    /usr/bin/time -v -o $MAP_DIR/bwa_memory.txt \
        bwa mem -t 4 \
        $GENOME_REF \
        $MAP_DIR/CD4_ATACseq_Day1_merged_filtered.sort_swapped_alleles_r1.fq \
        $MAP_DIR/CD4_ATACseq_Day1_merged_filtered.sort_swapped_alleles_r2.fq \
        | samtools view -b -o $MAP_DIR/remapped.unsorted.bam -

    # Sort and index remapped BAM
    samtools sort -o $MAP_DIR/remapped.bam $MAP_DIR/remapped.unsorted.bam
    samtools index $MAP_DIR/remapped.bam
    rm $MAP_DIR/remapped.unsorted.bam

    END_REMAP=$(date +%s)
    REMAP_TIME=$((END_REMAP - START_REMAP))
    BWA_MEM=$(grep "Maximum resident set size" $MAP_DIR/bwa_memory.txt | awk '{print $6}')
    BWA_MEM_MB=$((BWA_MEM / 1024))
    echo "✓ Remapping complete (Time: ${REMAP_TIME}s, Memory: ${BWA_MEM_MB} MB)"

    # Step 1c: Filter remapped reads
    echo "1c. Filtering remapped reads..."
    /usr/bin/time -v -o $MAP_DIR/filter_memory.txt \
        python src/mapping/__main__.py filter-remapped \
        $MAP_DIR/remapped.bam \
        --json $MAP_DIR/wasp_data.json \
        --out_bam $MAP_DIR/wasp_filtered.bam

    FILTER_MEM=$(grep "Maximum resident set size" $MAP_DIR/filter_memory.txt | awk '{print $6}')
    FILTER_MEM_MB=$((FILTER_MEM / 1024))

    END_MAP=$(date +%s)
    MAP_TIME=$((END_MAP - START_MAP))

    # Calculate total mapping memory (max of all steps)
    MAP_MEM=$MAKE_READS_MEM
    [ $BWA_MEM -gt $MAP_MEM ] && MAP_MEM=$BWA_MEM
    [ $FILTER_MEM -gt $MAP_MEM ] && MAP_MEM=$FILTER_MEM
    MAP_MEM_MB=$((MAP_MEM / 1024))

    echo "✓ WASP mapping complete"
    echo "  Total Time: ${MAP_TIME}s"
    echo "  Peak Memory: ${MAP_MEM_MB} MB (max across all steps)"
    echo ""

    # Use WASP-filtered BAM for counting
    WASP_BAM="$MAP_DIR/wasp_filtered.bam"

    # Get read counts
    ORIGINAL_READS=$(samtools view -c $TEST_BAM)
    WASP_READS=$(samtools view -c $WASP_BAM)
    KEPT_PCT=$((100 * WASP_READS / ORIGINAL_READS))
    echo "  Original reads: $ORIGINAL_READS"
    echo "  WASP filtered reads: $WASP_READS ($KEPT_PCT%)"
fi

# ===========================================
# STEP 2: Counting Alleles
# ===========================================

echo ""
echo "=========================================="
echo " Step 2: Counting Alleles"
echo "=========================================="
START_COUNT=$(date +%s)

/usr/bin/time -v -o $COUNT_DIR/memory_profile.txt \
    python src/counting/__main__.py count-variants \
    $WASP_BAM \
    $TEST_VCF \
    --samples $SAMPLE \
    --region $TEST_BED \
    --out $COUNT_DIR/counts.tsv \
    --temp $COUNT_DIR/temp

END_COUNT=$(date +%s)
COUNT_TIME=$((END_COUNT - START_COUNT))

COUNT_MEM=$(grep "Maximum resident set size" $COUNT_DIR/memory_profile.txt | awk '{print $6}')
COUNT_MEM_MB=$((COUNT_MEM / 1024))

if [ -f "$COUNT_DIR/counts.tsv" ]; then
    COUNT_ROWS=$(wc -l < $COUNT_DIR/counts.tsv)
    COUNT_MD5=$(md5sum $COUNT_DIR/counts.tsv | awk '{print $1}')

    # Count statistics
    COUNT_SNPS=$(tail -n +2 $COUNT_DIR/counts.tsv | wc -l)
    TOTAL_REF=$(tail -n +2 $COUNT_DIR/counts.tsv | awk -F'\t' '{sum+=$(NF-2)} END {print sum}')
    TOTAL_ALT=$(tail -n +2 $COUNT_DIR/counts.tsv | awk -F'\t' '{sum+=$(NF-1)} END {print sum}')
    TOTAL_OTHER=$(tail -n +2 $COUNT_DIR/counts.tsv | awk -F'\t' '{sum+=$NF} END {print sum}')
    TOTAL_COUNTS=$((TOTAL_REF + TOTAL_ALT + TOTAL_OTHER))

    echo "✓ Counting complete"
    echo "  Rows: $COUNT_ROWS (header + $COUNT_SNPS SNPs)"
    echo "  Total allele counts: $TOTAL_COUNTS (ref: $TOTAL_REF, alt: $TOTAL_ALT, other: $TOTAL_OTHER)"
    echo "  MD5: $COUNT_MD5"
    echo "  Time: ${COUNT_TIME}s, Peak Memory: ${COUNT_MEM_MB} MB"
else
    echo "✗ Counting failed - output file not created!"
    exit 1
fi

# ===========================================
# STEP 3: Analyzing Allelic Imbalance
# ===========================================

echo ""
echo "=========================================="
echo " Step 3: Analyzing Allelic Imbalance"
echo "=========================================="
START_ANALYSIS=$(date +%s)

/usr/bin/time -v -o $ANALYSIS_DIR/memory_profile.txt \
    python src/analysis/__main__.py find-imbalance \
    $COUNT_DIR/counts.tsv \
    --out $ANALYSIS_DIR/ai_results.tsv \
    --min 10 \
    --pseudocount 1

END_ANALYSIS=$(date +%s)
ANALYSIS_TIME=$((END_ANALYSIS - START_ANALYSIS))

ANALYSIS_MEM=$(grep "Maximum resident set size" $ANALYSIS_DIR/memory_profile.txt | awk '{print $6}')
ANALYSIS_MEM_MB=$((ANALYSIS_MEM / 1024))

if [ -f "$ANALYSIS_DIR/ai_results.tsv" ]; then
    ANALYSIS_ROWS=$(wc -l < $ANALYSIS_DIR/ai_results.tsv)
    ANALYSIS_MD5=$(md5sum $ANALYSIS_DIR/ai_results.tsv | awk '{print $1}')

    # Analysis statistics
    AI_REGIONS=$(tail -n +2 $ANALYSIS_DIR/ai_results.tsv | wc -l)
    SIG_REGIONS=$(tail -n +2 $ANALYSIS_DIR/ai_results.tsv | awk -F'\t' '$NF < 0.05 {count++} END {print count+0}')

    echo "✓ Analysis complete"
    echo "  Rows: $ANALYSIS_ROWS (header + $AI_REGIONS regions)"
    echo "  Significant regions (FDR < 0.05): $SIG_REGIONS"
    echo "  MD5: $ANALYSIS_MD5"
    echo "  Time: ${ANALYSIS_TIME}s, Peak Memory: ${ANALYSIS_MEM_MB} MB"
else
    echo "✗ Analysis failed - output file not created!"
    exit 1
fi

# ===========================================
# Save Comprehensive Metadata
# ===========================================

echo ""
echo "=========================================="
echo " Saving Pipeline Metadata"
echo "=========================================="

# Calculate totals
if [ "$SKIP_MAPPING" = true ]; then
    TOTAL_TIME=$((COUNT_TIME + ANALYSIS_TIME))
    PEAK_MEM=$((COUNT_MEM > ANALYSIS_MEM ? COUNT_MEM : ANALYSIS_MEM))
else
    TOTAL_TIME=$((MAP_TIME + COUNT_TIME + ANALYSIS_TIME))
    PEAK_MEM=$MAP_MEM
    [ $COUNT_MEM -gt $PEAK_MEM ] && PEAK_MEM=$COUNT_MEM
    [ $ANALYSIS_MEM -gt $PEAK_MEM ] && PEAK_MEM=$ANALYSIS_MEM
fi
PEAK_MEM_MB=$((PEAK_MEM / 1024))

cat > $BASELINE_DIR/pipeline_metadata.txt <<EOF
WASP2 Full Pipeline Baseline
=============================
Date: $(date)
Host: $(hostname)
User: $(whoami)

Dependencies:
-------------
bcftools: $(bcftools --version | head -1)
bedtools: $(bedtools --version)
samtools: $(samtools --version | head -1)
python: $(python --version)
EOF

if [ "$SKIP_MAPPING" != true ]; then
cat >> $BASELINE_DIR/pipeline_metadata.txt <<EOF
bwa: $(bwa 2>&1 | grep "Version" | head -1)
EOF
fi

cat >> $BASELINE_DIR/pipeline_metadata.txt <<EOF

Input Data:
-----------
BAM: $TEST_BAM ($(stat -c%s $TEST_BAM) bytes)
VCF: $TEST_VCF ($(stat -c%s $TEST_VCF) bytes)
BED: $TEST_BED ($(stat -c%s $TEST_BED) bytes)
Sample: $SAMPLE
EOF

if [ "$SKIP_MAPPING" != true ]; then
cat >> $BASELINE_DIR/pipeline_metadata.txt <<EOF
Genome: $GENOME_REF

Step 1: WASP Mapping
--------------------
Original reads: $ORIGINAL_READS
WASP filtered reads: $WASP_READS ($KEPT_PCT% kept)
Output: $MAP_DIR/wasp_filtered.bam
Time: ${MAP_TIME}s
Peak Memory: ${MAP_MEM_MB} MB (${MAP_MEM} KB)

  1a. Make reads - Memory: ${MAKE_READS_MEM_MB} MB
  1b. BWA remap - Memory: ${BWA_MEM_MB} MB, Time: ${REMAP_TIME}s
  1c. Filter - Memory: ${FILTER_MEM_MB} MB
EOF
else
cat >> $BASELINE_DIR/pipeline_metadata.txt <<EOF

Step 1: WASP Mapping
--------------------
Status: SKIPPED (no genome reference)
Using original BAM: $TEST_BAM
EOF
fi

cat >> $BASELINE_DIR/pipeline_metadata.txt <<EOF

Step 2: Counting
----------------
Output: $COUNT_DIR/counts.tsv
Rows: $COUNT_ROWS (header + $COUNT_SNPS SNPs)
Total allele counts: $TOTAL_COUNTS
  Ref alleles: $TOTAL_REF
  Alt alleles: $TOTAL_ALT
  Other alleles: $TOTAL_OTHER
MD5: $COUNT_MD5
Time: ${COUNT_TIME}s
Peak Memory: ${COUNT_MEM_MB} MB (${COUNT_MEM} KB)

Step 3: Analysis
----------------
Output: $ANALYSIS_DIR/ai_results.tsv
Rows: $ANALYSIS_ROWS (header + $AI_REGIONS regions tested)
Significant regions (FDR < 0.05): $SIG_REGIONS
MD5: $ANALYSIS_MD5
Time: ${ANALYSIS_TIME}s
Peak Memory: ${ANALYSIS_MEM_MB} MB (${ANALYSIS_MEM} KB)

Performance Summary:
--------------------
Total Pipeline Time: ${TOTAL_TIME}s
Peak Memory Usage: ${PEAK_MEM_MB} MB
Memory Efficiency: $(awk "BEGIN {printf \"%.2f\", $PEAK_MEM_MB / $TOTAL_TIME}") MB/s
EOF

echo "✓ Metadata saved to $BASELINE_DIR/pipeline_metadata.txt"

# Save checksums for regression testing
if [ "$SKIP_MAPPING" != true ]; then
    md5sum $MAP_DIR/wasp_filtered.bam > $BASELINE_DIR/mapping_baseline.md5
fi
md5sum $COUNT_DIR/counts.tsv > $BASELINE_DIR/counting_baseline.md5
md5sum $ANALYSIS_DIR/ai_results.tsv > $BASELINE_DIR/analysis_baseline.md5

# Save sample outputs for manual inspection
head -20 $COUNT_DIR/counts.tsv > $COUNT_DIR/counts_head20.txt
head -20 $ANALYSIS_DIR/ai_results.tsv > $ANALYSIS_DIR/ai_results_head20.txt

echo ""
echo "=========================================="
echo " Pipeline Execution Complete!"
echo "=========================================="
echo "Finished: $(date)"
echo ""
echo "Pipeline Summary:"
if [ "$SKIP_MAPPING" != true ]; then
    echo "  1. Mapping: $ORIGINAL_READS → $WASP_READS reads (${MAP_TIME}s, ${MAP_MEM_MB} MB)"
fi
echo "  2. Counting: $COUNT_SNPS SNPs, $TOTAL_COUNTS alleles (${COUNT_TIME}s, ${COUNT_MEM_MB} MB)"
echo "  3. Analysis: $AI_REGIONS regions, $SIG_REGIONS significant (${ANALYSIS_TIME}s, ${ANALYSIS_MEM_MB} MB)"
echo ""
echo "  Total time: ${TOTAL_TIME}s"
echo "  Peak memory: ${PEAK_MEM_MB} MB"
echo ""
echo "Baseline files saved in: baselines/"
echo "Memory profiles: baselines/{mapping,counting,analysis}/*memory*.txt"
echo "Metadata: baselines/pipeline_metadata.txt"
