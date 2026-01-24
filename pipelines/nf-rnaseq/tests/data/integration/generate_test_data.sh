#!/bin/bash
# =============================================================================
# WASP2 RNA-seq Integration Test Data Generation Script
# =============================================================================
# This script generates all test data files needed for integration testing:
# - FASTA index (.fai)
# - Compressed and indexed VCF (.vcf.gz, .vcf.gz.tbi)
# - Simulated paired-end FASTQ reads (sample1_R1.fq.gz, sample1_R2.fq.gz)
# - STAR genome index
#
# Prerequisites:
#   - samtools (for faidx)
#   - bcftools or bgzip/tabix (for VCF compression/indexing)
#   - wgsim (for read simulation, part of samtools)
#   - STAR (for genome indexing)
#
# Usage:
#   cd pipelines/nf-rnaseq/tests/data/integration
#   bash generate_test_data.sh
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "==================================================================="
echo " WASP2 Integration Test Data Generator"
echo "==================================================================="
echo "Working directory: $SCRIPT_DIR"
echo ""

# -----------------------------------------------------------------------------
# Check prerequisites
# -----------------------------------------------------------------------------
echo "[1/6] Checking prerequisites..."

check_tool() {
    if ! command -v "$1" &> /dev/null; then
        echo "ERROR: $1 is required but not found in PATH"
        exit 1
    fi
    echo "  ✓ $1 found"
}

check_tool samtools
check_tool bgzip
check_tool tabix
check_tool wgsim
check_tool STAR

echo ""

# -----------------------------------------------------------------------------
# Index FASTA reference
# -----------------------------------------------------------------------------
echo "[2/6] Indexing FASTA reference..."

if [[ -f "chr_test.fa.fai" ]]; then
    echo "  chr_test.fa.fai already exists, skipping"
else
    samtools faidx chr_test.fa
    echo "  ✓ Created chr_test.fa.fai"
fi

echo ""

# -----------------------------------------------------------------------------
# Compress and index VCF
# -----------------------------------------------------------------------------
echo "[3/6] Compressing and indexing VCF..."

if [[ -f "integration.vcf.gz" && -f "integration.vcf.gz.tbi" ]]; then
    echo "  integration.vcf.gz and .tbi already exist, skipping"
else
    # Remove any existing files to ensure clean state
    rm -f integration.vcf.gz integration.vcf.gz.tbi

    # Compress with bgzip (required for tabix)
    bgzip -c integration.vcf > integration.vcf.gz
    echo "  ✓ Created integration.vcf.gz"

    # Index with tabix
    tabix -p vcf integration.vcf.gz
    echo "  ✓ Created integration.vcf.gz.tbi"
fi

echo ""

# -----------------------------------------------------------------------------
# Simulate paired-end reads
# -----------------------------------------------------------------------------
echo "[4/6] Simulating paired-end FASTQ reads..."

# Parameters for wgsim:
#   -N 500    : Number of read pairs
#   -1 100    : Length of read 1
#   -2 100    : Length of read 2
#   -r 0      : Rate of mutations (0 = perfect reads)
#   -R 0      : Fraction of indels (0 = no indels)
#   -X 0      : Probability of extension
#   -e 0.001  : Base error rate
#   -S 42     : Random seed for reproducibility
#   -d 300    : Outer distance (fragment size)
#   -s 50     : Standard deviation of fragment size

NUM_READS=500
READ_LEN=100
FRAG_SIZE=300
FRAG_STD=50
ERROR_RATE=0.001
SEED=42

if [[ -f "sample1_R1.fq.gz" && -f "sample1_R2.fq.gz" ]]; then
    echo "  sample1_R1.fq.gz and sample1_R2.fq.gz already exist"
    echo "  To regenerate, delete these files and re-run"
else
    echo "  Generating ${NUM_READS} read pairs (${READ_LEN}bp, ${FRAG_SIZE}bp fragments)..."

    # Generate reads (wgsim outputs uncompressed FASTQ)
    wgsim -N $NUM_READS \
          -1 $READ_LEN \
          -2 $READ_LEN \
          -r 0 \
          -R 0 \
          -X 0 \
          -e $ERROR_RATE \
          -S $SEED \
          -d $FRAG_SIZE \
          -s $FRAG_STD \
          chr_test.fa \
          sample1_R1.fq \
          sample1_R2.fq \
          > /dev/null 2>&1

    # Compress the FASTQ files
    gzip -f sample1_R1.fq
    gzip -f sample1_R2.fq

    echo "  ✓ Created sample1_R1.fq.gz ($(du -h sample1_R1.fq.gz | cut -f1))"
    echo "  ✓ Created sample1_R2.fq.gz ($(du -h sample1_R2.fq.gz | cut -f1))"
fi

echo ""

# -----------------------------------------------------------------------------
# Build STAR index
# -----------------------------------------------------------------------------
echo "[5/6] Building STAR genome index..."

# Parameters for small genome:
#   --genomeSAindexNbases : Must be scaled for small genomes
#                           Formula: min(14, log2(genomeLength)/2 - 1)
#                           For ~20kb: min(14, log2(20000)/2 - 1) ≈ 6

STAR_INDEX_DIR="star_index"
GENOME_SIZE=$(grep -v "^>" chr_test.fa | tr -d '\n' | wc -c)

# Calculate appropriate --genomeSAindexNbases for small genome
# log2(20000) ≈ 14.3, so 14.3/2 - 1 ≈ 6
SA_INDEX_BASES=6

if [[ -f "${STAR_INDEX_DIR}/SAindex" ]]; then
    echo "  STAR index already exists in ${STAR_INDEX_DIR}/"
    echo "  To regenerate, delete the directory and re-run"
else
    echo "  Building index for ${GENOME_SIZE}bp genome (genomeSAindexNbases=${SA_INDEX_BASES})..."

    # Create fresh index directory
    rm -rf "${STAR_INDEX_DIR}"
    mkdir -p "${STAR_INDEX_DIR}"

    # Build STAR index
    STAR --runMode genomeGenerate \
         --genomeDir "${STAR_INDEX_DIR}" \
         --genomeFastaFiles chr_test.fa \
         --genomeSAindexNbases ${SA_INDEX_BASES} \
         --runThreadN 2 \
         --outFileNamePrefix "${STAR_INDEX_DIR}/" \
         2>&1 | tail -5

    echo "  ✓ Created STAR index ($(du -sh ${STAR_INDEX_DIR} | cut -f1))"
fi

echo ""

# -----------------------------------------------------------------------------
# Validate generated files
# -----------------------------------------------------------------------------
echo "[6/6] Validating generated files..."

validate_file() {
    if [[ -f "$1" ]]; then
        size=$(du -h "$1" | cut -f1)
        echo "  ✓ $1 ($size)"
    else
        echo "  ✗ $1 NOT FOUND"
        exit 1
    fi
}

validate_file "chr_test.fa"
validate_file "chr_test.fa.fai"
validate_file "integration.gtf"
validate_file "integration.vcf"
validate_file "integration.vcf.gz"
validate_file "integration.vcf.gz.tbi"
validate_file "sample1_R1.fq.gz"
validate_file "sample1_R2.fq.gz"
validate_file "samplesheet_integration.csv"
validate_file "star_index/SAindex"

echo ""
echo "==================================================================="
echo " SUCCESS! All integration test data generated."
echo "==================================================================="
echo ""
echo "Total disk usage: $(du -sh . | cut -f1)"
echo ""
echo "To run integration tests:"
echo "  cd pipelines/nf-rnaseq"
echo "  nf-test test tests/integration.nf.test --profile test_integration,conda"
echo ""
