#!/bin/bash
# =============================================================================
# WASP2 Shared Core Test Data Generator
# =============================================================================
# Generates the unified test dataset used by all WASP2 pipelines, Galaxy tools,
# CLI smoke tests, and container validation.
#
# Outputs (all committed to git, ~700K total):
#   chr_test.fa + .fai          - 20kb synthetic reference genome (2 gene regions)
#   variants.vcf + .gz + .tbi   - 10 het SNPs across 2 samples
#   annotation.gtf              - 2 genes, 6 exons
#   regions.bed                 - Peak/region file from exon coordinates
#   sample{1,2,3}.bam + .bai   - Aligned reads (wgsim + bwa)
#   bwa_index/                  - BWA index for chr_test.fa
#   expected_counts.tsv         - WASP2 counting output baseline
#   expected_analysis.tsv       - WASP2 analysis output baseline (placeholder)
#
# Prerequisites: samtools, bgzip, tabix, wgsim, bwa, bcftools
# Conda env:     conda activate WASP2_dev2
#
# Usage:
#   cd tests/shared_data
#   bash generate_core_data.sh
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "==================================================================="
echo " WASP2 Shared Core Test Data Generator"
echo "==================================================================="
echo "Working directory: $SCRIPT_DIR"
echo ""

# -----------------------------------------------------------------------------
# Check prerequisites
# -----------------------------------------------------------------------------
echo "[1/8] Checking prerequisites..."

check_tool() {
    if ! command -v "$1" &> /dev/null; then
        echo "ERROR: $1 is required but not found in PATH"
        echo "  Try: conda activate WASP2_dev2"
        exit 1
    fi
    echo "  ✓ $1 found: $(which $1)"
}

check_tool samtools
check_tool bgzip
check_tool tabix
check_tool wgsim
check_tool bwa
check_tool bcftools

echo ""

# -----------------------------------------------------------------------------
# Reference genome (reuse nf-rnaseq integration chr_test.fa)
# -----------------------------------------------------------------------------
echo "[2/8] Creating reference genome..."

INTEGRATION_FA="../../pipelines/nf-rnaseq/tests/data/integration/chr_test.fa"

if [[ -f "chr_test.fa" ]]; then
    echo "  chr_test.fa already exists, skipping"
else
    if [[ -f "$INTEGRATION_FA" ]]; then
        cp "$INTEGRATION_FA" chr_test.fa
        echo "  ✓ Copied chr_test.fa from nf-rnaseq integration ($(du -h chr_test.fa | cut -f1))"
    else
        echo "ERROR: Could not find source genome at $INTEGRATION_FA"
        exit 1
    fi
fi

# Index FASTA
if [[ ! -f "chr_test.fa.fai" ]]; then
    samtools faidx chr_test.fa
    echo "  ✓ Created chr_test.fa.fai"
fi

echo ""

# -----------------------------------------------------------------------------
# Annotation GTF (reuse from nf-rnaseq integration)
# -----------------------------------------------------------------------------
echo "[3/8] Creating annotation GTF..."

INTEGRATION_GTF="../../pipelines/nf-rnaseq/tests/data/integration/integration.gtf"

if [[ -f "annotation.gtf" ]]; then
    echo "  annotation.gtf already exists, skipping"
else
    if [[ -f "$INTEGRATION_GTF" ]]; then
        cp "$INTEGRATION_GTF" annotation.gtf
        echo "  ✓ Copied annotation.gtf from nf-rnaseq integration"
    else
        echo "ERROR: Could not find source GTF at $INTEGRATION_GTF"
        exit 1
    fi
fi

echo ""

# -----------------------------------------------------------------------------
# VCF with 10 het SNPs across 2 samples
# -----------------------------------------------------------------------------
echo "[4/8] Creating VCF with 10 het SNPs..."

if [[ -f "variants.vcf" ]]; then
    echo "  variants.vcf already exists, skipping"
else
    # Gene 1 (INTGENE001): exons at 500-1500, 2500-3500, 4500-5500 (+ strand)
    # Gene 2 (INTGENE002): exons at 10500-11500, 12500-13500, 14500-15500 (- strand)
    # Place 5 SNPs in each gene's exonic regions
    #
    # Two samples: SAMPLE1 has all 10 het, SAMPLE2 has 8 het + 2 hom-ref
    # This allows testing multi-sample handling and differential allele ratios

    cat > variants.vcf << 'EOVCF'
##fileformat=VCFv4.2
##fileDate=20260218
##source=WASP2SharedTestData
##reference=chr_test.fa
##contig=<ID=chr_test,length=19800>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
chr_test	750	snp001	C	T	100	PASS	DP=50	GT:DP	0/1:50	0/1:50
chr_test	1200	snp002	T	G	100	PASS	DP=50	GT:DP	0/1:50	0/1:50
chr_test	2800	snp003	A	C	100	PASS	DP=50	GT:DP	0/1:50	0/1:50
chr_test	3200	snp004	G	A	100	PASS	DP=50	GT:DP	0/1:50	0/0:50
chr_test	5000	snp005	G	T	100	PASS	DP=50	GT:DP	0/1:50	0/1:50
chr_test	10800	snp006	T	C	100	PASS	DP=50	GT:DP	0/1:50	0/1:50
chr_test	11200	snp007	A	G	100	PASS	DP=50	GT:DP	0/1:50	0/1:50
chr_test	12800	snp008	C	A	100	PASS	DP=50	GT:DP	0/1:50	0/0:50
chr_test	13200	snp009	G	T	100	PASS	DP=50	GT:DP	0/1:50	0/1:50
chr_test	15000	snp010	A	C	100	PASS	DP=50	GT:DP	0/1:50	0/1:50
EOVCF
    echo "  ✓ Created variants.vcf (10 het SNPs, 2 samples)"
fi

# Compress and index VCF
if [[ ! -f "variants.vcf.gz" || ! -f "variants.vcf.gz.tbi" ]]; then
    rm -f variants.vcf.gz variants.vcf.gz.tbi
    bgzip -c variants.vcf > variants.vcf.gz
    tabix -p vcf variants.vcf.gz
    echo "  ✓ Created variants.vcf.gz + .tbi"
fi

echo ""

# -----------------------------------------------------------------------------
# Regions BED from GTF exon coordinates
# -----------------------------------------------------------------------------
echo "[5/8] Creating regions BED..."

if [[ -f "regions.bed" ]]; then
    echo "  regions.bed already exists, skipping"
else
    # Extract exon coordinates from GTF → BED format
    # GTF exons from annotation.gtf:
    #   chr_test 500-1500 (exon 1, gene 1)
    #   chr_test 2500-3500 (exon 2, gene 1)
    #   chr_test 4500-5500 (exon 3, gene 1)
    #   chr_test 10500-11500 (exon 1, gene 2)
    #   chr_test 12500-13500 (exon 2, gene 2)
    #   chr_test 14500-15500 (exon 3, gene 2)
    cat > regions.bed << 'EOBED'
chr_test	499	1500	INTEXON001
chr_test	2499	3500	INTEXON002
chr_test	4499	5500	INTEXON003
chr_test	10499	11500	INTEXON004
chr_test	12499	13500	INTEXON005
chr_test	14499	15500	INTEXON006
EOBED
    echo "  ✓ Created regions.bed (6 exonic regions)"
fi

echo ""

# -----------------------------------------------------------------------------
# BWA index
# -----------------------------------------------------------------------------
echo "[6/8] Building BWA index..."

BWA_INDEX_DIR="bwa_index"
if [[ -f "${BWA_INDEX_DIR}/chr_test.fa.bwt" ]]; then
    echo "  BWA index already exists, skipping"
else
    mkdir -p "$BWA_INDEX_DIR"
    cp chr_test.fa "$BWA_INDEX_DIR/"
    bwa index "$BWA_INDEX_DIR/chr_test.fa" 2>&1 | tail -3
    echo "  ✓ Created BWA index ($(du -sh $BWA_INDEX_DIR | cut -f1))"
fi

echo ""

# -----------------------------------------------------------------------------
# Simulate reads for 3 samples and align with BWA
# -----------------------------------------------------------------------------
echo "[7/8] Simulating and aligning reads for 3 samples..."

NUM_READS=500
READ_LEN=100
ERROR_RATE=0.001

simulate_and_align() {
    local sample_name=$1
    local seed=$2
    local frag_size=$3
    local frag_std=$4

    if [[ -f "${sample_name}.bam" && -f "${sample_name}.bam.bai" ]]; then
        echo "  ${sample_name}.bam already exists, skipping"
        return
    fi

    echo "  Simulating ${sample_name} (seed=${seed}, frags=${frag_size}bp)..."

    # Simulate reads with wgsim
    wgsim -N $NUM_READS \
          -1 $READ_LEN \
          -2 $READ_LEN \
          -r 0 \
          -R 0 \
          -X 0 \
          -e $ERROR_RATE \
          -S $seed \
          -d $frag_size \
          -s $frag_std \
          chr_test.fa \
          "${sample_name}_R1.fq" \
          "${sample_name}_R2.fq" \
          > /dev/null 2>&1

    # Align with bwa mem
    bwa mem -t 2 \
        -R "@RG\tID:${sample_name}\tSM:${sample_name}\tPL:ILLUMINA\tLB:lib1" \
        "$BWA_INDEX_DIR/chr_test.fa" \
        "${sample_name}_R1.fq" \
        "${sample_name}_R2.fq" \
        2>/dev/null \
    | samtools sort -o "${sample_name}.bam" -

    samtools index "${sample_name}.bam"

    # Compress FASTQs for pipeline use
    gzip -f "${sample_name}_R1.fq"
    gzip -f "${sample_name}_R2.fq"

    local read_count=$(samtools view -c "${sample_name}.bam")
    echo "  ✓ ${sample_name}.bam: ${read_count} aligned reads ($(du -h ${sample_name}.bam | cut -f1))"
}

# Sample1: standard RNA-seq-like fragments (seed 42)
simulate_and_align "sample1" 42 300 50

# Sample2: slightly different fragment size (seed 43)
simulate_and_align "sample2" 43 250 40

# Sample3: third sample for OUTRIDER (seed 44)
simulate_and_align "sample3" 44 280 45

echo ""

# -----------------------------------------------------------------------------
# Validate all generated files
# -----------------------------------------------------------------------------
echo "[8/8] Validating generated files..."

ERRORS=0

validate_file() {
    local filepath=$1
    local min_size=${2:-1}  # minimum expected size in bytes

    if [[ -f "$filepath" ]]; then
        local size=$(stat -c%s "$filepath" 2>/dev/null || stat -f%z "$filepath" 2>/dev/null)
        if [[ $size -ge $min_size ]]; then
            echo "  ✓ $filepath ($(du -h "$filepath" | cut -f1))"
        else
            echo "  ✗ $filepath exists but too small (${size} bytes, expected >= ${min_size})"
            ERRORS=$((ERRORS + 1))
        fi
    else
        echo "  ✗ $filepath NOT FOUND"
        ERRORS=$((ERRORS + 1))
    fi
}

validate_bam() {
    local bam=$1
    if samtools quickcheck "$bam" 2>/dev/null; then
        local count=$(samtools view -c "$bam")
        echo "  ✓ $bam passes quickcheck (${count} reads)"
    else
        echo "  ✗ $bam FAILS quickcheck"
        ERRORS=$((ERRORS + 1))
    fi
}

echo ""
echo "  --- Reference files ---"
validate_file "chr_test.fa" 19000
validate_file "chr_test.fa.fai" 10

echo ""
echo "  --- Variant files ---"
validate_file "variants.vcf" 500
validate_file "variants.vcf.gz" 100
validate_file "variants.vcf.gz.tbi" 50

echo ""
echo "  --- Annotation files ---"
validate_file "annotation.gtf" 500
validate_file "regions.bed" 100

echo ""
echo "  --- BWA index ---"
validate_file "bwa_index/chr_test.fa.bwt" 1000

echo ""
echo "  --- BAM files ---"
validate_file "sample1.bam" 100
validate_file "sample1.bam.bai" 50
validate_bam "sample1.bam"
validate_file "sample2.bam" 100
validate_file "sample2.bam.bai" 50
validate_bam "sample2.bam"
validate_file "sample3.bam" 100
validate_file "sample3.bam.bai" 50
validate_bam "sample3.bam"

echo ""
echo "  --- FASTQ files ---"
validate_file "sample1_R1.fq.gz" 1000
validate_file "sample1_R2.fq.gz" 1000
validate_file "sample2_R1.fq.gz" 1000
validate_file "sample2_R2.fq.gz" 1000
validate_file "sample3_R1.fq.gz" 1000
validate_file "sample3_R2.fq.gz" 1000

echo ""
if [[ $ERRORS -eq 0 ]]; then
    echo "==================================================================="
    echo " SUCCESS! All shared core test data generated."
    echo "==================================================================="
else
    echo "==================================================================="
    echo " WARNING: $ERRORS validation errors found."
    echo "==================================================================="
    exit 1
fi

echo ""
echo "Total disk usage: $(du -sh . | cut -f1)"
echo ""
echo "Files ready for per-pipeline generators:"
echo "  pipelines/nf-atacseq/tests/data/generate_test_data.sh"
echo "  pipelines/nf-scatac/tests/data/generate_test_data.sh"
echo "  pipelines/nf-outrider/tests/data/generate_test_data.sh"
echo "  galaxy/tools/wasp2/test-data/generate_test_data.sh"
echo ""
