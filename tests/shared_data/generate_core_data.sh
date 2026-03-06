#!/bin/bash
# =============================================================================
# WASP2 Shared Core Test Data Generator
# =============================================================================
# Generates the unified test dataset used by all WASP2 pipelines, Galaxy tools,
# CLI smoke tests, and container validation.
#
# Outputs (all committed to git, ~700K total):
#   chr_test.fa + .fai          - 20kb random reference genome (high complexity)
#   variants.vcf + .gz + .tbi   - 30 het SNPs across 3 samples (phased)
#   annotation.gtf              - 12 genes, 16 exons
#   regions.bed                 - Peak/region file from exon coordinates
#   sample{1,2,3}.bam + .bai   - Aligned reads (wgsim + bwa)
#   bwa_index/                  - BWA index for chr_test.fa
#   expected_counts.tsv         - WASP2 counting output baseline
#   expected_analysis.tsv       - WASP2 analysis output baseline (placeholder)
#
# Prerequisites: python3, samtools, bgzip, tabix, wgsim, bwa, bcftools
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
    echo "  OK $1 found: $(which $1)"
}

check_tool python3
check_tool samtools
check_tool bgzip
check_tool tabix
check_tool wgsim
check_tool bwa
check_tool bcftools

echo ""

# -----------------------------------------------------------------------------
# Clean up old generated files (force full regeneration)
# -----------------------------------------------------------------------------
echo "[2/8] Generating random reference genome via Python..."

# Remove old reference and index to force regeneration
rm -f chr_test.fa chr_test.fa.fai

# Generate a realistic random reference with the Python script
python3 generate_reference.py > chr_test.fa

# Index FASTA
samtools faidx chr_test.fa
echo "  OK Created chr_test.fa ($(du -h chr_test.fa | cut -f1)) + .fai"

echo ""

# -----------------------------------------------------------------------------
# Annotation GTF (keep existing if present, otherwise copy from integration)
# -----------------------------------------------------------------------------
echo "[3/8] Checking annotation GTF..."

INTEGRATION_GTF="../../pipelines/nf-rnaseq/tests/data/integration/integration.gtf"

if [[ -f "annotation.gtf" ]]; then
    echo "  annotation.gtf already exists, keeping"
else
    if [[ -f "$INTEGRATION_GTF" ]]; then
        cp "$INTEGRATION_GTF" annotation.gtf
        echo "  OK Copied annotation.gtf from nf-rnaseq integration"
    else
        echo "ERROR: Could not find source GTF at $INTEGRATION_GTF"
        exit 1
    fi
fi

echo ""

# -----------------------------------------------------------------------------
# VCF with 30 het SNPs across 3 samples (REF alleles from actual reference)
# -----------------------------------------------------------------------------
echo "[4/8] Creating VCF with 30 phased het SNPs..."

# Read actual reference bases at SNP positions using Python
# This ensures REF alleles match the generated reference exactly
rm -f variants.vcf variants.vcf.gz variants.vcf.gz.tbi

python3 - << 'PYEOF'
import sys

# Read the reference sequence
with open("chr_test.fa") as f:
    lines = f.readlines()
seq = ''.join(line.strip() for line in lines if not line.startswith('>'))

# Deterministic ALT allele mapping (always different from REF)
alt_map = {'A': 'C', 'T': 'G', 'G': 'T', 'C': 'A'}

# SNP positions (1-based) spread across gene regions in annotation.gtf:
#   Gene1 (500-5500):   750, 1200, 2800, 3200, 5000
#   Gene3 (5800-6300):  6000, 6100
#   Gene4 (6500-7000):  6700, 6800
#   Gene5 (7200-7700):  7400, 7500
#   Gene6 (7900-8400):  8100, 8200
#   Gene7 (8600-9100):  8800, 8900
#   Gene2 (10500-15500): 10800, 11200, 12800, 13200, 15000
#   Gene8 (15800-16300): 16000, 16100
#   Gene9 (16500-17000): 16700, 16800
#   Gene10 (17200-17700): 17400, 17500
#   Gene11 (17900-18400): 18100, 18200
#   Gene12 (18600-19100): 18800, 18900
snps = [
    (750,   "snp001"),
    (1200,  "snp002"),
    (2800,  "snp003"),
    (3200,  "snp004"),
    (5000,  "snp005"),
    (6000,  "snp011"),
    (6100,  "snp012"),
    (6700,  "snp013"),
    (6800,  "snp014"),
    (7400,  "snp015"),
    (7500,  "snp016"),
    (8100,  "snp017"),
    (8200,  "snp018"),
    (8800,  "snp019"),
    (8900,  "snp020"),
    (10800, "snp006"),
    (11200, "snp007"),
    (12800, "snp008"),
    (13200, "snp009"),
    (15000, "snp010"),
    (16000, "snp021"),
    (16100, "snp022"),
    (16700, "snp023"),
    (16800, "snp024"),
    (17400, "snp025"),
    (17500, "snp026"),
    (18100, "snp027"),
    (18200, "snp028"),
    (18800, "snp029"),
    (18900, "snp030"),
]

# Genotype patterns for 3 samples (phased, 0|1)
# sample1: 28 het, 2 hom-ref   (high het)
# sample2: 22 het, 8 hom-ref   (medium het)
# sample3: 22 het, 8 hom-ref   (medium het, different SNPs from sample2)
genotypes = [
    # snp001-005 (Gene1)
    ("0|1", "0|1", "0|1"),
    ("0|1", "0|1", "0|1"),
    ("0|1", "0|1", "0|0"),
    ("0|1", "0|0", "0|0"),
    ("0|1", "0|1", "0|1"),
    # snp011-012 (Gene3)
    ("0|1", "0|1", "0|0"),
    ("0|1", "0|0", "0|1"),
    # snp013-014 (Gene4)
    ("0|1", "0|1", "0|1"),
    ("0|0", "0|1", "0|1"),
    # snp015-016 (Gene5)
    ("0|1", "0|1", "0|0"),
    ("0|1", "0|0", "0|1"),
    # snp017-018 (Gene6)
    ("0|1", "0|1", "0|1"),
    ("0|0", "0|1", "0|1"),
    # snp019-020 (Gene7)
    ("0|1", "0|1", "0|1"),
    ("0|1", "0|0", "0|1"),
    # snp006-010 (Gene2)
    ("0|1", "0|1", "0|0"),
    ("0|1", "0|1", "0|1"),
    ("0|1", "0|0", "0|0"),
    ("0|1", "0|1", "0|1"),
    ("0|1", "0|1", "0|0"),
    # snp021-022 (Gene8)
    ("0|1", "0|1", "0|1"),
    ("0|1", "0|0", "0|1"),
    # snp023-024 (Gene9)
    ("0|1", "0|1", "0|0"),
    ("0|0", "0|1", "0|1"),
    # snp025-026 (Gene10)
    ("0|1", "0|1", "0|1"),
    ("0|1", "0|0", "0|1"),
    # snp027-028 (Gene11)
    ("0|1", "0|1", "0|0"),
    ("0|0", "0|1", "0|1"),
    # snp029-030 (Gene12)
    ("0|1", "0|1", "0|1"),
    ("0|1", "0|0", "0|1"),
]

with open("variants.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n")
    f.write("##fileDate=20260306\n")
    f.write("##source=WASP2SharedTestData\n")
    f.write("##reference=chr_test.fa\n")
    f.write(f"##contig=<ID=chr_test,length={len(seq)}>\n")
    f.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n')
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2\tsample3\n")

    for i, (pos, snp_id) in enumerate(snps):
        ref = seq[pos - 1]  # 0-based index
        alt = alt_map[ref]
        gt1, gt2, gt3 = genotypes[i]
        f.write(f"chr_test\t{pos}\t{snp_id}\t{ref}\t{alt}\t100\tPASS\tDP=50\tGT:DP\t{gt1}:50\t{gt2}:50\t{gt3}:50\n")

print(f"Created variants.vcf with {len(snps)} SNPs, REF alleles verified against reference", file=sys.stderr)
PYEOF

echo "  OK Created variants.vcf (30 het SNPs, 3 samples, phased)"

# Compress and index VCF
bgzip -c variants.vcf > variants.vcf.gz
tabix -p vcf variants.vcf.gz
echo "  OK Created variants.vcf.gz + .tbi"

echo ""

# -----------------------------------------------------------------------------
# Regions BED from GTF exon coordinates
# -----------------------------------------------------------------------------
echo "[5/8] Creating regions BED..."

if [[ -f "regions.bed" ]]; then
    echo "  regions.bed already exists, keeping"
else
    cat > regions.bed << 'EOBED'
chr_test	499	1500	INTEXON001
chr_test	2499	3500	INTEXON002
chr_test	4499	5500	INTEXON003
chr_test	10499	11500	INTEXON004
chr_test	12499	13500	INTEXON005
chr_test	14499	15500	INTEXON006
EOBED
    echo "  OK Created regions.bed (6 exonic regions)"
fi

echo ""

# -----------------------------------------------------------------------------
# BWA index
# -----------------------------------------------------------------------------
echo "[6/8] Building BWA index..."

BWA_INDEX_DIR="bwa_index"
# Always rebuild BWA index when reference changes
rm -rf "$BWA_INDEX_DIR"
mkdir -p "$BWA_INDEX_DIR"
cp chr_test.fa "$BWA_INDEX_DIR/"
bwa index "$BWA_INDEX_DIR/chr_test.fa" 2>&1 | tail -3
echo "  OK Created BWA index ($(du -sh $BWA_INDEX_DIR | cut -f1))"

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

    # Always regenerate when reference changes
    rm -f "${sample_name}.bam" "${sample_name}.bam.bai"
    rm -f "${sample_name}_R1.fq" "${sample_name}_R2.fq"
    rm -f "${sample_name}_R1.fq.gz" "${sample_name}_R2.fq.gz"

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
    echo "  OK ${sample_name}.bam: ${read_count} aligned reads ($(du -h ${sample_name}.bam | cut -f1))"
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
            echo "  OK $filepath ($(du -h "$filepath" | cut -f1))"
        else
            echo "  FAIL $filepath exists but too small (${size} bytes, expected >= ${min_size})"
            ERRORS=$((ERRORS + 1))
        fi
    else
        echo "  FAIL $filepath NOT FOUND"
        ERRORS=$((ERRORS + 1))
    fi
}

validate_bam() {
    local bam=$1
    if samtools quickcheck "$bam" 2>/dev/null; then
        local count=$(samtools view -c "$bam")
        echo "  OK $bam passes quickcheck (${count} reads)"
    else
        echo "  FAIL $bam FAILS quickcheck"
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

# -----------------------------------------------------------------------------
# Quality report: MAPQ + proper pairing
# -----------------------------------------------------------------------------
echo ""
echo "  --- Alignment quality report ---"
for sample in sample1 sample2 sample3; do
    total=$(samtools view -c "${sample}.bam")
    mapq_gt0=$(samtools view -c -q 1 "${sample}.bam")
    mapq_pct=$(python3 -c "print(f'{100*${mapq_gt0}/${total}:.1f}')")
    proper=$(samtools view -c -f 2 "${sample}.bam")
    proper_pct=$(python3 -c "print(f'{100*${proper}/${total}:.1f}')")
    echo "  ${sample}: ${total} total, ${mapq_gt0} MAPQ>0 (${mapq_pct}%), ${proper} properly paired (${proper_pct}%)"
done

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
