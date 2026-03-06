#!/bin/bash
# =============================================================================
# WASP2 nf-atacseq Test Data Generator (v2 — realistic reference)
# =============================================================================
# Generates self-contained ATAC-seq test data with a non-repetitive reference
# so BWA alignment produces meaningful mapping rates (>80%).
#
# Previous version used the shared chr_test.fa which is a repetitive ATGC
# pattern yielding ~0% mapping. This version generates its own reference.
#
# To produce non-zero allele counts, reads are simulated from BOTH haplotypes:
# half from the REF haplotype, half from an ALT haplotype with het SNPs applied.
#
# Prerequisites: python3, samtools, bgzip, tabix, wgsim, bwa
#   (all available in WASP2_dev2 conda env or WASP2 micromamba env)
#
# Usage:
#   cd pipelines/nf-atacseq/tests/data
#   bash generate_test_data.sh
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# BWA may not be in PATH; check common conda/micromamba locations
if ! command -v bwa &>/dev/null; then
    for candidate in \
        /usr/local/Cellar/micromamba/*/envs/WASP2/bin/bwa \
        /usr/local/Cellar/micromamba/*/envs/WASP2_dev2/bin/bwa \
        "${HOME}/miniforge3/envs/WASP2/bin/bwa" \
        "${HOME}/miniconda3/envs/WASP2/bin/bwa"; do
        if [[ -x "$candidate" ]]; then
            export PATH="$(dirname "$candidate"):$PATH"
            break
        fi
    done
fi

echo "==================================================================="
echo " WASP2 nf-atacseq Test Data Generator (v2)"
echo "==================================================================="

# -----------------------------------------------------------------------------
# Check prerequisites
# -----------------------------------------------------------------------------
echo "[0/7] Checking prerequisites..."

check_tool() {
    if ! command -v "$1" &>/dev/null; then
        echo "ERROR: $1 is required but not found in PATH"
        echo "  Try: conda activate WASP2_dev2"
        exit 1
    fi
    echo "  OK: $1"
}

check_tool python3
check_tool samtools
check_tool bwa
check_tool wgsim
check_tool bgzip
check_tool tabix
echo ""

# -----------------------------------------------------------------------------
# Clean stale symlinks and old data (one-time migration from v1)
# -----------------------------------------------------------------------------
echo "[1/7] Cleaning stale data..."
for f in chr_test.fa chr_test.fa.fai variants.vcf.gz variants.vcf.gz.tbi annotation.gtf regions.bed; do
    if [[ -L "$f" ]]; then
        rm -f "$f"
        echo "  Removed symlink: $f"
    fi
done
rm -rf bwa_index
rm -f sample1_R1.fq.gz sample1_R2.fq.gz
rm -f chr_test.fa chr_test.fa.fai variants.vcf variants.vcf.gz variants.vcf.gz.tbi regions.bed
echo "  Cleaned previous outputs"
echo ""

# -----------------------------------------------------------------------------
# Generate realistic non-repetitive reference
# -----------------------------------------------------------------------------
echo "[2/7] Generating realistic reference genome..."
python3 "${SCRIPT_DIR}/generate_realistic_reference.py" chr_test.fa
samtools faidx chr_test.fa
echo "  Created chr_test.fa + .fai"
echo ""

# -----------------------------------------------------------------------------
# Generate VCF with ~30 het SNPs + ALT haplotype reference
# -----------------------------------------------------------------------------
echo "[3/7] Creating VCF with 30 het SNPs and ALT haplotype..."

python3 - <<'PYEOF'
import random

# Read reference
with open("chr_test.fa") as f:
    lines = f.readlines()
seq = "".join(l.strip() for l in lines[1:])

# Deterministic SNP positions spread across the reference
rng = random.Random(99)
positions = sorted(rng.sample(range(200, 19800), 30))

# Transition mapping for plausible variants
transitions = {"A": "G", "G": "A", "T": "C", "C": "T"}

# --- Write VCF ---
vcf_lines = []
vcf_lines.append("##fileformat=VCFv4.2")
vcf_lines.append("##source=WASP2_nf_atacseq_test_data_v2")
vcf_lines.append("##reference=chr_test.fa")
vcf_lines.append("##contig=<ID=chr_test,length=20000>")
vcf_lines.append('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">')
vcf_lines.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
vcf_lines.append('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">')
vcf_lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1")

snp_map = {}  # pos -> (ref, alt)
for i, pos in enumerate(positions):
    ref = seq[pos - 1]  # 1-based
    alt = transitions[ref]
    snp_id = f"snp{i+1:03d}"
    vcf_lines.append(
        f"chr_test\t{pos}\t{snp_id}\t{ref}\t{alt}\t100\tPASS\tDP=50\tGT:DP\t0|1:50"
    )
    snp_map[pos] = (ref, alt)

with open("variants.vcf", "w") as f:
    f.write("\n".join(vcf_lines) + "\n")

print(f"  Created variants.vcf with {len(positions)} het SNPs")

# --- Write ALT haplotype reference ---
alt_seq = list(seq)
for pos, (ref, alt) in snp_map.items():
    assert alt_seq[pos - 1] == ref, f"Mismatch at {pos}: expected {ref}, got {alt_seq[pos-1]}"
    alt_seq[pos - 1] = alt

with open("chr_test_alt.fa", "w") as f:
    f.write(">chr_test\n")
    alt_str = "".join(alt_seq)
    for i in range(0, len(alt_str), 60):
        f.write(alt_str[i:i+60] + "\n")

print(f"  Created chr_test_alt.fa (ALT haplotype with {len(snp_map)} substitutions)")
PYEOF

# Compress and index
rm -f variants.vcf.gz variants.vcf.gz.tbi
bgzip -c variants.vcf > variants.vcf.gz
tabix -p vcf variants.vcf.gz
echo "  Created variants.vcf.gz + .tbi"
echo ""

# -----------------------------------------------------------------------------
# Create regions BED covering all SNP positions
# -----------------------------------------------------------------------------
echo "[4/7] Creating regions BED..."

python3 - <<'PYEOF'
import random

rng = random.Random(99)
positions = sorted(rng.sample(range(200, 19800), 30))

# Create ~500bp regions centered on each SNP, merge overlapping
regions = []
for pos in positions:
    start = max(0, pos - 250)
    end = min(20000, pos + 250)
    regions.append((start, end))

# Merge overlapping regions
merged = [regions[0]]
for start, end in regions[1:]:
    if start <= merged[-1][1]:
        merged[-1] = (merged[-1][0], max(merged[-1][1], end))
    else:
        merged.append((start, end))

with open("regions.bed", "w") as f:
    for i, (start, end) in enumerate(merged):
        f.write(f"chr_test\t{start}\t{end}\tpeak_{i+1}\n")

print(f"  Created regions.bed with {len(merged)} peaks covering {len(positions)} SNPs")
PYEOF
echo ""

# -----------------------------------------------------------------------------
# Simulate ATAC-seq reads from BOTH haplotypes (REF + ALT)
# -----------------------------------------------------------------------------
echo "[5/7] Simulating ATAC-seq paired-end reads (dual haplotype)..."

# 20kb genome, 75bp reads, ~20x total coverage
# Split: ~1350 pairs from REF, ~1350 pairs from ALT
NUM_READS_PER_HAP=1350
READ_LEN=75
FRAG_SIZE=180
FRAG_STD=30
ERROR_RATE=0.001

# Simulate from REF haplotype
wgsim -N $NUM_READS_PER_HAP \
      -1 $READ_LEN \
      -2 $READ_LEN \
      -r 0 -R 0 -X 0 \
      -e $ERROR_RATE \
      -S 100 \
      -d $FRAG_SIZE \
      -s $FRAG_STD \
      chr_test.fa \
      ref_R1.fq \
      ref_R2.fq \
      > /dev/null 2>&1
echo "  Simulated ${NUM_READS_PER_HAP} pairs from REF haplotype"

# Simulate from ALT haplotype
wgsim -N $NUM_READS_PER_HAP \
      -1 $READ_LEN \
      -2 $READ_LEN \
      -r 0 -R 0 -X 0 \
      -e $ERROR_RATE \
      -S 200 \
      -d $FRAG_SIZE \
      -s $FRAG_STD \
      chr_test_alt.fa \
      alt_R1.fq \
      alt_R2.fq \
      > /dev/null 2>&1
echo "  Simulated ${NUM_READS_PER_HAP} pairs from ALT haplotype"

# Combine and compress
cat ref_R1.fq alt_R1.fq | gzip -c > sample1_R1.fq.gz
cat ref_R2.fq alt_R2.fq | gzip -c > sample1_R2.fq.gz
echo "  Combined into sample1_R{1,2}.fq.gz ($((NUM_READS_PER_HAP * 2)) total pairs)"

# Clean up temporary files
rm -f ref_R1.fq ref_R2.fq alt_R1.fq alt_R2.fq chr_test_alt.fa
echo ""

# -----------------------------------------------------------------------------
# Build BWA index
# -----------------------------------------------------------------------------
echo "[6/7] Building BWA index..."

BWA_INDEX_DIR="bwa_index"
mkdir -p "$BWA_INDEX_DIR"
cp chr_test.fa "$BWA_INDEX_DIR/"
bwa index "$BWA_INDEX_DIR/chr_test.fa" 2>&1 | tail -2
echo "  Created BWA index"
echo ""

# -----------------------------------------------------------------------------
# Create samplesheets (both test and local variants)
# -----------------------------------------------------------------------------
echo "[7/7] Creating samplesheets..."

# test samplesheet uses absolute paths
cat > samplesheet_test.csv << EOF
sample,fastq_1,fastq_2,sample_name
test_sample1,${SCRIPT_DIR}/sample1_R1.fq.gz,${SCRIPT_DIR}/sample1_R2.fq.gz,sample1
EOF
echo "  Created samplesheet_test.csv"

# local samplesheet uses ${projectDir} relative paths (for nextflow)
cat > samplesheet_local.csv << 'EOF'
sample,fastq_1,fastq_2,sample_name
test_sample1,${projectDir}/tests/data/sample1_R1.fq.gz,${projectDir}/tests/data/sample1_R2.fq.gz,sample1
EOF
echo "  Created samplesheet_local.csv"

# -----------------------------------------------------------------------------
# Quick validation
# -----------------------------------------------------------------------------
echo ""
echo "==================================================================="
echo " Validation"
echo "==================================================================="

# Check BWA alignment quality
echo ""
echo "--- Quick alignment test (first 100 pairs) ---"
bwa mem -t 2 \
    -R "@RG\tID:sample1\tSM:sample1\tPL:ILLUMINA\tLB:lib1" \
    "$BWA_INDEX_DIR/chr_test.fa" \
    <(gunzip -c sample1_R1.fq.gz | head -400) \
    <(gunzip -c sample1_R2.fq.gz | head -400) \
    2>/dev/null \
| samtools flagstat - 2>/dev/null

echo ""

# Check VCF REF alleles match reference
echo "--- VCF REF allele validation ---"
python3 - <<'PYEOF'
seq_lines = open("chr_test.fa").readlines()
seq = "".join(l.strip() for l in seq_lines[1:])

errors = 0
total = 0
with open("variants.vcf") as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        pos = int(fields[1])
        ref = fields[3]
        actual = seq[pos - 1]
        total += 1
        if ref != actual:
            print(f"  MISMATCH at pos {pos}: VCF REF={ref}, actual={actual}")
            errors += 1

if errors == 0:
    print(f"  All {total} REF alleles match reference")
else:
    print(f"  {errors}/{total} mismatches found!")
PYEOF

echo ""
echo "==================================================================="
echo " SUCCESS! nf-atacseq test data generated (v2)."
echo "==================================================================="
echo "Total: $(du -sh . | cut -f1)"
echo ""
