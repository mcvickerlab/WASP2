#!/bin/bash
# =============================================================================
# Generate enhanced test BAMs for WASP2 nf-modules nf-test
# =============================================================================
# Creates paired-end reads overlapping heterozygous variants from sample.vcf
# with both REF and ALT alleles for proper allele counting validation.
#
# Variants in sample.vcf (sample1 heterozygous sites):
#   chr1:100 (rs1, A>G, 0/1)
#   chr1:400 (rs4, T>C, 0/1)
#   chr2:100 (rs5, A>T, 0/1)
#
# Each het site gets 4 reads: 2 with REF allele, 2 with ALT allele
# Plus 2 read pairs not overlapping any variant
# Total: 14 read pairs (28 alignments)
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

if ! command -v samtools &> /dev/null; then
    echo "ERROR: samtools is required but not found"
    exit 1
fi

echo "Creating enhanced test SAM file..."

# Use Python to generate SAM with exact-length sequences to avoid SEQ/QUAL mismatch
python3 << 'PYEOF'
READ_LEN = 50
QUAL = "I" * READ_LEN
BASE = "G" * READ_LEN

# Header
header = [
    "@HD\tVN:1.6\tSO:coordinate",
    "@SQ\tSN:chr1\tLN:248956422",
    "@SQ\tSN:chr2\tLN:242193529",
    "@RG\tID:test\tSM:sample1\tPL:ILLUMINA",
    "@PG\tID:bwa\tPN:bwa\tVN:0.7.17",
]

reads = []

def make_read(name, flag, chrom, pos, mate_pos, tlen, seq=None, nm=0):
    """Generate a SAM record. seq defaults to all G's if not provided."""
    if seq is None:
        seq = BASE
    assert len(seq) == READ_LEN, f"{name}: seq len {len(seq)} != {READ_LEN}"
    return f"{name}\t{flag}\t{chrom}\t{pos}\t60\t{READ_LEN}M\t=\t{mate_pos}\t{tlen}\t{seq}\t{QUAL}\tRG:Z:test\tNM:i:{nm}"

def inject_base(seq, offset, base):
    """Replace base at offset in sequence."""
    s = list(seq)
    s[offset] = base
    return "".join(s)

# --- chr1:100 (A>G) - offset from read start ---
# REF reads: A at position 100 (offset = 100 - start)
for suffix, start in [("ref1", 80), ("ref2", 75)]:
    offset = 100 - start  # 20 or 25
    seq = inject_base(BASE, offset, "A")
    reads.append(make_read(f"c1p100_{suffix}", 99, "chr1", start, start+100, 150, seq, 0))
    reads.append(make_read(f"c1p100_{suffix}", 147, "chr1", start+100, start, -150, BASE, 0))

# ALT reads: G at position 100
for suffix, start in [("alt1", 80), ("alt2", 75)]:
    offset = 100 - start
    seq = inject_base(BASE, offset, "G")
    reads.append(make_read(f"c1p100_{suffix}", 99, "chr1", start, start+100, 150, seq, 1))
    reads.append(make_read(f"c1p100_{suffix}", 147, "chr1", start+100, start, -150, BASE, 0))

# --- chr1:400 (T>C) ---
for suffix, start in [("ref1", 380), ("ref2", 375)]:
    offset = 400 - start
    seq = inject_base(BASE, offset, "T")
    reads.append(make_read(f"c1p400_{suffix}", 99, "chr1", start, start+100, 150, seq, 0))
    reads.append(make_read(f"c1p400_{suffix}", 147, "chr1", start+100, start, -150, BASE, 0))

for suffix, start in [("alt1", 380), ("alt2", 375)]:
    offset = 400 - start
    seq = inject_base(BASE, offset, "C")
    reads.append(make_read(f"c1p400_{suffix}", 99, "chr1", start, start+100, 150, seq, 1))
    reads.append(make_read(f"c1p400_{suffix}", 147, "chr1", start+100, start, -150, BASE, 0))

# --- chr2:100 (A>T) ---
for suffix, start in [("ref1", 80), ("ref2", 75)]:
    offset = 100 - start
    seq = inject_base(BASE, offset, "A")
    reads.append(make_read(f"c2p100_{suffix}", 99, "chr2", start, start+100, 150, seq, 0))
    reads.append(make_read(f"c2p100_{suffix}", 147, "chr2", start+100, start, -150, BASE, 0))

for suffix, start in [("alt1", 80), ("alt2", 75)]:
    offset = 100 - start
    seq = inject_base(BASE, offset, "T")
    reads.append(make_read(f"c2p100_{suffix}", 99, "chr2", start, start+100, 150, seq, 1))
    reads.append(make_read(f"c2p100_{suffix}", 147, "chr2", start+100, start, -150, BASE, 0))

# --- Non-variant reads ---
reads.append(make_read("novar1", 99, "chr1", 500, 600, 150, BASE, 0))
reads.append(make_read("novar1", 147, "chr1", 600, 500, -150, BASE, 0))
reads.append(make_read("novar2", 99, "chr2", 500, 600, 150, BASE, 0))
reads.append(make_read("novar2", 147, "chr2", 600, 500, -150, BASE, 0))

with open("minimal.sam", "w") as f:
    for line in header:
        f.write(line + "\n")
    for line in reads:
        f.write(line + "\n")

print(f"Generated {len(reads)} SAM records ({len(reads)//2} read pairs)")
PYEOF

echo "Converting SAM to BAM..."
samtools view -bS minimal.sam > minimal_unsorted.bam

echo "Sorting BAM..."
samtools sort -o minimal.bam minimal_unsorted.bam

echo "Indexing BAM..."
samtools index minimal.bam

# Create remapped BAM (simulates user remapping — same reads, same positions)
cp minimal.bam minimal_remap.bam
samtools index minimal_remap.bam

# Update wasp_data.json with all read mappings
cat > sample.wasp_data.json << 'EOJSON'
{
  "bam_file": "minimal.bam",
  "variant_file": "sample.vcf.gz",
  "sample": "sample1",
  "phased": false,
  "read_mappings": {
    "c1p100_ref1": {
      "original_pos": 80,
      "chrom": "chr1",
      "variants": [{"pos": 100, "ref": "A", "alt": "G"}]
    },
    "c1p100_ref2": {
      "original_pos": 75,
      "chrom": "chr1",
      "variants": [{"pos": 100, "ref": "A", "alt": "G"}]
    },
    "c1p100_alt1": {
      "original_pos": 80,
      "chrom": "chr1",
      "variants": [{"pos": 100, "ref": "A", "alt": "G"}]
    },
    "c1p100_alt2": {
      "original_pos": 75,
      "chrom": "chr1",
      "variants": [{"pos": 100, "ref": "A", "alt": "G"}]
    },
    "c1p400_ref1": {
      "original_pos": 380,
      "chrom": "chr1",
      "variants": [{"pos": 400, "ref": "T", "alt": "C"}]
    },
    "c1p400_ref2": {
      "original_pos": 375,
      "chrom": "chr1",
      "variants": [{"pos": 400, "ref": "T", "alt": "C"}]
    },
    "c1p400_alt1": {
      "original_pos": 380,
      "chrom": "chr1",
      "variants": [{"pos": 400, "ref": "T", "alt": "C"}]
    },
    "c1p400_alt2": {
      "original_pos": 375,
      "chrom": "chr1",
      "variants": [{"pos": 400, "ref": "T", "alt": "C"}]
    },
    "c2p100_ref1": {
      "original_pos": 80,
      "chrom": "chr2",
      "variants": [{"pos": 100, "ref": "A", "alt": "T"}]
    },
    "c2p100_ref2": {
      "original_pos": 75,
      "chrom": "chr2",
      "variants": [{"pos": 100, "ref": "A", "alt": "T"}]
    },
    "c2p100_alt1": {
      "original_pos": 80,
      "chrom": "chr2",
      "variants": [{"pos": 100, "ref": "A", "alt": "T"}]
    },
    "c2p100_alt2": {
      "original_pos": 75,
      "chrom": "chr2",
      "variants": [{"pos": 100, "ref": "A", "alt": "T"}]
    }
  }
}
EOJSON

# Cleanup
rm -f minimal.sam minimal_unsorted.bam

echo ""
echo "Created test files:"
ls -la minimal*.bam* sample.wasp_data.json

echo ""
echo "BAM statistics:"
samtools flagstat minimal.bam
