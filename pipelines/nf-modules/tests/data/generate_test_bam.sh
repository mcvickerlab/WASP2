#!/bin/bash
# Generate minimal synthetic BAM for WASP2 nf-test
# Creates paired-end reads overlapping heterozygous variants from sample.vcf
#
# Variants in sample.vcf (sample1 heterozygous sites):
#   chr1:100 (rs1, A>G, 0/1)
#   chr1:400 (rs4, T>C, 0/1)
#   chr2:100 (rs5, A>T, 0/1)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Check for samtools
if ! command -v samtools &> /dev/null; then
    echo "ERROR: samtools is required but not found"
    exit 1
fi

echo "Creating minimal test SAM file..."

# Create SAM with 50bp reads (CIGAR 50M must match 50bp sequence)
cat > minimal.sam << 'EOSAM'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:248956422
@SQ	SN:chr2	LN:242193529
@RG	ID:test	SM:sample1	PL:ILLUMINA
@PG	ID:bwa	PN:bwa	VN:0.7.17
read001	99	chr1	90	60	50M	=	140	100	GGGGGGGGGGAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	RG:Z:test	NM:i:0
read001	147	chr1	140	60	50M	=	90	-100	GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	RG:Z:test	NM:i:0
read002	99	chr1	390	60	50M	=	440	100	GGGGGGGGGGTGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	RG:Z:test	NM:i:0
read002	147	chr1	440	60	50M	=	390	-100	GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	RG:Z:test	NM:i:0
read003	99	chr2	90	60	50M	=	140	100	GGGGGGGGGGAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	RG:Z:test	NM:i:0
read003	147	chr2	140	60	50M	=	90	-100	GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	RG:Z:test	NM:i:0
EOSAM

echo "Converting SAM to BAM..."
samtools view -bS minimal.sam > minimal_unsorted.bam

echo "Sorting BAM..."
samtools sort -o minimal.bam minimal_unsorted.bam

echo "Indexing BAM..."
samtools index minimal.bam

# Create a copy for "remapped" BAM (simulates user remapping step)
cp minimal.bam minimal_remap.bam
samtools index minimal_remap.bam

# Cleanup
rm -f minimal.sam minimal_unsorted.bam

echo "Created test files:"
ls -la minimal*.bam*

echo "BAM statistics:"
samtools flagstat minimal.bam
