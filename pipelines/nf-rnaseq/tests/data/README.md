# Test Data for nf-rnaseq Pipeline

This directory contains test data for the WASP2 RNA-seq ASE pipeline tests.

## Directory Structure

```
data/
├── README.md             # This file
├── test_snps.vcf.gz      # Minimal VCF with test SNPs (6 het SNPs, 2 samples)
├── test_snps.vcf.gz.tbi  # Tabix index
├── test.gtf              # Minimal GTF annotation (2 genes, 6 exons)
├── sample1_R1.fq.gz      # Placeholder FASTQ for SAMPLE1 R1
├── sample1_R2.fq.gz      # Placeholder FASTQ for SAMPLE1 R2
├── sample2_R1.fq.gz      # Placeholder FASTQ for SAMPLE2 R1
├── sample2_R2.fq.gz      # Placeholder FASTQ for SAMPLE2 R2
└── star_index/           # Placeholder STAR index directory
    └── .gitkeep
```

## For Stub Tests

Stub tests (`-stub` mode) don't require real data files - they use process stub
blocks that generate placeholder outputs. The file paths just need to exist.

## For Integration Tests

For full integration tests with real data, use the HG00731 RNA-seq test dataset:

```bash
# BAM file location
benchmarking/star_wasp_comparison/results/unified_2025-12-04_00-29-39/A_sorted.bam

# VCF file location
benchmarking/star_wasp_comparison/data/HG00731_het_only_chr.vcf.gz
```

## Generating Test Data

To regenerate minimal test data:

```bash
# Create minimal VCF
echo -e "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=248956422>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1" | bgzip > test_snps.vcf.gz
tabix -p vcf test_snps.vcf.gz

# Create minimal GTF
echo -e 'chr1\ttest\texon\t1000\t2000\t.\t+\t.\tgene_id "TEST001"; transcript_id "TEST001.1";' > test.gtf
```
