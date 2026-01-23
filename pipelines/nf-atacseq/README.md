# wasp2-nf-atacseq

ATAC-seq Allelic Imbalance (AI) Pipeline

## Features
- BWA-MEM2 alignment
- WASP2 mapping bias correction
- Peak-level allele counting
- caQTL integration support

## Quick Start
```bash
nextflow run . -profile singularity \
  --input samplesheet.csv \
  --vcf variants.vcf.gz \
  --fasta reference.fa \
  --peaks peaks.bed
```

## Output
- `results/counts/` - Allele counts at peaks
- `results/analysis/` - AI statistical results
- `results/bigwig/` - Signal tracks
