# wasp2-nf-rnaseq

RNA-seq Allele-Specific Expression (ASE) Pipeline

## Features
- STAR alignment with WASP mode
- WASP2 bias-corrected allele counting
- Beta-binomial statistical testing
- eQTL integration support

## Quick Start
```bash
nextflow run . -profile docker \
  --input samplesheet.csv \
  --vcf variants.vcf.gz \
  --fasta reference.fa \
  --gtf genes.gtf
```

## Samplesheet Format
```csv
sample,fastq_1,fastq_2
SAMPLE1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz
```

## Output
- `results/counts/` - Allele counts per sample
- `results/analysis/` - Statistical test results
- `results/multiqc/` - QC report
