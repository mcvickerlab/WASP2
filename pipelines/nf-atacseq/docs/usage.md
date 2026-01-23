# nf-atacseq: Usage

## Introduction

**nf-atacseq** is a Nextflow DSL2 pipeline for ATAC-seq allelic imbalance (AI) analysis. It uses WASP2 for mapping bias correction and performs beta-binomial statistical testing to identify regions with significant allelic imbalance.

## Pipeline Summary

1. Read QC ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Adapter trimming ([fastp](https://github.com/OpenGene/fastp))
3. Alignment ([BWA-MEM](https://github.com/lh3/bwa) or [Bowtie2](https://github.com/BenLangmead/bowtie2))
4. Duplicate marking ([Picard MarkDuplicates](https://broadinstitute.github.io/picard/))
5. Peak calling ([MACS2](https://github.com/macs3-project/MACS))
6. WASP2 mapping bias correction
7. Allele counting at peaks
8. Allelic imbalance statistical analysis
9. MultiQC report

## Quick Start

```bash
nextflow run nf-atacseq \
    --input samplesheet.csv \
    --vcf variants.vcf.gz \
    --fasta genome.fa \
    --outdir results \
    -profile docker
```

## Samplesheet Input

The pipeline requires a samplesheet CSV file with the following columns:

| Column | Description |
|--------|-------------|
| `sample` | Unique sample identifier |
| `fastq_1` | Path to R1 FASTQ file |
| `fastq_2` | Path to R2 FASTQ file |
| `sample_name` | (Optional) Sample name in VCF for het filtering |

Example samplesheet:

```csv
sample,fastq_1,fastq_2,sample_name
ATAC_sample1,/data/sample1_R1.fastq.gz,/data/sample1_R2.fastq.gz,NA12878
ATAC_sample2,/data/sample2_R1.fastq.gz,/data/sample2_R2.fastq.gz,HG00096
```

## Required Parameters

| Parameter | Description |
|-----------|-------------|
| `--input` | Path to samplesheet CSV |
| `--vcf` | Phased VCF/BCF/PGEN file with variants |
| `--fasta` | Reference genome FASTA |

## Optional Parameters

### Reference Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--bwa_index` | null | Pre-built BWA index directory |
| `--bowtie2_index` | null | Pre-built Bowtie2 index directory |
| `--peaks` | null | Pre-called peaks BED file |

### WASP2 Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--wasp_min_count` | 10 | Minimum allele count for AI analysis |
| `--wasp_pseudocount` | 1 | Pseudocount for beta-binomial model |
| `--wasp_phased` | false | Use phased haplotype model |
| `--wasp_include_indels` | false | Include indels in analysis |

### Processing Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--aligner` | 'bwa' | Aligner to use: 'bwa' or 'bowtie2' |
| `--macs_gsize` | 'hs' | MACS2 effective genome size |
| `--skip_trimming` | false | Skip adapter trimming |
| `--skip_dedup` | false | Skip duplicate marking |
| `--skip_wasp` | false | Skip WASP filtering |
| `--skip_peak_calling` | false | Skip peak calling (requires --peaks) |

### Output Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outdir` | './results' | Output directory |
| `--publish_dir_mode` | 'copy' | Publishing mode: 'copy', 'symlink', 'link' |

## Running with Profiles

### Docker

```bash
nextflow run nf-atacseq -profile docker --input samplesheet.csv ...
```

### Singularity

```bash
nextflow run nf-atacseq -profile singularity --input samplesheet.csv ...
```

### Conda

```bash
nextflow run nf-atacseq -profile conda --input samplesheet.csv ...
```

### Test Profile

Run with minimal test data:

```bash
nextflow run nf-atacseq -profile test,docker
```

## Example Commands

### Full Analysis with WASP

```bash
nextflow run nf-atacseq \
    --input samplesheet.csv \
    --vcf phased_variants.vcf.gz \
    --fasta hg38.fa \
    --bwa_index /ref/bwa_index \
    --macs_gsize hs \
    --outdir results \
    -profile docker
```

### Using Pre-called Peaks

```bash
nextflow run nf-atacseq \
    --input samplesheet.csv \
    --vcf variants.vcf.gz \
    --fasta hg38.fa \
    --peaks consensus_peaks.bed \
    --skip_peak_calling \
    --outdir results \
    -profile singularity
```

### Skip WASP (Basic ATAC-seq)

```bash
nextflow run nf-atacseq \
    --input samplesheet.csv \
    --fasta hg38.fa \
    --skip_wasp \
    --outdir results \
    -profile conda
```

## Resource Requirements

Typical resource usage per sample (30M paired-end reads):

| Process | CPUs | Memory | Time |
|---------|------|--------|------|
| Alignment (BWA-MEM) | 8 | 16 GB | 30-45 min |
| WASP make-reads | 4 | 8 GB | 5-10 min |
| WASP remapping | 8 | 16 GB | 15-20 min |
| WASP filtering | 4 | 8 GB | 5 min |
| Peak calling (MACS2) | 4 | 8 GB | 10 min |
| Allele counting | 4 | 8 GB | 10 min |

## Troubleshooting

### Common Issues

1. **Missing VCF index**: Ensure your VCF is bgzipped and indexed with tabix
2. **Memory errors**: Increase `--max_memory` or use a profile with more resources
3. **No peaks found**: Check that MACS2 `--gsize` matches your genome

### Resume Failed Runs

```bash
nextflow run nf-atacseq ... -resume
```

## Citation

If you use nf-atacseq, please cite:

- WASP2: [GitHub](https://github.com/your-org/WASP2)
- Nextflow: [Nextflow](https://www.nextflow.io/)
