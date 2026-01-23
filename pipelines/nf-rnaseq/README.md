# wasp2-nf-rnaseq

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.0-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg)](https://sylabs.io/docs/)

RNA-seq Allele-Specific Expression (ASE) Pipeline using WASP2 for mapping bias correction.

## Features

- **STAR alignment** with two-pass mode for optimal splice-aware mapping
- **WASP2 bias correction** using the remap-filter approach
- **Rust-accelerated** allele counting for high performance
- **Beta-binomial statistical testing** for allelic imbalance detection
- **eQTL integration support** for combining ASE with population genetics

## Quick Start

```bash
# Minimal example
nextflow run pipelines/nf-rnaseq -profile docker \
    --input samplesheet.csv \
    --vcf variants.vcf.gz \
    --star_index /path/to/star_index

# With GTF annotation for gene-level analysis
nextflow run pipelines/nf-rnaseq -profile docker \
    --input samplesheet.csv \
    --vcf variants.vcf.gz \
    --star_index /path/to/star_index \
    --gtf genes.gtf \
    --outdir my_results
```

## Test the Pipeline

```bash
# Stub run (validates structure, no real data needed)
nextflow run pipelines/nf-rnaseq -profile test_stub,docker

# Full test with minimal data
nextflow run pipelines/nf-rnaseq -profile test,docker
```

## Samplesheet Format

```csv
sample,fastq_1,fastq_2
SAMPLE1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz
SAMPLE2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz
```

## Pipeline Overview

```
FASTQ → STAR align → WASP2 make-reads → STAR remap → WASP2 filter → count → analyze
```

1. **STAR Alignment**: Initial splice-aware alignment
2. **WASP2 make-reads**: Generate allele-swapped reads for bias detection
3. **STAR Remap**: Re-align swapped reads
4. **WASP2 filter**: Remove reads with mapping bias
5. **Count alleles**: Count reads at heterozygous SNPs
6. **Analyze**: Statistical testing for allelic imbalance

## Output

```
results/
├── wasp_filtered/      # Bias-corrected BAM files
│   ├── {sample}_wasp_filt.bam
│   └── {sample}.filter_stats.txt
├── counts/             # Allele counts at het SNPs
│   └── {sample}_counts.tsv
├── analysis/           # Statistical test results
│   └── {sample}_ai_results.tsv
└── pipeline_info/      # Execution reports
    ├── execution_timeline.html
    └── software_versions.yml
```

## Requirements

- Nextflow >= 22.10.0
- Docker, Singularity, or Conda
- STAR genome index
- Indexed VCF with heterozygous variants

## Example Commands

### Basic Analysis

```bash
nextflow run pipelines/nf-rnaseq -profile docker \
    --input samplesheet.csv \
    --vcf het_snps.vcf.gz \
    --star_index /ref/star_index \
    --gtf /ref/genes.gtf
```

### HPC with Singularity

```bash
nextflow run pipelines/nf-rnaseq -profile singularity \
    --input samplesheet.csv \
    --vcf het_snps.vcf.gz \
    --star_index /ref/star_index \
    --max_cpus 32 \
    --max_memory 128.GB
```

### Custom WASP2 Parameters

```bash
nextflow run pipelines/nf-rnaseq -profile docker \
    --input samplesheet.csv \
    --vcf het_snps.vcf.gz \
    --star_index /ref/star_index \
    --min_count 20 \
    --pseudocount 0.5
```

### Resume Failed Run

```bash
nextflow run pipelines/nf-rnaseq -profile docker \
    --input samplesheet.csv \
    --vcf het_snps.vcf.gz \
    --star_index /ref/star_index \
    -resume
```

## Documentation

- [Usage Guide](docs/usage.md) - Detailed parameter documentation
- [Output Description](docs/output.md) - Output file formats and interpretation

## Testing

Run nf-test suite:

```bash
cd pipelines/nf-rnaseq
nf-test test --tag pipeline
```

## Citation

If you use this pipeline, please cite:

```
WASP2: Allele-specific analysis toolkit
https://github.com/Jaureguy760/WASP2-exp
```

## License

MIT License - see repository root for details.
