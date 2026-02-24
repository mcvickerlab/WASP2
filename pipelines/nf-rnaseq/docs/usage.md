# wasp2-nf-rnaseq: Usage

## Introduction

**wasp2-nf-rnaseq** is a Nextflow pipeline for RNA-seq Allele-Specific Expression (ASE) analysis using WASP2 for mapping bias correction.

The pipeline performs:
1. STAR alignment of RNA-seq reads
2. WASP2 mapping bias correction (remap-filter approach)
3. Allele counting at heterozygous SNPs
4. Statistical testing for allelic imbalance

## Quick Start

```bash
nextflow run pipelines/nf-rnaseq \
    --input samplesheet.csv \
    --vcf variants.vcf.gz \
    --star_index /path/to/star_index \
    --gtf genes.gtf \
    --outdir results \
    -profile docker
```

## Samplesheet Format

The pipeline requires a CSV samplesheet with the following columns:

| Column     | Description                                |
|------------|--------------------------------------------|
| `sample`   | Sample identifier (unique)                 |
| `fastq_1`  | Path to R1 FASTQ file (gzipped)           |
| `fastq_2`  | Path to R2 FASTQ file (optional, gzipped) |

### Example Samplesheet

```csv
sample,fastq_1,fastq_2
SAMPLE1,/data/sample1_R1.fastq.gz,/data/sample1_R2.fastq.gz
SAMPLE2,/data/sample2_R1.fastq.gz,/data/sample2_R2.fastq.gz
SAMPLE3,/data/sample3_R1.fastq.gz,
```

**Note:** Leave `fastq_2` empty for single-end data.

## Required Parameters

| Parameter     | Description                              |
|---------------|------------------------------------------|
| `--input`     | Path to samplesheet CSV                  |
| `--vcf`       | Path to VCF file with heterozygous SNPs (indexed) |
| `--star_index`| Path to STAR genome index directory      |

## Optional Parameters

### Reference Files

| Parameter | Description                               | Default |
|-----------|-------------------------------------------|---------|
| `--gtf`   | Path to GTF annotation file               | `null`  |

### WASP2 Options

| Parameter         | Description                          | Default |
|-------------------|--------------------------------------|---------|
| `--min_count`     | Minimum read count for testing       | `10`    |
| `--pseudocount`   | Pseudocount for ratio calculations   | `1`     |

### Output Options

| Parameter            | Description                      | Default    |
|----------------------|----------------------------------|------------|
| `--outdir`           | Output directory                 | `./results`|
| `--publish_dir_mode` | How to publish files             | `copy`     |

### Resource Limits

| Parameter      | Description              | Default    |
|----------------|--------------------------|------------|
| `--max_cpus`   | Maximum CPUs per process | `16`       |
| `--max_memory` | Maximum memory           | `128.GB`   |
| `--max_time`   | Maximum wall time        | `240.h`    |

## Running with Different Profiles

### Docker (recommended)

```bash
nextflow run pipelines/nf-rnaseq -profile docker \
    --input samplesheet.csv \
    --vcf variants.vcf.gz \
    --star_index /path/to/star_index
```

### Singularity (HPC environments)

```bash
nextflow run pipelines/nf-rnaseq -profile singularity \
    --input samplesheet.csv \
    --vcf variants.vcf.gz \
    --star_index /path/to/star_index
```

### Conda

```bash
nextflow run pipelines/nf-rnaseq -profile conda \
    --input samplesheet.csv \
    --vcf variants.vcf.gz \
    --star_index /path/to/star_index
```

## Test Profile

Run with minimal test data to validate installation:

```bash
nextflow run pipelines/nf-rnaseq -profile test,docker
```

For stub runs (no real data, validates workflow structure):

```bash
nextflow run pipelines/nf-rnaseq -profile test_stub,docker
```

## VCF Preparation

The VCF file should contain heterozygous SNPs for the samples being analyzed:

1. **Filter for heterozygous variants:**
   ```bash
   bcftools view -g het input.vcf.gz -Oz -o het_only.vcf.gz
   ```

2. **Index the VCF:**
   ```bash
   tabix -p vcf het_only.vcf.gz
   ```

3. **Optional - Filter for exonic regions:**
   ```bash
   bedtools intersect -a het_only.vcf.gz -b exons.bed -header | \
       bgzip > het_exonic.vcf.gz
   tabix -p vcf het_exonic.vcf.gz
   ```

## STAR Index Generation

If you don't have a STAR index:

```bash
STAR --runMode genomeGenerate \
    --runThreadN 8 \
    --genomeDir star_index \
    --genomeFastaFiles genome.fa \
    --sjdbGTFfile genes.gtf \
    --sjdbOverhang 100
```

## Pipeline Workflow

```
┌─────────────┐
│ Input FASTQ │
└──────┬──────┘
       │
       ▼
┌─────────────────┐
│ STAR Alignment  │ → Initial BAM
└──────┬──────────┘
       │
       ▼
┌─────────────────────┐
│ WASP2 make-reads    │ → Swapped FASTQs + keep.bam + to_remap.bam
└──────┬──────────────┘
       │
       ▼
┌─────────────────┐
│ STAR Re-align   │ → Remapped BAM
└──────┬──────────┘
       │
       ▼
┌──────────────────────┐
│ WASP2 filter-remapped│ → Bias-corrected BAM
└──────┬───────────────┘
       │
       ▼
┌─────────────────────┐
│ WASP2 count-variants│ → Allele counts TSV
└──────┬──────────────┘
       │
       ▼
┌───────────────────────┐
│ WASP2 find-imbalance  │ → Statistical results
└───────────────────────┘
```

## Troubleshooting

### Common Issues

**"VCF index not found"**
- Ensure your VCF is indexed with `tabix -p vcf file.vcf.gz`
- Both `.tbi` and `.csi` index formats are supported

**Out of memory errors**
- Reduce `--max_memory` or increase resources
- STAR alignment requires significant memory (~32GB for human genome)

**Missing samples in VCF**
- Ensure VCF sample names match samplesheet sample IDs
- Use `bcftools query -l variants.vcf.gz` to list samples

### Debug Mode

Enable debug output:

```bash
nextflow run pipelines/nf-rnaseq -profile debug,docker \
    --input samplesheet.csv \
    --vcf variants.vcf.gz \
    --star_index /path/to/star_index
```

## Citation

If you use this pipeline, please cite:

```
WASP2: Allele-specific analysis toolkit
https://github.com/mcvickerlab/WASP2
```
