# nf-scatac: Usage

## Introduction

**nf-scatac** is a Nextflow DSL2 pipeline for single-cell ATAC-seq allelic imbalance (AI) analysis. It processes 10x Genomics scATAC-seq data (fragments files or CellRanger ATAC output) and uses WASP2 for per-cell allele counting and statistical testing to identify regions with significant allelic imbalance.

## Pipeline Summary

1. Input validation (samplesheet, VCF)
2. VCF to BED conversion for SNP positions
3. Per-cell allele counting at heterozygous SNPs
4. Allelic imbalance statistical analysis
5. Results aggregation

## Quick Start

```bash
nextflow run nf-scatac \
    --input samplesheet.csv \
    --vcf variants.vcf.gz \
    --outdir results \
    -profile docker
```

## Samplesheet Input

The pipeline requires a samplesheet CSV file with the following columns:

| Column | Required | Description |
|--------|----------|-------------|
| `sample` | Yes | Unique sample identifier |
| `fragments` | Yes* | Path to 10x fragments.tsv.gz file |
| `cellranger_dir` | Yes* | Path to CellRanger ATAC output directory |
| `barcode_tag` | No | BAM tag for cell barcodes (default: CB) |
| `chemistry` | No | Library chemistry (default: 10x-atac-v2) |

*Either `fragments` or `cellranger_dir` is required.

### Example samplesheet (fragments):

```csv
sample,fragments,cellranger_dir,barcode_tag,chemistry
GM12878_rep1,/data/GM12878/fragments.tsv.gz,,CB,10x-atac-v2
GM12878_rep2,/data/GM12878_rep2/fragments.tsv.gz,,CB,10x-atac-v2
```

### Example samplesheet (CellRanger ATAC output):

```csv
sample,fragments,cellranger_dir,barcode_tag,chemistry
GM12878_cellranger,,/data/cellranger_atac/GM12878,CB,10x-atac-v2
```

## Required Parameters

| Parameter | Description |
|-----------|-------------|
| `--input` | Path to samplesheet CSV |
| `--vcf` | Phased VCF/BCF file with variants (must be bgzipped and indexed) |

## Optional Parameters

### WASP2 Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--wasp_min_count` | 10 | Minimum allele count for AI analysis |
| `--wasp_pseudocount` | 1 | Pseudocount for beta-binomial model |
| `--wasp_phased` | false | Use phased haplotype model |

### Single-Cell Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--barcode_tag` | 'CB' | Default BAM tag for cell barcodes |
| `--chemistry` | '10x-atac-v2' | Default library chemistry |

### Output Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outdir` | './results' | Output directory |
| `--publish_dir_mode` | 'copy' | Publishing mode: 'copy', 'symlink', 'link' |

## Running with Profiles

### Docker

```bash
nextflow run nf-scatac -profile docker --input samplesheet.csv ...
```

### Singularity

```bash
nextflow run nf-scatac -profile singularity --input samplesheet.csv ...
```

### Stub Test (CI/CD)

Run fast stub tests that validate workflow structure:

```bash
nextflow run nf-scatac -profile test_stub -stub-run
```

### Integration Test

Run with real GM12878 scATAC-seq data:

```bash
nextflow run nf-scatac -profile test_real,singularity
```

## Example Commands

### Full Single-Cell Analysis

```bash
nextflow run nf-scatac \
    --input samplesheet.csv \
    --vcf phased_variants.vcf.gz \
    --wasp_min_count 10 \
    --outdir results \
    -profile docker
```

### Using Phased Haplotype Model

```bash
nextflow run nf-scatac \
    --input samplesheet.csv \
    --vcf phased_variants.vcf.gz \
    --wasp_phased true \
    --outdir results \
    -profile singularity
```

## Supported Chemistries

| Chemistry | Description |
|-----------|-------------|
| `10x-atac-v1` | 10x Genomics Single Cell ATAC v1 |
| `10x-atac-v2` | 10x Genomics Single Cell ATAC v2 (default) |
| `custom` | Custom scATAC-seq library prep |

## Resource Requirements

Typical resource usage per sample (5000 cells):

| Process | CPUs | Memory | Time |
|---------|------|--------|------|
| VCF to BED conversion | 2 | 4 GB | 2-5 min |
| Per-cell allele counting | 4 | 8 GB | 15-30 min |
| Imbalance analysis | 2 | 4 GB | 5-10 min |

## Troubleshooting

### Common Issues

1. **Missing VCF index**: Ensure your VCF is bgzipped and indexed with tabix
   ```bash
   bgzip variants.vcf
   tabix -p vcf variants.vcf.gz
   ```

2. **Fragments file not indexed**: Index with tabix
   ```bash
   tabix -p bed fragments.tsv.gz
   ```

3. **Memory errors**: Increase `--max_memory` or use a profile with more resources

4. **No cells found**: Check that barcode_tag matches your BAM (default: CB)

### Resume Failed Runs

```bash
nextflow run nf-scatac ... -resume
```

## Citation

If you use nf-scatac, please cite:

- WASP2: [GitHub](https://github.com/your-org/WASP2)
- Nextflow: [Nextflow](https://www.nextflow.io/)
