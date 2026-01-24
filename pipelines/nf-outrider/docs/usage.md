# nf-outrider: Usage

## Introduction

**nf-outrider** is a Nextflow DSL2 pipeline that integrates WASP2's high-performance allele counting with OUTRIDER's autoencoder-based outlier detection. It identifies:

1. **Aberrant Expression**: Gene-level expression outliers using OUTRIDER's neural network approach
2. **Mono-allelic Expression (MAE)**: Allele-specific expression with binomial statistics

## Pipeline Summary

```
RNA-seq BAMs → WASP2 Count → Gene Aggregation → OUTRIDER → Outlier Calls
                    ↓                              ↓
              61× faster than              Autoencoder-based
              GATK ASEReadCounter          outlier detection
                    ↓
              MAE Detection (binomial)
```

## Quick Start

```bash
nextflow run nf-outrider \
    -profile docker \
    --input samplesheet.csv \
    --vcf variants.vcf.gz \
    --gtf annotation.gtf \
    --outdir results
```

## Samplesheet Input

The pipeline requires a samplesheet CSV file with aligned RNA-seq BAMs:

| Column | Required | Description |
|--------|----------|-------------|
| `sample` | Yes | Unique sample identifier |
| `bam` | Yes | Path to aligned BAM file |
| `bai` | No | Path to BAM index (auto-detected if not provided) |

### Example samplesheet:

```csv
sample,bam,bai
patient1,/data/patient1.bam,/data/patient1.bam.bai
patient2,/data/patient2.bam,/data/patient2.bam.bai
control1,/data/control1.bam,/data/control1.bam.bai
control2,/data/control2.bam,/data/control2.bam.bai
```

**Note**: OUTRIDER requires multiple samples (ideally 10+) for reliable outlier detection. The autoencoder learns "normal" expression patterns from the cohort.

## Required Parameters

| Parameter | Description |
|-----------|-------------|
| `--input` | Path to samplesheet CSV |
| `--vcf` | Phased VCF/BCF file with heterozygous variants |
| `--gtf` | Gene annotation GTF (Gencode, Ensembl, or RefSeq format) |

## Optional Parameters

### OUTRIDER Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outrider_padj` | 0.05 | Adjusted p-value cutoff for significance |
| `--outrider_zScore` | 2 | Z-score cutoff for outlier calls |
| `--outrider_q` | auto | Encoding dimension (auto-estimated if not set) |
| `--outrider_iterations` | 15 | Maximum iterations for model fitting |
| `--outrider_convergence` | 1e-5 | Convergence threshold |
| `--outrider_min_samples` | 3 | Minimum samples required |

### MAE Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--skip_mae` | false | Skip mono-allelic expression analysis |
| `--mae_min_count` | 10 | Minimum allele count for testing |
| `--mae_padj` | 0.05 | P-value cutoff for significance |
| `--mae_alt_ratio` | 0.8 | Alternative allele ratio threshold |

### WASP2 Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--wasp_min_count` | 10 | Minimum allele count for counting |
| `--gene_feature` | 'exon' | GTF feature type for gene aggregation |
| `--gene_attribute` | 'gene_id' | GTF attribute for gene ID |

### Output Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outdir` | './results' | Output directory |
| `--publish_dir_mode` | 'copy' | Publishing mode: 'copy', 'symlink', 'link' |

## Running with Profiles

### Docker (Recommended for local)

```bash
nextflow run nf-outrider -profile docker --input samplesheet.csv ...
```

### Singularity (Recommended for HPC)

```bash
nextflow run nf-outrider -profile singularity --input samplesheet.csv ...
```

### Conda

```bash
nextflow run nf-outrider -profile conda --input samplesheet.csv ...
```

### Test Profile

Run with minimal test data:

```bash
nextflow run nf-outrider -profile test,docker
```

## Example Commands

### Full Analysis

```bash
nextflow run nf-outrider \
    --input samplesheet.csv \
    --vcf phased_variants.vcf.gz \
    --gtf gencode.v38.annotation.gtf \
    --outrider_padj 0.05 \
    --outrider_zScore 2 \
    --outdir results \
    -profile docker
```

### Skip MAE Analysis

```bash
nextflow run nf-outrider \
    --input samplesheet.csv \
    --vcf variants.vcf.gz \
    --gtf annotation.gtf \
    --skip_mae \
    --outdir results \
    -profile singularity
```

### Custom OUTRIDER Parameters

```bash
nextflow run nf-outrider \
    --input samplesheet.csv \
    --vcf variants.vcf.gz \
    --gtf annotation.gtf \
    --outrider_q 10 \
    --outrider_iterations 20 \
    --outrider_padj 0.01 \
    --outdir results \
    -profile docker
```

## OUTRIDER Algorithm

OUTRIDER uses an autoencoder neural network approach:

1. **Encoding**: Learn a low-dimensional representation of expression patterns
2. **Decoding**: Reconstruct expected expression from the encoding
3. **Outlier Detection**: Compare observed vs expected, flag significant deviations

The encoding dimension `q` controls model complexity:
- **Auto mode** (default): Estimated from eigenvalue analysis
- **Manual**: Set via `--outrider_q` parameter

## MAE Detection

Mono-allelic expression is detected using:

1. **Binomial Test**: Statistical test for allelic imbalance at each heterozygous site
2. **Alt Ratio Threshold**: Identifies extreme imbalance (default ≥0.8 alt ratio)
3. **FDR Correction**: Benjamini-Hochberg multiple testing correction

## Resource Requirements

Typical resource usage (30 samples):

| Process | CPUs | Memory | Time |
|---------|------|--------|------|
| WASP2 counting (per sample) | 4 | 8 GB | 5-10 min |
| Gene aggregation | 2 | 4 GB | 2-5 min |
| Count merging | 2 | 8 GB | 5 min |
| OUTRIDER fitting | 4 | 16 GB | 15-30 min |
| MAE detection (per sample) | 2 | 4 GB | 2-5 min |

## Troubleshooting

### Common Issues

1. **Missing VCF index**: Ensure VCF is bgzipped and tabix-indexed
   ```bash
   bgzip variants.vcf
   tabix -p vcf variants.vcf.gz
   ```

2. **GTF format issues**: Ensure GTF has required attributes (gene_id, exon feature)

3. **Too few samples**: OUTRIDER needs multiple samples for reliable detection (10+ recommended)

4. **Memory errors during OUTRIDER**: Increase `--max_memory` or reduce sample count

### Resume Failed Runs

```bash
nextflow run nf-outrider ... -resume
```

## Citation

If you use nf-outrider, please cite:

- **WASP2**: High-performance allele-specific analysis
- **OUTRIDER**: Brechtmann et al. "OUTRIDER: A Statistical Method for Detecting Aberrantly Expressed Genes in RNA Sequencing Data." Am J Hum Genet (2018). doi:10.1016/j.ajhg.2018.10.025
