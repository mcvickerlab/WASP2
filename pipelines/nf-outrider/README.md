# nf-outrider: WASP2 + OUTRIDER Pipeline

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

**WASP2 + OUTRIDER for aberrant expression and mono-allelic expression detection.**

## Overview

nf-outrider integrates WASP2's high-performance allele counting with OUTRIDER's autoencoder-based outlier detection to identify:

1. **Aberrant Expression**: Gene-level expression outliers using OUTRIDER's neural network approach
2. **Mono-allelic Expression (MAE)**: Allele-specific expression with binomial statistics

## Workflow

```
RNA-seq BAMs → WASP2 Count → Gene Aggregation → OUTRIDER → Outlier Calls
                    ↓                              ↓
              61× faster than              Autoencoder-based
              GATK ASEReadCounter          outlier detection
                    ↓
              MAE Detection (binomial)
```

## Why WASP2 + OUTRIDER (vs nf-core/drop)?

| Feature | nf-core/drop | nf-outrider |
|---------|--------------|-------------|
| Allele counting | GATK ASEReadCounter (1600s) | **WASP2 (26s, 61× faster)** |
| Acceleration | None | **Rust-accelerated** |
| Statistics | Binomial | **Binomial with FDR correction** |
| MAE detection | Standard | **Enhanced with FDR-corrected binomial test** |
| Bias correction | Basic | **WASP mapping filter available** |

## Quick Start

```bash
# Run with Docker
nextflow run nf-outrider \
    -profile docker \
    --input samplesheet.csv \
    --vcf variants.vcf.gz \
    --gtf annotation.gtf \
    --outdir results

# Run with Singularity (HPC)
nextflow run nf-outrider \
    -profile singularity \
    --input samplesheet.csv \
    --vcf variants.vcf.gz \
    --gtf annotation.gtf
```

## Input

### Samplesheet (CSV)

```csv
sample,bam,bai
patient1,/path/to/patient1.bam,/path/to/patient1.bam.bai
patient2,/path/to/patient2.bam,/path/to/patient2.bam.bai
control1,/path/to/control1.bam,/path/to/control1.bam.bai
```

### Required Files

- **VCF**: Heterozygous variants (bgzip-compressed, tabix-indexed)
- **GTF**: Gene annotation (Gencode, Ensembl, or RefSeq format)
- **BAMs**: Aligned RNA-seq reads (coordinate-sorted, indexed)

## Output

```
results/
├── wasp2/
│   └── allele_counts/          # Per-sample, per-variant allele counts
├── aggregated/
│   └── *.gene_counts.tsv       # Gene-level expression for each sample
├── outrider/
│   ├── outrider_results.tsv    # Aberrant expression calls
│   ├── outrider_model.rds      # Trained OUTRIDER model
│   └── outrider_summary.html   # Interactive summary report
├── mae/
│   └── *.mae_results.tsv       # Mono-allelic expression calls
└── pipeline_info/
    └── execution_*.html        # Nextflow reports
```

## Parameters

### Core Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | required | Samplesheet CSV |
| `--vcf` | required | VCF with heterozygous variants |
| `--gtf` | required | Gene annotation GTF |
| `--outdir` | `./results` | Output directory |

### OUTRIDER Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outrider_padj` | 0.05 | Adjusted p-value cutoff |
| `--outrider_zScore` | 2 | Z-score cutoff |
| `--outrider_q` | auto | Encoding dimension |
| `--outrider_iterations` | 15 | Max iterations |
| `--outrider_convergence` | 1e-5 | Convergence threshold |

### MAE Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--skip_mae` | false | Skip MAE analysis |
| `--mae_min_count` | 10 | Min allele count |
| `--mae_padj` | 0.05 | P-value cutoff |
| `--mae_alt_ratio` | 0.8 | Alt ratio threshold |

## Profiles

```bash
# Docker (recommended for local)
-profile docker

# Singularity (recommended for HPC)
-profile singularity

# Conda
-profile conda

# Test with minimal data
-profile test
```

## OUTRIDER Algorithm

OUTRIDER uses an autoencoder neural network to:

1. **Learn latent representation**: Encode expression patterns into lower dimension
2. **Reconstruct expected counts**: Decode to predict "normal" expression
3. **Identify outliers**: Flag genes where observed >> expected

The encoding dimension `q` controls model complexity - auto-estimated by analyzing eigenvalue distribution.

## MAE Detection

Mono-allelic expression is detected using:

1. **Binomial test**: Statistical test for allelic imbalance
2. **Alt ratio threshold**: Identifies sites with extreme allelic ratios (default ≥0.8)
3. **Multiple testing correction**: Benjamini-Hochberg FDR control

## Citation

If you use nf-outrider, please cite:

```bibtex
@article{WASP2,
    title={WASP2: High-performance allele-specific analysis},
    author={...},
    journal={...},
    year={2024}
}

@article{OUTRIDER,
    title={OUTRIDER: A Statistical Method for Detecting Aberrantly Expressed Genes in RNA Sequencing Data},
    author={Brechtmann et al.},
    journal={Am J Hum Genet},
    year={2018},
    doi={10.1016/j.ajhg.2018.10.025}
}
```

## References

- [OUTRIDER Paper](https://doi.org/10.1016/j.ajhg.2018.10.025)
- [OUTRIDER GitHub](https://github.com/gagneurlab/OUTRIDER)
- [nf-core/drop](https://nf-co.re/drop/dev/) (reference implementation)
- [WASP2 Documentation](https://github.com/your-org/WASP2)

## License

MIT License - see [LICENSE](../../LICENSE) for details.

## Issue Tracking

- Issue: #35
