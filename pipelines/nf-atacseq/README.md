# nf-atacseq

[![nf-atacseq Tests](https://github.com/your-org/WASP2/actions/workflows/nf-atacseq-tests.yml/badge.svg)](https://github.com/your-org/WASP2/actions/workflows/nf-atacseq-tests.yml)

ATAC-seq Allelic Imbalance (AI) Pipeline with WASP2 mapping bias correction.

## Overview

**nf-atacseq** is a Nextflow DSL2 pipeline that performs allelic imbalance analysis on ATAC-seq data. It integrates WASP2 for mapping bias correction, ensuring accurate quantification of chromatin accessibility at heterozygous sites.

## Features

- **Dual aligner support**: BWA-MEM or Bowtie2
- **WASP2 mapping bias correction**: Eliminates reference bias
- **Peak calling**: MACS2 or use pre-called peaks
- **Allele counting**: Count reads at heterozygous SNPs within peaks
- **Statistical testing**: Beta-binomial model with FDR correction
- **Comprehensive QC**: FastQC, fastp, MultiQC reports

## Quick Start

### Minimal Example

```bash
nextflow run pipelines/nf-atacseq \
    --input samplesheet.csv \
    --vcf variants.vcf.gz \
    --fasta genome.fa \
    -profile docker
```

### With Pre-built Index

```bash
nextflow run pipelines/nf-atacseq \
    --input samplesheet.csv \
    --vcf phased_variants.vcf.gz \
    --fasta hg38.fa \
    --bwa_index /ref/bwa_index \
    --outdir results \
    -profile singularity
```

### Using Pre-called Peaks

```bash
nextflow run pipelines/nf-atacseq \
    --input samplesheet.csv \
    --vcf variants.vcf.gz \
    --fasta hg38.fa \
    --peaks consensus_peaks.bed \
    --skip_peak_calling \
    -profile docker
```

### Test Run

```bash
nextflow run pipelines/nf-atacseq -profile test,docker
nextflow run pipelines/nf-atacseq -profile test,docker -stub-run  # Workflow validation only
```

## Samplesheet Format

```csv
sample,fastq_1,fastq_2,sample_name
ATAC_sample1,/data/sample1_R1.fastq.gz,/data/sample1_R2.fastq.gz,NA12878
ATAC_sample2,/data/sample2_R1.fastq.gz,/data/sample2_R2.fastq.gz,HG00096
```

| Column | Description |
|--------|-------------|
| `sample` | Unique sample identifier |
| `fastq_1` | Path to R1 FASTQ file |
| `fastq_2` | Path to R2 FASTQ file (optional for single-end) |
| `sample_name` | Sample name in VCF (for het variant filtering) |

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | required | Samplesheet CSV path |
| `--vcf` | required | VCF/BCF with variants |
| `--fasta` | required | Reference genome FASTA |
| `--aligner` | 'bwa' | Aligner: 'bwa' or 'bowtie2' |
| `--peaks` | null | Pre-called peaks BED file |
| `--skip_wasp` | false | Skip WASP bias correction |
| `--wasp_min_count` | 10 | Min reads for AI testing |
| `--outdir` | './results' | Output directory |

See [docs/usage.md](docs/usage.md) for complete parameter reference.

## Output

```
results/
├── fastqc/          # Raw read QC
├── alignment/       # BAMs, stats, dup metrics
├── peaks/           # MACS2 narrowPeak files
├── wasp2/           # WASP-filtered BAMs
├── counts/          # Allele count tables
├── analysis/        # AI statistical results
├── multiqc/         # Aggregated QC report
└── pipeline_info/   # Execution reports
```

### Key Output Files

- **`*_counts.tsv`**: Per-SNP allele counts at peaks
- **`*_ai_results.tsv`**: Statistical test results with FDR p-values
- **`*_wasp_filt.bam`**: WASP-corrected BAM files

See [docs/output.md](docs/output.md) for detailed output descriptions.

## Testing

### Run nf-test Suite

```bash
cd pipelines/nf-atacseq

# Install nf-test
curl -fsSL https://code.askimed.com/install/nf-test | bash

# Run all tests
nf-test test

# Run stub tests only (fast)
nf-test test --tag ci_stub

# Run specific module tests
nf-test test --tag wasp2
```

### Manual Stub Run

Validate workflow structure without data:

```bash
nextflow run . -profile test -stub-run
```

## Profiles

| Profile | Description |
|---------|-------------|
| `docker` | Run with Docker containers |
| `singularity` | Run with Singularity containers |
| `conda` | Run with Conda environments |
| `test` | Minimal test configuration |
| `test_full` | Full test with real data |

## Pipeline DAG

```
FASTQ → FastQC → Fastp → BWA/Bowtie2 → Samtools → Picard → MACS2 → WASP2 → Counts → AI Analysis
                                                                    ↓
                                                              MultiQC Report
```

## Requirements

- Nextflow >= 23.04.0
- Java 11+
- Docker, Singularity, or Conda

## Citation

If you use nf-atacseq, please cite:

- **WASP2**: [GitHub Repository](https://github.com/your-org/WASP2)
- **Nextflow**: Di Tommaso, P., et al. (2017). Nextflow enables reproducible computational workflows. *Nature Biotechnology*.

## License

MIT License - see [LICENSE](../../LICENSE) for details.

## Support

- [Issues](https://github.com/your-org/WASP2/issues)
- [Documentation](docs/)
