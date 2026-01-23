# nf-atacseq: Output

## Introduction

This document describes the output files and directory structure produced by the nf-atacseq pipeline.

## Pipeline Output

The pipeline outputs are organized in the following directory structure:

```
results/
├── fastqc/                    # FastQC reports
├── fastp/                     # Trimming reports
├── alignment/                 # Aligned BAMs and statistics
│   ├── stats/                 # samtools stats/flagstat
│   └── picard_metrics/        # Duplicate metrics
├── peaks/                     # MACS2 peak calls
├── wasp2/                     # WASP2 outputs
│   ├── remap/                 # Intermediate remap files
│   └── filtered/              # WASP-filtered BAMs
├── counts/                    # Allele count tables
├── analysis/                  # Allelic imbalance results
├── multiqc/                   # MultiQC report
└── pipeline_info/             # Execution reports
```

## Output Files

### FastQC

- `*.html`: FastQC HTML report
- `*.zip`: FastQC data archive

### Fastp

- `*.json`: Trimming statistics in JSON format
- `*.html`: Trimming report
- `*.log`: Processing log

### Alignment

- `*.sorted.bam`: Coordinate-sorted BAM file
- `*.sorted.bam.bai`: BAM index
- `*.stats`: samtools stats output
- `*.flagstat`: samtools flagstat output
- `*.metrics.txt`: Picard MarkDuplicates metrics

### Peaks

- `*_peaks.narrowPeak`: MACS2 peak calls in narrowPeak format
- `*_peaks.xls`: MACS2 peak statistics
- `*_summits.bed`: Peak summit positions

### WASP2 Filtered BAMs

- `*_wasp_filt.bam`: WASP-filtered BAM (mapping bias corrected)
- `*_wasp_filt.bam.bai`: Index
- `*_wasp_stats.txt`: WASP filtering statistics

### Counts

**File**: `*_counts.tsv`

Allele counts at heterozygous SNPs within peaks:

| Column | Description |
|--------|-------------|
| chrom | Chromosome |
| pos | Position (1-based) |
| ref | Reference allele |
| alt | Alternate allele |
| GT | Genotype |
| region | Peak/region ID |
| ref_count | Reference allele read count |
| alt_count | Alternate allele read count |
| other_count | Other allele read count |

### Analysis

**File**: `*_ai_results.tsv`

Allelic imbalance statistical results:

| Column | Description |
|--------|-------------|
| region | Peak/region ID |
| ref_count | Total reference reads |
| alt_count | Total alternate reads |
| total_count | Total reads |
| pval | Beta-binomial p-value |
| fdr_pval | FDR-corrected p-value |
| log2_ratio | log2(ref/alt) ratio |
| dispersion | Estimated dispersion |

### MultiQC

- `multiqc_report.html`: Aggregated QC report
- `multiqc_data/`: MultiQC data files

### Pipeline Info

- `execution_report_*.html`: Nextflow execution report
- `execution_timeline_*.html`: Timeline visualization
- `execution_trace_*.txt`: Process trace file
- `pipeline_dag_*.html`: Pipeline DAG visualization

## Interpreting Results

### Allelic Imbalance

Regions with significant allelic imbalance (AI) have:
- `fdr_pval < 0.05`: Statistically significant after FDR correction
- `|log2_ratio| > 0.5`: At least 1.4-fold difference between alleles

### WASP Filtering Statistics

The WASP statistics file shows:
- Total reads processed
- Reads passing WASP filter
- Reads removed due to mapping bias

### Quality Metrics

Check the MultiQC report for:
- Read quality scores
- Adapter contamination
- Alignment rates
- Duplication rates
- Peak calling statistics

## Downstream Analysis

### Loading Results in R

```r
library(readr)

# Load allelic imbalance results
ai_results <- read_tsv("results/analysis/sample_ai_results.tsv")

# Filter significant regions
sig_ai <- ai_results %>%
  filter(fdr_pval < 0.05, abs(log2_ratio) > 0.5)
```

### Loading Results in Python

```python
import pandas as pd

# Load counts
counts = pd.read_csv("results/counts/sample_counts.tsv", sep="\t")

# Load AI results
ai_results = pd.read_csv("results/analysis/sample_ai_results.tsv", sep="\t")

# Filter significant
sig_ai = ai_results[
    (ai_results['fdr_pval'] < 0.05) &
    (abs(ai_results['log2_ratio']) > 0.5)
]
```

## Troubleshooting

### Empty Counts File

- Check that your VCF contains heterozygous variants for the sample
- Ensure peaks overlap with variants
- Verify BAM has reads at variant positions

### No Significant AI Results

- Increase sequencing depth
- Check variant calling quality
- Consider adjusting `--wasp_min_count`
