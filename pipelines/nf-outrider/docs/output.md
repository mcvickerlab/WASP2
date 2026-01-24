# nf-outrider: Output

## Introduction

This document describes the output files and directory structure produced by the nf-outrider pipeline for aberrant expression and mono-allelic expression detection.

## Pipeline Output

The pipeline outputs are organized in the following directory structure:

```
results/
├── wasp2/
│   └── allele_counts/           # Per-sample, per-variant allele counts
│       └── {sample}_counts.tsv
├── aggregated/
│   └── {sample}_gene_counts.tsv # Gene-level aggregated counts
├── merged/
│   └── gene_count_matrix.tsv    # Sample × Gene count matrix
├── outrider/
│   ├── outrider_results.tsv     # Aberrant expression calls
│   ├── outrider_model.rds       # Trained OUTRIDER model
│   └── outrider_summary.html    # Interactive summary report
├── mae/
│   └── {sample}_mae_results.tsv # Mono-allelic expression calls
└── pipeline_info/
    ├── timeline.html
    ├── report.html
    └── trace.txt
```

## Output Files

### WASP2 Allele Counts

**Directory**: `wasp2/allele_counts/`

**File**: `{sample}_counts.tsv`

Per-variant allele counts for each sample:

| Column | Description |
|--------|-------------|
| chrom | Chromosome |
| pos | Position (1-based) |
| ref | Reference allele |
| alt | Alternate allele |
| GT | Sample genotype |
| ref_count | Reference allele read count |
| alt_count | Alternate allele read count |
| total_count | Total reads (ref + alt) |

### Gene-Level Aggregated Counts

**Directory**: `aggregated/`

**File**: `{sample}_gene_counts.tsv`

Allele counts aggregated to gene level:

| Column | Description |
|--------|-------------|
| gene_id | Gene identifier (from GTF) |
| gene_name | Gene symbol (if available) |
| n_variants | Number of heterozygous variants in gene |
| ref_count | Total reference reads across gene |
| alt_count | Total alternate reads across gene |
| total_count | Total reads |
| mean_ref_ratio | Mean reference allele ratio |

### Merged Count Matrix

**Directory**: `merged/`

**File**: `gene_count_matrix.tsv`

Combined gene × sample matrix for OUTRIDER input:

| Column | Description |
|--------|-------------|
| gene_id | Gene identifier |
| sample1 | Total counts for sample1 |
| sample2 | Total counts for sample2 |
| ... | Additional samples |

### OUTRIDER Results

**Directory**: `outrider/`

#### Main Results: `outrider_results.tsv`

Aberrant expression outlier calls:

| Column | Description |
|--------|-------------|
| sampleID | Sample identifier |
| geneID | Gene identifier |
| observed | Observed expression (normalized counts) |
| expected | Expected expression from autoencoder |
| pValue | Raw p-value |
| padjust | FDR-adjusted p-value |
| zScore | Z-score of deviation |
| l2fc | Log2 fold change (observed/expected) |
| aberrant | Boolean: TRUE if significant outlier |

#### Model File: `outrider_model.rds`

Serialized OUTRIDER model (R object) containing:
- Trained autoencoder weights
- Normalization factors
- Sample/gene metadata
- Can be loaded in R for further analysis

#### Summary Report: `outrider_summary.html`

Interactive HTML report with:
- PCA plot of samples
- Outlier heatmap
- Sample-level aberrant gene counts
- Gene-level aberration frequencies

### MAE Results

**Directory**: `mae/`

**File**: `{sample}_mae_results.tsv`

Mono-allelic expression calls per sample:

| Column | Description |
|--------|-------------|
| gene_id | Gene identifier |
| chrom | Chromosome |
| pos | Position |
| ref | Reference allele |
| alt | Alternate allele |
| ref_count | Reference allele reads |
| alt_count | Alternate allele reads |
| total_count | Total reads |
| alt_ratio | Alternate allele ratio |
| pval | Binomial test p-value |
| padj | FDR-adjusted p-value |
| mae_call | Boolean: TRUE if significant MAE |

### Pipeline Info

**Directory**: `pipeline_info/`

- `execution_report_*.html`: Nextflow execution report
- `execution_timeline_*.html`: Timeline visualization
- `execution_trace_*.txt`: Process trace file

## Interpreting Results

### Aberrant Expression (OUTRIDER)

Genes with aberrant expression have:
- `aberrant == TRUE`: Flagged as significant outlier
- `padjust < 0.05`: Statistically significant after FDR correction
- `|zScore| > 2`: Substantial deviation from expected

**Direction of effect**:
- `l2fc > 0`: Over-expression (observed > expected)
- `l2fc < 0`: Under-expression (observed < expected)

### Mono-allelic Expression

Sites with MAE have:
- `mae_call == TRUE`: Significant allelic imbalance
- `padj < 0.05`: Statistically significant
- `alt_ratio >= 0.8` or `alt_ratio <= 0.2`: Extreme allelic ratio

## Downstream Analysis

### Loading OUTRIDER Results in R

```r
library(readr)
library(OUTRIDER)

# Load outlier calls
outliers <- read_tsv("results/outrider/outrider_results.tsv")

# Filter significant aberrations
sig_outliers <- outliers %>%
  filter(aberrant == TRUE, padjust < 0.05)

# Load trained model for further analysis
ods <- readRDS("results/outrider/outrider_model.rds")

# Plot specific gene
plotExpressedGenes(ods)
plotQQ(ods, gene = "ENSG00000123456")
```

### Loading Results in Python

```python
import pandas as pd

# Load OUTRIDER results
outliers = pd.read_csv(
    "results/outrider/outrider_results.tsv",
    sep="\t"
)

# Filter significant
sig_outliers = outliers[
    (outliers['aberrant'] == True) &
    (outliers['padjust'] < 0.05)
]

# Load MAE results
mae_results = pd.read_csv(
    "results/mae/patient1_mae_results.tsv",
    sep="\t"
)

sig_mae = mae_results[mae_results['mae_call'] == True]
```

### Integrating with Variant Data

```r
library(VariantAnnotation)

# Load VCF
vcf <- readVcf("variants.vcf.gz")

# Join with MAE results
mae <- read_tsv("results/mae/patient1_mae_results.tsv")

# Annotate with variant info
mae_annotated <- mae %>%
  mutate(variant_id = paste(chrom, pos, sep = ":")) %>%
  left_join(variant_annotations, by = "variant_id")
```

## Quality Control

### OUTRIDER Diagnostics

Check the HTML summary report for:
- **PCA plot**: Samples should cluster by expected groups (batch, condition)
- **Encoding dimension**: Auto-estimated q should be reasonable (typically 5-20)
- **Convergence**: Model should converge within iterations

### Sample Quality

Flag samples with:
- Very high outlier counts (potential technical issues)
- Unusual PCA positions (batch effects, contamination)
- Low variant coverage

## Troubleshooting

### Few or No Outliers Detected

- Check sample size (need 10+ samples for reliable detection)
- Verify GTF annotation matches genome build
- Consider adjusting p-value threshold

### Too Many Outliers

- Check for batch effects in PCA
- Verify sample quality
- Consider increasing z-score threshold

### MAE Results Empty

- Check VCF contains heterozygous variants for samples
- Verify variants overlap gene annotations
- Increase coverage by adjusting `--mae_min_count`
