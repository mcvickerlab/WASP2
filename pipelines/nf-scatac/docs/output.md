# nf-scatac: Output

## Introduction

This document describes the output files and directory structure produced by the nf-scatac pipeline for single-cell ATAC-seq allelic imbalance analysis.

## Pipeline Output

The pipeline outputs are organized in the following directory structure:

```
results/
├── allele_counts/           # Per-cell allele counts at het SNPs
│   └── {sample}_allele_counts.tsv
├── imbalance/               # Allelic imbalance analysis
│   └── {sample}_imbalance.tsv
├── variants/                # VCF-derived files
│   └── {sample}_het_snps.bed
└── pipeline_info/           # Execution reports
    ├── timeline.html
    ├── report.html
    └── trace.txt
```

## Output Files

### Allele Counts

**Directory**: `allele_counts/`

**File**: `{sample}_allele_counts.tsv`

Per-cell allele counts at heterozygous SNPs. Each row represents a cell-SNP combination:

| Column | Description |
|--------|-------------|
| cell_barcode | Cell barcode (from BAM CB tag) |
| chrom | Chromosome |
| pos | Position (1-based) |
| ref | Reference allele |
| alt | Alternate allele |
| ref_count | Reference allele read count for this cell |
| alt_count | Alternate allele read count for this cell |
| total_count | Total reads (ref + alt) for this cell |

### Allelic Imbalance Results

**Directory**: `imbalance/`

**File**: `{sample}_imbalance.tsv`

Statistical analysis of allelic imbalance at each SNP, aggregated across cells:

| Column | Description |
|--------|-------------|
| chrom | Chromosome |
| pos | Position (1-based) |
| ref | Reference allele |
| alt | Alternate allele |
| n_cells | Number of cells with coverage at this SNP |
| total_ref | Total reference reads across all cells |
| total_alt | Total alternate reads across all cells |
| total_count | Total reads (ref + alt) |
| pval | Beta-binomial p-value |
| fdr_pval | FDR-corrected p-value (Benjamini-Hochberg) |
| log2_ratio | log2(ref/alt) ratio |
| dispersion | Estimated overdispersion parameter |

### Variant Files

**Directory**: `variants/`

**File**: `{sample}_het_snps.bed`

BED file of heterozygous SNP positions extracted from the VCF:

| Column | Description |
|--------|-------------|
| chrom | Chromosome |
| start | Start position (0-based) |
| end | End position |
| name | SNP ID (chr:pos:ref>alt) |

### Pipeline Info

**Directory**: `pipeline_info/`

- `execution_report_*.html`: Nextflow execution report with process statistics
- `execution_timeline_*.html`: Timeline visualization of process execution
- `execution_trace_*.txt`: Detailed trace file for each process
- `pipeline_dag_*.html`: Pipeline DAG visualization

## Interpreting Results

### Allelic Imbalance Significance

Variants with significant allelic imbalance (AI) have:
- `fdr_pval < 0.05`: Statistically significant after FDR correction
- `|log2_ratio| > 0.5`: At least 1.4-fold difference between alleles
- `n_cells >= 10`: Sufficient cell coverage (recommended)

### Cell-Level vs Pseudo-bulk

The pipeline provides:
1. **Cell-level counts** (`allele_counts/`): For exploring cell-to-cell heterogeneity
2. **Pseudo-bulk analysis** (`imbalance/`): Aggregated statistical testing

### Quality Metrics

Check the execution report for:
- Number of cells processed
- Total variants analyzed
- Processing time per sample

## Downstream Analysis

### Loading Results in R

```r
library(readr)
library(dplyr)

# Load allelic imbalance results
ai_results <- read_tsv("results/imbalance/GM12878_imbalance.tsv")

# Filter significant variants
sig_ai <- ai_results %>%
  filter(fdr_pval < 0.05, abs(log2_ratio) > 0.5, n_cells >= 10)

# Load cell-level counts
cell_counts <- read_tsv("results/allele_counts/GM12878_allele_counts.tsv")

# Analyze per-cell heterogeneity at a specific variant
variant_cells <- cell_counts %>%
  filter(chrom == "chr1", pos == 12345)
```

### Loading Results in Python

```python
import pandas as pd

# Load counts
cell_counts = pd.read_csv(
    "results/allele_counts/GM12878_allele_counts.tsv",
    sep="\t"
)

# Load AI results
ai_results = pd.read_csv(
    "results/imbalance/GM12878_imbalance.tsv",
    sep="\t"
)

# Filter significant
sig_ai = ai_results[
    (ai_results['fdr_pval'] < 0.05) &
    (abs(ai_results['log2_ratio']) > 0.5) &
    (ai_results['n_cells'] >= 10)
]

# Cell-level analysis
pivot_counts = cell_counts.pivot_table(
    index='cell_barcode',
    columns=['chrom', 'pos'],
    values='ref_count',
    fill_value=0
)
```

### Integration with Scanpy/AnnData

```python
import scanpy as sc
import pandas as pd

# Load your scATAC AnnData
adata = sc.read_h5ad("scatac_data.h5ad")

# Load WASP2 cell counts
cell_counts = pd.read_csv(
    "results/allele_counts/GM12878_allele_counts.tsv",
    sep="\t"
)

# Add as additional layer or obsm
# ... depends on your analysis goals
```

## Troubleshooting

### Empty Counts File

- Check that your VCF contains heterozygous variants for the sample
- Ensure fragments file has the correct index (.tbi)
- Verify cell barcodes in fragments match expected format
- Check that the barcode_tag parameter matches your BAM format

### Low Cell Coverage

- Increase sequencing depth
- Check that cell barcodes are being parsed correctly
- Verify the chemistry parameter matches your library prep

### No Significant AI Results

- Increase sequencing depth per cell
- Pool related samples for more power
- Consider adjusting `--wasp_min_count` threshold
- Check that variants are truly heterozygous in your sample

### Memory Issues

For large datasets (>10,000 cells):
- Increase `--max_memory` parameter
- Consider running on a cluster with more resources
- Split input files by chromosome if necessary
