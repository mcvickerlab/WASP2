# wasp2-nf-rnaseq: Output

This document describes the output produced by the pipeline.

## Pipeline Overview

The pipeline produces the following outputs:

```
results/
├── wasp_filtered/          # WASP2 bias-corrected BAM files
├── counts/                 # Allele counts at heterozygous SNPs
├── analysis/               # Statistical test results
└── pipeline_info/          # Execution reports and logs
```

## Output Directories

### `wasp_filtered/`

Contains WASP2 bias-corrected BAM files ready for downstream analysis.

| File | Description |
|------|-------------|
| `{sample}_wasp_filt.bam` | Bias-corrected aligned reads |
| `{sample}_wasp_filt.bam.bai` | BAM index |
| `{sample}.filter_stats.txt` | Filtering statistics |

**Filter Statistics Example:**
```
Sample: SAMPLE1
Total reads in remapped BAM: 1500000
Total reads in keep BAM: 8500000
Total reads after WASP filter: 9200000
```

The WASP filter removes reads that map to different locations when alleles are swapped, indicating potential mapping bias.

### `counts/`

Allele counts at heterozygous SNPs for each sample.

| File | Description |
|------|-------------|
| `{sample}_counts.tsv` | Tab-separated allele counts |

**Counts File Format:**

| Column | Description |
|--------|-------------|
| `chrom` | Chromosome |
| `pos` | Position (1-based) |
| `ref` | Reference allele |
| `alt` | Alternate allele |
| `region` | Gene/region ID (from GTF) |
| `ref_count` | Reference allele read count |
| `alt_count` | Alternate allele read count |
| `other_count` | Reads with neither allele |
| `N` | Total informative reads |

**Example:**
```tsv
chrom	pos	ref	alt	region	ref_count	alt_count	other_count	N
chr1	1250	A	G	ENSG00000001	25	18	2	43
chr1	2200	C	T	ENSG00000001	30	35	0	65
chr1	4500	G	A	ENSG00000002	12	8	1	20
```

### `analysis/`

Statistical test results for allelic imbalance.

| File | Description |
|------|-------------|
| `{sample}_ai_results.tsv` | Allelic imbalance test results |

**Results File Format:**

| Column | Description |
|--------|-------------|
| `region` | Gene/region identifier |
| `snp_count` | Number of heterozygous SNPs in region |
| `ref_sum` | Total reference allele counts |
| `alt_sum` | Total alternate allele counts |
| `mu` | Estimated mean allelic ratio |
| `null_ll` | Null hypothesis log-likelihood |
| `alt_ll` | Alternative hypothesis log-likelihood |
| `LRT` | Likelihood ratio test statistic |
| `pvalue` | P-value from chi-squared test |
| `fdr` | Benjamini-Hochberg FDR-adjusted p-value |

**Example:**
```tsv
region	snp_count	ref_sum	alt_sum	mu	null_ll	alt_ll	LRT	pvalue	fdr
ENSG00000001	5	125	180	0.41	-450.2	-445.8	8.8	0.003	0.015
ENSG00000002	3	45	48	0.48	-120.5	-120.4	0.2	0.65	0.85
ENSG00000003	8	210	95	0.69	-380.1	-360.2	39.8	2.7e-10	5.4e-09
```

**Interpreting Results:**
- **mu < 0.5**: Reference allele underexpressed
- **mu > 0.5**: Reference allele overexpressed
- **mu ≈ 0.5**: Balanced expression
- **FDR < 0.05**: Significant allelic imbalance

### `pipeline_info/`

Execution information and reports.

| File | Description |
|------|-------------|
| `execution_timeline_{timestamp}.html` | Visual timeline of process execution |
| `execution_report_{timestamp}.html` | Detailed execution report |
| `execution_trace_{timestamp}.txt` | Tab-separated execution metrics |
| `pipeline_dag_{timestamp}.html` | Pipeline DAG visualization |
| `software_versions.yml` | Software versions used |

**Software Versions Example:**
```yaml
STAR_ALIGN_INITIAL:
    star: 2.7.11a
    samtools: 1.18
WASP2_UNIFIED_MAKE_READS:
    wasp2: 1.2.0
WASP2_FILTER_REMAPPED:
    wasp2: 1.2.0
    samtools: 1.18
WASP2_COUNT_ALLELES:
    wasp2: 1.2.0
WASP2_ANALYZE_IMBALANCE:
    wasp2: 1.2.0
```

## Workflow Outputs

The pipeline emits the following Nextflow channels for integration with other workflows:

| Channel | Description |
|---------|-------------|
| `wasp_bam` | Tuple of (meta, bam, bai) for filtered BAMs |
| `counts` | Tuple of (meta, counts_tsv) for allele counts |
| `results` | Tuple of (meta, results_tsv) for AI results |
| `versions` | Collected software versions |

## File Sizes

Approximate output sizes per sample (human whole transcriptome):

| Output | Size |
|--------|------|
| WASP-filtered BAM | 2-5 GB |
| Allele counts | 1-10 MB |
| AI results | 100 KB - 1 MB |
| Pipeline reports | 1-5 MB total |

## Quality Control

### Filtering Statistics

Check the filter statistics to ensure reasonable filtering rates:

```bash
cat results/wasp_filtered/*filter_stats.txt
```

**Expected behavior:**
- Most reads should be in "keep" BAM (no overlapping variants)
- 5-15% of reads typically need remapping
- 80-95% of remapped reads should pass WASP filter

### Coverage Check

Verify sufficient coverage at heterozygous sites:

```bash
awk -F'\t' '$9 >= 10' results/counts/*_counts.tsv | wc -l
```

### Significant Imbalance

Count significant genes:

```bash
awk -F'\t' 'NR>1 && $10 < 0.05' results/analysis/*_ai_results.tsv | wc -l
```

## Downstream Analysis

### Combining Samples

Merge allele counts across samples for multi-sample analysis:

```bash
# Concatenate with sample ID
for f in results/counts/*_counts.tsv; do
    sample=$(basename $f _counts.tsv)
    awk -v s=$sample 'NR>1 {print s"\t"$0}' $f
done > combined_counts.tsv
```

### Visualization

Load results in R for visualization:

```r
library(tidyverse)

results <- read_tsv("results/analysis/SAMPLE1_ai_results.tsv")

# Volcano plot
ggplot(results, aes(x = log2(mu/(1-mu)), y = -log10(pvalue))) +
    geom_point(aes(color = fdr < 0.05)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    labs(x = "log2(Allelic Ratio)", y = "-log10(P-value)")
```

### Integration with eQTL Data

The allele counts can be integrated with eQTL analysis:

```bash
# Join with eQTL results by gene
join -t $'\t' -1 1 -2 1 \
    <(sort -k1 results/analysis/*_ai_results.tsv) \
    <(sort -k1 eqtl_results.tsv) \
    > integrated_results.tsv
```
