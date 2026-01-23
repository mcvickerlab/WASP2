# nf-outrider: OUTRIDER Integration Pipeline

WASP2 + OUTRIDER for aberrant expression and mono-allelic expression detection.

## Workflow

```
RNA-seq BAMs → WASP2 Filter → WASP2 Count → Gene Aggregation → OUTRIDER → Outlier Calls
                    ↓                              ↓
              Bias-corrected             Autoencoder-based
                  reads                  outlier detection
```

## Why WASP2 + OUTRIDER (vs nf-core/drop)?

| Feature | nf-core/drop | nf-outrider |
|---------|--------------|-------------|
| Allele counting | GATK ASEReadCounter | WASP2 (61× faster) |
| Bias correction | Basic | WASP2 mapping filter |
| Statistics | Binomial | Beta-binomial (overdispersion-aware) |
| MAE detection | Standard | Enhanced with WASP2 AI calls |

## Inputs

```
params {
    bams           = "data/*.bam"
    vcf            = "variants.vcf.gz"
    gtf            = "annotation.gtf"
    sample_sheet   = "samples.csv"
}
```

## Outputs

```
results/
├── wasp2/
│   ├── filtered_bams/      # Bias-corrected BAMs
│   └── allele_counts/      # Per-sample, per-variant counts
├── aggregated/
│   └── gene_counts.tsv     # Gene-level expression matrix
├── outrider/
│   ├── outliers.tsv        # Aberrant expression calls
│   └── model.rds           # Trained OUTRIDER model
└── mae/
    └── mae_results.tsv     # Mono-allelic expression
```

## Pipeline Modules

```
include { WASP2_FILTER } from '../nf-modules/wasp2_filter'
include { WASP2_COUNT } from '../nf-modules/wasp2_count'
include { AGGREGATE_GENES } from './modules/aggregate'
include { OUTRIDER_FIT } from './modules/outrider'
include { MAE_DETECT } from './modules/mae'
```

## OUTRIDER Integration

```groovy
process OUTRIDER_FIT {
    container 'ghcr.io/gagneurlab/outrider:latest'

    input:
    path(count_matrix)

    output:
    path("outrider_results.tsv"), emit: outliers
    path("outrider_model.rds"), emit: model

    script:
    """
    Rscript -e "
    library(OUTRIDER)
    ods <- OutriderDataSet(countData = as.matrix(read.csv('${count_matrix}')))
    ods <- OUTRIDER(ods)
    res <- results(ods, padjCutoff = 0.05)
    write.csv(res, 'outrider_results.tsv')
    saveRDS(ods, 'outrider_model.rds')
    "
    """
}
```

## References

- [OUTRIDER Paper](https://doi.org/10.1016/j.ajhg.2018.10.025)
- [OUTRIDER GitHub](https://github.com/gagneurlab/OUTRIDER)
- [nf-core/drop](https://nf-co.re/drop/dev/) (reference implementation)
- Issue: #35
