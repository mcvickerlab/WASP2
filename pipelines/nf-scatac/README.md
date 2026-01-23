# nf-scatac

Single-Cell ATAC-seq Allelic Imbalance Pipeline

## Features

- 10x Genomics scATAC fragment support
- Allelic imbalance analysis at heterozygous SNPs
- Cell barcode propagation through pipeline stages
- Pseudo-bulk aggregation per sample
- Support for CellRanger ATAC output and custom fragments
- nf-core compliant subworkflow architecture (Issue #57)

## Architecture

```
nf-scatac/
├── main.nf                          # Entry point
├── nextflow.config                  # Pipeline configuration
├── workflows/
│   └── scatac.nf                    # Main workflow
├── subworkflows/
│   ├── local/
│   │   ├── utils_nfscatac_pipeline.nf  # Pipeline utilities
│   │   ├── wasp_allelic_sc/         # WASP2 single-cell integration
│   │   └── generate_fragments/      # Fragment file generation
│   └── nf-core/
│       └── bam_stats_samtools/      # BAM QC stats
├── modules/
│   └── local/
│       └── scatac_count_alleles/    # Per-cell allele counting
├── conf/
│   ├── base.config
│   ├── modules.config
│   ├── test_stub.config
│   └── test_real.config
└── tests/
    ├── main.nf.test
    └── subworkflows/
        ├── wasp_allelic_sc.nf.test
        └── generate_fragments.nf.test
```

## Quick Start

```bash
nextflow run . -profile docker \
  --input samplesheet.csv \
  --vcf variants.vcf.gz
```

## Samplesheet Format

| Column | Required | Description |
|--------|----------|-------------|
| sample | Yes | Sample identifier |
| fragments | Yes* | Path to fragments.tsv.gz |
| cellranger_dir | Yes* | Path to CellRanger ATAC output |
| barcode_tag | No | BAM tag for cell barcodes (default: CB) |
| chemistry | No | Library chemistry (default: 10x-atac-v2) |

*Either `fragments` or `cellranger_dir` is required

Example:
```csv
sample,fragments,cellranger_dir,barcode_tag,chemistry
GM12878_rep1,/path/to/fragments.tsv.gz,,CB,10x-atac-v2
GM12878_rep2,,/path/to/cellranger/output,CB,10x-atac-v2
```

## Single-Cell Meta Map

The pipeline propagates scATAC-specific metadata through all stages:

```groovy
[
    id: 'sample1',
    single_end: false,
    cell_barcode_tag: 'CB',
    umi_tag: null,  // ATAC typically doesn't have UMI
    chemistry: '10x-atac-v2'
]
```

## Subworkflows

### WASP_ALLELIC_SC

Single-cell WASP2 allelic imbalance analysis:

```groovy
include { WASP_ALLELIC_SC } from './subworkflows/local/wasp_allelic_sc/main'

WASP_ALLELIC_SC (
    ch_fragments,   // [ val(meta), path(fragments), path(tbi) ]
    ch_vcf          // [ val(meta), path(vcf), path(tbi) ]
)

// Outputs:
// - cell_counts: Per-cell allele counts at SNPs
// - imbalance: Allelic imbalance analysis results
```

### GENERATE_FRAGMENTS

Generate 10x-compatible fragments from BAM:

```groovy
include { GENERATE_FRAGMENTS } from './subworkflows/local/generate_fragments/main'

GENERATE_FRAGMENTS ( ch_bam )  // [ val(meta), path(bam), path(bai) ]

// Outputs:
// - fragments: [ val(meta), path(fragments.tsv.gz), path(tbi) ]
```

## Testing

### Stub Tests (CI/CD)

Run fast stub tests that validate workflow structure without real computation:

```bash
# Using nf-test
cd pipelines/nf-scatac
nf-test test --profile test_stub

# Or direct Nextflow stub run
nextflow run . -profile test_stub -stub-run
```

### Subworkflow Tests

```bash
# Test specific subworkflow
nf-test test tests/subworkflows/wasp_allelic_sc.nf.test
nf-test test tests/subworkflows/generate_fragments.nf.test
```

### Integration Tests (Real Data)

Run full pipeline with GM12878 scATAC-seq data:

```bash
nextflow run . -profile test_real,singularity
```

Test data locations:
- **BAM**: `/iblm/netapp/data3/aho/project_data/wasp2/10x_cellranger_atac/gm12878_el4/`
- **VCF**: `/iblm/netapp/data1/aho/variants/NA12878.vcf.gz`

### Test Structure

```
tests/
├── main.nf.test              # Pipeline & workflow tests
├── tags.yml                  # Test tags for filtering
├── stub/                     # Minimal stub test data
│   ├── samplesheet.csv
│   ├── fragments.tsv.gz
│   ├── variants.vcf.gz
│   └── variants.bed
├── real/                     # Integration test samplesheet
│   └── samplesheet.csv
└── subworkflows/             # Subworkflow-level tests
    ├── wasp_allelic_sc.nf.test
    └── generate_fragments.nf.test
```

## Output

```
results/
├── allele_counts/           # Per-cell allele counts at het SNPs
│   └── {sample}_allele_counts.tsv
├── imbalance/               # Allelic imbalance analysis
│   └── {sample}_imbalance.tsv
└── pipeline_info/           # Execution reports
    ├── timeline.html
    ├── report.html
    └── trace.txt
```

## Supported Chemistries

| Chemistry | Description |
|-----------|-------------|
| 10x-atac-v1 | 10x Genomics Single Cell ATAC v1 |
| 10x-atac-v2 | 10x Genomics Single Cell ATAC v2 (default) |
| custom | Custom scATAC-seq library prep |

## References

- Issue [#57](https://github.com/Jaureguy760/WASP2-final/issues/57) - nf-core Subworkflow Pattern Compliance
- Issue [#48](https://github.com/Jaureguy760/WASP2-final/issues/48) - Validation & Test Suite
- Issue [#32](https://github.com/Jaureguy760/WASP2-final/issues/32) - scATAC Pipeline
- [nf-test docs](https://code.askimed.com/nf-test/)
- [nf-core guidelines](https://nf-co.re/docs/guidelines)
