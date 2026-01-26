# nf-scatac

Single-Cell ATAC-seq Allelic Imbalance Pipeline

## Features

- **10x Genomics scATAC fragment support** - Direct input from CellRanger ATAC output
- **Allelic imbalance analysis** - At heterozygous SNPs using WASP2
- **Cell barcode filtering** - Optional whitelist for quality-filtered cells
- **Peak region filtering** - Restrict analysis to accessible regions
- **AnnData/H5AD output** - For scverse ecosystem (Scanpy, ArchR, Signac)
- **Zarr output** - For GenVarLoader integration
- **Pseudo-bulk aggregation** - Sample-level aggregation for statistical power
- **nf-core compliant** - Subworkflow architecture (Issue #57)

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
│       ├── scatac_add_haplotype_layers/  # Hap1/hap2 layer creation from phased VCF
│       ├── scatac_count_alleles/         # Per-cell allele counting (fragment-based)
│       ├── scatac_create_anndata/        # AnnData H5AD output
│       └── scatac_pseudobulk/            # Pseudo-bulk aggregation
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
| sample | Yes | Sample identifier (must match VCF sample name for BAM mode) |
| fragments | Yes* | Path to fragments.tsv.gz |
| cellranger_dir | Yes* | Path to CellRanger ATAC output |
| bam | No | Path to BAM file (enables allele-specific counting with ref/alt/hap1/hap2 layers) |
| barcode_tag | No | BAM tag for cell barcodes (default: CB) |
| chemistry | No | Library chemistry (default: 10x-atac-v2) |
| barcodes | No | File with valid cell barcodes (one per line) |
| peaks | No | BED file with peak regions to restrict analysis |

*Either `fragments` or `cellranger_dir` is required for fragment-based counting

**Note on counting modes:**
- **Fragment-based** (fragments only): Counts fragment overlaps at SNP positions. Output has only `X` layer (total overlaps).
- **BAM-based** (bam column provided): True allele-specific counting. Output has `X`, `ref`, `alt`, `hap1`, `hap2` layers.

Example:
```csv
sample,fragments,cellranger_dir,bam,barcode_tag,chemistry,barcodes,peaks
GM12878_rep1,/path/to/fragments.tsv.gz,,,CB,10x-atac-v2,/path/to/barcodes.txt,/path/to/peaks.bed
GM12878_rep2,,/path/to/cellranger/output,,CB,10x-atac-v2,,
NA12878_bam,,,/path/to/possorted_bam.bam,CB,10x-atac-v2,/path/to/barcodes.txt,
```

## Parameters

### Required
- `--input` - Samplesheet CSV (see format above)
- `--vcf` - Indexed VCF/BCF with heterozygous SNPs

### Processing Options
- `--min_fragments_per_cell` - Minimum fragments per cell to include [default: 1000]
- `--min_cells_per_snp` - Minimum cells per SNP for pseudo-bulk [default: 3]
- `--min_count` - Minimum count for imbalance testing [default: 10]

### Output Options
- `--outdir` - Output directory [default: ./results]
- `--create_zarr` - Also output Zarr format for GenVarLoader [default: false]
- `--skip_anndata` - Skip AnnData H5AD creation [default: false]
- `--skip_pseudobulk` - Skip pseudo-bulk aggregation and analysis [default: false]

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
    ch_fragments,   // [ val(meta), path(fragments), path(tbi), path(barcodes), path(peaks) ]
    ch_vcf          // [ val(meta), path(vcf), path(tbi) ]
)

// Outputs:
// - cell_counts: Per-cell allele counts at SNPs
// - anndata: AnnData H5AD files
// - zarr: Zarr directories (if enabled)
// - cell_qc: Cell QC metrics
// - pseudobulk: Aggregated counts
// - imbalance: Allelic imbalance results
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

## Output

```
results/
├── allele_counts/           # Per-cell allele counts at het SNPs
│   └── {sample}_allele_counts.tsv
├── count_stats/             # Counting statistics
│   └── {sample}_count_stats.tsv
├── anndata/                 # AnnData H5AD files for scverse
│   └── {sample}_allelic.h5ad
├── zarr/                    # Zarr directories (if --create_zarr)
│   └── {sample}_allelic.zarr/
├── cell_qc/                 # Cell QC metrics
│   └── {sample}_cell_qc.tsv
├── pseudobulk/              # Pseudo-bulk aggregated counts
│   ├── {sample}_pseudobulk_counts.tsv
│   └── {sample}_aggregation_stats.tsv
├── imbalance/               # Allelic imbalance analysis
│   └── {sample}_ai_results.tsv
├── variants/                # Processed variant BED
│   └── variants.variants.bed
└── pipeline_info/           # Execution reports
    ├── timeline.html
    ├── report.html
    └── trace.txt
```

## AnnData Output Format

The H5AD file contains different layers depending on input type:

### Fragment-based input (overlap counting)
- **X**: Sparse matrix of total fragment overlaps (cells × SNPs)

### BAM-based input (allele-specific counting)
- **X**: Sparse matrix of total counts (cells × SNPs)
- **layers**:
  - `ref`: Reference allele counts per cell/SNP
  - `alt`: Alternate allele counts per cell/SNP
  - `hap1`: Haplotype 1 counts (from phased VCF)
  - `hap2`: Haplotype 2 counts (from phased VCF)

### Common metadata
- **obs**: Cell metadata
  - `n_snps`: Number of SNPs with overlaps
  - `total_counts`: Total counts
  - `ref_counts`: Total reference allele counts (BAM input only)
  - `alt_counts`: Total alternate allele counts (BAM input only)
  - `hap1_counts`: Total haplotype 1 counts (BAM input only)
  - `hap2_counts`: Total haplotype 2 counts (BAM input only)
  - `chemistry`: Library chemistry
  - `sample_id`: Sample identifier
- **var**: SNP metadata
  - `chrom`: Chromosome
  - `pos`: Position
  - `ref`: Reference allele
  - `alt`: Alternate allele
- **uns**: Unstructured metadata
  - `sample_id`: Sample identifier
  - `pipeline`: 'nf-scatac'
  - `data_type`: 'scATAC_allelic_counts' or 'scATAC_allelic_counts_phased'
  - `phased_snps`: Number of phased SNPs (BAM input only)
  - `phasing_rate`: Fraction of SNPs with phasing info (BAM input only)

## Supported Chemistries

| Chemistry | Description |
|-----------|-------------|
| 10x-atac-v1 | 10x Genomics Single Cell ATAC v1 |
| 10x-atac-v2 | 10x Genomics Single Cell ATAC v2 (default) |
| custom | Custom scATAC-seq library prep |

## References

- Issue [#32](https://github.com/Jaureguy760/WASP2-final/issues/32) - scATAC Pipeline
- Issue [#57](https://github.com/Jaureguy760/WASP2-final/issues/57) - nf-core Subworkflow Pattern Compliance
- Issue [#48](https://github.com/Jaureguy760/WASP2-final/issues/48) - Validation & Test Suite
- [ArchR](https://www.archrproject.com/) - scATAC-seq analysis
- [Signac](https://satijalab.org/signac/) - scATAC-seq toolkit
- [AnnData](https://anndata.readthedocs.io/) - Annotated data format
- [nf-test docs](https://code.askimed.com/nf-test/)
- [nf-core guidelines](https://nf-co.re/docs/guidelines)
