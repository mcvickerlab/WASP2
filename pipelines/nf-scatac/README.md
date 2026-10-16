# nf-scatac

Single-cell ATAC-seq counting and coverage workflow

Fragment files contain coordinates, not observed alleles. Their coverage cannot support allelic imbalance testing. This workflow does not run statistical inference: BAM-to-inference aggregation is also not yet implemented. Use the WASP2 single-cell analysis CLI on validated BAM-derived allele counts for downstream testing.

## Features

- **10x Genomics scATAC fragment support** - Direct input from CellRanger ATAC output
- **BAM allele counting** - Reference and alternate counts at heterozygous SNPs
- **Cell barcode filtering** - Optional whitelist for quality-filtered cells
- **Peak region filtering** - Restrict analysis to accessible regions
- **AnnData/H5AD output** - For scverse ecosystem (Scanpy, ArchR, Signac)
- **Zarr output** - For GenVarLoader integration
- **Coverage aggregation** - Sample-level fragment overlap totals, without allele inference
- **Nextflow DSL2** - Reusable modules and subworkflows

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
│       ├── scatac_count_alleles/         # Per-cell fragment overlap coverage
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
- `--min_count` - Reserved bulk-analysis threshold; not used by this workflow

### Output Options
- `--outdir` - Output directory [default: ./results]
- `--create_zarr` - Also output Zarr format for GenVarLoader [default: false]
- `--skip_anndata` - Skip AnnData H5AD creation [default: false]
- `--skip_pseudobulk` - Skip fragment coverage aggregation [default: false]

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

Single-cell counting and coverage aggregation:

```groovy
include { WASP_ALLELIC_SC } from './subworkflows/local/wasp_allelic_sc/main'

WASP_ALLELIC_SC (
    ch_input,       // [ meta, fragments, tbi, barcodes, peaks, bam, bai ]
    ch_vcf          // [ val(meta), path(vcf), path(tbi) ]
)

// Outputs:
// - cell_counts: Per-cell fragment overlap coverage at SNPs
// - anndata: AnnData H5AD files
// - zarr: Zarr directories (if enabled)
// - cell_qc: Cell QC metrics
// - pseudobulk: Aggregated fragment coverage
// - imbalance: Empty (inference is not implemented)
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

Test data locations are configurable via Nextflow params (`--vcf`, `--input`)
or through environment variables consumed by `conf/test_real.config` and
`conf/test_full.config`. The reference runs use GM12878 scATAC-seq
(10x Genomics) with a matched NA12878 phased VCF; adapt to local paths.

## Output

```
results/
├── fragment_coverage/      # Per-cell fragment overlap coverage
│   └── {sample}_fragment_coverage.tsv
├── count_stats/             # Counting statistics
│   └── {sample}_count_stats.tsv
├── anndata/                 # AnnData H5AD files for scverse
│   ├── {sample}_coverage.h5ad          # Fragment input
│   └── {sample}_with_haplotypes.h5ad   # BAM input
├── zarr/                    # Zarr directories (if --create_zarr)
│   ├── {sample}_coverage.zarr/         # Fragment input
│   └── {sample}_with_haplotypes.zarr/  # BAM input
├── cell_qc/                 # Cell QC metrics
│   └── {sample}_cell_qc.tsv
├── pseudobulk/              # Aggregated fragment coverage
│   ├── {sample}_pseudobulk_coverage.tsv
│   └── {sample}_aggregation_stats.tsv
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
  - `data_type`: 'scATAC_fragment_coverage' for fragment input; 'scATAC_allelic_counts_phased' for phased BAM input
  - `phased_snps`: Number of phased SNPs (BAM input only)
  - `phasing_rate`: Fraction of SNPs with phasing info (BAM input only)

## Supported Chemistries

| Chemistry | Description |
|-----------|-------------|
| 10x-atac-v1 | 10x Genomics Single Cell ATAC v1 |
| 10x-atac-v2 | 10x Genomics Single Cell ATAC v2 (default) |
| custom | Custom scATAC-seq library prep |

## References

- Issue [#32](https://github.com/mcvickerlab/WASP2/issues/32) - scATAC Pipeline
- Issue [#57](https://github.com/mcvickerlab/WASP2/issues/57) - nf-core Subworkflow Pattern Compliance
- Issue [#48](https://github.com/mcvickerlab/WASP2/issues/48) - Validation & Test Suite
- [ArchR](https://www.archrproject.com/) - scATAC-seq analysis
- [Signac](https://satijalab.org/signac/) - scATAC-seq toolkit
- [AnnData](https://anndata.readthedocs.io/) - Annotated data format
- [nf-test docs](https://code.askimed.com/nf-test/)
- [nf-core guidelines](https://nf-co.re/docs/guidelines)
