# nf-scatac: Usage

## Introduction

**nf-scatac** is a Nextflow DSL2 workflow for single-cell ATAC-seq counting and coverage. BAM inputs support allele counting; fragment inputs support overlap coverage only. Fragment coordinates cannot identify observed alleles. Statistical inference is not run for either path; BAM-to-inference aggregation remains unimplemented. See the [output documentation](output.md) and use the WASP2 single-cell CLI separately for inference on validated BAM-derived counts.

## Pipeline Summary

1. Input validation (samplesheet, VCF)
2. VCF to BED conversion for SNP positions
3. BAM allele counting or fragment overlap coverage at heterozygous SNPs
4. AnnData conversion and optional fragment coverage aggregation

## Quick Start

```bash
nextflow run nf-scatac \
    --input samplesheet.csv \
    --vcf variants.vcf.gz \
    --outdir results \
    -profile docker
```

## Samplesheet Input

The pipeline requires a samplesheet CSV file with the following columns:

| Column | Required | Description |
|--------|----------|-------------|
| `sample` | Yes | Unique sample identifier |
| `fragments` | Yes* | Path to 10x fragments.tsv.gz file |
| `cellranger_dir` | Yes* | Path to CellRanger ATAC output directory |
| `barcode_tag` | No | BAM tag for cell barcodes (default: CB) |
| `chemistry` | No | Library chemistry (default: 10x-atac-v2) |

*Either `fragments` or `cellranger_dir` is required.

### Example samplesheet (fragments):

```csv
sample,fragments,cellranger_dir,barcode_tag,chemistry
GM12878_rep1,/data/GM12878/fragments.tsv.gz,,CB,10x-atac-v2
GM12878_rep2,/data/GM12878_rep2/fragments.tsv.gz,,CB,10x-atac-v2
```

### Example samplesheet (CellRanger ATAC output):

```csv
sample,fragments,cellranger_dir,barcode_tag,chemistry
GM12878_cellranger,,/data/cellranger_atac/GM12878,CB,10x-atac-v2
```

## Required Parameters

| Parameter | Description |
|-----------|-------------|
| `--input` | Path to samplesheet CSV |
| `--vcf` | Phased VCF/BCF file with variants (must be bgzipped and indexed) |

## Optional Parameters

### Coverage options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min_fragments_per_cell` | 1000 | Minimum overlap total per cell for the fragment path |
| `--min_cells_per_snp` | 3 | Minimum covered cells per SNP in aggregated fragment coverage |
| `--skip_pseudobulk` | false | Skip fragment coverage aggregation |
| `--skip_anndata` | false | Skip AnnData output |
| `--create_zarr` | false | Also write Zarr output |

The legacy `--min_count` and `--pseudocount` parameters do not run inference in this workflow.

### Single-Cell Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--barcode_tag` | 'CB' | Default BAM tag for cell barcodes |
| `--chemistry` | '10x-atac-v2' | Default library chemistry |

### Output Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outdir` | './results' | Output directory |
| `--publish_dir_mode` | 'copy' | Publishing mode: 'copy', 'symlink', 'link' |

## Running with Profiles

### Docker

```bash
nextflow run nf-scatac -profile docker --input samplesheet.csv ...
```

### Singularity

```bash
nextflow run nf-scatac -profile singularity --input samplesheet.csv ...
```

### Stub Test (CI/CD)

Run fast stub tests that validate workflow structure:

```bash
nextflow run nf-scatac -profile test_stub -stub-run
```

### Integration Test

Run with real GM12878 scATAC-seq data:

```bash
nextflow run nf-scatac -profile test_real,singularity
```

## Example Commands

### Cell-level counting and coverage

```bash
nextflow run pipelines/nf-scatac \
    --input samplesheet.csv \
    --vcf phased_variants.vcf.gz \
    --outdir results \
    -profile docker
```

Use `--skip_pseudobulk` to omit fragment coverage aggregation. Phase information is used for BAM haplotype layers; no phased or unphased inference is performed by this workflow.

## Supported Chemistries

| Chemistry | Description |
|-----------|-------------|
| `10x-atac-v1` | 10x Genomics Single Cell ATAC v1 |
| `10x-atac-v2` | 10x Genomics Single Cell ATAC v2 (default) |
| `custom` | Custom scATAC-seq library prep |

## Resources

Configure `--max_cpus`, `--max_memory`, and `--max_time` for the input size and execution environment. Check the execution report for observed resource use.

## Troubleshooting

### Common Issues

1. **Missing VCF index**: Ensure your VCF is bgzipped and indexed with tabix
   ```bash
   bgzip variants.vcf
   tabix -p vcf variants.vcf.gz
   ```

2. **Fragments file not indexed**: Index with tabix
   ```bash
   tabix -p bed fragments.tsv.gz
   ```

3. **Memory errors**: Increase `--max_memory` or use a profile with more resources

4. **No cells found**: Check that barcode_tag matches your BAM (default: CB)

### Resume Failed Runs

```bash
nextflow run nf-scatac ... -resume
```

## Citation

If you use nf-scatac, please cite:

- WASP2: [GitHub](https://github.com/your-org/WASP2)
- Nextflow: [Nextflow](https://www.nextflow.io/)
