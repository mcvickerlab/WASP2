<p align="center">
  <img src="https://raw.githubusercontent.com/mcvickerlab/WASP2/master/docs/source/_static/banner.svg" alt="WASP2 - Allele-Specific Analysis Pipeline" width="100%">
</p>

<p align="center">
  <a href="https://pypi.org/project/wasp2/">
    <img src="https://img.shields.io/pypi/v/wasp2" alt="PyPI">
  </a>
  <a href="https://bioconda.github.io/recipes/wasp2/README.html">
    <img src="https://img.shields.io/conda/vn/bioconda/wasp2" alt="Bioconda">
  </a>
  <a href="https://github.com/mcvickerlab/WASP2/actions/workflows/ci.yml">
    <img src="https://github.com/mcvickerlab/WASP2/actions/workflows/ci.yml/badge.svg" alt="CI">
  </a>
  <a href="https://mcvickerlab.github.io/WASP2/">
    <img src="https://img.shields.io/badge/docs-GitHub%20Pages-blue" alt="Documentation">
  </a>
</p>

## Installation

### Bioconda

```bash
mamba install -c conda-forge -c bioconda wasp2
```

### PyPI

```bash
pip install wasp2
```

The PyPI package does not install external tools such as `samtools`,
`bcftools`, or `bedtools`; install those separately.

### Docker

```bash
docker pull ghcr.io/mcvickerlab/wasp2:1.4.1
docker run --rm ghcr.io/mcvickerlab/wasp2:1.4.1 wasp2-count --help
```

### Singularity/Apptainer

```bash
singularity pull wasp2.sif docker://ghcr.io/mcvickerlab/wasp2:1.4.1
singularity exec wasp2.sif wasp2-count --help
```

## CLI Tools

WASP2 installs four command-line entry points:

- `wasp2-map`
- `wasp2-count`
- `wasp2-analyze`
- `***REMOVED***`

## Quick Start

### 1. Correct mapping bias

```bash
wasp2-map make-reads input.bam variants.vcf.gz -s sample1 -o remap_dir

# Realign remap_dir/*_swapped_alleles_r1.fq and r2.fq with the same aligner
# and settings used for the original BAM, then:

wasp2-map filter-remapped remapped.bam \
  -j remap_dir/input_wasp_data_files.json \
  -o filtered.bam
```

### 2. Count alleles

```bash
wasp2-count count-variants filtered.bam variants.vcf.gz -s sample1 -o counts.tsv
```

### 3. Test for imbalance

```bash
wasp2-analyze find-imbalance counts.tsv -o ai_results.tsv
```

## Single-Cell Example

```bash
wasp2-count count-variants-sc \
  cellranger.bam \
  variants.vcf.gz \
  barcodes.tsv \
  --samples sample1 \
  --feature genes.gtf \
  --out_file allele_counts.h5ad

wasp2-analyze find-imbalance-sc \
  allele_counts.h5ad \
  barcode_groups.tsv \
  --sample sample1 \
  -o ai_results.tsv
```

## CVPC Utilities

```bash
***REMOVED*** inventory --output inventory.tsv
***REMOVED*** manifest --output manifest.csv
***REMOVED*** validate
```

See the [documentation](https://mcvickerlab.github.io/WASP2/) for complete
usage, tutorials, and API details.
