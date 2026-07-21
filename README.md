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
docker pull ghcr.io/mcvickerlab/wasp2:1.5.0
docker run --rm ghcr.io/mcvickerlab/wasp2:1.5.0 wasp2-count --help
```

### Singularity/Apptainer

```bash
singularity pull wasp2.sif docker://ghcr.io/mcvickerlab/wasp2:1.5.0
singularity exec wasp2.sif wasp2-count --help
```

## CLI Tools

WASP2 installs three command-line entry points:

- `wasp2-map`
- `wasp2-count`
- `wasp2-analyze`

## Quick Start

> **Prerequisite — normalize your VCF first.** WASP2's allelic-imbalance model is biallelic.
> Split multi-allelic records and left-align before running:
>
> ```bash
> bcftools norm -m- -f reference.fa input.vcf.gz -Oz -o input.norm.vcf.gz
> bcftools index -t input.norm.vcf.gz
> ```
>
> By default WASP2 keeps **biallelic SNVs only** (multi-allelic sites are dropped, equivalent to
> `bcftools view -m2 -M2`). Pass `--include-multiallelic` to `count-variants`, `count-variants-sc`,
> or `make-reads` to instead emit one row per ALT at any remaining multi-allelic sites.

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

For a bulk cohort, use one final indexed BAM per donor and map donor identifiers
to VCF sample names explicitly:

```text
donor_id	vcf_sample	bam
donor1	VCF_sample_1	/path/to/donor1.bam
donor2	VCF_sample_2	/path/to/donor2.bam
```

```bash
wasp2-count count-cohort donors.tsv variants.vcf.gz cohort_counts \
  --unit snv --regions atac_peaks.bed
```

This creates a new locked directory containing `counts.tsv.gz`, normalized
`donors.tsv`, and `count_manifest.json`. In SNV mode, regions restrict included
sites but do not combine their tests. Use `--unit feature --regions peaks.bed`
to retain peak identifiers for feature-level analysis.

### 3. Test for imbalance

```bash
wasp2-analyze find-imbalance counts.tsv -o ai_results.tsv
```

### Donor-local cohort analysis

Analyze the locked count directory (or its `count_manifest.json` or
`counts.tsv.gz`) with an explicit unit. SNV analysis rejects `--phased` and tests
each variant independently within each donor:

```bash
wasp2-analyze find-imbalance cohort_counts \
  --scope per-donor --unit snv --min-donor-observations 50 \
  -o atac_snv_results.tsv
```

`--dispersion-scope per-donor` is the default. `--dispersion-scope global`
fits one nuisance model over unique eligible donor-SNV rows and reuses that exact
fit for every included donor; testing and BH correction remain donor-local. Both
`single` and `linear` dispersion models are supported.

Feature analysis preserves one SNV in every overlapping feature. `peak` remains
an alias for `feature`:

```bash
wasp2-count count-cohort donors.tsv variants.vcf.gz feature_counts \
  --unit feature --regions peaks.bed

wasp2-analyze find-imbalance feature_counts \
  --scope per-donor --unit feature --model linear --phased \
  -o atac_feature_results.tsv
```

The donor-local route defaults to pseudocount 0. It writes
`atac_feature_results.tsv`, `atac_feature_results.dispersion.tsv`,
`atac_feature_results.qc.tsv`, and `atac_feature_results.provenance.json` and
refuses to overwrite any of them. For standalone legacy count tables only,
`sample` is accepted and normalized to `donor_id`; providing both columns is an
error. Use `--expected-manifest-sha256` to require an externally trusted locked
bundle manifest.

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

See the [documentation](https://mcvickerlab.github.io/WASP2/) for complete
usage, tutorials, and API details.
