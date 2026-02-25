<p align="center">
  <img src="docs/source/_static/banner.svg" alt="WASP2 - Allele-Specific Analysis Pipeline" width="100%">
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
  <a href="https://github.com/mcvickerlab/WASP2/blob/master/LICENSE">
    <img src="https://img.shields.io/badge/license-MIT-green" alt="License">
  </a>
</p>

<p align="center">
  <a href="https://mcvickerlab.github.io/WASP2/">Documentation</a> •
  <a href="https://mcvicker.salk.edu/">McVicker Lab</a> •
  <a href="https://github.com/bmvdgeijn/WASP">Original WASP</a>
</p>

---

## Installation

### Recommended: Bioconda

```bash
mamba install -c conda-forge -c bioconda wasp2
```

Installs WASP2 and all dependencies (samtools, bcftools, bedtools, htslib) automatically. Available for Linux (x86_64, aarch64) and macOS (Intel, Apple Silicon). Requires [miniforge](https://github.com/conda-forge/miniforge).

### Via PyPI

```bash
pip install wasp2
```

Pre-built wheels for Linux (x86_64, aarch64) and macOS (Intel, Apple Silicon) with Python 3.10-3.13. The Rust extension and htslib are bundled in the wheel. Requires samtools, bcftools, and bedtools installed separately.

### For development

```bash
git clone https://github.com/mcvickerlab/WASP2.git
cd WASP2
pixi install        # resolves all dependencies including Rust toolchain
pixi run verify     # build + test
```

See the [documentation](https://mcvickerlab.github.io/WASP2/) for Docker and Singularity install options.

## Quick Start

WASP2 has three steps that run in order:

**Step 1: Remap reads** to correct mapping bias

```bash
wasp2-map make-reads input.bam variants.vcf.gz -s sample1 -o remap_dir/
# Realign the swapped-allele reads with your aligner, then:
wasp2-map filter-remapped remapped.bam -j remap_dir/sample1_wasp_data_files.json -o filtered.bam
```

**Step 2: Count alleles** at heterozygous SNPs

```bash
wasp2-count count-variants filtered.bam variants.vcf.gz -s sample1
```

**Step 3: Test for allelic imbalance**

```bash
wasp2-analyze find-imbalance counts.tsv -o results.tsv
```

See the [documentation](https://mcvickerlab.github.io/WASP2/) for detailed usage, single-cell workflows, and supported variant formats (VCF, BCF, PGEN).

## Authors

- **Aaron Ho** — Creator of WASP2
- **Jeff Jaureguy** — Developer and maintainer
- **[McVicker Lab](https://mcvicker.salk.edu/)**, Salk Institute

## Citation

If you use WASP2 in your research, please cite our paper (coming soon).
