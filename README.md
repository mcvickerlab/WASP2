<h1 align="center">
<img src="doc/wasp2_hex_logo_v1.png" width="300">
</h1>

<p align="center">
  <a href="https://github.com/Jaureguy760/WASP2-final/actions/workflows/ci.yml">
    <img src="https://github.com/Jaureguy760/WASP2-final/actions/workflows/ci.yml/badge.svg" alt="CI">
  </a>
  <a href="https://jaureguy760.github.io/WASP2-final/">
    <img src="https://img.shields.io/badge/docs-GitHub%20Pages-blue" alt="Documentation">
  </a>
  <a href="https://pypi.org/project/wasp2/">
    <img src="https://img.shields.io/pypi/v/wasp2?logo=pypi&logoColor=white" alt="PyPI">
  </a>
  <a href="https://github.com/Jaureguy760/WASP2-final/blob/main/LICENSE">
    <img src="https://img.shields.io/badge/license-MIT-green" alt="License">
  </a>
  <img src="https://img.shields.io/badge/python-3.10+-blue?logo=python&logoColor=white" alt="Python">
  <img src="https://img.shields.io/badge/rust-1.70+-orange?logo=rust&logoColor=white" alt="Rust">
</p>

# WASP2: Allele-specific pipeline for unbiased read mapping and allelic imbalance analysis

WASP2 provides high-performance tools for allele-specific analysis of genomic data, correcting mapping biases and detecting allelic imbalance in bulk and single-cell experiments.

## Features

- **Allele-specific counting** - Count reads at heterozygous variants in RNA-seq and ATAC-seq data
- **Mapping bias correction** - WASP filtering to remove reads with allelic mapping bias
- **Statistical analysis** - Beta-binomial model for detecting allelic imbalance
- **Single-cell support** - Per-cell allelic counts and cell-type-specific imbalance analysis
- **High-performance variant I/O** - Support for VCF, BCF (~7x faster), and PGEN (~25x faster) formats
- **Rust-accelerated core** - Critical paths implemented in Rust via PyO3
- **Container support** - Docker, Singularity, and Nextflow-ready

## Installation

```bash
pip install wasp2
```

Pre-built wheels are available for Linux (x86_64, aarch64) and macOS (Intel, Apple Silicon).

### Optional dependencies

```bash
# Faster VCF parsing (~7x speedup)
pip install wasp2[cyvcf2]

# PLINK2 PGEN format support (~25x speedup)
pip install wasp2[plink]
```

### Container

```bash
docker pull ghcr.io/jaureguy760/wasp2-final:latest
# or for HPC:
singularity pull wasp2.sif docker://ghcr.io/jaureguy760/wasp2-final:latest
```

## Quick Start

```bash
# Count allele-specific reads
wasp2-count count-variants reads.bam variants.vcf.gz -s sample1 -r regions.bed

# WASP mapping bias correction
wasp2-map make-reads reads.bam variants.vcf.gz
# ... remap with your aligner ...
wasp2-map filter-remapped remapped.bam --json wasp_data_files.json

# Detect allelic imbalance
wasp2-analyze find-imbalance counts.tsv -o results.tsv
```

## Documentation

Full documentation with tutorials and API reference: **[jaureguy760.github.io/WASP2-final](https://jaureguy760.github.io/WASP2-final/)**

## Contributing

```bash
# Clone and set up development environment
git clone https://github.com/Jaureguy760/WASP2-final.git
cd WASP2-final
conda env create -f environment.yml
conda activate WASP2

# Build Rust extension
maturin develop --release -m rust/Cargo.toml

# Run tests
pytest tests/ -v
```

## Citation

If you use WASP2 in your research, please cite our paper (coming soon).
