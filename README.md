<p align="center">
  <img src="docs/source/_static/banner.svg" alt="WASP2 - Allele-Specific Analysis Pipeline" width="100%">
</p>

<p align="center">
  <a href="https://github.com/mcvickerlab/WASP2/actions/workflows/ci.yml">
    <img src="https://github.com/mcvickerlab/WASP2/actions/workflows/ci.yml/badge.svg" alt="CI">
  </a>
  <a href="https://mcvickerlab.github.io/WASP2/">
    <img src="https://img.shields.io/badge/docs-GitHub%20Pages-blue" alt="Documentation">
  </a>
  <a href="https://github.com/mcvickerlab/WASP2/blob/main/LICENSE">
    <img src="https://img.shields.io/badge/license-MIT-green" alt="License">
  </a>
</p>

<p align="center">
  <a href="https://mcvickerlab.github.io/WASP2/">Documentation</a> •
  <a href="https://mcvicker.salk.edu/">McVicker Lab</a> •
  <a href="https://github.com/bmvdgeijn/WASP">Original WASP</a>
</p>

---

## Quick Start

```bash
pip install wasp2
```

WASP2 provides three pipelines that run in order:

```bash
# 1. Remap reads to correct mapping bias
wasp2-map make-reads input.bam variants.vcf.gz -o remap_dir/
wasp2-map filter-remapped remap_dir/ -o filtered.bam

# 2. Count allele-specific reads at heterozygous SNPs
wasp2-count count-variants filtered.bam variants.vcf.gz -s sample1

# 3. Test for allelic imbalance
wasp2-analyze find-imbalance counts.tsv -o results.tsv
```

## Authors

- **Aaron Ho** — Creator of WASP2
- **Jeff Jaureguy** — Developer and maintainer
- **[McVicker Lab](https://mcvicker.salk.edu/)**, Salk Institute

## Citation

If you use WASP2 in your research, please cite our paper (coming soon).
