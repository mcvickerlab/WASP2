# WASP2 Rust Acceleration - Summary

## Overview
High-performance Rust implementations of WASP2 computational bottlenecks, providing **2.5-5× speedups** while maintaining 100% output equivalence with Python.

## Performance at a Glance

| Pipeline Stage | Speedup | Threading | Status |
|----------------|---------|-----------|--------|
| Allele Counting | **2.5×** | Rayon (per-chromosome) | ✅ Production |
| Mapping Filter | **4.9×** | BGZF parallel (16 threads) | ✅ Production |
| AI Analysis | **3.4×** | Rayon (per-region) | ✅ Production |

**All outputs verified bit-for-bit identical to Python reference implementation.**

## Quick Start

```bash
# Build Rust extension (one-time setup)
cd rust/
maturin develop --release

# Use automatically (default if extension available)
python -m src.counting count-variants sample.bam variants.vcf

# Enable multi-threading for best performance
export WASP2_RUST_THREADS=16
python -m src.counting count-variants sample.bam variants.vcf
```

## Documentation

- **`rust/README.md`** - Comprehensive Rust implementation guide
- **`docs/RUST_BENCHMARKS.md`** - Detailed performance benchmarks
- **`analysis/publication_benchmarks.pdf`** - Publication-ready figure

## Key Technologies

- **rust-htslib** - Fast BAM/CRAM I/O via htslib C bindings
- **PyO3** - Zero-copy Python ↔ Rust integration  
- **Rayon** - Work-stealing parallelism (CPU-bound tasks)
- **rv** - Statistical distributions (BetaBinomial)
- **FxHash** - High-performance hash maps

## Citation

```bibtex
@article{vandegeijn2015wasp,
  title={WASP: allele-specific software for robust molecular quantitative trait locus discovery},
  author={van de Geijn, Bryce and McVicker, Graham and Gilad, Yoav and Pritchard, Jonathan K},
  journal={Nature Methods},
  volume={12},
  number={11},
  pages={1061--1063},
  year={2015}
}
```

## License
Same as WASP2 main project.
