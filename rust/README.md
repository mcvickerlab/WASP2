# WASP2 Rust Acceleration

High-performance BAM allele counting using Rust + PyO3.

## Overview

This Rust extension provides 4-15x speedup for BAM I/O operations in WASP2's counting pipeline by:
- Using `noodles` pure-Rust BAM parser (no htslib dependency)
- Releasing Python GIL for parallel processing
- Zero-copy data transfer via PyO3 bindings

## Performance

**Target:** 7.84s → 0.5-2s for 111K SNP counting (4-15x speedup)

Benchmark from OPTIMIZATION_ANALYSIS.md:
```
Python baseline: 7.84s BAM I/O (98% of 8.03s total)
Rust target: 0.5-2s (conservative estimate)
```

## Installation

### Quick build inside the WASP2 conda env

Prereqs (once per env):
```bash
mamba install -n WASP2 -y clang libclang llvmdev llvm-openmp gcc_linux-64 gxx_linux-64
```

Build/test:
```bash
conda run -n WASP2 bash -lc '
  export LIBCLANG_PATH="$CONDA_PREFIX_2/lib"
  export CC="$CONDA_PREFIX_2/bin/clang"
  export CFLAGS=""    # clear inherited flags that can break hts-sys
  export AR="$CONDA_PREFIX_2/bin/llvm-ar"
  export CPATH="$CONDA_PREFIX_2/lib/clang/21/include:$CONDA_PREFIX_2/include:$CONDA_PREFIX_2/x86_64-conda-linux-gnu/sysroot/usr/include"
  cargo test -q       # or cargo test --tests
'
```

### Prerequisites

```bash
# Install Rust toolchain
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# Install Python build dependencies
pip install maturin>=1.0
```

### Build & Install

```bash
# From project root
cd rust/

# Development build (fast compilation)
maturin develop

# Release build (optimized for performance)
maturin build --release
pip install target/wheels/wasp2_rust-*.whl
```

## Usage

### Python API

```python
from wasp2_rust import BamCounter

# Create counter
counter = BamCounter("/path/to/file.bam")

# Count alleles at SNP positions
# regions: [(chrom, pos, ref, alt), ...]
regions = [("chr1", 12345, "A", "G"), ("chr1", 67890, "C", "T")]
counts = counter.count_alleles(regions, min_qual=20)

# Returns: [(ref_count, alt_count, other_count), ...]
# Example: [(15, 8, 1), (22, 18, 0)]
```

### CLI Integration

The `wasp2-count` command automatically uses Rust acceleration when available:

```bash
# Use Rust (default if extension installed)
wasp2-count sample.bam variants.vcf --regions peaks.bed

# Force Python implementation
wasp2-count sample.bam variants.vcf --no-rust

# Check if Rust is available
python3 -c "from counting.count_alleles import RUST_AVAILABLE; print(f'Rust: {RUST_AVAILABLE}')"
```

## Architecture

### Components

**`rust/src/bam_counter.rs`**
Core BAM counting logic:
- BAM file opening via noodles
- Region-based read iteration
- Base extraction with quality filtering
- Allele classification (ref/alt/other)

**`rust/src/lib.rs`**
PyO3 module definition:
- Exposes `BamCounter` class to Python
- Handles Python ↔ Rust type conversions

**Python integration (`src/counting/count_alleles.py`)**
- `count_snp_alleles_rust()` wrapper
- Automatic fallback to Python if Rust unavailable
- `make_count_df(use_rust=True)` parameter

### Dependencies

```toml
[dependencies]
pyo3 = { version = "0.20", features = ["extension-module"] }
noodles = { version = "0.76", features = ["bam", "sam", "core"] }
rayon = "1.8"
anyhow = "1.0"
```

## Development

### Build for Testing

```bash
# Fast debug build
maturin develop

# Run Python tests with Rust
pytest tests/counting/ -v
```

### Benchmarking

```bash
# Compare Rust vs Python
python3 -c "
from counting.run_counting import run_count_variants

# Rust
run_count_variants('test.bam', 'test.vcf', use_rust=True)

# Python baseline
run_count_variants('test.bam', 'test.vcf', use_rust=False)
"
```

### Debugging

```bash
# Check compilation warnings
maturin build --release

# Test Rust imports
python3 -c "import wasp2_rust; print(dir(wasp2_rust))"

# Enable Rust backtraces
RUST_BACKTRACE=1 python3 script.py
```

## Limitations

1. **No indexed queries yet** - Current implementation iterates all reads
   - Future: Add BAI/CSI index support for chromosome filtering

2. **Chromosome filtering disabled** - All reads processed regardless of chromosome
   - Requires header lookup to map ref_id → chromosome name

3. **No parallelization** - Single-threaded BAM reading
   - Future: Parallel chromosome processing with rayon

## Future Optimizations

From OPTIMIZATION_MASTER_PLAN.md Phase 2:

1. **Indexed BAM queries** (2-3x additional speedup)
   ```rust
   let index = bam::bai::read("file.bam.bai")?;
   let query = reader.query(&header, &index, &region)?;
   ```

2. **Parallel chromosome processing** (2-4x on 8+ cores)
   ```rust
   use rayon::prelude::*;
   chrom_list.par_iter().map(|chrom| count_chromosome(chrom))
   ```

3. **SIMD vectorization** for base quality checks

4. **Memory-mapped BAM** for large files

## Troubleshooting

**Import Error: "No module named 'wasp2_rust'"**
```bash
# Rebuild and install
cd rust/ && maturin develop
```

**Compilation Error: "linker 'cc' not found"**
```bash
# Install build tools
sudo apt-get install build-essential  # Debian/Ubuntu
sudo yum groupinstall "Development Tools"  # RHEL/CentOS
```

**Runtime Error: "BAM file not found"**
```python
# Check path exists
from pathlib import Path
assert Path("file.bam").exists()
```

## References

- **PyO3 Book**: https://pyo3.rs/
- **noodles docs**: https://docs.rs/noodles/
- **OPTIMIZATION_MASTER_PLAN.md**: Full implementation strategy
- **OPTIMIZATION_ANALYSIS.md**: Performance profiling results
