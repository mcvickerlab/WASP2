# WASP2 Rust Acceleration

High-performance implementations of WASP2 bottlenecks using Rust + PyO3.

## Performance Summary

| Module | Implementation | Speedup | Notes |
|--------|---------------|---------|-------|
| **Counting** | ✅ Complete | **2.5×** single-thread<br>**2.5×** multi-thread* | *Chr10 test has 1 chromosome |
| **Mapping Filter** | ✅ Complete | **5×** (16 threads) | BGZF parallel decompression |
| **Analysis** | ✅ Complete | **3.4×** (parallel regions) | Beta-binomial model |

### Benchmarks (chr10 test dataset)

**Counting** (111K SNPs, with precomputed intermediates):
- Python: 6.86s → Rust: 2.70s (**2.5× faster**)
- Core counting loop: Python ~3.9s → Rust ~0.2s (**18× faster**)

**Mapping Filter** (20K read pairs):
- Python (1 thread): 0.44s → Rust (16 threads): 0.09s (**4.9× faster**)

**Analysis** (43 regions, beta-binomial AI detection):
- Python: 0.28s → Rust: 0.08s (**3.4× faster**)

All implementations produce **identical outputs** to Python (zero mismatches verified).

## Installation

### Prerequisites

```bash
# Install Rust toolchain
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# Install maturin for building Python extensions
pip install maturin>=1.0
```

### Build & Install

```bash
# From rust/ directory
cd rust/

# Development build (fast compilation)
maturin develop

# Release build (optimized for performance)
maturin develop --release
```

The extension will be automatically available to Python as `wasp2_rust`.

## Usage

### Automatic Integration

WASP2 automatically uses Rust acceleration when available. No code changes needed:

```bash
# Rust is used automatically if extension is built
python -m src.counting count-variants sample.bam variants.vcf

# Force Python implementation
python -m src.counting count-variants sample.bam variants.vcf --no-rust
```

### Threading

Enable multi-threading via environment variable:

```bash
# Use 16 threads for parallel processing
export WASP2_RUST_THREADS=16
python -m src.counting count-variants sample.bam variants.vcf
```

### Python API

#### Counting

```python
from wasp2_rust import BamCounter

# Create counter
counter = BamCounter("/path/to/file.bam")

# Count alleles at SNP positions
regions = [("chr1", 12345, "A", "G"), ("chr1", 67890, "C", "T")]
counts = counter.count_alleles(regions, min_qual=0, threads=1)

# Returns: [(ref_count, alt_count, other_count), ...]
```

#### Mapping Filter

```python
from wasp2_rust import filter_bam_wasp

# Filter remapped BAM to remove reads that moved
kept, removed_moved, removed_missing = filter_bam_wasp(
    to_remap_bam="original.bam",
    remapped_bam="remapped.bam",
    remap_keep_bam="filtered.bam",
    threads=16
)
```

#### Analysis

```python
from wasp2_rust import analyze_imbalance

# Run beta-binomial allelic imbalance analysis
results = analyze_imbalance(
    tsv_path="counts.tsv",
    min_count=10,
    pseudocount=1,
    method="single"
)

# Results is a list of dicts with: region, pval, fdr_pval, mu, etc.
```

## Architecture

### Key Technologies

- **rust-htslib**: Fast BAM/CRAM I/O (bindings to htslib C library)
- **PyO3**: Zero-copy Python ↔ Rust bindings
- **Rayon**: Work-stealing parallelism for multi-threading
- **rv**: Statistical distributions (BetaBinomial for analysis)
- **FxHash**: Fast hash maps for read deduplication

### Modules

**`bam_counter.rs`** - Allele counting from BAM files
- Per-chromosome parallel processing with Rayon
- Read deduplication via FxHashSet
- Encounter-order SNP assignment (matches Python behavior exactly)

**`mapping_filter.rs`** - WASP remap filter
- Parallel BGZF block decompression (htslib threading)
- Keeps reads that return to original position with all haplotype copies

**`analysis.rs`** - Beta-binomial allelic imbalance detection
- Parallel per-region likelihood optimization
- Uses rv crate for BetaBinomial distribution
- Golden section search for parameter optimization
- Benjamini-Hochberg FDR correction

**`bam_remapper.rs`** - Allele swapping for WASP remapping
- Parse bedtools intersect output
- Generate haplotype FASTQs for realignment
- Support for phased/unphased variants

## Development

### Running Benchmarks

```bash
# Counting benchmark (requires precomputed intermediates)
python analysis/run_counting_bench.py \
  --vcf-bed /tmp/counting_precomputed/filter_chr10.bed \
  --intersect-bed /tmp/counting_precomputed/filter_chr10_intersect.bed

# Mapping filter benchmark
python analysis/run_mapping_filter_bench.py

# Analysis benchmark
python analysis/run_analysis_bench.py
```

### Testing

```bash
# Run Rust unit tests
cd rust/
cargo test

# Run Python integration tests with Rust
cd ..
pytest tests/ -v
```

## Dependencies

```toml
[dependencies]
pyo3 = { version = "0.20", features = ["extension-module"] }
rust-htslib = "0.47"  # BAM/CRAM I/O
rayon = "1.8"          # Parallelism
rv = "0.19"            # Beta-binomial distribution
statrs = "0.17"        # Chi-squared distribution
rustc-hash = "1.1"     # Fast FxHash
anyhow = "1.0"         # Error handling
```

## Troubleshooting

**Import Error: "No module named 'wasp2_rust'"**
```bash
cd rust/ && maturin develop --release
```

**Build Error: "stddef.h not found"**
```bash
# Set C include path
export C_INCLUDE_PATH=/usr/lib/gcc/x86_64-redhat-linux/11/include
export LIBCLANG_PATH=/path/to/conda/envs/WASP2/lib
maturin develop --release
```

**Performance Not Improving**
- Ensure release build: `maturin develop --release`
- Enable threading: `export WASP2_RUST_THREADS=16`
- Check CPU count: `nproc` (use 8-16 threads typically)

## Citation

If you use this Rust acceleration in your research, please cite:

```
van de Geijn B, McVicker G, Gilad Y, Pritchard JK (2015)
WASP: allele-specific software for robust molecular quantitative trait locus discovery.
Nature Methods 12(11):1061-1063
```

## License

Same as WASP2 main project (see LICENSE in root directory).
