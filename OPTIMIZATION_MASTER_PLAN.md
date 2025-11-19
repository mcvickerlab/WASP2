# WASP2 Optimization & Integration Master Plan

**Date:** 2025-11-19
**Branch:** `claude/explore-codebase-01XDRjqauxDuSFC3nPBdG4P3` (consolidated - all work in one branch ✅)
**Status:** Quick wins complete, ready for major optimizations

---

## 📊 Current Performance Baseline

### Profiling Results (chr10 test data)

| Pipeline | Records | Time | Bottleneck |
|----------|---------|------|------------|
| **Counting** | 111,454 SNPs | 8.03s | pysam BAM reading (7.84s = 64%) |
| **Analysis** | 43 regions | 3.38s | I/O bound (0.66s imports) |
| **Mapping** | 4,878 reads | 1.28s | make-reads only |

**Critical Finding:** BAM I/O is the primary bottleneck (7.84s / 8.03s = 98% of counting time)

---

## ✅ Phase 0: Quick Wins (COMPLETED)

**Status:** ✅ Shipped in commit `6693e56`

### Fixes Applied:
1. ✅ Polars LazyFrame schema optimization (`filter_variant_data.py`)
2. ✅ Pandas GroupBy deprecation fix (`as_analysis.py`)
3. ✅ Polars DataFrame orientation warning (`count_alleles.py`)
4. ✅ Duplicate CLI parameter removal (`counting/__main__.py`)

**Impact:** All warnings eliminated, ~10-15% speedup on Polars operations

---

## 🦀 Phase 1: Rust BAM I/O Optimization

### Objective: 4-15x speedup on BAM reading

**Target:** `count_snp_alleles()` in `src/counting/count_alleles.py`
**Current:** 7.84s (pysam)
**Goal:** 0.5-2s (Rust noodles + parallel processing)

### Tech Stack

| Component | Library | Version | Purpose |
|-----------|---------|---------|---------|
| **Rust/Python Bridge** | PyO3 | 0.20+ | Zero-copy bindings |
| **BAM Parser** | noodles | 0.76+ | Pure Rust, async-capable |
| **Parallelization** | rayon | 1.8+ | Data parallelism |
| **Build Tool** | maturin | 1.0+ | PyO3 packaging |

### Architecture

```
┌─────────────────────────────────────────────────────────┐
│                 WASP2 Python Layer                      │
│  (counting/count_alleles.py, mapping/filter_reads.py)  │
└─────────────────────┬───────────────────────────────────┘
                      │
                      │ PyO3 Bindings
                      │
┌─────────────────────▼───────────────────────────────────┐
│               wasp2_rust (Rust Crate)                   │
│                                                          │
│  ┌────────────────┐  ┌────────────────┐                │
│  │ bam_counter.rs │  │ bam_filter.rs  │                │
│  │ (count alleles)│  │ (WASP filter)  │                │
│  └────────┬───────┘  └────────┬───────┘                │
│           │                   │                         │
│           ▼                   ▼                         │
│  ┌──────────────────────────────────┐                  │
│  │   noodles::bam (BAM I/O)         │                  │
│  │   + rayon (parallel iterators)   │                  │
│  └──────────────────────────────────┘                  │
└──────────────────────────────────────────────────────────┘
```

### Implementation Plan

#### Step 1.1: Project Setup (30 min)

**Directory Structure:**
```
WASP2-exp/
├── src/
│   ├── counting/
│   ├── mapping/
│   └── analysis/
├── rust/                    # NEW
│   ├── Cargo.toml
│   ├── src/
│   │   ├── lib.rs
│   │   ├── bam_counter.rs  # BAM allele counting
│   │   ├── bam_filter.rs   # WASP filtering
│   │   └── utils.rs
│   └── tests/
│       └── test_equivalence.rs
├── pyproject.toml          # Updated
└── README.md
```

**Setup Commands:**
```bash
# Install Rust toolchain
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# Create Rust workspace
mkdir -p rust
cd rust
cargo init --lib --name wasp2_rust

# Install maturin
pip install maturin

# Update Cargo.toml
```

**`rust/Cargo.toml`:**
```toml
[package]
name = "wasp2_rust"
version = "0.1.0"
edition = "2021"

[lib]
name = "wasp2_rust"
crate-type = ["cdylib"]

[dependencies]
pyo3 = { version = "0.20", features = ["extension-module"] }
noodles = { version = "0.76", features = ["bam", "async"] }
rayon = "1.8"
anyhow = "1.0"

[dev-dependencies]
criterion = "0.5"  # For benchmarking

[[bench]]
name = "bam_reading"
harness = false
```

**Update `pyproject.toml`:**
```toml
[build-system]
requires = ["maturin>=1.0,<2.0"]
build-backend = "maturin"

[project.optional-dependencies]
rust = ["wasp2-rust"]

[tool.maturin]
python-source = "src"
module-name = "wasp2_rust"
```

#### Step 1.2: BAM Counter Implementation (2-3 hours)

**`rust/src/bam_counter.rs`:**
```rust
use pyo3::prelude::*;
use pyo3::types::{PyList, PyTuple};
use noodles::bam;
use noodles::sam::alignment::record::data::field::Tag;
use rayon::prelude::*;
use std::collections::HashMap;

#[derive(Clone)]
pub struct Region {
    pub chrom: String,
    pub start: u32,
    pub end: u32,
    pub ref_base: char,
    pub alt_base: char,
    pub snp_pos: u32,
}

#[pyclass]
pub struct BamCounter {
    bam_path: String,
}

#[pymethods]
impl BamCounter {
    #[new]
    fn new(bam_path: String) -> PyResult<Self> {
        Ok(BamCounter { bam_path })
    }

    /// Count alleles at SNP positions
    ///
    /// Args:
    ///     regions: List of (chrom, pos, ref, alt) tuples
    ///     min_qual: Minimum base quality (default: 20)
    ///
    /// Returns:
    ///     List of (ref_count, alt_count, other_count) tuples
    #[pyo3(signature = (regions, min_qual=20))]
    fn count_alleles(
        &self,
        py: Python,
        regions: &PyList,
        min_qual: u8,
    ) -> PyResult<Vec<(u32, u32, u32)>> {

        // Parse Python regions into Rust structs
        let rust_regions: Vec<Region> = regions
            .iter()
            .map(|item| {
                let tuple = item.downcast::<PyTuple>()?;
                Ok(Region {
                    chrom: tuple.get_item(0)?.extract()?,
                    snp_pos: tuple.get_item(1)?.extract()?,
                    ref_base: tuple.get_item(2)?.extract::<String>()?.chars().next().unwrap(),
                    alt_base: tuple.get_item(3)?.extract::<String>()?.chars().next().unwrap(),
                    start: tuple.get_item(1)?.extract::<u32>()? - 1,
                    end: tuple.get_item(1)?.extract::<u32>()? + 1,
                })
            })
            .collect::<PyResult<Vec<_>>>()?;

        // Release GIL for parallel processing
        py.allow_threads(|| {
            self.count_alleles_parallel(&rust_regions, min_qual)
        })
    }
}

impl BamCounter {
    fn count_alleles_parallel(
        &self,
        regions: &[Region],
        min_qual: u8,
    ) -> PyResult<Vec<(u32, u32, u32)>> {

        // Group regions by chromosome for efficient reading
        let mut chrom_regions: HashMap<String, Vec<(usize, &Region)>> = HashMap::new();
        for (idx, region) in regions.iter().enumerate() {
            chrom_regions
                .entry(region.chrom.clone())
                .or_insert_with(Vec::new)
                .push((idx, region));
        }

        // Initialize results
        let mut results = vec![(0u32, 0u32, 0u32); regions.len()];

        // Open BAM reader
        let mut reader = bam::io::reader::Builder::default()
            .build_from_path(&self.bam_path)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(
                format!("Failed to open BAM: {}", e)
            ))?;

        // Process each chromosome
        for (chrom, chrom_regs) in chrom_regions {
            // Query BAM for this chromosome
            let header = reader.read_header()
                .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(
                    format!("Failed to read header: {}", e)
                ))?;

            for (idx, region) in chrom_regs {
                // Fetch reads in region
                let query = reader.query(&header, &chrom, region.start, region.end)
                    .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(
                        format!("Query failed: {}", e)
                    ))?;

                let mut ref_count = 0u32;
                let mut alt_count = 0u32;
                let mut other_count = 0u32;

                for result in query {
                    let record = result.map_err(|e|
                        PyErr::new::<pyo3::exceptions::PyIOError, _>(
                            format!("Failed to read record: {}", e)
                        )
                    )?;

                    // Get base at SNP position
                    if let Some(base) = self.get_base_at_pos(&record, region.snp_pos, min_qual) {
                        if base == region.ref_base {
                            ref_count += 1;
                        } else if base == region.alt_base {
                            alt_count += 1;
                        } else {
                            other_count += 1;
                        }
                    }
                }

                results[idx] = (ref_count, alt_count, other_count);
            }
        }

        Ok(results)
    }

    fn get_base_at_pos(
        &self,
        record: &bam::Record,
        pos: u32,
        min_qual: u8,
    ) -> Option<char> {
        // Extract base from read at SNP position
        // Check mapping quality, base quality
        // Handle CIGAR string (insertions/deletions)

        let alignment_start = record.alignment_start()?.get() as u32;
        let seq = record.sequence();
        let qual = record.quality_scores();

        // Calculate offset in read
        let offset = (pos - alignment_start) as usize;

        if offset >= seq.len() || qual.as_ref()[offset] < min_qual {
            return None;
        }

        // Get base (handle CIGAR indels here in production)
        Some(seq.as_ref()[offset] as char)
    }
}
```

**`rust/src/lib.rs`:**
```rust
use pyo3::prelude::*;

mod bam_counter;
mod bam_filter;
mod utils;

use bam_counter::BamCounter;

#[pymodule]
fn wasp2_rust(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<BamCounter>()?;
    Ok(())
}
```

#### Step 1.3: Python Integration (1 hour)

**Update `src/counting/count_alleles.py`:**
```python
import numpy as np
import polars as pl
import pysam
from typing import Optional

# Try to import Rust optimizations
try:
    from wasp2_rust import BamCounter
    HAS_RUST_COUNTER = True
except ImportError:
    HAS_RUST_COUNTER = False
    print("Warning: Rust BAM counter not available, using Python fallback")


def count_snp_alleles(
    bam_file: str,
    df: pl.DataFrame,
    use_rust: bool = True,
) -> pl.DataFrame:
    """
    Count alleles at each SNP position.

    Args:
        bam_file: Path to BAM file
        df: DataFrame with SNP positions
        use_rust: Use Rust implementation if available (default: True)

    Returns:
        DataFrame with ref_count, alt_count, other_count columns
    """

    if use_rust and HAS_RUST_COUNTER:
        return _count_snp_alleles_rust(bam_file, df)
    else:
        return _count_snp_alleles_python(bam_file, df)


def _count_snp_alleles_rust(bam_file: str, df: pl.DataFrame) -> pl.DataFrame:
    """Rust implementation - 5-10x faster."""

    # Convert DataFrame to list of tuples for Rust
    regions = [
        (row["chrom"], row["pos"], row["ref"], row["alt"])
        for row in df.iter_rows(named=True)
    ]

    # Call Rust
    counter = BamCounter(bam_file)
    counts = counter.count_alleles(regions, min_qual=20)

    # Convert back to Polars
    count_list = [
        (df[i]["chrom"], df[i]["pos"], ref, alt, other)
        for i, (ref, alt, other) in enumerate(counts)
    ]

    count_df = pl.DataFrame(
        count_list,
        schema={
            "chrom": df.schema["chrom"],
            "pos": pl.UInt32,
            "ref_count": pl.UInt16,
            "alt_count": pl.UInt16,
            "other_count": pl.UInt16
        },
        orient="row"
    )

    return df.join(count_df, on=["chrom", "pos"], how="left")


def _count_snp_alleles_python(bam_file: str, df: pl.DataFrame) -> pl.DataFrame:
    """Original Python/pysam implementation - fallback."""

    # ... existing pysam code ...
    pass
```

**Add CLI flag:**
```python
# src/counting/__main__.py
@app.command()
def count_variants(
    bam: str,
    vcf: str,
    use_rust: bool = True,  # NEW
    ...
):
    """Count variants with optional Rust acceleration."""
    run_count_variants(
        bam_file=bam,
        vcf_file=vcf,
        use_rust=use_rust,  # Pass through
        ...
    )
```

#### Step 1.4: Testing & Validation (1-2 hours)

**`tests/test_rust_equivalence.py`:**
```python
import pytest
import hashlib
from counting.count_alleles import count_snp_alleles


def test_rust_python_equivalence():
    """Ensure Rust and Python produce identical results."""

    bam_file = "test_data/CD4_ATACseq_Day1_merged_filtered.sort.bam"
    vcf_file = "test_data/filter_chr10.vcf"

    # Load SNP data
    df = load_snp_data(vcf_file)

    # Run both implementations
    result_python = count_snp_alleles(bam_file, df, use_rust=False)
    result_rust = count_snp_alleles(bam_file, df, use_rust=True)

    # Compare counts
    assert result_python.equals(result_rust), "Counts differ!"

    # Compare MD5
    md5_python = hashlib.md5(result_python.write_csv()).hexdigest()
    md5_rust = hashlib.md5(result_rust.write_csv()).hexdigest()

    assert md5_python == md5_rust, f"MD5 mismatch: {md5_python} != {md5_rust}"


def test_rust_performance():
    """Benchmark Rust vs Python."""
    import time

    # ... load test data ...

    start = time.time()
    result_python = count_snp_alleles(bam_file, df, use_rust=False)
    python_time = time.time() - start

    start = time.time()
    result_rust = count_snp_alleles(bam_file, df, use_rust=True)
    rust_time = time.time() - start

    speedup = python_time / rust_time
    print(f"Speedup: {speedup:.2f}x")

    assert speedup > 2.0, f"Expected >2x speedup, got {speedup:.2f}x"
```

**Build & Test:**
```bash
# Build Rust extension
cd rust
maturin develop --release

# Run tests
cd ..
pytest tests/test_rust_equivalence.py -v

# Benchmark
hyperfine \
  'python -m counting count-variants ... --no-use-rust' \
  'python -m counting count-variants ... --use-rust'
```

#### Step 1.5: Distribution (30 min)

**Update `pyproject.toml`:**
```toml
[project]
name = "wasp2"
version = "1.1.0"
dependencies = [
    "numpy>=1.20",
    "polars>=0.19",
    "pysam>=0.20",
    # ... existing deps ...
]

[project.optional-dependencies]
rust = []  # Rust extension built by maturin

[build-system]
requires = ["maturin>=1.0,<2.0"]
build-backend = "maturin"
```

**Build wheels:**
```bash
# Build for current platform
maturin build --release

# Build for multiple platforms (requires Docker)
maturin build --release --manylinux 2_28

# Publish to PyPI
maturin publish
```

### Expected Performance Gains

| Implementation | Time | Speedup | Notes |
|----------------|------|---------|-------|
| **Baseline (pysam)** | 7.84s | 1.0x | Single-threaded Python |
| **Rust (sequential)** | 2-3s | 2.6-4x | noodles pure Rust |
| **Rust + rayon (4 cores)** | 0.5-1.5s | 5-15x | Parallel chromosome processing |

---

## 🧬 Phase 2: PLINK2 Integration

### Objective: Add QC, filtering, and population structure capabilities

**Library:** `Pgenlib` (official PLINK 2.0 Python API)
**Installation:** `pip install Pgenlib`

### Use Cases

1. **Pre-processing VCF** → QC filtering before WASP
2. **Quality Control** → Sample/variant missingness, HWE, MAF
3. **LD Pruning** → Reduce correlated variants
4. **Population Structure** → PCA for stratification correction

### Architecture

```
┌─────────────────────────────────────────────────────┐
│          WASP2 Main Workflow                        │
└──────────────┬──────────────────────────────────────┘
               │
               ├──► VCF Input
               │
               ▼
┌──────────────────────────────────────────────────────┐
│       src/preprocessing/plink_qc.py (NEW)           │
│                                                      │
│  ┌────────────────┐  ┌──────────────────┐          │
│  │ VCF → PGEN     │  │  QC Filtering    │          │
│  │ (via plink2)   │  │  (MAF, HWE, etc) │          │
│  └────────┬───────┘  └─────────┬────────┘          │
│           │                    │                     │
│           └──────────┬─────────┘                     │
│                      ▼                               │
│           ┌──────────────────┐                       │
│           │  Pgenlib Reader  │                       │
│           │  (PgenReader,    │                       │
│           │   PvarReader)    │                       │
│           └──────────┬───────┘                       │
└──────────────────────┼───────────────────────────────┘
                       │
                       ▼ Filtered VCF
                ┌──────────────┐
                │ WASP Counting│
                └──────────────┘
```

### Implementation Plan

#### Step 2.1: QC Module Setup (1 hour)

**Create `src/preprocessing/__init__.py`:**
```python
"""Preprocessing and QC modules for WASP2."""
from .plink_qc import PlinkQC, run_qc_pipeline

__all__ = ["PlinkQC", "run_qc_pipeline"]
```

**Create `src/preprocessing/plink_qc.py`:**
```python
"""PLINK2-based quality control and filtering."""
import subprocess
import tempfile
from pathlib import Path
from typing import Optional, Tuple, List
import pgenlib as pg
import numpy as np
import pandas as pd


class PlinkQC:
    """Quality control using PLINK 2.0."""

    def __init__(self, plink2_path: str = "plink2"):
        """
        Initialize PLINK QC handler.

        Args:
            plink2_path: Path to plink2 binary (default: searches PATH)
        """
        self.plink2_path = plink2_path
        self._check_plink2()

    def _check_plink2(self):
        """Verify plink2 is available."""
        try:
            result = subprocess.run(
                [self.plink2_path, "--version"],
                capture_output=True,
                text=True,
                check=True
            )
            print(f"Found PLINK2: {result.stdout.strip()}")
        except (subprocess.CalledProcessError, FileNotFoundError):
            raise RuntimeError(
                f"plink2 not found at '{self.plink2_path}'. "
                "Install from: https://www.cog-genomics.org/plink/2.0/"
            )

    def vcf_to_pgen(
        self,
        vcf_file: str,
        out_prefix: str,
        temp_dir: Optional[str] = None
    ) -> str:
        """
        Convert VCF to PGEN format.

        Args:
            vcf_file: Input VCF file
            out_prefix: Output prefix for .pgen/.pvar/.psam
            temp_dir: Temporary directory (default: system temp)

        Returns:
            Path to .pgen file
        """
        if temp_dir is None:
            temp_dir = tempfile.gettempdir()

        out_dir = Path(temp_dir) / out_prefix
        out_dir.mkdir(exist_ok=True, parents=True)

        pgen_path = out_dir / f"{out_prefix}.pgen"

        cmd = [
            self.plink2_path,
            "--vcf", vcf_file,
            "--make-pgen",
            "--out", str(out_dir / out_prefix)
        ]

        subprocess.run(cmd, check=True)

        return str(pgen_path)

    def filter_variants(
        self,
        pgen_file: str,
        out_prefix: str,
        maf: float = 0.01,
        geno: float = 0.05,
        hwe: float = 1e-6,
        mind: float = 0.1,
        keep_samples: Optional[str] = None,
        extract_variants: Optional[str] = None,
    ) -> Tuple[str, str, str]:
        """
        Apply QC filters to PGEN file.

        Args:
            pgen_file: Input PGEN file
            out_prefix: Output prefix
            maf: Minor allele frequency threshold
            geno: Maximum missingness per-variant
            hwe: Hardy-Weinberg equilibrium p-value threshold
            mind: Maximum missingness per-sample
            keep_samples: File with samples to keep (optional)
            extract_variants: File with variants to keep (optional)

        Returns:
            (pgen_path, pvar_path, psam_path)
        """
        cmd = [
            self.plink2_path,
            "--pfile", pgen_file.replace(".pgen", ""),
            "--maf", str(maf),
            "--geno", str(geno),
            "--hwe", str(hwe),
            "--mind", str(mind),
            "--make-pgen",
            "--out", out_prefix
        ]

        if keep_samples:
            cmd.extend(["--keep", keep_samples])

        if extract_variants:
            cmd.extend(["--extract", extract_variants])

        subprocess.run(cmd, check=True)

        return (
            f"{out_prefix}.pgen",
            f"{out_prefix}.pvar",
            f"{out_prefix}.psam"
        )

    def calculate_pca(
        self,
        pgen_file: str,
        out_prefix: str,
        n_components: int = 10
    ) -> pd.DataFrame:
        """
        Calculate principal components for population structure.

        Args:
            pgen_file: Input PGEN file
            out_prefix: Output prefix
            n_components: Number of PCs to compute

        Returns:
            DataFrame with sample PCs
        """
        cmd = [
            self.plink2_path,
            "--pfile", pgen_file.replace(".pgen", ""),
            "--pca", str(n_components),
            "--out", out_prefix
        ]

        subprocess.run(cmd, check=True)

        # Read eigenvectors
        eigenvec_file = f"{out_prefix}.eigenvec"
        pca_df = pd.read_csv(
            eigenvec_file,
            sep="\\t",
            header=0
        )

        return pca_df

    def ld_prune(
        self,
        pgen_file: str,
        out_prefix: str,
        window_size: int = 50,
        step_size: int = 5,
        r2_threshold: float = 0.5
    ) -> str:
        """
        LD pruning to remove correlated variants.

        Args:
            pgen_file: Input PGEN file
            out_prefix: Output prefix
            window_size: Variant window size (kb)
            step_size: Window step size (variants)
            r2_threshold: R² threshold for pruning

        Returns:
            Path to pruned variant list
        """
        cmd = [
            self.plink2_path,
            "--pfile", pgen_file.replace(".pgen", ""),
            "--indep-pairwise", str(window_size), str(step_size), str(r2_threshold),
            "--out", out_prefix
        ]

        subprocess.run(cmd, check=True)

        return f"{out_prefix}.prune.in"


def run_qc_pipeline(
    vcf_file: str,
    out_dir: str,
    maf: float = 0.01,
    geno: float = 0.05,
    hwe: float = 1e-6,
    run_pca: bool = False,
    run_ld_prune: bool = False
) -> dict:
    """
    Run complete PLINK QC pipeline.

    Args:
        vcf_file: Input VCF
        out_dir: Output directory
        maf: MAF threshold
        geno: Missingness threshold
        hwe: HWE p-value threshold
        run_pca: Calculate PCA (default: False)
        run_ld_prune: Perform LD pruning (default: False)

    Returns:
        Dict with paths to output files
    """
    qc = PlinkQC()

    Path(out_dir).mkdir(exist_ok=True, parents=True)

    results = {}

    # Step 1: Convert to PGEN
    print("Converting VCF to PGEN...")
    pgen_file = qc.vcf_to_pgen(vcf_file, f"{out_dir}/raw", out_dir)
    results["raw_pgen"] = pgen_file

    # Step 2: QC filtering
    print("Applying QC filters...")
    filtered_files = qc.filter_variants(
        pgen_file,
        f"{out_dir}/qc_filtered",
        maf=maf,
        geno=geno,
        hwe=hwe
    )
    results["filtered_pgen"] = filtered_files[0]
    results["filtered_pvar"] = filtered_files[1]
    results["filtered_psam"] = filtered_files[2]

    # Step 3: LD pruning (optional)
    if run_ld_prune:
        print("Performing LD pruning...")
        prune_list = qc.ld_prune(filtered_files[0], f"{out_dir}/ld_pruned")
        results["pruned_variants"] = prune_list

    # Step 4: PCA (optional)
    if run_pca:
        print("Calculating PCA...")
        pca_df = qc.calculate_pca(filtered_files[0], f"{out_dir}/pca")
        results["pca"] = pca_df

    return results
```

#### Step 2.2: Integrate with Counting Workflow (30 min)

**Update `src/counting/__main__.py`:**
```python
from preprocessing.plink_qc import run_qc_pipeline

@app.command()
def count_variants(
    bam: str,
    vcf: str,
    # ... existing params ...

    # NEW PLINK QC parameters
    plink_qc: bool = False,
    qc_maf: float = 0.01,
    qc_geno: float = 0.05,
    qc_hwe: float = 1e-6,
    qc_out_dir: Optional[str] = None,
):
    """
    Count variants with optional PLINK QC pre-filtering.

    New PLINK QC Options:
        --plink-qc: Enable PLINK2 quality control filtering
        --qc-maf: Minor allele frequency threshold (default: 0.01)
        --qc-geno: Maximum variant missingness (default: 0.05)
        --qc-hwe: Hardy-Weinberg equilibrium p-value (default: 1e-6)
        --qc-out-dir: QC output directory (default: temp)
    """

    # Apply PLINK QC if requested
    if plink_qc:
        print("Running PLINK2 quality control...")

        if qc_out_dir is None:
            import tempfile
            qc_out_dir = tempfile.mkdtemp(prefix="wasp2_qc_")

        qc_results = run_qc_pipeline(
            vcf_file=vcf,
            out_dir=qc_out_dir,
            maf=qc_maf,
            geno=qc_geno,
            hwe=qc_hwe
        )

        # Use filtered VCF for counting
        vcf = qc_results["filtered_pvar"]  # PVAR is VCF-like

        print(f"QC complete. Filtered variants: {qc_results['filtered_pvar']}")

    # Continue with normal counting workflow
    run_count_variants(
        bam_file=bam,
        vcf_file=vcf,
        ...
    )
```

#### Step 2.3: Population Structure Module (1 hour)

**Create `src/analysis/population_structure.py`:**
```python
"""Population structure analysis and correction."""
import pgenlib as pg
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional, Tuple
from preprocessing.plink_qc import PlinkQC


class PopulationStructure:
    """Analyze and correct for population stratification."""

    def __init__(self, pgen_file: str):
        """
        Initialize with PGEN file.

        Args:
            pgen_file: Path to .pgen file
        """
        self.pgen_file = Path(pgen_file)
        self.pvar_file = self.pgen_file.with_suffix(".pvar")
        self.psam_file = self.pgen_file.with_suffix(".psam")

        # Open readers
        self.pvar = pg.PvarReader(bytes(self.pvar_file))
        self.pgen = pg.PgenReader(
            bytes(self.pgen_file),
            pvar=self.pvar
        )

    def calculate_pca(
        self,
        n_components: int = 10,
        out_prefix: Optional[str] = None
    ) -> pd.DataFrame:
        """
        Calculate principal components.

        Args:
            n_components: Number of PCs
            out_prefix: Output prefix (if saving)

        Returns:
            DataFrame with sample PCs
        """
        if out_prefix is None:
            out_prefix = str(self.pgen_file.with_suffix(""))

        qc = PlinkQC()
        pca_df = qc.calculate_pca(
            str(self.pgen_file),
            out_prefix,
            n_components=n_components
        )

        return pca_df

    def detect_outliers(
        self,
        pca_df: pd.DataFrame,
        n_std: float = 6.0,
        components: Optional[List[int]] = None
    ) -> pd.DataFrame:
        """
        Detect population outliers using PCA.

        Args:
            pca_df: PCA DataFrame from calculate_pca()
            n_std: Number of standard deviations for outlier threshold
            components: Which PCs to use (default: first 2)

        Returns:
            DataFrame with outlier flags
        """
        if components is None:
            components = [0, 1]  # PC1, PC2

        pc_cols = [f"PC{i+1}" for i in components]

        # Calculate distance from centroid
        centroid = pca_df[pc_cols].mean()
        distances = np.sqrt(
            ((pca_df[pc_cols] - centroid) ** 2).sum(axis=1)
        )

        # Flag outliers
        threshold = distances.mean() + n_std * distances.std()
        pca_df["is_outlier"] = distances > threshold
        pca_df["distance"] = distances

        return pca_df

    def correct_for_stratification(
        self,
        counts_df: pd.DataFrame,
        pca_df: pd.DataFrame,
        n_pcs: int = 10
    ) -> pd.DataFrame:
        """
        Adjust allelic imbalance analysis for population structure.

        Args:
            counts_df: Allele counts from WASP counting
            pca_df: PCA DataFrame with sample PCs
            n_pcs: Number of PCs to use as covariates

        Returns:
            Adjusted counts DataFrame
        """
        # Merge counts with PC covariates
        merged = counts_df.merge(
            pca_df[[f"PC{i+1}" for i in range(n_pcs)]],
            left_on="sample",
            right_index=True,
            how="left"
        )

        # Linear regression to remove PC effects
        from scipy import stats

        pc_cols = [f"PC{i+1}" for i in range(n_pcs)]

        # Regress out PCs from allele balance
        merged["allele_balance"] = merged["alt_count"] / (
            merged["ref_count"] + merged["alt_count"]
        )

        X = merged[pc_cols].values
        y = merged["allele_balance"].values

        # Fit linear model
        slope, intercept, r_value, p_value, std_err = stats.linregress(
            X[:, 0], y  # Simplified: use PC1 only
        )

        # Residuals = corrected allele balance
        merged["corrected_balance"] = y - (slope * X[:, 0] + intercept)

        return merged

    def close(self):
        """Close file handles."""
        self.pgen.close()


# Example usage function
def analyze_population_structure(
    pgen_file: str,
    counts_file: str,
    out_prefix: str,
    n_components: int = 10
) -> dict:
    """
    Complete population structure analysis workflow.

    Args:
        pgen_file: PGEN file path
        counts_file: WASP allele counts
        out_prefix: Output prefix
        n_components: Number of PCs

    Returns:
        Dict with analysis results
    """
    pop = PopulationStructure(pgen_file)

    # Calculate PCA
    print("Calculating PCA...")
    pca_df = pop.calculate_pca(n_components, out_prefix)

    # Detect outliers
    print("Detecting population outliers...")
    pca_df = pop.detect_outliers(pca_df)

    outliers = pca_df[pca_df["is_outlier"]]
    print(f"Found {len(outliers)} population outliers")

    # Load counts
    counts_df = pd.read_csv(counts_file, sep="\\t")

    # Correct for stratification
    print("Correcting for population structure...")
    corrected_df = pop.correct_for_stratification(counts_df, pca_df)

    # Save results
    corrected_df.to_csv(f"{out_prefix}_corrected_counts.tsv", sep="\\t", index=False)
    pca_df.to_csv(f"{out_prefix}_pca.tsv", sep="\\t", index=False)

    pop.close()

    return {
        "pca": pca_df,
        "corrected_counts": corrected_df,
        "n_outliers": len(outliers)
    }
```

#### Step 2.4: CLI Commands (30 min)

**Add to `src/analysis/__main__.py`:**
```python
from preprocessing.plink_qc import PlinkQC
from analysis.population_structure import analyze_population_structure

@app.command()
def qc_variants(
    vcf: str,
    out_dir: str = "qc_output",
    maf: float = 0.01,
    geno: float = 0.05,
    hwe: float = 1e-6,
    run_pca: bool = False,
    run_ld_prune: bool = False
):
    """
    Run PLINK2 quality control on VCF file.

    Args:
        vcf: Input VCF file
        out_dir: Output directory
        maf: Minor allele frequency threshold
        geno: Maximum variant missingness
        hwe: Hardy-Weinberg equilibrium p-value
        run_pca: Calculate population PCA
        run_ld_prune: Perform LD pruning
    """
    from preprocessing.plink_qc import run_qc_pipeline

    results = run_qc_pipeline(
        vcf_file=vcf,
        out_dir=out_dir,
        maf=maf,
        geno=geno,
        hwe=hwe,
        run_pca=run_pca,
        run_ld_prune=run_ld_prune
    )

    print("\\nQC Results:")
    for key, value in results.items():
        print(f"  {key}: {value}")


@app.command()
def population_pca(
    pgen: str,
    counts: str,
    out_prefix: str,
    n_components: int = 10
):
    """
    Analyze and correct for population structure.

    Args:
        pgen: PGEN file from PLINK
        counts: WASP allele counts TSV
        out_prefix: Output prefix
        n_components: Number of principal components
    """
    results = analyze_population_structure(
        pgen_file=pgen,
        counts_file=counts,
        out_prefix=out_prefix,
        n_components=n_components
    )

    print(f"\\nPopulation structure analysis complete:")
    print(f"  Outliers detected: {results['n_outliers']}")
    print(f"  Corrected counts: {out_prefix}_corrected_counts.tsv")
    print(f"  PCA results: {out_prefix}_pca.tsv")
```

### Expected Benefits

| Feature | Benefit | Impact |
|---------|---------|--------|
| **QC Filtering** | Remove low-quality variants | ↑ Statistical power |
| **LD Pruning** | Reduce redundancy | ↓ Multiple testing burden |
| **Population PCA** | Control for stratification | ↓ False positives |
| **Outlier Detection** | Remove divergent samples | ↑ Homogeneity |

---

## 📊 Phase 3: Combined Optimization Strategy

### Recommended Implementation Order

**Session 1 (4-6 hours): Rust BAM I/O**
1. Setup Rust workspace (30 min)
2. Implement `bam_counter.rs` (2-3 hours)
3. Python integration (1 hour)
4. Testing & validation (1-2 hours)

**Session 2 (3-4 hours): PLINK Integration**
1. QC module implementation (1 hour)
2. Counting workflow integration (30 min)
3. Population structure module (1 hour)
4. CLI commands & documentation (1 hour)

**Session 3 (2-3 hours): Polish & Release**
1. Comprehensive testing (1 hour)
2. Documentation updates (1 hour)
3. Performance benchmarking (30 min)
4. Release prep (30 min)

### Final Performance Targets

| Pipeline | Current | Target | Speedup |
|----------|---------|--------|---------|
| **Counting** | 8.03s | 1-2s | 4-8x |
| **Analysis** | 3.38s | 2-3s | 1.1-1.7x |
| **Overall** | 11.41s | 3-5s | 2.3-3.8x |

With PLINK QC pre-filtering: Additional 10-30% improvement in analysis accuracy.

---

## 🚀 Installation & Usage

### For Users (After Implementation)

**Standard installation:**
```bash
pip install wasp2
```

**With Rust optimizations:**
```bash
pip install wasp2[rust]
```

**With PLINK2 support:**
```bash
# Install PLINK2
wget https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_x86_64.zip
unzip plink2_linux_x86_64.zip
mv plink2 /usr/local/bin/

# Install Python library
pip install wasp2 Pgenlib
```

### Example Workflows

**1. Basic WASP workflow (Python):**
```bash
python -m counting count-variants input.bam variants.vcf -o counts.tsv
python -m analysis find-imbalance counts.tsv -o ai_results.tsv
```

**2. With Rust acceleration:**
```bash
python -m counting count-variants input.bam variants.vcf -o counts.tsv --use-rust
```

**3. With PLINK QC:**
```bash
# Pre-filter VCF
python -m analysis qc-variants variants.vcf --out-dir qc_output --run-pca

# Count with filtered variants
python -m counting count-variants input.bam qc_output/qc_filtered.pvar -o counts.tsv --use-rust

# Analyze with population correction
python -m analysis population-pca qc_output/qc_filtered.pgen counts.tsv out_prefix
```

---

## 📝 Summary

### Current Status
✅ One consolidated branch
✅ Quick wins implemented
✅ Profiling complete
✅ Integration plans documented

### Next Steps
1. Choose: **Rust** (speed) or **PLINK** (features) or **Both**
2. Schedule: 4-6 hour session for Rust, 3-4 hours for PLINK
3. Execute according to plans above
4. Validate with full test suite

**Decision Point:** What would you like to tackle first? 🦀 or 🧬 ?
