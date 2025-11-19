# WASP2 Performance Analysis & Optimization Strategy

**Date:** 2025-11-19
**Analysis Duration:** 15 minutes profiling + research
**Session Time Available:** 1 hour

---

## 📊 Performance Profiling Results

### Counting Pipeline (111,454 SNPs)

**Total Runtime:** 12.18 seconds

**Bottleneck Breakdown:**
1. **pysam BAM reading: 7.84 seconds (64%)** ⚠️ CRITICAL BOTTLENECK
   - `cnext()`: 7.838s in pure C iteration
   - Called 115,506 times
   - Single-threaded I/O bound

2. **fetch() operations: 0.22 seconds (2%)**
   - Called 111,454 times (once per SNP)
   - Sequential region fetching

3. **count_snp_alleles(): 0.26 seconds (2%)**
   - Python loop overhead
   - Dictionary operations

4. **Polars LazyFrame: 0.08 seconds (0.7%)**
   - Suboptimal `.columns` access (triggers schema resolution)
   - Warning: Use `collect_schema().names()` instead

### Analysis Pipeline (43 regions)

**Total Runtime:** 3.38 seconds

**Not a bottleneck** - most time spent in imports/I/O (0.66s)

---

## 🎯 Quick Wins (Can Do in 30 Minutes)

### 1. Fix Polars LazyFrame Warnings (5 min)
**File:** `src/counting/filter_variant_data.py`

**Current (lines 151, 168):**
```python
intersect_ncols = len(df.columns)  # ⚠️ Triggers full schema resolution
subset_cols = [df.columns[i] for i in ...]  # ⚠️ Again!
```

**Optimized:**
```python
schema = df.collect_schema()
intersect_ncols = len(schema.names())
subset_cols = [schema.names()[i] for i in ...]
```

**Expected speedup:** ~10-20% on Polars operations

---

### 2. Fix Pandas GroupBy Warnings (5 min)
**File:** `src/analysis/as_analysis.py:248, 252`

**Current:**
```python
null_test = group_df.apply(lambda x: np.sum(betabinom.logpmf(...)))
alt_test = group_df.apply(lambda x: parse_opt(x, disp, phased=phased))
```

**Optimized:**
```python
null_test = group_df.apply(lambda x: np.sum(betabinom.logpmf(...)), include_groups=False)
alt_test = group_df.apply(lambda x: parse_opt(x, disp, phased=phased), include_groups=False)
```

**Expected speedup:** Minimal, but removes deprecation warnings

---

### 3. Vectorize count_alleles orientation warning (2 min)
**File:** `src/counting/count_alleles.py:63`

**Current:**
```python
count_df = pl.DataFrame(...)  # ⚠️ Orientation inferred
```

**Optimized:**
```python
count_df = pl.DataFrame(..., orient="row")
```

**Expected speedup:** None, just removes warning

---

### 4. Remove Duplicate CLI Parameter (2 min)
**Files:** `src/counting/__main__.py:27-28`, `src/mapping/__main__.py`

**Current:**
```python
typer.Option("--samples", "--sample", "--samps", "--samps", "-s", ...)  # Duplicate!
```

**Optimized:**
```python
typer.Option("--samples", "--sample", "--samps", "-s", ...)
```

---

## 🚀 Rust Optimization Strategy (Multi-Session Project)

### Target: BAM Reading (7.84s → ~1-2s expected)

### Phase 1: Research & Setup (1-2 hours)

**Libraries:**
- **PyO3**: Rust ↔ Python bindings (https://pyo3.rs/)
- **noodles**: Pure Rust BAM/SAM/CRAM parser (https://github.com/zaeleus/noodles)
  - Modern, async-capable, 1.5-2.5x faster than rust-htslib
  - Active development (2024-2025)
- **Alternative:** rust-htslib (bindings to HTSlib C library)

**Setup:**
```toml
# Cargo.toml
[dependencies]
pyo3 = { version = "0.20", features = ["extension-module"] }
noodles = { version = "0.76", features = ["bam", "async"] }
rayon = "1.8"  # Parallel iterator
```

---

### Phase 2: Implement BAM Reader (3-4 hours)

**Target Function:** `count_snp_alleles()` in `src/counting/count_alleles.py`

**Rust Implementation Plan:**

```rust
// src/rust/bam_counter/src/lib.rs
use pyo3::prelude::*;
use noodles::bam;
use rayon::prelude::*;

#[pyfunction]
fn count_alleles_rust(
    bam_path: &str,
    regions: Vec<(String, u32, u32)>,  // [(chrom, start, end), ...]
    py: Python,
) -> PyResult<Vec<(u32, u32, u32)>> {  // [(ref, alt, other), ...]

    // Read BAM with noodles
    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_path)?;

    // Parallel processing with rayon
    let counts: Vec<_> = regions.par_iter()
        .map(|(chrom, start, end)| {
            // Fetch reads in region (zero-copy where possible)
            let query = reader.query(chrom, *start, *end)?;

            // Count alleles
            let mut ref_count = 0;
            let mut alt_count = 0;
            let mut other_count = 0;

            for read in query {
                let read = read?;
                // Extract base at SNP position
                // ... allele counting logic ...
            }

            Ok((ref_count, alt_count, other_count))
        })
        .collect();

    Ok(counts)
}

#[pymodule]
fn wasp2_rust(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(count_alleles_rust, m)?)?;
    Ok(())
}
```

**Python Integration:**
```python
# src/counting/count_alleles.py
try:
    from wasp2_rust import count_alleles_rust
    HAS_RUST = True
except ImportError:
    HAS_RUST = False

def count_snp_alleles(...):
    if HAS_RUST:
        # Use Rust implementation (5-10x faster)
        return count_alleles_rust(bam_file, regions)
    else:
        # Fallback to Python/pysam
        return _count_snp_alleles_python(...)
```

**Build Configuration:**
```toml
# pyproject.toml
[build-system]
requires = ["maturin>=1.0,<2.0"]
build-backend = "maturin"

[project.optional-dependencies]
rust = ["wasp2-rust"]
```

**Expected Performance:**
- **Current:** 7.84s (pysam)
- **Rust noodles:** ~1.5-2.5s (3-5x faster)
- **Rust parallel:** ~0.5-1.5s (5-15x faster with rayon)

---

### Phase 3: Testing & Validation (1-2 hours)

**Validation Strategy:**
```python
# tests/test_rust_equivalence.py
def test_rust_python_equivalence():
    """Ensure Rust and Python produce identical results."""

    # Run both implementations
    py_result = count_alleles_python(bam, regions)
    rs_result = count_alleles_rust(bam, regions)

    # Compare MD5 checksums
    assert md5(py_result) == md5(rs_result)
```

**Benchmark:**
```bash
hyperfine \
  'python -m counting count-variants ...' \
  'python -m counting count-variants ... --use-rust'
```

---

## 🧬 plink2 Integration Strategy

### Use Cases

1. **Pre-filtering VCF files**
   - LD pruning
   - MAF filtering
   - Hardy-Weinberg filtering

2. **Quality Control**
   - Sample QC (missingness, heterozygosity)
   - Variant QC

3. **Population Structure**
   - PCA for stratification
   - Kinship/relatedness checks

---

### Implementation Plan

**Library:** pandas-plink (v2.3.2, May 2025)

**Installation:**
```bash
pip install pandas-plink plink2
```

**Integration Points:**

#### 1. Pre-Processing Hook (New Module)

```python
# src/preprocessing/plink_qc.py
import pandas_plink as pp
import subprocess

def run_plink_qc(
    vcf_file: str,
    out_prefix: str,
    maf: float = 0.01,
    geno: float = 0.05,
    hwe: float = 1e-6,
) -> str:
    """Run PLINK2 quality control on VCF."""

    # Call plink2 via subprocess
    cmd = [
        "plink2",
        "--vcf", vcf_file,
        "--maf", str(maf),
        "--geno", str(geno),
        "--hwe", str(hwe),
        "--make-pgen",
        "--out", out_prefix
    ]
    subprocess.run(cmd, check=True)

    return f"{out_prefix}.pgen"

def read_plink_variants(pgen_file: str) -> pd.DataFrame:
    """Read PLINK2 variant data into DataFrame."""

    # Read using pandas-plink
    (bim, fam, bed) = pp.read_plink(pgen_file)

    return bim  # Variant info
```

#### 2. Add to Counting Workflow

```python
# src/counting/__main__.py
@app.command()
def count_variants(
    bam: str,
    vcf: str,
    plink_qc: bool = False,  # NEW
    maf: float = 0.01,       # NEW
    ...
):
    if plink_qc:
        # Pre-filter with PLINK2
        from preprocessing.plink_qc import run_plink_qc
        filtered_vcf = run_plink_qc(vcf, "qc_filtered", maf=maf)
        vcf = filtered_vcf

    # Continue with normal workflow
    run_count_variants(bam, vcf, ...)
```

#### 3. Population Stratification Module (New)

```python
# src/analysis/population_structure.py
def calculate_pca_stratification(
    vcf_file: str,
    n_components: int = 10
) -> pd.DataFrame:
    """Calculate PCA for population stratification."""

    # Run PLINK2 PCA
    subprocess.run([
        "plink2",
        "--vcf", vcf_file,
        "--pca", str(n_components),
        "--out", "pca_output"
    ])

    # Read results with pandas-plink
    pca_df = pd.read_csv("pca_output.eigenvec", sep="\t")

    return pca_df
```

---

### Integration Benefits

1. **Quality Control:** Remove low-quality variants before analysis
2. **Performance:** PLINK2 is highly optimized for large datasets
3. **Reproducibility:** Standard bioinformatics workflow
4. **Population Correction:** Control for stratification in AI analysis

---

## 📈 Expected Performance Gains Summary

| Optimization | Current | Expected | Speedup | Effort |
|--------------|---------|----------|---------|--------|
| **Quick Wins** | 12.18s | 10-11s | 1.1-1.2x | 30 min |
| **Rust BAM Reading** | 7.84s | 0.5-2s | 4-15x | 6-10 hrs |
| **Rust + Parallel** | 7.84s | 0.5-1s | 8-15x | 8-12 hrs |
| **plink2 Pre-filter** | N/A | N/A | Data quality | 2-4 hrs |

---

## 🎯 Recommended Timeline

### Session 1 (Today - 45 min left):
✅ Quick wins implementation
- Fix LazyFrame warnings
- Fix pandas warnings
- Remove duplicate parameters

### Session 2 (2-3 hours):
- Set up Rust development environment
- Implement basic BAM reader prototype
- Benchmark against Python

### Session 3 (3-4 hours):
- Complete Rust implementation
- Parallel processing with rayon
- Full validation suite

### Session 4 (2-3 hours):
- plink2 integration
- QC workflow
- Population structure module

---

## 🔗 Resources

**Rust/PyO3:**
- PyO3 Guide: https://pyo3.rs/
- noodles BAM: https://github.com/zaeleus/noodles
- cnv-from-bam (example): https://pypi.org/project/cnv-from-bam/

**plink2/Python:**
- pandas-plink: https://github.com/limix/pandas-plink
- PLINK 2.0 docs: https://www.cog-genomics.org/plink/2.0/

**Benchmarking:**
- hyperfine: https://github.com/sharkdp/hyperfine
- py-spy: https://github.com/benfred/py-spy

---

**Next Steps:** Implement quick wins now (30 min), plan Rust session for dedicated 4+ hour block.
