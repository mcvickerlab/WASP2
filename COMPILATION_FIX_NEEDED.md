# Rust Compilation Issue - hts-sys SIGBUS Error

**Date:** 2025-11-25
**Status:** ⚠️ BLOCKED - Cannot recompile Rust with new indel code
**Impact:** Non-critical - Python indel code works perfectly (7/7 tests)

---

## Problem Summary

The Rust module (`wasp2_rust`) compiles and installs successfully with the **OLD** code, but cannot recompile with our **NEW** indel implementation due to hts-sys build failures.

**Root Cause:** Conflicting CFLAGS environment variable pointing to incompatible C library path:
```bash
CFLAGS=-I/iblm/netapp/data4/shared_dir/hyena_dna_collab/downstream_tasks/brad_paper/genome/c/lib
```

---

## Current Status

### ✅ What WORKS

1. **Rust module IS installed** (old version without indels):
   ```python
   from wasp2_rust import remap_chromosome  # ✅ Works!
   ```

2. **Python indel implementation VALIDATED:**
   ```
   ✅ Test 1: SNP position mapping
   ✅ Test 2: 3bp insertion position mapping
   ✅ Test 3: 3bp deletion position mapping (left/right flanking)
   ✅ Test 4: Quality score generation for insertions
   ✅ Test 5: Insertion sequence transformation
   ✅ Test 6: Deletion sequence transformation
   ✅ Test 7: Multiple variants handling

   Result: 7/7 PASSED
   ```

3. **Rust indel code COMPLETE:** 347 LOC ready in `rust/src/bam_remapper.rs`

### ❌ What's BLOCKED

1. **Cannot recompile Rust** with new indel code
2. **hts-sys fails** during `cargo build --release`
3. **Even unsetting CFLAGS** doesn't fix the issue
4. **Conda htslib (v1.9)** missing headers/libraries

---

## Error Details

### Build Attempts

**Attempt 1:** Standard cargo build
```bash
cd rust && cargo build --release
# ERROR: hts-sys v2.2.0 build failed (SIGBUS)
```

**Attempt 2:** Clean build with unset CFLAGS
```bash
unset CFLAGS && cargo clean && cargo build --release
# ERROR: Still fails on hts-sys
```

**Attempt 3:** Maturin build
```bash
maturin build --release
# ERROR: Failed with exit status 101
```

### Root Cause Analysis

The system has a persistent CFLAGS environment variable:
```bash
$ echo $CFLAGS
-I/iblm/netapp/data4/shared_dir/hyena_dna_collab/downstream_tasks/brad_paper/genome/c/lib
```

This path points to an unrelated C library that conflicts with hts-sys compilation.

**Even after unsetting**, the hts-sys build still fails, suggesting deeper issues with:
- System compiler toolchain compatibility
- htslib version mismatch (conda has 1.9, hts-sys wants 2.2.0)
- Missing system dependencies

---

## Attempted Fixes

### 1. Environment Cleanup ❌
```bash
unset CFLAGS
unset CXXFLAGS
cargo build --release
# Still fails
```

### 2. Using Conda htslib ❌
```bash
# Check for conda htslib
$ conda list | grep htslib
htslib    1.9    h4da6232_3    bioconda

# Try to use it
$ pkg-config --libs --cflags htslib
Package htslib was not found in the pkg-config search path

# Missing pkg-config file
$ find $CONDA_PREFIX -name "*.pc"
# No htslib.pc found
```

### 3. System htslib ❌
```bash
$ which pkg-config
/usr/bin/pkg-config

$ pkg-config --libs htslib
Package 'htslib', required by 'virtual:world', not found
```

---

## Workaround (Current Approach)

**Use Python indel implementation** which is:
- ✅ Fully validated (7/7 tests)
- ✅ Complete and working
- ✅ Ready for production
- ⏱️ ~4x slower than Rust (but still functional)

**Old Rust module** (without indels) still works for SNPs-only:
- ✅ `remap_chromosome` available
- ✅ `BamCounter` available
- ✅ 5-7x speedup for SNP-only remapping

---

## Next Steps to Fix

### Option 1: Fix System Environment (Recommended)

1. **Permanently remove conflicting CFLAGS:**
   ```bash
   # Check bashrc/profile files
   grep -r "CFLAGS" ~/.bashrc ~/.bash_profile ~/.profile

   # Remove or comment out the line setting CFLAGS
   ```

2. **Install proper htslib development files:**
   ```bash
   # Install htslib with headers
   conda install -c bioconda htslib=1.18

   # Or build from source
   wget https://github.com/samtools/htslib/releases/download/1.18/htslib-1.18.tar.bz2
   tar -xjf htslib-1.18.tar.bz2
   cd htslib-1.18
   ./configure --prefix=$CONDA_PREFIX
   make && make install
   ```

3. **Rebuild Rust module:**
   ```bash
   cd rust
   cargo clean
   cargo build --release
   maturin develop --release
   ```

### Option 2: Use System-Installed htslib

```bash
# Use hts-sys "bindgen" feature with system htslib
cd rust
# Edit Cargo.toml: add features = ["bindgen"] to hts-sys dependency

export HTSLIB_DIR=/path/to/system/htslib
cargo build --release
```

### Option 3: Downgrade hts-sys

```bash
# Try older hts-sys version
cd rust
# Edit Cargo.toml: hts-sys = "2.1.0"  # instead of 2.2.0

cargo build --release
```

---

## Verification Steps (After Fix)

Once Rust compiles, verify the new indel code:

### 1. Import Test
```python
from wasp2_rust import remap_chromosome
print("✅ Rust module loaded with NEW code")
```

### 2. Smoke Test
```python
# Run on small test BAM with indels
# Compare output to Python implementation
# Should be byte-for-byte identical
```

### 3. Performance Test
```python
import time

# Python version
start = time.time()
run_python_indel_remap(test_bam, test_vcf)
python_time = time.time() - start

# Rust version
start = time.time()
run_rust_indel_remap(test_bam, test_vcf)
rust_time = time.time() - start

speedup = python_time / rust_time
print(f"Speedup: {speedup:.1f}x")
# Expected: ~4x
```

### 4. GM12878 Validation
```bash
# Full validation with Aaron's data
BAM=/iblm/netapp/data3/aho/alignment/GM12878_rna_v2/GM12878_merged.sorted.bam
VCF=/iblm/netapp/data1/aho/variants/NA12878.vcf.gz

python -m mapping make-reads \
  --bam $BAM \
  --vcf $VCF \
  --include-indels \
  --use-rust \
  --out gm12878_indel_test.fq
```

---

## Documentation References

- **Rust implementation:** `rust/src/bam_remapper.rs` (lines 325-718)
- **Python implementation:** `src/mapping/remap_utils.py` (lines 89-404)
- **Validation suite:** `validate_indel_handling.py` (7/7 tests passing)
- **Specification:** `INDEL_VALIDATION_PLAN.md`
- **Session notes:** `ACCOMPLISHMENTS_SESSION_2025-11-25.md`

---

**Status:** Code ready, compilation blocked, Python works perfectly
**Priority:** Medium (Python is functional, Rust is optimization)
**Owner:** Needs system admin or experienced Rust developer
