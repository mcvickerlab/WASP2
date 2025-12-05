# WASP2 Rust Threading Optimization Plan

**Date:** 2025-12-04
**Target:** Improve thread scaling efficiency across all WASP2 Rust modules
**Affected Modules:** 3 (bam_counter, bam_remapper, multi_sample)

---

## Executive Summary

Audit identified **3 modules** with missing BGZF decompression threading. All use the same
per-chromosome parallel pattern via Rayon but fail to enable htslib's internal thread pool
for BAM decompression. Fix is **3 one-line additions**.

| Module | Function | Line | Pipeline Stage | Impact |
|--------|----------|------|----------------|--------|
| `bam_counter.rs` | `process_chromosome_reads()` | 147 | Allele counting | HIGH |
| `bam_remapper.rs` | `swap_alleles_for_chrom()` | 352 | Read remapping | HIGH |
| `multi_sample.rs` | `swap_alleles_for_chrom_multi()` | 586 | Multi-sample | MEDIUM |

---

## 1. Architecture Overview

### 1.1 WASP2 Pipeline Flow
```
┌─────────────────────────────────────────────────────────────────────────────┐
│                         WASP2 Rust Pipeline                                  │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  Input BAM ──┬──► bam_filter.rs ──► Filter reads by variant overlap         │
│              │    (set_threads ✓)                                           │
│              │                                                              │
│              ├──► bam_intersect.rs ──► Find read-variant intersections      │
│              │    (set_threads ✓)                                           │
│              │                                                              │
│              ├──► bam_remapper.rs ──► Generate haplotype reads      ❌ FIX  │
│              │    (set_threads MISSING)                                     │
│              │                                                              │
│              ├──► bam_counter.rs ──► Count alleles at variants      ❌ FIX  │
│              │    (set_threads MISSING)                                     │
│              │                                                              │
│              └──► multi_sample.rs ──► Multi-sample processing       ❌ FIX  │
│                   (set_threads MISSING)                                     │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

### 1.2 Common Threading Pattern (All 3 Modules)
```rust
// Pattern used in all 3 modules:
fn parallel_processor() {
    let chromosomes = get_chromosomes();

    chromosomes.par_iter()  // Rayon parallel iterator
        .map(|chrom| {
            // Each thread opens its own BAM reader
            let mut bam = bam::IndexedReader::from_path(path)?;
            // ❌ MISSING: bam.set_threads(2).ok();

            bam.fetch(chrom)?;
            process_reads(&mut bam)
        })
        .collect()
}
```

### 1.3 Why This Matters
```
WITHOUT set_threads():
┌────────────────────────────────────────────────────────────┐
│  Rayon Thread 1: chr1  [decompress]─[process]─[decompress]─│ ◄─ Serial BGZF
│  Rayon Thread 2: chr2  [decompress]─[process]─[decompress]─│ ◄─ Serial BGZF
│  ...                                                        │
│  Rayon Thread 8: chrY  [decompress]─[process]              │ ◄─ Serial BGZF
└────────────────────────────────────────────────────────────┘
Bottleneck: BGZF decompression is ~30% of BAM read time

WITH set_threads(2):
┌────────────────────────────────────────────────────────────┐
│  Rayon Thread 1: chr1  [decomp║decomp]─[process]─...       │ ◄─ Parallel BGZF
│  Rayon Thread 2: chr2  [decomp║decomp]─[process]─...       │ ◄─ Parallel BGZF
│  ...                                                        │
└────────────────────────────────────────────────────────────┘
Result: 20-30% faster per-thread BAM reading
```

---

## 2. Module-by-Module Analysis

### 2.1 bam_counter.rs (Allele Counting)

**Purpose:** Count REF/ALT alleles at heterozygous SNP positions

**Current Code (lines 140-155):**
```rust
fn process_chromosome_reads(
    &self,
    chrom: &str,
    regions: &[(usize, Region)],
    min_qual: u8,
    debug_sites: &FxHashMap<(String, u32), usize>,
) -> PyResult<FxHashMap<usize, (u32, u32, u32)>> {
    let mut bam = bam::IndexedReader::from_path(&self.bam_path)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(
            format!("Failed to open BAM: {}", e)
        ))?;
    // ❌ MISSING: bam.set_threads(2).ok();

    let mut seen_reads: FxHashSet<Vec<u8>> = FxHashSet::default();
    // ... rest of function
}
```

**Called From:** `count_alleles_impl()` line 106 via `grouped.par_iter().map()`

**Benchmark Impact:**
- Current: 70% efficiency at 8 threads (5.6x speedup)
- Expected: 75-80% efficiency (6.0-6.4x speedup)

---

### 2.2 bam_remapper.rs (Read Remapping)

**Purpose:** Generate haplotype read sequences by swapping alleles

**Current Code (lines 346-365):**
```rust
pub fn swap_alleles_for_chrom(
    bam_path: &str,
    variants: &FxHashMap<Vec<u8>, Vec<VariantSpan>>,
    chrom: &str,
    config: &RemapConfig,
) -> Result<(Vec<HaplotypeRead>, RemapStats)> {
    let mut bam = bam::IndexedReader::from_path(bam_path)
        .context("Failed to open BAM file")?;
    // ❌ MISSING: bam.set_threads(2).ok();

    let mut results = Vec::new();
    let mut stats = RemapStats::default();

    // Fetch reads for this chromosome
    let header = bam.header().clone();
    let tid = header.tid(chrom.as_bytes())...
    bam.fetch(tid as i32)?;
    // ... rest of function
}
```

**Called From:**
- `process_all_chromosomes_parallel()` line 814 via `chromosomes.par_iter().map()`
- `process_and_write_parallel()` (streaming version)

**Performance Notes:**
- Documented as "7x faster than Python" in comments
- This is the main bottleneck in `make_remap_reads` step
- With fix: potentially 8-9x faster than Python

---

### 2.3 multi_sample.rs (Multi-Sample Processing)

**Purpose:** Process variants across multiple samples simultaneously

**Current Code (lines 576-602):**
```rust
pub fn swap_alleles_for_chrom_multi(
    bam_path: &str,
    variants: &FxHashMap<Vec<u8>, Vec<MultiSampleVariantSpan>>,
    chrom: &str,
    out_r1: &str,
    out_r2: &str,
    max_seqs: usize,
) -> Result<MultiSampleRemapStats> {
    use rustc_hash::FxHashMap;

    let mut bam = bam::IndexedReader::from_path(bam_path)
        .context("Failed to open BAM file")?;
    // ❌ MISSING: bam.set_threads(2).ok();

    let mut stats = MultiSampleRemapStats::default();

    // Get chromosome tid
    let header = bam.header().clone();
    let tid = match header.tid(chrom.as_bytes()) {
        Some(t) => t,
        None => { return Ok(stats); }
    };

    bam.fetch(tid as i32)?;
    // ... rest of function
}
```

**Called From:** Multi-sample parallel processing pipeline

---

## 3. Implementation Plan

### 3.1 Phase 1: Quick Fixes (3 one-liners)

**Estimated Time:** 10 minutes total

```diff
# Fix 1: bam_counter.rs (after line 150)
  let mut bam = bam::IndexedReader::from_path(&self.bam_path)
      .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(
          format!("Failed to open BAM: {}", e)
      ))?;
+ // Enable parallel BGZF decompression (2 threads per chromosome worker)
+ bam.set_threads(2).ok();

# Fix 2: bam_remapper.rs (after line 353)
  let mut bam = bam::IndexedReader::from_path(bam_path)
      .context("Failed to open BAM file")?;
+ // Enable parallel BGZF decompression (2 threads per chromosome worker)
+ bam.set_threads(2).ok();

# Fix 3: multi_sample.rs (after line 587)
  let mut bam = bam::IndexedReader::from_path(bam_path)
      .context("Failed to open BAM file")?;
+ // Enable parallel BGZF decompression (2 threads per chromosome worker)
+ bam.set_threads(2).ok();
```

### 3.2 Build and Test

```bash
# Build
cd rust
cargo build --release

# Quick smoke test
python -c "
from wasp2_rust import BamCounter
import time

counter = BamCounter('test.bam')
regions = [('chr1', 100, 'A', 'G'), ...]  # Load test regions

for threads in [1, 4, 8]:
    start = time.time()
    results = counter.count_alleles(regions, min_qual=20, threads=threads)
    print(f'{threads}t: {time.time()-start:.2f}s')
"
```

### 3.3 Phase 2: Advanced Optimizations (Optional)

If Phase 1 doesn't achieve >80% efficiency, consider:

1. **Variant-based chunking** (replaces per-chromosome)
2. **Shared thread pool** (prevents thread explosion)
3. **Single-pass streaming** (max efficiency)

See `BAMCOUNTER_THREADING_OPTIMIZATION.md` for detailed Level 2-4 plans.

---

## 4. Code Diffs

### 4.1 bam_counter.rs

```diff
--- a/rust/src/bam_counter.rs
+++ b/rust/src/bam_counter.rs
@@ -147,6 +147,9 @@ impl BamCounter {
             format!("Failed to open BAM: {}", e)
         ))?;

+        // Enable parallel BGZF decompression (2 threads per chromosome worker)
+        bam.set_threads(2).ok();
+
         let mut seen_reads: FxHashSet<Vec<u8>> = FxHashSet::default();
         let total_snps: usize = regions.len();
```

### 4.2 bam_remapper.rs

```diff
--- a/rust/src/bam_remapper.rs
+++ b/rust/src/bam_remapper.rs
@@ -352,6 +352,9 @@ pub fn swap_alleles_for_chrom(
     let mut bam = bam::IndexedReader::from_path(bam_path)
         .context("Failed to open BAM file")?;

+    // Enable parallel BGZF decompression (2 threads per chromosome worker)
+    bam.set_threads(2).ok();
+
     let mut results = Vec::new();
     let mut stats = RemapStats::default();
```

### 4.3 multi_sample.rs

```diff
--- a/rust/src/multi_sample.rs
+++ b/rust/src/multi_sample.rs
@@ -586,6 +586,9 @@ pub fn swap_alleles_for_chrom_multi(
     let mut bam = bam::IndexedReader::from_path(bam_path)
         .context("Failed to open BAM file")?;

+    // Enable parallel BGZF decompression (2 threads per chromosome worker)
+    bam.set_threads(2).ok();
+
     let mut stats = MultiSampleRemapStats::default();
```

---

## 5. Expected Results

### 5.1 BamCounter (Allele Counting)

| Threads | Before | After | Improvement |
|---------|--------|-------|-------------|
| 1 | 84.6s | ~70s | ~17% |
| 4 | 23.9s | ~20s | ~16% |
| 8 | 15.1s | ~12s | ~20% |
| **Efficiency** | **70%** | **~78%** | **+8%** |

### 5.2 BamRemapper (Read Remapping)

| Metric | Before | After |
|--------|--------|-------|
| Per-chrom time | ~0.15s | ~0.12s |
| Total 8-thread | Unknown | ~20% faster |
| vs Python | 7x faster | ~9x faster |

### 5.3 Overall Pipeline

```
WASP2 Full Pipeline (56M reads):

Before:
  Filter:    ~100s (set_threads ✓)
  Intersect: ~50s  (set_threads ✓)
  Remap:     ~200s (set_threads ❌)
  Count:     ~15s  (set_threads ❌)
  ─────────────────────────────────
  Total:     ~365s

After (estimated):
  Filter:    ~100s (no change)
  Intersect: ~50s  (no change)
  Remap:     ~160s (20% faster)  ✓ FIXED
  Count:     ~12s  (20% faster)  ✓ FIXED
  ─────────────────────────────────
  Total:     ~322s (~12% overall speedup)
```

---

## 6. Testing Plan

### 6.1 Unit Tests

```rust
#[cfg(test)]
mod threading_tests {
    use super::*;

    #[test]
    fn test_bam_counter_threads() {
        // Verify set_threads is called
        let counter = BamCounter::new("test.bam".into()).unwrap();
        // Run with 1 vs 4 threads, expect >3x speedup
    }

    #[test]
    fn test_remapper_threads() {
        // Similar test for bam_remapper
    }
}
```

### 6.2 Integration Benchmark

```python
#!/usr/bin/env python3
"""Benchmark threading improvements across all modules."""

import subprocess
import time

BAM = "benchmarking/comprehensive_results/A_sorted_rg.bam"
BED = "benchmarking/comprehensive_results/variants.bed"

# Test BamCounter
print("=== BamCounter ===")
for threads in [1, 4, 8]:
    # Run benchmark
    pass

# Test BamRemapper (via remap_all_chromosomes)
print("=== BamRemapper ===")
# Run make_remap_reads with timing

# Test unified pipeline
print("=== Unified Pipeline ===")
# Run full pipeline with timing
```

---

## 7. Risks and Mitigations

| Risk | Probability | Impact | Mitigation |
|------|-------------|--------|------------|
| Thread explosion (8 Rayon × 2 BGZF = 16 threads) | Medium | Performance regression on <16 core systems | Cap total threads, use shared pool |
| Memory increase | Low | OOM on large BAMs | Monitor, reduce if needed |
| Correctness regression | Very Low | Wrong results | Extensive testing vs GATK baseline |
| No improvement on SSD | Medium | Wasted effort | BGZF threading helps even on fast storage |

---

## 8. Checklist

- [ ] Apply fix to `bam_counter.rs:150`
- [ ] Apply fix to `bam_remapper.rs:353`
- [ ] Apply fix to `multi_sample.rs:587`
- [ ] Run `cargo build --release`
- [ ] Run `cargo test`
- [ ] Benchmark BamCounter (1, 4, 8 threads)
- [ ] Benchmark full pipeline
- [ ] Update documentation
- [ ] Commit with message: "perf: enable BGZF threading in per-chromosome BAM readers"

---

## Appendix: Reference Implementation

See `unified_pipeline.rs:700-707` for correct pattern:

```rust
// CRITICAL: Open a fresh IndexedReader for this thread
let mut bam = bam::IndexedReader::from_path(bam_path)
    .context("Failed to open indexed BAM")?;

// Fetch reads for this chromosome
bam.fetch(chrom).context("Failed to fetch chromosome")?;

// Use a few threads for BAM decompression within this worker
bam.set_threads(2).ok();  // ✓ CORRECT - enables parallel BGZF
```
