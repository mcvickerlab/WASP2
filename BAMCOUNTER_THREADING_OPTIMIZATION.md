# BamCounter Threading Optimization Plan

**Target:** Improve thread scaling efficiency from 70% to 85-90%
**Current:** 5.6x speedup with 8 threads (70% efficiency)
**Goal:** 7.0x+ speedup with 8 threads (87%+ efficiency)

---

## 1. Codebase Architecture

### 1.1 File Structure
```
rust/src/
├── lib.rs              # PyO3 module entry point, exposes BamCounter
├── bam_counter.rs      # ★ TARGET FILE - allele counting implementation
├── bam_filter.rs       # Reference: uses set_threads() correctly
├── bam_intersect.rs    # Reference: uses set_threads() correctly
├── unified_pipeline.rs # Reference: uses set_threads() correctly
└── ... (other modules)
```

### 1.2 Dependencies (Cargo.toml)
```toml
rust-htslib = "0.44"  # BAM I/O with threading support
rayon = "1.8"         # Work-stealing parallelism
rustc-hash = "1.1"    # Fast hashing (FxHashMap/FxHashSet)
```

---

## 2. Current Implementation Analysis

### 2.1 BamCounter Class (`bam_counter.rs`)
```rust
// Lines 48-55: Entry point from Python
#[pyo3(signature = (regions, min_qual=0, threads=1))]
fn count_alleles(&self, py: Python, regions: &PyList, min_qual: u8, threads: usize)

// Lines 94-123: Parallel dispatch by chromosome
if threads > 1 {
    rayon::ThreadPoolBuilder::new().num_threads(threads).build()?.install(|| {
        grouped.par_iter().map(|(chrom, chrom_regions)| {
            self.process_chromosome_reads(chrom, chrom_regions, min_qual, &debug_sites)
        }).collect()
    });
}

// Lines 140-151: Per-chromosome processing (MISSING set_threads!)
fn process_chromosome_reads(&self, chrom: &str, regions: &[(usize, Region)], ...) {
    let mut bam = bam::IndexedReader::from_path(&self.bam_path)?;
    // ⚠️ NO set_threads() call - each reader uses single-threaded decompression
    ...
}
```

### 2.2 Threading Model (Current)
```
User calls count_alleles(threads=8)
         │
         ▼
    ┌────────────────────────────────────────┐
    │       Rayon Thread Pool (8 threads)     │
    │                                         │
    │  Thread 1: chr1 (large)    → BAM reader │◄── No BGZF threading
    │  Thread 2: chr2 (large)    → BAM reader │◄── No BGZF threading
    │  Thread 3: chr3            → BAM reader │◄── No BGZF threading
    │  ...                                    │
    │  Thread 8: chrY (small)    → BAM reader │◄── No BGZF threading
    └────────────────────────────────────────┘
```

**Problems:**
1. **No BGZF threading** - Each BAM reader decompresses single-threaded
2. **Load imbalance** - chr1 has ~3x more variants than chr22
3. **Coarse granularity** - Only 24 work units (chromosomes) for 8 threads

### 2.3 Reference: bam_filter.rs (Correct Pattern)
```rust
// Lines 96-100: Proper use of set_threads()
fn phase2_collect_remap_names(bam_path: &str, ...) {
    let mut bam = bam::Reader::from_path(bam_path)?;
    let num_threads = config.read_threads.min(rayon::current_num_threads());
    bam.set_threads(num_threads).ok();  // ✓ Enables parallel BGZF decompression
    ...
}
```

---

## 3. Optimization Levels

### Level 1: BGZF Decompression Threading (Quick Win)

**Effort:** 5 minutes
**Impact:** +20-30% speedup
**Risk:** Low

#### 3.1.1 Code Change
```rust
// File: rust/src/bam_counter.rs
// Location: process_chromosome_reads(), after line 150

fn process_chromosome_reads(
    &self,
    chrom: &str,
    regions: &[(usize, Region)],
    min_qual: u8,
    debug_sites: &FxHashMap<(String, u32), usize>,
+   decompression_threads: usize,  // NEW PARAMETER
) -> PyResult<FxHashMap<usize, (u32, u32, u32)>> {
    let mut bam = bam::IndexedReader::from_path(&self.bam_path)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(
            format!("Failed to open BAM: {}", e)
        ))?;

+   // Enable parallel BGZF decompression
+   bam.set_threads(decompression_threads).ok();

    // ... rest unchanged
}
```

#### 3.1.2 Thread Allocation Strategy
```rust
// In count_alleles_impl(), calculate optimal thread split:
let total_threads = threads;
let num_chromosomes = grouped.len();

// Split threads between Rayon workers and BGZF decompression
// Rule: Use sqrt(threads) for BGZF if threads > 4
let decompression_threads = if total_threads > 4 {
    (total_threads as f64).sqrt().ceil() as usize
} else {
    1
};
let worker_threads = total_threads;  // Rayon handles work-stealing

// Example: 8 threads → 3 BGZF threads per reader, 8 Rayon workers
```

---

### Level 2: Variant-Based Chunking (Medium Effort)

**Effort:** 1-2 hours
**Impact:** +15-20% additional speedup
**Risk:** Medium (algorithm change)

#### 3.2.1 Problem
Current: Per-chromosome parallelism with ~24 chromosomes
- chr1: ~180K variants (takes longest)
- chr22: ~30K variants (finishes quickly)
- 8 threads waiting for chr1 to finish

#### 3.2.2 Solution
Chunk by variant count, not chromosome:
```rust
// NEW: Replace group_regions_by_chrom with chunk_regions_evenly

fn chunk_regions_evenly(
    regions: &[Region],
    num_chunks: usize,
) -> Vec<Vec<(usize, Region)>> {
    let chunk_size = (regions.len() + num_chunks - 1) / num_chunks;

    regions.iter()
        .enumerate()
        .collect::<Vec<_>>()
        .chunks(chunk_size)
        .map(|chunk| {
            chunk.iter()
                .map(|(idx, r)| (*idx, r.clone()))
                .collect()
        })
        .collect()
}
```

#### 3.2.3 New Processing Model
```rust
fn count_alleles_impl_v2(&self, regions: &[Region], min_qual: u8, threads: usize) {
    // Sort regions by (chrom, pos) for sequential BAM access
    let mut sorted_regions: Vec<(usize, Region)> = regions.iter()
        .enumerate()
        .map(|(i, r)| (i, r.clone()))
        .collect();
    sorted_regions.sort_by(|a, b| {
        (&a.1.chrom, a.1.pos).cmp(&(&b.1.chrom, b.1.pos))
    });

    // Chunk into equal-sized batches
    let chunks = chunk_regions_evenly(&sorted_regions, threads * 4);  // 4x oversubscription

    // Process with work-stealing
    let results: Vec<_> = chunks.par_iter()
        .map(|chunk| self.process_chunk(chunk, min_qual))
        .collect();
}

fn process_chunk(&self, chunk: &[(usize, Region)], min_qual: u8) {
    let mut bam = bam::IndexedReader::from_path(&self.bam_path)?;
    bam.set_threads(2).ok();

    // Group by chromosome within chunk for efficient fetching
    let by_chrom = group_by_chrom(chunk);

    for (chrom, regions) in by_chrom {
        // Fetch span once, then query positions
        self.process_chromosome_regions(&mut bam, &chrom, &regions, min_qual)
    }
}
```

---

### Level 3: Shared Thread Pool (Polish)

**Effort:** 30 minutes
**Impact:** Resource efficiency (prevents thread explosion)
**Risk:** Low

#### 3.3.1 Problem
With Level 1+2, each Rayon worker creates BGZF threads:
- 8 Rayon workers × 2 BGZF threads = 16 htslib threads
- Can cause thread contention on systems with fewer cores

#### 3.3.2 Solution
Use rust-htslib's shared thread pool:
```rust
use rust_htslib::tpool::ThreadPool;

fn count_alleles_impl_v3(&self, regions: &[Region], min_qual: u8, threads: usize) {
    // Create shared thread pool for all BAM readers
    let pool = ThreadPool::new(threads)?;

    // Store pool in Arc for thread-safe sharing
    let pool_arc = Arc::new(pool);

    chunks.par_iter().map(|chunk| {
        let mut bam = bam::IndexedReader::from_path(&self.bam_path)?;
        bam.set_thread_pool(&pool_arc)?;  // Share pool across readers
        self.process_chunk(&mut bam, chunk, min_qual)
    }).collect()
}
```

---

## 4. Implementation Plan

### Phase 1: Level 1 Only (Immediate)
```bash
# 1. Edit bam_counter.rs (5 min)
# 2. Build and test
cd rust && cargo build --release

# 3. Run benchmark
python -c "
from wasp2_rust import BamCounter
counter = BamCounter('test.bam')
# Test with 1, 2, 4, 8 threads
"
```

**Expected Results:**
| Threads | Before (s) | After (s) | Improvement |
|---------|------------|-----------|-------------|
| 1 | 84.65 | ~70 | ~17% |
| 8 | 15.07 | ~12 | ~20% |

### Phase 2: Level 2 (Next Sprint)
1. Implement `chunk_regions_evenly()`
2. Modify `count_alleles_impl()` to use chunking
3. Add `--chunk-size` parameter to Python API
4. Benchmark with various chunk sizes

### Phase 3: Level 3 (Polish)
1. Add ThreadPool support
2. Clean up thread allocation logic
3. Document optimal settings

---

## 5. Code Diff: Level 1 Implementation

```diff
--- a/rust/src/bam_counter.rs
+++ b/rust/src/bam_counter.rs
@@ -1,6 +1,7 @@
 use pyo3::prelude::*;
 use pyo3::types::PyList;
 use rayon::prelude::*;
+use rust_htslib::bam::Read;  // For set_threads trait
 use rust_htslib::{bam, bam::Read as BamRead, bam::ext::BamRecordExtensions};
 use rustc_hash::{FxHashMap, FxHashSet};
 use std::path::Path;
@@ -93,6 +94,10 @@ impl BamCounter {
         let grouped = self.group_regions_by_chrom(regions);
         let debug_sites = parse_debug_sites();

+        // Calculate decompression threads (sqrt heuristic for balance)
+        let decompression_threads = if threads > 4 {
+            ((threads as f64).sqrt().ceil() as usize).max(2)
+        } else { 1 };
+
         // Process chromosomes in parallel if threads > 1
         if threads > 1 {
             rayon::ThreadPoolBuilder::new()
@@ -104,7 +109,7 @@ impl BamCounter {
                     // Process chromosomes in parallel
                     let partial_results: Result<Vec<_>, _> = grouped
                         .par_iter()
-                        .map(|(chrom, chrom_regions)| {
+                        .map(|(chrom, chrom_regions)| {
                             self.process_chromosome_reads(chrom, chrom_regions, min_qual, &debug_sites)
                         })
                         .collect();
@@ -144,6 +149,9 @@ impl BamCounter {
             format!("Failed to open BAM: {}", e)
         ))?;

+    // Enable parallel BGZF decompression for faster BAM reading
+    bam.set_threads(2).ok();
+
     let mut seen_reads: FxHashSet<Vec<u8>> = FxHashSet::default();
     let total_snps: usize = regions.len();
```

---

## 6. Testing Plan

### 6.1 Unit Tests
```rust
#[cfg(test)]
mod threading_tests {
    #[test]
    fn test_scaling_efficiency() {
        let counter = BamCounter::new("test_data/sample.bam".into()).unwrap();
        let regions = load_test_regions();

        let t1 = time(|| counter.count_alleles(&regions, 0, 1));
        let t8 = time(|| counter.count_alleles(&regions, 0, 8));

        let efficiency = (t1 / t8) / 8.0;
        assert!(efficiency > 0.80, "Efficiency {:.1}% below 80% target", efficiency * 100.0);
    }
}
```

### 6.2 Benchmark Script
```python
#!/usr/bin/env python3
"""Benchmark BamCounter threading improvements."""

import time
from wasp2_rust import BamCounter

BAM = "benchmarking/comprehensive_results/A_sorted_rg.bam"
REGIONS = [...]  # Load 2.2M variants

counter = BamCounter(BAM)

for threads in [1, 2, 4, 8]:
    start = time.time()
    results = counter.count_alleles(REGIONS, min_qual=20, threads=threads)
    elapsed = time.time() - start

    speedup = baseline / elapsed if threads > 1 else 1.0
    efficiency = speedup / threads * 100

    print(f"{threads}t: {elapsed:.2f}s, {speedup:.2f}x speedup, {efficiency:.0f}% efficiency")
```

---

## 7. Expected Outcomes

### After Level 1 (BGZF threading)
| Threads | Time (s) | Speedup | Efficiency |
|---------|----------|---------|------------|
| 1 | ~70 | 1.0x | 100% |
| 2 | ~40 | 1.75x | 87% |
| 4 | ~22 | 3.2x | 80% |
| 8 | ~12 | 5.8x | 73% |

### After Level 1+2 (variant chunking)
| Threads | Time (s) | Speedup | Efficiency |
|---------|----------|---------|------------|
| 1 | ~70 | 1.0x | 100% |
| 2 | ~37 | 1.9x | 95% |
| 4 | ~19 | 3.7x | 92% |
| 8 | ~10 | 7.0x | 87% |

---

## 8. Risks and Mitigations

| Risk | Impact | Mitigation |
|------|--------|------------|
| Thread contention on NFS | Performance regression | Use shared thread pool (Level 3) |
| Memory increase per thread | OOM on large BAMs | Cap BGZF threads at 2-3 per reader |
| API breaking change | User code fails | Keep `threads` param semantics |
| Correctness regression | Wrong counts | Extensive testing against GATK |

---

## 9. Future Optimizations (Out of Scope)

1. **Single-pass streaming** - Read BAM once, check all variants per read
2. **Memory-mapped I/O** - Avoid kernel-user copies for BAM data
3. **SIMD base comparison** - Vectorized allele matching
4. **GPU acceleration** - Offload counting to CUDA/Metal

---

## Appendix: Useful References

- [rust-htslib set_threads](https://docs.rs/rust-htslib/latest/rust_htslib/bam/trait.Read.html#method.set_threads)
- [htslib BGZF threading](https://github.com/samtools/htslib/blob/develop/htslib/bgzf.h)
- [Rayon work-stealing](https://docs.rs/rayon/latest/rayon/)
- [Sambamba threading model](https://pmc.ncbi.nlm.nih.gov/articles/PMC4765878/)
