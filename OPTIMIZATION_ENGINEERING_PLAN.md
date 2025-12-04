# WASP2 Rust Optimization Engineering Plan

**Created:** 2025-12-03
**Branch:** ropc-indels
**Current State:** Unified pipeline at 99.8% baseline match (3,076 baseline-only, 344 unified-only)

---

## Executive Summary

This plan outlines optimizations for the WASP2 Rust acceleration module, targeting additional 2-5x speedup on top of the unified pipeline. Each optimization includes research findings, implementation approach, validation strategy, and rollback plan.

---

## Phase 1: Quick Wins (Low Risk, High Impact)

### 1.1 Upgrade rust-htslib to v0.51.0

**Research Findings:**
- Current: v0.44.0, Latest: v0.51.0 (Oct 2025)
- Key changes since v0.44:
  - v0.49.0: Fixed memory leak in faidx sequence fetching
  - v0.50.0: Added `bam::Record::set_cigar`, non-diploid genotype support
  - v0.51.0: Support for updating aux tags in-place
  - Various clippy/formatting fixes

**Implementation:**
```toml
# Cargo.toml
rust-htslib = { version = "0.51", default-features = false }
```

**Validation:**
```bash
cd rust && cargo build --release
cd .. && pytest tests/test_rust_python_match.py -v
```

**Risk:** Low - backwards compatible changes
**Expected Impact:** Bug fixes, minor performance improvements

---

### 1.2 Pre-allocated Record Reading

**Research Findings:**
- Current: Using `.records()` iterator (allocates new Record per read)
- Docs: "Using the iterator is less efficient than pre-allocating a Record"
- Expected: ~10-20% speedup on BAM streaming

**Current Code (unified_pipeline.rs:466):**
```rust
for result in bam.records() {
    let read = result?;
    // ...
}
```

**Optimized Code:**
```rust
let mut record = bam::Record::new();
while let Some(result) = bam.read(&mut record) {
    result?;
    // Process record (need to clone if buffering)
    // ...
}
```

**Caveat:** Our pair buffer stores records, so we need to clone when buffering:
```rust
// When buffering first mate:
pair_buffer.insert(read_name.clone(), record.clone());

// Reset for next iteration
record = bam::Record::new();  // Or reuse after pair completion
```

**Alternative - rc_records():**
```rust
for result in bam.rc_records() {
    let record = result?;
    // record is Rc<Record>, cheaper than full allocation
}
```

**Validation:**
```bash
# Run unified pipeline on chr21 test data
python -c "
import wasp2_rust
stats = wasp2_rust.unified_make_reads_py(
    'benchmarking/star_wasp_comparison/results/wasp2_run/A_sorted_chr21.bam',
    '/tmp/test_variants.bed',
    '/tmp/r1.fq.gz', '/tmp/r2.fq.gz'
)
print(stats)
"
```

**Risk:** Medium - requires careful handling of record ownership
**Expected Impact:** 10-20% speedup on BAM streaming phase

---

### 1.3 Shared ThreadPool for BAM I/O

**Research Findings:**
- rust-htslib supports `ThreadPool::new(n_threads)` for shared pools
- Allows controlling total thread count across multiple readers/writers
- Currently each call to `set_threads()` creates independent thread pool

**Implementation:**
```rust
use rust_htslib::tpool::ThreadPool;

// Create shared pool once
let tpool = ThreadPool::new(config.read_threads)?;

// Use for BAM reader
bam.set_thread_pool(&tpool)?;

// If we add BAM writer later, share the pool
writer.set_thread_pool(&tpool)?;
```

**Validation:** Same as 1.2
**Risk:** Low
**Expected Impact:** Better thread utilization, 5-10% in multi-file scenarios

---

## Phase 2: Parallel Processing (Medium Risk, High Impact)

### 2.1 Parallel FASTQ Compression with gzp

**Research Findings:**
- [gzp crate](https://github.com/sstadick/gzp): Multi-threaded compression, pigz-like
- Near drop-in replacement for `Write` trait
- Compression time scales proportionally with threads
- Supports BGZF format (bioinformatics standard)

**Current Code (unified_pipeline.rs:365-376):**
```rust
let mut r1_writer = BufWriter::with_capacity(
    1024 * 1024,
    GzEncoder::new(r1_file, Compression::fast()),
);
```

**Optimized Code:**
```rust
use gzp::{deflate::Gzip, ZBuilder};

let mut r1_writer = ZBuilder::<Gzip, _>::new()
    .num_threads(4)  // Parallel compression threads
    .compression_level(flate2::Compression::fast())
    .from_writer(BufWriter::with_capacity(1024 * 1024, r1_file));
```

**Cargo.toml Addition:**
```toml
gzp = { version = "0.11", default-features = false, features = ["deflate_default"] }
```

**Validation:**
```bash
# Compare output file sizes and content
zcat /tmp/r1_baseline.fq.gz | md5sum
zcat /tmp/r1_optimized.fq.gz | md5sum
```

**Risk:** Low - output format identical (concatenated gzip blocks)
**Expected Impact:** 2-4x speedup on FASTQ writing (currently ~30% of runtime)

---

### 2.2 Parallel Chromosome Processing with IndexedReader

**Research Findings:**
- rust-htslib `IndexedReader` supports region fetching via `.fetch()`
- Each chromosome can be processed independently with separate reader
- Current unified pipeline is single-threaded sequential BAM read

**Architecture:**
```
┌─────────────────────────────────────────────────────┐
│                    Main Thread                       │
│  1. Build VariantStore (coitrees)                   │
│  2. Get chromosome list from BAM header             │
│  3. Spawn per-chromosome workers                    │
└─────────────────────────────────────────────────────┘
                         │
         ┌───────────────┼───────────────┐
         ▼               ▼               ▼
┌─────────────┐  ┌─────────────┐  ┌─────────────┐
│  Thread 1   │  │  Thread 2   │  │  Thread N   │
│   chr1      │  │   chr2      │  │   chrN      │
│ IndexReader │  │ IndexReader │  │ IndexReader │
│ → channel   │  │ → channel   │  │ → channel   │
└─────────────┘  └─────────────┘  └─────────────┘
         │               │               │
         └───────────────┼───────────────┘
                         ▼
              ┌─────────────────────┐
              │   Writer Thread     │
              │  FASTQ output       │
              └─────────────────────┘
```

**Implementation Sketch:**
```rust
pub fn unified_make_reads_parallel(
    bam_path: &str,
    bed_path: &str,
    r1_path: &str,
    r2_path: &str,
    config: &UnifiedConfig,
) -> Result<UnifiedStats> {
    // 1. Build variant store (shared, read-only)
    let store = Arc::new(build_variant_store(bed_path)?);

    // 2. Get chromosomes from BAM header
    let bam = bam::IndexedReader::from_path(bam_path)?;
    let chroms: Vec<String> = (0..bam.header().target_count())
        .map(|tid| String::from_utf8_lossy(bam.header().tid2name(tid)).to_string())
        .collect();

    // 3. Setup output channel
    let (tx, rx) = bounded(config.channel_buffer);

    // 4. Spawn writer thread
    let writer_handle = spawn_writer_thread(rx, r1_path, r2_path);

    // 5. Process chromosomes in parallel
    let stats: Vec<UnifiedStats> = chroms.par_iter()
        .map(|chrom| {
            // Each thread opens its own IndexedReader
            let mut bam = bam::IndexedReader::from_path(bam_path)?;
            bam.fetch(chrom)?;

            process_chromosome(&bam, &store, chrom, &tx, config)
        })
        .collect::<Result<Vec<_>>>()?;

    // 6. Aggregate stats
    Ok(stats.into_iter().fold(UnifiedStats::default(), |a, b| a.merge(b)))
}
```

**Validation:**
```bash
# Compare read counts between sequential and parallel
python -c "
import wasp2_rust
# Run both versions, compare haplotypes_written
"
```

**Risk:** Medium - requires BAM index, different traversal order
**Expected Impact:** 3-8x speedup depending on chromosome count and thread availability

---

## Phase 3: Advanced Optimizations (Higher Risk)

### 3.1 bcf_reader for Parallel VCF Decompression

**Research Findings:**
- [bcf_reader crate](https://docs.rs/bcf_reader): Parallel BGZF decompression with rayon
- `IndexedBcfReader` for random access with CSI index
- Early stage development, but addresses decompression bottleneck

**Use Case:** Only beneficial if VCF→BED conversion is a bottleneck (likely not - we cache BED files)

**Risk:** High - early stage crate
**Expected Impact:** Minimal for current workflow (BED is pre-generated)

---

### 3.2 noodles Async I/O

**Research Findings:**
- [noodles](https://github.com/zaeleus/noodles): Pure Rust bioinformatics I/O
- Async support with Tokio for BAM, BCF, VCF
- Benchmarks show 1.5x speedup for async BAM reading
- `lazy_records()` for reduced allocation

**Use Case:** Could replace rust-htslib for BAM reading

**Risk:** High - major rewrite, different API
**Expected Impact:** 1.5x speedup on BAM reading, but high implementation cost

---

## Test Data and Validation Framework

### Available Test Data

| File | Size | Purpose |
|------|------|---------|
| `benchmarking/.../A_sorted_chr21.bam` | 32MB | Quick validation |
| `benchmarking/.../A_sorted.bam` | Full | Performance benchmark |
| `tests/data/sample.vcf` | <1KB | Unit tests |

### Validation Script Template

```python
#!/usr/bin/env python3
"""Validate optimization changes against baseline."""

import subprocess
import tempfile
import time
from pathlib import Path

def run_baseline(bam_path, bed_path, out_dir):
    """Run multi-pass baseline pipeline."""
    # filter -> intersect -> remap
    pass

def run_optimized(bam_path, bed_path, out_dir):
    """Run optimized unified pipeline."""
    import wasp2_rust
    return wasp2_rust.unified_make_reads_py(
        bam_path, bed_path,
        str(out_dir / "r1.fq.gz"),
        str(out_dir / "r2.fq.gz"),
    )

def compare_outputs(baseline_dir, optimized_dir):
    """Compare FASTQ outputs for correctness."""
    # Sort and compare read names/sequences
    pass

def benchmark(func, *args, n_runs=3):
    """Benchmark a function."""
    times = []
    for _ in range(n_runs):
        start = time.time()
        result = func(*args)
        times.append(time.time() - start)
    return min(times), result

if __name__ == "__main__":
    BAM = "benchmarking/star_wasp_comparison/results/wasp2_run/A_sorted_chr21.bam"
    BED = "/tmp/test_variants.bed"

    with tempfile.TemporaryDirectory() as tmpdir:
        baseline_time, _ = benchmark(run_baseline, BAM, BED, Path(tmpdir) / "baseline")
        optimized_time, stats = benchmark(run_optimized, BAM, BED, Path(tmpdir) / "optimized")

        print(f"Baseline: {baseline_time:.2f}s")
        print(f"Optimized: {optimized_time:.2f}s")
        print(f"Speedup: {baseline_time / optimized_time:.2f}x")
        print(f"Stats: {stats}")
```

---

## Implementation Order

| Priority | Optimization | Effort | Risk | Impact |
|----------|-------------|--------|------|--------|
| 1 | rust-htslib upgrade | 30min | Low | Bug fixes |
| 2 | Pre-allocated Record | 2hr | Medium | 10-20% |
| 3 | gzp parallel compression | 2hr | Low | 2-4x on writes |
| 4 | Shared ThreadPool | 1hr | Low | 5-10% |
| 5 | Parallel chromosome | 4hr | Medium | 3-8x |

---

## Rollback Strategy

Each optimization should be:
1. Implemented in a feature branch
2. Gated behind a feature flag or config option
3. Validated against baseline before merge
4. Tagged for easy rollback

```rust
// Example feature flag
pub struct UnifiedConfig {
    pub use_parallel_chromosomes: bool,  // Default: false
    pub use_gzp_compression: bool,       // Default: false
    // ...
}
```

---

## Success Metrics

| Metric | Current | Target |
|--------|---------|--------|
| 150M reads processing time | ~500s baseline | <100s |
| Baseline match rate | 99.8% | 99.9%+ |
| Memory usage | ~1.3GB | <2GB |
| Thread utilization | 1 core | All available |

---

## References

- [rust-htslib GitHub](https://github.com/rust-bio/rust-htslib)
- [rust-htslib CHANGELOG](https://github.com/rust-bio/rust-htslib/blob/master/CHANGELOG.md)
- [gzp - Multi-threaded Compression](https://github.com/sstadick/gzp)
- [coitrees - Interval Trees](https://github.com/dcjones/coitrees)
- [noodles - Bioinformatics I/O](https://github.com/zaeleus/noodles)
- [bcf_reader - Parallel BCF](https://docs.rs/bcf_reader)
