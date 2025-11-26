#!/usr/bin/env python3
"""
Profile the Rust indel processing code to identify bottlenecks.

This script analyzes the mapping_filter.rs code to understand where time
is spent during WASP filtering of remapped reads.
"""

import time
import sys
from pathlib import Path
from collections import defaultdict

# Analysis of mapping_filter.rs code structure
PROFILING_REPORT = """
================================================================================
RUST INDEL PROCESSING PROFILING REPORT
================================================================================

CODE ANALYSIS: rust/src/mapping_filter.rs

The filter_bam_wasp() function performs the following operations:

1. FIRST PASS - Read remapped BAM and build keep set:
   - Parse BAM records (I/O bound - via rust-htslib)
   - Parse WASP-encoded read names (CPU: string operations)
   - Extract position/count info from suffix
   - Build HashMap of expected positions
   - Track remaining copies to see
   - Position matching with optional slop tolerance

2. SECOND PASS - Filter original BAM:
   - Read original to_remap BAM (I/O bound)
   - Check if read name is in keep_set (CPU: HashMap lookup)
   - Write kept reads to output BAM (I/O bound)

================================================================================
BOTTLENECK ANALYSIS (based on code review)
================================================================================

TOP 3 POTENTIAL BOTTLENECKS:

1. STRING ALLOCATION IN HOT LOOP (Lines 72-80)
   -----------------------------------------------
   Code:
        let name = match std::str::from_utf8(orig_name) {
            Ok(s) => s.to_owned(),  // <-- ALLOCATION on EVERY read
            Err(_) => continue,
        };

   Impact: MEDIUM-HIGH
   - Allocates owned String for every read with variants
   - Could use Cow<str> or cache byte slices instead
   - Estimated overhead: ~10-20% for high-variant datasets

2. QNAME PARSING (Lines 41-75)
   ------------------------------
   Code:
        let split_idx = qname.windows(6).position(|w| w == b"_WASP_");
        let parts: Vec<&[u8]> = suffix.split(|b| *b == b'_').collect();
        parse_i64(parts[0])  // Parse positions

   Impact: MEDIUM
   - Linear scan for "_WASP_" marker on every read
   - Multiple allocations for split operation
   - UTF-8 validation and parsing for positions
   - Estimated overhead: ~15-25% of total time

3. HASHMAP OPERATIONS (Lines 77-88, 171)
   ----------------------------------------
   Code:
        if !pos_map.contains_key(&name) {
            pos_map.insert(name.clone(), (pos1, pos2));
            remaining.insert(name.clone(), total);
            keep_set.insert(name.clone());  // <-- Multiple clones
        }

   Impact: LOW-MEDIUM
   - Using FxHashSet (good choice!)
   - Multiple HashMap lookups and inserts
   - String cloning for keys (see #1)
   - Estimated overhead: ~5-10% for large datasets

================================================================================
SECONDARY BOTTLENECKS:
================================================================================

4. POSITION MATCHING LOGIC (Lines 95-108)
   - Branching logic for slop tolerance
   - Multiple abs() operations for indel matching
   - Impact: LOW (~5%)

5. BAM I/O (Lines 26-36, 149-177)
   - rust-htslib is already highly optimized
   - Cannot optimize further without changing library
   - Impact: Not a bottleneck (efficient C library)

================================================================================
MEMORY USAGE ANALYSIS
================================================================================

HashMaps grow linearly with number of reads:
- keep_set: String → O(n) where n = reads with variants
- pos_map: String → (i64, i64) → ~40 bytes per entry
- remaining: String → i64 → ~32 bytes per entry

For 1M reads with variants:
- keep_set: ~50 MB (assuming avg 32 bytes per name)
- pos_map: ~70 MB
- remaining: ~60 MB
Total: ~180 MB (reasonable, not a concern)

================================================================================
OPTIMIZATION RECOMMENDATIONS
================================================================================

PHASE 1: LOW-HANGING FRUIT (If profiling confirms bottleneck)

1. REDUCE STRING ALLOCATIONS:

   Current:
   ```rust
   let name = match std::str::from_utf8(orig_name) {
       Ok(s) => s.to_owned(),  // Allocates!
       Err(_) => continue,
   };
   ```

   Optimized:
   ```rust
   // Use byte slices as keys, only convert when writing output
   use rustc_hash::FxHashMap;
   let mut keep_set: FxHashSet<Vec<u8>> = FxHashSet::default();

   // Or use Cow for zero-copy when possible
   use std::borrow::Cow;
   let name: Cow<str> = match std::str::from_utf8(orig_name) {
       Ok(s) => Cow::Borrowed(s),
       Err(_) => continue,
   };
   ```

   Expected speedup: 10-20%

2. OPTIMIZE QNAME PARSING:

   Current: Linear scan + split + collect

   Optimized:
   ```rust
   // Use memchr for faster searching
   use memchr::memmem;
   let finder = memmem::Finder::new(b"_WASP_");
   let split_idx = finder.find(qname);

   // Avoid Vec allocation - parse directly
   let mut parts = suffix.splitn(5, |b| *b == b'_');
   let pos1 = parse_i64(parts.next()?)?;
   let pos2 = parse_i64(parts.next()?)?;
   // ... etc
   ```

   Expected speedup: 5-15%

3. PRE-ALLOCATE HASHMAPS:

   Current:
   ```rust
   let mut keep_set: FxHashSet<String> = FxHashSet::default();
   ```

   Optimized:
   ```rust
   // If you know approximate read count from BAM header
   let mut keep_set: FxHashSet<String> = FxHashSet::with_capacity_and_hasher(
       estimated_count, Default::default()
   );
   ```

   Expected speedup: 2-5%

PHASE 2: AGGRESSIVE OPTIMIZATION (If needed)

4. UNSAFE OPTIMIZATIONS (Use with caution!):

   - Use unsafe byte-to-string conversions if input is trusted
   - Skip UTF-8 validation for read names (usually ASCII)
   - Direct memory access for position parsing

   Expected speedup: 5-10% (but adds risk)

5. PARALLELIZATION:

   - Process chromosomes in parallel using rayon
   - Requires refactoring to work chromosome-by-chromosome

   Expected speedup: 2-4x on multi-core systems

================================================================================
PROFILING METHODOLOGY
================================================================================

RECOMMENDED APPROACH:

1. Create synthetic test data:
   - 100K reads with WASP suffixes
   - Mix of matching/non-matching positions
   - Varying numbers of variants per read

2. Use `cargo bench` with criterion.rs:
   - Micro-benchmarks for individual functions
   - Full integration benchmark
   - Compare before/after optimizations

3. Use `cargo flamegraph` for visual profiling:
   ```bash
   cargo install flamegraph
   cargo flamegraph --bench mapping_filter_bench
   ```

   This will show exactly where CPU time is spent.

4. Memory profiling with valgrind/massif:
   ```bash
   valgrind --tool=massif target/release/your_test
   ms_print massif.out.XXX
   ```

================================================================================
EXPECTED RESULTS
================================================================================

Based on code analysis, the CURRENT implementation is already quite good:

✓ Uses FxHashSet (fast hash function)
✓ Uses rust-htslib (optimized BAM I/O)
✓ Minimal allocations in critical path
✓ Simple, readable code

PREDICTION:

Without profiling on real data, I estimate:
- 60% of time: BAM I/O (can't optimize)
- 20% of time: String operations (optimizable)
- 10% of time: HashMap operations (slightly optimizable)
- 10% of time: Position matching (not worth optimizing)

LIKELY OUTCOME:
→ Code is probably ALREADY FAST ENOUGH
→ Optimizations would yield 10-30% speedup at most
→ Only optimize if profiling shows it's actually slow

================================================================================
NEXT STEPS
================================================================================

1. [ ] Create realistic test dataset (100K+ reads)
2. [ ] Measure baseline performance (reads/sec)
3. [ ] Run cargo bench for micro-benchmarks
4. [ ] Generate flamegraph to visualize bottlenecks
5. [ ] IF bottlenecks found: Implement optimizations
6. [ ] Measure speedup vs baseline
7. [ ] Validate correctness (results must match!)

DECISION CRITERIA:

- If throughput > 100K reads/sec → Don't optimize (fast enough)
- If throughput < 50K reads/sec → Optimize string allocations
- If throughput < 10K reads/sec → Major refactor needed

================================================================================
CONCLUSION
================================================================================

The Rust code in mapping_filter.rs is already well-written with good
algorithmic choices. The main opportunities for optimization are:

1. Reducing string allocations (use byte slices or Cow)
2. Faster WASP suffix parsing (memchr + avoid Vec allocation)
3. HashMap pre-allocation if read count is known

However, these optimizations are premature without profiling data.

RECOMMENDATION: Profile first, optimize only if needed.
"""


def main():
    print(PROFILING_REPORT)

    # Save report to file
    report_path = Path(__file__).parent / "RUST_PROFILING_REPORT.md"
    with open(report_path, "w") as f:
        f.write(PROFILING_REPORT)

    print(f"\nReport saved to: {report_path}")

    # Code metrics
    print("\n" + "=" * 80)
    print("CODE METRICS")
    print("=" * 80)

    mapping_filter_path = Path(__file__).parent / "rust" / "src" / "mapping_filter.rs"

    if mapping_filter_path.exists():
        with open(mapping_filter_path) as f:
            lines = f.readlines()

        # Count operations in main loop
        in_loop = False
        string_ops = 0
        hashmap_ops = 0
        parse_ops = 0

        for i, line in enumerate(lines, 1):
            if "for rec_res in remapped_reader.records()" in line:
                in_loop = True
            if in_loop and "for rec_res in to_reader.records()" in line:
                break

            if in_loop:
                if "to_owned()" in line or ".clone()" in line:
                    string_ops += 1
                    print(f"  Line {i}: String allocation - {line.strip()}")
                if ".insert(" in line or ".get(" in line or ".contains_key(" in line:
                    hashmap_ops += 1
                if ".parse::<i64>()" in line:
                    parse_ops += 1

        print(f"\nIn main processing loop:")
        print(f"  String allocations: {string_ops}")
        print(f"  HashMap operations: {hashmap_ops}")
        print(f"  Parse operations: {parse_ops}")
    else:
        print(f"Warning: Could not find {mapping_filter_path}")


if __name__ == "__main__":
    main()
