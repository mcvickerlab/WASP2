# Rust Mapping Optimization - Skeleton Summary

**Created:** 2025-11-20
**Status:** Skeleton complete, ready for implementation
**Expected Impact:** 7-20x speedup for mapping stage

---

## What Was Created

### ✅ Completed Skeleton Structure

```
rust/
├── src/
│   ├── lib.rs              # PyO3 bindings (UPDATED with new functions)
│   ├── bam_counter.rs      # Allele counting (EXISTING - already works)
│   ├── bam_remapper.rs     # ✨ NEW: Allele swapping skeleton
│   └── read_pairer.rs      # ✨ NEW: Read pairing utilities
└── README.md               # Updated with remapping docs
```

### New Files

#### 1. `rust/src/bam_remapper.rs` (367 lines)

**Data Structures:**
- `Variant` - Single SNP with haplotype alleles
- `RemapConfig` - Configuration (max_seqs, is_phased)
- `HaplotypeRead` - Generated read with swapped alleles
- `RemapStats` - Statistics tracking

**Core Functions (with TODOs):**
- `parse_intersect_bed()` - Parse BED → HashMap (replaces Polars)
- `swap_alleles_for_chrom()` - Process one chromosome
- `generate_haplotype_seqs()` - Core allele swapping
- `write_fastq_pair()` - Write FASTQ output
- `process_all_chromosomes_parallel()` - Parallel processing

**Expected Speedups:**
- `parse_intersect_bed`: 0.316s → 0.040s (8x)
- `swap_alleles_for_chrom`: 0.147s → 0.020s (7x)
- Overall: 0.500s → 0.070s sequential, 0.025s parallel (7-20x)

#### 2. `rust/src/read_pairer.rs` (218 lines)

**Data Structures:**
- `PairingStats` - Track discarded reads
- `ReadPairer` - Iterator for paired reads

**Core Functionality (with TODOs):**
- `ReadPairer::new()` - Create pairer for BAM file
- `ReadPairer::for_chromosome()` - Chromosome-specific pairing
- `passes_filters()` - Filter unmapped/secondary/supplementary
- `process_read()` - Pair reads by name

**Expected Speedup:** 2-3x over Python dict-based pairing

#### 3. `rust/src/lib.rs` (Updated)

**New PyO3 Functions:**
```python
# Single chromosome
wasp2_rust.remap_chromosome(
    bam_path, intersect_bed, chrom, out_r1, out_r2, max_seqs=64
) -> (pairs_processed, haplotypes_generated)

# All chromosomes (parallel)
wasp2_rust.remap_all_chromosomes(
    bam_path, intersect_bed, out_r1, out_r2, max_seqs=64
) -> (pairs_processed, haplotypes_generated)
```

Currently returns error: "not yet implemented - skeleton only"

---

## Profiling Data

### Baseline (Python)

From `baselines/mapping/allele_swap_profile.txt`:

| Component | Time | % |
|-----------|------|---|
| Polars DataFrame ops | 0.363s | 72% |
| make_intersect_df | 0.316s | 63% |
| swap_chrom_alleles | 0.147s | 29% |
| CSV parsing | 0.059s | 12% |
| partition_by | 0.040s | 8% |
| **TOTAL** | **0.500s** | **100%** |

### Target (Rust)

| Component | Rust Time | Speedup |
|-----------|-----------|---------|
| parse_intersect_bed (streaming) | 0.040s | 8x |
| swap_alleles (byte ops) | 0.020s | 7x |
| Parallel (8 cores) | 0.025s | 20x |

---

## Implementation Roadmap

### Phase 1: Core Infrastructure (Week 1-2)
- [x] Create module skeletons
- [x] Define data structures
- [x] Add PyO3 bindings
- [ ] Implement `ReadPairer`
- [ ] Add unit tests

### Phase 2: Intersection Parser (Week 2-3)
- [ ] Implement `parse_intersect_bed()`
- [ ] Benchmark vs Polars
- [ ] Integration test with real data

### Phase 3: Allele Swapping (Week 3-5)
- [ ] Implement `build_alignment_map()`
- [ ] Implement `generate_haplotype_seqs()`
- [ ] Implement `swap_alleles_for_chrom()`
- [ ] Validate output matches Python exactly

### Phase 4: I/O (Week 5-6)
- [ ] Implement `write_fastq_pair()`
- [ ] Handle paired-end FASTQ format
- [ ] Optimize buffer sizes

### Phase 5: Parallelization (Week 6-7)
- [ ] Implement `process_all_chromosomes_parallel()`
- [ ] Add rayon thread pool configuration
- [ ] Benchmark scaling on multiple cores

### Phase 6: Integration (Week 7-8)
- [ ] Update `src/mapping/make_remap_reads.py` to use Rust
- [ ] Add feature flag: `USE_RUST_REMAPPER`
- [ ] Full pipeline regression testing
- [ ] Performance benchmarking vs Python
- [ ] Documentation and examples

---

## Key Design Decisions

### 1. FxHashMap Instead of Python Dict
- **Why:** 2-3x faster lookups for byte keys
- **Where:** Read name → variants mapping

### 2. Byte Slices (`&[u8]`) Instead of Strings
- **Why:** Zero UTF-8 validation overhead
- **Where:** Read names, sequences, quality scores

### 3. In-Place Modification
- **Why:** Avoid Python's immutable string allocations
- **Where:** Allele swapping in sequences

### 4. Streaming BED Parser
- **Why:** Avoid full Polars DataFrame materialization
- **Where:** `parse_intersect_bed()`

### 5. Rayon for Parallelism
- **Why:** Easy chromosome-level parallelism
- **Where:** `process_all_chromosomes_parallel()`

---

## Testing Strategy

### Unit Tests (in each module)
```rust
#[cfg(test)]
mod tests {
    #[test]
    #[ignore]  // Remove when implementing
    fn test_function_name() {
        // TODO: Test implementation
    }
}
```

### Integration Tests
1. Create small test BAM + intersect.bed
2. Run Rust implementation
3. Compare output with Python (exact match)
4. Verify statistics match

### Regression Tests
1. Run full pipeline with Rust
2. Compare final counts with Python baseline
3. Verify no SNPs differ

---

## How to Continue Implementation

### Option A: Start with Intersection Parser (Easiest Win)

```rust
// In rust/src/bam_remapper.rs, replace:
pub fn parse_intersect_bed<P: AsRef<Path>>(
    intersect_bed: P,
) -> Result<FxHashMap<Vec<u8>, Vec<Variant>>> {
    // TODO: Implement streaming BED parser
    unimplemented!()
}

// With actual implementation:
pub fn parse_intersect_bed<P: AsRef<Path>>(
    intersect_bed: P,
) -> Result<FxHashMap<Vec<u8>, Vec<Variant>>> {
    let file = File::open(intersect_bed.as_ref())?;
    let reader = BufReader::new(file);
    let mut variants: FxHashMap<Vec<u8>, Vec<Variant>> = FxHashMap::default();

    for line in reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();

        // Parse BED format:
        // 0: chrom, 1: start, 2: end, 3: read/mate, ...
        // 10: vcf_chrom, 11: vcf_start, 12: vcf_end,
        // 13: ref, 14: alt, 15: GT (hap1|hap2)

        let read_name = fields[3].split('/').next().unwrap().as_bytes().to_vec();
        let chrom = fields[10].to_string();
        let pos = fields[11].parse::<u32>()?;
        let ref_allele = fields[13].as_bytes().to_vec();
        let alt_allele = fields[14].as_bytes().to_vec();

        // Parse phased genotype
        let gt_parts: Vec<&str> = fields[15].split('|').collect();
        let hap1 = gt_parts[0].as_bytes().to_vec();
        let hap2 = gt_parts[1].as_bytes().to_vec();

        let variant = Variant {
            chrom,
            pos,
            ref_allele,
            alt_allele,
            hap1,
            hap2,
        };

        variants.entry(read_name)
               .or_insert_with(Vec::new)
               .push(variant);
    }

    Ok(variants)
}
```

### Option B: Start with Read Pairer (More Complex)

See `rust/src/read_pairer.rs` TODOs - needs BAM iteration logic.

### Option C: Start with Tests (TDD Approach)

Write tests first, then implement to pass tests.

---

## Validation

To ensure correctness:

1. **Exact output match:** Rust FASTQ must match Python FASTQ byte-for-byte
2. **Statistics match:** Pair counts, discard counts, haplotype counts
3. **Downstream compatibility:** Remapped BAM produces identical counting results

---

## Documentation

- ✅ **Profiling:** `MAPPING_OPTIMIZATION_ANALYSIS.md`
- ✅ **Skeleton:** This document
- ✅ **API:** Rust module doc comments
- ✅ **Usage:** `rust/README.md`
- 📝 **Next:** Implementation notes as you code

---

## Questions / Decisions Needed

1. **BED format:** Verify the exact column indices for intersection file
2. **FASTQ names:** Confirm WASP naming format: `{name}_WASP_{pos1}_{pos2}_{num}_{total}`
3. **Error handling:** Should failures be fatal or continue with warnings?
4. **Memory limits:** Should we cap the HashMap size?
5. **Python integration:** Feature flag or auto-detect?

---

## Success Metrics

✅ **Skeleton Complete** - All data structures and function signatures defined
⏳ **Phase 1:** Rust code compiles with basic implementations
⏳ **Phase 2:** Unit tests pass
⏳ **Phase 3:** Integration tests match Python output exactly
⏳ **Phase 4:** Performance tests show 7-20x speedup
⏳ **Phase 5:** Full pipeline regression tests pass
⏳ **Phase 6:** Production ready

---

## Contact / Next Steps

**Ready to implement?** Pick a starting point:

1. **Easiest win:** Intersection parser (8x speedup alone)
2. **Most impactful:** Allele swapping (core bottleneck)
3. **Most fun:** Parallel processing (rayon magic)

**Need help?** See:
- `MAPPING_OPTIMIZATION_ANALYSIS.md` for detailed profiling
- `rust/src/bam_remapper.rs` for TODOs
- `rust/README.md` for build instructions
