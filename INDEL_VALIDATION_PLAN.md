# Indel Validation Plan - Sanity Check

**Date:** 2025-11-25
**Purpose:** Define expected indel behavior before Rust implementation
**Branch:** rust-optimization-plink2

---

## Executive Summary

This document defines **ground truth** for how WASP2 should handle insertions and deletions during read remapping. It serves as:

1. **Sanity check** for current Python implementation
2. **Acceptance criteria** for Rust implementation
3. **Test specification** for validation scripts

**Key Principle:** WASP2 swaps alleles to create alternate haplotypes for bias testing. Indels change sequence length, requiring special handling for position mapping and quality scores.

---

## 1. What Are Indels in WASP2?

### Insertion (Read has MORE bases than reference)
```
Reference:  ACGT----TGCA
Read:       ACGTAAAATGCA
            ^^^^    ^^^^
            Match  4bp INS
```

**VCF Representation:**
```
#CHROM  POS  REF  ALT
chr1    100  T    TAAAA   # 4bp insertion
```

### Deletion (Read has FEWER bases than reference)
```
Reference:  ACGTAAAATGCA
Read:       ACGT----TGCA
            ^^^^    ^^^^
            Match  4bp DEL
```

**VCF Representation:**
```
#CHROM  POS  REF     ALT
chr1    100  TAAAA   T      # 4bp deletion
```

---

## 2. Expected Behavior: Position Mapping

### Problem: Reference Position → Read Position

For **SNPs**, mapping is 1:1:
```
Ref pos:  100 101 102 103
Read pos:  50  51  52  53
```

For **Insertions**, read has extra bases:
```
Ref:  T G C
      | | |
Read: T G A A A C
Pos:  50 51 52 53 54 55

Query: "Where is ref pos 102 (C) in the read?"
Answer: Position 55 (must skip over insertion at 51-54)
```

For **Deletions**, read is missing bases:
```
Ref:  T G A A A C
      | |       |
Read: T G - - - C
Pos:  50 51    52

Query: "Where is ref pos 102 (A, deleted) in the read?"
Answer: Not in read - use flanking positions (51 or 52)
```

### Solution: Left/Right Flanking Maps

Python implementation uses **two dictionaries**:

```python
ref2q_left[ref_pos] = query_pos   # Nearest left position
ref2q_right[ref_pos] = query_pos  # Nearest right position
```

**For deletions:**
- `ref2q_left[102]` = 51 (position before deletion)
- `ref2q_right[102]` = 52 (position after deletion)

**For insertions:**
- Both maps point to the same position

---

## 3. Expected Behavior: Sequence Transformation

### Test Case 1: Simple SNP (baseline)

**Original Read:**
```
Sequence: ACGTACGTACGT
Position: 000111222333
          012345678901
```

**Variant:** chr1:5 G→T (SNP)

**Expected Output:**
```
Haplotype 1: ACGT[T]CGTACGT  (swap G→T at pos 5)
Haplotype 2: ACGT[G]CGTACGT  (keep reference G)
```

**Quality Scores:** Unchanged (same length)

---

### Test Case 2: 1bp Insertion

**Original Read:**
```
Sequence: ACGTACGTACGT
Ref Pos:  100-103 (matches ACGT)
```

**Variant:** chr1:102 T→TA (1bp insertion of A)

**Expected Transformation:**
```
Read segment split:
  [0]: "ACGT"      (before variant, ref pos 100-101)
  [1]: "A"         (variant site, ref pos 102)
  [2]: "CGTACGT"   (after variant, ref pos 103+)

Haplotype 1 (insert A):
  Allele: "TA" (reference T + inserted A)
  Sequence: "ACGT" + "TA" + "CGTACGT" = "ACGTTACGTACGT" (13 bases, +1)

Haplotype 2 (no insertion):
  Allele: "T" (reference only)
  Sequence: "ACGT" + "T" + "CGTACGT" = "ACGTTCGTACGT" (12 bases, original)
```

**Quality Scores (Haplotype 1):**
```
Original: [30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30]
Expected: [30, 30, 30, 30, 30, Q, 30, 30, 30, 30, 30, 30, 30]
                              ↑
                      Inserted base quality
```

Where `Q` = average of flanking qualities (or default Q30)

---

### Test Case 3: 3bp Insertion

**Original Read:**
```
Sequence: ACGTACGTACGT
Ref Pos:  100-111
```

**Variant:** chr1:105 A→AGGG (3bp insertion of GGG)

**Expected Transformation:**
```
Haplotype 1 (insert GGG):
  Split: "ACGTA" + "AGGG" + "CGTACGT"
  Result: "ACGTAAGGGCGTACGT" (16 bases, +3)

Haplotype 2 (no insertion):
  Split: "ACGTA" + "A" + "CGTACGT"
  Result: "ACGTAACGTACGT" (13 bases, original)
```

**Quality Scores (Haplotype 1):**
```
Original: [30] * 12
Expected: [30, 30, 30, 30, 30, Q1, Q2, Q3, 30, 30, 30, 30, ...]
                                ^^^^^^^^
                           3 inserted qualities
```

---

### Test Case 4: 1bp Deletion

**Original Read:**
```
Sequence: ACGTAACGTACGT
Ref Pos:  100-112
```

**Variant:** chr1:105 TA→T (1bp deletion of A)

**Expected Transformation:**
```
Haplotype 1 (delete A):
  Split: "ACGT" + "T" + "ACGTACGT"
  Result: "ACGTTACGTACGT" (12 bases, -1)

Haplotype 2 (keep A):
  Split: "ACGT" + "TA" + "ACGTACGT"
  Result: "ACGTAACGTACGT" (13 bases, original)
```

**Quality Scores (Haplotype 1):**
```
Original: [30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30]
Expected: [30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30]
                              ^^
                        Deleted quality removed
```

---

### Test Case 5: 3bp Deletion

**Original Read:**
```
Sequence: ACGTAAACGTACGT
Ref Pos:  100-113
```

**Variant:** chr1:105 TAAA→T (3bp deletion of AAA)

**Expected Transformation:**
```
Haplotype 1 (delete AAA):
  Split: "ACGT" + "T" + "CGTACGT"
  Result: "ACGTTCGTACGT" (12 bases, -3)

Haplotype 2 (keep AAA):
  Split: "ACGT" + "TAAA" + "CGTACGT"
  Result: "ACGTAAACGTACGT" (14 bases, original)
```

**Quality Scores (Haplotype 1):**
```
Original: [30] * 14
Expected: [30] * 12  (3 qualities removed)
```

---

## 4. Quality Score Handling

### Insertion: Add Quality Scores

**Strategy:** Average flanking regions

```python
def _fill_insertion_quals(insert_len, left_qual, right_qual, insert_qual=30):
    if len(left_qual) == 0 and len(right_qual) == 0:
        # No flanking data - use default Q30
        return np.full(insert_len, 30, dtype=np.uint8)

    # Average left and right flanking qualities
    flank_quals = np.concatenate([left_qual, right_qual])
    mean_qual = int(np.mean(flank_quals))
    return np.full(insert_len, mean_qual, dtype=np.uint8)
```

**Example:**
```
Left flank:  [35, 30, 28]
Right flank: [32, 30]
Average: (35+30+28+32+30) / 5 = 31
Insert 3 bases: [31, 31, 31]
```

### Deletion: Remove Quality Scores

**Strategy:** Truncate to new length

```python
# Original quality: [30, 30, 30, 30, 30]
# Delete 2 bases from allele
# New quality: [30, 30, 30]
new_qual = qual_part[:allele_len]
```

---

## 5. Edge Cases

### Edge Case 1: Deletion at Read Start
```
Read: ----ACGTACGT
Ref:  AAAAACGTACGT

Variant: chr1:100 AAAA→A (3bp deletion at start)

Expected:
  ref2q_left[100] = 0 (no position before)
  ref2q_right[100] = 0 (first available position)

Result: May fail mapping if variant overlaps unmapped region → return None
```

### Edge Case 2: Insertion at Read End
```
Read: ACGTACGT----
Ref:  ACGTACGTAAAA

Variant: chr1:108 T→TAAAA (3bp insertion at end)

Expected:
  Haplotype with insertion: "ACGTACGTTAAAA"
  Quality: Original + [Q, Q, Q] where Q = mean of left flank only
```

### Edge Case 3: Multiple Overlapping Indels
```
Read: ACGTACGTACGT
Variant 1: chr1:102 T→TA (insertion)
Variant 2: chr1:103 A→- (deletion)

Expected:
  Process sequentially
  Split positions account for cumulative length changes
```

**Current Implementation:** Not fully tested - may need special handling

---

## 6. Validation Criteria

### Python Implementation Must:

✅ **Correctness:**
1. Sequence length changes by expected amount (insertions +N, deletions -N)
2. Haplotype sequences differ only at variant sites
3. Non-variant segments are identical to original read
4. Quality scores have correct length (match sequence length)

✅ **Position Mapping:**
1. `ref2q_left` maps all reference positions (including deletions)
2. `ref2q_right` maps all reference positions (including deletions)
3. For deletions: left and right maps differ (flanking positions)
4. For insertions: both maps point to insertion site

✅ **Quality Scores:**
1. Inserted bases have quality ≥ 20 (acceptable Phred score)
2. Inserted bases use flanking average if available
3. Deleted bases properly removed (no quality array overflow)
4. Non-variant regions preserve original qualities

---

## 7. Test Data Generation

### Minimal Test Dataset

Create synthetic BAM + VCF with known variants:

```python
# test_indel_ground_truth.py

test_cases = [
    {
        "name": "1bp_insertion",
        "ref_seq": "ACGTACGTACGT",
        "variant": ("chr1", 105, "A", "AG"),  # Insert G
        "expected_hap1": "ACGTAGCGTACGT",
        "expected_hap2": "ACGTACGTACGT",
    },
    {
        "name": "3bp_insertion",
        "ref_seq": "ACGTACGTACGT",
        "variant": ("chr1", 105, "A", "AGGG"),
        "expected_hap1": "ACGTAAGGGCGTACGT",
        "expected_hap2": "ACGTAACGTACGT",
    },
    {
        "name": "1bp_deletion",
        "ref_seq": "ACGTAACGTACGT",
        "variant": ("chr1", 105, "TA", "T"),
        "expected_hap1": "ACGTTACGTACGT",
        "expected_hap2": "ACGTAACGTACGT",
    },
    {
        "name": "3bp_deletion",
        "ref_seq": "ACGTAAACGTACGT",
        "variant": ("chr1", 105, "TAAA", "T"),
        "expected_hap1": "ACGTTCGTACGT",
        "expected_hap2": "ACGTAAACGTACGT",
    },
]
```

---

## 8. Rust Implementation Acceptance

### Before Accepting Rust Code:

1. **Output Comparison:**
   - Run Python implementation on test cases
   - Run Rust implementation on same test cases
   - Sequences must match **byte-for-byte**
   - Quality scores must match **exactly**

2. **Performance Validation:**
   - Rust must be ≥ 3x faster than Python
   - Memory usage must be ≤ Python
   - Handle ≥ 100K reads without issues

3. **Edge Case Handling:**
   - Deletions at read boundaries
   - Insertions at read boundaries
   - Multiple variants per read
   - Reads with no variants (passthrough)

---

## 9. Known Issues and Limitations

### Current Python Implementation:

✅ **Working:**
- Single insertions (tested in `remap_utils.py:get_read_het_data()`)
- Single deletions (tested with `_build_ref2read_maps()`)
- Quality score generation for insertions
- Quality score truncation for deletions

⚠️ **Uncertain:**
- Multiple overlapping indels on same read
- Very large indels (>50bp)
- Complex indels (combination of ins + del)

❌ **Not Supported:**
- Structural variants (>1kb)
- MNPs (multi-nucleotide polymorphisms)

### Rust Implementation (to be built):

Target: Match Python behavior exactly for known-working cases

---

## 10. Validation Script Specification

### Script: `validate_indel_handling.py`

**Purpose:** Test Python implementation against ground truth

**Steps:**
1. Create synthetic BAM with test reads
2. Create synthetic VCF with test variants
3. Run `write_remap_bam()` with `include_indels=True`
4. Parse output FASTQ files
5. Compare to expected sequences
6. Report pass/fail for each test case

**Output:**
```
Testing indel handling...
✅ 1bp_insertion: PASS (sequences match)
✅ 3bp_insertion: PASS (sequences match)
✅ 1bp_deletion: PASS (sequences match)
✅ 3bp_deletion: PASS (sequences match)
❌ 10bp_insertion: FAIL (quality score mismatch)

Summary: 4/5 tests passed
```

---

## 11. Success Criteria

Before proceeding with Rust implementation:

1. ✅ Python implementation passes all basic test cases (1bp, 3bp, 10bp indels)
2. ✅ Quality score handling validated
3. ✅ Edge cases documented (even if not fully handled)
4. ✅ Ground truth sequences confirmed

After Rust implementation:

1. ✅ Rust output matches Python byte-for-byte
2. ✅ Performance improvement ≥ 3x
3. ✅ All edge cases handled identically
4. ✅ Integration tests pass

---

## 12. References

### Key Files:
- **Python Implementation:** `src/mapping/remap_utils.py` (lines 89-404)
  - `_build_ref2read_maps()`: Position mapping
  - `_fill_insertion_quals()`: Quality score generation
  - `make_phased_seqs_with_qual()`: Sequence + quality construction

- **Rust Target:** `rust/src/bam_remapper.rs`
  - Will add similar functions in Rust

### Algorithms:
- **Two-pass CIGAR walking:** Forward for left map, backward for right map
- **Quality averaging:** Mean of flanking regions for insertions
- **Quality truncation:** Simple slice for deletions

---

## Next Steps

1. **Create validation script** (`validate_indel_handling.py`)
2. **Run Python validation** to confirm ground truth
3. **Document any failures** or unexpected behavior
4. **Fix Python issues** if found
5. **Proceed to Rust implementation** once Python validated

---

**Status:** ✅ SPECIFICATION COMPLETE
**Ready for:** Python validation testing
**Blocks:** Rust indel implementation (Phase 1-4)
