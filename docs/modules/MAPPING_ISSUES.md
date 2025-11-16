# Mapping Module - Technical Debt & Issues

**Phase 1.4 Issue Catalog**
**Total Issues:** 12
**Critical:** 2 | **High:** 4 | **Medium:** 4 | **Low:** 2

---

## Critical Issues

### M1: Same None-Check Bug
**Severity:** Critical
**File:** __main__.py:86
**Impact:** TypeError crash

**Evidence:**
```python
# Line 86
if len(samples) > 0:  # TypeError if samples=None!
    samples=samples[0]
else:
    samples=None
```

**Fix:** `if samples is not None and len(samples) > 0:`
**Note:** Same bug found in counting and analysis modules

---

### M2: Single-End and Unphased Not Implemented
**Severity:** Critical
**File:** run_mapping.py:75-82
**Impact:** Limits usability

**Evidence:**
```python
if not wasp_files.is_paired:
    raise ValueError("Single-End not Implemented")

if not wasp_files.is_phased:
    raise ValueError("Unphased not Implemented")

if wasp_files.samples is None:
    raise ValueError("Zero samples not supported yet")
```

**Fix:** Implement unphased support (uses DP like analysis module) or document as permanent limitation

---

## High Priority Issues

### M3: Dead Commented Code (200+ LOC)
**Severity:** High
**Files:** make_remap_reads.py, intersect_variant_data.py
**Impact:** Code bloat, maintainability

**Evidence:**
```python
# make_remap_reads.py:396-499 (104 lines of old swap_chrom_alleles)
# intersect_variant_data.py:100-236 (137 lines of old make_intersect_df/process_bam)
```

**Fix:** Remove all commented-out code (~241 lines total)

---

### M4: Subprocess Usage Instead of Pysam
**Severity:** High
**Files:** make_remap_reads.py, intersect_variant_data.py
**Impact:** Inconsistency, error handling

**Evidence:**
```python
# make_remap_reads.py:238-245
collate_cmd = ["samtools", "collate", "-u", "-O", out_bam]
fastq_cmd = ["samtools", "fastq", "-1", r1_out, "-2", r2_out]
collate_process = subprocess.run(collate_cmd, stdout=PIPE, check=True)
fastq_process = subprocess.run(fastq_cmd, input=collate_process.stdout, check=True)

# intersect_variant_data.py - Uses bcftools via subprocess
```

**Fix:** Standardize on either subprocess or pysam wrappers, not mixed

---

### M5: No Input Validation
**Severity:** High
**Files:** run_mapping.py, wasp_data_files.py
**Impact:** Silent failures, confusing errors

**Evidence:**
```python
# wasp_data_files.py:75-76 - Regex split could fail
vcf_prefix = re.split(r'.vcf|.bcf', Path(self.vcf_file).name)[0]
bam_prefix = Path(self.bam_file).name.rsplit(".bam")[0]

# No check if files exist before processing
```

**Fix:** Add file existence checks, validate formats

---

### M6: Memory Risk - Loading All Read Names
**Severity:** High
**File:** intersect_variant_data.py:84-86
**Impact:** Memory usage on large BAMs

**Evidence:**
```python
unique_reads = np.unique(
    [read.query_name for read in bam.fetch(until_eof=True)])
file.write("\n".join(unique_reads))
```

**Fix:** Stream read names instead of loading all into memory

---

## Medium Priority Issues

### M7: TODOs in Production
**Severity:** Medium
**Files:** Multiple
**Impact:** Incomplete features

**Evidence:**
```python
# wasp_data_files.py:11
# TODO, GOTTA INCLUDE ALL POSSIBLE DATA COMBOS

# wasp_data_files.py:51
# TODO GOTTA FIX THIS TO CHECK IF PHASED

# wasp_data_files.py:60
# TODO GOTTA WARN UNPHASED BAD

# wasp_data_files.py:61
# TODO WARN SOME UNPHASED WHILE OTHERS PHASED

# wasp_data_files.py:67
# TODO handle temp loc, maybe make default if temp not made?

# remap_utils.py:78
# TODO MULTISAMP AND MAX SEQS

# intersect_variant_data.py:68
# TODO FIX ALL OF THESE TO USE A CLASS
```

**Fix:** Complete TODOs or create GitHub issues

---

### M8: Hardcoded Column Indices
**Severity:** Medium
**File:** intersect_variant_data.py:251
**Impact:** Fragile code

**Evidence:**
```python
subset_cols = [df.columns[i] for i in np.r_[0, 3, 1, 2, -num_samps:0]]
```

**Fix:** Use named constants or column names

---

### M9: Temporary Workaround in Production
**Severity:** Medium
**File:** wasp_data_files.py:68-70
**Impact:** Temp files not cleaned up

**Evidence:**
```python
# Temporary workaround until figure out temp dir options
if self.temp_loc is None:
    self.temp_loc = self.out_dir  # Uses output dir instead of temp!
```

**Fix:** Implement proper temp directory handling

---

### M10: Bug Workaround Instead of Fix
**Severity:** Medium
**File:** make_remap_reads.py:299-302
**Impact:** Silent data loss possible

**Evidence:**
```python
else:
    # TEMPORARY FIX FOR BUG????
    # NOT SURE WHY SOME READS WOULD SHOW UP BUT NOT OVERLAP A SNP
    continue
```

**Fix:** Investigate root cause instead of silently skipping

---

## Low Priority Issues

### M11: Commented Debug Prints
**Severity:** Low
**Files:** Multiple
**Impact:** Code clutter

**Evidence:**
```python
# run_mapping.py:72 - print(*vars(wasp_files).items(), sep="\n")
# make_remap_reads.py:221 - print(f"{chrom}: Processed...")
# filter_remap_reads.py:62-66 - Multiple print statements
# intersect_variant_data.py - Many commented prints
```

**Fix:** Remove or convert to proper logging

---

### M12: Inconsistent Error Messages
**Severity:** Low
**Files:** run_mapping.py
**Impact:** User experience

**Evidence:**
```python
raise ValueError("Single-End not Implemented")
raise ValueError("Unphased not Implemented")
raise ValueError("Zero samples not supported yet")
```

**Fix:** Standardize error messages with helpful suggestions

---

## Summary by Category

| Category | Count | Example |
|----------|-------|---------|
| Correctness | 2 | M1 (None bug), M10 (bug workaround) |
| Code Quality | 3 | M3 (dead code), M7 (TODOs), M11 (prints) |
| Implementation Gaps | 1 | M2 (single-end/unphased unsupported) |
| Performance | 1 | M6 (memory risk) |
| Maintainability | 3 | M4 (subprocess), M8 (hardcoded), M9 (workaround) |
| Validation | 1 | M5 (no input validation) |
| User Experience | 1 | M12 (error messages) |

---

## Estimated Fix Times

| Issue | Effort | Risk |
|-------|--------|------|
| M1 | 15 min | Low |
| M2 | 2-3 days | High (requires algorithm implementation) |
| M3 | 30 min | Low (delete code) |
| M4 | 2 hours | Medium (standardization) |
| M5 | 2 hours | Low |
| M6 | 1 hour | Low (streaming) |
| M7 | 3 hours | Medium (depends on TODOs) |
| M8 | 1 hour | Low |
| M9 | 1 hour | Low |
| M10 | 2 hours | Medium (investigation) |
| M11 | 30 min | Low |
| M12 | 30 min | Low |

**Total Effort (excluding M2):** ~14 hours
**With M2 (unphased support):** ~3-4 days

---

## Recommended Fix Order

**Week 1: Critical Bugs**
1. M1 - None-check bug (15 min)
2. M5 - Input validation (2 hours)

**Week 2: Code Quality**
3. M3 - Delete dead code (30 min)
4. M11 - Remove debug prints (30 min)
5. M7 - Resolve TODOs (3 hours)

**Week 3: Improvements**
6. M4 - Standardize subprocess usage (2 hours)
7. M6 - Stream read names (1 hour)
8. M8 - Named constants (1 hour)
9. M9 - Proper temp handling (1 hour)
10. M10 - Investigate bug workaround (2 hours)
11. M12 - Better error messages (30 min)

**Future: Feature Implementation**
12. M2 - Unphased/single-end support (separate epic)

---

## Key Insight

Mapping module is well-structured with decorator patterns and clear separation of concerns, but has **3 major limitations**:
1. Paired-end only
2. Phased VCF only
3. Sample required

These are documented but limit usability compared to original WASP tool.
