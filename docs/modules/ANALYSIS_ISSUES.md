# Analysis Module - Technical Debt & Issues

**Phase 1.3 Issue Catalog**
**Total Issues:** 18
**Critical:** 3 | **High:** 6 | **Medium:** 7 | **Low:** 2

---

## Critical Issues (Fix Immediately)

### A1: Dead Code Bloat (300+ LOC)
**Severity:** Critical
**Files:** as_analysis.py, compare_ai.py
**Impact:** Code maintainability, confusion

**Description:**
Massive amounts of commented-out old implementations still in codebase.

**Evidence:**
```python
# as_analysis.py:205-243 (39 lines of old single_model)
# as_analysis.py:340-370 (31 lines of old linear_model)
# as_analysis.py:496-556 (61 lines of old get_imbalance_sc)
# compare_ai.py:316-516 (200 lines of get_compared_imbalance_diff_snps v0)
# run_analysis.py:109-166 (58 lines of old WaspAnalysisData class)
```

**Fix:** Delete all commented-out code. Total removal: ~389 lines

---

### A2: Debug Print Statements in Production
**Severity:** Critical
**Files:** as_analysis_sc.py, run_compare_ai.py
**Impact:** Production logs, performance

**Evidence:**
```python
# as_analysis_sc.py:89
print(disp) # DEEBUG BY SHOWING DISP

# run_compare_ai.py:39-40
print(*vars(ai_files).items(), sep="\n") # For debugging
print(adata_inputs) # For debugging
```

**Fix:** Remove debug prints or use proper logging framework

---

### A3: Same Bug Pattern as Counting Module
**Severity:** Critical
**Files:** __main__.py:133, 171, 209
**Impact:** TypeError crashes

**Description:**
Same `if len(groups) > 0:` bug without None check found in counting module.

**Evidence:**
```python
# Line 133 (find_imbalance_sc)
if len(groups) > 0:  # TypeError if groups=None!
    groups=groups[0]

# Line 171 (compare_imbalance)
if len(groups) > 0:  # Same bug
    groups=groups[0]

# Line 209 (run_analysis main wrapper)
if len(samples) > 0:  # Same bug
    samples=samples[0]
```

**Fix:** Add None check: `if groups is not None and len(groups) > 0:`

---

## High Priority Issues

### A4: Duplicate Allele Counting Code
**Severity:** High
**Files:** count_alleles.py, count_alleles_sc.py
**Impact:** Code duplication, maintenance burden

**Description:**
Analysis module has its own copies of counting functions (244 LOC duplicated from counting module).

**Fix:** Import from counting module instead of duplicating

---

### A5: Inconsistent FDR Correction
**Severity:** High
**Files:** as_analysis.py, as_analysis_sc.py, compare_ai.py
**Impact:** Statistical correctness

**Evidence:**
```python
# as_analysis.py:410 - Uses custom bh_correction()
df["fdr_pval"] = bh_correction(df)

# as_analysis_sc.py:256 - Uses scipy directly
df["fdr_pval"] = false_discovery_control(df["pval"], method="bh")

# compare_ai.py:311 - Uses scipy
df["fdr_pval"] = false_discovery_control(df["pval"], method="bh")

# compare_ai.py:513 - Uses custom bh_correction()
df = bh_correction(df)
```

**Fix:** Standardize on scipy.stats.false_discovery_control everywhere

---

### A6: No Input Validation
**Severity:** High
**Files:** run_analysis.py, filter_data.py
**Impact:** Silent failures, confusing errors

**Evidence:**
```python
# run_analysis.py:50-52 - Accepts invalid model silently
if ((self.model is None) or
    (self.model not in {"single", "linear"})):
    self.model = "single"  # No error, just defaults!

# filter_data.py:24-36 - No VCF validation before processing
```

**Fix:** Raise ValueError for invalid inputs instead of silent defaults

---

### A7: Hardcoded Magic Numbers
**Severity:** High
**Files:** filter_data.py
**Impact:** Maintainability, flexibility

**Evidence:**
```python
# filter_data.py:77
usecols=[0, 1, 3, 4, 10, 11, 12]  # What do these mean?

# filter_data.py:95
usecols=[0, 1, 3, 4, 12, 18]  # Different magic numbers
```

**Fix:** Use named constants or column names

---

### A8: Unused Alternative Implementation
**Severity:** High
**Files:** compare_ai.py
**Impact:** Code bloat (200 LOC unused)

**Description:**
`get_compared_imbalance_diff_snps()` and `compare_imbalance_between_groups_diff_snps()` are complete implementations that are never called.

**Fix:** Remove if truly unused, or document why kept as alternative

---

### A9: Memory Inefficiency - No Chunking
**Severity:** High
**Files:** run_analysis_sc.py, run_compare_ai.py
**Impact:** Memory usage on large datasets

**Evidence:**
```python
# Loads entire h5ad into memory
adata = ad.read_h5ad(ai_files.adata_file)  # Could be GB-sized
```

**Fix:** Use AnnData backed mode for large files

---

## Medium Priority Issues

### A10: Inconsistent Parameter Naming
**Severity:** Medium
**Files:** run_analysis_sc.py, run_compare_ai.py
**Impact:** API confusion

**Evidence:**
```python
# run_analysis_sc.py:208
def run_ai_analysis_sc(..., phase=None, ...):  # "phase"

# WaspAnalysisSC.__init__:39
self.phased = phased  # "phased"
```

**Fix:** Standardize on "phased" everywhere

---

### A11: TODO Comments in Production
**Severity:** Medium
**Files:** Multiple
**Impact:** Incomplete features

**Evidence:**
```python
# run_analysis.py:18
# TODO GOTTA IMPLEMENT THIS

# run_analysis.py:44
# TODO parse vcf for phased instead of default unphased

# run_analysis_sc.py:43
# TODO ADD GROUP DISP and other model types

# as_analysis_sc.py:51
# TODO add options to identify and filter GT errors
```

**Fix:** Complete TODOs or create GitHub issues

---

### A12: No Error Handling for File I/O
**Severity:** Medium
**Files:** run_analysis.py, run_analysis_sc.py
**Impact:** Crashes on missing files

**Evidence:**
```python
# run_analysis.py:64 - No try/except
with open(self.count_file) as f:
    count_cols = next(reader(f, delimiter = "\t"))

# run_analysis_sc.py:228 - No validation
adata = ad.read_h5ad(ai_files.adata_file)
```

**Fix:** Add try/except with helpful error messages

---

### A13: Commented-Out Optimizations
**Severity:** Medium
**Files:** as_analysis_sc.py, count_alleles_sc.py
**Impact:** Performance

**Evidence:**
```python
# as_analysis_sc.py:95, 104, 108 - SparseArray commented out
# allele_counts.append(SparseArray(np.zeros(num_cols), fill_value=0))
allele_counts.append(np.zeros(num_cols, dtype=np.int32))

# count_alleles_sc.py:183 - Converts to Sparse later anyway
df = df.astype({group: "Sparse[int]" for group in cols})
```

**Fix:** Either use SparseArray throughout or remove commented code

---

### A14: Mutable Default Arguments
**Severity:** Medium
**Files:** Multiple
**Impact:** Potential bugs

**Note:** Actually handled correctly with `if param is None:` checks, but still anti-pattern

---

### A15: Inconsistent Return Types
**Severity:** Medium
**Files:** as_analysis_sc.py, compare_ai.py
**Impact:** API confusion

**Evidence:**
```python
# as_analysis_sc.py:157 - Returns dict of DataFrames
return df_dict

# compare_ai.py:194 - Also returns dict but with tuple keys
df_dict[(group1, group2)] = df
```

**Fix:** Document return types clearly in docstrings

---

### A16: Z-Score Cutoff Has No Default
**Severity:** Medium
**Files:** run_analysis_sc.py
**Impact:** Confusing behavior

**Evidence:**
```python
# WaspAnalysisSC.__init__:40
self.z_cutoff = z_cutoff # Should i default to something like 4 or 5?

# adata_count_qc:22
if z_cutoff is None and gt_error is None:
    return adata  # No QC if not specified
```

**Fix:** Set sensible default (z=3 or 4) with option to disable

---

## Low Priority Issues

### A17: Pandas SettingWithCopyWarning Risk
**Severity:** Low
**Files:** as_analysis_sc.py
**Impact:** Warnings in logs

**Evidence:**
```python
# as_analysis_sc.py:144 - Uses .copy() to avoid warning
adata = adata[adata.obs[sample].isin(['1|0', '0|1', '1/0', '0/1'])].copy()
```

**Fix:** Already using .copy(), but pattern should be consistent

---

### A18: Cryptic Variable Names
**Severity:** Low
**Files:** as_analysis.py, as_analysis_sc.py
**Impact:** Readability

**Evidence:**
```python
mu, rho, disp, ll, lrt  # Statistical notation is fine
a_counter, pile_tup     # Less clear
```

**Fix:** Add inline comments explaining abbreviations

---

## Summary by Category

| Category | Count | Example |
|----------|-------|---------|
| Code Quality | 5 | A1 (dead code), A4 (duplication) |
| Correctness | 3 | A3 (None bug), A5 (FDR inconsistency) |
| Performance | 2 | A9 (memory), A13 (sparse arrays) |
| Maintainability | 4 | A7 (magic numbers), A11 (TODOs) |
| API Design | 3 | A6 (validation), A10 (naming) |
| Documentation | 1 | A15 (return types) |

---

## Estimated Fix Times

| Issue | Effort | Risk |
|-------|--------|------|
| A1 | 30 min | Low (just delete) |
| A2 | 15 min | Low |
| A3 | 30 min | Low |
| A4 | 2 hours | Medium (import refactor) |
| A5 | 1 hour | Low |
| A6 | 2 hours | Low |
| A7 | 1 hour | Low |
| A8 | 30 min | Low (delete or document) |
| A9 | 3 hours | Medium (backed mode) |
| A10 | 1 hour | Low |
| A11 | 2 hours | Medium (depends on TODOs) |
| A12 | 2 hours | Low |
| A13 | 1 hour | Low |
| A14-A18 | 3 hours | Low |

**Total Estimated Effort:** ~20 hours

---

## Recommended Fix Order (4-Week Sprint)

**Week 1: Critical Bugs**
1. A3 - None check bug (30 min)
2. A2 - Debug prints (15 min)
3. A5 - FDR standardization (1 hour)

**Week 2: Code Quality**
4. A1 - Delete dead code (30 min)
5. A8 - Remove unused alt implementation (30 min)
6. A11 - Resolve TODOs (2 hours)

**Week 3: Improvements**
7. A6 - Input validation (2 hours)
8. A7 - Named constants (1 hour)
9. A10 - Consistent naming (1 hour)
10. A12 - Error handling (2 hours)

**Week 4: Optimizations**
11. A4 - Deduplicate counting code (2 hours)
12. A9 - Memory optimization (3 hours)
13. A13-A18 - Remaining issues (4 hours)
