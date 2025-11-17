# mypy Type Errors Analysis

Complete analysis of all type errors found after adding type hints to counting and mapping modules.

## Summary

- **Total Errors**: 16
- **Counting Module**: 5 errors in 2 files
- **Mapping Module**: 11 errors in 3 files

## Severity Classification

### 🔴 CRITICAL - Will Crash at Runtime (1 error)
These cause actual crashes if the code path is executed.

### 🟡 WARNING - Type Safety Issues (13 errors)
Won't crash but indicate incorrect type usage that could fail with certain inputs.

### 🟢 INFO - Missing Annotations (2 errors)
Just incomplete type hints, no functional impact.

---

## Counting Module Errors (5 total)

### File: `src/counting/filter_variant_data.py`

#### 🔴 CRITICAL: Line 210 - Undefined Variable
```python
raise ValueError(f"Could not recognize BED format. Expected 3-6 columns, got {n_cols} columns")
```

**Error**: `Name "n_cols" is not defined`

**Issue**: Variable `n_cols` is used but never defined in the function.

**Impact**: **WILL CRASH** with `NameError` if this error path is hit.

**Fix**: Should be `len(df.columns)`:
```python
raise ValueError(f"Could not recognize BED format. Expected 3-6 columns, got {len(df.columns)} columns")
```

**Root Cause**: Copy-paste error or incomplete refactoring.

---

### File: `src/counting/__main__.py`

#### 🟢 INFO: Line 19 - Missing Return Type
```python
def count_variants(
    bam: Annotated[str, typer.Argument(help="Bam File")],
    ...
```

**Error**: `Function is missing a return type annotation`

**Fix**: Add `-> None`:
```python
def count_variants(...) -> None:
```

**Impact**: None - just incomplete type hints.

---

#### 🟡 WARNING: Line 119 - Type Mismatch
```python
samples: Optional[List[str]]  # Declared type

if samples is not None and len(samples) > 0:
    samples = samples[0]  # ❌ Assigning str to List[str] | None
```

**Error**: `Incompatible types in assignment (expression has type "str", variable has type "list[str] | None")`

**Issue**: `samples` is typed as `List[str] | None` but code assigns a single string to it.

**Impact**: Type checker confusion. Code works at runtime because Python is dynamic, but violates type contract.

**Fix Option 1** - Use separate variable:
```python
samples_list: Optional[List[str]] = samples
sample_str: Optional[str] = None

if samples_list is not None and len(samples_list) > 0:
    sample_str = samples_list[0]
```

**Fix Option 2** - Change parameter type:
```python
samples: Optional[str]  # If it's always a single sample
```

**Root Cause**: Typer CLI returns `List[str]` but internal logic treats it as single string.

---

#### 🟢 INFO: Line 143 - Missing Return Type
Same as line 19 - missing `-> None`.

---

#### 🟡 WARNING: Line 204 - Type Mismatch
Same issue as line 119 - different function, same pattern.

---

## Mapping Module Errors (11 total)

### File: `src/mapping/wasp_data_files.py`

#### 🟡 WARNING: Line 64 - None Check Missing
```python
self.samples: Optional[List[str]]  # Can be None

# Later in code:
samps_phased = [vcf_samps[s].phased for s in self.samples]  # ❌ Might be None!
```

**Error**: `Item "None" of "list[str] | None" has no attribute "__iter__" (not iterable)`

**Issue**: Code iterates over `self.samples` without checking if it's None.

**Impact**: Would crash with `TypeError: 'NoneType' object is not iterable` if samples is None.

**Fix**: Add None check:
```python
if self.samples is not None:
    samps_phased = [vcf_samps[s].phased for s in self.samples]
```

---

#### 🟡 WARNING: Line 103 - Optional Assignment
```python
self.remap_fq2: str  # Declared as str (not optional)

# Later:
self.remap_fq2 = None  # ❌ Assigning None to str
```

**Error**: `Incompatible types in assignment (expression has type "None", variable has type "str")`

**Issue**: For single-end reads, `remap_fq2` is set to None, but type says it's always a string.

**Fix**: Change type to Optional:
```python
self.remap_fq2: Optional[str]
```

---

#### 🟡 WARNING: Line 114 - Path with None
```python
out_dir: Optional[Union[str, Path]]  # Can be None

Path(self.out_dir)  # ❌ Path() doesn't accept None
```

**Error**: `Argument 1 to "Path" has incompatible type "str | Path | None"; expected "str | PathLike[str]"`

**Issue**: `out_dir` could be None, but Path() constructor requires a valid path.

**Fix**: Ensure out_dir is set before use, or provide default:
```python
if self.out_dir is None:
    self.out_dir = "."
Path(self.out_dir)
```

---

### File: `src/mapping/run_mapping.py`

#### 🟡 WARNING: Line 94 - Same Path(None) Issue
Same as wasp_data_files.py:114

---

#### 🟡 WARNING: Line 98 - Union[str, Path] vs str
```python
wasp_files.vcf_file: Union[str, Path]

vcf_to_bed(vcf_file=wasp_files.vcf_file, ...)  # Function expects str only
```

**Error**: `Argument "vcf_file" to "vcf_to_bed" has incompatible type "str | Path"; expected "str"`

**Issue**: Function signatures don't match - caller has Union type, callee expects only str.

**Fix Option 1** - Convert to str:
```python
vcf_to_bed(vcf_file=str(wasp_files.vcf_file), ...)
```

**Fix Option 2** - Update function signature to accept Union:
```python
def vcf_to_bed(vcf_file: Union[str, Path], ...) -> str:
```

---

#### 🟡 WARNING: Line 100 - samples Type Mismatch
```python
samples: Union[str, List[str]]  # Can be single string or list

vcf_to_bed(..., samples=wasp_files.samples)  # Expects List[str] | None
```

**Error**: `Argument "samples" to "vcf_to_bed" has incompatible type "str | list[str]"; expected "list[str] | None"`

**Issue**: Caller allows `str`, callee expects only `List[str]` or None.

**Fix**: Normalize to list in WaspDataFiles:
```python
if isinstance(samples, str):
    self.samples = [samples]
```

---

#### 🟡 WARNING: Line 103 - Same Union[str, Path] Issue
Same as line 98.

---

#### 🟡 WARNING: Line 125 - Same samples Issue
Same as line 100.

---

### File: `src/mapping/__main__.py`

#### 🟡 WARNING: Line 87 - len() on Optional
```python
samples: Optional[List[str]]

if len(samples) > 1:  # ❌ samples might be None!
```

**Error**: `Argument 1 to "len" has incompatible type "list[str] | None"; expected "Sized"`

**Fix**: Add None check:
```python
if samples is not None and len(samples) > 1:
```

---

#### 🟡 WARNING: Line 88 - Indexing Optional
```python
samples: Optional[List[str]]

samples = [samples[0]]  # ❌ Can't index None!
```

**Error**: `Value of type "list[str] | None" is not indexable`

**Fix**: Add None check:
```python
if samples is not None:
    samples = [samples[0]]
```

---

#### 🟡 WARNING: Line 88 - Type Assignment Mismatch
```python
samples: Optional[List[str]]

samples = [samples[0]]  # Type of samples[0] is str | Any, but samples expects List[str] | None
```

**Error**: `Incompatible types in assignment (expression has type "str | Any", variable has type "list[str] | None")`

**Issue**: Complex type inference issue with the reassignment.

---

## Summary by Category

### 🔴 Critical Bugs (Will Crash): 1
- `filter_variant_data.py:210` - Undefined variable `n_cols`

### 🟡 Type Safety Issues (Could Crash): 13
- Missing None checks: 4 errors
- Union type mismatches: 7 errors
- Optional type mismatches: 2 errors

### 🟢 Missing Annotations (Cosmetic): 2
- Missing return types: 2 errors

## Recommended Fix Priority

### High Priority (Fix Now)
1. **filter_variant_data.py:210** - Fix undefined `n_cols` (actual bug!)
2. **wasp_data_files.py:64** - Add None check for samples iteration
3. **mapping/__main__.py:87-88** - Add None checks before len() and indexing

### Medium Priority (Improves Type Safety)
4. Normalize Union[str, List[str]] → List[str] in constructors
5. Fix Optional[str] assignments (remap_fq2)
6. Convert Union[str, Path] to str when passing to functions

### Low Priority (Cosmetic)
7. Add missing `-> None` return type annotations

## Test Coverage Analysis

**Why don't tests catch these?**

Most errors are in edge cases not covered by current tests:
- `n_cols` error is in an uncommon BED format branch
- None checks fail when optional parameters are actually None
- Union type issues work in Python's dynamic runtime but fail static checks

**Recommendation**: These type errors reveal gaps in test coverage. Consider adding tests for:
- Edge case BED formats
- None values for optional parameters
- Mixed type inputs (str vs List[str])
