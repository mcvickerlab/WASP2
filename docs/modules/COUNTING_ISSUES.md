# Counting Module - Technical Debt & Issues Inventory

**Module**: `src/counting/`
**Generated**: 2025-11-15 (Phase 1.2)
**Total Issues**: 24 catalogued

---

## Issue Summary

| Severity | Count | Description |
|----------|-------|-------------|
| 🔴 **Critical** | 3 | Bugs that affect functionality or performance |
| 🟡 **High** | 5 | Code quality issues, major duplication |
| 🟢 **Medium** | 8 | Error handling, minor improvements |
| ⚪ **Low** | 8 | TODOs, cleanup, documentation |

---

## 🔴 Critical Issues

### C1: Binary Search Not Used (Performance Bug)

**File**: `src/counting/count_alleles.py`
**Lines**: 107-120
**Severity**: 🔴 **CRITICAL** (Performance)

**Problem**:
```python
# Line 109: TODO Update with binary search
for qpos, refpos in read.get_aligned_pairs(True):  # O(n) linear search!
    if refpos == pos-1:
        if seq[qpos] == ref:
            ref_count+=1
        ...
```

**Description**:
- Uses O(n) linear search through all aligned pairs
- Binary search function `find_read_aln_pos()` EXISTS in same file but NOT USED!
- Single-cell version (count_alleles_sc.py:209) DOES use binary search
- Inconsistency between bulk and single-cell implementations

**Impact**:
- Slower performance on reads with many indels
- Unnecessary CPU cycles

**Fix**:
```python
# Use existing binary search function:
qpos = find_read_aln_pos(read, pos-1)
if qpos is not None:
    if seq[qpos] == ref:
        ref_count += 1
    elif seq[qpos] == alt:
        alt_count += 1
    else:
        other_count += 1
```

**Effort**: 🟢 Low (10 lines of code)
**Priority**: High (easy performance win)

---

### C2: AnnData Dimensions Transposed

**File**: `src/counting/count_alleles_sc.py`
**Lines**: 105-117, 128-145
**Severity**: 🔴 **CRITICAL** (Compatibility)

**Problem**:
```python
# Line 108: Creates (SNPs × Cells) matrix
sparse_ref = csr_matrix(
    ...,
    shape=(snp_df.shape[0], len(bc_dict)),  # (n_snps, n_cells)
    dtype=np.uint8
)

# Line 139-145: Assigns incorrectly
adata.obs = snp_df.to_pandas()      # SNPs as observations (WRONG!)
adata.var_names = bc_dict.keys()    # Cells as variables (WRONG!)
```

**Standard Convention**:
- Scanpy expects: `adata.X` shape = `(cells, genes)`
- SnapATAC2 expects: `adata.X` shape = `(cells, peaks/bins)`
- **Current**: `adata.X` shape = `(SNPs, cells)` ← TRANSPOSED!

**Impact**:
- Breaks compatibility with scanpy/SnapATAC2
- Downstream analysis tools may fail or give wrong results
- Users would need to manually transpose

**Fix**:
```python
# Transpose all matrices:
sparse_ref = csr_matrix(..., shape=(len(bc_dict), snp_df.shape[0]))
# (n_cells, n_snps) ✓

# Correct assignment:
adata.obs_names = bc_dict.keys()  # Cells as observations ✓
adata.var = snp_df.to_pandas()     # SNPs as variables ✓
```

**Effort**: 🟡 Medium (need to verify downstream analysis still works)
**Priority**: **CRITICAL** (breaks standard workflows)

**Validation Needed**: Test with scanpy and SnapATAC2

---

### C3: Sample Parsing Bug (TypeError)

**File**: `src/counting/__main__.py`
**Lines**: 118-121, 203-206
**Severity**: 🔴 **CRITICAL** (Runtime Error)

**Problem**:
```python
# Line 118: TypeError if samples is None!
if len(samples) > 0:  # ← Crashes when samples=None
    samples=samples[0]
else:
    samples=None
```

**Description**:
- Calls `len()` on potentially `None` value
- Raises `TypeError: object of type 'NoneType' has no len()`
- Same bug in both `count_variants()` and `count_variants_sc()`

**Trigger**:
```bash
python -m src.counting count-variants data.bam variants.vcf
# ← No --samples flag → samples=None → CRASH
```

**Fix**:
```python
if samples is not None and len(samples) > 0:
    samples = samples[0]
else:
    samples = None
```

**Effort**: 🟢 Trivial (one-line fix, two locations)
**Priority**: High (breaks basic usage)

**Note**: Unclear why only first sample is used - may indicate design issue

---

## 🟡 High Priority Issues

### H1: Massive Code Duplication (80%)

**Files**: `run_counting.py` vs `run_counting_sc.py`
**Lines**: `WaspCountFiles` (16-104) vs `WaspCountSC` (19-116)
**Severity**: 🟡 **HIGH** (Maintainability)

**Duplication Examples**:

| Code | WaspCountFiles | WaspCountSC | Identical? |
|------|----------------|-------------|------------|
| Sample parsing | Lines 39-48 | Lines 51-60 | ✅ 100% |
| VCF prefix parsing | Lines 61-62 | Lines 73-74 | ✅ 100% |
| File extension regex | Lines 73-93 | Lines 86-106 | ✅ 95% |
| Temp directory handling | Lines 57-58 | Lines 68-69 | ✅ 100% |

**Only Differences**:
1. Parameter name: `region_file` vs `feature_file`
2. Parameter name: `barcode_file` (new in SC)
3. Default output: `.tsv` vs `.h5ad`

**Impact**:
- Bug fixes must be applied twice
- Increases maintenance burden
- Violates DRY principle

**Fix**: Extract common base class
```python
class WaspCountFilesBase:
    def __init__(self, bam_file, vcf_file, ...):
        # Common logic

class WaspCountFiles(WaspCountFilesBase):
    # Bulk-specific

class WaspCountSC(WaspCountFilesBase):
    def __init__(self, ..., barcode_file, ...):
        super().__init__(...)
        self.barcode_file = barcode_file
```

**Effort**: 🟡 Medium (refactor + test)
**Priority**: High (reduces future bugs)

---

### H2: Dead Code (Never Called)

**Files**: `filter_variant_data.py`, `parse_gene_data.py`
**Severity**: 🟡 **HIGH** (Code Bloat)

**Dead Functions**:

| Function | File | Lines | Reason |
|----------|------|-------|--------|
| `parse_intersect_region()` | filter_variant_data.py | 180-236 | Old version, replaced by `_new()` |
| `parse_intersect_genes()` | parse_gene_data.py | 156-179 | Old version, replaced by `_new()` |
| `gtf_to_bed()` | filter_variant_data.py | 75-105 | Never called, GTF handled in parse_gene_data |

**Evidence**:
```bash
# grep shows no callers:
$ grep -r "parse_intersect_region(" src/
# → Only definition, no calls!

$ grep -r "gtf_to_bed(" src/
# → Only definition, no calls!
```

**Impact**:
- 123 lines of unused code
- Confuses developers ("Which version to use?")
- Increases cognitive load

**Fix**: Delete dead functions

**Effort**: 🟢 Low (delete + verify tests pass)
**Priority**: Medium (cleanup)

---

### H3: No Error Handling (Print Instead of Raise)

**Files**: Multiple
**Severity**: 🟡 **HIGH** (User Experience)

**Examples**:

```python
# run_counting.py:93
print("invalid ftype")  # Should raise ValueError!

# parse_gene_data.py:63
print(f"exon, gene or transcript not found...")  # Should raise!

# parse_gene_data.py:80
print(f"No 'ID', '{feature}_id' or 'Name' attribute found...")  # Should raise!

# filter_variant_data.py:211
print("COULD NOT RECOGNIZE FORMAT OR WRONG NUMBER OF COLS")  # Should raise!

# count_alleles.py:52
print(f"Skipping {chrom}: Contig not found\n")  # Warning is OK
```

**Impact**:
- Errors are silent (just print, continue running)
- Users don't know why output is wrong
- Hard to debug failures

**Fix**: Raise appropriate exceptions
```python
# Instead of:
print("invalid ftype")

# Use:
raise ValueError(f"Unsupported file type: {f_ext}. "
                 f"Expected: .bed, .narrowPeak, .gtf, .gff, .gff3")
```

**Effort**: 🟢 Low (4-5 locations)
**Priority**: High (better UX)

---

### H4: Unclear Read Set Behavior

**Files**: `count_alleles.py`, `count_alleles_sc.py`
**Lines**: count_alleles.py:88,100-103 / count_alleles_sc.py:181,190-192
**Severity**: 🟡 **HIGH** (Logic Ambiguity)

**Problem**:
```python
def count_snp_alleles(bam, chrom, snp_list):
    read_set = set()  # ← OUTSIDE SNP loop!

    for pos, ref, alt in snp_list:  # Multiple SNPs
        for read in bam.fetch(chrom, pos-1, pos):
            if read.query_name in read_set:
                continue  # ← Skip read at ALL subsequent SNPs!
            read_set.add(read.query_name)
```

**Behavior**: A read spanning 2 SNPs is only counted at the FIRST SNP encountered

**Questions**:
1. Is this intentional?
2. Should reads be counted at each SNP they overlap?
3. Or only once per chromosome?

**Evidence of Confusion**:
```python
# Line 93 (commented out):
# read_set = set()
# ↑ Was this meant to be INSIDE the SNP loop?
```

**Impact**:
- Unclear specification
- Potentially incorrect counts

**Fix Needed**:
1. Document intended behavior
2. Add test case for reads spanning multiple SNPs
3. Either:
   - Move `read_set` inside loop (count at each SNP), OR
   - Keep as-is but document (count once per chromosome)

**Effort**: 🟡 Medium (needs design decision + tests)
**Priority**: High (affects correctness)

---

### H5: Debug Code in Production

**File**: `src/counting/run_counting_sc.py`
**Line**: 140
**Severity**: 🟡 **HIGH** (Code Quality)

**Problem**:
```python
print(*vars(count_files).items(), sep="\n")  # For debugging
```

**Description**:
- Prints all class attributes to stdout
- Not behind any debug flag
- Runs in production

**Output**:
```
('bam_file', '/path/to/data.bam')
('vcf_file', '/path/to/variants.vcf')
('barcode_file', '/path/to/barcodes.txt')
...
```

**Impact**:
- Clutters output
- Unprofessional
- May leak paths in logs

**Fix**: Remove or use proper logging
```python
import logging
logger = logging.getLogger(__name__)
logger.debug(f"WaspCountSC config: {vars(count_files)}")
```

**Effort**: 🟢 Trivial
**Priority**: Medium (polish)

---

## 🟢 Medium Priority Issues

### M1: No Type Hints

**Files**: All files
**Severity**: 🟢 **MEDIUM** (Developer Experience)

**Problem**: No function signatures have type hints

**Example**:
```python
# Current:
def count_snp_alleles(bam, chrom, snp_list):
    ...

# Should be:
def count_snp_alleles(
    bam: AlignmentFile,
    chrom: str,
    snp_list: Iterator[Tuple[int, str, str]]
) -> List[Tuple[str, int, int, int, int]]:
    ...
```

**Impact**:
- No IDE autocomplete
- No static type checking
- Harder to understand API

**Fix**: Add type hints throughout

**Effort**: 🟡 Medium (1,430 lines to annotate)
**Priority**: Medium (Phase 2 task)

---

### M2: Missing Docstrings

**Files**: All files
**Severity**: 🟢 **MEDIUM** (Documentation)

**Problem**: Many functions lack docstrings

**Examples**:
```python
# Has docstring:
def make_count_df(bam_file, df):
    """
    Make DF containing all intersections and allele counts
    ...
    """

# Missing docstring:
def tempdir_decorator(func):  # ← No docstring!

def find_read_aln_pos(read, pos):  # ← No docstring!

class WaspCountFiles:  # ← No class docstring!
```

**Impact**:
- Harder to understand API
- No auto-generated docs

**Fix**: Add Google/NumPy style docstrings

**Effort**: 🟡 Medium
**Priority**: Medium (Phase 2)

---

### M3: No Progress Indicators

**Files**: `count_alleles.py`, `count_alleles_sc.py`
**Severity**: 🟢 **MEDIUM** (User Experience)

**Problem**: Long-running operations have no progress bar

**Current**:
```python
for chrom in chrom_list:  # Could be 24 chromosomes!
    count_snp_alleles(...)
    # User sees nothing until chrom finishes
```

**Fix**: Add tqdm progress bar
```python
from tqdm import tqdm

for chrom in tqdm(chrom_list, desc="Counting chromosomes"):
    count_snp_alleles(...)
```

**Effort**: 🟢 Low (add tqdm dependency + 2 locations)
**Priority**: Medium (nice UX improvement)

---

### M4: Hard-Coded Constants

**Files**: Multiple
**Severity**: 🟢 **MEDIUM** (Flexibility)

**Examples**:
```python
# filter_variant_data.py:20-21
view_cmd = ["bcftools", "view", str(vcf_file),
            "-m2", "-M2", "-v", "snps", "-Ou"]  # ← Hard-coded filters

# count_alleles_sc.py:109, 116, 123
dtype=np.uint8  # ← Max count = 255, hard-coded
```

**Impact**:
- Cannot change SNP filters without code edit
- uint8 limits counts to 255 (may overflow for high-coverage)

**Fix**: Make configurable
```python
# Configuration file or parameters:
MIN_ALLELES = 2  # -m2
MAX_ALLELES = 2  # -M2
VARIANT_TYPE = "snps"  # -v
COUNT_DTYPE = np.uint16  # Support up to 65,535 counts
```

**Effort**: 🟡 Medium
**Priority**: Low (current defaults are reasonable)

---

### M5: Commented-Out Code

**Files**: Multiple
**Severity**: 🟢 **MEDIUM** (Code Cleanliness)

**Examples**:
```python
# __main__.py:12-13
# app = typer.Typer()
# app = typer.Typer(pretty_exceptions_show_locals=False)

# count_alleles.py:93
# read_set = set()

# run_counting.py:151-153
# temporarily disable for ASE
# if not count_files.is_gene_file:
#     with_gt = True

# run_counting.py:216-219
# df = parse_intersect_region(...)  # Old version

# count_alleles_sc.py:63
# chrom_list = chrom_list[:3] # Testing purposes

# count_alleles_sc.py:164-169
# region_snp_dict = dict(...)  # Old approach
```

**Impact**:
- Clutters code
- Unclear if needed
- Git history already preserves old code

**Fix**: Delete commented code

**Effort**: 🟢 Trivial
**Priority**: Low (cleanup)

---

### M6: Inconsistent Subprocess Error Handling

**Files**: `filter_variant_data.py`
**Severity**: 🟢 **MEDIUM** (Robustness)

**Problem**:
```python
# Some calls use check=True:
subprocess.run(intersect_cmd, stdout=file, check=True)  # ← Good!

# Others don't capture stderr:
subprocess.run(view_cmd, stdout=subprocess.PIPE, check=True)
# ← If bcftools writes error to stderr, user doesn't see it!
```

**Fix**: Consistent error handling
```python
try:
    result = subprocess.run(
        cmd, check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
except subprocess.CalledProcessError as e:
    raise RuntimeError(
        f"Command failed: {' '.join(cmd)}\n"
        f"Return code: {e.returncode}\n"
        f"Error: {e.stderr}"
    ) from e
```

**Effort**: 🟡 Medium (wrap all subprocess calls)
**Priority**: Medium

---

### M7: No Input Validation

**Files**: All orchestrators
**Severity**: 🟢 **MEDIUM** (User Experience)

**Problem**: No validation that input files exist/are readable

**Example**:
```python
def run_count_variants(bam_file, vcf_file, ...):
    # No check if files exist!
    # Will fail later with cryptic pysam error
```

**Fix**: Validate early
```python
def run_count_variants(bam_file, vcf_file, ...):
    # Validate inputs
    bam_path = Path(bam_file)
    if not bam_path.exists():
        raise FileNotFoundError(f"BAM file not found: {bam_file}")
    if not bam_path.is_file():
        raise ValueError(f"BAM path is not a file: {bam_file}")

    # Check BAM is indexed
    bai_path = Path(str(bam_file) + ".bai")
    if not bai_path.exists():
        raise FileNotFoundError(f"BAM index not found: {bai_path}. "
                                f"Run: samtools index {bam_file}")
    # ... validate VCF, etc.
```

**Effort**: 🟡 Medium (add validation for all inputs)
**Priority**: Medium (better error messages)

---

### M8: Regex Duplication

**Files**: `run_counting.py` vs `run_counting_sc.py`
**Severity**: 🟢 **MEDIUM** (Duplication)

**Problem**: Identical regex patterns defined twice

**Example**:
```python
# Both files have:
if re.search(r'\.(.*Peak|bed)(?:\.gz)?$', f_ext, re.I):
    ...
elif re.search(r'\.g[tf]f(?:\.gz)?$', f_ext, re.I):
    ...
elif re.search(r'\.gff3(?:\.gz)?$', f_ext, re.I):
    ...
```

**Fix**: Extract to constants or utility function
```python
# utils.py
FILE_TYPE_PATTERNS = {
    'bed': r'\.(.*Peak|bed)(?:\.gz)?$',
    'gtf': r'\.g[tf]f(?:\.gz)?$',
    'gff3': r'\.gff3(?:\.gz)?$',
}

def detect_file_type(file_path):
    ext = "".join(Path(file_path).suffixes)
    for ftype, pattern in FILE_TYPE_PATTERNS.items():
        if re.search(pattern, ext, re.I):
            return ftype
    return None
```

**Effort**: 🟢 Low
**Priority**: Low (part of H1 refactor)

---

## ⚪ Low Priority Issues (TODOs & Cleanup)

### L1-L13: TODO Comments

**Severity**: ⚪ **LOW** (Documentation)

| File | Line | TODO | Priority |
|------|------|------|----------|
| `__main__.py` | 16 | `# TODO GOTTA TEST THIS` | Medium (implies untested) |
| `__main__.py` | 138-139 | Unit tests for Typer | Medium |
| `__main__.py` | 175 | Implement GTF/GFF for SC | Low (feature request) |
| `run_counting.py` | 15 | Move class to separate file? | Low (preference) |
| `run_counting.py` | 98, 174 | Add auto-parsers & validators | Medium |
| `run_counting.py` | 164 | Parse gene features/attrs | Low (done) |
| `run_counting.py` | 196 | Validate dataframe | Medium |
| `run_counting.py` | 221 | Include filtered BAM step? | Low (design decision) |
| `run_counting.py` | 230 | Return for analysis pipeline? | Low (design decision) |
| `run_counting_sc.py` | 31 | Accept .h5 barcode file | Low (feature) |
| `run_counting_sc.py` | 164 | Handle multi-column barcodes | Low (feature) |
| `run_counting_sc.py` | 177 | Output options (MTX, dense) | Low (feature) |
| `count_alleles_sc.py` | 152 | Other features besides 'region' | Low (feature) |
| `parse_gene_data.py` | 78, 99 | Error handling | Medium |
| `filter_variant_data.py` | 100 | Extra validation | Medium |

**Fix**: Either implement or remove TODOs

**Effort**: Varies
**Priority**: Low to Medium

---

## Issue Prioritization for Phase 2

### Phase 2.1: Critical Fixes (Week 1)

1. **C3**: Fix sample parsing bug (1 line × 2 locations) - 🟢 Easy
2. **C1**: Use binary search in bulk counting (10 lines) - 🟢 Easy
3. **C2**: Fix AnnData transpose (needs testing) - 🟡 Medium
4. **H3**: Add proper error handling (4-5 locations) - 🟢 Easy

### Phase 2.2: Code Quality (Week 2)

5. **H1**: Refactor WaspCountFiles duplication - 🟡 Medium
6. **H2**: Remove dead code - 🟢 Easy
7. **H5**: Remove debug print - 🟢 Trivial
8. **H4**: Document/fix read_set behavior - 🟡 Medium

### Phase 2.3: Polish (Week 3)

9. **M1**: Add type hints - 🟡 Medium
10. **M2**: Add docstrings - 🟡 Medium
11. **M3**: Add progress bars - 🟢 Easy
12. **M7**: Add input validation - 🟡 Medium

### Phase 2.4: Cleanup (Week 4)

13. **M5**: Remove commented code - 🟢 Trivial
14. **M6**: Consistent subprocess handling - 🟡 Medium
15. **L1-L13**: Address TODOs - Varies

---

## Testing Requirements

For each fix, need:

1. **Unit tests**: Test individual functions
2. **Integration tests**: Test full workflow
3. **Regression tests**: Compare against baseline outputs

**Baseline Testing**:
```bash
# Before fix:
./scripts/run_baseline.sh  # Establish baseline

# After fix:
./scripts/validate_against_baseline.sh  # Must pass!
```

---

## Metrics

### Current State

- **Total LOC**: 1,430
- **Functions**: ~30
- **Classes**: 3
- **TODO comments**: 13
- **Dead code**: 123 lines (8.6%)
- **Duplicated code**: ~200 lines (14%)

### Target State (Post Phase 2)

- **Dead code**: 0 lines
- **Duplicated code**: <50 lines (<3.5%)
- **Type hints**: 100%
- **Docstrings**: 100%
- **Error handling**: Consistent raise-based
- **Test coverage**: >80%

---

**Document Version**: 1.0
**Next**: Phase 2 refactoring using this inventory
