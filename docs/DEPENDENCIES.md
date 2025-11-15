# WASP2 Dependency Analysis

**Generated**: 2025-11-15
**Purpose**: Comprehensive analysis of external dependencies, their usage, versions, and potential issues

---

## Overview

WASP2 has **9 core dependencies** specified in `environment.yml`:

| Dependency | Version Constraint | Type | Usage Count |
|------------|-------------------|------|-------------|
| `python` | 3.9.* | Runtime | - |
| `numpy` | Latest | Data processing | 8 files |
| `pandas` | Latest | Data processing | 9 files |
| `polars` | Latest | Data processing | 8 files |
| `scipy` | Latest | Statistics | 8 files |
| `pysam` | Latest | Bioinformatics | 6 files |
| `pybedtools` | Latest | Bioinformatics | 1 file |
| `bedtools` | Latest | System tool | Indirect (via pybedtools) |
| `typer` | Latest | CLI | 3 files |
| `anndata` | Latest | Single-cell | 6 files |

**Total Dependencies**: 9 direct + 1 indirect (bedtools)

---

## Python Version

**Specified**: `python=3.9.*`

**Status**: ✅ Good choice
- Python 3.9 is stable and widely supported
- Released October 2020, maintenance through October 2025
- Good balance of stability and modern features

**Considerations**:
- Python 3.9 reaches end-of-life in October 2025
- Should plan migration to 3.10+ in next 6-12 months
- Type hints improvements in 3.10+ would benefit codebase

**Recommendation**:
- Keep 3.9 for now (stable)
- Add 3.10 and 3.11 compatibility testing
- Plan migration to 3.10+ by Q3 2025

---

## Conda Channels

**Order** (priority):
1. `bioconda` - Bioinformatics-specific packages
2. `conda-forge` - Community-driven packages
3. `defaults` - Anaconda official packages

**Status**: ✅ Correct channel priority
- `bioconda` first ensures correct bioinformatics tool versions
- `conda-forge` provides newer/better-maintained packages
- `defaults` as fallback

**Considerations**:
- Channel order matters for package resolution
- Some packages may have conflicts between channels
- `bioconda` depends on `conda-forge`, so both are needed

---

## Dependency Deep Dive

### 1. NumPy - Numerical Computing

**Version**: Latest (from conda-forge)

**Purpose**: Array operations, numerical computations

**Usage Locations**:
- `analysis/as_analysis.py` - Statistical calculations, array operations
- `analysis/as_analysis_sc.py` - Single-cell numerical operations
- `analysis/compare_ai.py` - Comparison calculations
- `analysis/run_analysis_sc.py` - Data aggregation
- `counting/filter_variant_data.py` - Variant filtering
- `counting/count_alleles_sc.py` - Sparse matrix operations
- `mapping/intersect_variant_data.py` - Position calculations

**Key Functions Used**:
- Array operations: `np.array`, `np.zeros`, `np.ones`
- Mathematical functions: `np.log`, `np.exp`, `np.sum`, `np.mean`
- Statistical functions: `np.median`, `np.std`
- Indexing and filtering

**Risk**: 🟢 LOW
- Very stable, mature library
- Core dependency for scientific Python
- Breaking changes are rare

**Version Constraint**: NONE specified
- ⚠️ **Issue**: No version pinning could lead to compatibility issues
- **Recommendation**: Pin to compatible range (e.g., `numpy>=1.20,<2.0`)

---

### 2. Pandas - Data Manipulation

**Version**: Latest (from conda-forge)

**Purpose**: DataFrame operations, CSV I/O, data aggregation

**Usage Locations**:
- `analysis/as_analysis.py` - Data processing, groupby operations
- `analysis/as_analysis_sc.py` - Single-cell data frames
- `analysis/compare_ai.py` - Comparative analysis data
- `analysis/run_analysis.py` - Input/output handling
- `analysis/run_analysis_sc.py` - Cell type grouping
- `analysis/run_compare_ai.py` - Comparison orchestration
- `counting/count_alleles_sc.py` - Cell barcode mapping
- `analysis/count_alleles_sc.py` - (duplicate module name)

**Key Functions Used**:
- DataFrame creation and manipulation
- `pd.read_csv`, `pd.DataFrame`
- GroupBy operations: `df.groupby().agg()`
- Merging and joining dataframes
- Statistical aggregations

**Risk**: 🟡 MEDIUM
- Pandas 2.0+ introduces breaking changes
- Performance differences between versions
- Memory usage patterns vary

**Version Constraint**: NONE specified
- ⚠️ **Issue**: Pandas 2.0+ may break existing code
- **Recommendation**: Pin to `pandas>=1.3,<2.1` and test compatibility

**Pandas vs Polars Redundancy**:
- ⚠️ **Both pandas AND polars are used**
- This is intentional but adds complexity
- Polars for performance-critical operations
- Pandas for compatibility and convenience

---

### 3. Polars - High-Performance DataFrames

**Version**: Latest (from conda-forge)

**Purpose**: High-performance DataFrame operations, large file handling

**Usage Locations**:
- `counting/count_alleles.py` - Fast allele counting
- `counting/count_alleles_sc.py` - Single-cell counting
- `counting/filter_variant_data.py` - VCF filtering
- `mapping/intersect_variant_data.py` - VCF/BAM intersection

**Key Functions Used**:
- `pl.read_csv` - Fast CSV reading
- `pl.DataFrame` - DataFrame creation
- Lazy evaluation for large files
- Fast filtering and aggregation
- Schema validation

**Risk**: 🟡 MEDIUM
- Relatively new library (1.0 released 2023)
- API still evolving rapidly
- Breaking changes more frequent than Pandas

**Version Constraint**: NONE specified
- ⚠️ **Critical Issue**: Polars API changes frequently
- **Recommendation**: Pin to specific minor version (e.g., `polars>=0.19,<0.20`)

**Why Both Pandas and Polars?**:
- **Polars**: Performance-critical operations (VCF parsing, large BAM files)
- **Pandas**: Statistical analysis, compatibility with scipy/scikit-learn
- **Trade-off**: Complexity vs performance

**Recommendation**: Document when to use each, consider consolidation in Phase 2

---

### 4. SciPy - Scientific Computing

**Version**: Latest (from conda-forge)

**Purpose**: Statistical distributions, optimization, scientific functions

**Usage Locations** (8 files - CRITICAL DEPENDENCY):
- `analysis/as_analysis.py` - **CORE STATISTICAL LOGIC**
- `analysis/as_analysis_sc.py` - Single-cell statistics
- `analysis/compare_ai.py` - Differential analysis
- `counting/count_alleles_sc.py` - Sparse matrices

**Key Modules Used**:

#### `scipy.stats`:
- `betabinom` - Beta-binomial distribution (core AI model)
- `chi2` - Chi-squared distribution (p-values)
- `binom` - Binomial distribution
- `zscore` - Z-score normalization
- `rankdata` - Rank calculations
- `false_discovery_control` - FDR correction

#### `scipy.optimize`:
- `minimize_scalar` - Dispersion parameter optimization
- `minimize` - Multi-parameter optimization

#### `scipy.special`:
- `expit` - Logistic sigmoid function

#### `scipy.sparse`:
- `csr_matrix` - Compressed sparse row matrices (single-cell data)

**Risk**: 🟢 LOW to 🟡 MEDIUM
- Stable API for core functions
- Statistical functions well-tested
- Optimization algorithms robust

**Version Constraint**: NONE specified
- ⚠️ **Issue**: Statistical API changes could affect results
- **Recommendation**: Pin to `scipy>=1.7,<1.12` for stability

**Critical for Correctness**:
- Beta-binomial implementation is core to scientific validity
- Any version changes should be tested for numerical consistency
- Optimization results may vary between versions

---

### 5. pysam - SAM/BAM/VCF File Handling

**Version**: Latest (from bioconda)

**Purpose**: Read and write SAM/BAM/CRAM/VCF files

**Usage Locations**:
- `analysis/count_alleles.py` - Read BAM files
- `counting/count_alleles.py` - Read BAM files
- `counting/count_alleles_sc.py` - Read BAM with cell barcodes
- `mapping/intersect_variant_data.py` - VCF and BAM intersection
- `mapping/wasp_data_files.py` - VCF reading, BAM I/O
- `mapping/make_remap_reads.py` - (likely, not in sample)

**Key Classes Used**:
- `AlignmentFile` - BAM file reading/writing
- `VariantFile` - VCF file reading

**Key Operations**:
- Reading alignments from BAM
- Fetching reads by region
- Reading VCF variant records
- Writing BAM files

**Risk**: 🟢 LOW
- Mature bioinformatics standard library
- Wrapper around htslib (C library)
- Well-tested, stable API

**Version Constraint**: NONE specified
- ⚠️ **Issue**: htslib updates could change behavior
- **Recommendation**: Pin to `pysam>=0.19,<0.22` for stability

**Performance Notes**:
- BAM reading is I/O bound
- Indexing (BAI files) is critical for performance
- Random access vs sequential reading patterns

---

### 6. pybedtools - Genomic Interval Operations

**Version**: Latest (from bioconda)

**Purpose**: Wrapper around bedtools for genomic interval operations

**Usage Locations**:
- `mapping/intersect_variant_data.py` - BedTool import (1 file only!)

**Key Functions Used**:
- `BedTool` - Genomic interval representation
- Intersection operations
- Coordinate manipulation

**Risk**: 🟡 MEDIUM
- Wrapper around bedtools binary (system dependency)
- Requires bedtools to be installed
- Less frequently maintained than pysam

**Version Constraint**: NONE specified
- **Recommendation**: Pin to `pybedtools>=0.9,<1.0`

**Dependency Chain**:
- `pybedtools` (Python) → `bedtools` (system binary)
- Both must be installed via conda

**⚠️ Limited Usage**:
- Only used in 1 file (intersect_variant_data.py)
- Consider replacing with pure pysam operations
- Or polars interval operations
- **Phase 2 opportunity**: Remove dependency if possible

---

### 7. bedtools - Genomic Toolset (System Binary)

**Version**: Latest (from bioconda)

**Purpose**: System dependency for pybedtools

**Risk**: 🟢 LOW
- Industry standard tool
- Very stable

**Note**: Only needed because of pybedtools
- If pybedtools is removed, this can be removed too

---

### 8. Typer - CLI Framework

**Version**: Latest (from conda-forge)

**Purpose**: Command-line interface framework

**Usage Locations**:
- `counting/__main__.py` - Counting CLI commands
- `analysis/__main__.py` - Analysis CLI commands
- `mapping/__main__.py` - Mapping CLI commands

**Key Features Used**:
- `@app.command()` decorators
- Type annotations for arguments
- `Annotated` types for parameter descriptions
- Automatic help generation

**Dependencies**:
- `typing_extensions` - Extended type hints (imported separately)

**Risk**: 🟢 LOW
- Well-maintained by tiangolo (FastAPI creator)
- Stable API
- Good documentation

**Version Constraint**: NONE specified
- **Recommendation**: Pin to `typer>=0.9,<1.0`

**Quality Notes**:
- All 3 modules have TODO comments: "GOTTA TEST THIS"
- CLI interfaces may not be fully tested
- Phase 1 should document CLI contract
- Phase 2 should add CLI tests

---

### 9. anndata - Annotated Data Matrices

**Version**: Latest (from conda-forge)

**Purpose**: Single-cell data format (h5ad files)

**Usage Locations**:
- `analysis/as_analysis_sc.py` - Single-cell analysis
- `analysis/run_compare_ai.py` - Celltype comparison
- `analysis/run_analysis_sc.py` - SC orchestration
- `counting/count_alleles_sc.py` - SC count output

**Key Classes Used**:
- `AnnData` - Annotated data matrix
- Cell × SNP matrices
- Obs/var annotations
- H5AD file I/O

**Risk**: 🟡 MEDIUM
- Evolving API for single-cell standards
- Large files can have memory issues
- Compatibility with scanpy ecosystem important

**Version Constraint**: NONE specified
- **Recommendation**: Pin to `anndata>=0.8,<0.10`

**Single-Cell Workflow**:
```
Counting → AnnData (cells × SNPs) → H5AD file
↓
Analysis → Per-celltype aggregation → Results
```

---

## Missing Dependencies

### Detected in Code but Not in environment.yml

1. **`typing_extensions`** - Used in 3 `__main__.py` files
   - Provides `Annotated` type for Typer
   - Should be added to environment.yml
   - Or rely on Typer to pull it in

**Recommendation**: Add `typing_extensions` explicitly

---

## Security Vulnerabilities

**⚠️ No Version Pinning = Security Risk**

**Issues**:
1. **No version constraints** on any package
2. `conda install` will pull latest versions
3. Security vulnerabilities in older versions not protected against
4. Supply chain attacks possible with auto-updates

**Recommendations**:
1. Pin all major versions (e.g., `numpy>=1.20,<2.0`)
2. Use `conda list --export > requirements.txt` to freeze exact versions
3. Regularly update and test with newer versions
4. Subscribe to security advisories for key packages

**High-Risk Packages** (due to popularity):
- numpy, pandas, scipy - Widely targeted
- pysam - Handles untrusted binary files (BAM/VCF)

---

## Performance Considerations

### Memory Usage

**High Memory Dependencies**:
1. **pandas** - DataFrames can be memory-intensive
   - Entire DataFrame loaded into memory
   - Recommendation: Use chunking for large files

2. **anndata** - H5AD files can be very large (single-cell)
   - Backed mode available for large datasets
   - Consider streaming operations

3. **pysam** - BAM files read into memory per region
   - Use indexed fetching to reduce memory
   - Process chromosomes sequentially

### CPU Performance

**Optimized Libraries**:
- **polars** - Rust-based, very fast
- **numpy** - Compiled C/Fortran code
- **pysam** - C library (htslib)

**Optimization Opportunities**:
- Parallelize scipy.optimize calls (currently sequential)
- Use polars lazy evaluation more extensively
- Batch BAM file operations

---

## Dependency Redundancies

### Pandas vs Polars

**Current Usage**:
- Both used throughout codebase
- Polars for: VCF parsing, large file operations
- Pandas for: Statistical operations, groupby, compatibility

**Issues**:
- Increases package size
- Requires converting between formats
- Developers need to know both APIs

**Recommendations**:
1. **Keep both** for now (different strengths)
2. Document usage guidelines:
   - Use Polars for: File I/O, filtering, large datasets
   - Use Pandas for: Statistics, scipy integration, small datasets
3. **Phase 2**: Evaluate if consolidation is possible

---

## Dependency Graph

```
WASP2
├── Python 3.9
├── Data Processing
│   ├── numpy (arrays, math)
│   ├── pandas (DataFrames, stats)
│   └── polars (fast DataFrames)
├── Statistics
│   └── scipy (betabinom, optimize)
├── Bioinformatics
│   ├── pysam (BAM/VCF I/O)
│   ├── pybedtools (intervals)
│   └── bedtools (system binary)
├── Single-Cell
│   └── anndata (h5ad format)
└── CLI
    ├── typer (CLI framework)
    └── typing_extensions (type hints)
```

---

## Version Compatibility Matrix

| Package | Current | Tested | Min Recommended | Max Recommended |
|---------|---------|--------|-----------------|-----------------|
| python | 3.9.* | 3.9 | 3.9 | 3.11 |
| numpy | latest | ? | 1.20 | <2.0 |
| pandas | latest | ? | 1.3 | <2.1 |
| polars | latest | ? | 0.18 | <1.0 |
| scipy | latest | ? | 1.7 | <1.12 |
| pysam | latest | ? | 0.19 | <0.22 |
| pybedtools | latest | ? | 0.9 | <1.0 |
| bedtools | latest | ? | 2.30 | <3.0 |
| typer | latest | ? | 0.9 | <1.0 |
| anndata | latest | ? | 0.8 | <0.10 |

**⚠️ Testing Status**: Unknown - no test suite found

---

## Recommendations Summary

### Critical (Phase 2, Task 2.1)

1. **Add version constraints to environment.yml**:
```yaml
dependencies:
  - python=3.9.*
  - numpy>=1.20,<2.0
  - pandas>=1.3,<2.1
  - polars>=0.18,<1.0
  - scipy>=1.7,<1.12
  - pysam>=0.19,<0.22
  - pybedtools>=0.9,<1.0
  - bedtools>=2.30,<3.0
  - typer>=0.9,<1.0
  - anndata>=0.8,<0.10
  - typing_extensions
```

2. **Create requirements.txt with frozen versions**:
```bash
conda list --export > requirements-frozen.txt
```

3. **Security audit**:
   - Check for known vulnerabilities
   - Set up Dependabot or Renovate
   - Regular dependency updates

### High Priority (Phase 2)

4. **Evaluate pybedtools removal**:
   - Only used in 1 file
   - Can likely replace with pysam or polars
   - Reduces dependency count

5. **Document pandas vs polars usage guidelines**:
   - When to use each
   - Conversion patterns
   - Performance considerations

6. **Add compatibility tests**:
   - Test with Python 3.10, 3.11
   - Test with pandas 2.0+
   - Test with latest polars

### Medium Priority

7. **Optimize scipy usage**:
   - Profile optimization calls
   - Consider parallelization
   - Validate numerical stability

8. **Monitor anndata evolution**:
   - Single-cell standards changing
   - Keep compatible with scanpy ecosystem

---

## Phase 1 Action Items

- ✅ **Task 1.1.2**: Complete (this document)
- ⏭️ **Next**: Create architecture diagram (Task 1.1.3)

**No code changes in Phase 1** - documentation only!

Version pinning and dependency updates scheduled for **Phase 2, Task 2.1**.
