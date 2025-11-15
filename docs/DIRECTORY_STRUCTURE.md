# WASP2 Directory Structure

**Generated**: 2025-11-15
**Purpose**: Complete mapping of repository structure and file purposes

---

## Repository Overview

```
WASP2-exp/
├── bin/                    # Executables
├── doc/                    # Documentation assets (logos)
├── src/                    # Source code (3 main modules)
│   ├── counting/          # Module 1: Allele counting
│   ├── analysis/          # Module 2: Allelic imbalance analysis
│   └── mapping/           # Module 3: Mapping bias correction
├── test_data/             # Test datasets (on test-data-bundle branch)
├── docs/                  # Project documentation (Phase 1 deliverables)
├── CHANGELOG.md           # Version history
├── LICENSE                # License information
├── README.md              # User documentation
├── ENGINEERING_PLAN.md    # Phase 1 & 2 work plan
└── environment.yml        # Conda environment specification
```

**Code Statistics**:
- **Total Python LOC**: 5,768 lines
- **Modules**: 3 (counting, analysis, mapping)
- **Python files**: 24
- **Entry points**: 4 (bin/WASP2 + 3 module __main__.py files)

---

## `/bin` - Executables

| File | Lines | Purpose | Status |
|------|-------|---------|--------|
| `WASP2` | 38 | Main executable entry point | ⚠️ INCOMPLETE |

**Issues Identified**:
- Contains `TODO MAIN EXECUTABLE` comment
- Only supports "count" and "analysis" commands (not "mapping")
- Hardcoded to import from analysis module for both commands (BUG)
- Does not route to counting module correctly

**Current Implementation**:
```python
# bin/WASP2 (simplified)
if (sys.argv[1] == "count") or (sys.argv[1] == "analysis"):
    sys.path.append(str(root_dir / "src" / "analysis"))  # BUG: Always analysis!
    from run_analysis import parse_cmd, run
```

**Expected Invocation** (per README, using modules directly):
- `python WASP2/src/counting count-variants ...`
- `python WASP2/src/analysis find-imbalance ...`
- `python WASP2/src/mapping make-reads ...`

---

## `/doc` - Documentation Assets

| File | Size | Purpose |
|------|------|---------|
| `wasp2_hex_logo_v1.png` | 60 KB | WASP2 logo (hexagon style) |
| `wasp_hex_logo_v1.png` | 94 KB | Original WASP logo |

---

## `/src/counting` - Counting Module

**Purpose**: Process allele-specific read counts from BAM files overlapping heterozygous SNPs

| File | Lines | Purpose |
|------|-------|---------|
| `__main__.py` | 221 | Typer CLI interface (count-variants, count-variants-sc) |
| `run_counting.py` | 230 | Orchestrator for bulk counting workflow |
| `count_alleles.py` | 123 | Core allele counting logic (BAM + VCF processing) |
| `filter_variant_data.py` | 237 | VCF filtering by sample heterozygosity and regions |
| `parse_gene_data.py` | 213 | GTF/GFF3 annotation parsing |
| `run_counting_sc.py` | 178 | Orchestrator for single-cell counting |
| `count_alleles_sc.py` | 228 | Single-cell allele counting with cell barcodes |
| **Total** | **1,430** | |

### Entry Point Commands

**1. Bulk Counting**: `count-variants`
```bash
python -m src.counting count-variants [BAM] [VCF] {OPTIONS}
```

**Options**:
- `-s/--samples`: Filter for heterozygous SNPs in samples
- `-r/--region`: Filter SNPs overlapping regions (BED/narrowPeak/GTF/GFF3)
- `-o/--out_file`: Output file (default: counts.tsv)
- `-t/--temp_loc`: Keep temporary files for debugging
- `--use_region_names`: Use BED name column instead of coordinates
- RNA-seq specific:
  - `--gene_feature`: Feature type for counting (default: 'exon')
  - `--gene_attribute`: Attribute for ID (default: feature_id/ID)
  - `--gene_parent`: Parent attribute (default: transcript_id/Parent)

**2. Single-Cell Counting**: `count-variants-sc`
```bash
python -m src.counting count-variants-sc [BAM] [VCF] [BARCODES] {OPTIONS}
```

**Output**: H5AD file (AnnData format) with cell x SNP count matrix

### Data Flow

```
BAM + VCF → filter_variant_data → parse_gene_data (if regions)
         → count_alleles → counts.tsv (or .h5ad for sc)
```

---

## `/src/analysis` - Analysis Module

**Purpose**: Statistical analysis of allelic imbalance using beta-binomial models

| File | Lines | Purpose |
|------|-------|---------|
| `__main__.py` | 351 | Typer CLI (find-imbalance, find-imbalance-sc, compare-imbalance) |
| `run_analysis.py` | 204 | Orchestrator for bulk allelic imbalance analysis |
| `as_analysis.py` | 674 | **CORE**: Beta-binomial model, dispersion estimation, optimization |
| `filter_data.py` | 124 | Data filtering and validation |
| `count_alleles.py` | 121 | Analysis-specific counting utilities |
| `run_analysis_sc.py` | 266 | Orchestrator for single-cell per-celltype analysis |
| `as_analysis_sc.py` | 258 | Single-cell statistical methods |
| `run_compare_ai.py` | 77 | Orchestrator for celltype comparison |
| `compare_ai.py` | 516 | Differential allelic imbalance between celltypes |
| `count_alleles_sc.py` | 185 | Single-cell counting support |
| **Total** | **2,776** | |

### Entry Point Commands

**1. Bulk Analysis**: `find-imbalance`
```bash
python -m src.analysis find-imbalance [COUNTS] {OPTIONS}
```

**Options**:
- `-o/--out_file`: Output file (default: ai_results.tsv)
- `--min`: Minimum allele count (default: 10)
- `-p/--pseudocount`: Pseudocount for AI calculation (default: 1)
- `--phased`: Use phased haplotype model (vs unphased)
- `--region_col`: Name of region column (auto-parses by default)
- `--groupby`: Report by parent group instead of feature level (RNA-seq)

**2. Single-Cell Analysis**: `find-imbalance-sc`
```bash
python -m src.analysis find-imbalance-sc [COUNTS.h5ad] [BARCODE_MAP] {OPTIONS}
```

**3. Comparative Analysis**: `compare-imbalance`
```bash
python -m src.analysis compare-imbalance [COUNTS.h5ad] [BARCODE_MAP] {OPTIONS}
```

### Key Statistical Components

**`as_analysis.py`** (674 lines - CRITICAL):
- Beta-binomial distribution modeling
- Dispersion parameter estimation (single and linear models)
- scipy.optimize for maximum likelihood estimation
- Phased vs unphased genotype support
- Z-score filtering for outliers

### Data Flow

```
counts.tsv → filter_data → as_analysis (beta-binomial) → ai_results.tsv
```

For single-cell:
```
counts.h5ad + barcode_map → group by celltype → as_analysis_sc → ai_results_[GROUP].tsv
```

---

## `/src/mapping` - Mapping Module

**Purpose**: Correct for mapping bias via allele swapping and remapping (3-step WASP pipeline)

| File | Lines | Purpose |
|------|-------|---------|
| `__main__.py` | 179 | Typer CLI (make-reads, filter-remapped) |
| `run_mapping.py` | 239 | Orchestrator for both mapping steps |
| `make_remap_reads.py` | 498 | **CORE**: Identify reads, swap alleles, generate FASTQs |
| `filter_remap_reads.py` | 96 | Filter reads that fail to remap to same position |
| `wasp_data_files.py` | 110 | File I/O and metadata management (JSON output) |
| `intersect_variant_data.py` | 305 | VCF/BAM intersection logic |
| `remap_utils.py` | 135 | Utility functions for mapping operations |
| **Total** | **1,562** | |

### Entry Point Commands

**Step 1**: `make-reads` - Create swapped allele reads
```bash
python -m src.mapping make-reads [BAM] [VCF] {OPTIONS}
```

**Options**:
- `--threads`: Thread allocation
- `-s/--samples`: Filter polymorphic SNPs in samples
- `-o/--out_dir`: Output directory for remapped data
- `-t/--temp_loc`: Keep temporary files
- `-j/--out_json`: JSON file with WASP metadata (default: [PREFIX]_wasp_data_files.json)

**Output**:
- `[PREFIX]_swapped_alleles_r1.fq` and `_r2.fq`: Reads with swapped alleles
- `[PREFIX]_to_remap.bam`: Reads to be remapped
- `[PREFIX]_keep.bam`: Reads not overlapping SNPs (kept)
- `[PREFIX]_wasp_data_files.json`: Metadata for step 3

**Step 2**: External remapping (user's choice of aligner)
```bash
# Example with BWA
bwa mem genome.fa swapped_r1.fq swapped_r2.fq | samtools view -b > remapped.bam
samtools sort -o remapped.bam remapped.bam
samtools index remapped.bam
```

**Step 3**: `filter-remapped` - Filter inconsistent reads
```bash
python -m src.mapping filter-remapped [REMAPPED.bam] --json [WASP_JSON] {OPTIONS}
```

**Output**:
- `[PREFIX]_wasp_filt.bam`: Unbiased BAM (only reads that remapped correctly)

### Data Flow

```
Step 1: BAM + VCF → intersect_variant_data → make_remap_reads
        → swapped FASTQs + to_remap.bam + keep.bam + metadata.json

Step 2: [USER REMAPS] swapped FASTQs → remapped.bam

Step 3: remapped.bam + metadata.json → filter_remap_reads
        → wasp_filt.bam (unbiased)
```

---

## `/test_data` - Test Datasets

**Location**: Available on `test-data-bundle` branch only

| File | Size | Description |
|------|------|-------------|
| `wasp2_test_bundle.tar.gz` | 9.3 MB | Compressed test bundle |
| `filter_chr10.vcf` | 12 MB | NA12878 chr10 SNPs (111,454 variants) |
| `NA12878_snps_chr10.bed` | 2.6 MB | Matching intervals |
| `CD4_ATACseq_Day1_merged_filtered.sort.bam` | 7.6 MB | Small ATAC-seq BAM subset |
| `CD4_ATACseq_Day1_merged_filtered.sort.bam.bai` | 2.0 MB | BAM index |
| `as_counts.txt` | 274 bytes | Toy allele-specific counts |
| `README.md` | 345 bytes | Test data description |

**Source**: Aaron Ho's WASP test data for lightweight regression testing

---

## `/docs` - Project Documentation

**Created during Phase 1**: Documentation deliverables

Expected structure:
```
docs/
├── ARCHITECTURE.md                  # High-level system overview
├── DEPENDENCIES.md                  # Dependency analysis
├── DIRECTORY_STRUCTURE.md          # This file
├── modules/
│   ├── COUNTING_MODULE.md          # Counting module deep dive
│   ├── COUNTING_ISSUES.md          # Technical debt inventory
│   ├── ANALYSIS_MODULE.md          # Analysis module deep dive
│   ├── ANALYSIS_ISSUES.md          # Technical debt inventory
│   ├── STATISTICAL_METHODS.md      # Mathematical documentation
│   ├── MAPPING_MODULE.md           # Mapping module deep dive
│   └── MAPPING_ISSUES.md           # Technical debt inventory
├── testing/
│   ├── TEST_STRATEGY.md            # Comprehensive test plan
│   ├── TEST_DATA_REQUIREMENTS.md   # Test data specs
│   ├── CURRENT_COVERAGE.md         # Existing test inventory
│   └── RUNNING_TESTS.md            # Test execution guide
├── security/
│   └── SECURITY_AUDIT.md           # Security findings
├── performance/
│   ├── PERFORMANCE_ANALYSIS.md     # Bottleneck analysis
│   └── BENCHMARKS.md               # Performance metrics
├── CROSS_CUTTING_ISSUES.md         # Common problems
├── REFACTORING_OPPORTUNITIES.md    # Consolidation opportunities
├── PHASE1_FINDINGS.md              # Executive summary
├── PHASE2_ROADMAP.md               # Prioritized work plan
└── CONFIGURATION.md                # Configuration guide
```

---

## Configuration Files

### `environment.yml`

**Purpose**: Conda environment specification

**Python Version**: 3.9.*

**Dependencies**:
- Data processing: `numpy`, `pandas`, `polars`
- Bioinformatics: `pysam`, `pybedtools`, `bedtools`, `anndata`
- Statistics: `scipy`
- CLI: `typer`

**Channels**:
1. `bioconda` - Bioinformatics packages
2. `conda-forge` - Community packages
3. `defaults` - Anaconda defaults

---

## Entry Points Summary

### Module Invocation (Current Working Method)

```bash
# Counting module
python -m src.counting count-variants [BAM] [VCF] {OPTIONS}
python -m src.counting count-variants-sc [BAM] [VCF] [BARCODES] {OPTIONS}

# Analysis module
python -m src.analysis find-imbalance [COUNTS] {OPTIONS}
python -m src.analysis find-imbalance-sc [COUNTS] [BARCODE_MAP] {OPTIONS}
python -m src.analysis compare-imbalance [COUNTS] [BARCODE_MAP] {OPTIONS}

# Mapping module
python -m src.mapping make-reads [BAM] [VCF] {OPTIONS}
python -m src.mapping filter-remapped [REMAPPED_BAM] --json [JSON] {OPTIONS}
```

### Main Executable (bin/WASP2) - BROKEN

**Status**: ⚠️ Incomplete implementation, has bugs

**TODO**: Fix routing to proper modules or document as deprecated

---

## File Type Breakdown

| Type | Count | Purpose |
|------|-------|---------|
| Python modules (`*.py`) | 24 | Source code |
| CLI entry points (`__main__.py`) | 3 | Module invocation |
| Configuration files | 1 | `environment.yml` |
| Documentation (`*.md`) | 3 | README, CHANGELOG, ENGINEERING_PLAN |
| Images | 2 | Logos |
| Shell scripts | 1 | `bin/WASP2` |

---

## Code Distribution by Module

| Module | Files | Lines | % of Code | Complexity |
|--------|-------|-------|-----------|------------|
| **Analysis** | 10 | 2,776 | 48.1% | HIGH (statistical) |
| **Mapping** | 7 | 1,562 | 27.1% | MEDIUM |
| **Counting** | 7 | 1,430 | 24.8% | MEDIUM |
| **TOTAL** | 24 | 5,768 | 100% | |

**Largest Files** (by LOC):
1. `analysis/as_analysis.py` - 674 lines (beta-binomial models)
2. `analysis/compare_ai.py` - 516 lines (differential AI)
3. `mapping/make_remap_reads.py` - 498 lines (allele swapping)
4. `analysis/__main__.py` - 351 lines (CLI definitions)
5. `mapping/intersect_variant_data.py` - 305 lines (VCF/BAM intersection)

---

## Critical Observations

### ⚠️ Issues Found

1. **bin/WASP2 is broken**:
   - Contains TODO comment
   - Only supports 2 of 3 modules
   - Hardcoded bug in module routing
   - Doesn't match README documentation

2. **No tests directory**:
   - No unit tests found
   - No integration tests
   - Test data exists but no test harness

3. **No shared utilities**:
   - Code duplication likely between modules
   - No common library for shared functions

4. **Missing documentation**:
   - No API documentation
   - No developer guide
   - No architectural diagrams

### ✅ Strengths

1. **Modular design**: Clear separation of concerns
2. **CLI framework**: Typer provides good UX
3. **Single-cell support**: Well-integrated sc workflows
4. **Test data available**: Lightweight regression dataset exists

---

## Next Steps

As per ENGINEERING_PLAN.md:
- ✅ **Task 1.1.1**: Complete (this document)
- ⏭️ **Task 1.1.2**: Analyze external dependencies → `docs/DEPENDENCIES.md`
- ⏭️ **Task 1.1.3**: Create architecture diagram → `docs/ARCHITECTURE.md`
