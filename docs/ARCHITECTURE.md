# WASP2 Architecture Documentation

**Generated**: 2025-11-15
**Purpose**: High-level system architecture, design patterns, and data flow

---

## Table of Contents

1. [System Overview](#system-overview)
2. [Architectural Principles](#architectural-principles)
3. [Module Architecture](#module-architecture)
4. [Data Flow](#data-flow)
5. [Design Patterns](#design-patterns)
6. [Technology Stack](#technology-stack)
7. [Integration Points](#integration-points)
8. [Architectural Issues](#architectural-issues)

---

## System Overview

WASP2 is a **modular bioinformatics pipeline** for detecting and analyzing allele-specific patterns in genomics data.

### Core Capabilities

1. **Allele Counting** - Quantify read counts for each allele at heterozygous SNPs
2. **Imbalance Analysis** - Statistically measure allelic imbalance
3. **Mapping Bias Correction** - Remove reads with mapping bias

### Target Data Types

- Bulk RNA-seq and ATAC-seq
- Single-cell RNA-seq (scRNA-seq)
- Single-cell ATAC-seq (scATAC-seq)

### High-Level Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                         WASP2 SYSTEM                            │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐         │
│  │   COUNTING   │  │   ANALYSIS   │  │   MAPPING    │         │
│  │    MODULE    │  │    MODULE    │  │    MODULE    │         │
│  │              │  │              │  │              │         │
│  │ BAM + VCF    │  │   Counts     │  │ BAM + VCF    │         │
│  │   ↓          │  │     ↓        │  │   ↓          │         │
│  │ Allele       │  │ Beta-        │  │ Allele       │         │
│  │ Counts       │  │ Binomial     │  │ Swapping     │         │
│  │              │  │ Model        │  │   ↓          │         │
│  │              │  │   ↓          │  │ Remap        │         │
│  │              │  │ AI Results   │  │   ↓          │         │
│  │              │  │              │  │ Filter       │         │
│  └──────────────┘  └──────────────┘  └──────────────┘         │
│                                                                 │
│  Each module: Independent | Typer CLI | Python Package          │
└─────────────────────────────────────────────────────────────────┘
```

---

## Architectural Principles

### 1. **Modularity**

**Design**: Three independent modules with clear boundaries

```
src/
├── counting/    # Input: BAM + VCF → Output: Counts
├── analysis/    # Input: Counts → Output: AI Statistics
└── mapping/     # Input: BAM + VCF → Output: Filtered BAM
```

**Benefits**:
- Modules can be used independently
- Clear separation of concerns
- Easier testing and maintenance

**Trade-offs**:
- Some code duplication between modules
- Data format dependencies between modules

### 2. **CLI-First Design**

**Framework**: Typer (Python CLI framework)

**Pattern**: Each module = standalone CLI application

```python
# Example: src/counting/__main__.py
app = typer.Typer()

@app.command()
def count_variants(bam: str, vcf: str, ...):
    run_count_variants(bam, vcf, ...)
```

**Invocation**:
```bash
python -m src.counting count-variants [args]
python -m src.analysis find-imbalance [args]
python -m src.mapping make-reads [args]
```

**Benefits**:
- Self-documenting with `--help`
- Type-safe argument parsing
- Automatic validation

**Issues**:
- bin/WASP2 wrapper is incomplete/broken
- No Python API for programmatic use

### 3. **Data Pipeline Architecture**

**Pattern**: Multi-stage pipelines with intermediate files

```
Stage 1: Input → Processing → Intermediate Files
Stage 2: Intermediate Files → Analysis → Output
```

**Rationale**:
- Allows pipeline restarts at any stage
- Facilitates debugging (inspect intermediate files)
- Enables manual intervention between stages

**Trade-offs**:
- More I/O overhead
- Requires temp file management
- Disk space requirements

---

## Module Architecture

### Counting Module

**Purpose**: Extract allele-specific read counts from sequencing data

#### Component Structure

```
counting/
│
├── __main__.py          ← CLI Entry Point (Typer commands)
│   ├── count_variants()       → Bulk counting
│   └── count_variants_sc()    → Single-cell counting
│
├── run_counting.py      ← Orchestrator (Bulk)
│   └── run_count_variants()
│       ├── Validate inputs
│       ├── Filter variants
│       ├── Parse gene data (if regions provided)
│       ├── Count alleles
│       └── Write output
│
├── run_counting_sc.py   ← Orchestrator (Single-Cell)
│   └── run_count_variants_sc()
│       └── Similar flow + cell barcode handling
│
├── filter_variant_data.py
│   └── filter_vcf_samples()
│       └── Filter VCF for heterozygous SNPs in samples
│
├── parse_gene_data.py
│   └── parse_genes()
│       └── Parse GTF/GFF3 annotations
│
├── count_alleles.py     ← Core Logic (Bulk)
│   └── count_alleles()
│       ├── Read BAM with pysam
│       ├── Intersect reads with SNP positions
│       ├── Count ref/alt alleles
│       └── Return polars DataFrame
│
└── count_alleles_sc.py  ← Core Logic (Single-Cell)
    └── count_alleles_sc()
        ├── Extract cell barcodes
        ├── Build cell × SNP matrix
        └── Return AnnData object
```

#### Data Flow

```
INPUT: BAM + VCF + [optional: regions, samples]
  │
  ├─→ filter_variant_data.py
  │     └─→ Filtered VCF (heterozygous SNPs only)
  │
  ├─→ parse_gene_data.py (if regions)
  │     └─→ Region definitions (BED format)
  │
  └─→ count_alleles.py
        ├─ Read BAM alignments (pysam)
        ├─ For each read overlapping SNP:
        │    ├─ Extract base at SNP position
        │    └─ Increment ref/alt count
        └─→ OUTPUT: counts.tsv or counts.h5ad

Format (TSV):
chrom  pos    ref  alt  region   ref_count  alt_count  other_count
chr1   12345  A    G    peak123  15         8          0
```

#### Key Algorithms

**Binary Search for Position Lookup**:
- SNPs sorted by position
- `bisect_left()` for O(log n) position finding
- Critical for performance with many SNPs

**Cell Barcode Handling** (Single-Cell):
- Extract CB tag from BAM
- Build sparse matrix (scipy.sparse.csr_matrix)
- Efficient for sparse single-cell data

---

### Analysis Module

**Purpose**: Statistical analysis of allelic imbalance

#### Component Structure

```
analysis/
│
├── __main__.py          ← CLI Entry Point
│   ├── find_imbalance()        → Bulk analysis
│   ├── find_imbalance_sc()     → Single-cell per-celltype
│   └── compare_imbalance()     → Differential AI
│
├── run_analysis.py      ← Orchestrator (Bulk)
│   └── run_ai_analysis()
│       ├── Load count data
│       ├── Filter by min counts
│       ├── Group by region/gene
│       └── Call as_analysis.get_imbalance()
│
├── as_analysis.py       ← CORE STATISTICAL ENGINE ⭐
│   ├── get_imbalance()
│   │   ├── Aggregate counts by region
│   │   ├── Estimate dispersion (beta-binomial)
│   │   ├── Calculate imbalance metrics
│   │   └── Compute p-values and FDR
│   │
│   ├── opt_prob() / opt_unphased_dp() / opt_phased_new()
│   │   └── Optimization functions for dispersion
│   │
│   └── bh_correction()
│       └── Benjamini-Hochberg FDR correction
│
├── as_analysis_sc.py    ← Single-Cell Statistics
│   ├── get_imbalance_sc()
│   └── adata_count_qc()
│       └── Z-score filtering for QC
│
├── compare_ai.py        ← Differential Imbalance
│   └── get_compared_imbalance()
│       └── Compare AI between celltypes
│
├── run_analysis_sc.py   ← SC Orchestrator
├── run_compare_ai.py    ← Comparison Orchestrator
├── filter_data.py       ← Data validation
└── count_alleles*.py    ← Utilities (some duplication with counting/)
```

#### Core Statistical Model

**Beta-Binomial Distribution**:

```
Read counts ~ BetaBinomial(n, α, β)

Where:
- n = total reads (ref + alt)
- α, β = shape parameters of beta prior
- Dispersion ρ = 1 / (α + β + 1)
```

**Null Hypothesis**: No allelic imbalance (50/50 ref/alt)

**Alternative Models**:

1. **Unphased** (default):
   - Equal probability for either haplotype
   - H0: p = 0.5 (no imbalance)

2. **Phased**:
   - Known maternal/paternal assignment
   - Can test directional imbalance

#### Optimization Pipeline

```
┌────────────────────────────────────────────────┐
│  1. Aggregate Counts by Region                 │
│     ref_sum, alt_sum, total                    │
└─────────────────┬──────────────────────────────┘
                  │
┌─────────────────▼──────────────────────────────┐
│  2. Estimate Dispersion Parameter (ρ)          │
│     Method: Maximum Likelihood Estimation      │
│     Optimizer: scipy.optimize.minimize_scalar  │
│     Search: [0, 1] with bounds                 │
└─────────────────┬──────────────────────────────┘
                  │
┌─────────────────▼──────────────────────────────┐
│  3. Calculate Allelic Imbalance Metrics        │
│     - Log2(alt/ref) fold change                │
│     - Beta-binomial p-value                    │
│     - Effect size estimates                    │
└─────────────────┬──────────────────────────────┘
                  │
┌─────────────────▼──────────────────────────────┐
│  4. Multiple Testing Correction                │
│     Method: Benjamini-Hochberg FDR             │
└─────────────────┬──────────────────────────────┘
                  │
                OUTPUT: AI statistics per region
```

#### Data Flow

```
INPUT: counts.tsv
  │
  ├─→ filter_data.py
  │     └─→ Filter by min_count threshold
  │
  ├─→ Group by region/gene
  │
  ├─→ as_analysis.py
  │     ├─ For each region:
  │     │   ├─ Aggregate ref + alt counts
  │     │   ├─ Optimize dispersion parameter
  │     │   ├─ Calculate log2(alt/ref)
  │     │   ├─ Compute p-value (betabinom test)
  │     │   └─ Estimate effect size
  │     │
  │     └─ FDR correction across all regions
  │
  └─→ OUTPUT: ai_results.tsv

Format:
region   ref_total  alt_total  log2_fc  pvalue   fdr      dispersion
peak123  150        80         -0.91    0.0023   0.045    0.012
```

---

### Mapping Module

**Purpose**: Remove allele-specific mapping bias

#### Component Structure

```
mapping/
│
├── __main__.py          ← CLI Entry Point
│   ├── make_reads()         → Step 1: Identify + swap
│   └── filter_remapped()    → Step 3: Filter remapped
│
├── run_mapping.py       ← Orchestrators
│   ├── run_make_remap_reads()   → Step 1 orchestration
│   └── run_wasp_filt()          → Step 3 orchestration
│
├── make_remap_reads.py  ← CORE: Allele Swapping ⭐
│   └── make_remap_reads()
│       ├── Intersect BAM with VCF
│       ├── Identify reads overlapping SNPs
│       ├── Swap alleles at SNP positions
│       ├── Write FASTQ with swapped reads
│       └── Track metadata for filtering
│
├── filter_remap_reads.py ← Step 3: Filter
│   └── filter_remap_reads()
│       ├── Compare original vs remapped positions
│       ├── Keep only identically mapped reads
│       └── Write filtered BAM
│
├── intersect_variant_data.py  ← VCF/BAM Intersection
│   └── intersect_bam_vcf()
│       └── Find reads overlapping SNP positions
│
├── wasp_data_files.py   ← File Management
│   └── WaspDataFiles (class)
│       ├── Manage temp files
│       ├── Track file paths
│       └── Write/read metadata JSON
│
└── remap_utils.py       ← Utilities
    └── Helper functions
```

#### Three-Step WASP Pipeline

```
┌──────────────────────────────────────────────────────────────┐
│  STEP 1: make-reads (WASP Tool)                              │
├──────────────────────────────────────────────────────────────┤
│                                                              │
│  INPUT: original.bam + variants.vcf                          │
│     │                                                        │
│     ├─→ intersect_variant_data                              │
│     │     └─ Find reads overlapping SNPs                    │
│     │                                                        │
│     ├─→ make_remap_reads                                    │
│     │     ├─ For each read with SNP:                        │
│     │     │   ├─ Extract read sequence                      │
│     │     │   ├─ Swap allele at SNP position(s)            │
│     │     │   └─ Write to FASTQ                            │
│     │     │                                                 │
│     │     └─ Split reads into 3 categories:                │
│     │         ├─ to_remap.bam (reads with SNPs)           │
│     │         ├─ keep.bam (reads without SNPs)            │
│     │         └─ swapped.fq (FASTQs for remapping)        │
│     │                                                        │
│     └─→ wasp_data_files                                     │
│           └─ Write metadata JSON                            │
│                                                              │
│  OUTPUT: swapped_r1.fq, swapped_r2.fq, metadata.json        │
└──────────────────────────────────────────────────────────────┘
                           │
                           ▼
┌──────────────────────────────────────────────────────────────┐
│  STEP 2: Remap (USER'S ALIGNER - BWA/STAR/etc)              │
├──────────────────────────────────────────────────────────────┤
│                                                              │
│  Example:                                                    │
│  $ bwa mem genome.fa swapped_r1.fq swapped_r2.fq \          │
│      | samtools view -b > remapped.bam                       │
│  $ samtools sort -o remapped.bam remapped.bam               │
│  $ samtools index remapped.bam                              │
│                                                              │
│  OUTPUT: remapped.bam                                        │
└──────────────────────────────────────────────────────────────┘
                           │
                           ▼
┌──────────────────────────────────────────────────────────────┐
│  STEP 3: filter-remapped (WASP Tool)                         │
├──────────────────────────────────────────────────────────────┤
│                                                              │
│  INPUT: remapped.bam + metadata.json                         │
│     │                                                        │
│     ├─→ filter_remap_reads                                  │
│     │     ├─ For each read in remapped.bam:                │
│     │     │   ├─ Find corresponding read in to_remap.bam   │
│     │     │   ├─ Compare genomic positions                 │
│     │     │   └─ KEEP if: same chr, pos, strand, CIGAR    │
│     │     │       DISCARD if: different mapping            │
│     │     │                                                 │
│     │     └─ Merge kept reads with keep.bam               │
│     │                                                        │
│     └─→ OUTPUT: wasp_filt.bam (unbiased BAM)               │
│                                                              │
└──────────────────────────────────────────────────────────────┘
```

#### Allele Swapping Algorithm

```python
# Pseudocode
for read in bam:
    if read overlaps SNPs:
        for snp in overlapping_snps:
            # Get read base at SNP position
            read_base = read.seq[snp.read_offset]

            # Determine swap
            if read_base == snp.ref:
                swapped_base = snp.alt
            elif read_base == snp.alt:
                swapped_base = snp.ref
            else:
                continue  # Doesn't match ref or alt

            # Create swapped read
            new_seq = swap_base(read.seq, snp.read_offset, swapped_base)
            new_qual = read.qual  # Quality scores unchanged

        # Write swapped read to FASTQ
        write_fastq(read.name, new_seq, new_qual)
```

#### Metadata Management

**JSON Structure** (`wasp_data_files.json`):
```json
{
  "bam_prefix": "sample123",
  "to_remap_bam": "sample123_to_remap.bam",
  "keep_bam": "sample123_keep.bam",
  "swapped_r1": "sample123_swapped_r1.fq",
  "swapped_r2": "sample123_swapped_r2.fq",
  "read_count": 1234567,
  "snp_overlaps": 45678
}
```

**Purpose**: Links files across 3-step pipeline

---

## Data Flow

### End-to-End Workflow Example: ATAC-seq Analysis

```
┌────────────────────────────────────────────────────────────────┐
│  1. MAP READS (External: BWA/Bowtie2)                          │
│     FASTQ → Aligned BAM                                         │
└─────────────────────┬──────────────────────────────────────────┘
                      │
┌─────────────────────▼──────────────────────────────────────────┐
│  2. REMOVE MAPPING BIAS (WASP Mapping Module)                  │
│     python -m src.mapping make-reads bam vcf                   │
│     [User remaps]                                               │
│     python -m src.mapping filter-remapped remapped.bam         │
│     → wasp_filt.bam (unbiased)                                 │
└─────────────────────┬──────────────────────────────────────────┘
                      │
┌─────────────────────▼──────────────────────────────────────────┐
│  3. COUNT ALLELES (WASP Counting Module)                       │
│     python -m src.counting count-variants \                    │
│       wasp_filt.bam variants.vcf \                             │
│       --samples SAMPLE123 \                                     │
│       --region atac_peaks.bed                                   │
│     → counts.tsv                                                │
└─────────────────────┬──────────────────────────────────────────┘
                      │
┌─────────────────────▼──────────────────────────────────────────┐
│  4. ANALYZE IMBALANCE (WASP Analysis Module)                   │
│     python -m src.analysis find-imbalance counts.tsv           │
│     → ai_results.tsv                                            │
└─────────────────────┬──────────────────────────────────────────┘
                      │
                   RESULTS
```

### Single-Cell Workflow

```
scATAC-seq BAM (with CB tags) + VCF
          │
          ├─→ [Optional] WASP filtering (mapping module)
          │
          ├─→ count-variants-sc (counting module)
          │     → allele_counts.h5ad (AnnData)
          │          cells × SNPs matrix
          │
          ├─→ find-imbalance-sc (analysis module)
          │     + barcode_map.tsv (cell → celltype)
          │     → ai_results_Tcell.tsv
          │     → ai_results_Bcell.tsv
          │        (one file per celltype)
          │
          └─→ compare-imbalance (analysis module)
                → ai_comparison_Tcell_vs_Bcell.tsv
                   (differential AI between celltypes)
```

---

## Design Patterns

### 1. **Orchestrator Pattern**

**Structure**:
```
__main__.py (CLI)
    ↓
run_*.py (Orchestrator)
    ↓
core_logic.py (Worker functions)
```

**Example** (Counting):
```python
# __main__.py - CLI interface
@app.command()
def count_variants(bam, vcf, ...):
    run_count_variants(bam, vcf, ...)

# run_counting.py - Orchestrator
def run_count_variants(bam, vcf, ...):
    # 1. Validate inputs
    # 2. Filter VCF
    # 3. Parse genes
    # 4. Count alleles
    # 5. Write output

# count_alleles.py - Worker
def count_alleles(bam, vcf, ...):
    # Core counting logic
```

**Benefits**:
- Separation of concerns
- CLI logic separate from business logic
- Easier testing (test orchestrator independently)

### 2. **File-Based Data Passing**

**Pattern**: Modules communicate via files, not in-memory objects

```
Module A → writes file.txt → Module B reads file.txt
```

**Rationale**:
- Language-agnostic (could rewrite modules in other languages)
- Checkpointing (restart pipeline mid-way)
- Debugging (inspect intermediate files)

**Trade-offs**:
- I/O overhead
- Type safety only at file format level
- No compiler checking of interfaces

### 3. **Temporary File Management**

**Pattern**: Decorator for temp file cleanup

```python
@tempfile_wrapper(temp_loc, keep_temp)
def process_data(...):
    # Use temp files
    # Decorator handles cleanup
```

**Issues Found**:
- Not consistently used across modules
- Some temp files may leak
- Error paths may not clean up

### 4. **Lazy Loading for Large Files**

**Pattern**: Use polars lazy evaluation

```python
# Don't load entire file
df = pl.scan_csv("huge.vcf")
      .filter(pl.col("QUAL") > 30)
      .collect()  # Only now loads filtered data
```

**Benefits**: Memory efficiency

**Usage**: Inconsistent - not used everywhere it could be

---

## Technology Stack

### Layer Architecture

```
┌────────────────────────────────────────────────┐
│             CLI Layer (Typer)                  │
│  User-facing commands, argument parsing        │
└─────────────────┬──────────────────────────────┘
                  │
┌─────────────────▼──────────────────────────────┐
│         Orchestration Layer (run_*.py)         │
│  Workflow management, validation, I/O          │
└─────────────────┬──────────────────────────────┘
                  │
┌─────────────────▼──────────────────────────────┐
│         Business Logic Layer                   │
│  Counting, statistics, allele swapping         │
└─────────────────┬──────────────────────────────┘
                  │
┌─────────────────▼──────────────────────────────┐
│         Data Access Layer                      │
│  pysam: BAM/VCF I/O                           │
│  polars/pandas: DataFrame operations           │
│  scipy: Statistical functions                  │
└────────────────────────────────────────────────┘
```

### Dependency Usage by Layer

| Layer | Dependencies | Purpose |
|-------|--------------|---------|
| CLI | typer, typing_extensions | Command interface |
| Orchestration | pathlib, tempfile | File management |
| Business Logic | numpy, polars, pandas | Data processing |
| Statistics | scipy, numpy | Beta-binomial models |
| Bioinformatics | pysam, pybedtools | File I/O |
| Single-Cell | anndata, scipy.sparse | H5AD format |

---

## Integration Points

### Module Interfaces

#### Counting → Analysis

**Contract**: TSV file with specific columns

```
Required columns:
- chrom: str
- pos: int
- ref: str (single base)
- alt: str (single base)
- ref_count: int
- alt_count: int

Optional columns:
- region: str (for ATAC peaks)
- gene_id: str (for RNA-seq)
- transcript_id: str (for RNA-seq)
```

**Issue**: ⚠️ No formal schema validation

#### Mapping → Counting

**Contract**: Filtered BAM file

- Must be sorted
- Must be indexed (.bai file)
- Compatible with pysam

**Issue**: ⚠️ No verification that BAM is WASP-filtered

#### Single-Cell Format

**Contract**: AnnData H5AD file

```python
adata structure:
    .X: sparse matrix (cells × SNPs) with allele counts
    .obs: cell metadata
    .var: SNP metadata (chrom, pos, ref, alt)
    .uns: unstructured metadata (samples, etc.)
```

### External Tool Integration

**Mapping Module Step 2**: User must run external aligner

```bash
# WASP provides input FASTQs
# User runs their aligner of choice:
bwa mem ...
# or
STAR ...
# or
bowtie2 ...
```

**Issue**: ⚠️ No validation of user's aligner output

---

## Architectural Issues

### 🔴 Critical Issues

1. **bin/WASP2 Main Executable is Broken**
   - Only supports 2 of 3 modules
   - Hardcoded to wrong module
   - Contains TODO comment
   - **Impact**: Users can't use main executable
   - **Recommendation**: Fix or deprecate

2. **No Test Suite**
   - Zero unit tests found
   - Zero integration tests
   - No CI/CD
   - **Impact**: Regression risk, hard to refactor
   - **Recommendation**: Phase 2 priority

3. **No Shared Utilities Library**
   - Code duplication between modules
   - E.g., `count_alleles.py` exists in both counting/ and analysis/
   - **Impact**: Bug fixes must be applied multiple times
   - **Recommendation**: Extract common code to shared lib

### 🟡 Medium Issues

4. **Inconsistent Error Handling**
   - No custom exception hierarchy
   - Inconsistent try/except patterns
   - **Impact**: Poor user experience, hard debugging
   - **Recommendation**: Standardize in Phase 2

5. **No Type Hints**
   - Functions lack type annotations
   - **Impact**: IDE support poor, runtime errors
   - **Recommendation**: Add type hints in Phase 2

6. **Pandas + Polars Redundancy**
   - Both used throughout
   - Unclear when to use which
   - **Impact**: Larger dependency footprint, confusion
   - **Recommendation**: Document guidelines or consolidate

7. **Limited Configuration**
   - Hard-coded constants
   - No config file support
   - **Impact**: Users can't tune parameters easily
   - **Recommendation**: Add config system

### 🟢 Minor Issues

8. **TODO Comments in Production Code**
   - "TODO GOTTA TEST THIS" in all 3 `__main__.py` files
   - **Impact**: Unclear if CLIs are tested
   - **Recommendation**: Remove or implement

9. **No Logging Framework**
   - Print statements instead of logging
   - No log levels
   - **Impact**: Production debugging difficult
   - **Recommendation**: Implement logging in Phase 2

10. **Temporary File Leaks Possible**
    - Cleanup not guaranteed in error paths
    - **Impact**: Disk space issues
    - **Recommendation**: Use context managers consistently

---

## Architectural Strengths

### ✅ Good Design Choices

1. **Modular Structure**
   - Clean separation between counting/analysis/mapping
   - Can use modules independently

2. **File-Based Interfaces**
   - Language-agnostic
   - Easy to checkpoint and debug
   - Users can inspect intermediate files

3. **Typer CLI Framework**
   - Modern, type-safe CLI
   - Auto-generated help
   - Good user experience

4. **Single-Cell Support**
   - Well-integrated scRNA/scATAC workflows
   - Uses standard formats (H5AD/AnnData)
   - Compatible with scanpy ecosystem

5. **Statistical Rigor**
   - Beta-binomial model is appropriate
   - FDR correction included
   - Phased/unphased support

---

## Recommended Architectural Improvements

### Phase 2 Priorities

1. **Create Shared Utilities Package**
   ```
   src/
   ├── common/
   │   ├── __init__.py
   │   ├── bam_utils.py
   │   ├── vcf_utils.py
   │   ├── file_utils.py
   │   └── validation.py
   ├── counting/
   ├── analysis/
   └── mapping/
   ```

2. **Add Schema Validation**
   - Validate TSV files before processing
   - Check BAM/VCF format compliance
   - Clear error messages for invalid inputs

3. **Implement Logging**
   ```python
   import logging
   logger = logging.getLogger(__name__)
   logger.info("Processing BAM file...")
   ```

4. **Add Type Hints**
   ```python
   def count_alleles(
       bam_path: Path,
       vcf_path: Path,
       min_qual: int = 20
   ) -> pl.DataFrame:
       ...
   ```

5. **Fix bin/WASP2**
   - Properly route to modules
   - Or document as deprecated, use `python -m src.module` instead

---

## Future Architecture Considerations

### Scalability

**Current Limitations**:
- Single-threaded processing (mostly)
- Limited parallelization
- Memory-intensive for large files

**Potential Improvements**:
- Parallelize chromosome processing
- Use Dask for out-of-core computing
- Implement streaming for BAM files

### Python API

**Current**: CLI-only interface

**Future**: Provide programmatic API
```python
from wasp2 import count_variants, analyze_imbalance

counts = count_variants(bam="data.bam", vcf="variants.vcf")
results = analyze_imbalance(counts, min_count=10)
```

### Plugin Architecture

**Future**: Allow custom statistical models
```python
from wasp2.analysis import register_model

@register_model("custom-beta-binomial")
def my_model(counts, params):
    ...
```

---

## Summary

### Architecture Type
**Modular Pipeline Architecture** with file-based interfaces

### Strengths
- Clear module boundaries
- Flexible workflows
- Standard file formats
- Statistical rigor

### Weaknesses
- Code duplication
- No tests
- Broken main executable
- Limited configuration

### Phase 1 Status
- ✅ **Task 1.1.1**: Directory structure documented
- ✅ **Task 1.1.2**: Dependencies analyzed
- ✅ **Task 1.1.3**: Architecture documented (this file)

### Next Steps
- Commit Phase 1.1 deliverables
- Begin Phase 1.2: Module deep dives (starting with Counting)

---

**Document Version**: 1.0
**Last Updated**: 2025-11-15
**Next Review**: After Phase 1.2 (Counting module deep dive)
