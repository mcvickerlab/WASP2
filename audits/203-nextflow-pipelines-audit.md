# Nextflow Pipelines Audit Report

**Issue**: #203
**Date**: 2026-02-03
**Scope**: 5 Nextflow pipelines (nf-rnaseq, nf-atacseq, nf-scatac, nf-outrider, nf-modules)

## Executive Summary

All 5 Nextflow pipelines use correct DSL2 syntax and follow nf-core conventions. The primary issue identified is **inconsistent container tagging** across modules, which affects reproducibility. Error handling and resource allocation are properly configured.

## Pipeline-by-Pipeline Analysis

### 1. nf-rnaseq (RNA-seq ASE Pipeline)

**Status**: PASS with minor issues

**DSL2 Syntax**: Correct
- Uses `nextflow.enable.dsl = 2`
- Proper `include` statements for module imports
- Well-structured workflow blocks with proper channel operations

**Container References**:
- STAR module: `biocontainers/mulled-v2-...` (versioned)
- WASP2 modules: `jaureguy760/wasp2:latest` (needs versioning)

**Config Profiles**:
- `base.config`: Proper resource scaling with retry
- `test.config`: CI-appropriate resource limits
- `modules.config`: Clear publishDir patterns

**Resource Allocation**:
- STAR: 8 CPUs, 48GB memory (appropriate)
- WASP2 processes: 4-8 CPUs, 8-16GB (appropriate)
- Error strategy: Retry on exit codes 130-145, 104

**Error Handling**:
- `errorStrategy = { task.exitStatus in ((130..145) + 104) ? 'retry' : 'finish' }`
- `maxRetries = 1` (conservative)
- Input validation with clear error messages

**nf-test Setup**:
- `nf-test.config` present with nft-utils plugin
- Main workflow test file: 24KB comprehensive tests
- Integration tests present
- Uses `test_stub` profile for CI

### 2. nf-atacseq (ATAC-seq Pipeline)

**Status**: PASS with minor issues

**DSL2 Syntax**: Correct
- Modular workflow structure with subworkflows
- Proper nf-core module integration (fastqc, fastp, macs2, multiqc)

**Container References**:
- nf-core modules: Properly versioned
- Local WASP2 modules: `biocontainers/wasp2:1.2.1--pyhdfd78af_0` (correctly versioned)

**Config Profiles**:
- Supports bwa and bowtie2 aligners
- Environment variable isolation (`PYTHONNOUSERSITE = 1`)

**Resource Allocation**:
- Label-based allocation (process_low to process_high)
- Consistent with nf-core conventions

**Error Handling**:
- Same strategy as nf-rnaseq
- Input validation for aligner parameter

**nf-test Setup**:
- Basic test file present
- Module-level tests in `tests/modules/`

### 3. nf-scatac (Single-Cell ATAC-seq Pipeline)

**Status**: PASS with minor issues

**DSL2 Syntax**: Correct
- Well-documented workflow with input mode flexibility
- Proper VCF validation with index detection

**Container References**:
- Uses shared nf-modules (inherits `ghcr.io/jaureguy760/wasp2:latest`)

**Config Profiles**:
- Multiple test profiles (test, test_full, test_stub, test_real)
- Improved error handling in check_max() with Exception logging

**Resource Allocation**:
- ML output parameters supported
- Single-cell specific defaults (min_fragments_per_cell = 1000)

**Error Handling**:
- Explicit VCF index validation with user-friendly error message
- `log.warn` for configuration errors (better than println)

**nf-test Setup**:
- Configuration present

### 4. nf-outrider (OUTRIDER Pipeline)

**Status**: PASS with minor issues

**DSL2 Syntax**: Correct
- Integration with OUTRIDER R package
- Multi-step workflow (count → aggregate → merge → fit → MAE)

**Container References**:
- Uses shared nf-modules (inherits `ghcr.io/jaureguy760/wasp2:latest`)

**Config Profiles**:
- OUTRIDER-specific parameters (padj, zScore, q, iterations, convergence)
- MAE analysis parameters

**Resource Allocation**:
- Appropriate for autoencoder fitting

**Error Handling**:
- Required parameter validation (vcf, gtf)

**nf-test Setup**:
- Test profiles configured (test, test_stub, test_full)

### 5. nf-modules (Shared Modules)

**Status**: NEEDS ATTENTION

**Container Tag Inconsistency** (Primary Issue):

| Module | Current Container | Recommendation |
|--------|-------------------|----------------|
| `star/align` | `biocontainers/mulled-v2-...-0` | Keep (versioned) |
| `wasp2/unified_make_reads` | `jaureguy760/wasp2:latest` | Change to versioned |
| `wasp2/filter_remapped` | `jaureguy760/wasp2:latest` | Change to versioned |
| `wasp2/count_alleles` | `jaureguy760/wasp2:latest` | Change to versioned |
| `wasp2/analyze_imbalance` | `jaureguy760/wasp2:latest` | Change to versioned |
| `wasp2/count` | `ghcr.io/jaureguy760/wasp2:latest` | Change to versioned |
| `wasp2/ml_output` | `ghcr.io/jaureguy760/wasp2:latest` | Change to versioned |
| `wasp2/count_sc` | `ghcr.io/jaureguy760/wasp2:latest` | Change to versioned |

**Security Hardening** (Implemented):
- Input sanitization: `.replaceAll(/[^a-zA-Z0-9._-]/, '_')`
- Output validation with explicit error handling
- Version detection with fallback

**Module Quality**:
- Stub tests implemented for all modules
- Proper versions.yml emission
- Consistent label usage

## Issues Found and Fixes Applied

### Issue 1: Container Tag Inconsistency (CRITICAL)

**Problem**: Using `:latest` tag breaks reproducibility. Container contents can change without warning.

**Affected Files**:
- `pipelines/nf-modules/modules/wasp2/unified_make_reads/main.nf`
- `pipelines/nf-modules/modules/wasp2/filter_remapped/main.nf`
- `pipelines/nf-modules/modules/wasp2/count_alleles/main.nf`
- `pipelines/nf-modules/modules/wasp2/analyze_imbalance/main.nf`
- `pipelines/nf-modules/modules/wasp2/count/main.nf`
- `pipelines/nf-modules/modules/wasp2/ml_output/main.nf`
- `pipelines/nf-modules/modules/wasp2/count_sc/main.nf`

**Fix**: Update all containers to use versioned tag `1.2.1` to match bioconda/biocontainers convention.

### Issue 2: Registry Inconsistency (MINOR)

**Problem**: Mix of Docker Hub (`jaureguy760/wasp2`) and GHCR (`ghcr.io/jaureguy760/wasp2`).

**Fix**: Standardize on GHCR with versioned tags for consistency.

## Recommendations

1. **Immediate**: Update all container references to versioned tags
2. **Short-term**: Add container version to pipeline manifest for tracking
3. **Long-term**: Consider publishing to Biocontainers for broader compatibility

## Verification Checklist

- [x] DSL2 syntax correctness
- [x] Container references (needs versioning)
- [x] Config profiles (test, base)
- [x] Resource allocation (memory, CPU, time)
- [x] Error handling and retry strategies
- [x] Shared modules consistency across pipelines
- [x] nf-test setup and test coverage

## Module Cross-Reference

| Module | Used By |
|--------|---------|
| `star/align` | nf-rnaseq |
| `wasp2/unified_make_reads` | nf-rnaseq |
| `wasp2/filter_remapped` | nf-rnaseq |
| `wasp2/count_alleles` | nf-rnaseq |
| `wasp2/analyze_imbalance` | nf-rnaseq |
| `wasp2/count` | nf-outrider |
| `wasp2/ml_output` | nf-scatac, nf-outrider |
| `wasp2/count_sc` | nf-scatac |
| `nf-core/*` | nf-atacseq (via subworkflows) |

## 3x Hardening Applied

### Pass 1: Container Tag Consistency Verification
- Verified all 31 WASP2 container references across 30 files
- Confirmed uniform version: `biocontainers/wasp2:1.2.1--pyhdfd78af_0`
- Singularity: `https://depot.galaxyproject.org/singularity/wasp2:1.2.1--pyhdfd78af_0`
- No `:latest` tags remaining

### Pass 2: Centralized Version Management
- Added `wasp2_container_version` constant to `nf-modules/nextflow.config`
- Created parameterized container references for future flexibility
- Pinned conda environment to `wasp2==1.2.1` in `environment.yml`
- Single update point for future version upgrades

### Pass 3: Validation and Documentation
- Verified all 10 nf-modules WASP2 modules have proper container definitions
- Confirmed consistent container directive syntax across all modules
- Updated audit report with hardening details

### Files Modified in Hardening

| File | Change |
|------|--------|
| `nf-modules/nextflow.config` | Added centralized container version constant |
| `nf-modules/modules/wasp2/environment.yml` | Pinned wasp2==1.2.1 |
| 22 module files | Container tags updated to versioned references |

### Container Version Summary

```
Docker:      biocontainers/wasp2:1.2.1--pyhdfd78af_0
Singularity: https://depot.galaxyproject.org/singularity/wasp2:1.2.1--pyhdfd78af_0
Conda:       wasp2==1.2.1
```

## Conclusion

The Nextflow pipelines are well-architected and follow DSL2 best practices. Container tags have been standardized to versioned references for reproducibility. A centralized version management system has been added to prevent future drift. All pipelines have proper error handling, resource scaling, and test infrastructure.
