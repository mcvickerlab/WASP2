# WASP2-final v1.3.0 Pre-Merge Validation Report

**Date:** 2026-02-19
**Repository:** Jaureguy760/WASP2-final
**Target upstream:** mcvickerlab/WASP2
**Environment:** RHEL 9 HPC, Rust 1.91.1, Python 3.10.10 (WASP2_dev2 mamba env)

---

## Executive Summary

| Phase | Status | Blocking? |
|-------|--------|-----------|
| 1A. Rust Tests | PASS (97/97) | -- |
| 1B. Rust Clippy | FAIL (style/perf) | No (correctness OK) |
| 1C. Python Lint (ruff) | FAIL (6 errors) | Low |
| 1D. Python Format (ruff) | FAIL (16 files) | Low |
| 1E. Rust Format (cargo fmt) | FAIL (7 files) | Low |
| 1F. Pre-commit Hooks | FAIL (multiple) | Low |
| 1G. Bandit (Python security) | PASS (0 issues) | -- |
| 1H. Cargo Audit (Rust security) | **FAIL (2 vulns)** | **YES** |
| 1I. Sanity Tests (chr21) | PASS (8/8) | -- |
| 1J. Docs Build (Sphinx) | PASS (12 warnings) | No |
| 2A. Container Configs | PASS (v1.3.0 consistent) | -- |
| 2B. Singularity.def | PASS (correct tag) | -- |
| 3. Nextflow Pipelines | **FAIL (all 4)** | **YES** |
| 4A. Galaxy Tools | FAIL (version mismatch) | Medium |
| 4B. Bioconda Recipes | PASS (review only) | -- |
| Python Test Suite (prior) | PASS (80/80, 84 skip) | -- |

**Verdict: NOT ready for merge.** Two blocking issues (security vulns, Nextflow pipelines) plus several medium-priority items must be resolved first.

---

## Phase 1: Local Validation

### 1A. Rust Tests (`cargo test`) -- PASS

- **97 passed, 0 failed, 7 ignored** (+ 3 ignored doc-tests)
- Modules tested: analysis (6), bam_filter (3), bam_remapper (30), bam_counter (1), bam_intersect (4), cigar_utils (5), mapping_filter (10), multi_sample (9), seq_decode (2), unified_pipeline (6), vcf_to_bed (4)
- 4 compiler warnings: unused imports in test code, 1 deprecated function use (`apply_allele_substitutions` -> `apply_allele_substitutions_cigar_aware`)

### 1B. Rust Clippy (`cargo clippy -- -D warnings`) -- FAIL

Exit code 101. All issues are **style/performance**, not correctness:

| Lint | Count | Files |
|------|-------|-------|
| `too_many_arguments` (15/7) | 6+ | lib.rs, bam_remapper.rs, unified_pipeline.rs |
| `needless_range_loop` | multiple | bam_remapper.rs, multi_sample.rs |
| `assign_op_pattern` | 2 | analysis.rs:235,242 |
| `cast_abs_to_unsigned` | 2 | vcf_to_bed.rs:245,309 |
| `manual_contains` | 2 | vcf_to_bed.rs:282,283 |
| `field_reassign_with_default` | 1 | unified_pipeline.rs:1772 |
| `needless_borrows_for_generic_args` | 1 | lib.rs:774 |
| `manual_pattern_char_comparison` | multiple | vcf_to_bed.rs, bam_intersect.rs |
| Other (collapsible_if, unnecessary_cast, etc.) | many | various |

**Note:** The `too_many_arguments` in `lib.rs:545,620` are PyO3 function bindings mirroring Python signatures -- these are hard to refactor without breaking the Python API. Consider `#[allow(clippy::too_many_arguments)]` for those.

### 1C-1E. Python & Rust Formatting -- FAIL

| Check | Result |
|-------|--------|
| `ruff check src/ tests/` | 6 errors (4 auto-fixable: F401 unused import, F541 f-string no placeholders, B007 unused loop vars) |
| `ruff format --check` | 16 Python files need reformatting |
| `cargo fmt --check` | 7 Rust files need formatting |

**Fix:** `ruff check --fix src/ tests/ && ruff format src/ tests/ && cd rust && cargo fmt`

### 1F. Pre-commit Hooks -- FAIL

| Hook | Status |
|------|--------|
| ruff | FAIL (15 errors incl. benchmarking/docs E402, B904) |
| ruff-format | FAIL (29 files) |
| trailing whitespace | PASS |
| end-of-file-fixer | FAIL (9 files, notebooks) |
| check-yaml | FAIL (bioconda meta.yaml Jinja templates) |
| large files / merge conflicts | PASS |
| private key / bandit / gitleaks | PASS |
| basedpyright | FAIL (60 errors, 3 warnings) |

The **basedpyright** errors are typical for scientific Python (numpy/pandas/polars type inference, PyO3 binding stubs). Not blocking.

### 1G. Bandit (Python Security) -- PASS

Zero issues across 7,908 lines. Config excludes B101, B603, B607, B404 per `pyproject.toml`.

### 1H. Cargo Audit (Rust Security) -- **FAIL (BLOCKING)**

**2 NEW vulnerabilities requiring immediate action:**

| Crate | Version | Advisory | Fix |
|-------|---------|----------|-----|
| **bytes** | 1.11.0 | RUSTSEC-2026-0007 (integer overflow in BytesMut::reserve) | Upgrade to >= 1.11.1 |
| **pyo3** | 0.28.1 | RUSTSEC-2026-0013 (type confusion with abi3 + Python 3.12+) | Upgrade to >= 0.28.2 |

**Action:** Update `rust/Cargo.toml`:
```toml
pyo3 = "0.28.2"  # was "0.28"
# bytes is transitive -- update Cargo.lock with `cargo update -p bytes`
```

3 known informational warnings (custom_derive RUSTSEC-2025-0058, paste RUSTSEC-2024-0436, lru RUSTSEC-2026-0002) -- all transitive, no upstream fixes available.

### 1I. Sanity Tests (chr21 real data) -- PASS

**8/8 passed** in 8.63s against HG00731 chr21 data:
- TestAlleleCounts (2 tests)
- TestFastqGeneration (2 tests)
- TestAnalysis (3 tests)
- TestPipelineIntegration (1 test)

### 1J. Documentation Build (Sphinx) -- PASS

`docs/build/html/index.html` generated (55KB, 27 pages). 12 non-blocking warnings:
- 1 malformed RST table in `quickstart_mapping.rst:122`
- 6 unreferenced citations
- 3 unknown 'nextflow' Pygments lexer in `seqera_ai_integration.md`
- 1 missing cross-reference to `WASP2_ECOSYSTEM.md`

---

## Phase 2: Container Validation (Review Only)

### 2A. Container Configs -- PASS

All version strings consistent at **v1.3.0**:

| File | Version Source | Consistent? |
|------|--------------|-------------|
| `rust/Cargo.toml` | `version = "1.3.0"` | Baseline |
| `Dockerfile` | `ARG VERSION=1.3.0` | Yes |
| `Dockerfile.optimized` | `ARG VERSION=1.3.0` | Yes |
| `Singularity.def` | `From: jaureguy760/wasp2:1.3.0` | Yes |

Base images pinned: `rust:1.87-bookworm`, `python:3.11-slim-bookworm`, `uv:0.9.26`.
Health checks present. Multi-stage builds correct.

**Minor observations:**
- `Dockerfile.optimized` has `debug = true` in release profile (intentional for profiling?)
- `scripts/check-version-consistency.sh` referenced but missing
- typer version skew between Dockerfile and Dockerfile.optimized

### 2B. Singularity.def -- PASS

Correctly bootstraps from `docker://jaureguy760/wasp2:1.3.0`. Test section runs `wasp2-count --version`.

---

## Phase 3: Nextflow Pipeline Tests -- **FAIL (BLOCKING)**

Nextflow v25.10.4 installed via separate conda env (nextflow_env with Java 17).

| Pipeline | Status | Error |
|----------|--------|-------|
| **nf-atacseq** | FAIL | `fromSamplesheet` not defined by nf-schema@2.6.1 (API breaking change from nf-validation) |
| **nf-rnaseq** | BLOCKED | Test data is stubs (36-byte FASTQs, empty STAR index dir) |
| **nf-outrider** | FAIL | `params.outrider_q` evaluates to null -- null channel error at `outrider.nf:122` |
| **nf-scatac** | BLOCKED | Stub test data only |

**Root causes:**
1. **nf-atacseq:** Uses `fromSamplesheet` channel factory from the deprecated `nf-validation` plugin, but `nf-schema@2.6.1` (which Nextflow auto-downloads) removed this API. Need to migrate to `nf-schema` 2.x `samplesheetToList()` or pin nf-schema to 1.x.
2. **nf-atacseq test config:** References incorrect nf-core test-datasets URLs (wrong filenames/paths).
3. **nf-rnaseq/nf-scatac:** Test data directories contain only stub/placeholder files -- need real minimal test data.
4. **nf-outrider:** Pipeline code bug -- `params.outrider_q` defaults to `null` (for auto-estimation) but the `val` input declaration rejects null values. Also has `first` operator warnings on value channels.

---

## Phase 4: Packaging Validation

### 4A. Galaxy Tools -- FAIL (Version Mismatch)

`planemo lint` output: All 4 Galaxy XMLs have warnings (XMLOrder, BioToolsValid). These are non-blocking.

**Critical finding:** `galaxy/tools/wasp2/macros.xml` defines:
```xml
<token name="@TOOL_VERSION@">1.2.0</token>
```
This **must be updated to 1.3.0** before release.

Galaxy tool XMLs affected:
- `wasp2_count_variants.xml`
- `wasp2_filter_remapped.xml`
- `wasp2_find_imbalance.xml`
- `wasp2_make_reads.xml`

### 4B. Bioconda Recipes -- PASS (Review)

Two recipe variants exist:
- `bioconda/meta.yaml` -- Simpler, v1.3.0 correct, SHA256 placeholder
- `bioconda-recipe/meta.yaml` -- Comprehensive, v1.3.0 correct, SHA256 placeholder, bio.tools identifier

**Recommendation:** Consolidate to single recipe (`bioconda-recipe/meta.yaml` is more complete). SHA256 hash must be filled after PyPI release.

---

## Priority Action Items

### Blocking (must fix before merge)

1. **Cargo dependency security patches** -- Update `rust/Cargo.toml` for pyo3 >= 0.28.2 and run `cargo update -p bytes` for bytes >= 1.11.1. Rebuild and re-audit.

2. **Nextflow nf-atacseq migration** -- Migrate from `fromSamplesheet` (nf-validation) to `samplesheetToList` (nf-schema 2.x) or pin nf-schema to 1.x in `nextflow.config`.

3. **Nextflow nf-outrider null channel** -- Fix `params.outrider_q` handling in `workflows/outrider.nf:122` to accept null values (use `params.outrider_q ?: 'auto'` pattern).

4. **Nextflow test data** -- Populate real minimal test data for nf-rnaseq and nf-scatac pipelines (or convert to `-stub-run` compatible).

### High Priority (should fix before merge)

5. **Galaxy version bump** -- Update `macros.xml` @TOOL_VERSION@ from 1.2.0 to 1.3.0.

6. **Run formatters** -- `ruff format && cargo fmt` to clean up all formatting.

7. **Fix ruff lint errors** -- 6 errors in sanity tests, 9 in benchmarking/docs.

### Low Priority (can fix post-merge)

8. **Clippy cleanup** -- Address style lints across Rust code, allow `too_many_arguments` for PyO3 bindings.

9. **basedpyright configuration** -- Add type stubs or pyright overrides for numpy/pandas/polars/PyO3.

10. **Docs warnings** -- Fix RST table in `quickstart_mapping.rst`, add nextflow Pygments lexer.

11. **Consolidate bioconda recipes** -- Choose one of the two meta.yaml variants.

---

## Tests Summary

| Test Suite | Passed | Failed | Skipped |
|-----------|--------|--------|---------|
| Python (pytest) | 80 | 0 | 84 |
| Rust (cargo test) | 97 | 0 | 10 |
| Sanity (chr21) | 8 | 0 | 0 |
| **Total** | **185** | **0** | **94** |

All functional tests pass. Zero test failures across the entire codebase.
