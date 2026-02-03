# Audit: Test Suite Quality and Coverage Analysis

**Issue:** #200
**Date:** 2026-02-02
**Scope:** `tests/`, `rust/src/` (test modules), CI/CD test configuration

---

## Executive Summary

The WASP2 test suite contains **254 test items** (154 Python, ~100 Rust) across 16 Python
test files and 12 Rust modules. Overall structure is sound, but several significant quality
issues need attention: two test files bypass pytest entirely, test isolation is weak in
places, coverage gaps exist in critical paths (statistical analysis, single-cell, pipeline
CLI runners), and the Rust side lacks integration tests.

| Metric | Value | Assessment |
|--------|-------|------------|
| Python test files | 16 | Adequate |
| Python test items | 154 | Moderate |
| Rust test items | ~100 | Good for unit tests |
| Rust integration tests | 0 | Gap |
| `conftest.py` size | 227 lines (7.3 KB) | Reasonable — not the 7435-line file feared |
| `test_indel_correctness.py` size | 340 lines (12 KB) | Reasonable — not generated bloat |
| CI test execution | pytest + cargo test | Good |

---

## 1. `conftest.py` Review (227 lines, 7.3 KB)

**Verdict: Well-structured, no splitting needed.**

The issue flagged `conftest.py` as "7435 lines" — this is incorrect. The actual file is
227 lines containing:

- 5 session-scoped fixtures (test data paths, VCF/PGEN file generation)
- 3 function-scoped fixtures (temp dirs, expected variants)
- 4 custom marker registrations
- 3 helper functions (`has_command`, `skip_without_*`)

**Findings:**

| ID | Finding | Severity |
|----|---------|----------|
| C-1 | `sample_vcf_gz` fixture opens a file without `with` statement (line 91: `stdout=open(vcf_gz_path, "wb")`) — resource leak on failure | Low |
| C-2 | Markers registered in both `conftest.py` and `pyproject.toml` (redundant) | Informational |
| C-3 | `benchmarks/conftest.py` (14 KB) is larger and contains synthetic data generators — appropriate for its purpose | OK |

---

## 2. `test_indel_correctness.py` Review (340 lines, 12 KB)

**Verdict: Valid tests, but structural issues.**

The issue flagged this as "12128 lines" — incorrect. The actual file is 340 lines with
10 well-crafted correctness tests for INDEL handling (position mapping, quality filling,
phased sequence building, multi-sample).

**Findings:**

| ID | Finding | Severity |
|----|---------|----------|
| IC-1 | Contains `run_all_tests()` manual runner with a **typo** on line 316: `except AssertionError` (missing 'r') — this dead code silently catches `NameError` instead of `AssertionError`, so failures would show as "ERROR" not "FAIL" | Medium |
| IC-2 | Tests use `print()` statements with emoji output (✅/❌) — noise in pytest output, non-functional | Low |
| IC-3 | Uses `sys.path.insert(0, ...)` instead of proper package installation — fragile | Low |
| IC-4 | Tests are well-written with meaningful assertions and clear documentation | Strength |

---

## 3. Test Files That Bypass pytest

**Two files are scripts, not pytest test modules:**

### `test_rust_python_match.py` (203 lines)

| ID | Finding | Severity |
|----|---------|----------|
| RP-1 | **Zero pytest-collectible tests.** File executes comparison code at module import time (lines 17-202 run as top-level statements). pytest would collect 0 tests from this file. | High |
| RP-2 | Uses `global passed, failed` counters instead of assertions — non-standard test pattern | High |
| RP-3 | `test_validation_quick.py` wraps this in a subprocess call and handles the `returncode == 5` (no tests collected) case by skipping — confirms the issue is known | Medium |

### `test_rust_bam_filter.py` (126 lines)

| ID | Finding | Severity |
|----|---------|----------|
| RB-1 | `test_rust_filter_matches_samtools()` uses bare `return` instead of `pytest.skip()` when data is missing (lines 39-40) — pytest collects it but the test silently passes with no assertions when data is absent | High |
| RB-2 | Uses `print()` for pass/fail instead of `assert` for comparison results (line 92-93 check equality but line 113 returns `False` instead of asserting) | High |
| RB-3 | Depends on external benchmark data (`benchmarking/star_wasp_comparison/`) that likely doesn't exist in most environments — effectively a dead test | Medium |

---

## 4. Test Isolation Issues

| ID | Finding | Severity |
|----|---------|----------|
| TI-1 | `test_validation_quick.py::test_rust_python_parity` spawns a subprocess to run `test_rust_python_match.py` — tests-calling-tests pattern creates hidden dependencies and unclear failure attribution | Medium |
| TI-2 | `test_validation_quick.py::test_indel_correctness` similarly spawns a subprocess — if the inner tests fail, the outer test shows subprocess stderr, not the actual assertion failure | Medium |
| TI-3 | Session-scoped fixtures (`sample_vcf`, `sample_vcf_gz`) write to `tests/data/` which is version-controlled — test runs modify tracked files | Medium |
| TI-4 | `sample_vcf_gz` fixture depends on external tools (bcftools/bgzip) — tests silently skip on systems without these tools, creating environment-dependent coverage | Low |

---

## 5. Coverage Gaps in `src/`

### Python Coverage Analysis

Mapping test files to `src/` modules reveals significant gaps:

| `src/` Module | Test Coverage | Gap Severity |
|---------------|--------------|--------------|
| `wasp2/io/variant_source.py` | `tests/io/test_variant_source.py` (37 tests) | Good |
| `wasp2/io/vcf_source.py` | `tests/io/test_vcf_source.py` (18 tests) | Good |
| `wasp2/io/cyvcf2_source.py` | `tests/io/test_cyvcf2_source.py` (21 tests) | Good |
| `wasp2/io/compat.py` | `tests/io/test_compat.py` (7 tests) | Good |
| `mapping/remap_utils.py` | `test_indel_correctness.py` (10 tests) | Adequate |
| `wasp2/cli.py` | **No tests** | Medium |
| `analysis/as_analysis.py` | **No dedicated tests** | **Critical** |
| `analysis/compare_ai.py` | **No tests** | **Critical** |
| `analysis/as_analysis_sc.py` | **No tests** | **High** |
| `analysis/filter_data.py` | **No tests** | High |
| `counting/count_alleles.py` | **No dedicated Python tests** | High |
| `counting/count_alleles_sc.py` | **No tests** | High |
| `counting/filter_variant_data.py` | **No tests** | Medium |
| `counting/parse_gene_data.py` | **No tests** | Medium |
| `mapping/intersect_variant_data.py` | Tested indirectly via regression | Low |
| `mapping/make_remap_reads.py` | Tested indirectly via regression | Low |
| `mapping/filter_remap_reads.py` | **No tests** | Medium |
| `mapping/wasp_data_files.py` | **No tests** | Low |
| `mapping/run_mapping.py` (CLI) | **No tests** | Low |
| All `run_*.py` CLI runners | **No tests** | Low |

**Key gap:** The statistical analysis core (`as_analysis.py`, `compare_ai.py`) has zero
dedicated tests. These contain the beta-binomial optimization, FDR correction, and
likelihood ratio test logic — the scientific heart of WASP2.

---

## 6. Rust Test Coverage Gaps

| Module | Tests | Assessment |
|--------|-------|------------|
| `bam_remapper.rs` | 42 | Excellent |
| `mapping_filter.rs` | 14 | Very good |
| `multi_sample.rs` | 11 | Very good |
| `cigar_utils.rs` | 6 | Good |
| `unified_pipeline.rs` | 6 | Good |
| `bam_intersect.rs` | 4 | Adequate |
| `vcf_to_bed.rs` | 4 | Adequate |
| `read_pairer.rs` | 4 | Adequate |
| `bam_filter.rs` | 3 | Under-tested |
| `analysis.rs` | 3 | Under-tested |
| `seq_decode.rs` | 2 | Minimal |
| `bam_counter.rs` | 1 | **Critical gap** — main `count_alleles()` API untested |
| `lib.rs` | 0 | **No tests** — PyO3 FFI layer completely untested |

**No integration tests exist** (`rust/tests/` directory absent). All Rust tests are
in-module unit tests.

---

## 7. Regression Tests Assessment

**Verdict: Well-designed but dependent on external baseline data.**

| ID | Finding | Severity |
|----|---------|----------|
| R-1 | `test_pipeline_regression.py` tests all skip when baseline data is absent (every test has `pytest.skip` guard) — CI likely skips all 10 tests | Medium |
| R-2 | Quickbench parity tests (`test_quickbench_snv_parity.py`, `test_quickbench_indel_parity.py`, `test_quickbench_indel_trim_invariants.py`) are well-structured and meaningful — they compare unified Rust path against multi-pass Python baseline | Strength |
| R-3 | Quickbench tests depend on `benchmarking.quickbench` module — skip when not available | Low |
| R-4 | Performance regression thresholds (30% time, 20% memory) are reasonable | Strength |
| R-5 | MD5 checksums for output validation is a strong approach | Strength |

---

## 8. Benchmark Tests Assessment

**Verdict: Comprehensive synthetic benchmarks, well-organized.**

The `tests/benchmarks/` directory contains 4 test files with 25 parametrized tests across
scales (100 to 1M items). `benchmarks/conftest.py` provides synthetic data generators.
`benchmarks/utils/visualization.py` (26 KB) generates performance plots.

| ID | Finding | Severity |
|----|---------|----------|
| B-1 | Benchmarks are correctly excluded from CI (`--ignore=tests/benchmarks/`) | OK |
| B-2 | Dedicated workflow (`benchmarks.yml`) runs on manual dispatch, releases, and weekly schedule | Strength |
| B-3 | Results produce actionable output (timing, memory, scaling curves) | Strength |

---

## 9. Test Interdependencies

| ID | Finding | Severity |
|----|---------|----------|
| D-1 | `test_validation_quick.py` depends on `test_indel_correctness.py` and `test_rust_python_match.py` via subprocess — if either is renamed/moved, the wrapper silently skips | Medium |
| D-2 | Session-scoped fixtures create cascading dependencies: `sample_vcf_gz` depends on `sample_vcf` depends on `test_data_dir` — failure in early fixture cascades to all VCF-dependent tests | Low (expected) |
| D-3 | No circular dependencies detected | OK |

---

## 10. CI/CD Test Configuration

| Aspect | Status |
|--------|--------|
| Python tests in CI | `pytest tests/ -v --tb=short -x --ignore=tests/benchmarks/` |
| Rust tests in CI | `cargo test` |
| Linting | ruff + mypy |
| Security scanning | bandit + cargo-audit |
| Benchmarks | Separate workflow (manual/weekly/release) |
| Coverage reporting | Configured in `pyproject.toml` but **not enforced in CI** |

| ID | Finding | Severity |
|----|---------|----------|
| CI-1 | No coverage threshold enforced in CI — coverage can regress silently | Medium |
| CI-2 | `pytest -x` (fail-fast) means only the first failure is reported per CI run | Low |

---

## Summary of Findings by Severity

### Critical (2)
- **No tests for statistical analysis core** (`as_analysis.py`, `compare_ai.py`) — the scientific heart of WASP2
- **Rust `bam_counter.rs::count_alleles()` has zero tests** for its main public API

### High (4)
- `test_rust_python_match.py` has zero pytest-collectible tests (script, not test module)
- `test_rust_bam_filter.py` silently passes when data is missing (no `pytest.skip`)
- No tests for single-cell analysis (`as_analysis_sc.py`, `count_alleles_sc.py`)
- No Rust integration tests (`rust/tests/` absent)

### Medium (9)
- `test_indel_correctness.py` `run_all_tests()` has typo: `AssertionError` (line 316)
- Tests-calling-tests via subprocess pattern (validation_quick → indel/parity)
- Session fixtures write to version-controlled `tests/data/`
- Regression tests all skip when baseline data absent (CI likely runs 0 regression tests)
- No coverage threshold in CI
- Missing tests for `filter_data.py`, `filter_remap_reads.py`, `filter_variant_data.py`
- Redundant marker registration (conftest.py + pyproject.toml)
- `conftest.py` size was misreported in issue (227 lines, not 7435)
- `test_indel_correctness.py` size was misreported (340 lines, not 12128)

### Low / Informational (5)
- Resource leak in `sample_vcf_gz` fixture (unclosed file handle)
- `print()` noise in test output
- `sys.path.insert` usage instead of proper packaging
- Environment-dependent test skipping (bcftools, plink2)
- CLI runners untested (low risk — thin wrappers)

---

## Recommendations

### Immediate (address in this PR or next)
1. Convert `test_rust_python_match.py` to proper pytest functions
2. Fix `test_rust_bam_filter.py` to use `pytest.skip()` instead of bare `return`
3. Fix typo `AssertionError` → `AssertionError` in `test_indel_correctness.py` line 316
4. Add `pytest-cov` minimum threshold to CI (even 30% to start)

### Short-term
5. Add unit tests for `as_analysis.py` core functions (`opt_prob`, `single_model`, `linear_model`, `compare_results`)
6. Add tests for `bam_counter.rs::count_alleles()`
7. Create `rust/tests/` integration test directory
8. Stop writing session-scoped fixture output to `tests/data/` (use `tmp_path_factory`)

### Long-term
9. Add single-cell analysis test coverage
10. Add property-based tests for statistical functions (hypothesis library)
11. Enforce coverage threshold increase over time (ratchet pattern)
