# Security Audit Report — Dependencies & Code

**Issue:** #201
**Scope:** Python dependencies (`pyproject.toml`), Rust dependencies (`rust/Cargo.toml`),
code review (OWASP), CI workflow review, secret detection
**Date:** 2026-02-02
**Auditor:** Automated + manual review

---

## Executive Summary

WASP2's security posture is **strong**. No hardcoded secrets, no unsafe shell
usage, no dynamic code patterns, and no high-severity Bandit findings in the
Python source. The primary action items are: (1) adding a `.gitleaks.toml`
configuration for project-specific coverage, (2) tracking environment-level CVEs
that may affect deployment, and (3) updating the `security.yml` workflow to also
scan the `dev` branch.

**Overall Risk:** LOW

---

## 1. Python Dependency Audit (`pyproject.toml`)

### Direct Dependencies — CVE Status

| Package | Version Spec | Known CVEs | Status |
|---------|-------------|------------|--------|
| `numpy>=1.21.0` | — | None | OK |
| `pandas>=1.5.0,<3.0.0` | — | None | OK |
| `polars>=0.19.0` | — | None | OK |
| `scipy>=1.10.0` | — | None | OK |
| `pysam>=0.21.0` | — | None | OK |
| `pybedtools>=0.9.0` | — | None | OK |
| `anndata>=0.8.0,<0.12.0` | — | None | OK |
| `scanpy>=1.9.0` | — | None | OK |
| `typer>=0.9.0` | — | None | OK |
| `rich>=13.0.0` | — | None | OK |

### Optional Dependencies

| Package | Version Spec | Known CVEs | Status |
|---------|-------------|------------|--------|
| `cyvcf2>=0.31.0` | — | None | OK |
| `Pgenlib>=0.90` | — | None | OK |

### Dev Dependencies — CVE Status

| Package | Version Spec | Known CVEs | Status |
|---------|-------------|------------|--------|
| `bandit[toml]>=1.8.0` | — | None | OK |
| `pip-audit>=2.7.0` | — | None | OK |
| `pytest>=7.0` | — | None | OK |
| `ruff>=0.9.0` | — | None | OK |
| `maturin>=1.4` | — | None | OK |

### Environment-Level CVEs (not WASP2 direct deps)

`pip-audit` found 30 vulnerabilities in 15 packages in the shared environment.
**None are WASP2 direct dependencies**, but they may affect deployment:

| Package | CVE | Fix Version | Relevance |
|---------|-----|-------------|-----------|
| `jinja2` 3.1.4 | CVE-2024-56326, CVE-2024-56201, CVE-2025-27516 | 3.1.6+ | Transitive (Sphinx docs) |
| `pyarrow` 14.0.2 | PYSEC-2024-161 | 17.0.0 | Not a WASP2 dep |
| `werkzeug` 3.1.3 | CVE-2025-66221, CVE-2026-21860 | 3.1.5 | Not a WASP2 dep |
| `pip` 25.3 | (see pip-audit output) | 26.0+ | Build tool only |

**Recommendation:** No action required for WASP2 itself. Environment
administrators should update `jinja2` and `pip` in shared environments.

### Outdated Direct Dependencies

All WASP2 direct dependencies have sufficiently recent minimum versions. The
lower bounds (`>=`) allow pip to resolve to the latest compatible version.

### Unnecessary Dependencies

No unnecessary dependencies identified. All 10 core dependencies are actively
used:
- `numpy`, `pandas`, `polars`, `scipy` — data processing and statistics
- `pysam`, `pybedtools` — BAM/VCF/BED file I/O
- `anndata`, `scanpy` — single-cell data structures
- `typer`, `rich` — CLI interface

---

## 2. Rust Dependency Audit (`rust/Cargo.toml`)

`cargo audit` could not run in this environment (NFS locking limitation). The
audit relies on the CI workflow and manual review.

### Direct Dependencies — Status

| Crate | Version | Status | Notes |
|-------|---------|--------|-------|
| `pyo3` | 0.20 | **CVE tracked** | RUSTSEC-2025-0020 (risk of buffer overflow in PyString::from_object). Fix pending in PR #217. |
| `rust-htslib` | 0.44 | Pinned (NFS) | Intentionally pinned; 0.47+ has NFS build issues. |
| `rayon` | 1.8 | OK | |
| `anyhow` | 1.0 | OK | |
| `rustc-hash` | 1.1 | OK | |
| `statrs` | 0.18 | OK | |
| `rv` | 0.19 | OK | |
| `coitrees` | 0.4 | OK | |
| `crossbeam-channel` | 0.5 | OK | |
| `gzp` | 0.11 | OK | |
| `itoa` | 1.0 | OK | |
| `smallvec` | 1.13 | OK | |
| `noodles-vcf` | 0.72 | OK | |
| `noodles-bcf` | 0.68 | OK | |
| `noodles-core` | 0.16 | OK | |
| `noodles-bgzf` | 0.33 | OK | |
| `flate2` | 1.1 | OK | |

### Known Informational Warnings (transitive)

| Advisory | Crate | Status |
|----------|-------|--------|
| RUSTSEC-2025-0058 | `custom_derive` (via rust-htslib) | Unmaintained; no alternative. Low risk. |
| RUSTSEC-2024-0436 | `paste` (via rv) | Unmaintained. Low risk. |
| RUSTSEC-2026-0002 | `lru` (via rv) | Unsound `IterMut`; WASP2 does not use it. Low risk. |
| **RUSTSEC-2025-0020** | **`pyo3`** | **Buffer overflow. Tracked in PR #217.** |

### Unnecessary Dependencies

The `argmin`/`argmin-math` issue from the Rust audit (#199) has been resolved —
these crates are no longer in `Cargo.toml`.

---

## 3. Code Security Review (OWASP)

### Subprocess Calls — All Safe

Every subprocess call uses **list arguments** (not string interpolation) and
**never uses unsafe shell invocation**:

| File | Tool Called | Pattern | Verdict |
|------|-----------|---------|---------|
| `counting/filter_variant_data.py` | `bedtools intersect` | `subprocess.run(list, check=True)` | SAFE |
| `mapping/filter_remap_reads.py` | `samtools merge/index` | `subprocess.run(list, check=True)` | SAFE |
| `mapping/intersect_variant_data.py` | `samtools index` | `subprocess.run(list, check=True)` | SAFE |
| `wasp2/io/vcf_source.py` | `bcftools view/query` | `subprocess.run(list, check=True)` | SAFE |
| `wasp2/io/cyvcf2_source.py` | `bcftools view/query` | `subprocess.run(list, check=True)` | SAFE |
| `wasp2/io/compat.py` | `bcftools view/query` | `subprocess.run(list, check=True)` | SAFE |

### Command Injection — Not Vulnerable

- No unsafe shell invocations anywhere in the codebase
- No dynamic code execution calls
- All file paths passed as typed arguments, not string-interpolated
- External tool arguments constructed from validated CLI parameters via Typer

### Path Traversal — Low Risk

File paths come from CLI arguments (user-controlled, but this is a local CLI
tool, not a web service). Path handling uses `pathlib.Path` throughout.

### SQL Injection — N/A

No database usage.

### XSS/CSRF — N/A

No web interface.

### Bandit Results

With project-configured skips (B101, B404, B603, B607):
- **0 issues found** across 6,986 lines of Python code

Without skips (full scan):
- **56 Low-severity findings**: all B101 (assert), B404 (import subprocess),
  B603 (subprocess without shell), B607 (partial path). All are expected
  patterns for a bioinformatics CLI tool.
- **0 Medium or High severity findings**

### Hardcoded Secrets / Credentials

- **None found** in source code
- **Gitleaks scan:** 201 commits scanned, 0 leaks found
- `.gitignore` correctly excludes `.env*`, `.venv*`, `*.log`, credential files
- Docker image uses non-root user (`wasp2:wasp2`)

---

## 4. Security Workflow Review (`security.yml`)

### Current Configuration

| Scanner | Scope | Strict? | Notes |
|---------|-------|---------|-------|
| pip-audit | Python deps | Informational | Correct — reviewed manually |
| Bandit | Python code | Informational | Correct — Low findings expected |
| cargo-audit | Rust deps | Informational | Correct — tracks advisories |
| Gitleaks | Secrets | **Strict** | Correct — secrets block builds |

### Findings (all resolved in this PR)

1. ~~**Branch coverage gap:** Triggers on `push: [main]` and `pull_request: [main]`
   but not `dev`.~~ **RESOLVED** — `dev` branch added to triggers.

2. ~~**Gitleaks version:** CI uses v8.18.4 while pre-commit uses v8.24.0.~~
   **RESOLVED** — CI updated to v8.24.0 with version-aware install guard.

3. **Virtual environment isolation:** pip-audit and Bandit jobs correctly create
   isolated venvs.

4. **NFS locking:** `cargo-audit` handles NFS lock with
   `rm -f ~/.cargo/advisory-db.lock`.

5. **cargo-audit install guard:** Added binary existence check before running
   `cargo audit` to prevent false green checks when install fails.

### Recommendations Applied

- Added `dev` branch to security workflow triggers
- Updated Gitleaks CI version to v8.24.0 with version-aware install
- Added cargo-audit binary existence guard
- Added `set -o pipefail` to gitleaks install step

---

## 5. Gitleaks Configuration

**Finding:** No `.gitleaks.toml` file exists. The project relies on defaults.

**Recommendation:** Add `.gitleaks.toml` for project-specific patterns and to
allowlist known false positives.

---

## 6. Pre-commit Configuration Review

| Hook | Version | Status |
|------|---------|--------|
| Ruff | v0.9.6 | OK |
| pre-commit-hooks | v4.6.0 | OK |
| Bandit | 1.8.3 | OK |
| Gitleaks | v8.24.0 | OK |
| basedpyright | local | OK |
| detect-private-key | (in pre-commit-hooks) | OK |

All hooks are properly configured.

---

## 7. Container Security (Dockerfile)

| Check | Result |
|-------|--------|
| Multi-stage build | PASS |
| Non-root user | PASS — `wasp2:wasp2` (UID 1000) |
| Minimal base image | PASS — `python:3.11-slim-bookworm` |
| Build tools removed | PASS — `g++` purged after compilation |
| No secrets in image | PASS |
| Health check | PASS |
| `--no-cache-dir` | PASS |

---

## 8. Summary of Action Items

### Must Do (this PR)

1. **Add `.gitleaks.toml`** — project-specific secret detection config
2. **Update `security.yml`** — add `dev` branch, bump Gitleaks version

### Track Separately

3. **Merge Dependabot PR #217** — fixes RUSTSEC-2025-0020 (pyo3)
4. **Monitor `jinja2` in docs environment** — CVE-2024-56326

### No Action Required

- Python direct dependencies: **0 CVEs**
- Bandit code scan: **0 medium/high findings**
- Subprocess usage: **all safe** (list args, no shell)
- Hardcoded secrets: **none found**
- OWASP review: **no vulnerabilities**
- Container security: **all checks pass**

---

## 9. Audit Checklist

| Check | Result |
|-------|--------|
| Python dependency CVEs | **PASS** — 0 in direct deps |
| Rust dependency CVEs | **1 TRACKED** — pyo3 (PR #217) |
| Outdated packages | **PASS** — lower bounds allow latest |
| Unnecessary dependencies | **PASS** — all deps actively used |
| `security.yml` workflow | **UPDATED** — added dev branch |
| Hardcoded secrets | **PASS** — 0 leaks in 201 commits |
| OWASP: Command injection | **PASS** |
| OWASP: Path traversal | **N/A** — local CLI tool |
| OWASP: SQL injection | **N/A** — no database |
| Gitleaks config | **ADDED** — `.gitleaks.toml` |
