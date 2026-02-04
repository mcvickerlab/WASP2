# Packaging Audit: Version Consistency

**Date:** 2026-02-03
**Issue:** #205
**Scope:** `pyproject.toml`, `rust/Cargo.toml`, `CHANGELOG.md`, `Dockerfile`, `bioconda-recipe/meta.yaml`, `Singularity.def`

---

## Executive Summary

All packaging files are consistent at version **1.3.0**. Entry points resolve correctly, optional dependency groups match expectations, and PyPI publishing configuration via maturin is properly set up. No issues found.

---

## Version Consistency

**Single source of truth:** `rust/Cargo.toml` (version `1.3.0`)

| File | Version | Method | Status |
|------|---------|--------|--------|
| `rust/Cargo.toml` | 1.3.0 | Hardcoded (source of truth) | PASS |
| `pyproject.toml` | dynamic | `dynamic = ["version"]` via maturin | PASS |
| `Dockerfile` | 1.3.0 | `ARG VERSION=1.3.0` | PASS |
| `Singularity.def` | 1.3.0 | `From:` tag + `Version` label | PASS |
| `bioconda-recipe/meta.yaml` | 1.3.0 | `{% set version = "1.3.0" %}` | PASS |
| `CHANGELOG.md` | 1.3.0 | `## [1.3.0] - 2025-01-29` (manual) | PASS |

**Automated verification:** `scripts/check-version-consistency.sh` passes with exit code 0. Note: the script checks Dockerfile, Singularity.def, meta.yaml, and pyproject.toml against Cargo.toml. CHANGELOG.md was verified manually.

---

## Entry Points

| Console Script | Target | Module Exists | `app` Defined |
|----------------|--------|---------------|---------------|
| `wasp2-count` | `counting.__main__:app` | Yes | Yes (`typer.Typer`) |
| `wasp2-map` | `mapping.__main__:app` | Yes | Yes (`typer.Typer`) |
| `wasp2-analyze` | `analysis.__main__:app` | Yes | Yes (`typer.Typer`) |

All entry points reference valid modules under `src/` with correctly defined `app` objects.

---

## Optional Dependency Groups

| Group | Purpose | Package Count | Status |
|-------|---------|---------------|--------|
| `dev` | Testing, linting, type checking, security tools | 13 | PASS |
| `benchmark` | Performance profiling and visualization | 4 | PASS |
| `docs` | Sphinx documentation generation | 8 | PASS |
| `rust` | Rust extension building via maturin | 1 | PASS |
| `plink` | PLINK2 `.pgen` format support via Pgenlib | 1 | PASS |
| `cyvcf2` | Fast VCF parsing via cyvcf2 | 1 | PASS |

All groups declared in the issue (`cyvcf2`, `plink`, `docs`, `dev`, `benchmark`) are present. The `rust` group is an additional valid group for building the Rust extension.

---

## PyPI Publishing Configuration

| Setting | Value | Status |
|---------|-------|--------|
| Build backend | `maturin>=1.6,<2.0` | PASS |
| Bindings | `pyo3` | PASS |
| Python source | `src` | PASS |
| Python packages | `counting`, `mapping`, `analysis`, `wasp2` | PASS |
| Manifest path | `rust/Cargo.toml` | PASS |
| Strip binaries | `true` | PASS |
| Includes | `LICENSE`, `README.md` | PASS |
| Classifiers | Production/Stable, Python 3.10-3.12 | PASS |
| License | MIT | PASS |

---

## Findings

No issues found. All packaging metadata is consistent and correctly configured.

### Existing Safeguards

- `scripts/check-version-consistency.sh` automates version drift detection across key packaging files (Dockerfile, Singularity.def, meta.yaml, pyproject.toml)
- `pyproject.toml` uses `dynamic = ["version"]` to avoid manual synchronization with `Cargo.toml`
- Dockerfile and Singularity.def include comments pointing to `Cargo.toml` as the source of truth
