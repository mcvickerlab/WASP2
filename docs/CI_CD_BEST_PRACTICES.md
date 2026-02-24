# WASP2 CI/CD Best Practices Guide

> Based on analysis of GenVarLoader, pysam, rust-bio, polars, and uv projects

## Architecture Overview

### Runner Configuration

WASP2 uses **3 specialized self-hosted runners** on Mac M3 Max for optimal parallelization:

| Runner | Labels | Purpose |
|--------|--------|---------|
| `wasp2-python-runner` | `python, testing, lint, fast` | Fast Python tests, linting |
| `wasp2-rust-runner` | `rust, build, maturin` | Rust builds, wheel building |
| `wasp2-analysis-runner` | `analysis, bioinformatics, docker, slow` | Heavy analysis, Docker |

### Workflow Selection

```yaml
# Fast Python tasks
runs-on: [self-hosted, macOS, ARM64, python]

# Rust compilation
runs-on: [self-hosted, macOS, ARM64, rust]

# Heavy analysis/Docker
runs-on: [self-hosted, macOS, ARM64, analysis]
```

## Workflow Schedule

All times in UTC. Staggered to prevent resource contention.

| Workflow | Schedule | Purpose |
|----------|----------|---------|
| Security | Monday 2:05am | Weekly security scans |
| CodeQL | Monday 2:15am | Static analysis |
| Nightly | Daily 3:00am | Extended tests, benchmarks |
| Dependabot | Monday 3:00am | Dependency updates |

## Caching Strategy

### Rust (using Swatinem/rust-cache)

```yaml
- uses: Swatinem/rust-cache@v2
  with:
    workspaces: rust
    save-if: github.ref == 'refs/heads/main'
```

**Benefits:**
- 50-70% faster incremental builds
- Only saves cache on main branch (prevents cache pollution)

### Python

```yaml
- uses: actions/cache@v4
  with:
    path: ~/.cache/pip
    key: pip-${{ runner.os }}-${{ hashFiles('pyproject.toml') }}
    restore-keys: |
      pip-${{ runner.os }}-
```

### Docker (BuildKit)

```yaml
- uses: docker/build-push-action@v6
  with:
    cache-from: type=gha,scope=wasp2
    cache-to: type=gha,mode=max,scope=wasp2
```

## Matrix Testing Strategy

### Python Version Matrix

```yaml
strategy:
  fail-fast: false
  matrix:
    python-version: ['3.10', '3.11', '3.12']
```

### Platform Matrix (Release)

```yaml
matrix:
  include:
    - os: ubuntu-latest
      target: x86_64
    - os: ubuntu-latest
      target: aarch64
    - os: macos-13
      target: x86_64-apple-darwin
    - os: macos-14
      target: aarch64-apple-darwin
```

## Concurrency Control

Prevent parallel runs on same branch:

```yaml
concurrency:
  group: ${{ github.workflow }}-${{ github.ref }}
  cancel-in-progress: true
```

## Dependabot Automation

### Auto-approve + Auto-merge for Patches

The `dependabot-auto-merge.yml` workflow:
- **Auto-approves** patch and minor updates
- **Auto-merges** patch updates only (safer)
- **Labels** major updates for manual review

### Dependabot Configuration

```yaml
# .github/dependabot.yml
version: 2
updates:
  - package-ecosystem: "pip"
    schedule:
      interval: "weekly"
      day: "monday"
      time: "03:00"
    commit-message:
      prefix: "chore(deps)"
```

## Security Scanning

### Tools

| Tool | Language | Purpose |
|------|----------|---------|
| `pip-audit` | Python | Dependency vulnerabilities |
| `bandit` | Python | Code security issues |
| `cargo-audit` | Rust | Dependency vulnerabilities |
| `gitleaks` | All | Secret detection |
| `CodeQL` | All | Static analysis |

### Blocking vs Informational

- **Blocking:** gitleaks (secrets must never be committed)
- **Informational:** pip-audit, bandit, cargo-audit (logged but don't fail PR)

## Release Workflow

### Multi-Platform Wheel Building

Using `PyO3/maturin-action@v1`:

```yaml
- uses: PyO3/maturin-action@v1
  with:
    target: ${{ matrix.target }}
    args: --release --out dist -m rust/Cargo.toml
    manylinux: "2014"
    before-script-linux: |
      yum install -y bzip2-devel xz-devel zlib-devel
```

### Trusted Publishing (OIDC)

No PyPI tokens needed! Configure in PyPI settings:

```yaml
permissions:
  id-token: write  # OIDC

- uses: pypa/gh-action-pypi-publish@release/v1
  # No token required - uses OIDC
```

## Performance Optimization Tips

### M3 Max Specific

1. **6-8 parallel jobs** maximum (leaving cores for OS)
2. **Separate runners** for different workload types
3. **sccache** for Rust incremental compilation
4. **Cache benchmark data** to avoid regeneration

### Expected Performance Gains

| Stage | Without Optimization | With Optimization |
|-------|---------------------|-------------------|
| Lint/Format | 2-3 min | 30-45 sec |
| Unit Tests | 5-8 min | 2-3 min |
| Rust Build (cold) | 8-12 min | 4-6 min |
| Rust Build (warm) | 2-3 min | 30-60 sec |
| Full Pipeline | 35-45 min | 8-12 min |

## Nightly Testing

The `nightly.yml` workflow runs:
- Extended unit tests (all Python versions)
- Integration tests
- Performance benchmarks
- Nextflow pipeline tests (optional)

### Running Manually

```bash
gh workflow run nightly.yml -f test_set=benchmarks
```

## Runner Management

### Setup Multi-Runners

```bash
./scripts/setup-multi-runners.sh
```

### Management Commands

```bash
# Check all runners
for d in ~/wasp2-runners/*/; do (cd "$d" && ./svc.sh status); done

# Stop all runners
for d in ~/wasp2-runners/*/; do (cd "$d" && ./svc.sh stop); done

# Start all runners
for d in ~/wasp2-runners/*/; do (cd "$d" && ./svc.sh start); done
```

### View on GitHub

https://github.com/mcvickerlab/WASP2/settings/actions/runners

## References

- [GenVarLoader CI/CD](https://github.com/mcvickerlab/GenVarLoader) - Pixi, commitizen, multi-platform
- [uv workflows](https://github.com/astral-sh/uv) - Planning jobs, caching
- [pysam CI](https://github.com/pysam-developers/pysam) - maturin, multi-platform wheels
- [rust-bio CI](https://github.com/rust-bio/rust-bio) - Rust best practices
- [GitHub Actions Best Practices](https://docs.github.com/en/actions/security-for-github-actions/security-guides/security-hardening-for-github-actions)
