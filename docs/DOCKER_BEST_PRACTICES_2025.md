# Docker Best Practices for Python Applications (2025-2026)

**Research Date:** January 2025
**Target Project:** WASP2 (Python + Rust bioinformatics)
**Author:** Performance Analysis

---

## Table of Contents

1. [Base Image Selection](#1-base-image-selection)
2. [Package Manager Comparison](#2-package-manager-comparison)
3. [Multi-Stage Build Patterns](#3-multi-stage-build-patterns)
4. [Image Size Reduction](#4-image-size-reduction)
5. [Build Cache Optimization](#5-build-cache-optimization)
6. [Runtime Performance](#6-runtime-performance)
7. [Security Hardening](#7-security-hardening)
8. [WASP2 Specific Recommendations](#8-wasp2-specific-recommendations)

---

## 1. Base Image Selection

### Comparison Matrix

| Image | Size | glibc | Shells | CVEs* | Scientific Computing | Bioinformatics |
|-------|------|-------|--------|-------|---------------------|----------------|
| `python:3.11-slim-bookworm` | 221MB | Yes | bash | Medium | Excellent | Excellent |
| `python:3.11-alpine` | 89MB | No (musl) | sh | Low | Poor | Poor |
| `gcr.io/distroless/python3-debian12` | 95MB | Yes | None | Very Low | Limited | Poor |
| `cgr.dev/chainguard/python` | ~80MB | Yes (Wolfi) | None | Near Zero | Good | Limited |
| `ghcr.io/astral-sh/uv:python3.11-bookworm-slim` | 288MB | Yes | bash | Medium | Excellent | Excellent |

*CVEs = Common Vulnerabilities and Exposures (typical count)

### Recommendations by Use Case

**Scientific/Bioinformatics (WASP2):**
- **Use:** `python:3.11-slim-bookworm`
- **Reason:** glibc compatibility for numpy/scipy wheels, apt for samtools/bcftools

**Minimal API Services:**
- **Use:** `gcr.io/distroless/python3-debian12`
- **Reason:** Smallest attack surface, no shell access

**Enterprise with SLA requirements:**
- **Use:** `cgr.dev/chainguard/python` (paid)
- **Reason:** Daily CVE patches, signed builds, SBOM included

**Development/CI:**
- **Use:** `ghcr.io/astral-sh/uv:python3.11-bookworm-slim`
- **Reason:** UV pre-installed, fast dependency resolution

### Why NOT Alpine for Scientific Python

Alpine uses musl libc instead of glibc, causing:

1. **No pre-built wheels:** numpy, scipy, pandas require compilation
2. **Build time:** 10-30 minutes vs 30 seconds with glibc wheels
3. **Binary compatibility:** Some C extensions fail or have subtle bugs
4. **Performance:** musl can be slower for numerical workloads (5-15%)

```dockerfile
# BAD for scientific Python
FROM python:3.11-alpine
RUN pip install numpy  # Compiles from source, may fail

# GOOD for scientific Python
FROM python:3.11-slim-bookworm
RUN pip install numpy  # Uses pre-built wheel, instant
```

---

## 2. Package Manager Comparison

### UV vs pip vs Poetry (2025)

| Feature | UV (v0.9.26) | pip (24.x) | Poetry (1.8.x) |
|---------|--------------|------------|----------------|
| **Install Speed** | 10-100x faster | Baseline | 0.3-0.5x |
| **Lockfile** | `uv.lock` (native) | None (pip-tools) | `poetry.lock` |
| **Resolution** | Parallel, incremental | Sequential | Sequential |
| **Docker Integration** | Excellent | Basic | Limited |
| **Cache Efficiency** | Excellent | Good | Moderate |
| **Binary Size** | 10MB | Built-in | ~50MB |
| **Rust Builds** | Native (maturin) | Manual | Manual |

### UV Docker Features

```dockerfile
# 1. Copy UV binary (10MB, pinned version)
COPY --from=ghcr.io/astral-sh/uv:0.9.26 /uv /uvx /bin/

# 2. Cache mounts for 10-100x faster rebuilds
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --locked

# 3. Bytecode compilation for 20-30% faster startup
ENV UV_COMPILE_BYTECODE=1

# 4. Non-editable install for production
RUN uv sync --locked --no-editable

# 5. Skip dev dependencies
ENV UV_NO_DEV=1
```

### Migration from pip to UV

```dockerfile
# BEFORE (pip)
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
COPY . .
RUN pip install --no-cache-dir .

# AFTER (UV) - 10-100x faster
COPY --from=ghcr.io/astral-sh/uv:0.9.26 /uv /bin/uv
COPY pyproject.toml uv.lock ./
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --locked --no-install-project
COPY . .
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --locked --no-editable
```

---

## 3. Multi-Stage Build Patterns

### Pattern for Python + Rust Projects (WASP2)

```dockerfile
# syntax=docker/dockerfile:1.7

# ============================================================
# Stage 1: Build Rust extension
# ============================================================
FROM rust:1.75-bookworm AS rust-builder

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3-dev libclang-dev libhts-dev libbz2-dev liblzma-dev \
    zlib1g-dev pkg-config && rm -rf /var/lib/apt/lists/*

RUN pip3 install --break-system-packages maturin>=1.4

WORKDIR /build
COPY rust/ rust/
COPY pyproject.toml README.md LICENSE ./

RUN --mount=type=cache,target=/root/.cargo/registry \
    --mount=type=cache,target=/root/.cargo/git \
    maturin build --release -m rust/Cargo.toml -o /wheels

# ============================================================
# Stage 2: Build Python dependencies
# ============================================================
FROM python:3.11-slim-bookworm AS python-builder

COPY --from=ghcr.io/astral-sh/uv:0.9.26 /uv /bin/uv
ENV UV_COMPILE_BYTECODE=1 UV_LINK_MODE=copy UV_NO_DEV=1

WORKDIR /app

# Install dependencies (cached layer)
COPY pyproject.toml uv.lock ./
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --locked --no-install-project

# Install Rust wheel
COPY --from=rust-builder /wheels/*.whl /tmp/
RUN --mount=type=cache,target=/root/.cache/uv \
    uv pip install /tmp/*.whl

# Install project (non-editable)
COPY src/ ./src/
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --locked --no-editable

# ============================================================
# Stage 3: Minimal runtime
# ============================================================
FROM python:3.11-slim-bookworm AS runtime

# Runtime dependencies only
RUN apt-get update && apt-get install -y --no-install-recommends \
    samtools bcftools bedtools tabix libhts3 procps \
    && rm -rf /var/lib/apt/lists/*

# Non-root user
RUN groupadd -g 1000 app && useradd -u 1000 -g app -m app
USER app

# Copy venv from builder
COPY --from=python-builder --chown=app:app /app/.venv /app/.venv
ENV PATH="/app/.venv/bin:$PATH"

WORKDIR /data
CMD ["python", "--version"]
```

### Key Benefits

| Aspect | Improvement |
|--------|-------------|
| Final image size | 50-70% smaller (no Rust toolchain) |
| Build cache hits | 80-95% on dependency changes |
| Security surface | Minimal (no compilers in runtime) |
| Rebuild time | 10-100x faster with UV caches |

---

## 4. Image Size Reduction

### Techniques Ranked by Impact

| Technique | Typical Savings | Implementation |
|-----------|-----------------|----------------|
| Multi-stage builds | 1-2GB | Separate build/runtime stages |
| Remove build tools | 500MB-1GB | Don't install gcc, make, etc. in runtime |
| Use slim base | 300-500MB | `-slim` variant instead of full |
| UV vs pip | 50-100MB | Better dependency resolution |
| `--no-cache-dir` | 50-200MB | `pip install --no-cache-dir` |
| Combine RUN layers | 20-50MB | Chain commands with `&&` |
| `.dockerignore` | Variable | Exclude tests, docs, .git |
| Strip binaries | 10-30MB | `strip` on compiled extensions |
| Remove apt lists | 30-50MB | `rm -rf /var/lib/apt/lists/*` |

### .dockerignore Best Practices

```gitignore
# Version control
.git
.gitignore

# Development
tests/
docs/
*.md
!README.md
.vscode/
.idea/

# Python artifacts
__pycache__/
*.pyc
.venv/
dist/
build/
*.egg-info/

# Rust artifacts (if not using multi-stage properly)
rust/target/

# Security sensitive
.env
.env.*
*.pem
*.key

# Large data files
*.bam
*.vcf
*.fastq
```

---

## 5. Build Cache Optimization

### BuildKit Cache Mounts

```dockerfile
# Python package cache
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --locked

# Rust crates cache
RUN --mount=type=cache,target=/root/.cargo/registry \
    --mount=type=cache,target=/root/.cargo/git \
    cargo build --release

# apt cache (Debian/Ubuntu)
RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    apt-get update && apt-get install -y package
```

### Layer Ordering Strategy

```dockerfile
# GOOD: Least-changing files first
COPY pyproject.toml uv.lock ./          # 1. Rarely changes
RUN uv sync --locked --no-install-project  # 2. Cached if deps unchanged
COPY src/ ./src/                         # 3. Changes frequently
RUN uv sync --locked                     # 4. Only rebuilds app

# BAD: Invalidates cache on any change
COPY . .
RUN uv sync --locked
```

### Bind Mounts for Config Files

```dockerfile
# Avoid COPY cache invalidation for config-only operations
RUN --mount=type=bind,source=pyproject.toml,target=pyproject.toml \
    --mount=type=bind,source=uv.lock,target=uv.lock \
    uv sync --locked --no-install-project
```

---

## 6. Runtime Performance

### Environment Variables

```dockerfile
# Bytecode compilation (20-30% faster startup)
ENV UV_COMPILE_BYTECODE=1

# Prevent runtime .pyc generation
ENV PYTHONDONTWRITEBYTECODE=1

# Unbuffered output for logging
ENV PYTHONUNBUFFERED=1

# Memory optimization for long-running processes
ENV MALLOC_ARENA_MAX=2

# Disable hash randomization (reproducible, slightly faster)
ENV PYTHONHASHSEED=0
```

### Startup Time Optimization

| Technique | Impact | Trade-off |
|-----------|--------|-----------|
| Bytecode compilation | -30% startup | +10% image size |
| Lazy imports | -50% startup | Code complexity |
| Module preloading | -20% startup | Memory usage |
| Static linking | -10% startup | Larger binaries |

### Memory Optimization

```dockerfile
# For memory-constrained environments
ENV MALLOC_ARENA_MAX=2
ENV PYTHONMALLOC=malloc

# For numpy/scipy
ENV OMP_NUM_THREADS=4
ENV MKL_NUM_THREADS=4
```

---

## 7. Security Hardening

### Non-Root User

```dockerfile
# Create user with specific UID/GID
RUN groupadd --system --gid 1000 app && \
    useradd --system --gid 1000 --uid 1000 \
    --create-home --shell /sbin/nologin app

# Switch before running application
USER app

# Verify permissions
RUN whoami && id
```

### Read-Only Filesystem

```bash
# Runtime flag (not in Dockerfile)
docker run --read-only --tmpfs /tmp --tmpfs /app/cache myimage
```

### Capability Dropping

```bash
# Runtime flag (not in Dockerfile)
docker run --cap-drop=ALL --cap-add=NET_BIND_SERVICE myimage
```

### Security Scanning

```bash
# Trivy (recommended)
trivy image --severity HIGH,CRITICAL wasp2:latest

# Grype
grype wasp2:latest

# Snyk
snyk container test wasp2:latest
```

### Image Signing and Verification

```bash
# Sign with cosign
cosign sign --key cosign.key wasp2:latest

# Verify
cosign verify --key cosign.pub wasp2:latest
```

### Security Checklist

- [ ] Non-root user (`USER` directive)
- [ ] No secrets in image (use runtime secrets)
- [ ] Base image pinned (digest or specific tag)
- [ ] Regular vulnerability scanning
- [ ] HEALTHCHECK defined
- [ ] Minimal installed packages
- [ ] No shell in production (consider distroless)
- [ ] Read-only filesystem where possible
- [ ] Capabilities dropped

---

## 8. WASP2 Specific Recommendations

### Summary

| Aspect | Current | Recommended |
|--------|---------|-------------|
| Base image | `python:3.11-slim-bookworm` | Keep (optimal for bioinformatics) |
| Package manager | pip | **UV** (10-100x faster) |
| Build stages | 2 | 3 (separate Python builder) |
| Cache mounts | None | **Add BuildKit mounts** |
| Non-root user | Yes | Keep |
| Bytecode compilation | No | **Enable** |
| Image signing | No | Consider for releases |

### Migration Path

1. **Phase 1 (Quick Win):** Add BuildKit cache mounts to existing Dockerfile
2. **Phase 2 (Medium):** Replace pip with UV, add bytecode compilation
3. **Phase 3 (Full):** Use optimized 3-stage Dockerfile

### Expected Improvements

| Metric | Current (est.) | Optimized | Improvement |
|--------|----------------|-----------|-------------|
| Image size | 800-1000MB | 500-600MB | 40% smaller |
| Fresh build | 10-15 min | 5-8 min | 50% faster |
| Cached rebuild | 3-5 min | 30-60 sec | 80% faster |
| Startup time | 2-3 sec | 1-2 sec | 40% faster |
| CVE count | Varies | Baseline-10% | Reduced |

### Files Created

- `/Users/jeffjaureguy/Projects/WASP2-final/Dockerfile.optimized` - New optimized Dockerfile
- `/Users/jeffjaureguy/Projects/WASP2-final/scripts/benchmark_docker_build.sh` - Build comparison script
- `/Users/jeffjaureguy/Projects/WASP2-final/.dockerignore` - Updated ignore patterns

---

## References

- [UV Docker Guide](https://docs.astral.sh/uv/guides/integration/docker/) - Official UV documentation
- [uv-docker-example](https://github.com/astral-sh/uv-docker-example) - Official example repository
- [Google Distroless](https://github.com/GoogleContainerTools/distroless) - Minimal base images
- [Chainguard Images](https://edu.chainguard.dev/chainguard/chainguard-images/getting-started/python/) - Zero-CVE images
- [Docker BuildKit](https://docs.docker.com/build/buildkit/) - Advanced build features
- [Maturin Docker](https://github.com/PyO3/maturin) - Rust/Python wheel building

---

*Document generated: January 2025*
*UV Version: 0.9.26*
*Python Target: 3.11*
