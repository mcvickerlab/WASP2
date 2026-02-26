# syntax=docker/dockerfile:1
# WASP2 Multi-stage Dockerfile
# Builds Rust extension and packages for Nextflow DSL2 modules
# Uses cargo-chef for Rust dependency caching and BuildKit cache mounts

# ============================================================================
# Stage 1: Install cargo-chef
# ============================================================================
FROM rust:1.88-bookworm AS chef
RUN cargo install cargo-chef
WORKDIR /build

# ============================================================================
# Stage 2: Plan — generate dependency recipe (only depends on Cargo.toml/lock)
# ============================================================================
FROM chef AS planner
WORKDIR /build/rust
COPY rust/Cargo.toml rust/Cargo.lock ./
COPY rust/src/ src/
RUN cargo chef prepare --recipe-path /build/recipe.json

# ============================================================================
# Stage 3: Cook + Build — compile deps (cached), then compile source
# ============================================================================
FROM chef AS rust-builder

# Install build dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3-dev \
    python3-pip \
    libclang-dev \
    libhts-dev \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    pkg-config \
    cmake \
    && rm -rf /var/lib/apt/lists/*

# Install maturin (before source copy so this layer is cached)
RUN pip3 install --break-system-packages --no-cache-dir "maturin>=1.4"

# Cook dependencies — this layer only rebuilds when Cargo.toml/lock changes
WORKDIR /build/rust
COPY --from=planner /build/recipe.json /build/recipe.json
COPY rust/Cargo.toml rust/Cargo.lock ./
RUN --mount=type=cache,target=/usr/local/cargo/registry \
    --mount=type=cache,target=/usr/local/cargo/git \
    --mount=type=cache,target=/build/rust/target \
    cargo chef cook --release --recipe-path /build/recipe.json

# Copy full source and build wheel
WORKDIR /build
COPY rust/ rust/
COPY src/ src/
COPY pyproject.toml .
COPY README.md .
COPY LICENSE .

RUN --mount=type=cache,target=/usr/local/cargo/registry \
    --mount=type=cache,target=/usr/local/cargo/git \
    --mount=type=cache,target=/build/rust/target \
    maturin build --release -m rust/Cargo.toml -o /wheels

# ============================================================================
# Stage 4: Runtime image
# ============================================================================
FROM python:3.11-slim-bookworm

# Version: keep in sync with rust/Cargo.toml (single source of truth)
# Run scripts/check-version-consistency.sh to verify
ARG VERSION=1.4.0

LABEL org.opencontainers.image.source="https://github.com/mcvickerlab/WASP2"
LABEL org.opencontainers.image.description="WASP2: Allele-specific analysis of NGS data with Rust acceleration"
LABEL org.opencontainers.image.licenses="MIT"
LABEL org.opencontainers.image.vendor="mcvickerlab"
LABEL org.opencontainers.image.title="WASP2"
LABEL org.opencontainers.image.version="${VERSION}"
LABEL maintainer="Jeff Jaureguy <jeffpjaureguy@gmail.com>"

# Install runtime deps + temporary build deps for pybedtools C++ extension
# Combined into one RUN to minimize layers; build tools purged at the end
RUN apt-get update && apt-get install -y --no-install-recommends \
    # Bioinformatics tools
    samtools \
    bcftools \
    bedtools \
    tabix \
    # htslib runtime libs
    libhts3 \
    libbz2-1.0 \
    liblzma5 \
    zlib1g \
    libcurl4 \
    # Procps for ps command (Nextflow needs it)
    procps \
    # Temporary: build tools for pybedtools C++ extension
    g++ \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy wheel from builder and install, then purge build tools in same layer
COPY --from=rust-builder /wheels/*.whl /tmp/
RUN --mount=type=cache,target=/root/.cache/pip \
    pip install /tmp/*.whl \
    && rm -rf /tmp/*.whl \
    && apt-get purge -y --auto-remove g++ zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Verify non-Python tools are available (Python tools skipped during build
# because Polars uses AVX2 instructions that fail under QEMU emulation
# on ARM64 CI runners building linux/amd64 images)
RUN samtools --version && bcftools --version && bedtools --version

# Create non-root user for security
RUN groupadd -g 1000 wasp2 && \
    useradd -u 1000 -g wasp2 -m -s /sbin/nologin wasp2 && \
    mkdir -p /data && chown wasp2:wasp2 /data

# Switch to non-root user
USER wasp2

# Bundle test data and smoke test for container validation (~300K)
COPY --chown=wasp2:wasp2 tests/shared_data/chr_test.fa \
     tests/shared_data/chr_test.fa.fai \
     tests/shared_data/variants.vcf \
     tests/shared_data/variants.vcf.gz \
     tests/shared_data/variants.vcf.gz.tbi \
     tests/shared_data/annotation.gtf \
     tests/shared_data/regions.bed \
     tests/shared_data/sample1.bam \
     tests/shared_data/sample1.bam.bai \
     /opt/wasp2/test-data/
COPY --chown=wasp2:wasp2 scripts/container_smoke_test.sh /opt/wasp2/scripts/

# Prevent Python from writing bytecode and ensure output is unbuffered
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# Set working directory for Nextflow
WORKDIR /data

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD wasp2-count --version || exit 1

# Default command
CMD ["wasp2-count", "--help"]
