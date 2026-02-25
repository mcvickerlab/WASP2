# WASP2 Multi-stage Dockerfile
# Builds Rust extension and packages for Nextflow DSL2 modules

# ============================================================================
# Stage 1: Build Rust extension
# ============================================================================
FROM rust:1.87-bookworm AS rust-builder

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

# Install maturin
RUN pip3 install --break-system-packages --no-cache-dir maturin>=1.4

# Copy source files needed for maturin build
WORKDIR /build
COPY rust/ rust/
COPY src/ src/
COPY pyproject.toml .
COPY README.md .
COPY LICENSE .

# Build wheels
RUN maturin build --release -m rust/Cargo.toml -o /wheels

# ============================================================================
# Stage 2: Runtime image
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

# Install runtime dependencies + temporary build deps for pybedtools (C++ extension)
RUN apt-get update && apt-get install -y --no-install-recommends \
    # Bioinformatics tools
    samtools \
    bcftools \
    bedtools \
    tabix \
    # For htslib
    libhts3 \
    libbz2-1.0 \
    liblzma5 \
    zlib1g \
    libcurl4 \
    # Procps for ps command (Nextflow needs it)
    procps \
    # Build tools needed to compile pybedtools C++ extension
    g++ \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy wheel from builder and install
COPY --from=rust-builder /wheels/*.whl /tmp/
RUN pip install --no-cache-dir /tmp/*.whl && rm -rf /tmp/*.whl

# Remove build tools to reduce image size
RUN apt-get purge -y --auto-remove g++ zlib1g-dev && rm -rf /var/lib/apt/lists/*

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
