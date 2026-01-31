# WASP2 Multi-stage Dockerfile
# Builds Rust extension and packages for Nextflow DSL2 modules

# ============================================================================
# Stage 1: Build Rust extension
# ============================================================================
FROM rust:1.87-bookworm AS rust-builder

# Install build dependencies
RUN apt-get update && apt-get install -y \
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
RUN pip3 install --break-system-packages maturin>=1.4

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

# Build arguments for versioning (can be overridden at build time)
ARG VERSION=1.2.0

LABEL org.opencontainers.image.source="https://github.com/Jaureguy760/WASP2-final"
LABEL org.opencontainers.image.description="WASP2: Allele-specific analysis of NGS data with Rust acceleration"
LABEL org.opencontainers.image.licenses="MIT"
LABEL org.opencontainers.image.vendor="Jaureguy760"
LABEL org.opencontainers.image.title="WASP2"
LABEL org.opencontainers.image.version="${VERSION}"
LABEL maintainer="Jeff Jaureguy <jeffpjaureguy@gmail.com>"

# Install runtime dependencies
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
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean

# Copy wheel from builder and install
COPY --from=rust-builder /wheels/*.whl /tmp/
RUN pip install --no-cache-dir /tmp/*.whl && rm -rf /tmp/*.whl

# Copy Python source and install (non-editable for production)
WORKDIR /app
COPY src/ src/
COPY pyproject.toml .
COPY README.md .
RUN pip install --no-cache-dir . --no-build-isolation

# Verify installation
RUN wasp2-count --help && \
    wasp2-map --help && \
    wasp2-analyze --help && \
    python -c "import wasp2_rust; print('Rust extension loaded successfully')" && \
    samtools --version && \
    bcftools --version

# Create non-root user for security
RUN groupadd -g 1000 wasp2 && \
    useradd -u 1000 -g wasp2 -m -s /bin/bash wasp2 && \
    mkdir -p /data && chown wasp2:wasp2 /data

# Switch to non-root user
USER wasp2

# Set working directory for Nextflow
WORKDIR /data

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD wasp2-count --version || exit 1

# Default command
CMD ["wasp2-count", "--help"]
