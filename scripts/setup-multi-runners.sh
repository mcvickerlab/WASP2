#!/bin/bash
# ==============================================================================
# WASP2 Multi-Runner Setup Script
# Sets up 3 specialized GitHub Actions runners on Mac M3 Max for parallelization
# Based on best practices from GenVarLoader, uv, pysam, and polars projects
# ==============================================================================

set -eo pipefail

# Configuration
REPO="Jaureguy760/WASP2-final"
RUNNERS_BASE="${HOME}/wasp2-runners"

# Runner configurations (3 specialized runners for M3 Max)
declare -A RUNNERS
RUNNERS["python"]="wasp2-python-runner:python,testing,lint,fast"
RUNNERS["rust"]="wasp2-rust-runner:rust,build,maturin"
RUNNERS["analysis"]="wasp2-analysis-runner:analysis,bioinformatics,docker,slow"

echo "=============================================="
echo "WASP2 Multi-Runner Setup (M3 Max Optimized)"
echo "=============================================="
echo ""
echo "This will set up 3 specialized runners:"
echo "  1. python-runner  - Fast Python tests, linting"
echo "  2. rust-runner    - Rust builds, maturin wheel building"
echo "  3. analysis-runner - Heavy analysis, Docker, slow tests"
echo ""

# ------------------------------------------------------------------------------
# Step 1: Check prerequisites
# ------------------------------------------------------------------------------
echo "[1/6] Checking prerequisites..."

# Check for gh CLI
if ! command -v gh &> /dev/null; then
    echo "  ❌ GitHub CLI not found. Install with: brew install gh"
    exit 1
fi
echo "  ✅ GitHub CLI installed"

# Check gh auth
if ! gh auth status &> /dev/null; then
    echo "  ❌ GitHub CLI not authenticated. Running 'gh auth login'..."
    gh auth login
fi
echo "  ✅ GitHub CLI authenticated"

# Check for Docker
if ! command -v docker &> /dev/null; then
    echo "  ⚠️  Docker not found. Analysis runner won't have Docker support."
else
    echo "  ✅ Docker installed"
fi

# Check for Python
if ! command -v python3 &> /dev/null; then
    echo "  ❌ Python3 not found. Install with: brew install python@3.11"
    exit 1
fi
echo "  ✅ Python3 installed ($(python3 --version))"

# Check for Rust
if ! command -v cargo &> /dev/null; then
    echo "  ⚠️  Rust not found. Installing..."
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
    source "$HOME/.cargo/env"
fi
echo "  ✅ Rust installed ($(cargo --version))"

# ------------------------------------------------------------------------------
# Step 2: Create runners base directory
# ------------------------------------------------------------------------------
echo ""
echo "[2/6] Creating runners directory structure..."
mkdir -p "$RUNNERS_BASE"
echo "  ✅ Created $RUNNERS_BASE"

# ------------------------------------------------------------------------------
# Step 3: Download runner binary (if not cached)
# ------------------------------------------------------------------------------
echo ""
echo "[3/6] Downloading GitHub Actions runner..."

RUNNER_VERSION=$(curl -s https://api.github.com/repos/actions/runner/releases/latest | grep '"tag_name":' | sed -E 's/.*"v([^"]+)".*/\1/')
if [[ -z "$RUNNER_VERSION" ]]; then
    echo "  ❌ Failed to determine latest runner version from GitHub API."
    echo "     This may be caused by network issues or GitHub API rate limiting."
    exit 1
fi
RUNNER_ARCHIVE="actions-runner-osx-arm64-${RUNNER_VERSION}.tar.gz"
RUNNER_CACHE="${RUNNERS_BASE}/.cache"

mkdir -p "$RUNNER_CACHE"

if [[ ! -f "${RUNNER_CACHE}/${RUNNER_ARCHIVE}" ]]; then
    echo "  ⏳ Downloading runner v${RUNNER_VERSION}..."
    curl -fSo "${RUNNER_CACHE}/${RUNNER_ARCHIVE}" -L \
        "https://github.com/actions/runner/releases/download/v${RUNNER_VERSION}/${RUNNER_ARCHIVE}"

    # Verify SHA256 checksum (supply chain protection)
    # GitHub embeds checksums in release notes body; extract via API
    echo "  ⏳ Verifying checksum..."
    EXPECTED_HASH=$(curl -sL https://api.github.com/repos/actions/runner/releases/latest \
        | grep -A1 "$RUNNER_ARCHIVE" \
        | grep -oE '[a-f0-9]{64}' \
        | head -1)
    if [[ -z "$EXPECTED_HASH" ]] || [[ ! "$EXPECTED_HASH" =~ ^[a-f0-9]{64}$ ]]; then
        echo "  ⚠️  Could not retrieve checksum from GitHub release notes."
        echo "     Skipping verification (manual verification recommended)."
    else
        ACTUAL_HASH=$(shasum -a 256 "${RUNNER_CACHE}/${RUNNER_ARCHIVE}" | awk '{print $1}')
        if [[ "$EXPECTED_HASH" != "$ACTUAL_HASH" ]]; then
            echo "  ❌ Checksum verification failed!"
            echo "     Expected: $EXPECTED_HASH"
            echo "     Actual:   $ACTUAL_HASH"
            rm -f "${RUNNER_CACHE}/${RUNNER_ARCHIVE}"
            exit 1
        fi
        echo "  ✅ Checksum verified"
    fi
    echo "  ✅ Downloaded runner"
else
    echo "  ✅ Runner v${RUNNER_VERSION} already cached"
fi

# ------------------------------------------------------------------------------
# Step 4: Set up each runner
# ------------------------------------------------------------------------------
echo ""
echo "[4/6] Setting up individual runners..."

for runner_type in "${!RUNNERS[@]}"; do
    IFS=':' read -r runner_name runner_labels <<< "${RUNNERS[$runner_type]}"
    runner_dir="${RUNNERS_BASE}/${runner_type}"

    echo ""
    echo "  Setting up: ${runner_name}"
    echo "    Labels: self-hosted, macOS, ARM64, ${runner_labels}"

    # Create runner directory
    mkdir -p "$runner_dir"

    # Extract runner if not already set up
    if [[ ! -f "${runner_dir}/config.sh" ]]; then
        echo "    ⏳ Extracting runner..."
        tar xzf "${RUNNER_CACHE}/${RUNNER_ARCHIVE}" -C "$runner_dir"
    fi

    # Get registration token (requires POST)
    echo "    ⏳ Getting registration token..."
    TOKEN=$(gh api -X POST "repos/${REPO}/actions/runners/registration-token" --jq '.token')

    if [[ -z "$TOKEN" ]]; then
        echo "    ❌ Failed to get registration token for ${runner_name}"
        continue
    fi

    # Remove old configuration if exists
    if [[ -f "${runner_dir}/.runner" ]]; then
        echo "    ⚠️  Removing old configuration..."
        cd "$runner_dir"
        ./config.sh remove --token "$TOKEN" 2>/dev/null || true
    fi

    # Configure runner
    echo "    ⏳ Configuring runner..."
    cd "$runner_dir"
    ./config.sh \
        --url "https://github.com/${REPO}" \
        --token "$TOKEN" \
        --name "$runner_name" \
        --labels "${runner_labels}" \
        --work "_work" \
        --replace \
        --unattended

    echo "    ✅ ${runner_name} configured"
done

# ------------------------------------------------------------------------------
# Step 5: Install as services
# ------------------------------------------------------------------------------
echo ""
echo "[5/6] Installing runners as services..."

for runner_type in "${!RUNNERS[@]}"; do
    IFS=':' read -r runner_name runner_labels <<< "${RUNNERS[$runner_type]}"
    runner_dir="${RUNNERS_BASE}/${runner_type}"

    echo "  Installing ${runner_name} service..."
    cd "$runner_dir"

    # Uninstall old service if exists
    ./svc.sh uninstall 2>/dev/null || true

    # Install as service
    ./svc.sh install
    echo "    ✅ ${runner_name} service installed"
done

# ------------------------------------------------------------------------------
# Step 6: Start all runners
# ------------------------------------------------------------------------------
echo ""
echo "[6/6] Starting all runners..."

for runner_type in "${!RUNNERS[@]}"; do
    IFS=':' read -r runner_name runner_labels <<< "${RUNNERS[$runner_type]}"
    runner_dir="${RUNNERS_BASE}/${runner_type}"

    echo "  Starting ${runner_name}..."
    cd "$runner_dir"
    ./svc.sh start
    echo "    ✅ ${runner_name} started"
done

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------
echo ""
echo "=============================================="
echo "✅ WASP2 Multi-Runner Setup Complete!"
echo "=============================================="
echo ""
echo "Runners configured:"
for runner_type in "${!RUNNERS[@]}"; do
    IFS=':' read -r runner_name runner_labels <<< "${RUNNERS[$runner_type]}"
    echo "  - ${runner_name}"
    echo "    Labels: self-hosted, macOS, ARM64, ${runner_labels}"
    echo "    Path: ${RUNNERS_BASE}/${runner_type}"
done
echo ""
echo "Management Commands:"
echo "  Check all:   for d in ${RUNNERS_BASE}/*/; do (cd \"\$d\" && ./svc.sh status); done"
echo "  Stop all:    for d in ${RUNNERS_BASE}/*/; do (cd \"\$d\" && ./svc.sh stop); done"
echo "  Start all:   for d in ${RUNNERS_BASE}/*/; do (cd \"\$d\" && ./svc.sh start); done"
echo ""
echo "Workflow labels to use:"
echo "  runs-on: [self-hosted, macOS, ARM64, python]     # Fast Python tests"
echo "  runs-on: [self-hosted, macOS, ARM64, rust]       # Rust builds"
echo "  runs-on: [self-hosted, macOS, ARM64, analysis]   # Heavy workloads"
echo "  runs-on: [self-hosted, macOS, ARM64, docker]     # Docker builds"
echo ""
echo "Verify on GitHub:"
echo "  https://github.com/${REPO}/settings/actions/runners"
echo ""
