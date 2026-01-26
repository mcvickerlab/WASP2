#!/bin/bash
# ==============================================================================
# WASP2 Mac Runner Setup Script
# Sets up a self-hosted GitHub Actions runner on macOS with Docker support
# ==============================================================================

set -e

# Configuration
RUNNER_DIR="${HOME}/actions-runner"
REPO="Jaureguy760/WASP2-final"
RUNNER_NAME="wasp2-mac-runner"
RUNNER_LABELS="macOS,ARM64,docker,wasp2"

echo "=============================================="
echo "WASP2 Mac Runner Setup"
echo "=============================================="
echo ""

# ------------------------------------------------------------------------------
# Step 1: Check prerequisites
# ------------------------------------------------------------------------------
echo "[1/6] Checking prerequisites..."

# Check for Homebrew
if ! command -v brew &> /dev/null; then
    echo "  ❌ Homebrew not found. Installing..."
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
else
    echo "  ✅ Homebrew installed"
fi

# Check for Docker
if ! command -v docker &> /dev/null; then
    echo "  ❌ Docker not found. Please install Docker Desktop for Mac first:"
    echo "     https://www.docker.com/products/docker-desktop/"
    exit 1
else
    echo "  ✅ Docker installed ($(docker --version))"
fi

# Check Docker is running
if ! docker info &> /dev/null; then
    echo "  ❌ Docker is not running. Please start Docker Desktop."
    exit 1
else
    echo "  ✅ Docker is running"
fi

# Check for gh CLI
if ! command -v gh &> /dev/null; then
    echo "  ⏳ Installing GitHub CLI..."
    brew install gh
else
    echo "  ✅ GitHub CLI installed"
fi

# Check gh auth
if ! gh auth status &> /dev/null; then
    echo "  ❌ GitHub CLI not authenticated. Running 'gh auth login'..."
    gh auth login
fi
echo "  ✅ GitHub CLI authenticated"

# ------------------------------------------------------------------------------
# Step 2: Install development tools
# ------------------------------------------------------------------------------
echo ""
echo "[2/6] Installing development tools..."

# Python
if ! command -v python3 &> /dev/null; then
    echo "  ⏳ Installing Python..."
    brew install python@3.11
else
    echo "  ✅ Python installed ($(python3 --version))"
fi

# Rust
if ! command -v cargo &> /dev/null; then
    echo "  ⏳ Installing Rust..."
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
    source "$HOME/.cargo/env"
else
    echo "  ✅ Rust installed ($(cargo --version))"
fi

# Python tools
echo "  ⏳ Installing Python tools..."
pip3 install --quiet --upgrade pip
pip3 install --quiet maturin pytest pytest-cov ruff bandit mypy

echo "  ✅ Development tools installed"

# ------------------------------------------------------------------------------
# Step 3: Download GitHub Actions runner
# ------------------------------------------------------------------------------
echo ""
echo "[3/6] Setting up GitHub Actions runner..."

mkdir -p "$RUNNER_DIR"
cd "$RUNNER_DIR"

# Get latest runner version
RUNNER_VERSION=$(curl -s https://api.github.com/repos/actions/runner/releases/latest | grep '"tag_name":' | sed -E 's/.*"v([^"]+)".*/\1/')
echo "  Latest runner version: v${RUNNER_VERSION}"

# Download if not present or outdated
if [[ ! -f "$RUNNER_DIR/run.sh" ]]; then
    echo "  ⏳ Downloading runner..."
    curl -o actions-runner.tar.gz -L "https://github.com/actions/runner/releases/download/v${RUNNER_VERSION}/actions-runner-osx-arm64-${RUNNER_VERSION}.tar.gz"
    tar xzf actions-runner.tar.gz
    rm actions-runner.tar.gz
    echo "  ✅ Runner downloaded"
else
    echo "  ✅ Runner already exists"
fi

# ------------------------------------------------------------------------------
# Step 4: Get registration token and configure
# ------------------------------------------------------------------------------
echo ""
echo "[4/6] Configuring runner..."

# Get registration token
TOKEN=$(gh api -X POST "repos/${REPO}/actions/runners/registration-token" --jq '.token')

if [[ -z "$TOKEN" ]]; then
    echo "  ❌ Failed to get registration token. Check your GitHub permissions."
    exit 1
fi

# Check if already configured
if [[ -f "$RUNNER_DIR/.runner" ]]; then
    echo "  ⚠️  Runner already configured. Removing old configuration..."
    ./config.sh remove --token "$TOKEN" || true
fi

# Configure runner
echo "  ⏳ Registering runner with GitHub..."
./config.sh \
    --url "https://github.com/${REPO}" \
    --token "$TOKEN" \
    --name "$RUNNER_NAME" \
    --labels "$RUNNER_LABELS" \
    --work "_work" \
    --replace

echo "  ✅ Runner configured"

# ------------------------------------------------------------------------------
# Step 5: Install as service
# ------------------------------------------------------------------------------
echo ""
echo "[5/6] Installing as launch service..."

# Install service
./svc.sh install

echo "  ✅ Service installed"

# ------------------------------------------------------------------------------
# Step 6: Start the runner
# ------------------------------------------------------------------------------
echo ""
echo "[6/6] Starting runner..."

./svc.sh start

echo "  ✅ Runner started"

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------
echo ""
echo "=============================================="
echo "✅ WASP2 Mac Runner Setup Complete!"
echo "=============================================="
echo ""
echo "Runner Details:"
echo "  Name:     $RUNNER_NAME"
echo "  Labels:   self-hosted, $RUNNER_LABELS"
echo "  Repo:     $REPO"
echo "  Location: $RUNNER_DIR"
echo ""
echo "Management Commands:"
echo "  Check status:  cd $RUNNER_DIR && ./svc.sh status"
echo "  Stop runner:   cd $RUNNER_DIR && ./svc.sh stop"
echo "  Start runner:  cd $RUNNER_DIR && ./svc.sh start"
echo "  View logs:     cd $RUNNER_DIR && cat _diag/*.log"
echo ""
echo "Verify on GitHub:"
echo "  https://github.com/${REPO}/settings/actions/runners"
echo ""
