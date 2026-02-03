#!/usr/bin/env bash
# Check that version strings across all packaging files match rust/Cargo.toml
# (the single source of truth for WASP2 version).
#
# Usage: ./scripts/check-version-consistency.sh
# Exit code: 0 if all consistent, 1 if mismatches found

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"

# Extract version from Cargo.toml (single source of truth)
CARGO_VERSION=$(grep '^version' "$REPO_ROOT/rust/Cargo.toml" | head -1 | sed 's/.*"\(.*\)".*/\1/')

if [ -z "$CARGO_VERSION" ]; then
    echo "ERROR: Could not extract version from rust/Cargo.toml"
    exit 1
fi

echo "Source of truth (rust/Cargo.toml): $CARGO_VERSION"
echo "---"

ERRORS=0

check_version() {
    local file="$1"
    local found="$2"

    if [ "$found" = "$CARGO_VERSION" ]; then
        echo "OK   $file ($found)"
    else
        echo "FAIL $file: expected $CARGO_VERSION, found $found"
        ERRORS=$((ERRORS + 1))
    fi
}

# Dockerfile: ARG VERSION=x.y.z
# pyproject.toml should use dynamic version, not a hardcoded one
if grep -q '^version\s*=' "$REPO_ROOT/pyproject.toml"; then
    echo "FAIL pyproject.toml: contains hardcoded version (should use dynamic = [\"version\"])"
    ERRORS=$((ERRORS + 1))
else
    echo "OK   pyproject.toml (dynamic version via maturin)"
fi

DOCKER_VERSION=$(grep '^ARG VERSION=' "$REPO_ROOT/Dockerfile" | sed 's/ARG VERSION=//' || true)
check_version "Dockerfile" "$DOCKER_VERSION"

SING_FROM_VERSION=$(grep '^From:' "$REPO_ROOT/Singularity.def" | sed 's/.*://' || true)
check_version "Singularity.def (From)" "$SING_FROM_VERSION"

SING_LABEL_VERSION=$(grep '^\s*Version' "$REPO_ROOT/Singularity.def" | awk '{print $2}' || true)
check_version "Singularity.def (Label)" "$SING_LABEL_VERSION"

META_VERSION=$(grep 'set version' "$REPO_ROOT/bioconda-recipe/meta.yaml" | sed 's/.*"\(.*\)".*/\1/' || true)
check_version "bioconda-recipe/meta.yaml" "$META_VERSION"

echo "---"
if [ "$ERRORS" -gt 0 ]; then
    echo "FAILED: $ERRORS version mismatch(es) found"
    echo "Update all files to match rust/Cargo.toml version: $CARGO_VERSION"
    exit 1
else
    echo "All versions consistent: $CARGO_VERSION"
    exit 0
fi
