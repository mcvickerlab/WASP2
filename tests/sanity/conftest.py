"""Fixtures and data helpers for sanity tests using real chr21 HG00731 data.

This module provides:
- Fixtures for loading sanity test data from a release tarball
- Helper to download sanity dataset from GitHub releases
- Markers for sanity tests
"""

from __future__ import annotations

import hashlib
import os
import shutil
import subprocess
import tarfile
from pathlib import Path
from typing import TYPE_CHECKING

import pytest

if TYPE_CHECKING:
    from collections.abc import Generator

# Sanity data version and paths
SANITY_DATA_VERSION = "v1"
SANITY_TARBALL_NAME = f"wasp2-sanity-chr21-{SANITY_DATA_VERSION}.tar.xz"
SANITY_DATA_DIR = Path(os.environ.get("WASP2_SANITY_DATA_DIR", Path(__file__).parent / "data"))

# No public archive is currently published. Use a reviewed local dataset or
# configure an archive URL and its independently verified SHA-256 explicitly.
GITHUB_RELEASE_URL = os.environ.get("WASP2_SANITY_DATA_URL")
ZENODO_DOI_URL: str | None = None

# Expected files in sanity dataset
SANITY_FILES = [
    "chr21.bam",
    "chr21.bam.bai",
    "chr21.vcf.gz",
    "chr21.vcf.gz.tbi",
    "expected_counts.tsv",
    "expected_r1.fq.gz",
    "expected_r2.fq.gz",
    "expected_analysis.tsv",
]


def pytest_configure(config: pytest.Config) -> None:
    """Register sanity test markers."""
    config.addinivalue_line(
        "markers",
        "sanity: marks tests as sanity tests using real chr21 data "
        "(deselect with '-m \"not sanity\"')",
    )


def is_sanity_data_available() -> bool:
    """Check if all sanity data files are present."""
    return all((SANITY_DATA_DIR / f).exists() for f in SANITY_FILES)


def download_sanity_data(
    release_url: str | None = None,
    force: bool = False,
) -> Path:
    """Download and extract sanity dataset from GitHub release or Zenodo.

    Data hosting strategy:
    - Primary: GitHub Releases (fast, CI-integrated, cached by Actions)
    - Fallback: Zenodo (DOI-backed archival for academic citation)

    Parameters
    ----------
    release_url : str | None
        URL to the release tarball. If None, attempts GitHub Releases first,
        then falls back to Zenodo if available.
    force : bool
        If True, re-download even if data exists.

    Returns
    -------
    Path
        Path to the extracted data directory.

    Raises
    ------
    RuntimeError
        If download or extraction fails from all sources.
    """
    if is_sanity_data_available() and not force:
        return SANITY_DATA_DIR

    SANITY_DATA_DIR.mkdir(parents=True, exist_ok=True)

    # Build list of URLs to try (primary first, then fallbacks)
    urls_to_try = []
    if release_url is not None:
        urls_to_try.append(release_url)
    else:
        if GITHUB_RELEASE_URL:
            urls_to_try.append(GITHUB_RELEASE_URL)
        if ZENODO_DOI_URL is not None:
            urls_to_try.append(ZENODO_DOI_URL)

    if not urls_to_try:
        raise RuntimeError(
            "No public sanity archive configured. Set WASP2_SANITY_DATA_DIR for local data, or WASP2_SANITY_DATA_URL and WASP2_SANITY_DATA_SHA256."
        )
    expected_hash = os.environ.get("WASP2_SANITY_DATA_SHA256")
    if not expected_hash or len(expected_hash) != 64:
        raise RuntimeError("Set WASP2_SANITY_DATA_SHA256 to the reviewed archive digest")
    tarball_path = SANITY_DATA_DIR / SANITY_TARBALL_NAME

    # Try each URL until one succeeds
    last_error = None
    for url in urls_to_try:
        try:
            if shutil.which("wget"):
                subprocess.run(
                    ["wget", "-q", "-O", str(tarball_path), url],
                    check=True,
                    capture_output=True,
                )
            elif shutil.which("curl"):
                subprocess.run(
                    ["curl", "--fail", "--show-error", "-sL", "-o", str(tarball_path), url],
                    check=True,
                    capture_output=True,
                )
            else:
                raise RuntimeError("Neither wget nor curl available for download")
            # Success - break out of loop
            break
        except subprocess.CalledProcessError as e:
            last_error = e
            tarball_path.unlink(missing_ok=True)  # Clean up partial download
            continue
    else:
        # All URLs failed
        raise RuntimeError(
            f"Failed to download sanity data from any source. Last error: {last_error}"
        )

    # Verify before opening the archive; only regular files/directories are allowed.
    if hashlib.sha256(tarball_path.read_bytes()).hexdigest() != expected_hash:
        tarball_path.unlink(missing_ok=True)
        raise RuntimeError("Sanity archive checksum mismatch")
    try:
        with tarfile.open(tarball_path, "r:xz") as tar:
            # Extract to data directory, stripping top-level dir
            for member in tar.getmembers():
                # Strip the top-level directory from paths
                parts = Path(member.name).parts
                if len(parts) > 1:
                    member.name = str(Path(*parts[1:]))
                    destination = (SANITY_DATA_DIR / member.name).resolve()
                    if not destination.is_relative_to(SANITY_DATA_DIR.resolve()):
                        raise RuntimeError("Sanity archive contains an unsafe path")
                    if member.isdir():
                        destination.mkdir(parents=True, exist_ok=True)
                    elif member.isfile():
                        destination.parent.mkdir(parents=True, exist_ok=True)
                        with tar.extractfile(member) as source, destination.open("wb") as target:
                            shutil.copyfileobj(source, target)
                    else:
                        raise RuntimeError("Sanity archive contains a non-regular file")
    except (tarfile.TarError, OSError) as e:
        raise RuntimeError(f"Failed to extract sanity data: {e}") from e
    finally:
        # Clean up tarball
        tarball_path.unlink(missing_ok=True)

    if not is_sanity_data_available():
        raise RuntimeError("Sanity archive is missing required expected outputs")
    return SANITY_DATA_DIR


@pytest.fixture(scope="session")
def sanity_data_dir() -> Path:
    """Return path to sanity test data directory.

    Skips test if data is not available.
    """
    if not is_sanity_data_available():
        pytest.skip(
            "Sanity data not available. Run 'make download-sanity-data' "
            "after configuring a reviewed archive URL and SHA-256."
        )
    return SANITY_DATA_DIR


@pytest.fixture(scope="session")
def sanity_data(sanity_data_dir: Path) -> dict[str, Path]:
    """Load sanity test data paths.

    Returns a dictionary with paths to all sanity data files:
    - bam: chr21.bam
    - bam_index: chr21.bam.bai
    - vcf: chr21.vcf.gz
    - vcf_index: chr21.vcf.gz.tbi
    - expected_counts: expected_counts.tsv
    - expected_r1: expected_r1.fq.gz
    - expected_r2: expected_r2.fq.gz
    - expected_analysis: expected_analysis.tsv
    """
    return {
        "bam": sanity_data_dir / "chr21.bam",
        "bam_index": sanity_data_dir / "chr21.bam.bai",
        "vcf": sanity_data_dir / "chr21.vcf.gz",
        "vcf_index": sanity_data_dir / "chr21.vcf.gz.tbi",
        "expected_counts": sanity_data_dir / "expected_counts.tsv",
        "expected_r1": sanity_data_dir / "expected_r1.fq.gz",
        "expected_r2": sanity_data_dir / "expected_r2.fq.gz",
        "expected_analysis": sanity_data_dir / "expected_analysis.tsv",
    }


@pytest.fixture(scope="function")
def sanity_tmp_dir(tmp_path: Path) -> Generator[Path, None, None]:
    """Provide a temporary directory for sanity test outputs."""
    output_dir = tmp_path / "sanity_output"
    output_dir.mkdir()
    yield output_dir
