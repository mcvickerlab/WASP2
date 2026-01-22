"""
Pytest configuration and shared fixtures for WASP2 tests.

This module provides:
- Test data fixtures (VCF, PGEN files)
- Temporary directory fixtures
- Mock objects for testing
"""

import gzip
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pytest

# Project root
ROOT = Path(__file__).parent.parent
TEST_DATA_DIR = ROOT / "tests" / "data"


# ============================================================================
# Session-scoped fixtures (created once per test session)
# ============================================================================

@pytest.fixture(scope="session")
def test_data_dir() -> Path:
    """Return path to test data directory, creating if needed."""
    TEST_DATA_DIR.mkdir(parents=True, exist_ok=True)
    return TEST_DATA_DIR


@pytest.fixture(scope="session")
def sample_vcf_content() -> str:
    """Generate minimal VCF content for testing."""
    return """\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2
chr1\t100\trs1\tA\tG\t30\tPASS\tDP=50\tGT\t0/1\t0/0
chr1\t200\trs2\tC\tT\t30\tPASS\tDP=45\tGT\t1/1\t0/1
chr1\t300\trs3\tG\tA\t30\tPASS\tDP=60\tGT\t0/0\t1/1
chr1\t400\trs4\tT\tC\t30\tPASS\tDP=55\tGT\t0/1\t0/1
chr2\t100\trs5\tA\tT\t30\tPASS\tDP=40\tGT\t0/1\t0/0
chr2\t200\trs6\tG\tC\t30\tPASS\tDP=35\tGT\t./.\t0/1
"""


@pytest.fixture(scope="session")
def sample_vcf(test_data_dir, sample_vcf_content) -> Path:
    """Create a sample VCF file for testing."""
    vcf_path = test_data_dir / "sample.vcf"
    vcf_path.write_text(sample_vcf_content)
    return vcf_path


@pytest.fixture(scope="session")
def sample_vcf_gz(test_data_dir, sample_vcf) -> Path:
    """Create a bgzipped and indexed VCF file for testing.

    Uses bcftools to properly bgzip the file (required for pysam/tabix).
    """
    vcf_gz_path = test_data_dir / "sample.vcf.gz"

    # Remove old file if exists (might be wrong format)
    if vcf_gz_path.exists():
        vcf_gz_path.unlink()
    tbi_path = Path(str(vcf_gz_path) + ".tbi")
    if tbi_path.exists():
        tbi_path.unlink()

    # Use bcftools to properly bgzip (required for pysam)
    try:
        subprocess.run(
            ["bcftools", "view", "-Oz", "-o", str(vcf_gz_path), str(sample_vcf)],
            check=True, capture_output=True
        )
        # Create tabix index
        subprocess.run(
            ["bcftools", "index", "-t", str(vcf_gz_path)],
            check=True, capture_output=True
        )
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        # Fall back to bgzip if bcftools fails
        try:
            subprocess.run(
                ["bgzip", "-c", str(sample_vcf)],
                stdout=open(vcf_gz_path, 'wb'),
                check=True
            )
            subprocess.run(
                ["tabix", "-p", "vcf", str(vcf_gz_path)],
                check=True, capture_output=True
            )
        except (subprocess.CalledProcessError, FileNotFoundError):
            pytest.skip(f"bcftools/bgzip not available for bgzip compression")

    return vcf_gz_path


@pytest.fixture(scope="session")
def sample_pgen_files(test_data_dir, sample_vcf) -> Dict[str, Path]:
    """Create sample PGEN/PVAR/PSAM files for testing.

    Returns dict with 'pgen', 'pvar', 'psam' keys.
    """
    pgen_prefix = test_data_dir / "sample"
    pgen_path = pgen_prefix.with_suffix('.pgen')
    pvar_path = pgen_prefix.with_suffix('.pvar')
    psam_path = pgen_prefix.with_suffix('.psam')

    # Try to convert VCF to PGEN using plink2
    try:
        subprocess.run([
            "plink2",
            "--vcf", str(sample_vcf),
            "--make-pgen",
            "--out", str(pgen_prefix),
            "--allow-extra-chr",
        ], check=True, capture_output=True)
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        pytest.skip(f"plink2 not available or conversion failed: {e}")

    return {
        'pgen': pgen_path,
        'pvar': pvar_path,
        'psam': psam_path,
        'prefix': pgen_prefix,
    }


# ============================================================================
# Function-scoped fixtures (created per test)
# ============================================================================

@pytest.fixture
def tmp_output_dir(tmp_path) -> Path:
    """Provide a temporary directory for test outputs."""
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    return output_dir


@pytest.fixture
def vcf_expected_variants() -> List[Dict]:
    """Expected variant data from sample VCF."""
    return [
        {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "G", "id": "rs1"},
        {"chrom": "chr1", "pos": 200, "ref": "C", "alt": "T", "id": "rs2"},
        {"chrom": "chr1", "pos": 300, "ref": "G", "alt": "A", "id": "rs3"},
        {"chrom": "chr1", "pos": 400, "ref": "T", "alt": "C", "id": "rs4"},
        {"chrom": "chr2", "pos": 100, "ref": "A", "alt": "T", "id": "rs5"},
        {"chrom": "chr2", "pos": 200, "ref": "G", "alt": "C", "id": "rs6"},
    ]


@pytest.fixture
def vcf_expected_het_sites_sample1() -> List[Dict]:
    """Expected heterozygous sites for sample1."""
    return [
        {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "G"},  # 0/1
        {"chrom": "chr1", "pos": 400, "ref": "T", "alt": "C"},  # 0/1
        {"chrom": "chr2", "pos": 100, "ref": "A", "alt": "T"},  # 0/1
    ]


@pytest.fixture
def vcf_expected_het_sites_sample2() -> List[Dict]:
    """Expected heterozygous sites for sample2."""
    return [
        {"chrom": "chr1", "pos": 200, "ref": "C", "alt": "T"},  # 0/1
        {"chrom": "chr1", "pos": 400, "ref": "T", "alt": "C"},  # 0/1
        {"chrom": "chr2", "pos": 200, "ref": "G", "alt": "C"},  # 0/1
    ]


# ============================================================================
# Markers
# ============================================================================

def pytest_configure(config):
    """Register custom markers not in pytest.ini."""
    # Note: slow, rust, integration, unit, benchmark markers are in pytest.ini
    config.addinivalue_line(
        "markers", "requires_plink2: marks tests that require plink2"
    )
    config.addinivalue_line(
        "markers", "requires_bcftools: marks tests that require bcftools"
    )


# ============================================================================
# Helper functions (not fixtures)
# ============================================================================

def has_command(cmd: str) -> bool:
    """Check if a command is available in PATH."""
    return shutil.which(cmd) is not None


def skip_without_plink2():
    """Skip test if plink2 is not available."""
    if not has_command("plink2"):
        pytest.skip("plink2 not available")


def skip_without_bcftools():
    """Skip test if bcftools is not available."""
    if not has_command("bcftools"):
        pytest.skip("bcftools not available")


def skip_without_pgenlib():
    """Skip test if pgenlib is not available."""
    try:
        import pgenlib
    except ImportError:
        pytest.skip("pgenlib not available")
