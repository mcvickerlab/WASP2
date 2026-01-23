"""
Pytest configuration and fixtures for WASP2 performance benchmarks.

This module provides:
- Synthetic data generation at various scales (1K, 10K, 100K, 1M variants)
- Memory profiling utilities
- Benchmark result collection and comparison
"""

import gc
import subprocess
from collections.abc import Generator
from pathlib import Path
from typing import Any

import numpy as np
import pytest

# Try to import memory_profiler for memory benchmarks
try:
    from memory_profiler import memory_usage

    HAS_MEMORY_PROFILER = True
except ImportError:
    HAS_MEMORY_PROFILER = False

# ============================================================================
# Benchmark configuration
# ============================================================================

# Scale levels for parametrized benchmarks
BENCHMARK_SCALES = {
    "tiny": 100,
    "small": 1_000,
    "medium": 10_000,
    "large": 100_000,
    "xlarge": 1_000_000,
}

# Default chromosomes for synthetic data
CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


# ============================================================================
# Benchmark fixtures
# ============================================================================


@pytest.fixture(scope="session")
def benchmark_data_dir(tmp_path_factory) -> Path:
    """Session-scoped temporary directory for benchmark data."""
    return tmp_path_factory.mktemp("benchmark_data")


@pytest.fixture(scope="session")
def rng() -> np.random.Generator:
    """Seeded random number generator for reproducible benchmarks."""
    return np.random.default_rng(seed=42)


# ============================================================================
# VCF data generation
# ============================================================================


def generate_synthetic_vcf(
    output_path: Path,
    n_variants: int,
    n_samples: int,
    rng: np.random.Generator,
    chromosomes: list[str] | None = None,
) -> Path:
    """Generate a synthetic VCF file with specified parameters."""
    if chromosomes is None:
        chromosomes = CHROMOSOMES

    bases = ["A", "C", "G", "T"]
    sample_names = [f"sample{i:04d}" for i in range(n_samples)]
    variants_per_chrom = n_variants // len(chromosomes)
    remainder = n_variants % len(chromosomes)

    with open(output_path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n')
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n')
        for chrom in chromosomes:
            f.write(f"##contig=<ID={chrom},length=250000000>\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
        for sample in sample_names:
            f.write(f"\t{sample}")
        f.write("\n")

        variant_id = 0
        for chrom_idx, chrom in enumerate(chromosomes):
            n_chrom_variants = variants_per_chrom + (1 if chrom_idx < remainder else 0)
            positions = np.sort(rng.integers(1, 249_000_000, size=n_chrom_variants))

            for pos in positions:
                variant_id += 1
                ref = rng.choice(bases)
                alt = rng.choice([b for b in bases if b != ref])
                gt_probs = [0.4, 0.35, 0.2, 0.05]
                gt_choices = ["0/0", "0/1", "1/1", "./."]
                genotypes = rng.choice(gt_choices, size=n_samples, p=gt_probs)
                depths = rng.integers(10, 100, size=n_samples)

                f.write(f"{chrom}\t{pos}\trs{variant_id}\t{ref}\t{alt}\t30\tPASS\tDP=50\tGT:DP")
                for gt, dp in zip(genotypes, depths, strict=False):
                    f.write(f"\t{gt}:{dp}")
                f.write("\n")

    return output_path


@pytest.fixture(scope="session")
def vcf_tiny(benchmark_data_dir: Path, rng: np.random.Generator) -> Path:
    """Generate tiny VCF (100 variants, 2 samples)."""
    return generate_synthetic_vcf(
        benchmark_data_dir / "tiny.vcf",
        n_variants=BENCHMARK_SCALES["tiny"],
        n_samples=2,
        rng=rng,
    )


@pytest.fixture(scope="session")
def vcf_small(benchmark_data_dir: Path, rng: np.random.Generator) -> Path:
    """Generate small VCF (1K variants, 10 samples)."""
    return generate_synthetic_vcf(
        benchmark_data_dir / "small.vcf",
        n_variants=BENCHMARK_SCALES["small"],
        n_samples=10,
        rng=rng,
    )


@pytest.fixture(scope="session")
def vcf_medium(benchmark_data_dir: Path, rng: np.random.Generator) -> Path:
    """Generate medium VCF (10K variants, 50 samples)."""
    return generate_synthetic_vcf(
        benchmark_data_dir / "medium.vcf",
        n_variants=BENCHMARK_SCALES["medium"],
        n_samples=50,
        rng=rng,
    )


@pytest.fixture(scope="session")
def vcf_large(benchmark_data_dir: Path, rng: np.random.Generator) -> Path:
    """Generate large VCF (100K variants, 100 samples)."""
    return generate_synthetic_vcf(
        benchmark_data_dir / "large.vcf",
        n_variants=BENCHMARK_SCALES["large"],
        n_samples=100,
        rng=rng,
    )


# ============================================================================
# Bgzipped/indexed VCF generation
# ============================================================================


def bgzip_and_index_vcf(vcf_path: Path) -> Path | None:
    """Bgzip and tabix-index a VCF file."""
    vcf_gz_path = Path(str(vcf_path) + ".gz")

    try:
        subprocess.run(
            ["bcftools", "view", "-Oz", "-o", str(vcf_gz_path), str(vcf_path)],
            check=True,
            capture_output=True,
        )
        subprocess.run(
            ["bcftools", "index", "-t", str(vcf_gz_path)],
            check=True,
            capture_output=True,
        )
        return vcf_gz_path
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None


@pytest.fixture(scope="session")
def vcf_small_gz(vcf_small: Path) -> Path | None:
    """Bgzipped small VCF."""
    return bgzip_and_index_vcf(vcf_small)


@pytest.fixture(scope="session")
def vcf_medium_gz(vcf_medium: Path) -> Path | None:
    """Bgzipped medium VCF."""
    return bgzip_and_index_vcf(vcf_medium)


# ============================================================================
# PGEN file generation
# ============================================================================


def vcf_to_pgen(vcf_path: Path, output_prefix: Path) -> dict[str, Path] | None:
    """Convert VCF to PGEN format using plink2."""
    try:
        subprocess.run(
            [
                "plink2",
                "--vcf",
                str(vcf_path),
                "--make-pgen",
                "--out",
                str(output_prefix),
                "--allow-extra-chr",
            ],
            check=True,
            capture_output=True,
        )
        return {
            "pgen": output_prefix.with_suffix(".pgen"),
            "pvar": output_prefix.with_suffix(".pvar"),
            "psam": output_prefix.with_suffix(".psam"),
            "prefix": output_prefix,
        }
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None


@pytest.fixture(scope="session")
def pgen_small(vcf_small: Path, benchmark_data_dir: Path) -> dict[str, Path] | None:
    """Convert small VCF to PGEN."""
    return vcf_to_pgen(vcf_small, benchmark_data_dir / "small_pgen")


@pytest.fixture(scope="session")
def pgen_medium(vcf_medium: Path, benchmark_data_dir: Path) -> dict[str, Path] | None:
    """Convert medium VCF to PGEN."""
    return vcf_to_pgen(vcf_medium, benchmark_data_dir / "medium_pgen")


# ============================================================================
# BAM file generation
# ============================================================================


def generate_synthetic_bam(
    output_path: Path,
    n_reads: int,
    rng: np.random.Generator,
    reference_length: int = 10_000_000,
) -> Path | None:
    """Generate a synthetic BAM file with specified parameters."""
    sam_path = output_path.with_suffix(".sam")
    bases = ["A", "C", "G", "T"]

    with open(sam_path, "w") as f:
        f.write("@HD\tVN:1.6\tSO:coordinate\n")
        f.write(f"@SQ\tSN:chr1\tLN:{reference_length}\n")
        f.write("@RG\tID:benchmark\tSM:sample1\n")

        read_length = 150
        for i in range(n_reads):
            pos = rng.integers(1, reference_length - read_length)
            seq = "".join(rng.choice(bases, size=read_length))
            qual = "".join(["I"] * read_length)
            flag = 99 if i % 2 == 0 else 147
            f.write(f"read{i:08d}\t{flag}\tchr1\t{pos}\t60\t{read_length}M\t=\t")
            f.write(f"{pos + 200}\t350\t{seq}\t{qual}\tRG:Z:benchmark\n")

    try:
        subprocess.run(
            ["samtools", "view", "-bS", "-o", str(output_path), str(sam_path)],
            check=True,
            capture_output=True,
        )
        subprocess.run(
            ["samtools", "index", str(output_path)],
            check=True,
            capture_output=True,
        )
        sam_path.unlink()
        return output_path
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None


@pytest.fixture(scope="session")
def bam_small(benchmark_data_dir: Path, rng: np.random.Generator) -> Path | None:
    """Generate small BAM (10K reads)."""
    return generate_synthetic_bam(
        benchmark_data_dir / "small.bam",
        n_reads=10_000,
        rng=rng,
    )


@pytest.fixture(scope="session")
def bam_medium(benchmark_data_dir: Path, rng: np.random.Generator) -> Path | None:
    """Generate medium BAM (100K reads)."""
    return generate_synthetic_bam(
        benchmark_data_dir / "medium.bam",
        n_reads=100_000,
        rng=rng,
    )


# ============================================================================
# Memory profiling utilities
# ============================================================================


class MemoryBenchmark:
    """Context manager for memory benchmarking."""

    def __init__(self):
        self.peak_memory: float = 0.0
        self.baseline_memory: float = 0.0

    def __enter__(self):
        gc.collect()
        if HAS_MEMORY_PROFILER:
            self.baseline_memory = memory_usage(-1, interval=0.1, timeout=1)[0]
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        gc.collect()
        return False

    def measure(self, func, *args, **kwargs) -> tuple[Any, float]:
        """Execute function and measure peak memory usage."""
        if not HAS_MEMORY_PROFILER:
            result = func(*args, **kwargs)
            return result, 0.0

        mem_usage, result = memory_usage(
            (func, args, kwargs),
            interval=0.1,
            timeout=None,
            retval=True,
            max_usage=True,
        )
        self.peak_memory = mem_usage - self.baseline_memory
        return result, self.peak_memory


@pytest.fixture
def memory_benchmark() -> Generator[MemoryBenchmark, None, None]:
    """Fixture providing memory benchmarking capability."""
    benchmark = MemoryBenchmark()
    with benchmark:
        yield benchmark


# ============================================================================
# Benchmark comparison utilities
# ============================================================================


@pytest.fixture(scope="session")
def benchmark_results_dir() -> Path:
    """Directory for storing benchmark results."""
    results_dir = Path(".benchmarks")
    results_dir.mkdir(exist_ok=True)
    return results_dir


def skip_if_no_tool(tool_name: str):
    """Skip test if external tool is not available."""
    import shutil

    if shutil.which(tool_name) is None:
        pytest.skip(f"{tool_name} not available")


@pytest.fixture
def skip_without_bcftools():
    """Skip if bcftools not available."""
    skip_if_no_tool("bcftools")


@pytest.fixture
def skip_without_samtools():
    """Skip if samtools not available."""
    skip_if_no_tool("samtools")


@pytest.fixture
def skip_without_plink2():
    """Skip if plink2 not available."""
    skip_if_no_tool("plink2")


# ============================================================================
# Parametrized scale fixtures
# ============================================================================


@pytest.fixture(params=["tiny", "small", "medium"])
def vcf_scale(
    request,
    vcf_tiny: Path,
    vcf_small: Path,
    vcf_medium: Path,
) -> tuple[str, Path]:
    """Parametrized fixture for multiple VCF scales."""
    scale_map = {
        "tiny": vcf_tiny,
        "small": vcf_small,
        "medium": vcf_medium,
    }
    return (request.param, scale_map[request.param])
