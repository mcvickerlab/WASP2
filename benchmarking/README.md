# WASP2 Benchmark Framework

Reproducible benchmarking framework for validating WASP2's performance claims.

## Performance Claims

WASP2 makes the following performance claims that this framework validates:

| Claim | Description | Benchmark |
|-------|-------------|-----------|
| **61x faster** | WASP filtering vs WASP v1 | `benchmark_mapping.py` |
| **6.4x faster** | Counting vs phASER | `benchmark_counting.py` |
| **r² > 0.99** | Concordance with GATK ASEReadCounter | `benchmark_concordance.py` |

## Quick Start

```bash
# Run all benchmarks with quick settings
python benchmarking/run_benchmarks.py --quick

# Run specific benchmark types
python benchmarking/run_benchmarks.py --counting
python benchmarking/run_benchmarks.py --mapping
python benchmarking/run_benchmarks.py --concordance

# Run comprehensive benchmarks
python benchmarking/run_benchmarks.py --all --n-variants 100000 --iterations 10
```

## Directory Structure

```
benchmarking/
├── __init__.py              # Package initialization
├── utils.py                 # Shared utilities (timing, reporting)
├── run_benchmarks.py        # Main CLI entry point
├── scripts/
│   ├── benchmark_counting.py    # Counting speed benchmarks
│   ├── benchmark_mapping.py     # Mapping filter benchmarks
│   └── benchmark_concordance.py # Accuracy validation
├── data/                    # Generated test data (gitignored)
└── results/                 # Benchmark results (JSON)
```

## Benchmark Descriptions

### 1. Counting Speed (`--counting`)

Compares WASP2's allelic imbalance analysis against phASER and GATK ASEReadCounter.

**What it measures:**
- Time to process allele count data through the analysis pipeline
- Statistical analysis performance (beta-binomial fitting)
- Region-level aggregation speed

**Expected result:** WASP2 should be approximately 6.4x faster than phASER.

```bash
python benchmarking/scripts/benchmark_counting.py --n-variants 10000
```

### 2. Mapping Filter Speed (`--mapping`)

Compares WASP2's Rust-accelerated read filtering against WASP v1's Python implementation.

**What it measures:**
- Time to filter reads that fail to remap to the same location
- BAM I/O performance
- Read position comparison efficiency

**Expected result:** WASP2 should be approximately 61x faster than WASP v1.

```bash
python benchmarking/scripts/benchmark_mapping.py --n-reads 10000
```

### 3. Concordance Validation (`--concordance`)

Validates that WASP2 produces accurate allele counts compared to GATK ASEReadCounter.

**What it measures:**
- Correlation (r²) between WASP2 and GATK allele counts
- Systematic biases in counting
- Edge case handling

**Expected result:** r² > 0.99 correlation.

```bash
python benchmarking/scripts/benchmark_concordance.py --n-variants 1000
```

## Reproducing Benchmarks

### Prerequisites

1. **Install WASP2 with Rust extension:**
   ```bash
   conda activate WASP2
   maturin develop --release -m rust/Cargo.toml
   pip install -e ".[dev]"
   ```

2. **Install benchmark dependencies:**
   ```bash
   pip install pytest-benchmark memory-profiler matplotlib seaborn
   ```

3. **Optional external tools (for real comparisons):**
   ```bash
   # phASER (for counting comparison)
   pip install phaser

   # GATK (for concordance validation)
   # Follow GATK installation instructions

   # WASP v1 (for mapping comparison)
   git clone https://github.com/bmvdgeijn/WASP.git
   ```

### Running Benchmarks

**Quick validation (CI default):**
```bash
python benchmarking/run_benchmarks.py --quick
```

**Full benchmark suite:**
```bash
python benchmarking/run_benchmarks.py --all \
    --n-variants 100000 \
    --n-reads 100000 \
    --iterations 10 \
    --output benchmarking/results/full_benchmark.json
```

**Pytest-benchmark integration:**
```bash
python tests/benchmarks/run_benchmarks.py --groups variant_scaling sample_scaling
```

### Interpreting Results

Results are saved as JSON with the following structure:

```json
{
  "timestamp": "2024-01-15T10:30:00",
  "benchmarks": [
    {
      "name": "wasp2_analysis_10000",
      "tool": "WASP2",
      "mean": 0.123,
      "std": 0.005,
      "min": 0.118,
      "max": 0.131,
      "iterations": 5,
      "parameters": {
        "n_variants": 10000,
        "n_regions": 1000
      }
    }
  ]
}
```

### CI Integration

Benchmarks run automatically via GitHub Actions:

- **Weekly:** Quick benchmarks to track performance regressions
- **Releases:** Full benchmark suite with baseline comparison
- **Manual:** Trigger via workflow_dispatch

See `.github/workflows/benchmarks.yml` for configuration.

## Methodology

### Statistical Rigor

- **Warmup:** 2 iterations discarded before measurement
- **Iterations:** Minimum 5 timed iterations
- **Garbage collection:** Forced between iterations
- **Metrics:** Mean, standard deviation, min, max, median

### Simulated Benchmarks

When external tools (phASER, GATK, WASP v1) are not installed, benchmarks use
simulated workloads based on published performance characteristics. These
simulations:

1. Reproduce the algorithmic complexity of each tool
2. Apply realistic overhead multipliers from published benchmarks
3. Enable CI testing without external dependencies

For definitive validation, install the actual tools.

### Data Generation

Synthetic data is generated with:
- Reproducible random seed (42)
- Realistic allele count distributions (beta-binomial)
- Representative genomic positions
- Configurable scale (100 to 1M+ variants)

## Extending the Framework

### Adding New Benchmarks

1. Create a new script in `benchmarking/scripts/`
2. Use `BenchmarkTimer` from `utils.py` for timing
3. Return `BenchmarkResult` objects
4. Add to `run_benchmarks.py` CLI

Example:
```python
from benchmarking.utils import BenchmarkTimer, BenchmarkResult

def benchmark_new_feature(n_items: int, iterations: int = 5) -> BenchmarkResult:
    timer = BenchmarkTimer("new_feature", iterations=iterations)

    for t in timer:
        with t:
            run_feature(n_items)

    timer.result.tool = "WASP2"
    timer.result.parameters = {"n_items": n_items}
    return timer.result
```

### Custom Comparisons

To add comparison with a new tool:

1. Implement a wrapper class following the pattern in `benchmark_counting.py`
2. Add tool availability check
3. Implement simulated fallback for CI

## Troubleshooting

### "Rust extension not available"

Build the Rust extension:
```bash
maturin develop --release -m rust/Cargo.toml
```

### "samtools not found"

Install bioinformatics tools:
```bash
conda install -c bioconda samtools bcftools
```

### Memory errors with large datasets

Reduce benchmark scale:
```bash
python benchmarking/run_benchmarks.py --n-variants 10000 --n-reads 10000
```

## License

MIT License - see LICENSE file in project root.
