# WASP2 Multi-Format Variant Support: Design Document

## Executive Summary

This document outlines the design for integrating PLINK2 (PGEN/PVAR/PSAM) format support into WASP2, alongside existing VCF support. The design follows software engineering best practices using the **Strategy + Factory + Registry** pattern to enable extensible, maintainable, and testable multi-format support.

---

## 1. Current State Analysis

### 1.1 Existing VCF Handling in WASP2-exp

| Module | File | VCF Handling | Issues |
|--------|------|--------------|--------|
| mapping | `intersect_variant_data.py` | `vcf_to_bed()` via bcftools subprocess | Duplicated in counting module |
| mapping | `make_remap_reads.py` | Uses BED output from above | Tightly coupled to VCF |
| counting | `filter_variant_data.py` | `vcf_to_bed()` (duplicate) | Code duplication |

### 1.2 Key Problems with Current Architecture

1. **Code Duplication**: `vcf_to_bed()` exists in both mapping and counting modules
2. **Format Lock-in**: Direct bcftools subprocess calls hardcode VCF format
3. **No Abstraction Layer**: Business logic mixed with file format handling
4. **Subprocess Dependency**: Relies on external bcftools binary
5. **No Format Auto-detection**: User must know and specify format

### 1.3 Existing PLINK2 Implementation (WASP2-improved-new)

The `WASP2-improved-new` repo has substantial PLINK2 support:

| File | Status | Quality |
|------|--------|---------|
| `pgen_utils.py` | Complete | Good - handles VCF→PGEN conversion, normalization |
| `pgen_genotype_reader.py` | Complete | Good - reads genotypes via pgenlib |
| `variant_reader.py` | Complete | Good - ABC pattern already implemented |

**What's Good:**
- Abstract `VariantReader` base class
- `VcfVariantReader` and `PgenVariantReader` implementations
- `open_variant_reader()` factory function
- Chunked reading for memory efficiency

**What Needs Improvement:**
- No registry pattern (can't easily add new formats)
- Missing `to_bed()` method for bedtools compatibility
- Not integrated with WASP2-exp's `WaspDataFiles`
- Lacks heterozygous site filtering at the source level

---

## 2. Proposed Architecture

### 2.1 Design Pattern: Strategy + Factory + Registry

```
┌─────────────────────────────────────────────────────────────────────┐
│                        User / CLI Layer                              │
│     wasp2 mapping --variants data.pgen --bam reads.bam              │
└─────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────┐
│                    VariantSourceFactory                              │
│  ┌─────────────────────────────────────────────────────────────┐   │
│  │  Registry: {'.vcf': VCFSource, '.pgen': PGENSource, ...}    │   │
│  └─────────────────────────────────────────────────────────────┘   │
│  • Auto-detect format from extension/magic bytes                    │
│  • Return appropriate VariantSource implementation                  │
│  • @register decorator for extensibility                            │
└─────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────┐
│                  VariantSource (Abstract Base Class)                 │
│  ═══════════════════════════════════════════════════════════════   │
│  Properties:                                                         │
│  • samples: List[str]                                               │
│  • variant_count: int                                               │
│  • sample_count: int                                                │
│                                                                      │
│  Abstract Methods:                                                   │
│  • iter_variants(samples?) -> Iterator[Variant]                     │
│  • get_het_sites(sample) -> Iterator[Variant]                       │
│  • get_genotype(sample, chrom, pos) -> Genotype                     │
│  • query_region(chrom, start, end) -> Iterator[Variant]             │
│  • to_bed(output, samples?, het_only?) -> Path                      │
│                                                                      │
│  Concrete Methods:                                                   │
│  • get_sample_idx(sample_id) -> int                                 │
│  • validate() -> bool                                               │
└─────────────────────────────────────────────────────────────────────┘
            │                       │                       │
            ▼                       ▼                       ▼
┌───────────────────┐   ┌───────────────────┐   ┌───────────────────┐
│    VCFSource      │   │   PGENSource      │   │  Future Formats   │
│   ─────────────   │   │   ────────────    │   │   ─────────────   │
│ • pysam/cyvcf2    │   │ • pgenlib         │   │ • BCF             │
│ • bcftools query  │   │ • Direct binary   │   │ • BGEN            │
│ • Indexed access  │   │ • Chunked read    │   │ • Zarr            │
└───────────────────┘   └───────────────────┘   └───────────────────┘
```

### 2.2 Core Data Structures

```python
from dataclasses import dataclass
from typing import Optional, Tuple
from enum import Enum

class Genotype(Enum):
    """Standardized genotype representation."""
    HOM_REF = 0      # 0/0
    HET = 1          # 0/1 or 1/0
    HOM_ALT = 2      # 1/1
    MISSING = -1     # ./.

@dataclass(frozen=True, slots=True)
class Variant:
    """Immutable variant representation."""
    chrom: str
    pos: int           # 1-based position
    ref: str
    alt: str
    id: Optional[str] = None

    @property
    def pos0(self) -> int:
        """0-based position for BED format."""
        return self.pos - 1

    def to_bed_line(self) -> str:
        """Convert to BED format line."""
        return f"{self.chrom}\t{self.pos0}\t{self.pos}\t{self.ref}\t{self.alt}"

@dataclass
class VariantGenotype:
    """Variant with genotype information."""
    variant: Variant
    genotype: Genotype
    allele1: Optional[str] = None  # For phased data
    allele2: Optional[str] = None

    @property
    def is_het(self) -> bool:
        return self.genotype == Genotype.HET
```

### 2.3 Abstract Base Class

```python
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Iterator, List, Optional, Dict, Any

class VariantSource(ABC):
    """
    Abstract interface for variant data sources.

    Implementations handle format-specific reading while exposing
    a unified API for WASP2's mapping and counting modules.
    """

    # Class-level registry for format handlers
    _registry: Dict[str, type] = {}

    @classmethod
    def register(cls, *extensions: str):
        """Decorator to register format handlers."""
        def decorator(subclass):
            for ext in extensions:
                cls._registry[ext.lower().lstrip('.')] = subclass
            return subclass
        return decorator

    @classmethod
    def open(cls, path: Path, **kwargs) -> 'VariantSource':
        """Factory method with auto-detection."""
        path = Path(path)
        ext = cls._detect_format(path)
        if ext not in cls._registry:
            raise ValueError(f"Unsupported format: {ext}. "
                           f"Supported: {list(cls._registry.keys())}")
        return cls._registry[ext](path, **kwargs)

    @classmethod
    def _detect_format(cls, path: Path) -> str:
        """Detect format from extension, handling compression."""
        suffixes = path.suffixes
        if suffixes[-1] in ('.gz', '.bgz', '.zst'):
            return suffixes[-2].lstrip('.') if len(suffixes) > 1 else ''
        return suffixes[-1].lstrip('.') if suffixes else ''

    # ─────────────────────────────────────────────────────────────
    # Abstract Properties
    # ─────────────────────────────────────────────────────────────

    @property
    @abstractmethod
    def samples(self) -> List[str]:
        """List of sample IDs in the file."""
        ...

    @property
    @abstractmethod
    def variant_count(self) -> int:
        """Total number of variants."""
        ...

    @property
    @abstractmethod
    def sample_count(self) -> int:
        """Total number of samples."""
        ...

    # ─────────────────────────────────────────────────────────────
    # Abstract Methods - Must be implemented by subclasses
    # ─────────────────────────────────────────────────────────────

    @abstractmethod
    def iter_variants(self,
                      samples: Optional[List[str]] = None,
                      het_only: bool = False) -> Iterator[VariantGenotype]:
        """
        Iterate over variants, optionally filtered by sample/het status.

        Args:
            samples: Sample IDs to include (None = all)
            het_only: If True, only yield heterozygous sites

        Yields:
            VariantGenotype objects
        """
        ...

    @abstractmethod
    def get_genotype(self, sample: str, chrom: str, pos: int) -> Genotype:
        """Get genotype for a specific sample at a position."""
        ...

    @abstractmethod
    def query_region(self,
                     chrom: str,
                     start: int,
                     end: int,
                     samples: Optional[List[str]] = None) -> Iterator[VariantGenotype]:
        """Query variants in a genomic region (1-based, inclusive)."""
        ...

    @abstractmethod
    def to_bed(self,
               output: Path,
               samples: Optional[List[str]] = None,
               het_only: bool = True,
               include_genotypes: bool = True) -> Path:
        """
        Export variants to BED format for bedtools intersection.

        This is the key method for WASP2 integration - it replaces
        the current vcf_to_bed() subprocess calls.

        Args:
            output: Output BED file path
            samples: Samples to include
            het_only: Only include heterozygous sites
            include_genotypes: Include genotype columns

        Returns:
            Path to output BED file
        """
        ...

    # ─────────────────────────────────────────────────────────────
    # Concrete Methods - Shared implementation
    # ─────────────────────────────────────────────────────────────

    def get_sample_idx(self, sample_id: str) -> int:
        """Get 0-based index for a sample ID."""
        try:
            return self.samples.index(sample_id)
        except ValueError:
            raise ValueError(f"Sample '{sample_id}' not found. "
                           f"Available: {self.samples[:5]}...")

    def validate(self) -> bool:
        """Validate the variant source is readable."""
        try:
            _ = self.variant_count
            _ = self.sample_count
            return True
        except Exception:
            return False

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def close(self):
        """Clean up resources. Override in subclasses if needed."""
        pass
```

### 2.4 VCF Implementation

```python
@VariantSource.register('vcf', 'vcf.gz', 'bcf')
class VCFSource(VariantSource):
    """VCF/BCF variant source using pysam."""

    def __init__(self, path: Path, **kwargs):
        import pysam
        self.path = Path(path)
        self._vcf = pysam.VariantFile(str(self.path))
        self._samples = list(self._vcf.header.samples)
        self._variant_count = None  # Lazy computation

    @property
    def samples(self) -> List[str]:
        return self._samples

    @property
    def variant_count(self) -> int:
        if self._variant_count is None:
            # Use tabix index if available
            if self.path.suffix == '.gz':
                try:
                    import subprocess
                    result = subprocess.run(
                        ['bcftools', 'index', '--nrecords', str(self.path)],
                        capture_output=True, text=True
                    )
                    self._variant_count = int(result.stdout.strip())
                except:
                    self._variant_count = sum(1 for _ in self._vcf)
                    self._vcf.reset()
            else:
                self._variant_count = sum(1 for _ in self._vcf)
                self._vcf.reset()
        return self._variant_count

    @property
    def sample_count(self) -> int:
        return len(self._samples)

    def iter_variants(self, samples=None, het_only=False):
        self._vcf.reset()
        sample_indices = None
        if samples:
            sample_indices = [self.get_sample_idx(s) for s in samples]

        for record in self._vcf:
            variant = Variant(
                chrom=record.contig,
                pos=record.pos,
                ref=record.ref,
                alt=record.alts[0] if record.alts else '.',
                id=record.id
            )

            # Get genotypes for requested samples
            for idx, sample in enumerate(samples or self._samples):
                gt = record.samples[sample].get('GT', (None, None))
                genotype = self._parse_gt(gt)

                if het_only and genotype != Genotype.HET:
                    continue

                alleles = self._get_alleles(record, gt)
                yield VariantGenotype(
                    variant=variant,
                    genotype=genotype,
                    allele1=alleles[0],
                    allele2=alleles[1]
                )

    def to_bed(self, output, samples=None, het_only=True, include_genotypes=True):
        """Export to BED using bcftools for efficiency."""
        import subprocess

        # Build bcftools pipeline
        view_cmd = ['bcftools', 'view', str(self.path),
                    '-m2', '-M2', '-v', 'snps', '-Ou']

        if samples:
            view_cmd.extend(['-s', ','.join(samples)])
            if het_only and len(samples) == 1:
                # Filter het genotypes
                view_proc = subprocess.run(view_cmd, capture_output=True)
                het_cmd = ['bcftools', 'view', '--genotype', 'het', '-Ou']
                view_proc = subprocess.run(het_cmd, input=view_proc.stdout,
                                          capture_output=True)
                view_output = view_proc.stdout
            else:
                view_proc = subprocess.run(view_cmd, capture_output=True)
                view_output = view_proc.stdout
        else:
            view_cmd.append('--drop-genotypes')
            view_proc = subprocess.run(view_cmd, capture_output=True)
            view_output = view_proc.stdout

        # Query to BED format
        fmt = '%CHROM\t%POS0\t%END\t%REF\t%ALT'
        if include_genotypes and samples:
            fmt += r'[\t%TGT]'
        fmt += '\n'

        query_cmd = ['bcftools', 'query', '-f', fmt, '-o', str(output)]
        subprocess.run(query_cmd, input=view_output, check=True)

        return Path(output)

    def _parse_gt(self, gt) -> Genotype:
        if None in gt:
            return Genotype.MISSING
        if sum(gt) == 0:
            return Genotype.HOM_REF
        if all(a == gt[0] for a in gt):
            return Genotype.HOM_ALT
        return Genotype.HET

    def close(self):
        if self._vcf:
            self._vcf.close()
```

### 2.5 PGEN Implementation

```python
@VariantSource.register('pgen')
class PGENSource(VariantSource):
    """PLINK2 PGEN variant source using pgenlib."""

    def __init__(self, path: Path, **kwargs):
        import pgenlib
        import pandas as pd

        self.path = Path(path)
        self.pvar_path = self.path.with_suffix('.pvar')
        self.psam_path = self.path.with_suffix('.psam')

        # Validate files exist
        for p in [self.path, self.pvar_path, self.psam_path]:
            if not p.exists():
                raise FileNotFoundError(f"Required file not found: {p}")

        # Read sample info
        self._psam_df = self._read_psam()
        self._samples = self._psam_df['IID'].tolist()

        # Read variant info
        self._pvar_df = self._read_pvar()

        # Initialize pgenlib reader with multiallelic support
        allele_counts = self._pvar_df['ALT'].str.count(',') + 2
        self._allele_idx_offsets = np.zeros(len(self._pvar_df) + 1, dtype=np.uintp)
        self._allele_idx_offsets[1:] = np.cumsum(allele_counts)

        self._reader = pgenlib.PgenReader(
            bytes(str(self.path), 'utf-8'),
            allele_idx_offsets=self._allele_idx_offsets
        )

    @property
    def samples(self) -> List[str]:
        return self._samples

    @property
    def variant_count(self) -> int:
        return self._reader.get_variant_ct()

    @property
    def sample_count(self) -> int:
        return self._reader.get_raw_sample_ct()

    def iter_variants(self, samples=None, het_only=False):
        sample_indices = None
        if samples:
            sample_indices = np.array([self.get_sample_idx(s) for s in samples],
                                     dtype=np.uint32)
            self._reader.change_sample_subset(sample_indices)

        genotype_buf = np.empty(2, dtype=np.int32)

        for var_idx in range(self.variant_count):
            row = self._pvar_df.iloc[var_idx]
            variant = Variant(
                chrom=str(row['CHROM']),
                pos=int(row['POS']),
                ref=row['REF'],
                alt=row['ALT'].split(',')[0],  # First alt for biallelic
                id=row.get('ID', '.')
            )

            # Read genotype
            self._reader.read_alleles(var_idx, genotype_buf)
            genotype = self._parse_alleles(genotype_buf)

            if het_only and genotype != Genotype.HET:
                continue

            yield VariantGenotype(
                variant=variant,
                genotype=genotype,
                allele1=self._allele_to_base(genotype_buf[0], variant),
                allele2=self._allele_to_base(genotype_buf[1], variant)
            )

    def to_bed(self, output, samples=None, het_only=True, include_genotypes=True):
        """Export to BED format directly from PGEN."""
        with open(output, 'w') as f:
            for vg in self.iter_variants(samples=samples, het_only=het_only):
                line = vg.variant.to_bed_line()
                if include_genotypes:
                    line += f"\t{vg.allele1}|{vg.allele2}"
                f.write(line + '\n')
        return Path(output)

    def _read_psam(self) -> pd.DataFrame:
        """Read PSAM file with standard column detection."""
        df = pd.read_csv(self.psam_path, sep='\t', dtype=str)
        df.columns = [c.lstrip('#') for c in df.columns]
        return df

    def _read_pvar(self) -> pd.DataFrame:
        """Read PVAR file skipping header comments."""
        return pd.read_csv(self.pvar_path, sep='\t', comment='#',
                          names=['CHROM', 'POS', 'ID', 'REF', 'ALT'],
                          dtype={'CHROM': str, 'POS': int, 'ID': str,
                                'REF': str, 'ALT': str})

    def _parse_alleles(self, buf) -> Genotype:
        if buf[0] < 0 or buf[1] < 0:
            return Genotype.MISSING
        if buf[0] == 0 and buf[1] == 0:
            return Genotype.HOM_REF
        if buf[0] == buf[1]:
            return Genotype.HOM_ALT
        return Genotype.HET

    def _allele_to_base(self, allele_idx: int, variant: Variant) -> str:
        if allele_idx < 0:
            return '.'
        if allele_idx == 0:
            return variant.ref
        alts = variant.alt.split(',')
        return alts[allele_idx - 1] if allele_idx <= len(alts) else '.'

    def close(self):
        if self._reader:
            self._reader.close()
```

---

## 3. Integration Plan

### 3.1 File Structure

```
src/
├── wasp2/
│   ├── __init__.py
│   ├── io/                          # NEW: I/O abstraction layer
│   │   ├── __init__.py
│   │   ├── variant_source.py        # ABC and factory
│   │   ├── vcf_source.py            # VCF implementation
│   │   ├── pgen_source.py           # PGEN implementation
│   │   └── formats/                 # Future formats
│   │       └── __init__.py
│   ├── mapping/
│   │   ├── intersect_variant_data.py  # UPDATED: Use VariantSource
│   │   ├── make_remap_reads.py
│   │   └── ...
│   └── counting/
│       ├── filter_variant_data.py     # UPDATED: Use VariantSource
│       └── ...
```

### 3.2 Migration Steps

| Phase | Task | Changes |
|-------|------|---------|
| 1 | Create `io/` module | New files, no breaking changes |
| 2 | Implement `VCFSource` | Port existing bcftools logic |
| 3 | Implement `PGENSource` | Port from WASP2-improved-new |
| 4 | Update `intersect_variant_data.py` | Replace `vcf_to_bed()` with `source.to_bed()` |
| 5 | Update `filter_variant_data.py` | Remove duplicate `vcf_to_bed()` |
| 6 | Update CLI | Add `--variant-format` auto-detection |
| 7 | Add tests | Unit + integration tests |

### 3.3 Backward Compatibility

```python
# Old code (still works):
from mapping.intersect_variant_data import vcf_to_bed
vcf_to_bed(vcf_file, out_bed, samples)

# New code:
from wasp2.io import VariantSource
with VariantSource.open(variant_file) as source:
    source.to_bed(out_bed, samples=samples, het_only=True)

# The old vcf_to_bed becomes a thin wrapper:
def vcf_to_bed(vcf_file, out_bed, samples=None):
    """Deprecated: Use VariantSource.to_bed() instead."""
    warnings.warn("vcf_to_bed is deprecated, use VariantSource", DeprecationWarning)
    with VariantSource.open(vcf_file) as source:
        return source.to_bed(out_bed, samples=samples, het_only=True)
```

---

## 4. Benchmarking Plan

### 4.1 Metrics to Measure

| Metric | Description | Tool |
|--------|-------------|------|
| **Wall time** | End-to-end execution time | `time` / `timeit` |
| **Peak memory** | Maximum RSS during execution | `/usr/bin/time -v` / `memory_profiler` |
| **I/O throughput** | Variants processed per second | Custom logging |
| **CPU utilization** | User vs system time | `time` |

### 4.2 Test Datasets

| Dataset | Size | Variants | Samples | Source |
|---------|------|----------|---------|--------|
| Small | ~10MB | 100K | 1 | Synthetic |
| Medium | ~500MB | 5M | 10 | 1000 Genomes subset |
| Large | ~5GB | 50M | 100 | iPSCORE subset |
| WGS | ~50GB | 500M | 1 | Full WGS sample |

### 4.3 Benchmark Scenarios

```python
# benchmark_config.py
BENCHMARKS = {
    "vcf_to_bed_single_sample": {
        "description": "Export het sites for single sample to BED",
        "formats": ["vcf", "vcf.gz", "pgen"],
        "samples": [1],
        "het_only": True,
    },
    "vcf_to_bed_multi_sample": {
        "description": "Export het sites for multiple samples",
        "formats": ["vcf", "vcf.gz", "pgen"],
        "samples": [1, 10, 100],
        "het_only": True,
    },
    "full_pipeline_mapping": {
        "description": "Complete WASP mapping pipeline",
        "formats": ["vcf.gz", "pgen"],
        "samples": [1],
        "include": ["vcf_to_bed", "intersect", "remap"],
    },
    "genotype_lookup": {
        "description": "Random access genotype queries",
        "formats": ["vcf.gz", "pgen"],
        "queries": [100, 1000, 10000],
    },
}
```

### 4.4 Benchmark Script Structure

```python
# benchmarks/run_benchmarks.py
import time
import tracemalloc
from pathlib import Path
from dataclasses import dataclass
from typing import List, Dict, Any
import json

@dataclass
class BenchmarkResult:
    name: str
    format: str
    dataset: str
    wall_time_sec: float
    peak_memory_mb: float
    variants_processed: int
    throughput_variants_per_sec: float

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

class VariantSourceBenchmark:
    """Benchmark suite for VariantSource implementations."""

    def __init__(self, output_dir: Path):
        self.output_dir = Path(output_dir)
        self.results: List[BenchmarkResult] = []

    def benchmark_to_bed(self,
                         source_path: Path,
                         samples: List[str],
                         het_only: bool = True,
                         n_runs: int = 3) -> BenchmarkResult:
        """Benchmark the to_bed() operation."""
        from wasp2.io import VariantSource

        times = []
        memories = []

        for _ in range(n_runs):
            tracemalloc.start()
            start = time.perf_counter()

            with VariantSource.open(source_path) as source:
                out_bed = self.output_dir / "bench_output.bed"
                source.to_bed(out_bed, samples=samples, het_only=het_only)
                variant_count = source.variant_count

            elapsed = time.perf_counter() - start
            current, peak = tracemalloc.get_traced_memory()
            tracemalloc.stop()

            times.append(elapsed)
            memories.append(peak / 1024 / 1024)  # MB

        avg_time = sum(times) / len(times)
        avg_memory = sum(memories) / len(memories)

        return BenchmarkResult(
            name="to_bed",
            format=source_path.suffix,
            dataset=source_path.stem,
            wall_time_sec=avg_time,
            peak_memory_mb=avg_memory,
            variants_processed=variant_count,
            throughput_variants_per_sec=variant_count / avg_time
        )

    def run_all(self, datasets: Dict[str, Path]) -> None:
        """Run all benchmarks on all datasets."""
        for name, path in datasets.items():
            # Test different scenarios
            for n_samples in [1, 10]:
                samples = [f"sample_{i}" for i in range(n_samples)]
                result = self.benchmark_to_bed(path, samples)
                self.results.append(result)

        # Save results
        with open(self.output_dir / "benchmark_results.json", "w") as f:
            json.dump([r.to_dict() for r in self.results], f, indent=2)

    def generate_report(self) -> str:
        """Generate markdown benchmark report."""
        # ... generate comparison tables and charts
```

### 4.5 Expected Performance Comparison

| Operation | VCF (bcftools) | VCF (pysam) | PGEN (pgenlib) | Expected Winner |
|-----------|----------------|-------------|----------------|-----------------|
| Load metadata | Fast | Medium | Fast | Tie |
| Single sample het export | Medium | Slow | Fast | PGEN (2-3x) |
| Multi-sample het export | Medium | Slow | Fast | PGEN (5-10x) |
| Random access query | Fast (indexed) | Fast | Fast | Tie |
| Memory (large file) | Low (streaming) | High | Low | VCF/PGEN |
| Full pipeline | Baseline | - | TBD | TBD |

### 4.6 Validation Tests

```python
def validate_output_equivalence(vcf_path: Path, pgen_path: Path, sample: str):
    """Ensure VCF and PGEN produce identical BED output."""
    from wasp2.io import VariantSource

    with VariantSource.open(vcf_path) as vcf_source:
        vcf_source.to_bed(Path("/tmp/vcf.bed"), samples=[sample])

    with VariantSource.open(pgen_path) as pgen_source:
        pgen_source.to_bed(Path("/tmp/pgen.bed"), samples=[sample])

    # Compare outputs
    import filecmp
    assert filecmp.cmp("/tmp/vcf.bed", "/tmp/pgen.bed"), \
        "VCF and PGEN outputs differ!"
```

---

## 5. Testing Strategy

### 5.1 Unit Tests

```python
# tests/test_variant_source.py
import pytest
from wasp2.io import VariantSource, VCFSource, PGENSource

class TestVariantSourceFactory:
    def test_auto_detect_vcf(self, vcf_file):
        source = VariantSource.open(vcf_file)
        assert isinstance(source, VCFSource)

    def test_auto_detect_pgen(self, pgen_file):
        source = VariantSource.open(pgen_file)
        assert isinstance(source, PGENSource)

    def test_unsupported_format(self, tmp_path):
        bad_file = tmp_path / "data.xyz"
        bad_file.touch()
        with pytest.raises(ValueError, match="Unsupported format"):
            VariantSource.open(bad_file)

class TestVCFSource:
    def test_samples(self, vcf_file):
        with VCFSource(vcf_file) as source:
            assert len(source.samples) > 0

    def test_iter_het_only(self, vcf_file):
        with VCFSource(vcf_file) as source:
            het_sites = list(source.iter_variants(het_only=True))
            for site in het_sites:
                assert site.genotype == Genotype.HET

class TestPGENSource:
    def test_samples(self, pgen_file):
        with PGENSource(pgen_file) as source:
            assert len(source.samples) > 0

    def test_to_bed_matches_vcf(self, vcf_file, pgen_file, tmp_path):
        """Ensure PGEN and VCF produce equivalent BED output."""
        # ... comparison test
```

### 5.2 Integration Tests

```python
# tests/test_integration.py
class TestMappingPipeline:
    def test_full_pipeline_vcf(self, vcf_file, bam_file):
        """Test complete mapping pipeline with VCF input."""
        # ... end-to-end test

    def test_full_pipeline_pgen(self, pgen_file, bam_file):
        """Test complete mapping pipeline with PGEN input."""
        # ... end-to-end test

    def test_pipeline_equivalence(self, vcf_file, pgen_file, bam_file):
        """Ensure VCF and PGEN produce identical WASP results."""
        # ... comparison test
```

---

## 6. Timeline and Milestones

| Week | Milestone | Deliverables |
|------|-----------|--------------|
| 1 | Core architecture | `VariantSource` ABC, factory, data classes |
| 2 | VCF implementation | `VCFSource` with full test coverage |
| 3 | PGEN implementation | `PGENSource` ported and tested |
| 4 | Integration | Update mapping/counting modules |
| 5 | Benchmarking | Run benchmarks, generate report |
| 6 | Documentation | Update docs, examples, migration guide |

---

## 7. Risks and Mitigations

| Risk | Impact | Mitigation |
|------|--------|------------|
| pgenlib API changes | High | Pin version, add compatibility layer |
| Performance regression | Medium | Benchmark at each phase |
| bcftools dependency | Low | Keep as fallback option |
| Memory issues with large files | Medium | Ensure streaming/chunked processing |

---

## 8. References

- [Stack Overflow: Design patterns for multiple file formats](https://stackoverflow.com/questions/35139016/which-design-pattern-to-use-to-process-different-files-in-java)
- [Hail Import/Export](https://hail.is/docs/0.2/methods/impex.html)
- [scikit-allel I/O utilities](https://scikit-allel.readthedocs.io/en/stable/io.html)
- [pgenlib Python API](https://github.com/chrchang/plink-ng/tree/master/2.0/Python)
- [PLINK2 file formats](https://www.cog-genomics.org/plink/2.0/formats)
