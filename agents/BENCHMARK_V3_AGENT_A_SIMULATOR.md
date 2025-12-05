# Agent A: ASE Simulator V3 - Core Engine

## Mission
Refactor `simulation/simulate_paired_end_ase.py` to fix critical design flaws and create publication-quality ground truth data.

---

## Repository Context

**GitHub:** https://github.com/Jaureguy760/WASP2-exp.git
**Branch:** `sim/benchmark-v3`
**Working Directory:** `/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp`
**Conda Environment:** `WASP2_dev2`

---

## Problem Statement

### Current Flaws (MUST FIX)

1. **Shared positions**: All replicates/coverages share the same genomic positions
   - Result: Counts are cumulative, not independent
   - Fix: Each test case at UNIQUE position

2. **No allele tracking**: Can't determine which reads came from which haplotype
   - Result: No perfect ground truth validation
   - Fix: Encode haplotype in read names

3. **Replicate design**: Multiple replicates at same position are meaningless
   - Result: Validation is broken
   - Fix: Each test is one variant at one unique position

---

## Design Specification

### Read Name Format (CRITICAL)

```
@{chrom}_{pos}_{read_idx}_{HAP1|HAP2}/1
@{chrom}_{pos}_{read_idx}_{HAP1|HAP2}/2

Examples:
@chr1_50000_001_HAP1/1    # R1 from REF haplotype at pos 50000
@chr1_50000_001_HAP1/2    # R2 (mate) from REF haplotype
@chr1_50000_002_HAP2/1    # R1 from ALT haplotype
@chr1_50000_002_HAP2/2    # R2 (mate) from ALT haplotype
```

This enables PERFECT ground truth validation by parsing read names.

### Unique Positions

```python
# BAD (current):
for replicate in range(10):
    for coverage in [20, 50, 100]:
        generate_reads(pos=50000, ...)  # SAME POSITION!

# GOOD (new):
position_counter = 50000
for vtype in ['SNP', 'INS', 'DEL']:
    for size in [1, 3, 5, 10]:
        for coverage in [20, 50, 100]:
            for ratio in [1.0, 2.0, 4.0]:
                generate_reads(pos=position_counter, ...)
                position_counter += 10000  # UNIQUE POSITION!
```

---

## Implementation

### File: `simulation/ase_simulator_v3.py`

```python
#!/usr/bin/env python3
"""
ASE Simulator V3 - Publication-Quality Ground Truth Generation

Key Features:
1. Unique genomic position per test case
2. Allele origin encoded in read names (@pos_idx_HAP1/HAP2)
3. Perfect ground truth tracking
4. Support for SNPs, insertions, deletions of various sizes
"""

import pysam
import random
import numpy as np
import pandas as pd
from pathlib import Path
from dataclasses import dataclass, asdict
from typing import List, Tuple, Dict, Optional
import gzip
import subprocess


@dataclass
class Variant:
    """A single variant to test."""
    chrom: str
    pos: int
    ref: str
    alt: str
    variant_type: str  # SNP, INS, DEL
    variant_size: int  # Size of INDEL (1 for SNP)


@dataclass
class TestCase:
    """A single test case with ground truth."""
    variant: Variant
    coverage: int
    true_ratio: float  # REF:ALT ratio (1.0 = balanced)
    expected_ref: int
    expected_alt: int
    seed: int


@dataclass
class SimulatedRead:
    """A simulated read with known origin."""
    name: str  # Format: {chrom}_{pos}_{idx}_{HAP1|HAP2}
    seq_r1: str
    qual_r1: str
    seq_r2: str
    qual_r2: str
    haplotype: str  # HAP1 (REF) or HAP2 (ALT)


class ASESimulatorV3:
    """
    Publication-quality ASE simulation with ground truth tracking.

    Usage:
        simulator = ASESimulatorV3(reference_fasta='ref.fa', seed=42)
        test_cases = simulator.generate_test_suite(output_dir='output/')
    """

    def __init__(
        self,
        reference_fasta: Optional[str] = None,
        genome_size: int = 10_000_000,  # 10Mb synthetic genome
        seed: int = 42
    ):
        """
        Initialize simulator.

        Args:
            reference_fasta: Path to reference FASTA (or None for synthetic)
            genome_size: Size of synthetic genome if no reference provided
            seed: Random seed for reproducibility
        """
        self.seed = seed
        random.seed(seed)
        np.random.seed(seed)

        if reference_fasta and Path(reference_fasta).exists():
            self.reference = pysam.FastaFile(reference_fasta)
            self.use_synthetic = False
        else:
            self.synthetic_genome = self._generate_synthetic_genome(genome_size)
            self.use_synthetic = True

    def _generate_synthetic_genome(self, size: int) -> str:
        """Generate random genome sequence."""
        return ''.join(random.choices('ACGT', k=size))

    def _get_sequence(self, chrom: str, start: int, end: int) -> str:
        """Get sequence from reference or synthetic genome."""
        if self.use_synthetic:
            return self.synthetic_genome[start:end]
        else:
            return self.reference.fetch(chrom, start, end)

    def generate_variant_set(
        self,
        n_snps: int = 100,
        n_insertions: int = 100,
        n_deletions: int = 100,
        min_spacing: int = 10000,
        chrom: str = 'chr1'
    ) -> List[Variant]:
        """
        Generate non-overlapping variants at unique positions.

        Each variant is at a unique position with minimum spacing.
        """
        variants = []
        current_pos = 50000  # Start position

        # SNPs
        for i in range(n_snps):
            ref_base = self._get_sequence(chrom, current_pos, current_pos + 1)
            alt_base = random.choice([b for b in 'ACGT' if b != ref_base])
            variants.append(Variant(
                chrom=chrom,
                pos=current_pos,
                ref=ref_base,
                alt=alt_base,
                variant_type='SNP',
                variant_size=1
            ))
            current_pos += min_spacing

        # Insertions (various sizes)
        insertion_sizes = [1, 2, 3, 4, 5, 6, 8, 10, 15, 20]
        for i in range(n_insertions):
            size = insertion_sizes[i % len(insertion_sizes)]
            ref_base = self._get_sequence(chrom, current_pos, current_pos + 1)
            insert_seq = ''.join(random.choices('ACGT', k=size))
            variants.append(Variant(
                chrom=chrom,
                pos=current_pos,
                ref=ref_base,
                alt=ref_base + insert_seq,
                variant_type='INS',
                variant_size=size
            ))
            current_pos += min_spacing

        # Deletions (various sizes)
        deletion_sizes = [1, 2, 3, 4, 5, 6, 8, 10, 15, 20]
        for i in range(n_deletions):
            size = deletion_sizes[i % len(deletion_sizes)]
            ref_seq = self._get_sequence(chrom, current_pos, current_pos + size + 1)
            variants.append(Variant(
                chrom=chrom,
                pos=current_pos,
                ref=ref_seq,
                alt=ref_seq[0],  # Keep first base, delete rest
                variant_type='DEL',
                variant_size=size
            ))
            current_pos += min_spacing

        return variants

    def generate_reads_for_test_case(
        self,
        test_case: TestCase,
        read_length: int = 150,
        insert_mean: int = 300,
        insert_std: int = 50,
        error_rate: float = 0.001
    ) -> List[SimulatedRead]:
        """
        Generate paired-end reads for a single test case.

        Read names encode the haplotype origin for perfect ground truth.
        """
        var = test_case.variant
        reads = []

        # Get reference context (1kb around variant)
        context_start = max(0, var.pos - 500)
        context_end = var.pos + 500
        ref_context = self._get_sequence(var.chrom, context_start, context_end)

        # Generate REF haplotype reads
        for i in range(test_case.expected_ref):
            read_name = f"{var.chrom}_{var.pos}_{i:04d}_HAP1"
            r1_seq, r1_qual, r2_seq, r2_qual = self._make_paired_read(
                ref_context,
                var.pos - context_start,  # Relative position
                var.ref,  # Use REF allele
                var.ref,
                read_length,
                insert_mean,
                insert_std,
                error_rate
            )
            reads.append(SimulatedRead(
                name=read_name,
                seq_r1=r1_seq,
                qual_r1=r1_qual,
                seq_r2=r2_seq,
                qual_r2=r2_qual,
                haplotype='HAP1'
            ))

        # Generate ALT haplotype reads
        for i in range(test_case.expected_alt):
            read_name = f"{var.chrom}_{var.pos}_{i:04d}_HAP2"
            r1_seq, r1_qual, r2_seq, r2_qual = self._make_paired_read(
                ref_context,
                var.pos - context_start,
                var.ref,
                var.alt,  # Use ALT allele
                read_length,
                insert_mean,
                insert_std,
                error_rate
            )
            reads.append(SimulatedRead(
                name=read_name,
                seq_r1=r1_seq,
                qual_r1=r1_qual,
                seq_r2=r2_seq,
                qual_r2=r2_qual,
                haplotype='HAP2'
            ))

        return reads

    def _make_paired_read(
        self,
        context_seq: str,
        var_pos_rel: int,
        ref_allele: str,
        use_allele: str,
        read_length: int,
        insert_mean: int,
        insert_std: int,
        error_rate: float
    ) -> Tuple[str, str, str, str]:
        """Generate a single paired-end read."""
        # Generate insert size
        insert_size = max(read_length + 50, int(np.random.normal(insert_mean, insert_std)))

        # Position fragment to cover variant
        offset = random.randint(-50, 50)
        frag_start = max(0, var_pos_rel - insert_size // 2 + offset)

        # Build fragment with allele
        left_seq = context_seq[frag_start:var_pos_rel]
        right_seq = context_seq[var_pos_rel + len(ref_allele):]
        fragment = left_seq + use_allele + right_seq

        # Ensure fragment is long enough
        if len(fragment) < insert_size:
            fragment = fragment + context_seq[len(fragment):insert_size]
        fragment = fragment[:insert_size]

        # R1: forward strand
        r1_seq = fragment[:read_length]
        r1_seq = self._add_errors(r1_seq, error_rate)
        r1_qual = self._generate_quality(read_length)

        # R2: reverse complement of end
        r2_seq = self._reverse_complement(fragment[-read_length:])
        r2_seq = self._add_errors(r2_seq, error_rate)
        r2_qual = self._generate_quality(read_length)

        return r1_seq, r1_qual, r2_seq, r2_qual

    def _reverse_complement(self, seq: str) -> str:
        """Return reverse complement."""
        comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        return ''.join(comp.get(b, 'N') for b in reversed(seq))

    def _add_errors(self, seq: str, rate: float) -> str:
        """Add sequencing errors."""
        seq_list = list(seq)
        for i in range(len(seq_list)):
            if random.random() < rate:
                seq_list[i] = random.choice([b for b in 'ACGT' if b != seq_list[i]])
        return ''.join(seq_list)

    def _generate_quality(self, length: int) -> str:
        """Generate quality string."""
        quals = np.random.normal(35, 5, length)
        quals = np.clip(quals, 10, 40).astype(int)
        return ''.join(chr(q + 33) for q in quals)

    def generate_test_suite(
        self,
        output_dir: str,
        n_snps: int = 100,
        n_insertions: int = 100,
        n_deletions: int = 100,
        coverages: List[int] = [20, 50, 100],
        ratios: List[float] = [1.0, 2.0, 4.0]
    ) -> pd.DataFrame:
        """
        Generate complete test suite with ground truth.

        Returns DataFrame with all test cases and expected values.
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Generate variants
        variants = self.generate_variant_set(n_snps, n_insertions, n_deletions)

        # Generate test cases
        test_cases = []
        all_reads = []

        var_idx = 0
        for coverage in coverages:
            for ratio in ratios:
                for var in variants:
                    # Calculate expected counts
                    alt_fraction = 1.0 / (1.0 + ratio)
                    n_alt = int(coverage * alt_fraction)
                    n_ref = coverage - n_alt

                    test = TestCase(
                        variant=var,
                        coverage=coverage,
                        true_ratio=ratio,
                        expected_ref=n_ref,
                        expected_alt=n_alt,
                        seed=self.seed + var_idx
                    )
                    test_cases.append(test)

                    # Generate reads
                    reads = self.generate_reads_for_test_case(test)
                    all_reads.extend(reads)

                    var_idx += 1

        # Write FASTQ files
        self._write_fastq(all_reads, output_dir / 'R1.fq.gz', output_dir / 'R2.fq.gz')

        # Write VCF
        self._write_vcf(variants, output_dir / 'variants.vcf.gz')

        # Write ground truth
        ground_truth = self._create_ground_truth_df(test_cases)
        ground_truth.to_csv(output_dir / 'ground_truth.csv', index=False)

        # Write reference if synthetic
        if self.use_synthetic:
            self._write_reference(output_dir / 'reference.fa')

        return ground_truth

    def _write_fastq(self, reads: List[SimulatedRead], r1_path: Path, r2_path: Path):
        """Write paired FASTQ files."""
        with gzip.open(r1_path, 'wt') as r1_fh, gzip.open(r2_path, 'wt') as r2_fh:
            for read in reads:
                # R1
                r1_fh.write(f"@{read.name}/1\n")
                r1_fh.write(f"{read.seq_r1}\n")
                r1_fh.write("+\n")
                r1_fh.write(f"{read.qual_r1}\n")

                # R2
                r2_fh.write(f"@{read.name}/2\n")
                r2_fh.write(f"{read.seq_r2}\n")
                r2_fh.write("+\n")
                r2_fh.write(f"{read.qual_r2}\n")

    def _write_vcf(self, variants: List[Variant], vcf_path: Path):
        """Write VCF file."""
        vcf_content = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSIMULATED
"""
        for var in variants:
            var_id = f"{var.variant_type}_{var.pos}_{var.variant_size}bp"
            vcf_content += f"{var.chrom}\t{var.pos}\t{var_id}\t{var.ref}\t{var.alt}\t60\tPASS\t.\tGT\t0|1\n"

        # Write and compress
        uncompressed = vcf_path.with_suffix('')
        with open(uncompressed, 'w') as f:
            f.write(vcf_content)

        subprocess.run(['bgzip', '-f', str(uncompressed)], check=True)
        subprocess.run(['tabix', '-p', 'vcf', str(vcf_path)], check=True)

    def _write_reference(self, ref_path: Path):
        """Write synthetic reference genome."""
        with open(ref_path, 'w') as f:
            f.write(">chr1\n")
            # Write in 80-char lines
            for i in range(0, len(self.synthetic_genome), 80):
                f.write(self.synthetic_genome[i:i+80] + "\n")

        # Index
        subprocess.run(['samtools', 'faidx', str(ref_path)], check=True)
        subprocess.run(['bwa', 'index', str(ref_path)], check=True)

    def _create_ground_truth_df(self, test_cases: List[TestCase]) -> pd.DataFrame:
        """Create ground truth DataFrame."""
        rows = []
        for tc in test_cases:
            rows.append({
                'chrom': tc.variant.chrom,
                'pos': tc.variant.pos,
                'ref': tc.variant.ref,
                'alt': tc.variant.alt,
                'variant_type': tc.variant.variant_type,
                'variant_size': tc.variant.variant_size,
                'coverage': tc.coverage,
                'true_ratio': tc.true_ratio,
                'expected_ref': tc.expected_ref,
                'expected_alt': tc.expected_alt,
                'seed': tc.seed
            })
        return pd.DataFrame(rows)


def main():
    """CLI entry point."""
    import argparse

    parser = argparse.ArgumentParser(description='ASE Simulator V3')
    parser.add_argument('--output', '-o', required=True, help='Output directory')
    parser.add_argument('--reference', '-r', help='Reference FASTA (optional)')
    parser.add_argument('--n-snps', type=int, default=100)
    parser.add_argument('--n-insertions', type=int, default=100)
    parser.add_argument('--n-deletions', type=int, default=100)
    parser.add_argument('--coverages', default='20,50,100')
    parser.add_argument('--ratios', default='1.0,2.0,4.0')
    parser.add_argument('--seed', type=int, default=42)

    args = parser.parse_args()

    coverages = [int(x) for x in args.coverages.split(',')]
    ratios = [float(x) for x in args.ratios.split(',')]

    simulator = ASESimulatorV3(
        reference_fasta=args.reference,
        seed=args.seed
    )

    ground_truth = simulator.generate_test_suite(
        output_dir=args.output,
        n_snps=args.n_snps,
        n_insertions=args.n_insertions,
        n_deletions=args.n_deletions,
        coverages=coverages,
        ratios=ratios
    )

    print(f"Generated {len(ground_truth)} test cases")
    print(f"Output: {args.output}")


if __name__ == '__main__':
    main()
```

---

## Tasks

### Task 1: Create `simulation/ase_simulator_v3.py`
- Implement the full class as specified above
- Test with minimal parameters first

### Task 2: Test the simulator
```bash
source /iblm/netapp/home/jjaureguy/mambaforge/etc/profile.d/conda.sh
conda activate WASP2_dev2
cd /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp

python simulation/ase_simulator_v3.py \
    --output /tmp/sim_v3_test \
    --n-snps 10 \
    --n-insertions 10 \
    --n-deletions 10 \
    --coverages 50 \
    --ratios 1.0
```

### Task 3: Validate read names contain haplotype info
```bash
zcat /tmp/sim_v3_test/R1.fq.gz | head -20
# Should see: @chr1_50000_0001_HAP1/1, @chr1_50000_0002_HAP2/1, etc.
```

### Task 4: Validate unique positions
```bash
zcat /tmp/sim_v3_test/variants.vcf.gz | grep -v "^#" | cut -f2 | sort -u | wc -l
# Should equal total number of variants (30 in test case)
```

---

## Success Criteria

- [ ] Each variant at UNIQUE genomic position
- [ ] Read names encode haplotype: `@{chrom}_{pos}_{idx}_{HAP1|HAP2}/1`
- [ ] Ground truth CSV has expected_ref and expected_alt columns
- [ ] VCF is valid (bgzip + tabix indexed)
- [ ] FASTQ files are paired correctly
- [ ] Reference is indexed (BWA + samtools)
- [ ] Test with 30 variants passes (10 SNP + 10 INS + 10 DEL)

---

## Commit Template

```bash
git add simulation/ase_simulator_v3.py
git commit -m "feat(sim): implement ASE simulator v3 with ground truth tracking

Key improvements over v2:
- Each variant at unique genomic position
- Haplotype origin encoded in read names (@pos_idx_HAP1/HAP2)
- Perfect ground truth for validation
- Support for SNPs and INDELs of various sizes

Test matrix: 300 variants × 3 coverages × 3 ratios = 2,700 tests
"
```
