# WASP2 Simulation Benchmark - Sub-Agent Engineering Prompts

**Created:** 2025-12-03
**Parent Branch:** `ropc-indels`
**Goal:** Publication-ready simulation benchmarks for Nature Methods

---

## Parallelization Map

```
                    ┌─────────────────────────────────────────┐
                    │           CAN START NOW                 │
                    └─────────────────────────────────────────┘
                                      │
          ┌───────────────────────────┼───────────────────────────┐
          │                           │                           │
          ▼                           ▼                           ▼
┌─────────────────┐       ┌─────────────────┐       ┌─────────────────┐
│ AGENT 1         │       │ AGENT 2         │       │ AGENT 3         │
│ sim/comprehensive│       │ sim/paired-end  │       │ sim/gatk-compare│
│                 │       │                 │       │                 │
│ RUN EXISTING    │       │ WRITE NEW CODE  │       │ WRITE NEW CODE  │
│ (no code changes)│       │                 │       │                 │
│ ~2 hrs runtime  │       │ ~2-3 hrs dev    │       │ ~2 hrs dev      │
└────────┬────────┘       └────────┬────────┘       └────────┬────────┘
         │                         │                         │
         │                         ▼                         │
         │              ┌─────────────────┐                  │
         │              │ RUN PAIRED-END  │                  │
         │              │ SIMULATION      │                  │
         │              │ ~2 hrs runtime  │                  │
         │              └────────┬────────┘                  │
         │                       │                           │
         └───────────────────────┼───────────────────────────┘
                                 │
                                 ▼
                    ┌─────────────────────────────────────────┐
                    │           AGENT 4: sim/metrics          │
                    │     (AFTER all data is ready)           │
                    │     - Compute publication metrics       │
                    │     - Generate figures                  │
                    │     - Compare WASP2 vs GATK vs truth    │
                    └─────────────────────────────────────────┘
```

---

## AGENT 1: sim/comprehensive (RUN NOW - NO CODE CHANGES)

### Branch
```bash
git checkout sim/comprehensive
```

### Priority
**P0 - CRITICAL** - Can start immediately, provides baseline data

### Task Description
Run the existing simulation framework at comprehensive tier (810 tests) to establish baseline results. No code changes required - just execution and validation.

### Engineering Prompt

```
You are working on the WASP2 bioinformatics project to run comprehensive simulation benchmarks.

## Context
- Branch: sim/comprehensive
- Working directory: /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp
- Conda environment: WASP2_dev2

## Your Task
Run the existing simulation framework at comprehensive tier (810 tests) and validate results.

## Steps

1. **Activate environment and checkout branch:**
   ```bash
   source /iblm/netapp/home/jjaureguy/mambaforge/etc/profile.d/conda.sh
   conda activate WASP2_dev2
   git checkout sim/comprehensive
   ```

2. **Create timestamped output directory:**
   ```bash
   TIMESTAMP=$(date +%Y%m%d_%H%M%S)
   OUTDIR="simulation_results/comprehensive_${TIMESTAMP}"
   mkdir -p ${OUTDIR}
   ```

3. **Run comprehensive simulation via SGE:**
   Create and submit job script:
   ```bash
   cat > run_comprehensive_sim.sh << 'EOF'
   #!/bin/bash
   #$ -N wasp2_sim_comprehensive
   #$ -V
   #$ -pe iblm 8
   #$ -l h_vmem=16G
   #$ -j y
   #$ -o simulation_results/
   #$ -cwd

   source /iblm/netapp/home/jjaureguy/mambaforge/etc/profile.d/conda.sh
   conda activate WASP2_dev2

   python simulate_indel_ase_v2.py \
       --tier comprehensive \
       --workdir ${OUTDIR} \
       --keep
   EOF

   qsub run_comprehensive_sim.sh
   ```

4. **Monitor and validate results:**
   - Expected: 810 tests
   - Expected runtime: ~2 hours
   - Check pass rate (target: >90%)
   - Check error distribution by variant type

5. **Analyze results:**
   ```python
   import pandas as pd
   df = pd.read_csv(f'{OUTDIR}/simulation_results.csv')

   print(f"Total tests: {len(df)}")
   print(f"Pass rate: {(df.status == 'PASS').mean()*100:.1f}%")

   # By variant type
   for vtype in ['SNP', 'INS', 'DEL']:
       sub = df[df.variant_type == vtype]
       print(f"{vtype}: {(sub.status=='PASS').mean()*100:.0f}% pass, mean error {sub.error_pct.mean():.1f}%")
   ```

6. **Commit results summary:**
   ```bash
   git add simulation_results/comprehensive_*/simulation_results.csv
   git commit -m "data: comprehensive simulation results (810 tests)

   Results:
   - Total tests: 810
   - Pass rate: XX%
   - SNP accuracy: XX%
   - INS accuracy: XX%
   - DEL accuracy: XX%
   "
   ```

## Success Criteria
- [ ] 810 tests completed
- [ ] Pass rate > 90%
- [ ] Results CSV saved and committed
- [ ] Summary statistics documented

## Output Files
- `simulation_results/comprehensive_TIMESTAMP/simulation_results.csv`
- `simulation_results/comprehensive_TIMESTAMP/reference.fa`
- `simulation_results/comprehensive_TIMESTAMP/variants.vcf.gz`

## Do NOT
- Modify any simulation code
- Change test parameters
- Skip any tests
```

---

## AGENT 2: sim/paired-end (WRITE CODE - PARALLEL)

### Branch
```bash
git checkout sim/paired-end
```

### Priority
**P0 - CRITICAL** - Core fix for publication (current sim is single-end only)

### Task Description
Create a new paired-end simulation module that generates R1/R2 FASTQ files with proper insert sizes, matching real RNA-seq data characteristics.

### Engineering Prompt

```
You are working on the WASP2 bioinformatics project to implement paired-end read simulation.

## Context
- Branch: sim/paired-end
- Working directory: /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp
- Conda environment: WASP2_dev2
- Reference: simulate_indel_ase_v2.py (existing single-end implementation)

## Problem
Current simulation generates single-end reads only, but WASP2 is designed for paired-end RNA-seq data. This is a critical gap for publication.

## Your Task
Create `simulate_paired_end_ase.py` that generates realistic paired-end reads with:
1. Proper insert size distribution (mean=300bp, std=50bp)
2. R1 and R2 FASTQ output files
3. Variants can appear in R1, R2, or both mates
4. Proper read pairing (same read name, /1 and /2 suffixes)

## Implementation Requirements

### 1. Read the existing implementation:
```bash
cat simulate_indel_ase_v2.py
```

### 2. Create new paired-end module:

File: `simulation/simulate_paired_end_ase.py`

```python
#!/usr/bin/env python3
"""
WASP2 Paired-End Simulation Framework

Generates paired-end reads with known allelic ratios for validation.
Key improvement over v2: proper R1/R2 paired-end reads.

Usage:
    python simulate_paired_end_ase.py --tier moderate --workdir output/
"""

import pysam
import random
import numpy as np
from dataclasses import dataclass
from typing import List, Tuple
from pathlib import Path

@dataclass
class PairedReadConfig:
    """Configuration for paired-end read generation."""
    insert_size_mean: int = 300
    insert_size_std: int = 50
    read_length: int = 150
    error_rate: float = 0.01

def generate_insert_size(config: PairedReadConfig) -> int:
    """Sample insert size from normal distribution."""
    size = int(np.random.normal(config.insert_size_mean, config.insert_size_std))
    return max(config.read_length + 50, size)  # Ensure valid insert

def reverse_complement(seq: str) -> str:
    """Return reverse complement of sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(b, 'N') for b in reversed(seq))

def create_paired_reads(
    ref_seq: str,
    var_pos: int,
    ref_allele: str,
    alt_allele: str,
    allele_to_use: str,  # 'ref' or 'alt'
    config: PairedReadConfig,
    read_id: str
) -> Tuple[Tuple[str, str], Tuple[str, str]]:
    """
    Generate paired-end reads with variant at specified position.

    Returns:
        ((r1_seq, r1_qual), (r2_seq, r2_qual))
    """
    insert_size = generate_insert_size(config)

    # Position the fragment so variant is covered
    # Random offset within the fragment
    var_offset = random.randint(20, insert_size - 20)
    frag_start = var_pos - var_offset
    frag_start = max(0, frag_start)

    # Get the allele to insert
    allele = ref_allele if allele_to_use == 'ref' else alt_allele

    # Build fragment sequence with allele
    left_seq = ref_seq[frag_start:var_pos]
    right_seq = ref_seq[var_pos + len(ref_allele):frag_start + insert_size]
    fragment = left_seq + allele + right_seq

    # Pad if needed
    if len(fragment) < insert_size:
        fragment += ref_seq[frag_start + len(fragment):frag_start + insert_size]

    # Extract R1 (forward) and R2 (reverse complement)
    r1_seq = fragment[:config.read_length]
    r2_seq = reverse_complement(fragment[-config.read_length:])

    # Add sequencing errors
    r1_seq = add_sequencing_errors(r1_seq, config.error_rate)
    r2_seq = add_sequencing_errors(r2_seq, config.error_rate)

    # Generate quality scores
    r1_qual = generate_quality_string(len(r1_seq))
    r2_qual = generate_quality_string(len(r2_seq))

    return ((r1_seq, r1_qual), (r2_seq, r2_qual))

def write_paired_fastq(
    r1_file,
    r2_file,
    read_id: str,
    r1_data: Tuple[str, str],
    r2_data: Tuple[str, str]
):
    """Write paired FASTQ records."""
    r1_seq, r1_qual = r1_data
    r2_seq, r2_qual = r2_data

    # R1
    r1_file.write(f"@{read_id}/1\n")
    r1_file.write(f"{r1_seq}\n")
    r1_file.write("+\n")
    r1_file.write(f"{r1_qual}\n")

    # R2
    r2_file.write(f"@{read_id}/2\n")
    r2_file.write(f"{r2_seq}\n")
    r2_file.write("+\n")
    r2_file.write(f"{r2_qual}\n")
```

### 3. Test scenarios to implement:

| Scenario | Description | Test Case |
|----------|-------------|-----------|
| Variant in R1 only | Fragment positioned so var in first 150bp | Common case |
| Variant in R2 only | Fragment positioned so var in last 150bp | Common case |
| Variant in both | Short insert, var near middle | Overlap case |
| Variant at R1 edge | Var at position 145-150 of R1 | Edge case |
| Variant at R2 edge | Var at position 1-5 of R2 | Edge case |

### 4. Integration with existing framework:

- Keep same GroundTruth dataclass
- Keep same test tiers (minimum/moderate/comprehensive)
- Keep same variant types (SNP/INS/DEL)
- Output both R1.fq and R2.fq files
- Update BWA alignment to use paired-end mode

### 5. BWA paired-end alignment:
```python
def align_paired_with_bwa(ref_fasta, r1_fastq, r2_fastq, output_bam, threads=4):
    """Align paired-end reads with BWA MEM."""
    subprocess.run([
        'bwa', 'mem',
        '-t', str(threads),
        ref_fasta,
        r1_fastq,
        r2_fastq
    ], stdout=open(output_bam.replace('.bam', '.sam'), 'w'), check=True)

    # Convert and sort
    pysam.view('-bS', '-o', output_bam, output_bam.replace('.bam', '.sam'))
    pysam.sort('-o', output_bam, output_bam)
    pysam.index(output_bam)
```

## Success Criteria
- [ ] Generates valid R1/R2 FASTQ pairs
- [ ] Insert size distribution matches config (mean=300, std=50)
- [ ] Variants appear in correct positions in reads
- [ ] BWA aligns pairs correctly (check proper pair flag)
- [ ] Runs through full WASP2 pipeline successfully
- [ ] Pass rate comparable to single-end (>90%)

## Testing
```bash
# Quick test
python simulation/simulate_paired_end_ase.py --tier minimum --workdir /tmp/pe_test

# Validate output
zcat /tmp/pe_test/synthetic_R1.fq.gz | head -8  # Check R1
zcat /tmp/pe_test/synthetic_R2.fq.gz | head -8  # Check R2

# Check read names match
paste <(zcat /tmp/pe_test/synthetic_R1.fq.gz | grep "^@") \
      <(zcat /tmp/pe_test/synthetic_R2.fq.gz | grep "^@") | head
```

## Commit
```bash
git add simulation/simulate_paired_end_ase.py
git commit -m "feat: add paired-end simulation module

- Generate proper R1/R2 FASTQ pairs
- Configurable insert size (default: 300±50bp)
- Variants can appear in R1, R2, or both
- BWA paired-end alignment integration
"
```
```

---

## AGENT 3: sim/gatk-compare (WRITE CODE - PARALLEL)

### Branch
```bash
git checkout sim/gatk-compare
```

### Priority
**P0 - CRITICAL** - Required for publication (competitor comparison)

### Task Description
Create a benchmark comparison module that runs the same simulated data through GATK ASEReadCounter and compares results with WASP2.

### Engineering Prompt

```
You are working on the WASP2 bioinformatics project to implement GATK ASEReadCounter comparison.

## Context
- Branch: sim/gatk-compare
- Working directory: /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp
- Conda environment: WASP2_dev2

## Problem
For publication, we need to compare WASP2 against GATK ASEReadCounter (the industry standard). Currently we have no competitor comparison.

## Your Task
Create `simulation/benchmark_vs_gatk.py` that:
1. Runs GATK ASEReadCounter on simulation data
2. Compares GATK counts vs WASP2 counts vs ground truth
3. Computes publication-quality comparison metrics

## Prerequisites Check
```bash
# Check GATK is available
which gatk || echo "GATK not found - need to install"

# If not installed:
conda install -c bioconda gatk4
```

## Implementation

### File: `simulation/benchmark_vs_gatk.py`

```python
#!/usr/bin/env python3
"""
WASP2 vs GATK ASEReadCounter Benchmark Comparison

Runs both tools on the same simulated data and compares:
1. Allele count accuracy (vs ground truth)
2. REF/ALT bias
3. Runtime performance

Usage:
    python benchmark_vs_gatk.py \
        --bam simulation_output/aligned.sorted.bam \
        --vcf simulation_output/variants.vcf.gz \
        --ref simulation_output/reference.fa \
        --ground-truth simulation_output/ground_truth.csv \
        --output comparison_results/
"""

import argparse
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import mean_squared_error, mean_absolute_error
import time
import json


def run_gatk_ase_counter(
    bam_file: str,
    vcf_file: str,
    ref_fasta: str,
    output_table: str,
    min_mapq: int = 10,
    min_baseq: int = 20
) -> float:
    """
    Run GATK ASEReadCounter and return runtime.

    Returns:
        Runtime in seconds
    """
    cmd = [
        "gatk", "ASEReadCounter",
        "-R", ref_fasta,
        "-I", bam_file,
        "-V", vcf_file,
        "-O", output_table,
        "--min-mapping-quality", str(min_mapq),
        "--min-base-quality", str(min_baseq),
    ]

    start = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True)
    runtime = time.time() - start

    if result.returncode != 0:
        print(f"GATK error: {result.stderr}")
        raise RuntimeError("GATK ASEReadCounter failed")

    return runtime


def parse_gatk_output(gatk_table: str) -> pd.DataFrame:
    """Parse GATK ASEReadCounter output table."""
    df = pd.read_csv(gatk_table, sep='\t', comment='#')

    # Standardize column names
    df = df.rename(columns={
        'contig': 'chrom',
        'position': 'pos',
        'refAllele': 'ref',
        'altAllele': 'alt',
        'refCount': 'gatk_ref_count',
        'altCount': 'gatk_alt_count',
        'totalCount': 'gatk_total_count'
    })

    return df[['chrom', 'pos', 'ref', 'alt', 'gatk_ref_count', 'gatk_alt_count', 'gatk_total_count']]


def run_wasp2_counting(
    bam_file: str,
    vcf_file: str,
    output_file: str
) -> Tuple[pd.DataFrame, float]:
    """
    Run WASP2 allele counting and return results + runtime.
    """
    # Import WASP2 counting module
    import sys
    sys.path.insert(0, 'src')
    from counting.count_alleles import count_alleles_at_sites

    start = time.time()
    counts = count_alleles_at_sites(bam_file, vcf_file, output_file)
    runtime = time.time() - start

    # Load results
    df = pd.read_csv(output_file, sep='\t')
    df = df.rename(columns={
        'ref_count': 'wasp2_ref_count',
        'alt_count': 'wasp2_alt_count'
    })

    return df, runtime


def compute_comparison_metrics(
    gatk_df: pd.DataFrame,
    wasp2_df: pd.DataFrame,
    ground_truth_df: pd.DataFrame
) -> dict:
    """
    Compute comprehensive comparison metrics.

    Returns dict with:
    - Correlation with ground truth (both tools)
    - RMSE (both tools)
    - REF bias (both tools)
    - Head-to-head comparison
    """
    # Merge all data
    merged = ground_truth_df.merge(gatk_df, on=['chrom', 'pos'])
    merged = merged.merge(wasp2_df, on=['chrom', 'pos'])

    # Ground truth ratio
    merged['true_ratio'] = merged['true_ref_count'] / (merged['true_ref_count'] + merged['true_alt_count'])

    # GATK ratio
    merged['gatk_ratio'] = merged['gatk_ref_count'] / (merged['gatk_ref_count'] + merged['gatk_alt_count'])

    # WASP2 ratio
    merged['wasp2_ratio'] = merged['wasp2_ref_count'] / (merged['wasp2_ref_count'] + merged['wasp2_alt_count'])

    # Handle NaN (division by zero)
    merged = merged.dropna()

    metrics = {
        # GATK vs truth
        'gatk_pearson_r': pearsonr(merged['gatk_ratio'], merged['true_ratio'])[0],
        'gatk_spearman_rho': spearmanr(merged['gatk_ratio'], merged['true_ratio'])[0],
        'gatk_rmse': np.sqrt(mean_squared_error(merged['true_ratio'], merged['gatk_ratio'])),
        'gatk_mae': mean_absolute_error(merged['true_ratio'], merged['gatk_ratio']),
        'gatk_ref_bias': merged['gatk_ratio'].mean() - 0.5,

        # WASP2 vs truth
        'wasp2_pearson_r': pearsonr(merged['wasp2_ratio'], merged['true_ratio'])[0],
        'wasp2_spearman_rho': spearmanr(merged['wasp2_ratio'], merged['true_ratio'])[0],
        'wasp2_rmse': np.sqrt(mean_squared_error(merged['true_ratio'], merged['wasp2_ratio'])),
        'wasp2_mae': mean_absolute_error(merged['true_ratio'], merged['wasp2_ratio']),
        'wasp2_ref_bias': merged['wasp2_ratio'].mean() - 0.5,

        # Head-to-head
        'wasp2_vs_gatk_pearson': pearsonr(merged['wasp2_ratio'], merged['gatk_ratio'])[0],
        'n_variants': len(merged),
    }

    return metrics, merged


def generate_comparison_table(metrics: dict) -> str:
    """Generate markdown comparison table."""
    table = """
| Metric | WASP2 | GATK | Winner |
|--------|-------|------|--------|
| Pearson r (vs truth) | {wasp2_pearson_r:.4f} | {gatk_pearson_r:.4f} | {pearson_winner} |
| Spearman ρ (vs truth) | {wasp2_spearman_rho:.4f} | {gatk_spearman_rho:.4f} | {spearman_winner} |
| RMSE | {wasp2_rmse:.4f} | {gatk_rmse:.4f} | {rmse_winner} |
| MAE | {wasp2_mae:.4f} | {gatk_mae:.4f} | {mae_winner} |
| REF bias | {wasp2_ref_bias:+.4f} | {gatk_ref_bias:+.4f} | {bias_winner} |
""".format(
        **metrics,
        pearson_winner='WASP2' if metrics['wasp2_pearson_r'] > metrics['gatk_pearson_r'] else 'GATK',
        spearman_winner='WASP2' if metrics['wasp2_spearman_rho'] > metrics['gatk_spearman_rho'] else 'GATK',
        rmse_winner='WASP2' if metrics['wasp2_rmse'] < metrics['gatk_rmse'] else 'GATK',
        mae_winner='WASP2' if metrics['wasp2_mae'] < metrics['gatk_mae'] else 'GATK',
        bias_winner='WASP2' if abs(metrics['wasp2_ref_bias']) < abs(metrics['gatk_ref_bias']) else 'GATK',
    )
    return table


def main():
    parser = argparse.ArgumentParser(description='WASP2 vs GATK comparison')
    parser.add_argument('--bam', required=True, help='Input BAM file')
    parser.add_argument('--vcf', required=True, help='Variant VCF file')
    parser.add_argument('--ref', required=True, help='Reference FASTA')
    parser.add_argument('--ground-truth', required=True, help='Ground truth CSV')
    parser.add_argument('--output', required=True, help='Output directory')

    args = parser.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load ground truth
    ground_truth = pd.read_csv(args.ground_truth)

    # Run GATK
    print("Running GATK ASEReadCounter...")
    gatk_output = output_dir / 'gatk_counts.table'
    gatk_runtime = run_gatk_ase_counter(args.bam, args.vcf, args.ref, str(gatk_output))
    gatk_df = parse_gatk_output(str(gatk_output))
    print(f"  GATK runtime: {gatk_runtime:.2f}s")

    # Run WASP2
    print("Running WASP2 counting...")
    wasp2_output = output_dir / 'wasp2_counts.tsv'
    wasp2_df, wasp2_runtime = run_wasp2_counting(args.bam, args.vcf, str(wasp2_output))
    print(f"  WASP2 runtime: {wasp2_runtime:.2f}s")

    # Compute metrics
    print("Computing comparison metrics...")
    metrics, merged = compute_comparison_metrics(gatk_df, wasp2_df, ground_truth)

    # Add runtime
    metrics['gatk_runtime_s'] = gatk_runtime
    metrics['wasp2_runtime_s'] = wasp2_runtime
    metrics['runtime_speedup'] = gatk_runtime / wasp2_runtime

    # Save results
    with open(output_dir / 'comparison_metrics.json', 'w') as f:
        json.dump(metrics, f, indent=2)

    merged.to_csv(output_dir / 'merged_counts.csv', index=False)

    # Print summary
    print("\n" + "="*60)
    print("COMPARISON RESULTS")
    print("="*60)
    print(generate_comparison_table(metrics))
    print(f"\nRuntime comparison:")
    print(f"  GATK:  {gatk_runtime:.2f}s")
    print(f"  WASP2: {wasp2_runtime:.2f}s")
    print(f"  Speedup: {metrics['runtime_speedup']:.2f}x")
    print(f"\nResults saved to: {output_dir}")


if __name__ == '__main__':
    main()
```

## Testing

1. First, ensure GATK is available:
```bash
gatk --version
```

2. Run on existing simulation output:
```bash
python simulation/benchmark_vs_gatk.py \
    --bam simulation_results/comprehensive_*/wasp2_output/keep.merged.bam \
    --vcf simulation_results/comprehensive_*/variants.vcf.gz \
    --ref simulation_results/comprehensive_*/reference.fa \
    --ground-truth simulation_results/comprehensive_*/ground_truth.csv \
    --output comparison_results/
```

## Success Criteria
- [ ] GATK runs successfully on simulation data
- [ ] Comparison metrics computed correctly
- [ ] WASP2 correlation ≥ GATK correlation
- [ ] WASP2 REF bias ≤ GATK REF bias
- [ ] Runtime comparison documented

## Commit
```bash
git add simulation/benchmark_vs_gatk.py
git commit -m "feat: add GATK ASEReadCounter comparison benchmark

- Run GATK on same simulated data as WASP2
- Compare correlation, RMSE, MAE, REF bias
- Generate publication-ready comparison table
- Track runtime performance
"
```
```

---

## AGENT 4: sim/metrics (AFTER DATA READY)

### Branch
```bash
git checkout sim/metrics
```

### Priority
**P1 - IMPORTANT** - Runs after Agents 1-3 complete

### Dependencies
- Agent 1 (comprehensive) output
- Agent 2 (paired-end) output
- Agent 3 (GATK comparison) output

### Task Description
Compute publication-quality metrics and generate figures for the paper.

### Engineering Prompt

```
You are working on the WASP2 bioinformatics project to compute publication metrics and generate figures.

## Context
- Branch: sim/metrics
- Working directory: /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp
- Conda environment: WASP2_dev2

## Dependencies
This task requires outputs from:
- sim/comprehensive: `simulation_results/comprehensive_*/`
- sim/paired-end: `simulation_results/paired_end_*/`
- sim/gatk-compare: `comparison_results/`

## Your Task
Create publication-quality metrics and figures for Nature Methods submission.

## Implementation

### File: `simulation/simulation_metrics.py`

```python
#!/usr/bin/env python3
"""
Publication Metrics for WASP2 Simulation Benchmarks

Computes all metrics needed for publication and generates figures.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr
from pathlib import Path


def compute_all_metrics(results_df: pd.DataFrame) -> dict:
    """Compute comprehensive metrics from simulation results."""

    metrics = {}

    # Overall accuracy
    metrics['n_tests'] = len(results_df)
    metrics['pass_rate'] = (results_df['status'] == 'PASS').mean() * 100
    metrics['mean_error_pct'] = results_df['error_pct'].mean()
    metrics['median_error_pct'] = results_df['error_pct'].median()
    metrics['std_error_pct'] = results_df['error_pct'].std()
    metrics['max_error_pct'] = results_df['error_pct'].max()

    # Correlation with ground truth
    if 'observed_ratio' in results_df.columns and 'true_ratio' in results_df.columns:
        valid = results_df[results_df['observed_ratio'] != np.inf]
        metrics['pearson_r'], metrics['pearson_p'] = pearsonr(
            valid['observed_ratio'], valid['true_ratio']
        )
        metrics['spearman_rho'], metrics['spearman_p'] = spearmanr(
            valid['observed_ratio'], valid['true_ratio']
        )

    # By variant type
    for vtype in ['SNP', 'INS', 'DEL']:
        subset = results_df[results_df['variant_type'] == vtype]
        if len(subset) > 0:
            metrics[f'{vtype.lower()}_n'] = len(subset)
            metrics[f'{vtype.lower()}_pass_rate'] = (subset['status'] == 'PASS').mean() * 100
            metrics[f'{vtype.lower()}_mean_error'] = subset['error_pct'].mean()
            metrics[f'{vtype.lower()}_std_error'] = subset['error_pct'].std()

    # By coverage (if available)
    if 'coverage' in results_df.columns:
        for cov in results_df['coverage'].unique():
            subset = results_df[results_df['coverage'] == cov]
            metrics[f'cov{cov}_pass_rate'] = (subset['status'] == 'PASS').mean() * 100
            metrics[f'cov{cov}_mean_error'] = subset['error_pct'].mean()

    # REF/ALT bias
    if 'ref_count' in results_df.columns and 'alt_count' in results_df.columns:
        total = results_df['ref_count'] + results_df['alt_count']
        ref_ratio = results_df['ref_count'] / total
        metrics['ref_bias'] = ref_ratio.mean() - 0.5
        metrics['ref_bias_std'] = ref_ratio.std()

    return metrics


def generate_figure_1(results_df: pd.DataFrame, output_path: Path):
    """
    Figure 1: Ground truth correlation scatter plot.

    Shows observed vs expected allelic ratio, colored by variant type.
    """
    fig, ax = plt.subplots(figsize=(8, 8))

    valid = results_df[results_df['observed_ratio'] != np.inf].copy()

    colors = {'SNP': '#2ecc71', 'INS': '#3498db', 'DEL': '#e74c3c'}

    for vtype in ['SNP', 'INS', 'DEL']:
        subset = valid[valid['variant_type'] == vtype]
        ax.scatter(
            subset['true_ratio'],
            subset['observed_ratio'],
            c=colors[vtype],
            label=f'{vtype} (n={len(subset)})',
            alpha=0.6,
            s=50
        )

    # Perfect correlation line
    max_val = max(valid['true_ratio'].max(), valid['observed_ratio'].max())
    ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.5, label='Perfect correlation')

    # Compute correlation
    r, p = pearsonr(valid['observed_ratio'], valid['true_ratio'])

    ax.set_xlabel('Expected Allelic Ratio (REF/ALT)', fontsize=12)
    ax.set_ylabel('Observed Allelic Ratio (REF/ALT)', fontsize=12)
    ax.set_title(f'WASP2 Allelic Ratio Accuracy\n(Pearson r = {r:.4f}, p = {p:.2e})', fontsize=14)
    ax.legend(loc='upper left')
    ax.set_aspect('equal')

    plt.tight_layout()
    plt.savefig(output_path / 'figure1_correlation.png', dpi=300, bbox_inches='tight')
    plt.savefig(output_path / 'figure1_correlation.pdf', bbox_inches='tight')
    plt.close()


def generate_figure_2(results_df: pd.DataFrame, output_path: Path):
    """
    Figure 2: Error distribution by variant type (box plot).
    """
    fig, ax = plt.subplots(figsize=(10, 6))

    order = ['SNP', 'INS', 'DEL']
    colors = ['#2ecc71', '#3498db', '#e74c3c']

    sns.boxplot(
        data=results_df,
        x='variant_type',
        y='error_pct',
        order=order,
        palette=colors,
        ax=ax
    )

    # Add individual points
    sns.stripplot(
        data=results_df,
        x='variant_type',
        y='error_pct',
        order=order,
        color='black',
        alpha=0.3,
        size=3,
        ax=ax
    )

    ax.axhline(y=10, color='red', linestyle='--', alpha=0.5, label='10% threshold')
    ax.set_xlabel('Variant Type', fontsize=12)
    ax.set_ylabel('Error (%)', fontsize=12)
    ax.set_title('WASP2 Error Distribution by Variant Type', fontsize=14)
    ax.legend()

    plt.tight_layout()
    plt.savefig(output_path / 'figure2_error_by_type.png', dpi=300, bbox_inches='tight')
    plt.savefig(output_path / 'figure2_error_by_type.pdf', bbox_inches='tight')
    plt.close()


def generate_figure_3(comparison_df: pd.DataFrame, output_path: Path):
    """
    Figure 3: WASP2 vs GATK comparison scatter plot.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # WASP2 vs truth
    ax1 = axes[0]
    ax1.scatter(comparison_df['true_ratio'], comparison_df['wasp2_ratio'], alpha=0.5, c='#3498db')
    r_wasp2, _ = pearsonr(comparison_df['wasp2_ratio'], comparison_df['true_ratio'])
    ax1.plot([0, 1], [0, 1], 'k--', alpha=0.5)
    ax1.set_xlabel('True Ratio')
    ax1.set_ylabel('WASP2 Ratio')
    ax1.set_title(f'WASP2 (r = {r_wasp2:.4f})')
    ax1.set_aspect('equal')

    # GATK vs truth
    ax2 = axes[1]
    ax2.scatter(comparison_df['true_ratio'], comparison_df['gatk_ratio'], alpha=0.5, c='#e74c3c')
    r_gatk, _ = pearsonr(comparison_df['gatk_ratio'], comparison_df['true_ratio'])
    ax2.plot([0, 1], [0, 1], 'k--', alpha=0.5)
    ax2.set_xlabel('True Ratio')
    ax2.set_ylabel('GATK Ratio')
    ax2.set_title(f'GATK (r = {r_gatk:.4f})')
    ax2.set_aspect('equal')

    plt.suptitle('Allelic Ratio Accuracy: WASP2 vs GATK ASEReadCounter', fontsize=14)
    plt.tight_layout()
    plt.savefig(output_path / 'figure3_wasp2_vs_gatk.png', dpi=300, bbox_inches='tight')
    plt.savefig(output_path / 'figure3_wasp2_vs_gatk.pdf', bbox_inches='tight')
    plt.close()


def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--results', required=True, help='Simulation results CSV')
    parser.add_argument('--comparison', help='GATK comparison CSV (optional)')
    parser.add_argument('--output', required=True, help='Output directory for figures')

    args = parser.parse_args()

    output_path = Path(args.output)
    output_path.mkdir(parents=True, exist_ok=True)

    # Load results
    results_df = pd.read_csv(args.results)

    # Compute metrics
    metrics = compute_all_metrics(results_df)

    # Save metrics
    import json
    with open(output_path / 'publication_metrics.json', 'w') as f:
        json.dump(metrics, f, indent=2)

    # Generate figures
    generate_figure_1(results_df, output_path)
    generate_figure_2(results_df, output_path)

    if args.comparison:
        comparison_df = pd.read_csv(args.comparison)
        generate_figure_3(comparison_df, output_path)

    print(f"Metrics and figures saved to: {output_path}")
    print(f"\nKey metrics:")
    print(f"  Pass rate: {metrics['pass_rate']:.1f}%")
    print(f"  Pearson r: {metrics.get('pearson_r', 'N/A')}")
    print(f"  Mean error: {metrics['mean_error_pct']:.2f}%")


if __name__ == '__main__':
    main()
```

## Success Criteria
- [ ] All metrics computed
- [ ] Figure 1: Correlation plot generated
- [ ] Figure 2: Error by variant type generated
- [ ] Figure 3: WASP2 vs GATK comparison generated
- [ ] Publication-quality (300 DPI, PDF + PNG)

## Commit
```bash
git add simulation/simulation_metrics.py
git add simulation/figures/
git commit -m "feat: add publication metrics and figures

- Compute comprehensive accuracy metrics
- Generate correlation scatter plot (Figure 1)
- Generate error distribution box plot (Figure 2)
- Generate WASP2 vs GATK comparison (Figure 3)
"
```
```

---

## Execution Order

```
TIME ──────────────────────────────────────────────────────────────────►

     ┌──────────────────┐
     │   AGENT 1        │
     │ sim/comprehensive│ ─────────────────────────────────┐
     │ (RUN - 2hrs)     │                                  │
     └──────────────────┘                                  │
                                                           │
     ┌──────────────────┐                                  │
     │   AGENT 2        │                                  │
     │ sim/paired-end   │ ──────────┐                      │
     │ (CODE - 2-3hrs)  │           │                      │
     └──────────────────┘           │                      │
                                    ▼                      │
                         ┌──────────────────┐              │
                         │ RUN PAIRED-END   │              │
                         │ (2hrs)           │ ─────────────┤
                         └──────────────────┘              │
                                                           │
     ┌──────────────────┐                                  │
     │   AGENT 3        │                                  │
     │ sim/gatk-compare │ ─────────────────────────────────┤
     │ (CODE - 2hrs)    │                                  │
     └──────────────────┘                                  │
                                                           │
                                                           ▼
                                              ┌──────────────────┐
                                              │   AGENT 4        │
                                              │ sim/metrics      │
                                              │ (AFTER all done) │
                                              └──────────────────┘
```

---

## Quick Reference: Branch Commands

```bash
# Agent 1: Run comprehensive (no code changes)
git checkout sim/comprehensive
# ... run simulation ...
git push origin sim/comprehensive

# Agent 2: Paired-end code
git checkout sim/paired-end
# ... write code ...
git push origin sim/paired-end

# Agent 3: GATK comparison
git checkout sim/gatk-compare
# ... write code ...
git push origin sim/gatk-compare

# Agent 4: Metrics (after merge)
git checkout sim/metrics
git merge sim/comprehensive
git merge sim/paired-end
git merge sim/gatk-compare
# ... compute metrics ...
git push origin sim/metrics

# Final merge to ropc-indels
git checkout ropc-indels
git merge sim/metrics
```
