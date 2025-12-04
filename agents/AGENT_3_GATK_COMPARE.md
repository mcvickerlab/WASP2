# Agent 3: GATK ASEReadCounter Comparison

## Mission
Create a benchmark comparison that runs GATK ASEReadCounter on the same simulated data as WASP2, enabling head-to-head accuracy and performance comparison. This is **required for publication**.

---

## Repository Context

**GitHub:** https://github.com/Jaureguy760/WASP2-exp.git
**Branch:** `sim/gatk-compare`
**Parent Branch:** `ropc-indels`

**Working Directory:**
```
/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp
```

**Conda Environment:** `WASP2_dev2`

---

## Background: Why This Matters

For Nature Methods publication, reviewers will ask:
> "How does WASP2 compare to GATK ASEReadCounter, the industry standard?"

We need to show:
1. **Accuracy**: WASP2 correlation with ground truth ≥ GATK
2. **Bias**: WASP2 REF bias ≤ GATK REF bias
3. **Speed**: WASP2 runtime vs GATK runtime
4. **Concordance**: Agreement between tools on same data

---

## GATK ASEReadCounter Reference

### What It Does
GATK ASEReadCounter counts reads supporting REF vs ALT alleles at heterozygous sites. It's the standard tool for allele-specific expression analysis.

### Installation Check
```bash
# Check if GATK is installed
which gatk
gatk --version

# If not installed:
conda install -c bioconda gatk4
```

### Command Syntax
```bash
gatk ASEReadCounter \
    -R reference.fa \           # Reference genome (required)
    -I input.bam \              # Input BAM file (required)
    -V variants.vcf \           # VCF with het sites (required)
    -O output.table \           # Output table (required)
    --min-mapping-quality 10 \  # Min MAPQ (default: 10)
    --min-base-quality 20       # Min base quality (default: 20)
```

### Output Format
Tab-separated table with columns:
```
contig  position  variantID  refAllele  altAllele  refCount  altCount  totalCount  lowMAPQDepth  lowBaseQDepth  rawDepth  otherBases  improperPairs
chr1    50000     .          A          G          25        24        49          0             0              49        0           0
```

### Key Columns for Comparison:
- `refCount`: Reads supporting REF allele
- `altCount`: Reads supporting ALT allele
- `totalCount`: Total reads at site

---

## WASP2 Counting API Reference

### File: src/counting/count_alleles.py

```python
from wasp2_rust import BamCounter as RustBamCounter

def count_snp_alleles_rust(bam_file, chrom, snp_list, threads=1):
    """
    Rust-accelerated allele counting.

    Args:
        bam_file: Path to BAM file
        chrom: Chromosome name
        snp_list: Iterator of (pos, ref, alt) tuples
        threads: Number of threads

    Returns:
        List of (chrom, pos, ref_count, alt_count, other_count) tuples
    """
    regions = [(chrom, pos, ref, alt) for pos, ref, alt in snp_list]
    counter = RustBamCounter(bam_file)
    counts = counter.count_alleles(regions, min_qual=0, threads=threads)
    return [
        (chrom, pos, ref_count, alt_count, other_count)
        for (_, pos, _, _), (ref_count, alt_count, other_count) in zip(regions, counts)
    ]
```

---

## Implementation Specification

### New File: `simulation/benchmark_vs_gatk.py`

```python
#!/usr/bin/env python3
"""
WASP2 vs GATK ASEReadCounter Benchmark Comparison

Runs both tools on simulated data with known ground truth.
Computes accuracy metrics and generates comparison table.

Usage:
    python benchmark_vs_gatk.py \
        --bam simulation_output/aligned.sorted.bam \
        --vcf simulation_output/variants.vcf.gz \
        --ref simulation_output/reference.fa \
        --ground-truth simulation_output/simulation_results.csv \
        --output comparison_results/
"""

import argparse
import subprocess
import json
import time
from pathlib import Path
from typing import Tuple, Dict

import pandas as pd
import numpy as np
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import mean_squared_error, mean_absolute_error


def check_gatk_available() -> bool:
    """Check if GATK is installed and accessible."""
    try:
        result = subprocess.run(['gatk', '--version'],
                                capture_output=True, text=True)
        return result.returncode == 0
    except FileNotFoundError:
        return False


def run_gatk_ase_counter(
    bam_file: str,
    vcf_file: str,
    ref_fasta: str,
    output_table: str,
    min_mapq: int = 10,
    min_baseq: int = 20
) -> Tuple[float, pd.DataFrame]:
    """
    Run GATK ASEReadCounter.

    Args:
        bam_file: Input BAM (must be indexed)
        vcf_file: VCF with variant sites
        ref_fasta: Reference FASTA (must be indexed with .fai and .dict)
        output_table: Output path for counts table
        min_mapq: Minimum mapping quality
        min_baseq: Minimum base quality

    Returns:
        (runtime_seconds, counts_dataframe)
    """
    # Ensure reference has dictionary
    dict_file = ref_fasta.replace('.fa', '.dict')
    if not Path(dict_file).exists():
        print("Creating sequence dictionary...")
        subprocess.run([
            'gatk', 'CreateSequenceDictionary',
            '-R', ref_fasta,
            '-O', dict_file
        ], check=True)

    # Ensure BAM is indexed
    if not Path(f"{bam_file}.bai").exists():
        subprocess.run(['samtools', 'index', bam_file], check=True)

    # Run ASEReadCounter
    cmd = [
        'gatk', 'ASEReadCounter',
        '-R', ref_fasta,
        '-I', bam_file,
        '-V', vcf_file,
        '-O', output_table,
        '--min-mapping-quality', str(min_mapq),
        '--min-base-quality', str(min_baseq),
    ]

    print(f"Running: {' '.join(cmd)}")
    start_time = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True)
    runtime = time.time() - start_time

    if result.returncode != 0:
        print(f"GATK stderr: {result.stderr}")
        raise RuntimeError(f"GATK ASEReadCounter failed: {result.stderr}")

    # Parse output
    df = pd.read_csv(output_table, sep='\t', comment='#')
    df = df.rename(columns={
        'contig': 'chrom',
        'position': 'pos',
        'refAllele': 'ref',
        'altAllele': 'alt',
        'refCount': 'gatk_ref',
        'altCount': 'gatk_alt',
        'totalCount': 'gatk_total'
    })

    return runtime, df[['chrom', 'pos', 'ref', 'alt', 'gatk_ref', 'gatk_alt', 'gatk_total']]


def run_wasp2_counting(
    bam_file: str,
    vcf_file: str
) -> Tuple[float, pd.DataFrame]:
    """
    Run WASP2 allele counting.

    Returns:
        (runtime_seconds, counts_dataframe)
    """
    import sys
    sys.path.insert(0, 'src')
    from wasp2_rust import BamCounter

    # Parse VCF to get sites
    import pysam
    vcf = pysam.VariantFile(vcf_file)
    sites = []
    for record in vcf:
        sites.append({
            'chrom': record.chrom,
            'pos': record.pos,
            'ref': record.ref,
            'alt': record.alts[0] if record.alts else '.'
        })
    vcf.close()

    # Count alleles
    start_time = time.time()
    counter = BamCounter(bam_file)

    results = []
    for site in sites:
        regions = [(site['chrom'], site['pos'], site['ref'], site['alt'])]
        counts = counter.count_alleles(regions, min_qual=0, threads=1)
        if counts:
            ref_count, alt_count, other_count = counts[0]
            results.append({
                'chrom': site['chrom'],
                'pos': site['pos'],
                'ref': site['ref'],
                'alt': site['alt'],
                'wasp2_ref': ref_count,
                'wasp2_alt': alt_count,
                'wasp2_total': ref_count + alt_count
            })

    runtime = time.time() - start_time
    df = pd.DataFrame(results)

    return runtime, df


def load_ground_truth(csv_path: str) -> pd.DataFrame:
    """Load ground truth from simulation results."""
    df = pd.read_csv(csv_path)

    # Aggregate by position (remove replicates)
    agg_df = df.groupby(['chrom', 'pos']).agg({
        'ref_count': 'sum',
        'alt_count': 'sum',
        'true_ratio': 'first',
        'variant_type': 'first'
    }).reset_index()

    agg_df = agg_df.rename(columns={
        'ref_count': 'truth_ref',
        'alt_count': 'truth_alt'
    })
    agg_df['truth_total'] = agg_df['truth_ref'] + agg_df['truth_alt']

    return agg_df


def compute_metrics(
    merged_df: pd.DataFrame,
    tool_prefix: str
) -> Dict[str, float]:
    """
    Compute accuracy metrics for a tool.

    Args:
        merged_df: DataFrame with truth and tool counts
        tool_prefix: Column prefix ('gatk' or 'wasp2')

    Returns:
        Dictionary of metrics
    """
    # Filter out zeros to avoid division issues
    df = merged_df[
        (merged_df['truth_total'] > 0) &
        (merged_df[f'{tool_prefix}_total'] > 0)
    ].copy()

    # Compute ratios
    df['truth_ratio'] = df['truth_ref'] / df['truth_total']
    df[f'{tool_prefix}_ratio'] = df[f'{tool_prefix}_ref'] / df[f'{tool_prefix}_total']

    metrics = {
        'n_sites': len(df),

        # Correlation with truth
        'pearson_r': pearsonr(df[f'{tool_prefix}_ratio'], df['truth_ratio'])[0],
        'pearson_p': pearsonr(df[f'{tool_prefix}_ratio'], df['truth_ratio'])[1],
        'spearman_rho': spearmanr(df[f'{tool_prefix}_ratio'], df['truth_ratio'])[0],

        # Error metrics
        'rmse': np.sqrt(mean_squared_error(df['truth_ratio'], df[f'{tool_prefix}_ratio'])),
        'mae': mean_absolute_error(df['truth_ratio'], df[f'{tool_prefix}_ratio']),

        # Bias
        'ref_bias': df[f'{tool_prefix}_ratio'].mean() - 0.5,
        'ref_bias_std': df[f'{tool_prefix}_ratio'].std(),

        # Count correlation
        'ref_count_r': pearsonr(df[f'{tool_prefix}_ref'], df['truth_ref'])[0],
        'alt_count_r': pearsonr(df[f'{tool_prefix}_alt'], df['truth_alt'])[0],
    }

    return metrics


def generate_comparison_table(
    gatk_metrics: Dict,
    wasp2_metrics: Dict
) -> str:
    """Generate markdown comparison table."""

    def winner(metric, gatk_val, wasp2_val, lower_better=False):
        if lower_better:
            return 'WASP2' if abs(wasp2_val) < abs(gatk_val) else 'GATK'
        else:
            return 'WASP2' if wasp2_val > gatk_val else 'GATK'

    table = f"""
## Accuracy Comparison: WASP2 vs GATK ASEReadCounter

| Metric | WASP2 | GATK | Better |
|--------|-------|------|--------|
| Pearson r (vs truth) | {wasp2_metrics['pearson_r']:.4f} | {gatk_metrics['pearson_r']:.4f} | {winner('r', gatk_metrics['pearson_r'], wasp2_metrics['pearson_r'])} |
| Spearman ρ (vs truth) | {wasp2_metrics['spearman_rho']:.4f} | {gatk_metrics['spearman_rho']:.4f} | {winner('rho', gatk_metrics['spearman_rho'], wasp2_metrics['spearman_rho'])} |
| RMSE | {wasp2_metrics['rmse']:.4f} | {gatk_metrics['rmse']:.4f} | {winner('rmse', gatk_metrics['rmse'], wasp2_metrics['rmse'], lower_better=True)} |
| MAE | {wasp2_metrics['mae']:.4f} | {gatk_metrics['mae']:.4f} | {winner('mae', gatk_metrics['mae'], wasp2_metrics['mae'], lower_better=True)} |
| REF bias | {wasp2_metrics['ref_bias']:+.4f} | {gatk_metrics['ref_bias']:+.4f} | {winner('bias', gatk_metrics['ref_bias'], wasp2_metrics['ref_bias'], lower_better=True)} |
| Sites compared | {wasp2_metrics['n_sites']} | {gatk_metrics['n_sites']} | - |
"""
    return table


def main():
    parser = argparse.ArgumentParser(
        description='Compare WASP2 vs GATK ASEReadCounter'
    )
    parser.add_argument('--bam', required=True, help='Input BAM file')
    parser.add_argument('--vcf', required=True, help='Variant VCF file')
    parser.add_argument('--ref', required=True, help='Reference FASTA')
    parser.add_argument('--ground-truth', required=True,
                        help='Ground truth CSV from simulation')
    parser.add_argument('--output', required=True, help='Output directory')

    args = parser.parse_args()

    # Check GATK
    if not check_gatk_available():
        print("ERROR: GATK not found. Install with: conda install -c bioconda gatk4")
        return 1

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("="*60)
    print("WASP2 vs GATK ASEReadCounter Benchmark")
    print("="*60)

    # Load ground truth
    print("\nLoading ground truth...")
    truth_df = load_ground_truth(args.ground_truth)
    print(f"  {len(truth_df)} variant sites")

    # Run GATK
    print("\nRunning GATK ASEReadCounter...")
    gatk_runtime, gatk_df = run_gatk_ase_counter(
        args.bam, args.vcf, args.ref,
        str(output_dir / 'gatk_counts.table')
    )
    print(f"  GATK runtime: {gatk_runtime:.2f}s")
    print(f"  Sites counted: {len(gatk_df)}")

    # Run WASP2
    print("\nRunning WASP2 counting...")
    wasp2_runtime, wasp2_df = run_wasp2_counting(args.bam, args.vcf)
    print(f"  WASP2 runtime: {wasp2_runtime:.2f}s")
    print(f"  Sites counted: {len(wasp2_df)}")

    # Merge all data
    print("\nMerging results...")
    merged = truth_df.merge(gatk_df, on=['chrom', 'pos'], how='inner')
    merged = merged.merge(wasp2_df, on=['chrom', 'pos'], how='inner')
    print(f"  Common sites: {len(merged)}")

    # Compute metrics
    print("\nComputing metrics...")
    gatk_metrics = compute_metrics(merged, 'gatk')
    wasp2_metrics = compute_metrics(merged, 'wasp2')

    # Add runtime
    gatk_metrics['runtime_s'] = gatk_runtime
    wasp2_metrics['runtime_s'] = wasp2_runtime

    # Save results
    merged.to_csv(output_dir / 'merged_counts.csv', index=False)

    with open(output_dir / 'gatk_metrics.json', 'w') as f:
        json.dump(gatk_metrics, f, indent=2)

    with open(output_dir / 'wasp2_metrics.json', 'w') as f:
        json.dump(wasp2_metrics, f, indent=2)

    # Generate report
    report = generate_comparison_table(gatk_metrics, wasp2_metrics)
    report += f"""

## Runtime Comparison

| Tool | Runtime (s) | Speedup |
|------|-------------|---------|
| GATK | {gatk_runtime:.2f} | 1.0x |
| WASP2 | {wasp2_runtime:.2f} | {gatk_runtime/wasp2_runtime:.1f}x |
"""

    with open(output_dir / 'comparison_report.md', 'w') as f:
        f.write(report)

    print("\n" + report)
    print(f"\nResults saved to: {output_dir}")

    return 0


if __name__ == '__main__':
    exit(main())
```

---

## Step-by-Step Execution

### Step 1: Setup
```bash
source /iblm/netapp/home/jjaureguy/mambaforge/etc/profile.d/conda.sh
conda activate WASP2_dev2
cd /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp

git checkout sim/gatk-compare
git pull origin sim/gatk-compare

mkdir -p simulation
```

### Step 2: Check/Install GATK
```bash
# Check if installed
gatk --version

# If not, install:
conda install -c bioconda gatk4 -y

# Verify
gatk ASEReadCounter --help | head -20
```

### Step 3: Create the Module
```bash
# Create simulation/benchmark_vs_gatk.py using the code above
```

### Step 4: Test on Existing Simulation Data

Wait for Agent 1 (sim/comprehensive) to complete, then:

```bash
# Find the simulation output
SIM_DIR=$(ls -td simulation_results/comprehensive_* | head -1)
echo "Using: ${SIM_DIR}"

# Run comparison
python simulation/benchmark_vs_gatk.py \
    --bam ${SIM_DIR}/aligned.sorted.bam \
    --vcf ${SIM_DIR}/variants.vcf.gz \
    --ref ${SIM_DIR}/reference.fa \
    --ground-truth ${SIM_DIR}/simulation_results.csv \
    --output comparison_results/
```

### Step 5: Validate Results
```bash
# Check output files
ls comparison_results/

# View comparison report
cat comparison_results/comparison_report.md

# Check metrics
cat comparison_results/wasp2_metrics.json
cat comparison_results/gatk_metrics.json
```

---

## Expected Output

### comparison_report.md
```markdown
## Accuracy Comparison: WASP2 vs GATK ASEReadCounter

| Metric | WASP2 | GATK | Better |
|--------|-------|------|--------|
| Pearson r (vs truth) | 0.9823 | 0.9801 | WASP2 |
| Spearman ρ (vs truth) | 0.9756 | 0.9734 | WASP2 |
| RMSE | 0.0234 | 0.0256 | WASP2 |
| MAE | 0.0189 | 0.0203 | WASP2 |
| REF bias | +0.0012 | +0.0034 | WASP2 |

## Runtime Comparison

| Tool | Runtime (s) | Speedup |
|------|-------------|---------|
| GATK | 45.23 | 1.0x |
| WASP2 | 3.12 | 14.5x |
```

---

## Success Criteria

- [ ] GATK ASEReadCounter runs successfully
- [ ] WASP2 counting runs successfully
- [ ] Metrics computed correctly
- [ ] WASP2 Pearson r ≥ 0.95 (comparable to GATK)
- [ ] WASP2 REF bias ≤ GATK REF bias
- [ ] Comparison report generated
- [ ] Code committed to sim/gatk-compare

---

## Commit Template

```bash
git add simulation/benchmark_vs_gatk.py
git add comparison_results/comparison_report.md
git add comparison_results/*.json

git commit -m "feat: add GATK ASEReadCounter comparison benchmark

Compares WASP2 vs GATK on simulated data with ground truth.

Results:
- WASP2 Pearson r: X.XXXX
- GATK Pearson r: X.XXXX
- WASP2 REF bias: +X.XXXX
- GATK REF bias: +X.XXXX
- WASP2 speedup: X.Xx faster

Conclusion: WASP2 achieves comparable/better accuracy with Xx speedup.
"

git push origin sim/gatk-compare
```

---

## Troubleshooting

### GATK Errors:

**"Sequence dictionary required":**
```bash
gatk CreateSequenceDictionary -R reference.fa
```

**"BAM not indexed":**
```bash
samtools index input.bam
```

**"VCF not indexed":**
```bash
gatk IndexFeatureFile -I variants.vcf.gz
```

---

## Handoff

When complete:
1. Push to `sim/gatk-compare`
2. Document WASP2 vs GATK metrics in commit message
3. Ready for Agent 4 to generate publication figures
