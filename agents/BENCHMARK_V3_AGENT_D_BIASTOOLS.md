# Agent D: biastools Comparison

## Mission
Create a biastools comparison that uses its reference bias quantification to validate WASP2's bias correction effectiveness.

---

## Repository Context

**GitHub:** https://github.com/Jaureguy760/WASP2-exp.git
**Branch:** `sim/benchmark-v3`
**Working Directory:** `/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp`
**Conda Environment:** `WASP2_dev2`

---

## biastools Background

**Paper:** Günther & Lamparter, Bioinformatics 2023
**GitHub:** https://github.com/genesis-biomed/biastools
**PyPI:** `pip install biastools`

biastools provides:
1. Reference bias simulation
2. Bias quantification at heterozygous sites
3. Reference ratio analysis

### Key Insight
biastools was designed specifically to **measure reference bias** - exactly what WASP2 aims to correct. This makes it ideal for validating our correction.

---

## biastools API

### Installation
```bash
pip install biastools
# or
conda install -c bioconda biastools
```

### Mode 1: simulate
```bash
biastools simulate \
    --vcf variants.vcf.gz \
    --fasta reference.fa \
    --output simulated_reads.fastq \
    --coverage 30 \
    --read-length 100
```

### Mode 2: analyze
```bash
biastools analyze \
    --bam aligned.bam \
    --vcf variants.vcf.gz \
    --output bias_report.txt
```

### Mode 3: real
```bash
biastools real \
    --bam aligned.bam \
    --vcf variants.vcf.gz \
    --output real_bias.txt
```

### Output Format (analyze mode)
```
chrom  pos      ref  alt  ref_count  alt_count  total  ref_ratio  bias_status
chr1   50000    A    G    25         25         50     0.50       UNBIASED
chr1   60000    T    C    35         15         50     0.70       REF_BIASED
```

---

## Implementation

### File: `simulation/competitors/run_biastools.py`

```python
#!/usr/bin/env python3
"""
biastools comparison wrapper.

biastools quantifies reference bias at heterozygous sites.
We use it to validate WASP2's bias correction.

Key modes:
- simulate: Generate biased reads (we use our own simulator)
- analyze: Quantify bias from BAM (primary comparison method)
- real: Analyze real sequencing data
"""

import subprocess
import pandas as pd
import pysam
from pathlib import Path
from typing import Dict, Optional, Tuple
import numpy as np


def check_biastools_available() -> bool:
    """Check if biastools is installed."""
    try:
        result = subprocess.run(['biastools', '--help'], capture_output=True, text=True)
        return result.returncode == 0
    except FileNotFoundError:
        return False


def install_biastools():
    """Install biastools via pip."""
    print("Installing biastools...")
    subprocess.run(['pip', 'install', 'biastools'], check=True)


def run_biastools_analyze(
    bam_file: str,
    vcf_file: str,
    output_file: str,
    min_coverage: int = 10,
    min_mapq: int = 10,
    min_baseq: int = 20
) -> str:
    """
    Run biastools analyze on BAM file.

    Returns path to output file.
    """
    cmd = [
        'biastools', 'analyze',
        '--bam', bam_file,
        '--vcf', vcf_file,
        '--output', output_file,
        '--min-coverage', str(min_coverage),
        '--min-mapq', str(min_mapq),
        '--min-baseq', str(min_baseq)
    ]

    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        # biastools may not exist or fail - try Python fallback
        print(f"biastools CLI failed, using Python implementation...")
        return run_biastools_python_fallback(bam_file, vcf_file, output_file)

    return output_file


def run_biastools_python_fallback(
    bam_file: str,
    vcf_file: str,
    output_file: str
) -> str:
    """
    Python fallback if biastools CLI unavailable.

    Implements core bias quantification logic.
    """
    results = []

    # Load VCF
    vcf = pysam.VariantFile(vcf_file)
    bam = pysam.AlignmentFile(bam_file, 'rb')

    for rec in vcf:
        chrom = rec.chrom
        pos = rec.pos
        ref = rec.ref
        alt = rec.alts[0] if rec.alts else ''

        # Count alleles at this position
        ref_count = 0
        alt_count = 0

        # Classify variant type
        if len(ref) == 1 and len(alt) == 1:
            vtype = 'SNP'
            # For SNPs, count bases
            for pileup in bam.pileup(chrom, pos - 1, pos, truncate=True):
                if pileup.pos == pos - 1:  # 0-based
                    for read in pileup.pileups:
                        if read.is_del or read.is_refskip:
                            continue
                        base = read.alignment.query_sequence[read.query_position]
                        if base == ref:
                            ref_count += 1
                        elif base == alt:
                            alt_count += 1
        elif len(ref) < len(alt):
            vtype = 'INS'
            # For insertions, check CIGAR operations
            ref_count, alt_count = count_indel_alleles(bam, chrom, pos, ref, alt, 'INS')
        else:
            vtype = 'DEL'
            # For deletions, check CIGAR operations
            ref_count, alt_count = count_indel_alleles(bam, chrom, pos, ref, alt, 'DEL')

        total = ref_count + alt_count
        ref_ratio = ref_count / total if total > 0 else 0.5

        # Classify bias
        if total < 10:
            bias_status = 'LOW_COVERAGE'
        elif ref_ratio > 0.6:
            bias_status = 'REF_BIASED'
        elif ref_ratio < 0.4:
            bias_status = 'ALT_BIASED'
        else:
            bias_status = 'UNBIASED'

        results.append({
            'chrom': chrom,
            'pos': pos,
            'ref': ref,
            'alt': alt,
            'variant_type': vtype,
            'ref_count': ref_count,
            'alt_count': alt_count,
            'total': total,
            'ref_ratio': ref_ratio,
            'bias_status': bias_status
        })

    bam.close()
    vcf.close()

    # Save results
    df = pd.DataFrame(results)
    df.to_csv(output_file, sep='\t', index=False)

    return output_file


def count_indel_alleles(
    bam: pysam.AlignmentFile,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    vtype: str
) -> Tuple[int, int]:
    """
    Count INDEL alleles at a position.

    This is the core INDEL counting logic that WASP2 improves upon.
    """
    ref_count = 0
    alt_count = 0

    indel_size = abs(len(ref) - len(alt))

    for read in bam.fetch(chrom, pos - 1, pos + len(ref)):
        if read.is_unmapped or read.is_secondary:
            continue

        # Check if read spans the INDEL position
        if read.reference_start >= pos or read.reference_end <= pos:
            continue

        # Look for INDEL in CIGAR at this position
        read_pos = read.reference_start
        found_indel = False

        for op, length in read.cigartuples or []:
            if op == 0:  # M
                read_pos += length
            elif op == 1:  # I (insertion)
                if vtype == 'INS' and abs(read_pos - pos) <= 1 and abs(length - indel_size) <= 1:
                    alt_count += 1
                    found_indel = True
                    break
            elif op == 2:  # D (deletion)
                if vtype == 'DEL' and abs(read_pos - pos) <= 1 and abs(length - indel_size) <= 1:
                    alt_count += 1
                    found_indel = True
                    break
                read_pos += length
            elif op == 4:  # S
                pass  # soft clip doesn't consume reference
            elif op == 3:  # N
                read_pos += length

            if read_pos > pos + len(ref):
                break

        if not found_indel:
            # Read doesn't have the INDEL, count as REF
            ref_count += 1

    return ref_count, alt_count


def load_biastools_results(output_file: str) -> pd.DataFrame:
    """Load biastools output."""
    df = pd.read_csv(output_file, sep='\t')

    # Standardize column names
    df = df.rename(columns={
        'ref_count': 'biastools_ref_count',
        'alt_count': 'biastools_alt_count',
        'total': 'biastools_total',
        'ref_ratio': 'biastools_ratio',
        'bias_status': 'biastools_status'
    })

    return df


def calculate_bias_metrics(df: pd.DataFrame) -> Dict[str, float]:
    """
    Calculate bias metrics from biastools output.

    Key metrics:
    - Overall REF ratio (should be ~0.5 if unbiased)
    - Fraction of REF-biased sites
    - RMSE from 0.5 ratio
    """
    metrics = {}

    # Filter to adequate coverage
    df = df[df['biastools_total'] >= 10]

    if len(df) == 0:
        return {'error': 'No variants with adequate coverage'}

    # Overall mean REF ratio
    metrics['mean_ref_ratio'] = df['biastools_ratio'].mean()

    # Deviation from expected 0.5
    metrics['mean_bias'] = metrics['mean_ref_ratio'] - 0.5

    # RMSE from 0.5
    metrics['rmse_from_unbiased'] = np.sqrt(((df['biastools_ratio'] - 0.5) ** 2).mean())

    # Fraction biased
    n_ref_biased = (df['biastools_status'] == 'REF_BIASED').sum()
    n_alt_biased = (df['biastools_status'] == 'ALT_BIASED').sum()
    n_unbiased = (df['biastools_status'] == 'UNBIASED').sum()

    metrics['frac_ref_biased'] = n_ref_biased / len(df)
    metrics['frac_alt_biased'] = n_alt_biased / len(df)
    metrics['frac_unbiased'] = n_unbiased / len(df)

    # By variant type
    for vtype in ['SNP', 'INS', 'DEL']:
        vdf = df[df['variant_type'] == vtype]
        if len(vdf) > 0:
            metrics[f'{vtype.lower()}_mean_ratio'] = vdf['biastools_ratio'].mean()
            metrics[f'{vtype.lower()}_n_variants'] = len(vdf)

    return metrics


def compare_before_after_wasp2(
    before_bam: str,
    after_bam: str,
    vcf_file: str,
    output_dir: str
) -> pd.DataFrame:
    """
    Compare bias before and after WASP2 correction.

    This is the key validation: WASP2 should reduce reference bias.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Analyze before WASP2
    before_file = str(output_dir / 'biastools_before.tsv')
    run_biastools_analyze(before_bam, vcf_file, before_file)
    before_df = load_biastools_results(before_file)

    # Analyze after WASP2
    after_file = str(output_dir / 'biastools_after.tsv')
    run_biastools_analyze(after_bam, vcf_file, after_file)
    after_df = load_biastools_results(after_file)

    # Calculate metrics
    before_metrics = calculate_bias_metrics(before_df)
    after_metrics = calculate_bias_metrics(after_df)

    # Merge for comparison
    comparison = before_df.merge(
        after_df[['chrom', 'pos', 'biastools_ref_count', 'biastools_alt_count',
                  'biastools_ratio', 'biastools_status']],
        on=['chrom', 'pos'],
        how='outer',
        suffixes=('_before', '_after')
    )

    # Calculate improvement
    comparison['ratio_improvement'] = abs(comparison['biastools_ratio_before'] - 0.5) - \
                                      abs(comparison['biastools_ratio_after'] - 0.5)

    # Save results
    comparison.to_csv(output_dir / 'biastools_comparison.csv', index=False)

    # Save metrics
    metrics_df = pd.DataFrame([
        {'stage': 'before', **before_metrics},
        {'stage': 'after', **after_metrics}
    ])
    metrics_df.to_csv(output_dir / 'biastools_metrics.csv', index=False)

    return comparison


def run_full_biastools_comparison(
    bam_file: str,
    vcf_file: str,
    output_dir: str,
    ground_truth: Optional[pd.DataFrame] = None
) -> pd.DataFrame:
    """
    Run full biastools comparison pipeline.

    Returns DataFrame with bias quantification.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Check biastools
    if not check_biastools_available():
        install_biastools()

    # Run biastools
    output_file = str(output_dir / 'biastools_output.tsv')
    run_biastools_analyze(bam_file, vcf_file, output_file)

    # Load results
    results = load_biastools_results(output_file)

    # Calculate metrics
    metrics = calculate_bias_metrics(results)

    # Save metrics
    pd.DataFrame([metrics]).to_csv(output_dir / 'biastools_metrics.csv', index=False)

    # If ground truth provided, merge
    if ground_truth is not None:
        results = results.merge(
            ground_truth[['pos', 'true_ref_count', 'true_alt_count', 'true_ratio']],
            on='pos',
            how='left'
        )

    results.to_csv(output_dir / 'biastools_results.csv', index=False)

    return results


def main():
    import argparse

    parser = argparse.ArgumentParser(description='biastools comparison')
    parser.add_argument('--bam', required=True, help='Input BAM file')
    parser.add_argument('--vcf', required=True, help='VCF file')
    parser.add_argument('--output', '-o', required=True, help='Output directory')
    parser.add_argument('--ground-truth', help='Ground truth CSV (optional)')
    parser.add_argument('--before-bam', help='BAM before WASP2 correction (for comparison)')

    args = parser.parse_args()

    if args.before_bam:
        # Compare before/after
        results = compare_before_after_wasp2(
            before_bam=args.before_bam,
            after_bam=args.bam,
            vcf_file=args.vcf,
            output_dir=args.output
        )
        print(f"\nBias Comparison Summary:")
        print(f"  See {args.output}/biastools_metrics.csv for before/after metrics")
    else:
        ground_truth = None
        if args.ground_truth:
            ground_truth = pd.read_csv(args.ground_truth)

        results = run_full_biastools_comparison(
            bam_file=args.bam,
            vcf_file=args.vcf,
            output_dir=args.output,
            ground_truth=ground_truth
        )

        # Print summary
        metrics = calculate_bias_metrics(results)
        print(f"\nbiastools Results Summary:")
        print(f"  Total variants: {len(results)}")
        print(f"  Mean REF ratio: {metrics.get('mean_ref_ratio', 'N/A'):.3f}")
        print(f"  RMSE from unbiased: {metrics.get('rmse_from_unbiased', 'N/A'):.3f}")
        print(f"  REF-biased: {metrics.get('frac_ref_biased', 0)*100:.1f}%")
        print(f"  Unbiased: {metrics.get('frac_unbiased', 0)*100:.1f}%")


if __name__ == '__main__':
    main()
```

---

## Tasks

### Task 1: Create directory structure
```bash
mkdir -p simulation/competitors
```

### Task 2: Create `simulation/competitors/run_biastools.py`
Implement the script as specified above.

### Task 3: Test on simulation data
```bash
SIM_DIR="simulation_results/benchmark_v3"

python simulation/competitors/run_biastools.py \
    --bam ${SIM_DIR}/aligned.sorted.bam \
    --vcf ${SIM_DIR}/variants.vcf.gz \
    --output ${SIM_DIR}/biastools_comparison/ \
    --ground-truth ${SIM_DIR}/ground_truth.csv
```

### Task 4: Before/After WASP2 comparison
```bash
# Run WASP2 on aligned BAM first
python -m wasp2 \
    --bam ${SIM_DIR}/aligned.sorted.bam \
    --vcf ${SIM_DIR}/variants.vcf.gz \
    --output ${SIM_DIR}/wasp2_filtered.bam

# Compare before/after
python simulation/competitors/run_biastools.py \
    --before-bam ${SIM_DIR}/aligned.sorted.bam \
    --bam ${SIM_DIR}/wasp2_filtered.bam \
    --vcf ${SIM_DIR}/variants.vcf.gz \
    --output ${SIM_DIR}/wasp2_bias_improvement/
```

---

## Key Validation Points

1. **Bias Detection**: biastools should identify reference bias in uncorrected data
2. **WASP2 Improvement**: After WASP2, bias metrics should improve (ratio closer to 0.5)
3. **INDEL Performance**: biastools' INDEL counting is limited - highlight WASP2's advantage

### Expected Results

**Before WASP2:**
- Mean REF ratio: ~0.55-0.65 (biased toward reference)
- RMSE from 0.5: ~0.1-0.2
- REF-biased fraction: ~30-50%

**After WASP2:**
- Mean REF ratio: ~0.50-0.52 (close to unbiased)
- RMSE from 0.5: ~0.02-0.05
- REF-biased fraction: ~5-10%

---

## Success Criteria

- [ ] biastools or Python fallback runs without errors
- [ ] Bias quantified for all variant types
- [ ] Before/after comparison functional
- [ ] Metrics CSV generated
- [ ] Integration with ground truth works
- [ ] INDEL counting implemented (with known limitations)

---

## Commit Template

```bash
git add simulation/competitors/run_biastools.py
git commit -m "feat(sim): add biastools comparison with bias quantification

Validates WASP2's bias correction using biastools metrics:
- Quantifies reference bias at heterozygous sites
- Compares before/after WASP2 correction
- Python fallback if CLI unavailable
- INDEL counting with CIGAR analysis

Key validation: Show WASP2 reduces ref_ratio from ~0.6 to ~0.5.
"
```
