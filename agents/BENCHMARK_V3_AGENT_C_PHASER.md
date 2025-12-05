# Agent C: phASER Comparison

## Mission
Create a phASER comparison that properly handles haplotype-level ASE and validates WASP2's approach against this production-grade tool.

---

## Repository Context

**GitHub:** https://github.com/Jaureguy760/WASP2-exp.git
**Branch:** `sim/benchmark-v3`
**Working Directory:** `/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp`
**Conda Environment:** `WASP2_dev2`

---

## phASER Background

**Paper:** Castel et al., Nature Communications 2016
**GitHub:** https://github.com/secastel/phaser

phASER performs:
1. Haplotype phasing from RNA-seq (read-backed phasing)
2. Allele-specific expression quantification
3. Gene-level ASE aggregation

### Key Insight
phASER is designed for **real RNA-seq data** where phasing is unknown. Our simulation has **known phasing**, so we can:
1. Provide phased VCF (saves computation)
2. Compare allele counts directly
3. Evaluate accuracy against ground truth

---

## phASER API

### Main Command: `phaser.py`
```bash
python phaser.py \
    --bam input.bam \
    --vcf variants.vcf.gz \
    --sample SAMPLE_ID \
    --baseq 10 \
    --mapq 10 \
    --o output_prefix \
    --threads 4
```

### Output Files
- `{prefix}.allelic_counts.txt` - Per-variant allele counts
- `{prefix}.haplotypic_counts.txt` - Haplotype-level counts
- `{prefix}.vcf.gz` - Phased VCF

### Allelic Counts Format
```
contig  position  variantID  refAllele  altAllele  refCount  altCount  totalCount  cA  cC  cG  cT  cN  refBias
chr1    50000     snp1       A          G          25        25        50          0   0   25  25  0   0.0
```

### Gene-Level: `phaser_gene_ae.py`
```bash
python phaser_gene_ae.py \
    --haplotypic_counts haplotypic_counts.txt \
    --features genes.bed \
    --o gene_ae.txt
```

---

## Implementation

### File: `simulation/competitors/run_phaser.py`

```python
#!/usr/bin/env python3
"""
phASER comparison wrapper.

phASER provides haplotype-aware ASE counting.
We compare its allelic counts against our ground truth.
"""

import subprocess
import pandas as pd
import pysam
from pathlib import Path
from typing import Dict, Optional
import tempfile
import os


def check_phaser_available() -> bool:
    """Check if phASER is installed."""
    try:
        # phASER is a Python script, check if importable
        result = subprocess.run(
            ['python', '-c', 'import phaser'],
            capture_output=True, text=True
        )
        if result.returncode == 0:
            return True

        # Try direct script
        phaser_path = os.environ.get('PHASER_PATH', '')
        if phaser_path and Path(phaser_path).exists():
            return True

        return False
    except Exception:
        return False


def install_phaser():
    """Install phASER from GitHub."""
    print("Installing phASER...")

    # Clone repository
    install_dir = Path.home() / '.local' / 'phaser'
    if not install_dir.exists():
        subprocess.run([
            'git', 'clone',
            'https://github.com/secastel/phaser.git',
            str(install_dir)
        ], check=True)

    # Set environment variable
    os.environ['PHASER_PATH'] = str(install_dir / 'phaser' / 'phaser.py')

    print(f"phASER installed to {install_dir}")
    print("Note: phASER requires IntervalTree: pip install intervaltree")
    subprocess.run(['pip', 'install', 'intervaltree'], check=True)


def get_phaser_script() -> str:
    """Get path to phaser.py script."""
    phaser_path = os.environ.get('PHASER_PATH', '')
    if phaser_path and Path(phaser_path).exists():
        return phaser_path

    # Default location
    default = Path.home() / '.local' / 'phaser' / 'phaser' / 'phaser.py'
    if default.exists():
        return str(default)

    raise RuntimeError("phASER not found. Set PHASER_PATH environment variable.")


def prepare_phased_vcf(input_vcf: str, output_vcf: str, sample_name: str = "SIMULATED"):
    """
    Prepare VCF for phASER with proper phasing.

    Our simulation VCFs have ground truth phasing.
    REF = HAP1, ALT = HAP2 (encoded in read names)
    """
    vcf_in = pysam.VariantFile(input_vcf)

    # Create new header with sample
    new_header = vcf_in.header.copy()
    if sample_name not in new_header.samples:
        new_header.add_sample(sample_name)

    # Add FORMAT fields if missing
    if 'GT' not in new_header.formats:
        new_header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')

    vcf_out = pysam.VariantFile(output_vcf, 'w', header=new_header)

    for rec in vcf_in:
        new_rec = vcf_out.new_record()
        new_rec.chrom = rec.chrom
        new_rec.pos = rec.pos
        new_rec.id = rec.id
        new_rec.alleles = rec.alleles
        new_rec.qual = rec.qual

        # Set phased genotype (0|1 means REF on hap1, ALT on hap2)
        new_rec.samples[sample_name]['GT'] = (0, 1)
        new_rec.samples[sample_name].phased = True

        vcf_out.write(new_rec)

    vcf_in.close()
    vcf_out.close()

    # Index
    pysam.tabix_index(output_vcf, preset='vcf', force=True)


def run_phaser(
    bam_file: str,
    vcf_file: str,
    output_prefix: str,
    sample_name: str = "SIMULATED",
    min_baseq: int = 10,
    min_mapq: int = 10,
    threads: int = 4
) -> str:
    """
    Run phASER on BAM file.

    Returns path to allelic counts file.
    """
    phaser_script = get_phaser_script()

    # Prepare phased VCF
    phased_vcf = output_prefix + '.phased.vcf.gz'
    prepare_phased_vcf(vcf_file, phased_vcf, sample_name)

    cmd = [
        'python', phaser_script,
        '--bam', bam_file,
        '--vcf', phased_vcf,
        '--sample', sample_name,
        '--baseq', str(min_baseq),
        '--mapq', str(min_mapq),
        '--o', output_prefix,
        '--threads', str(threads),
        '--pass_only'  # Only use PASS variants
    ]

    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"phASER stderr: {result.stderr}")
        # phASER may succeed but print warnings
        if not Path(f"{output_prefix}.allelic_counts.txt").exists():
            raise RuntimeError(f"phASER failed: {result.stderr}")

    return f"{output_prefix}.allelic_counts.txt"


def load_phaser_counts(phaser_output: str) -> pd.DataFrame:
    """Load phASER allelic counts."""
    df = pd.read_csv(phaser_output, sep='\t')

    # Standardize column names
    df = df.rename(columns={
        'contig': 'chrom',
        'position': 'pos',
        'variantID': 'variant_id',
        'refAllele': 'ref',
        'altAllele': 'alt',
        'refCount': 'phaser_ref_count',
        'altCount': 'phaser_alt_count',
        'totalCount': 'phaser_total'
    })

    # Calculate ratio
    df['phaser_ratio'] = df['phaser_ref_count'] / df['phaser_total']
    df['phaser_ratio'] = df['phaser_ratio'].fillna(0.5)  # Handle zero coverage

    return df


def load_vcf_variants(vcf_file: str) -> pd.DataFrame:
    """Load variants from VCF for comparison."""
    variants = []
    vcf = pysam.VariantFile(vcf_file)

    for rec in vcf:
        ref = rec.ref
        alt = rec.alts[0] if rec.alts else ''

        # Classify variant
        if len(ref) == 1 and len(alt) == 1:
            vtype = 'SNP'
        elif len(ref) < len(alt):
            vtype = 'INS'
        else:
            vtype = 'DEL'

        variants.append({
            'chrom': rec.chrom,
            'pos': rec.pos,
            'ref': ref,
            'alt': alt,
            'variant_type': vtype,
            'variant_id': rec.id or f"{rec.chrom}_{rec.pos}"
        })

    return pd.DataFrame(variants)


def compare_phaser_to_ground_truth(
    phaser_counts: pd.DataFrame,
    vcf_file: str,
    ground_truth: Optional[pd.DataFrame] = None
) -> pd.DataFrame:
    """
    Compare phASER counts to VCF variants and optionally ground truth.

    Returns merged DataFrame with comparison metrics.
    """
    vcf_variants = load_vcf_variants(vcf_file)

    # Merge phASER counts with VCF variants
    merged = vcf_variants.merge(
        phaser_counts[['chrom', 'pos', 'phaser_ref_count', 'phaser_alt_count',
                       'phaser_total', 'phaser_ratio']],
        on=['chrom', 'pos'],
        how='left'
    )

    # Fill missing values (phASER couldn't process)
    merged['phaser_ref_count'] = merged['phaser_ref_count'].fillna(0)
    merged['phaser_alt_count'] = merged['phaser_alt_count'].fillna(0)
    merged['phaser_total'] = merged['phaser_total'].fillna(0)
    merged['phaser_ratio'] = merged['phaser_ratio'].fillna(0.5)

    # Add status
    merged['phaser_status'] = 'OK'
    merged.loc[merged['phaser_total'] == 0, 'phaser_status'] = 'MISSING'

    # If ground truth provided, merge it
    if ground_truth is not None:
        merged = merged.merge(
            ground_truth[['pos', 'true_ref_count', 'true_alt_count', 'true_ratio']],
            on='pos',
            how='left'
        )

    return merged


def run_full_phaser_comparison(
    bam_file: str,
    vcf_file: str,
    output_dir: str,
    ground_truth: Optional[pd.DataFrame] = None
) -> pd.DataFrame:
    """
    Run full phASER comparison pipeline.

    Returns DataFrame with variant-level comparison.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Check phASER
    if not check_phaser_available():
        install_phaser()

    # Run phASER
    output_prefix = str(output_dir / 'phaser')
    phaser_file = run_phaser(bam_file, vcf_file, output_prefix)

    # Load results
    phaser_counts = load_phaser_counts(phaser_file)

    # Compare to ground truth
    comparison = compare_phaser_to_ground_truth(phaser_counts, vcf_file, ground_truth)

    # Save
    comparison.to_csv(output_dir / 'phaser_comparison.csv', index=False)

    return comparison


def main():
    import argparse

    parser = argparse.ArgumentParser(description='phASER comparison')
    parser.add_argument('--bam', required=True, help='Input BAM file')
    parser.add_argument('--vcf', required=True, help='VCF file')
    parser.add_argument('--output', '-o', required=True, help='Output directory')
    parser.add_argument('--ground-truth', help='Ground truth CSV (optional)')

    args = parser.parse_args()

    ground_truth = None
    if args.ground_truth:
        ground_truth = pd.read_csv(args.ground_truth)

    results = run_full_phaser_comparison(
        bam_file=args.bam,
        vcf_file=args.vcf,
        output_dir=args.output,
        ground_truth=ground_truth
    )

    print(f"\nphASER Results Summary:")
    print(f"  Total variants: {len(results)}")
    print(f"  Variants counted: {(results['phaser_status'] == 'OK').sum()}")
    print(f"  Missing: {(results['phaser_status'] == 'MISSING').sum()}")

    # By variant type
    for vtype in ['SNP', 'INS', 'DEL']:
        vtype_df = results[results['variant_type'] == vtype]
        ok = (vtype_df['phaser_status'] == 'OK').sum()
        print(f"  {vtype}: {ok}/{len(vtype_df)} counted")


if __name__ == '__main__':
    main()
```

---

## Tasks

### Task 1: Create directory structure
```bash
mkdir -p simulation/competitors
```

### Task 2: Create `simulation/competitors/run_phaser.py`
Implement the script as specified above.

### Task 3: Test on simulation data
```bash
# Use data from Agent A's simulation
SIM_DIR="simulation_results/benchmark_v3"

python simulation/competitors/run_phaser.py \
    --bam ${SIM_DIR}/aligned.sorted.bam \
    --vcf ${SIM_DIR}/variants.vcf.gz \
    --output ${SIM_DIR}/phaser_comparison/ \
    --ground-truth ${SIM_DIR}/ground_truth.csv
```

### Task 4: Handle phASER installation issues
phASER has several dependencies. Handle gracefully:
- intervaltree
- scipy
- pysam

---

## phASER Limitations for INDELs

phASER was primarily designed for **SNPs**. Its INDEL handling:
1. Can count alleles at INDEL positions
2. Uses read-backed phasing (may fail for isolated INDELs)
3. Less validated than SNP counting

**Key Comparison Point:**
> "phASER provides robust SNP-based ASE counting with read-backed phasing but was not specifically optimized for INDEL allele counting. WASP2 provides direct INDEL support."

---

## Success Criteria

- [ ] phASER runs without errors
- [ ] Allelic counts parsed correctly
- [ ] Comparison merged with VCF variants
- [ ] Missing variants tracked
- [ ] Ground truth comparison functional
- [ ] Works with Agent A's simulation output

---

## Commit Template

```bash
git add simulation/competitors/run_phaser.py
git commit -m "feat(sim): add phASER comparison with allelic count parsing

Compares WASP2 against phASER for ASE quantification:
- Prepares phased VCF from simulation ground truth
- Runs phASER allelic counting
- Merges with VCF variants for comparison
- Tracks missing variants and INDEL coverage

Key finding: phASER optimized for SNPs, limited INDEL validation.
"
```
