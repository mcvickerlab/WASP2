# Agent B: GATK ASEReadCounter Comparison

## Mission
Create a proper GATK ASEReadCounter comparison that handles the per-base vs per-variant mismatch and validates WASP2's INDEL advantage.

---

## Repository Context

**GitHub:** https://github.com/Jaureguy760/WASP2-exp.git
**Branch:** `sim/benchmark-v3`
**Working Directory:** `/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp`
**Conda Environment:** `WASP2_dev2`

---

## Problem Statement

### GATK Behavior
GATK ASEReadCounter outputs **per-base**, not **per-variant**:
- SNP at pos 50000 → 1 row
- 10bp deletion at pos 60000 → 10 rows (one per deleted base)
- Alleles in output don't match VCF alleles

### Our Previous Error
```
VCF: chr1  90000  DEL  GCGCGCGCGCGCGCGCGCGC  G
GATK output: 20 rows with positions 90000-90019, random alleles
```

This caused R² = nan in our figures.

---

## GATK ASEReadCounter API

### Command
```bash
gatk ASEReadCounter \
    -R reference.fa \
    -I input.bam \
    -V variants.vcf.gz \
    -O output.table \
    --min-mapping-quality 10 \
    --min-base-quality 20 \
    --min-depth-of-non-filtered-base 10
```

### Output Format
```
contig  position  variantID  refAllele  altAllele  refCount  altCount  totalCount  lowMAPQDepth  lowBaseQDepth  rawDepth  otherBases  improperPairs
chr1    50000     SNP_50000  A          G          25        25        50          0             0              52        2           0
```

### Key Limitations
1. Only processes **biallelic SNPs** well
2. INDELs: Outputs each position separately
3. **Does NOT count INDEL alleles correctly** - counts bases, not INDEL events

---

## Implementation

### File: `simulation/competitors/run_gatk.py`

```python
#!/usr/bin/env python3
"""
GATK ASEReadCounter wrapper with proper variant-level aggregation.

GATK outputs per-base counts. For INDELs, this means multiple rows.
We aggregate back to variant-level for fair comparison.
"""

import subprocess
import pandas as pd
import pysam
from pathlib import Path
from typing import Dict, Optional
import tempfile


def check_gatk_available() -> bool:
    """Check if GATK is installed."""
    try:
        result = subprocess.run(['gatk', '--version'], capture_output=True, text=True)
        return result.returncode == 0
    except FileNotFoundError:
        return False


def install_gatk():
    """Install GATK via conda."""
    print("Installing GATK...")
    subprocess.run(['conda', 'install', '-y', '-c', 'bioconda', 'gatk4'], check=True)


def prepare_bam_for_gatk(input_bam: str, output_bam: str):
    """
    Prepare BAM for GATK (add read groups if missing).

    GATK requires read groups. Simulation BAMs may not have them.
    """
    # Check if BAM has read groups
    bam = pysam.AlignmentFile(input_bam, 'rb')
    has_rg = 'RG' in bam.header

    if has_rg:
        # Just copy
        subprocess.run(['cp', input_bam, output_bam], check=True)
        subprocess.run(['samtools', 'index', output_bam], check=True)
    else:
        # Add read groups
        subprocess.run([
            'gatk', 'AddOrReplaceReadGroups',
            '-I', input_bam,
            '-O', output_bam,
            '-RGID', 'SIMULATED',
            '-RGLB', 'lib1',
            '-RGPL', 'ILLUMINA',
            '-RGPU', 'unit1',
            '-RGSM', 'SIMULATED'
        ], check=True)
        subprocess.run(['samtools', 'index', output_bam], check=True)


def prepare_reference_for_gatk(ref_fasta: str):
    """
    Prepare reference for GATK (create sequence dictionary if missing).
    """
    dict_file = ref_fasta.replace('.fa', '.dict')
    if not Path(dict_file).exists():
        subprocess.run([
            'gatk', 'CreateSequenceDictionary',
            '-R', ref_fasta,
            '-O', dict_file
        ], check=True)


def run_gatk_ase_counter(
    bam_file: str,
    vcf_file: str,
    ref_fasta: str,
    output_table: str,
    min_mapping_quality: int = 10,
    min_base_quality: int = 20,
    min_depth: int = 1
) -> str:
    """
    Run GATK ASEReadCounter.

    Returns path to output table.
    """
    # Prepare inputs
    with tempfile.NamedTemporaryFile(suffix='.bam', delete=False) as tmp:
        prepared_bam = tmp.name

    prepare_bam_for_gatk(bam_file, prepared_bam)
    prepare_reference_for_gatk(ref_fasta)

    # Run GATK
    cmd = [
        'gatk', 'ASEReadCounter',
        '-R', ref_fasta,
        '-I', prepared_bam,
        '-V', vcf_file,
        '-O', output_table,
        '--min-mapping-quality', str(min_mapping_quality),
        '--min-base-quality', str(min_base_quality),
        '--min-depth-of-non-filtered-base', str(min_depth)
    ]

    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"GATK error: {result.stderr}")
        raise RuntimeError(f"GATK failed: {result.stderr}")

    # Cleanup
    Path(prepared_bam).unlink(missing_ok=True)
    Path(prepared_bam + '.bai').unlink(missing_ok=True)

    return output_table


def load_vcf_variants(vcf_file: str) -> pd.DataFrame:
    """Load variants from VCF."""
    variants = []
    vcf = pysam.VariantFile(vcf_file)

    for rec in vcf:
        variants.append({
            'chrom': rec.chrom,
            'pos': rec.pos,
            'ref': rec.ref,
            'alt': rec.alts[0] if rec.alts else '',
            'variant_id': rec.id or f"{rec.chrom}_{rec.pos}",
            'variant_type': classify_variant(rec.ref, rec.alts[0] if rec.alts else '')
        })

    return pd.DataFrame(variants)


def classify_variant(ref: str, alt: str) -> str:
    """Classify variant type."""
    if len(ref) == 1 and len(alt) == 1:
        return 'SNP'
    elif len(ref) < len(alt):
        return 'INS'
    else:
        return 'DEL'


def aggregate_gatk_to_variant_level(
    gatk_table: str,
    vcf_file: str
) -> pd.DataFrame:
    """
    Aggregate GATK per-base output to variant level.

    Strategy:
    - SNPs: Use single row directly
    - INDELs: Use FIRST position (GATK can't count INDEL alleles properly)
    - Mark INDEL counts as approximate
    """
    # Load GATK output
    gatk = pd.read_csv(gatk_table, sep='\t', comment='#')

    # Load VCF variants
    vcf_variants = load_vcf_variants(vcf_file)

    results = []

    for _, var in vcf_variants.iterrows():
        pos = var['pos']
        ref = var['ref']
        alt = var['alt']
        vtype = var['variant_type']

        # Find GATK row(s) for this position
        gatk_rows = gatk[gatk['position'] == pos]

        if len(gatk_rows) == 0:
            # GATK missed this variant
            results.append({
                'chrom': var['chrom'],
                'pos': pos,
                'ref': ref,
                'alt': alt,
                'variant_type': vtype,
                'gatk_ref_count': 0,
                'gatk_alt_count': 0,
                'gatk_total': 0,
                'gatk_status': 'MISSING',
                'gatk_ratio': float('nan')
            })
        else:
            row = gatk_rows.iloc[0]

            # Calculate ratio (REF / (REF + ALT))
            total = row['refCount'] + row['altCount']
            ratio = row['refCount'] / total if total > 0 else float('nan')

            # Status depends on variant type
            if vtype == 'SNP':
                status = 'SNP_OK'
            else:
                # GATK doesn't count INDELs correctly
                # The counts are for individual bases, not INDEL events
                status = f'INDEL_APPROXIMATE'

            results.append({
                'chrom': var['chrom'],
                'pos': pos,
                'ref': ref,
                'alt': alt,
                'variant_type': vtype,
                'gatk_ref_count': row['refCount'],
                'gatk_alt_count': row['altCount'],
                'gatk_total': total,
                'gatk_status': status,
                'gatk_ratio': ratio
            })

    return pd.DataFrame(results)


def run_full_gatk_comparison(
    bam_file: str,
    vcf_file: str,
    ref_fasta: str,
    output_dir: str
) -> pd.DataFrame:
    """
    Run full GATK comparison pipeline.

    Returns DataFrame with aggregated variant-level counts.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Check GATK
    if not check_gatk_available():
        install_gatk()

    # Run GATK
    gatk_table = str(output_dir / 'gatk_raw.table')
    run_gatk_ase_counter(bam_file, vcf_file, ref_fasta, gatk_table)

    # Aggregate to variant level
    gatk_counts = aggregate_gatk_to_variant_level(gatk_table, vcf_file)

    # Save
    gatk_counts.to_csv(output_dir / 'gatk_variant_level.csv', index=False)

    return gatk_counts


def main():
    import argparse

    parser = argparse.ArgumentParser(description='GATK ASEReadCounter comparison')
    parser.add_argument('--bam', required=True, help='Input BAM file')
    parser.add_argument('--vcf', required=True, help='VCF file')
    parser.add_argument('--ref', required=True, help='Reference FASTA')
    parser.add_argument('--output', '-o', required=True, help='Output directory')

    args = parser.parse_args()

    results = run_full_gatk_comparison(
        bam_file=args.bam,
        vcf_file=args.vcf,
        ref_fasta=args.ref,
        output_dir=args.output
    )

    print(f"\nGATK Results Summary:")
    print(f"  Total variants: {len(results)}")
    print(f"  SNPs OK: {(results['gatk_status'] == 'SNP_OK').sum()}")
    print(f"  INDELs (approximate): {results['gatk_status'].str.contains('INDEL').sum()}")
    print(f"  Missing: {(results['gatk_status'] == 'MISSING').sum()}")


if __name__ == '__main__':
    main()
```

---

## Tasks

### Task 1: Create directory structure
```bash
mkdir -p simulation/competitors
```

### Task 2: Create `simulation/competitors/run_gatk.py`
Implement the script as specified above.

### Task 3: Test on existing simulation data
```bash
# Use data from previous simulation runs
SIM_DIR=$(ls -td simulation_results/paired_end_* 2>/dev/null | head -1)

python simulation/competitors/run_gatk.py \
    --bam ${SIM_DIR}/aligned.sorted.bam \
    --vcf ${SIM_DIR}/variants.vcf.gz \
    --ref ${SIM_DIR}/reference.fa \
    --output /tmp/gatk_test/
```

### Task 4: Validate INDEL handling
Check that INDELs are marked as approximate and explain the limitation.

---

## Success Criteria

- [ ] GATK runs without errors
- [ ] SNPs have status 'SNP_OK'
- [ ] INDELs have status 'INDEL_APPROXIMATE'
- [ ] Output CSV has one row per variant (not per base)
- [ ] gatk_ratio column calculated correctly
- [ ] Works with simulation data from Agent A

---

## Key Insight for Paper

**GATK ASEReadCounter cannot properly count INDEL alleles.**

It counts bases at each position, not INDEL events. This is a fundamental limitation that WASP2 addresses with proper INDEL-aware counting.

Quote for paper:
> "GATK ASEReadCounter outputs per-base allele counts, which provides accurate results for SNPs but cannot properly quantify allele-specific expression at INDEL sites. WASP2 addresses this limitation with INDEL-aware allele counting."

---

## Commit Template

```bash
git add simulation/competitors/run_gatk.py
git commit -m "feat(sim): add GATK ASEReadCounter comparison with variant-level aggregation

Addresses GATK's per-base output limitation:
- Aggregates per-base counts to variant level
- Marks INDEL counts as approximate
- Documents GATK's inability to count INDEL alleles

Key finding: GATK cannot properly count INDEL alleles,
providing evidence for WASP2's INDEL support value.
"
```
