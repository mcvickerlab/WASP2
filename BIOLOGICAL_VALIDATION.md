# Biological Validation: How to Know Indel Allelic Imbalance is Correct

**The Critical Question**: The code works, but are we detecting **real biology** or introducing **false positives/negatives**?

---

## Understanding the Problem

### What WASP2 Does

WASP2 removes **reference mapping bias** by:
1. Taking a read that overlaps heterozygous variants
2. Creating alternate version(s) with swapped alleles
3. Remapping both versions
4. **Keeping only reads that map to the same locus regardless of allele**

**Why this matters for indels**: Aligners have **stronger bias** for indels than SNPs because:
- Indels change read length → alignment score differences
- Indels in repetitive regions → micro-homology shifts
- Large indels → more likely to fail mapping

---

## How to Validate Biological Correctness

### Strategy 1: Positive Controls (Known True Cases)

Use variants with **known allelic imbalance** to verify WASP2 detects them correctly.

#### 1A. Simulated Data with Ground Truth

**Create reads with known allele ratios**:

```bash
# Create test script
cat > validate_simulation.py <<'EOF'
#!/usr/bin/env python3
"""
Simulate reads with known allelic imbalance to validate WASP2.

Ground truth: chr1:1000000 has 2bp insertion
- Haplotype 1 (REF): C
- Haplotype 2 (ALT): CAT
- Expression ratio: 70% REF, 30% ALT (true imbalance)
"""

import pysam
import random
import sys

def create_simulation():
    # Create reference genome segment
    ref_seq = "ATCG" * 1000  # 4000bp reference

    # Insert heterozygous indel at position 1000
    het_pos = 1000
    ref_allele = "C"
    alt_allele = "CAT"  # 2bp insertion

    # Create reads with 70:30 ratio (true allelic imbalance)
    ref_reads = []
    alt_reads = []

    for i in range(100):
        # 70% of reads have reference allele
        if random.random() < 0.7:
            read_seq = ref_seq[het_pos-50:het_pos] + ref_allele + ref_seq[het_pos+1:het_pos+101]
            ref_reads.append(read_seq)
        else:
            # 30% have alternate allele (insertion)
            read_seq = ref_seq[het_pos-50:het_pos] + alt_allele + ref_seq[het_pos+1:het_pos+99]
            alt_reads.append(read_seq)

    print(f"Generated simulation:")
    print(f"  Reference allele reads: {len(ref_reads)} (expected: ~70)")
    print(f"  Alternate allele reads: {len(alt_reads)} (expected: ~30)")
    print(f"  True ratio: {len(ref_reads)}/{len(alt_reads)} = {len(ref_reads)/len(alt_reads):.2f}")

    return ref_reads, alt_reads

if __name__ == "__main__":
    create_simulation()
EOF

python validate_simulation.py
```

**Expected validation**:
1. Before WASP: Raw counts may show bias (e.g., 80:20 due to mapping bias)
2. After WASP: Corrected counts should match ground truth (70:30)

---

#### 1B. Known Imprinted Genes

Use **imprinted genes** (expected ~100:0 or 0:100 ratio):

**Example genes**:
- **H19** (chr11:2016405-2019209) - maternally expressed
- **IGF2** (chr11:2150342-2170832) - paternally expressed
- **SNRPN** (chr15:25068793-25232196) - paternally expressed

**Validation**:
```bash
# Extract reads from imprinted region
python -m mapping make-reads \
    sample.bam \
    variants.vcf.gz \
    --samples SAMPLE_ID \
    --indels \
    --out_dir imprinted_test/

# Remap and filter
# ... (remap step)

# Check allelic ratio at heterozygous indels in imprinted region
# Should be ~100:0 or 0:100 for true imprinting
```

**Expected**: WASP-filtered reads should show **extreme bias** (>90:10) for imprinted genes.

---

### Strategy 2: Negative Controls (Should NOT Show Imbalance)

Use regions with **expected balanced expression** (50:50 ratio).

#### 2A. Housekeeping Genes

Test on constitutively expressed genes (GAPDH, ACTB):

```bash
# Extract housekeeping gene region
samtools view sample.bam chr12:6534405-6538375  # GAPDH

# Run WASP2
python -m mapping make-reads sample.bam variants.vcf.gz \
    --samples SAMPLE_ID \
    --indels \
    --out_dir housekeeping_test/

# Count alleles at heterozygous indels
# Should be ~50:50 (no imbalance)
```

**Expected**: After WASP filtering, allelic ratio should be **close to 1:1** (within statistical noise).

---

#### 2B. Technical Replicates

If you have **technical replicates** (same sample, sequenced twice):

```bash
# Run WASP on both replicates
python -m mapping make-reads replicate1.bam variants.vcf.gz --indels --out_dir rep1/
python -m mapping make-reads replicate2.bam variants.vcf.gz --indels --out_dir rep2/

# Compare allelic ratios between replicates
# Should be highly correlated (R² > 0.95)
```

**Expected**: If WASP is correct, **same biology should give same allelic ratios** across replicates.

---

### Strategy 3: Compare SNPs vs Indels

**Hypothesis**: Indels and nearby SNPs should show **consistent allelic imbalance** (they're on the same haplotype).

#### Test Haplotype Consistency

```python
# Check if indel and nearby SNP show same allelic ratio
# Example: chr1:1000000 (indel) and chr1:1000050 (SNP)
# Both heterozygous, on same haplotype
# Should have same allelic ratio after WASP filtering

import pysam

def check_haplotype_consistency(bam_file, indel_pos, snp_pos):
    """Verify indel and SNP on same haplotype show same imbalance."""

    bam = pysam.AlignmentFile(bam_file)

    # Count alleles at indel position
    indel_ref = 0
    indel_alt = 0

    # Count alleles at SNP position
    snp_ref = 0
    snp_alt = 0

    for read in bam.fetch(...):
        # ... count alleles
        pass

    indel_ratio = indel_ref / indel_alt
    snp_ratio = snp_ref / snp_alt

    print(f"Indel ratio: {indel_ratio:.2f}")
    print(f"SNP ratio:   {snp_ratio:.2f}")
    print(f"Difference:  {abs(indel_ratio - snp_ratio):.2f}")

    # Should be very similar (< 0.2 difference)
    assert abs(indel_ratio - snp_ratio) < 0.2, "Haplotype inconsistency!"
```

**Expected**: Variants on the **same haplotype** should show **correlated allelic ratios**.

---

### Strategy 4: Before/After WASP Comparison

**Key test**: WASP should **reduce bias**, not create it.

#### 4A. Measure Mapping Bias

```python
#!/usr/bin/env python3
"""
Measure reference allele bias before and after WASP filtering.

Expected: WASP should move ratio closer to 50:50 for balanced genes.
"""

import pysam

def measure_bias(bam_file, vcf_file, region):
    """Calculate reference allele frequency before/after WASP."""

    # Before WASP (raw BAM)
    raw_ref_count = 0
    raw_alt_count = 0

    bam = pysam.AlignmentFile(bam_file)
    for read in bam.fetch(region):
        # Count reference vs alternate alleles at heterozygous indels
        # ... (parse CIGAR for indels)
        pass

    raw_ratio = raw_ref_count / (raw_ref_count + raw_alt_count)

    # After WASP (filtered BAM)
    wasp_ref_count = 0
    wasp_alt_count = 0

    # ... similar counting on WASP-filtered BAM

    wasp_ratio = wasp_ref_count / (wasp_ref_count + wasp_alt_count)

    print(f"Before WASP: {raw_ratio:.2%} reference allele")
    print(f"After WASP:  {wasp_ratio:.2%} reference allele")
    print(f"Bias reduction: {abs(0.5 - wasp_ratio) / abs(0.5 - raw_ratio):.1%}")

    # For balanced genes, WASP should reduce bias
    # (move closer to 50%)
    return raw_ratio, wasp_ratio

# Test on housekeeping gene
measure_bias("sample.bam", "variants.vcf.gz", "chr12:6534405-6538375")
```

**Expected for balanced genes**:
- Before WASP: ~60% reference allele (reference bias)
- After WASP: ~50% reference allele (bias corrected)

**Expected for imprinted genes**:
- Before WASP: ~80% one allele (true + bias)
- After WASP: ~100% one allele (bias removed, true imbalance revealed)

---

### Strategy 5: Statistical Tests

#### 5A. Binomial Test for Allelic Imbalance

```python
from scipy.stats import binom_test

def test_allelic_imbalance(ref_count, alt_count, alpha=0.05):
    """Test if allelic ratio differs significantly from 50:50."""

    total = ref_count + alt_count
    p_value = binom_test(ref_count, total, p=0.5, alternative='two-sided')

    if p_value < alpha:
        print(f"Significant imbalance detected (p={p_value:.4f})")
        print(f"  Ratio: {ref_count}:{alt_count} ({ref_count/total:.1%})")
        return True
    else:
        print(f"No significant imbalance (p={p_value:.4f})")
        return False

# Example: Test indel at chr1:1000000
test_allelic_imbalance(ref_count=75, alt_count=25)  # Significant
test_allelic_imbalance(ref_count=52, alt_count=48)  # Not significant
```

**Use this to**:
- Identify which indels have **statistically significant** allelic imbalance
- Set thresholds for calling imbalanced indels

---

#### 5B. Compare to Known eQTL Databases

Cross-reference your indels with **known eQTLs** (expression QTLs):

**Databases**:
- **GTEx** (Genotype-Tissue Expression): https://gtexportal.org
- **eQTL Catalogue**: https://www.ebi.ac.uk/eqtl/

**Validation**:
```bash
# Download GTEx eQTLs for your tissue type
# Check if your indels with allelic imbalance are known eQTLs

# Expected: Indels showing imbalance should be enriched for eQTLs
```

---

### Strategy 6: Orthogonal Validation Methods

#### 6A. RNA-seq Allelic Counts (ASEReadCounter)

Use **GATK ASEReadCounter** as independent method:

```bash
# GATK method (different from WASP)
gatk ASEReadCounter \
    -R genome.fa \
    -I sample.bam \
    -V variants.vcf.gz \
    -O ase_counts.txt

# Compare WASP indel ratios vs GATK ratios
# Should be highly correlated (R² > 0.9)
```

**Expected**: Different methods should **agree on allelic ratios**.

---

#### 6B. Allele-Specific PCR Validation

For critical indels, validate with **wet-lab experiments**:

1. Design primers flanking the indel
2. PCR amplify from cDNA
3. Sanger sequence to measure allelic ratio
4. Compare to WASP-derived ratio

**Gold standard**: Lab validation confirms computational results.

---

## Comprehensive Validation Pipeline

### Step 1: Run WASP2 on Test Dataset

```bash
# Use dataset with known biology (e.g., lymphoblastoid cell lines with known ASE)
python -m mapping make-reads \
    GM12878.bam \
    1000genomes_indels.vcf.gz \
    --samples GM12878 \
    --indels \
    --max-indel-len 10 \
    --out_dir wasp_output/

# Remap
bwa mem genome.fa wasp_output/swapped_alleles_*.fq | samtools sort -o wasp_output/remapped.bam -

# Filter
python -m mapping filter-remapped wasp_output/remapped.bam \
    --json wasp_output/*_wasp_data_files.json \
    --same-locus-slop 2
```

---

### Step 2: Count Alleles at Indels

```python
#!/usr/bin/env python3
"""Count reference vs alternate alleles at heterozygous indels."""

import pysam
import polars as pl

def count_alleles_at_indels(bam_file, vcf_file, output_file):
    """Count alleles at all heterozygous indels."""

    results = []

    vcf = pysam.VariantFile(vcf_file)
    bam = pysam.AlignmentFile(bam_file)

    for variant in vcf:
        if not is_indel(variant):
            continue

        if not is_het(variant):
            continue

        # Count alleles at this position
        ref_count = 0
        alt_count = 0

        for pileupcolumn in bam.pileup(variant.chrom, variant.pos, variant.pos+1):
            if pileupcolumn.pos == variant.pos:
                for pileupread in pileupcolumn.pileups:
                    # Determine which allele this read supports
                    allele = get_allele_from_read(pileupread, variant)
                    if allele == "REF":
                        ref_count += 1
                    elif allele == "ALT":
                        alt_count += 1

        results.append({
            "chrom": variant.chrom,
            "pos": variant.pos,
            "ref": variant.ref,
            "alt": variant.alts[0],
            "ref_count": ref_count,
            "alt_count": alt_count,
            "ratio": ref_count / alt_count if alt_count > 0 else float('inf'),
            "p_value": binom_test(ref_count, ref_count + alt_count, p=0.5)
        })

    df = pl.DataFrame(results)
    df.write_csv(output_file)
    return df

# Run
df = count_alleles_at_indels(
    "wasp_output/remapped.keep.bam",
    "1000genomes_indels.vcf.gz",
    "indel_allelic_counts.csv"
)

# Summary
print(f"Total indels tested: {len(df)}")
print(f"Significant imbalance: {len(df.filter(pl.col('p_value') < 0.05))}")
print(f"Mean ratio: {df['ratio'].mean():.2f}")
```

---

### Step 3: Validate Against Known Truth

```python
# Compare to positive controls (imprinted genes)
imprinted_indels = df.filter(pl.col("gene").is_in(["H19", "IGF2", "SNRPN"]))
print(f"Imprinted indels with extreme ratio (>3:1): {len(imprinted_indels.filter(pl.col('ratio') > 3))}")

# Compare to negative controls (housekeeping genes)
housekeeping_indels = df.filter(pl.col("gene").is_in(["GAPDH", "ACTB"]))
print(f"Housekeeping indels with balanced ratio (0.67-1.5): {len(housekeeping_indels.filter((pl.col('ratio') > 0.67) & (pl.col('ratio') < 1.5)))}")
```

---

### Step 4: Cross-Validate with SNPs

```python
# Load SNP allelic counts from same sample
snp_df = pl.read_csv("snp_allelic_counts.csv")

# For each indel, find nearby SNPs (same haplotype)
for indel in df.iter_rows(named=True):
    nearby_snps = snp_df.filter(
        (pl.col("chrom") == indel["chrom"]) &
        (pl.col("pos").is_between(indel["pos"] - 1000, indel["pos"] + 1000))
    )

    if len(nearby_snps) > 0:
        # Check correlation
        indel_ratio = indel["ratio"]
        snp_ratios = nearby_snps["ratio"].to_list()
        correlation = calculate_correlation(indel_ratio, snp_ratios)

        print(f"Indel {indel['chrom']}:{indel['pos']} - SNP correlation: {correlation:.2f}")
```

---

## Red Flags: When to Suspect Incorrect Results

### 🚩 Red Flag 1: Housekeeping Genes Show Extreme Imbalance

If GAPDH/ACTB have >2:1 ratios → **Something is wrong**

**Possible causes**:
- Mapping bias not fully corrected
- Sample contamination
- Copy number variation

---

### 🚩 Red Flag 2: Indels and SNPs on Same Haplotype Disagree

If indel shows 3:1 but nearby SNP shows 1:1 → **Haplotype inconsistency**

**Possible causes**:
- Indel position mapping error
- Recombination (rare in cis)
- WASP filtering too aggressive

---

### 🚩 Red Flag 3: Technical Replicates Don't Correlate

If Rep1 and Rep2 show different ratios (R² < 0.8) → **Not reproducible**

**Possible causes**:
- Low read depth (stochastic noise)
- Batch effects
- Algorithm instability

---

### 🚩 Red Flag 4: All Indels Show Reference Bias

If 90% of indels have >60% reference allele → **Systematic bias**

**Possible causes**:
- WASP not correcting indel bias properly
- Aligner strongly prefers reference
- Need different `--same-locus-slop` value

---

## Success Criteria

Your indel allelic imbalance calls are **correct** if:

1. ✅ **Positive controls** (imprinted genes) show expected extreme ratios
2. ✅ **Negative controls** (housekeeping genes) show balanced ratios (~1:1)
3. ✅ **Haplotype consistency**: Indels and SNPs agree
4. ✅ **Reproducible**: Technical replicates correlate (R² > 0.9)
5. ✅ **Bias reduction**: WASP moves ratios closer to truth
6. ✅ **Statistical significance**: P-values make sense
7. ✅ **Orthogonal validation**: Agrees with GATK/wet-lab

---

## Recommended Validation Workflow

```bash
# 1. Run on known dataset
python -m mapping make-reads GM12878.bam known_variants.vcf.gz --indels --out_dir test/

# 2. Count alleles
python count_alleles.py test/remapped.keep.bam known_variants.vcf.gz -o counts.csv

# 3. Check controls
python validate_controls.py counts.csv --imprinted H19,IGF2 --housekeeping GAPDH,ACTB

# 4. Compare to SNPs
python compare_indels_vs_snps.py counts.csv snp_counts.csv

# 5. Statistical tests
python test_significance.py counts.csv --alpha 0.05

# 6. Generate report
python generate_validation_report.py counts.csv -o validation_report.html
```

---

## Bottom Line

**You know indel allelic imbalance is correct when**:

1. Known imprinted genes show expected imbalance ✅
2. Housekeeping genes show balanced expression ✅
3. Results are reproducible across replicates ✅
4. Indels and SNPs on same haplotype agree ✅
5. WASP reduces mapping bias (validated with simulations) ✅
6. Orthogonal methods agree (GATK, lab validation) ✅

**The code being correct ≠ biology being correct.** You need **experimental validation** with positive/negative controls!

---

**Next step**: I can help you create the validation scripts above. Which validation strategy would you like to implement first?
