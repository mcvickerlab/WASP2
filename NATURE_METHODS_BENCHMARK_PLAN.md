# WASP2 Nature Methods Benchmark Plan

## Executive Summary

This document outlines a publication-quality benchmark strategy for WASP2, based on analysis of existing ASE tool publications and current best practices in the field.

---

## Literature Analysis

### What Published ASE Tools Did

| Tool | Journal | Year | Simulation Approach | Real Data | Key Metrics |
|------|---------|------|---------------------|-----------|-------------|
| **WASP** | Nature Methods | 2015 | Simulated all overlapping reads at het sites in GM12878 | LCL RNA-seq, ChIP-seq | FDR, allelic balance |
| **STAR+WASP** | Genome Biology | 2024 | **None** - used proxy REF-ALT metric | 16 samples (1000G, ENCODE) | REF-ALT difference, speedup |
| **phASER** | Nature Comms | 2016 | Simulated haplotype connections with error rates | GTEx RNA-seq | Phasing accuracy, correlation with eQTL |
| **biastools** | Genome Biology | 2024 | Full simulation with Mason read simulator | NA12878 WGS | ALT fraction, bias categorization |
| **BEERS2** | Briefings Bioinf | 2024 | CAMPAREE diploid expression + full pipeline | Mouse liver | Ground truth alignments, ASE tracking |

### Critical Insight

**STAR+WASP (our direct competitor) did NOT use simulation** - they stated:
> "Since ground truth is unknown in real RNA-seq data, we use a proxy metric REF-ALT"

This is a **weakness we can exploit** by providing proper simulation benchmarks.

---

## Problems with Our Current Approach

### Issue 1: Simulation Design Flaw
```
Current: All replicates share the SAME position
         → Counts are cumulative, not per-replicate

Required: Each test case at UNIQUE position
          OR track reads by replicate ID
```

### Issue 2: GATK Comparison Mismatch
```
Current: GATK outputs per-BASE, not per-VARIANT
         → 20bp deletion = 20 rows of output
         → Alleles don't match VCF

Required: Proper variant-level aggregation
          OR use tools designed for INDEL ASE
```

### Issue 3: Missing Ground Truth Validation
```
Current: No tracking of which reads came from which allele

Required: Embed allele origin in read names/tags
          → Enables perfect ground truth validation
```

---

## Recommended Benchmark Strategy

### Tier 1: Simulation-Based Validation (Required)

#### Option A: Custom Simulation (Recommended - More Control)

Based on original WASP methodology + improvements:

```python
class ASESimulator:
    """
    Generate paired-end reads with known allele origin.

    Key features:
    1. Each variant at UNIQUE genomic position
    2. Allele origin encoded in read name: @read_001_HAP1 or @read_001_HAP2
    3. Known REF:ALT ratio per variant
    4. Realistic insert sizes (300±50bp)
    5. Sequencing errors (0.1-1%)
    6. Support for SNPs, insertions, deletions
    """

    def generate_test_case(
        self,
        chrom: str,
        position: int,
        ref_allele: str,
        alt_allele: str,
        true_ratio: float,  # REF:ALT ratio (1.0 = balanced)
        coverage: int,
        read_length: int = 150,
        insert_size: int = 300,
        insert_std: int = 50,
        error_rate: float = 0.001
    ) -> TestCase:
        """Generate one test case with ground truth."""

        # Calculate expected counts
        total_pairs = coverage
        alt_fraction = 1.0 / (1.0 + true_ratio)
        n_alt = int(total_pairs * alt_fraction)
        n_ref = total_pairs - n_alt

        reads = []
        for i in range(n_ref):
            # Generate read pair from REF haplotype
            r1, r2 = self.make_paired_reads(
                haplotype='REF',
                position=position,
                allele=ref_allele,
                read_id=f"read_{position}_{i}_HAP1"  # Encode origin!
            )
            reads.extend([r1, r2])

        for i in range(n_alt):
            # Generate read pair from ALT haplotype
            r1, r2 = self.make_paired_reads(
                haplotype='ALT',
                position=position,
                allele=alt_allele,
                read_id=f"read_{position}_{i}_HAP2"  # Encode origin!
            )
            reads.extend([r1, r2])

        return TestCase(
            chrom=chrom,
            position=position,
            ref=ref_allele,
            alt=alt_allele,
            true_ratio=true_ratio,
            expected_ref=n_ref,
            expected_alt=n_alt,
            reads=reads
        )
```

#### Option B: Use BEERS2 + CAMPAREE (Industry Standard)

```bash
# Install BEERS2 and CAMPAREE
pip install beers2
git clone https://github.com/itmat/CAMPAREE

# Generate diploid expression with known ASE
camparee \
    --genome hg38.fa \
    --vcf sample_hets.vcf \
    --gtf gencode.v38.gtf \
    --ase-ratios ase_ground_truth.tsv \
    --output camparee_output/

# Simulate reads with BEERS2
beers2 \
    --molecules camparee_output/molecules.tsv \
    --config beers2_config.yaml \
    --output beers2_reads/

# Ground truth is in read names: transcript_1 or transcript_2
```

### Tier 2: Real Data Validation (Required)

#### Dataset 1: GM12878 (Gold Standard)
- Fully phased genome (Genome in a Bottle)
- Extensive RNA-seq data available
- Used by original WASP paper

```bash
# Download GM12878 data
wget ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/...

# Download phased VCF
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/...
```

#### Dataset 2: GTEx Samples
- Use samples with known eQTLs
- Validate ASE correlates with eQTL effect sizes
- Compare to phASER haplotype-level data

#### Dataset 3: HG00731 (Current Benchmark)
- Already have this data
- 56M paired-end reads
- 2.2M heterozygous variants

### Tier 3: Tool Comparison (Required for Nature Methods)

#### Competitors to Benchmark Against:

| Tool | Type | GitHub/Source | Key Metric |
|------|------|---------------|------------|
| **WASP (original)** | Python | github.com/bmvdgeijn/WASP | Baseline |
| **STAR+WASP** | C++ | Built into STAR | Speed competitor |
| **phASER** | Python | github.com/secastel/phaser | Haplotype-level |
| **GATK ASEReadCounter** | Java | GATK toolkit | SNP counting |
| **biastools** | Python | github.com/maojanlin/biastools | Bias measurement |

#### Comparison Metrics:

```python
METRICS = {
    # Accuracy metrics
    'pearson_r': 'Correlation with ground truth ratio',
    'spearman_rho': 'Rank correlation with ground truth',
    'rmse': 'Root mean squared error',
    'mae': 'Mean absolute error',

    # Bias metrics
    'ref_bias': 'Mean(observed_ratio - 0.5) for balanced sites',
    'alt_fraction_deviation': 'Deviation from expected 50%',

    # Performance metrics
    'runtime_seconds': 'Wall clock time',
    'peak_memory_gb': 'Maximum RAM usage',
    'throughput_reads_per_sec': 'Processing speed',

    # Coverage metrics
    'variants_detected': 'Number of variants with counts',
    'variants_missed': 'Variants with zero counts',
    'indel_accuracy': 'Accuracy specifically for INDELs'
}
```

---

## Detailed Test Matrix

### Variant Types (Critical for INDEL Paper)

| Category | Variants | Purpose |
|----------|----------|---------|
| SNPs | A→G, T→C, G→A, C→T | Baseline |
| Small INS | 1-5bp | Standard INDEL |
| Medium INS | 6-15bp | Challenging |
| Large INS | 16-50bp | Edge case |
| Small DEL | 1-5bp | Standard INDEL |
| Medium DEL | 6-15bp | Challenging |
| Large DEL | 16-50bp | Edge case |
| Complex | MNP, delins | Advanced |

### Allelic Ratios to Test

| Ratio (REF:ALT) | Biological Meaning |
|-----------------|-------------------|
| 1:1 (0.5) | Balanced expression |
| 2:1 (0.67) | Moderate imbalance |
| 3:1 (0.75) | Strong imbalance |
| 4:1 (0.80) | Very strong imbalance |
| 10:1 (0.91) | Near monoallelic |
| 1:0 (1.0) | Monoallelic REF |
| 0:1 (0.0) | Monoallelic ALT |

### Coverage Levels

| Coverage | Use Case |
|----------|----------|
| 10x | Low coverage threshold |
| 20x | Minimum recommended |
| 50x | Standard RNA-seq |
| 100x | High coverage |
| 200x | Deep sequencing |

### Test Matrix Size

```
Variant types: 8 categories × 3 sizes = 24 variants
Ratios: 7 levels
Coverage: 5 levels
Replicates: 10 per condition

Total: 24 × 7 × 5 × 10 = 8,400 test cases
```

---

## Implementation Plan

### Phase 1: Fix Simulation Framework (Week 1)

```python
# New file: simulation/ase_simulator_v3.py

class ASESimulatorV3:
    """
    Publication-quality ASE simulation.

    Improvements over v2:
    1. Unique position per test case
    2. Allele origin in read names
    3. Proper INDEL handling
    4. Realistic error models
    """

    def __init__(self, reference_fasta: str):
        self.reference = pysam.FastaFile(reference_fasta)

    def generate_comprehensive_benchmark(
        self,
        output_dir: Path,
        n_variants_per_type: int = 100,
        coverages: List[int] = [20, 50, 100],
        ratios: List[float] = [1.0, 2.0, 4.0],
        seed: int = 42
    ) -> pd.DataFrame:
        """Generate full benchmark suite."""

        # Use DIFFERENT positions for each variant
        positions = self._select_unique_positions(
            n_variants_per_type * len(VARIANT_TYPES)
        )

        ground_truth = []

        for i, (vtype, vsize) in enumerate(VARIANT_CONFIGS):
            pos = positions[i]

            for cov in coverages:
                for ratio in ratios:
                    # Each combination is ONE test case
                    # with ONE unique position
                    test = self.generate_test_case(
                        position=pos,
                        variant_type=vtype,
                        variant_size=vsize,
                        coverage=cov,
                        true_ratio=ratio
                    )
                    ground_truth.append(test.to_dict())

        return pd.DataFrame(ground_truth)
```

### Phase 2: Implement Proper GATK Comparison (Week 1)

```python
# New file: simulation/gatk_comparison_v2.py

def aggregate_gatk_to_variant_level(
    gatk_table: pd.DataFrame,
    vcf_variants: pd.DataFrame
) -> pd.DataFrame:
    """
    Aggregate GATK per-base output to variant level.

    GATK outputs one row per base position.
    For a 10bp deletion, that's 10 rows.
    We need to aggregate to ONE row per variant.
    """

    results = []

    for _, var in vcf_variants.iterrows():
        pos = var['POS']
        ref = var['REF']
        alt = var['ALT']

        # Get all GATK rows overlapping this variant
        var_start = pos
        var_end = pos + len(ref)

        gatk_rows = gatk_table[
            (gatk_table['position'] >= var_start) &
            (gatk_table['position'] < var_end)
        ]

        if len(gatk_rows) == 0:
            # GATK missed this variant
            results.append({
                'chrom': var['CHROM'],
                'pos': pos,
                'ref': ref,
                'alt': alt,
                'gatk_ref_count': 0,
                'gatk_alt_count': 0,
                'gatk_status': 'MISSING'
            })
        else:
            # Aggregate: use the FIRST position for SNPs
            # For INDELs, this is complex...
            if len(ref) == 1 and len(alt) == 1:
                # SNP - simple
                row = gatk_rows.iloc[0]
                results.append({
                    'chrom': var['CHROM'],
                    'pos': pos,
                    'ref': ref,
                    'alt': alt,
                    'gatk_ref_count': row['refCount'],
                    'gatk_alt_count': row['altCount'],
                    'gatk_status': 'OK'
                })
            else:
                # INDEL - GATK doesn't handle well
                # Take counts from first position only
                row = gatk_rows.iloc[0]
                results.append({
                    'chrom': var['CHROM'],
                    'pos': pos,
                    'ref': ref,
                    'alt': alt,
                    'gatk_ref_count': row['refCount'],
                    'gatk_alt_count': row['altCount'],
                    'gatk_status': 'INDEL_APPROX'
                })

    return pd.DataFrame(results)
```

### Phase 3: Run Comprehensive Benchmark (Week 2)

```bash
#!/bin/bash
# run_nature_methods_benchmark.sh

# 1. Generate simulation data
python simulation/ase_simulator_v3.py \
    --reference hg38.fa \
    --output simulation_results/nature_methods_benchmark/ \
    --n-variants 100 \
    --coverages 20,50,100 \
    --ratios 1.0,2.0,4.0,10.0 \
    --seed 42

# 2. Run WASP2 on simulated data
python -m wasp2.pipeline \
    --bam simulation_results/aligned.bam \
    --vcf simulation_results/variants.vcf.gz \
    --output wasp2_results/

# 3. Run competitors
./run_wasp1.sh  # Original WASP
./run_star_wasp.sh  # STAR+WASP
./run_phaser.sh  # phASER
./run_gatk_ase.sh  # GATK ASEReadCounter

# 4. Compare results
python compare_tools.py \
    --ground-truth simulation_results/ground_truth.csv \
    --wasp2 wasp2_results/counts.tsv \
    --wasp1 wasp1_results/counts.tsv \
    --star-wasp star_wasp_results/counts.tsv \
    --gatk gatk_results/counts.tsv \
    --output comparison_results/
```

### Phase 4: Generate Publication Figures (Week 2)

#### Required Figures for Nature Methods:

1. **Figure 1: Method Overview**
   - Pipeline diagram
   - WASP2 algorithm illustration

2. **Figure 2: Simulation Accuracy**
   - Scatter: Observed vs True ratio (colored by variant type)
   - Separate panels for SNP, INS, DEL
   - R² and RMSE annotations

3. **Figure 3: Tool Comparison**
   - Bar chart: Accuracy metrics (Pearson r, RMSE)
   - Grouped by variant type
   - WASP2 vs WASP1 vs STAR+WASP vs GATK

4. **Figure 4: INDEL Performance**
   - Accuracy vs INDEL size
   - Show where GATK fails
   - Highlight WASP2 advantage

5. **Figure 5: Runtime Performance**
   - Already have this
   - 51x faster than WASP1

6. **Supplementary: Real Data Validation**
   - GM12878 allelic balance
   - Correlation with eQTL effects

---

## Success Criteria

### Minimum for Publication:

- [ ] 1000+ simulated test cases with ground truth
- [ ] Pearson r > 0.95 for SNPs
- [ ] Pearson r > 0.90 for INDELs
- [ ] REF bias < 5% (deviation from 0.5)
- [ ] INDEL accuracy >> GATK (statistically significant)
- [ ] Runtime 10x+ faster than WASP1
- [ ] Real data validation on GM12878

### Ideal for High-Impact Publication:

- [ ] 8000+ simulated test cases
- [ ] Multiple real datasets (GM12878, GTEx, HG00731)
- [ ] Comparison against 4+ tools
- [ ] Statistical significance tests (paired t-test, Wilcoxon)
- [ ] Reproducible analysis pipeline
- [ ] Public benchmark dataset release

---

## Timeline

| Week | Tasks |
|------|-------|
| 1 | Fix simulation framework, implement v3 |
| 1 | Fix GATK comparison, proper aggregation |
| 2 | Run comprehensive benchmark (8400 tests) |
| 2 | Generate publication figures |
| 3 | Real data validation (GM12878) |
| 3 | Write methods section |
| 4 | Review, iterate, finalize |

---

## References

1. van de Geijn B, et al. WASP: allele-specific software for robust molecular quantitative trait locus discovery. Nature Methods. 2015;12:1061-1063.

2. Dobin A, et al. STAR+WASP reduces reference bias in the allele-specific mapping of RNA-seq reads. Genome Biology. 2024.

3. Castel SE, et al. Rare variant phasing and haplotypic expression from RNA sequencing with phASER. Nature Communications. 2016;7:12817.

4. Lin MJ, et al. Measuring, visualizing, and diagnosing reference bias with biastools. Genome Biology. 2024;25:104.

5. Brooks TG, et al. BEERS2: RNA-Seq simulation through high fidelity in silico modeling. Briefings in Bioinformatics. 2024;25(3):bbae164.

---

## Appendix: Key URLs

- WASP: https://github.com/bmvdgeijn/WASP
- phASER: https://github.com/secastel/phaser
- biastools: https://github.com/maojanlin/biastools
- BEERS2: https://github.com/itmat/BEERS2
- CAMPAREE: https://github.com/itmat/CAMPAREE
- GATK: https://gatk.broadinstitute.org/hc/en-us/articles/360037054312-ASEReadCounter
