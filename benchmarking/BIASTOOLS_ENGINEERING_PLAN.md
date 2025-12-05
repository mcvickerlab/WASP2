# Biastools Validation Engineering Plan

## Overview

Use [biastools](https://github.com/maojanlin/biastools) (Genome Biology, 2024) to validate that WASP2 effectively removes reference mapping bias. This provides biological validation complementary to the GATK/phASER concordance benchmarks.

## Tool Information

- **Paper**: [Measuring, visualizing, and diagnosing reference bias with biastools](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03240-8) (Genome Biology, April 2024)
- **GitHub**: https://github.com/maojanlin/biastools
- **PyPI**: https://pypi.org/project/biastools/
- **Installation**: `pip install biastools`

## Validation Strategy

### Approach: Before/After Comparison

Compare reference bias metrics at heterozygous variant sites:

1. **Before WASP2**: Raw STAR-aligned BAM → measure baseline reference bias
2. **After WASP2**: WASP2-filtered BAM → measure reduced reference bias

### Key Metrics to Compare

| Metric | Description | Expected Change |
|--------|-------------|-----------------|
| **Biased sites** | # of variants showing significant ref bias | Decrease |
| **Mean ALT fraction** | Average ALT allele fraction at het sites | Increase toward 0.5 |
| **Bias distribution** | 25th/50th/75th percentile of ALT fractions | Tighter around 0.5 |

## Input Files

### Already Available

```bash
# Reference genome
REF="/iblm/netapp/data1/external/GRC38/combined/GDC_hg38/bwa_index/GRCh38.d1.vd1.fa"

# VCF with heterozygous variants (HG00731)
VCF="benchmarking/comprehensive_results/HG00731_biallelic_snps.vcf.gz"

# After WASP2 filtering (current benchmark BAM)
BAM_AFTER="benchmarking/comprehensive_results/A_sorted_rg.bam"
```

### Need to Locate/Generate

```bash
# Before WASP2 filtering (raw STAR output)
# Option A: Use original STAR BAM before WASP2 processing
# Option B: Re-run STAR alignment without WASP2
```

## Execution Plan

### Phase 1: Installation

```bash
source /iblm/netapp/home/jjaureguy/mambaforge/etc/profile.d/conda.sh
conda activate WASP2_dev2
pip install biastools
```

### Phase 2: Locate Pre-WASP2 BAM

Check for original STAR output before WASP2 filtering:

```bash
# Look for raw alignment BAM
ls -la benchmarking/star_wasp_comparison/results/*/
```

### Phase 3: Run Biastools Analysis

#### Option A: Real Data Analysis (Recommended)

```bash
# Analyze BEFORE WASP2 (raw alignment)
biastools --analyze --real -t 8 \
    -o benchmarking/biastools_results/before_wasp2 \
    -g $REF \
    -v $VCF \
    -s HG00731 \
    -r before_wasp2 \
    --bam $BAM_BEFORE

# Analyze AFTER WASP2 (filtered alignment)
biastools --analyze --real -t 8 \
    -o benchmarking/biastools_results/after_wasp2 \
    -g $REF \
    -v $VCF \
    -s HG00731 \
    -r after_wasp2 \
    --bam $BAM_AFTER
```

#### Option B: Comparative BAM Analysis

```bash
# Direct comparison of two BAMs
biastools_scan --compare_bam \
    -o benchmarking/biastools_results/comparison \
    -g $REF \
    -s HG00731 \
    -r wasp2_comparison \
    -i $BAM_BEFORE \
    -i2 $BAM_AFTER
```

### Phase 4: Generate Figures

Create publication-quality figures showing:

1. **Bias distribution histogram**: ALT fraction distribution before/after
2. **Biased sites count**: Bar chart comparing # of biased sites
3. **Violin/box plot**: REF ratio distribution at het sites

## Expected Output Files

```
benchmarking/biastools_results/
├── before_wasp2/
│   ├── HG00731.bias.all          # Per-variant bias report
│   ├── HG00731.bias.pdf          # Allele balance plots
│   └── HG00731.suspicious.bed    # Suspicious sites
├── after_wasp2/
│   ├── HG00731.bias.all
│   ├── HG00731.bias.pdf
│   └── HG00731.suspicious.bed
└── figures/
    ├── bias_reduction_histogram.png
    ├── bias_reduction_histogram.pdf
    ├── biased_sites_comparison.png
    └── biased_sites_comparison.pdf
```

## Publication Figure Concept

```
┌─────────────────────────────────────────────────────────────┐
│  Reference Bias Reduction by WASP2                          │
├─────────────────────────────────────────────────────────────┤
│                                                             │
│  A. ALT Fraction Distribution       B. Biased Sites Count  │
│  ┌─────────────────────┐           ┌─────────────────────┐ │
│  │    Before WASP2     │           │  ████████████ 1,234 │ │
│  │    (skewed left)    │           │  Before WASP2       │ │
│  │                     │           │                     │ │
│  │    After WASP2      │           │  ████ 312           │ │
│  │    (centered 0.5)   │           │  After WASP2        │ │
│  │                     │           │  (75% reduction)    │ │
│  └─────────────────────┘           └─────────────────────┘ │
│                                                             │
└─────────────────────────────────────────────────────────────┘
```

## Success Criteria

- [ ] ALT fraction distribution shifts toward 0.5 after WASP2
- [ ] Number of biased sites decreases significantly (>50%)
- [ ] Mean bias score decreases
- [ ] Results are statistically significant (p < 0.05)

## Dependencies

- biastools >= 0.1.0
- samtools >= 1.11
- bcftools >= 1.9
- Python >= 3.7

## Notes

1. biastools is designed for DNA-seq but works on RNA-seq at variant sites
2. Focus on heterozygous SNPs where we expect 50:50 allele ratio
3. The "real data" mode (`--real`) uses context-aware assignment without simulation
4. May need to filter VCF to only include heterozygous variants in expressed regions
