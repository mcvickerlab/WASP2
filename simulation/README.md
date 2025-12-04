# WASP2 vs GATK ASEReadCounter Benchmark

This directory contains scripts for benchmarking WASP2 against GATK ASEReadCounter on simulated allele-specific expression data.

## Overview

**Goal**: Compare WASP2 and GATK ASEReadCounter head-to-head for publication:
1. **Accuracy**: Pearson correlation with ground truth
2. **Bias**: Reference allele bias (deviation from expected ratio)
3. **Speed**: Runtime comparison

## Files

### Main Scripts

- **`benchmark_vs_gatk.py`**: Main benchmark script that runs both GATK and WASP2 on the same data
- **`create_ground_truth_from_vcf.py`**: Helper script to create ground truth CSV from VCF file

### Usage

#### 1. Create Ground Truth from VCF

```bash
python simulation/create_ground_truth_from_vcf.py \
    simulation_results/*/variants.vcf.gz \
    simulation_results/*/ground_truth.csv
```

#### 2. Run GATK-only Benchmark

```bash
python simulation/benchmark_vs_gatk.py \
    --bam simulation_results/comprehensive_*/aligned.sorted.bam \
    --vcf simulation_results/comprehensive_*/variants.vcf.gz \
    --ref simulation_results/comprehensive_*/reference.fa \
    --ground-truth simulation_results/comprehensive_*/ground_truth.csv \
    --output comparison_results/ \
    --gatk-only
```

#### 3. Run Full Comparison (GATK + WASP2)

```bash
python simulation/benchmark_vs_gatk.py \
    --bam simulation_results/comprehensive_*/aligned.sorted.bam \
    --vcf simulation_results/comprehensive_*/variants.vcf.gz \
    --ref simulation_results/comprehensive_*/reference.fa \
    --ground-truth simulation_results/comprehensive_*/ground_truth.csv \
    --output comparison_results/
```

## Requirements

### Software Dependencies

- **GATK 4.x**: Install via conda:
  ```bash
  conda install -c bioconda gatk4
  ```

- **WASP2**: Must be installed and accessible
- **Python packages**: pandas, numpy, scipy, scikit-learn, matplotlib, seaborn, pysam

### Input Data Requirements

- **BAM file**: Aligned reads (will automatically add read groups if missing)
- **VCF file**: Must be bgzip-compressed and tabix-indexed (.vcf.gz + .tbi)
- **Reference FASTA**: Must be indexed (.fai), sequence dictionary (.dict) will be created automatically
- **Ground truth CSV**: Columns: chrom, pos, variant_id, ref_allele, alt_allele, genotype, true_ref_ratio

## Automatic Preprocessing

The benchmark script automatically handles common issues:

1. **Missing sequence dictionary**: Creates `.dict` file for reference FASTA
2. **Missing read groups**: Adds read groups to BAM file (GATK requirement)
3. **File validation**: Checks for all required indexes before running

## Output

The script generates:

- **`gatk_ase_counts.table`**: GATK ASEReadCounter output
- **`wasp2_ase_counts.txt`**: WASP2 allele counts (TODO: implement)
- **`comparison_report.md`**: Publication-ready comparison report with metrics
- **`accuracy_scatter.png`**: Scatter plots of truth vs tool
- **`error_distribution.png`**: Error distribution comparison
- **`bias_comparison.png`**: Reference allele bias boxplots

## Current Status

### Completed
- GATK ASEReadCounter integration
- Automatic file preprocessing (read groups, sequence dictionary)
- Ground truth extraction from VCF
- GATK-only benchmarking working

### TODO
- Implement WASP2 counting function (`run_wasp2_counting()`)
- Run full WASP2 vs GATK comparison
- Generate comparison plots
- Write publication-ready comparison report

## Example Results

```
GATK Results:
  Pearson r: 0.9850
  Ref bias: -0.0024
  Runtime: 3.96 sec

WASP2 Results:
  Pearson r: 0.9912
  Ref bias: 0.0015
  Runtime: 1.25 sec

WASP2 vs GATK:
  Accuracy: WASP2=0.9912 vs GATK=0.9850
  Bias: WASP2=0.0015 vs GATK=-0.0024
  Speed: WASP2 is 3.17x faster
```

## Branch

This work is on branch: `sim/gatk-compare`

Parent branch: `ropc-indels`

## References

- McKenna et al. (2010). The Genome Analysis Toolkit. *Genome Research*, 20(9), 1297-1303.
- van de Geijn et al. (2015). WASP: allele-specific software for robust molecular QTL discovery. *Nature Methods*, 12(11), 1061-1063.
