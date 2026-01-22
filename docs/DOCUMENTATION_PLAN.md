# WASP2 Professional Documentation Plan

## Executive Summary

This document provides a comprehensive plan for creating professional, user-friendly documentation for WASP2, a bioinformatics CLI tool for allele-specific analysis. The plan draws on best practices from successful tools like STAR, salmon, cellranger, and bcftools.

**Current State**: WASP2 has solid foundation documentation (README, Sphinx API docs, user guides).

**Goal**: Elevate documentation to production-grade standards with progressive tutorials, comprehensive CLI help, improved discoverability, and accessibility for users at all skill levels.

---

## 1. README Best Practices

### Current Strengths
- Clear logo and branding
- Good badge coverage (Docs, PyPI, License, Python/Rust versions)
- Comprehensive CLI quick reference
- Performance documentation (VCF/PGEN formats)
- Installation instructions including conda and Rust build

### Recommended Improvements

#### 1.1 Enhanced Badge Section
```markdown
<p align="center">
  <!-- CI/CD Status -->
  <a href="https://github.com/Jaureguy760/WASP2-exp/actions/workflows/ci.yml">
    <img src="https://github.com/Jaureguy760/WASP2-exp/actions/workflows/ci.yml/badge.svg" alt="CI">
  </a>
  <a href="https://codecov.io/gh/Jaureguy760/WASP2-exp">
    <img src="https://codecov.io/gh/Jaureguy760/WASP2-exp/branch/main/graph/badge.svg" alt="Coverage">
  </a>

  <!-- Package Status -->
  <a href="https://pypi.org/project/wasp2/">
    <img src="https://img.shields.io/pypi/v/wasp2" alt="PyPI">
  </a>
  <a href="https://anaconda.org/bioconda/wasp2">
    <img src="https://img.shields.io/conda/vn/bioconda/wasp2" alt="Bioconda">
  </a>

  <!-- Documentation -->
  <a href="https://jaureguy760.github.io/WASP2-exp/">
    <img src="https://img.shields.io/badge/docs-GitHub%20Pages-blue" alt="Documentation">
  </a>

  <!-- Language & License -->
  <a href="https://github.com/Jaureguy760/WASP2-exp/blob/master/LICENSE">
    <img src="https://img.shields.io/badge/license-MIT-green" alt="License">
  </a>
  <img src="https://img.shields.io/badge/python-3.10+-blue" alt="Python">
  <img src="https://img.shields.io/badge/rust-1.70+-orange" alt="Rust">

  <!-- Usage Stats (Optional) -->
  <img src="https://img.shields.io/pypi/dm/wasp2" alt="Downloads">
  <img src="https://img.shields.io/github/stars/Jaureguy760/WASP2-exp?style=social" alt="Stars">
</p>
```

#### 1.2 Quick Start Section (Front and Center)
Place BEFORE installation for better UX. Users want to see *what* before *how*.

```markdown
## Quick Start (5 minutes)

Get started with WASP2 in three commands:

```bash
# 1. Install
pip install wasp2

# 2. Count allele-specific reads
wasp2-count count-variants sample.bam variants.vcf.gz -s sample1 -o counts.tsv

# 3. Detect allelic imbalance
wasp2-analyze find-imbalance counts.tsv -o results.tsv
```

**Output**: Statistical test results for allelic imbalance at heterozygous SNPs.

**Next**: See [Full Tutorial](#tutorial) or [Documentation](https://jaureguy760.github.io/WASP2-exp/)
```

#### 1.3 Feature Highlights Section
Use visual hierarchy and icons (plain text, not emoji):

```markdown
## Key Features

### Allele-Specific Analysis
- **Count Module**: Quantify ref/alt allele reads at heterozygous SNPs
- **Analysis Module**: Beta-binomial statistical testing for allelic imbalance
- **Mapping Module**: WASP algorithm for unbiased read mapping

### Performance
- **Rust Acceleration**: Core algorithms implemented in Rust (10-25x faster)
- **Multi-Format Support**: VCF, BCF, PGEN (up to 25x faster I/O)
- **High-Performance VCF**: Optional cyvcf2 backend (7x faster parsing)

### Applications
- RNA-seq allele-specific expression (ASE)
- ATAC-seq allelic chromatin accessibility
- Single-cell RNA-seq/ATAC-seq
- ChIP-seq allelic binding analysis

### Data Types
- Bulk RNA-seq, ATAC-seq, ChIP-seq
- Single-cell RNA-seq (10x Genomics, etc.)
- Paired-end and single-end reads
- Any organism with a reference genome
```

#### 1.4 Installation Options Section
Structured by user type:

```markdown
## Installation

### For Users (Recommended)

**Option 1: PyPI (Python package)**
```bash
pip install wasp2
```

**Option 2: Bioconda** (when available)
```bash
conda install -c bioconda wasp2
```

**Option 3: Install with performance enhancements**
```bash
# Install with cyvcf2 (7x faster VCF parsing)
pip install wasp2[cyvcf2]

# Install with PLINK2 support (25x faster variant I/O)
pip install wasp2[plink]

# Install all optional dependencies
pip install wasp2[all]
```

### For Developers

**From source with Rust acceleration:**
```bash
# Clone repository
git clone https://github.com/Jaureguy760/WASP2-exp.git
cd WASP2-exp

# Create conda environment
conda env create -f environment.yml
conda activate WASP2

# Build Rust extension
export LIBCLANG_PATH=$CONDA_PREFIX/lib
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
maturin develop --release -m rust/Cargo.toml

# Install development dependencies
pip install -e ".[dev]"
```

### Cloud Development

**GitHub Codespaces** (zero setup):
1. Click "Code" → "Codespaces" → "Create codespace"
2. Wait 2-3 minutes for automatic setup
3. Start using WASP2 immediately

See [.devcontainer/README.md](.devcontainer/README.md) for details.
```

#### 1.5 Citation Section
Essential for academic tools:

```markdown
## Citation

If you use WASP2 in your research, please cite:

```bibtex
@article{wasp2_2025,
  title={WASP2: High-performance allele-specific analysis of next-generation sequencing data},
  author={Ho, Aaron and Jaureguy, Jeff and McVicker, Graham},
  journal={Bioinformatics},
  year={2025},
  volume={XX},
  pages={XXX-XXX},
  doi={10.1093/bioinformatics/XXXXX}
}
```

**Original WASP paper:**
van de Geijn B, McVicker G, Gilad Y, Pritchard JK (2015). WASP: allele-specific software for robust molecular quantitative trait locus discovery. *Nature Methods* 12:1061-1063. [doi:10.1038/nmeth.3582](https://doi.org/10.1038/nmeth.3582)
```

#### 1.6 Comparison Table
Help users understand positioning:

```markdown
## Comparison with Other Tools

| Feature | WASP2 | GATK ASEReadCounter | phASER | MBASED |
|---------|-------|---------------------|---------|---------|
| **Mapping Bias Correction** | Yes (WASP) | No | No | No |
| **Statistical Testing** | Beta-binomial | No | Phasing only | Beta-binomial |
| **Single-Cell Support** | Yes | No | No | No |
| **Performance** | Fast (Rust) | Slow | Medium | Medium |
| **Variant Formats** | VCF/BCF/PGEN | VCF only | VCF only | VCF only |
| **Indel Support** | Yes | Yes | No | No |
| **License** | MIT | BSD | MIT | GPL |
```

#### 1.7 Learning Path Section
Guide users to appropriate resources:

```markdown
## Learning Resources

- **New to allele-specific analysis?** Start with [Concepts](docs/concepts.md)
- **Want to try WASP2 quickly?** Follow [Quick Start Tutorial](docs/tutorials/quickstart.md) (5 min)
- **Analyzing RNA-seq?** See [RNA-seq ASE Tutorial](docs/tutorials/rnaseq_ase.md) (30 min)
- **Working with ATAC-seq?** See [ATAC-seq Tutorial](docs/tutorials/atac_ase.md) (30 min)
- **Single-cell data?** See [Single-Cell Guide](docs/tutorials/single_cell.md) (45 min)
- **Need API reference?** Browse [API Documentation](https://jaureguy760.github.io/WASP2-exp/)
```

---

## 2. Tutorial Types and Structure

### 2.1 Tutorial Hierarchy

```
tutorials/
├── 00_concepts.md              # Background for newcomers
├── 01_quickstart.md            # 5-minute intro
├── 02_installation_guide.md    # Comprehensive setup
├── 03_basic_workflow.md        # Complete pipeline walkthrough
├── 04_rnaseq_ase.md           # RNA-seq specific
├── 05_atac_ase.md             # ATAC-seq specific
├── 06_single_cell.md          # Single-cell workflows
├── 07_advanced_options.md     # Power user features
├── 08_troubleshooting.md      # Common issues
└── 09_performance_tuning.md   # Optimization guide
```

### 2.2 Tutorial Template Structure

Each tutorial follows consistent structure (inspired by diataxis framework):

```markdown
# Tutorial Title

**Time**: X minutes
**Level**: Beginner/Intermediate/Advanced
**Prerequisites**: List of required knowledge/tools
**Data**: Link to example data

## Learning Objectives

By the end of this tutorial, you will:
- [ ] Objective 1
- [ ] Objective 2
- [ ] Objective 3

## Background

Brief context (2-3 paragraphs)

## Setup

```bash
# Download example data
wget https://example.com/data.tar.gz
tar -xzf data.tar.gz
cd tutorial_data/
```

## Step 1: [Action Verb]

**Goal**: What you'll accomplish in this step

**Command**:
```bash
wasp2-count count-variants sample.bam variants.vcf.gz \
  --samples NA12878 \
  --region genes.gtf \
  --out_file counts.tsv
```

**Explanation**: Line-by-line breakdown of flags

**Expected Output**:
```
Processing 45,283 variants...
Found 12,456 heterozygous SNPs in NA12878
Counted reads at 9,821 SNPs overlapping genes
Output written to counts.tsv
```

**Verification**:
```bash
head counts.tsv
wc -l counts.tsv  # Should be ~9,822 (header + 9,821 SNPs)
```

## Step 2: [Next Action]

[Same structure...]

## Interpreting Results

**Understanding the output**:
- Column A means...
- Column B means...

**Quality checks**:
1. Check total counts
2. Look for coverage distribution
3. Verify expected patterns

## Next Steps

- Try with your own data
- See [Advanced Tutorial] for more options
- Read about [Concept X] for deeper understanding

## Troubleshooting

**Problem**: Error message X
**Solution**: Do Y

**Problem**: Unexpected results
**Solution**: Check Z

## Summary

Quick recap of what was learned

## Further Reading

- Link to related tutorials
- Link to API docs
- Link to relevant papers
```

### 2.3 Specific Tutorial Content

#### Tutorial 0: Concepts (concepts.md)
```markdown
# Understanding Allele-Specific Analysis

## What is Allelic Imbalance?

In diploid organisms, each individual carries two copies (alleles) of most genes.
Normally, both alleles are expressed equally, but sometimes one allele is
preferentially expressed due to:

1. **Cis-regulatory variants**: SNPs affecting transcription factor binding
2. **Imprinting**: Parent-of-origin-specific expression
3. **X-inactivation**: Random inactivation of one X chromosome
4. **Allele-specific methylation**: Epigenetic regulation

## Why Does Reference Bias Matter?

Standard aligners preferentially map reads matching the reference genome:
- Reads with alternate alleles have more mismatches
- More mismatches = lower alignment scores
- Lower scores = more likely to be filtered

This creates artificial allelic imbalance favoring the reference allele.

## The WASP Solution

WASP corrects reference bias by:
1. Identifying reads overlapping variants
2. Swapping alleles in those reads
3. Re-mapping swapped reads
4. Keeping only reads that map to the same location

[Diagram illustrating WASP workflow]

## When to Use Each WASP2 Module

[Decision tree diagram]

**Counting Module**: Already have unbiased BAM? Just need allele counts?
**Mapping Module**: Have standard BAM? Need to correct reference bias first
**Analysis Module**: Have allele counts? Need statistical testing for imbalance?
```

#### Tutorial 1: Quick Start (quickstart.md)
```markdown
# WASP2 Quick Start (5 minutes)

**Level**: Beginner
**Time**: 5 minutes
**Prerequisites**: Python 3.10+

## 1. Install

```bash
pip install wasp2
```

## 2. Download Example Data

```bash
# Small test dataset (chr10, ~50MB)
wget https://github.com/Jaureguy760/WASP2-exp/raw/main/test_data/quickstart_bundle.tar.gz
tar -xzf quickstart_bundle.tar.gz
cd quickstart_data/
```

Contains:
- `sample.bam` - Aligned RNA-seq reads (chromosome 10)
- `variants.vcf.gz` - Heterozygous SNPs
- `genes.gtf` - Gene annotations

## 3. Count Allele-Specific Reads

```bash
wasp2-count count-variants \
  sample.bam \
  variants.vcf.gz \
  --samples NA12878 \
  --region genes.gtf \
  --out_file counts.tsv
```

**Output**: `counts.tsv` with ref/alt counts per SNP per gene

## 4. Detect Allelic Imbalance

```bash
wasp2-analyze find-imbalance \
  counts.tsv \
  --out_file results.tsv
```

**Output**: `results.tsv` with statistical tests for each gene

## 5. Inspect Results

```bash
# View significant genes (FDR < 0.05)
awk 'NR==1 || $8 < 0.05' results.tsv | column -t

# Count significant genes
awk 'NR>1 && $8 < 0.05' results.tsv | wc -l
```

## What's Next?

- **Understand the output**: See [Interpreting Results](interpreting_results.md)
- **Use your data**: See [Full Pipeline Tutorial](basic_workflow.md)
- **ATAC-seq analysis**: See [ATAC-seq Tutorial](atac_ase.md)
- **Single-cell data**: See [Single-Cell Guide](single_cell.md)
```

#### Tutorial 3: Basic Workflow (basic_workflow.md)
```markdown
# Complete WASP2 Pipeline Walkthrough

**Level**: Intermediate
**Time**: 30 minutes
**Prerequisites**: Basic command line, understanding of BAM/VCF formats

## Overview

This tutorial covers the complete WASP2 workflow:
1. Data preparation and QC
2. WASP mapping (bias correction)
3. Allele counting
4. Statistical analysis
5. Result interpretation

## Pipeline Diagram

```
Raw Reads (FASTQ)
    ↓
Standard Alignment (STAR/BWA/bowtie2)
    ↓
WASP Mapping Filter (wasp2-map)
    ├── make-reads: Generate swapped alleles
    ├── remap: Re-align swapped reads
    └── filter-remapped: Keep consistent mappings
    ↓
Unbiased BAM
    ↓
Allele Counting (wasp2-count)
    ↓
Statistical Analysis (wasp2-analyze)
    ↓
Allelic Imbalance Results
```

## Data Requirements

Before starting, ensure you have:
- [ ] Aligned BAM file (sorted, indexed)
- [ ] VCF file with genotypes (bgzipped, indexed)
- [ ] Optional: Gene/peak annotations (GTF/BED)
- [ ] Sample ID present in VCF

## Step 1: Quality Control

[Detailed QC steps...]

## Step 2: WASP Mapping Filter

[Complete mapping workflow...]

## Step 3: Allele Counting

[Counting with different options...]

## Step 4: Statistical Analysis

[Analysis and interpretation...]

## Step 5: Visualization

[Basic plotting in R/Python...]
```

#### Tutorial 4: RNA-seq ASE (rnaseq_ase.md)
```markdown
# RNA-seq Allele-Specific Expression Analysis

**Level**: Intermediate
**Time**: 45 minutes
**Data**: Download from [link]

## Use Case

You have RNA-seq data from a heterozygous individual and want to:
- Identify genes with allelic imbalance
- Detect potential cis-regulatory variants
- Find imprinted genes

## Biological Questions

1. Which genes show preferential expression of one allele?
2. Are there parent-of-origin effects (imprinting)?
3. Do allelic ratios differ between conditions/tissues?

## Dataset

- Sample: GM12878 (lymphoblastoid cell line)
- Sequencing: Paired-end 100bp RNA-seq
- Depth: ~30M reads
- Genome: GRCh38

## Workflow

### Part A: Standard RNA-seq Alignment

```bash
# Using STAR aligner
STAR --runThreadN 8 \
  --genomeDir /path/to/star_index \
  --readFilesIn sample_R1.fastq.gz sample_R2.fastq.gz \
  --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix sample_
```

### Part B: WASP Mapping Correction

[Detailed WASP steps...]

### Part C: Gene-Level Allele Counting

```bash
wasp2-count count-variants \
  sample_wasp_filtered.bam \
  genotypes.vcf.gz \
  --samples GM12878 \
  --region gencode.v38.gtf \
  --gene_feature exon \
  --gene_attribute gene_id \
  --out_file gene_counts.tsv
```

**Key options for RNA-seq**:
- `--gene_feature exon`: Count SNPs in exons
- `--gene_attribute gene_id`: Use Ensembl gene IDs
- `--gene_parent transcript_id`: Track which transcript

### Part D: Gene-Level Imbalance Analysis

```bash
wasp2-analyze find-imbalance \
  gene_counts.tsv \
  --min 10 \
  --groupby gene_id \
  --out_file gene_imbalance.tsv
```

**Key options**:
- `--min 10`: Require ≥10 total reads per gene
- `--groupby gene_id`: Aggregate by gene (not transcript)

### Part E: Interpretation

[How to interpret results...]

## Expected Results

- ~15,000 genes with sufficient coverage
- ~500-1,000 genes with significant allelic imbalance (FDR < 0.05)
- Known imprinted genes should show strong imbalance

## Validation

Compare your results to known imprinted genes:
[List of expected imprinted genes...]

## Troubleshooting RNA-seq Specific Issues

**Low coverage genes**: Use `--min 20` for stricter threshold
**Transcript ambiguity**: Add `--use_region_names` with transcript-level analysis
**Multi-mapping reads**: Consider `STAR --outFilterMultimapNmax 1`
```

#### Tutorial 5: ATAC-seq ASE (atac_ase.md)
```markdown
# ATAC-seq Allelic Chromatin Accessibility

**Level**: Intermediate
**Time**: 45 minutes

## Use Case

Measure allele-specific chromatin accessibility in ATAC-seq data to:
- Identify regulatory variants affecting accessibility
- Map allele-specific transcription factor binding
- Compare accessibility between conditions

## Key Differences from RNA-seq

| Aspect | RNA-seq | ATAC-seq |
|--------|---------|----------|
| **Features** | Genes/Transcripts | Peaks/Regions |
| **Annotation** | GTF/GFF | BED/narrowPeak |
| **Coverage** | Exons | Open chromatin |
| **Expected AI** | Imprinting, eQTLs | caQTLs, TF binding |

## Workflow

### Part A: Peak Calling

```bash
# Use MACS2 for peak calling
macs2 callpeak \
  -t sample.bam \
  -f BAMPE \
  -g hs \
  -n sample \
  --outdir peaks/ \
  -q 0.01
```

### Part B: WASP Mapping (Same as RNA-seq)

[WASP steps...]

### Part C: Peak-Level Allele Counting

```bash
wasp2-count count-variants \
  sample_wasp_filtered.bam \
  genotypes.vcf.gz \
  --samples NA12878 \
  --region peaks/sample_peaks.narrowPeak \
  --out_file peak_counts.tsv
```

**Key difference**: Use `narrowPeak` file instead of GTF

### Part D: Peak-Level Analysis

```bash
wasp2-analyze find-imbalance \
  peak_counts.tsv \
  --min 10 \
  --out_file peak_imbalance.tsv
```

### Part E: TF Binding Motif Enrichment

```bash
# Extract imbalanced peaks
awk 'NR==1 || $8 < 0.05' peak_imbalance.tsv > imbalanced_peaks.tsv

# Convert to BED for motif analysis
awk 'NR>1 {print $1"\t"$2-1"\t"$2}' imbalanced_peaks.tsv > imbalanced_peaks.bed

# Run motif enrichment (e.g., HOMER)
findMotifsGenome.pl imbalanced_peaks.bed hg38 motifs/ -size 200
```

## Interpretation

- Peaks with AI likely contain caQTLs
- Look for TF motifs disrupted by variants
- Compare accessibility between haplotypes

## Advanced: Footprinting Analysis

[Integration with footprinting tools...]
```

#### Tutorial 6: Single-Cell Analysis (single_cell.md)
```markdown
# Single-Cell Allele-Specific Analysis

**Level**: Advanced
**Time**: 60 minutes

## Overview

WASP2 provides specialized tools for single-cell data:
- `count-variants-sc`: Per-cell allele counting
- `find-imbalance-sc`: Cell-type-specific imbalance
- `compare-imbalance`: Differential AI between cell types

## Workflow

### Part A: Cell Barcode Preparation

```bash
# Extract cell barcodes from filtered cells (10x Genomics)
zcat filtered_feature_bc_matrix/barcodes.tsv.gz > cell_barcodes.txt
```

### Part B: Single-Cell Allele Counting

```bash
wasp2-count count-variants-sc \
  possorted_genome_bam.bam \
  genotypes.vcf.gz \
  cell_barcodes.txt \
  --samples donor1 \
  --feature peaks.bed \
  --out_file sc_allele_counts.h5ad
```

**Output**: AnnData object (h5ad) with:
- `.X`: Cell × SNP count matrix
- `.var`: SNP annotations
- `.obs`: Cell annotations

### Part C: Cell Type Annotation

Create barcode-to-celltype mapping:
```bash
# Format: BARCODE\tCELLTYPE
# Example:
AAACCTGAGAAACCAT-1    CD4_T
AAACCTGAGAAACCGC-1    CD4_T
AAACCTGAGAAACCTA-1    CD8_T
```

### Part D: Cell-Type-Specific Imbalance

```bash
wasp2-analyze find-imbalance-sc \
  sc_allele_counts.h5ad \
  barcode_celltype_map.tsv \
  --groups CD4_T,CD8_T,B_cell \
  --min 20 \
  --out_file celltype_imbalance.tsv
```

### Part E: Differential AI Between Cell Types

```bash
wasp2-analyze compare-imbalance \
  sc_allele_counts.h5ad \
  barcode_celltype_map.tsv \
  --groups CD4_T,CD8_T \
  --out_file CD4_vs_CD8_imbalance.tsv
```

## Interpretation

[How to interpret single-cell AI results...]

## Visualization in Python

```python
import scanpy as sc
import anndata as ad

# Load results
adata = ad.read_h5ad('sc_allele_counts.h5ad')

# Plot allelic ratio per cell type
sc.pl.violin(adata, 'allelic_ratio', groupby='celltype')
```

## Troubleshooting Single-Cell Issues

**Low SNP coverage**: Single cells have sparse data, use `--min 5` or aggregate
**Too many cells**: Subsample or analyze cell types separately
**Memory issues**: Process chromosomes separately
```

#### Tutorial 8: Troubleshooting Guide (troubleshooting.md)
```markdown
# WASP2 Troubleshooting Guide

Comprehensive guide organized by module and error type.

## General Issues

### Installation Problems

#### Problem: Rust extension fails to build
```
error: failed to run custom build command for `wasp2-rust`
```

**Causes**:
1. Missing Rust compiler
2. Missing libclang
3. Incompatible maturin version

**Solutions**:
```bash
# Install Rust
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# Install libclang (Ubuntu/Debian)
sudo apt-get install libclang-dev

# Install libclang (macOS)
brew install llvm
export LIBCLANG_PATH=$(brew --prefix llvm)/lib

# Update maturin
pip install --upgrade maturin

# Retry build
maturin develop --release -m rust/Cargo.toml
```

[... More troubleshooting sections ...]

## Module-Specific Issues

### Counting Module

#### No output SNPs

**Symptoms**: `counts.tsv` is empty or has only header

**Diagnostic**:
```bash
# Check VCF has heterozygous SNPs for your sample
bcftools view -s sample1 -g ^0/0,^1/1,^./. variants.vcf.gz | head -20

# Check BAM has reads
samtools view -c sample.bam

# Check coordinate overlap
samtools view sample.bam chr10:1000000-2000000 | head
bcftools view -r chr10:1000000-2000000 variants.vcf.gz | head
```

**Solutions**:
1. Verify sample name: `bcftools query -l variants.vcf.gz`
2. Check chromosome naming (chr10 vs 10)
3. Ensure BAM and VCF use same reference genome

[... More troubleshooting ...]

## Performance Issues

### Slow VCF Parsing

**Symptoms**: Counting takes >1 hour on large VCF

**Solutions**:
1. Install cyvcf2: `pip install wasp2[cyvcf2]` (7x speedup)
2. Convert to BCF: `bcftools view -O b variants.vcf.gz > variants.bcf` (5-8x speedup)
3. Convert to PGEN: `plink2 --vcf variants.vcf.gz --make-pgen` (25x speedup)

### High Memory Usage

**Symptoms**: Process killed with "Out of memory"

**Solutions**:
1. Process chromosomes separately: `--region chr10.bed`
2. Reduce threads: `--threads 1`
3. Use PGEN format instead of VCF (lower memory)
4. Filter VCF to heterozygous SNPs first:
   ```bash
   bcftools view -s sample1 -g ^0/0,^1/1 input.vcf.gz -O z -o het_only.vcf.gz
   ```

[... More performance tips ...]

## Error Messages Reference

| Error | Module | Cause | Solution |
|-------|--------|-------|----------|
| `FileNotFoundError: variants.vcf.gz.tbi` | count | Missing VCF index | Run `bcftools index variants.vcf.gz` |
| `ValueError: Sample not found in VCF` | count | Wrong sample name | Check with `bcftools query -l` |
| `RuntimeError: BAM file not sorted` | count | Unsorted BAM | Run `samtools sort` |
| `OSError: [Errno 28] No space left` | All | Disk full | Clean temp files or use `--temp_loc` |

[... Complete error reference ...]
```

#### Tutorial 9: Performance Tuning (performance_tuning.md)
```markdown
# WASP2 Performance Optimization

Get maximum performance from WASP2 for large-scale analyses.

## Variant Format Selection

### Performance Comparison

| Format | Read Speed | Memory | Recommendation |
|--------|------------|--------|----------------|
| VCF.gz (pysam) | 1x | Medium | Default, testing |
| VCF.gz (cyvcf2) | 7x | Medium | Production |
| BCF | 5-8x | Medium | Good balance |
| PGEN | 25x | Low | Large cohorts |

### When to Use Each Format

**VCF.gz + cyvcf2**:
- Best for most production workflows
- Preserves all VCF fields
- Compatible with all tools
- `pip install wasp2[cyvcf2]`

**BCF**:
- Binary VCF with no information loss
- Faster than VCF.gz
- Use when sharing with collaborators who have bcftools

**PGEN**:
- Best for genotype-only workflows
- Lowest memory usage
- 25x faster I/O
- Use for large cohorts (>1000 samples)

### Format Conversion

```bash
# VCF to BCF
bcftools view -O b variants.vcf.gz > variants.bcf
bcftools index variants.bcf

# VCF to PGEN
plink2 --vcf variants.vcf.gz \
  --make-pgen \
  --out variants

# PGEN back to VCF (if needed)
plink2 --pfile variants \
  --export vcf bgz \
  --out variants_from_pgen
```

## Threading and Parallelization

### Optimal Thread Counts

```bash
# Counting module (Rust-accelerated)
wasp2-count count-variants sample.bam variants.pgen --threads 4

# Mapping module
wasp2-map filter-remapped remap.bam --threads 4

# Analysis module (Python)
# Single-threaded optimization is sufficient
```

**Guidelines**:
- Use threads ≤ physical cores
- Diminishing returns beyond 8 threads
- I/O bottleneck often limits scaling

## Memory Optimization

### Large VCF Files

```bash
# Problem: 100GB VCF file causes OOM
# Solution 1: Convert to PGEN (lower memory)
plink2 --vcf huge.vcf.gz --make-pgen --out huge

# Solution 2: Process by chromosome
for chr in {1..22} X Y; do
  wasp2-count count-variants sample.bam huge.vcf.gz \
    --region chr${chr}.bed \
    --out_file counts_chr${chr}.tsv
done

# Combine results
head -1 counts_chr1.tsv > all_counts.tsv
tail -n +2 -q counts_chr*.tsv >> all_counts.tsv
```

### Large BAM Files

```bash
# Enable Rust acceleration (lower memory footprint)
export WASP2_USE_RUST=1

# Process regions separately
bedtools makewindows -g genome.txt -w 10000000 > windows.bed
parallel -j 4 wasp2-count count-variants sample.bam variants.vcf.gz \
  --region {} --out_file {/.}.tsv ::: windows_*.bed
```

## Disk I/O Optimization

### Temporary File Location

```bash
# Use fast local SSD instead of network storage
export TMPDIR=/scratch/local/tmp

# Or specify in command
wasp2-count count-variants sample.bam variants.vcf.gz \
  --temp_loc /scratch/local/tmp
```

### Pre-computed Intermediate Files

```bash
# Skip VCF-to-BED conversion on repeated runs
wasp2-count count-variants sample.bam variants.vcf.gz \
  --vcf-bed precomputed_vcf.bed \
  --intersect-bed precomputed_intersect.bed
```

## Pipeline Parallelization

### Processing Multiple Samples

```bash
# GNU parallel for multiple samples
parallel -j 4 \
  wasp2-count count-variants {}.bam variants.pgen \
    --samples {} \
    --out_file {}_counts.tsv \
  ::: sample1 sample2 sample3 sample4

# Nextflow pipeline (example)
process count_alleles {
  input:
    tuple val(sample_id), path(bam), path(bai)

  output:
    path("${sample_id}_counts.tsv")

  script:
  """
  wasp2-count count-variants ${bam} ${params.vcf} \
    --samples ${sample_id} \
    --out_file ${sample_id}_counts.tsv
  """
}
```

## Benchmark Results

### Real-World Performance

**Dataset**: 1000 Genomes, 30x WGS, ~100M variants

| Configuration | Time | Memory |
|---------------|------|--------|
| VCF.gz (pysam) | 45 min | 8 GB |
| VCF.gz (cyvcf2) | 6.5 min | 8 GB |
| BCF | 8 min | 8 GB |
| PGEN | 1.8 min | 4 GB |

**Recommendation**: Use PGEN for >10M variants, cyvcf2 otherwise

## Profiling Your Workflow

```bash
# Time each step
time wasp2-count count-variants sample.bam variants.vcf.gz

# Profile memory usage
/usr/bin/time -v wasp2-count count-variants sample.bam variants.vcf.gz

# Identify bottlenecks with Python profiler
python -m cProfile -o profile.stats count_script.py
python -c "import pstats; p = pstats.Stats('profile.stats'); p.sort_stats('cumulative').print_stats(20)"
```

## Cloud Computing Optimization

### AWS Batch / Google Cloud

```bash
# Use instance storage for temp files
export TMPDIR=/mnt/local-ssd

# Download data to local storage first
aws s3 cp s3://bucket/sample.bam /mnt/local-ssd/
wasp2-count count-variants /mnt/local-ssd/sample.bam variants.pgen

# Upload results
aws s3 cp counts.tsv s3://bucket/results/
```

### HPC Clusters

```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00

module load python/3.10
module load rust/1.70

wasp2-count count-variants sample.bam variants.pgen \
  --threads 8 \
  --temp_loc $TMPDIR \
  --out_file counts.tsv
```

## Summary Recommendations

1. **Always use cyvcf2 or PGEN** for production
2. **Process by chromosome** for very large files
3. **Use local SSD** for temp files
4. **Enable Rust acceleration** (default in v1.2+)
5. **Parallelize across samples**, not within sample
6. **Pre-filter VCF** to heterozygous SNPs only
```

---

## 3. CLI Documentation Best Practices

### 3.1 Enhanced --help Output

#### Current State
WASP2 uses Typer which generates decent help automatically.

#### Recommended Improvements

**Structure for each command**:

```
Usage: wasp2-count count-variants [OPTIONS] BAM VARIANTS

  Count allele-specific reads at heterozygous SNP positions.

  This command quantifies the number of reads supporting each allele (reference
  vs. alternate) at heterozygous SNPs. Results can be filtered by sample genotype
  and annotated with genomic regions (genes, peaks).

  Examples:
    # Basic counting
    wasp2-count count-variants sample.bam variants.vcf.gz

    # With sample filtering and gene annotation
    wasp2-count count-variants sample.bam variants.vcf.gz \
      --samples NA12878 \
      --region genes.gtf \
      --out_file counts.tsv

    # Using high-performance PGEN format
    wasp2-count count-variants sample.bam variants.pgen \
      --samples NA12878 \
      --out_file counts.tsv

Arguments:
  BAM        Path to aligned reads (BAM format, must be sorted and indexed)
  VARIANTS   Path to variants (VCF, BCF, or PGEN format)

Options:
  Input Filtering:
    -s, --samples TEXT         Sample ID(s) to filter heterozygous SNPs
                               Accepts: sample1,sample2 or file with one ID per line
    -r, --region PATH          Filter SNPs overlapping regions
                               Accepts: BED, GTF, GFF3, narrowPeak formats

  Output:
    -o, --out_file PATH        Output file path [default: counts.tsv]
    --temp_loc PATH            Directory for intermediate files [default: system temp]

  Region Annotation (for GTF/GFF3):
    --gene_feature TEXT        Feature type to count [default: exon]
    --gene_attribute TEXT      Attribute for feature ID [default: gene_id]
    --gene_parent TEXT         Parent attribute [default: transcript_id]
    --use_region_names         Use region names instead of coordinates

  Performance:
    --use-rust / --no-rust     Enable Rust acceleration [default: use-rust]
    --include-indels           Include indels in addition to SNPs

  Advanced:
    --vcf-bed PATH            Pre-computed VCF BED file (skip conversion)
    --intersect-bed PATH      Pre-computed intersect BED file (skip intersection)

  Other:
    -h, --help                Show this message and exit
    --version                 Show version and exit

Output Format:
  Tab-separated file with columns:
    chr, pos, ref, alt          - Variant coordinates and alleles
    ref_count, alt_count        - Reads supporting each allele
    other_count                 - Reads with other alleles
    total_count                 - Total overlapping reads
    region (if --region used)   - Overlapping gene/peak

Performance Tips:
  - Use PGEN format for large variant files (25x faster I/O)
  - Install cyvcf2 for faster VCF parsing: pip install wasp2[cyvcf2]
  - Process chromosomes separately for very large files

See Also:
  wasp2-analyze find-imbalance  - Detect allelic imbalance from counts
  wasp2-map make-reads          - Generate reads for WASP mapping

  Full documentation: https://jaureguy760.github.io/WASP2-exp/
```

#### Implementation

Enhance Typer command docstrings:

```python
@app.command()
def count_variants(
    bam: Annotated[str, typer.Argument(
        help="Path to aligned reads (BAM format, must be sorted and indexed)",
        metavar="BAM"
    )],
    variants: Annotated[str, typer.Argument(
        help="Path to variants (VCF, BCF, or PGEN format)",
        metavar="VARIANTS"
    )],
    # ... rest of parameters
) -> None:
    """
    Count allele-specific reads at heterozygous SNP positions.

    This command quantifies the number of reads supporting each allele (reference
    vs. alternate) at heterozygous SNPs. Results can be filtered by sample genotype
    and annotated with genomic regions (genes, peaks).

    \b
    Examples:
      # Basic counting
      wasp2-count count-variants sample.bam variants.vcf.gz

      # With sample filtering and gene annotation
      wasp2-count count-variants sample.bam variants.vcf.gz \\
        --samples NA12878 \\
        --region genes.gtf \\
        --out_file counts.tsv

    \b
    Output Format:
      Tab-separated file with columns:
        chr, pos, ref, alt - Variant coordinates
        ref_count, alt_count - Read counts per allele

    \b
    Performance Tips:
      - Use PGEN format for 25x faster I/O
      - Install cyvcf2: pip install wasp2[cyvcf2]

    See: https://jaureguy760.github.io/WASP2-exp/user_guide/counting.html
    """
```

### 3.2 Man Pages

Create traditional Unix man pages for each command.

#### File Structure
```
man/
├── man1/
│   ├── wasp2.1                    # Main command overview
│   ├── wasp2-count.1              # Count module overview
│   ├── wasp2-count-variants.1     # Specific command
│   ├── wasp2-count-variants-sc.1
│   ├── wasp2-map.1
│   ├── wasp2-map-make-reads.1
│   ├── wasp2-map-filter-remapped.1
│   ├── wasp2-analyze.1
│   └── wasp2-analyze-find-imbalance.1
```

#### Example Man Page (wasp2-count-variants.1)

```nroff
.TH WASP2-COUNT-VARIANTS 1 "January 2025" "WASP2 1.2.1" "WASP2 Manual"

.SH NAME
wasp2-count-variants \- count allele-specific reads at heterozygous SNPs

.SH SYNOPSIS
.B wasp2-count count-variants
.RI [ OPTIONS ]
.I BAM VARIANTS

.SH DESCRIPTION
.B wasp2-count count-variants
quantifies allele-specific read counts at heterozygous single nucleotide
polymorphism (SNP) positions. It processes aligned reads from a BAM file
and variant calls from a VCF/BCF/PGEN file to count reads supporting each
allele.

This is typically the first step in allelic imbalance analysis, followed
by statistical testing with
.BR wasp2-analyze (1).

.SH ARGUMENTS
.TP
.I BAM
Path to aligned reads in BAM format. Must be coordinate-sorted and indexed
(i.e., .bai file must exist).

.TP
.I VARIANTS
Path to variant calls. Supports VCF (.vcf, .vcf.gz), BCF (.bcf), and
PLINK2 PGEN (.pgen) formats. VCF/BCF files should be indexed (.tbi or .csi).

.SH OPTIONS
.SS Input Filtering
.TP
.BR \-s ", " \-\-samples =\fISAMPLE\fR
Filter SNPs to those heterozygous in the specified sample(s). Accepts
comma-separated sample IDs or a file with one sample per line.

.TP
.BR \-r ", " \-\-region =\fIPATH\fR
Filter SNPs overlapping genomic regions. Accepts BED, GTF, GFF3, or
narrowPeak format files.

.SS Output
.TP
.BR \-o ", " \-\-out_file =\fIPATH\fR
Output file path.
.I Default:
counts.tsv

.TP
.BR \-\-temp_loc =\fIDIR\fR
Directory for intermediate files. If not specified, uses system temporary
directory and removes files after completion.

.SS Region Annotation
.TP
.BR \-\-gene_feature =\fITYPE\fR
Feature type from GTF/GFF3 to use for counting.
.I Default:
exon

.TP
.BR \-\-gene_attribute =\fINAME\fR
Attribute name for feature identifier.
.I Default:
gene_id (GTF) or ID (GFF3)

.TP
.BR \-\-use_region_names
Use region names instead of coordinates in output. Names taken from
4th column of BED files.

.SS Performance
.TP
.BR \-\-use\-rust / \-\-no\-rust
Enable or disable Rust acceleration.
.I Default:
use-rust

.TP
.BR \-\-include\-indels
Include insertion/deletion variants in addition to SNPs.

.SH OUTPUT FORMAT
Tab-separated file with the following columns:

.TP
.B chr
Chromosome name

.TP
.B pos
SNP position (1-based)

.TP
.B ref
Reference allele

.TP
.B alt
Alternate allele

.TP
.B ref_count
Number of reads supporting reference allele

.TP
.B alt_count
Number of reads supporting alternate allele

.TP
.B other_count
Number of reads with other alleles

.TP
.B region
Overlapping genomic region (if --region specified)

.SH EXAMPLES
Basic counting:
.PP
.RS
.nf
wasp2-count count-variants sample.bam variants.vcf.gz
.fi
.RE

Count heterozygous SNPs for specific sample:
.PP
.RS
.nf
wasp2-count count-variants sample.bam variants.vcf.gz \\
  --samples NA12878 \\
  --out_file counts.tsv
.fi
.RE

Annotate with gene regions:
.PP
.RS
.nf
wasp2-count count-variants rnaseq.bam variants.pgen \\
  --samples NA12878 \\
  --region gencode.v38.gtf \\
  --out_file gene_counts.tsv
.fi
.RE

ATAC-seq with peak annotation:
.PP
.RS
.nf
wasp2-count count-variants atac.bam variants.bcf \\
  --samples NA12878 \\
  --region peaks.narrowPeak \\
  --out_file peak_counts.tsv
.fi
.RE

.SH EXIT STATUS
.TP
.B 0
Success

.TP
.B 1
General error (missing files, invalid arguments)

.TP
.B 2
Data processing error (empty output, incompatible formats)

.SH ENVIRONMENT
.TP
.B WASP2_DISABLE_RUST
Set to 1 to disable Rust acceleration (use Python fallback)

.TP
.B TMPDIR
Directory for temporary files if --temp_loc not specified

.SH FILES
.TP
.I counts.tsv
Default output filename if --out_file not specified

.SH NOTES
.SS Performance Optimization
For large variant files (>10M variants), use PGEN format for ~25x speedup:
.PP
.RS
.nf
plink2 --vcf variants.vcf.gz --make-pgen --out variants
wasp2-count count-variants sample.bam variants.pgen
.fi
.RE

Alternatively, install cyvcf2 for ~7x faster VCF parsing:
.PP
.RS
.nf
pip install wasp2[cyvcf2]
.fi
.RE

.SS Reference Genome Compatibility
Ensure BAM and VCF files use the same reference genome build (e.g., both
GRCh38 or both hg19). Chromosome naming (chr10 vs 10) must also match.

.SH BUGS
Report bugs at https://github.com/Jaureguy760/WASP2-exp/issues

.SH SEE ALSO
.BR wasp2 (1),
.BR wasp2-analyze (1),
.BR wasp2-map (1),
.BR samtools (1),
.BR bcftools (1)

Full documentation:
.UR https://jaureguy760.github.io/WASP2-exp/
.UE

.SH AUTHORS
Aaron Ho, Jeff Jaureguy, Graham McVicker

.SH COPYRIGHT
Copyright \(co 2025 Aaron Ho, Jeff Jaureguy, McVicker Lab
.br
License: MIT
```

#### Installation

Add to `setup.py` or `pyproject.toml`:

```toml
[tool.setuptools]
data_files = [
    ("share/man/man1", [
        "man/man1/wasp2.1",
        "man/man1/wasp2-count-variants.1",
        # ... other man pages
    ])
]
```

### 3.3 Shell Completion Scripts

Provide tab completion for bash, zsh, fish.

#### Generate with Typer

```python
# scripts/generate_completions.py
import typer
from counting.__main__ import app as count_app
from mapping.__main__ import app as map_app
from analysis.__main__ import app as analysis_app

def generate_all_completions():
    """Generate shell completions for all WASP2 commands"""

    # Create main app
    main_app = typer.Typer()
    main_app.add_typer(count_app, name="count")
    main_app.add_typer(map_app, name="map")
    main_app.add_typer(analysis_app, name="analyze")

    # Generate completions
    for shell in ["bash", "zsh", "fish"]:
        completion = typer.completion.get_completion(main_app, shell=shell)
        output_file = f"completions/wasp2.{shell}"
        with open(output_file, "w") as f:
            f.write(completion)
        print(f"Generated {output_file}")

if __name__ == "__main__":
    generate_all_completions()
```

#### Installation Instructions (in README)

```markdown
### Shell Completion (Optional)

Enable tab completion for WASP2 commands:

**Bash**:
```bash
# Add to ~/.bashrc
eval "$(wasp2-count --show-completion bash)"
eval "$(wasp2-map --show-completion bash)"
eval "$(wasp2-analyze --show-completion bash)"

# Or install completion script
sudo cp completions/wasp2.bash /etc/bash_completion.d/wasp2
```

**Zsh**:
```bash
# Add to ~/.zshrc
eval "$(wasp2-count --show-completion zsh)"
eval "$(wasp2-map --show-completion zsh)"
eval "$(wasp2-analyze --show-completion zsh)"
```

**Fish**:
```bash
wasp2-count --show-completion fish > ~/.config/fish/completions/wasp2-count.fish
wasp2-map --show-completion fish > ~/.config/fish/completions/wasp2-map.fish
wasp2-analyze --show-completion fish > ~/.config/fish/completions/wasp2-analyze.fish
```
```

### 3.4 Example Commands Reference

Create `examples/` directory with common use cases.

```
examples/
├── README.md                      # Overview of all examples
├── basic_rnaseq.sh               # Basic RNA-seq ASE
├── basic_atacseq.sh              # Basic ATAC-seq ASE
├── full_pipeline.sh              # Complete WASP pipeline
├── single_cell.sh                # Single-cell workflow
├── multiple_samples.sh           # Batch processing
├── performance_optimized.sh      # Performance tuning
└── data/                         # Small test datasets
    ├── sample.bam
    ├── variants.vcf.gz
    └── genes.gtf
```

#### Example: examples/basic_rnaseq.sh

```bash
#!/bin/bash
# WASP2 Example: Basic RNA-seq Allele-Specific Expression Analysis
#
# This script demonstrates a complete RNA-seq ASE workflow using WASP2.
# Expected runtime: ~5 minutes on test data

set -euo pipefail  # Exit on error, undefined variables, pipe failures

# ==============================================================================
# Configuration
# ==============================================================================

# Input files (update paths for your data)
BAM="data/rnaseq_sample.bam"
VCF="data/genotypes.vcf.gz"
GTF="data/genes.gtf"
SAMPLE_ID="NA12878"

# Output directory
OUTDIR="results/rnaseq_ase"
mkdir -p "$OUTDIR"

# ==============================================================================
# Step 1: Quality Control
# ==============================================================================

echo "==> Step 1: Quality Control"

# Check BAM alignment statistics
samtools flagstat "$BAM" > "$OUTDIR/alignment_stats.txt"

# Check variant file
echo "Total variants: $(bcftools view -H "$VCF" | wc -l)"
echo "Het SNPs for $SAMPLE_ID: $(bcftools view -s "$SAMPLE_ID" -g ^0/0,^1/1 "$VCF" | wc -l)"

# ==============================================================================
# Step 2: Count Allele-Specific Reads
# ==============================================================================

echo "==> Step 2: Counting allele-specific reads"

wasp2-count count-variants \
  "$BAM" \
  "$VCF" \
  --samples "$SAMPLE_ID" \
  --region "$GTF" \
  --gene_feature exon \
  --gene_attribute gene_id \
  --out_file "$OUTDIR/gene_counts.tsv"

# Inspect output
echo "Counted SNPs in $(tail -n +2 "$OUTDIR/gene_counts.tsv" | wc -l) genes"
head "$OUTDIR/gene_counts.tsv"

# ==============================================================================
# Step 3: Detect Allelic Imbalance
# ==============================================================================

echo "==> Step 3: Statistical analysis for allelic imbalance"

wasp2-analyze find-imbalance \
  "$OUTDIR/gene_counts.tsv" \
  --min 10 \
  --groupby gene_id \
  --out_file "$OUTDIR/gene_imbalance.tsv"

# Summary statistics
echo "Genes tested: $(tail -n +2 "$OUTDIR/gene_imbalance.tsv" | wc -l)"
echo "Significant genes (FDR < 0.05): $(awk 'NR>1 && $8 < 0.05' "$OUTDIR/gene_imbalance.tsv" | wc -l)"

# ==============================================================================
# Step 4: Extract Significant Results
# ==============================================================================

echo "==> Step 4: Extracting significant genes"

# Genes with significant allelic imbalance
awk 'NR==1 || $8 < 0.05' "$OUTDIR/gene_imbalance.tsv" \
  > "$OUTDIR/significant_genes.tsv"

# Sort by effect size
sort -t$'\t' -k6,6nr "$OUTDIR/significant_genes.tsv" \
  > "$OUTDIR/significant_genes_sorted.tsv"

echo "Top 10 genes with strongest allelic imbalance:"
head -11 "$OUTDIR/significant_genes_sorted.tsv" | column -t

# ==============================================================================
# Complete
# ==============================================================================

echo ""
echo "==> Analysis complete!"
echo "Results in: $OUTDIR/"
echo "  - gene_counts.tsv: Raw allele counts"
echo "  - gene_imbalance.tsv: Statistical test results"
echo "  - significant_genes.tsv: FDR < 0.05 genes"
echo ""
echo "Next steps:"
echo "  1. Visualize results (see examples/plot_results.R)"
echo "  2. Compare with known imprinted genes"
echo "  3. Perform gene set enrichment analysis"
```

---

## 4. API Documentation Best Practices

### 4.1 Docstring Standards

#### Recommendation: Google Style
WASP2's Sphinx is already configured for Google docstrings. This style is:
- More readable than NumPy style for shorter functions
- Well-supported by Sphinx with napoleon extension
- Popular in bioinformatics (used by scanpy, seaborn, etc.)

#### Comprehensive Docstring Template

```python
def run_count_variants(
    bam_file: str,
    variant_file: str,
    region_file: Optional[str] = None,
    samples: Optional[str] = None,
    use_region_names: bool = False,
    out_file: Optional[str] = None,
    temp_loc: Optional[str] = None,
    gene_feature: Optional[str] = None,
    gene_attribute: Optional[str] = None,
    gene_parent: Optional[str] = None,
    use_rust: bool = True,
    precomputed_vcf_bed: Optional[str] = None,
    precomputed_intersect: Optional[str] = None,
    include_indels: bool = False
) -> None:
    """Count allele-specific reads at heterozygous SNP positions.

    Quantifies the number of reads supporting reference vs. alternate alleles
    at heterozygous single nucleotide polymorphisms (SNPs). This is the first
    step in allelic imbalance analysis.

    The function processes aligned reads from a BAM file and variant calls from
    a VCF/BCF/PGEN file. Results can be filtered by sample genotype and annotated
    with genomic regions (genes, ATAC-seq peaks, etc.).

    Args:
        bam_file: Path to aligned reads (BAM format). Must be coordinate-sorted
            and indexed (.bai file required).
        variant_file: Path to variant calls. Supports VCF (.vcf, .vcf.gz),
            BCF (.bcf), and PLINK2 PGEN (.pgen) formats. VCF/BCF files should
            be indexed (.tbi or .csi).
        region_file: Path to genomic regions for SNP filtering and annotation.
            Accepts BED, GTF, GFF3, or narrowPeak formats. If provided, only
            SNPs overlapping these regions are counted. Optional.
        samples: Sample ID(s) to filter heterozygous SNPs. Accepts comma-separated
            IDs (e.g., "sample1,sample2") or path to file with one ID per line.
            If not provided, all variants are used regardless of genotype. Optional.
        use_region_names: If True, use region names (4th column of BED file) in
            output instead of genomic coordinates. Ignored if region_file is not
            BED format. Default: False.
        out_file: Output file path for allele counts. Tab-separated format with
            columns: chr, pos, ref, alt, ref_count, alt_count, other_count.
            Default: "counts.tsv".
        temp_loc: Directory for intermediate files. If None, uses system temporary
            directory and removes files after completion. Specify a path to preserve
            intermediate files for debugging. Optional.
        gene_feature: Feature type from GTF/GFF3 to use for SNP counting (e.g.,
            "exon", "CDS"). Only relevant if region_file is GTF/GFF3 format.
            Default: "exon".
        gene_attribute: Attribute name for feature identifier in GTF/GFF3 files
            (e.g., "gene_id", "transcript_id"). Default: "gene_id" for GTF,
            "ID" for GFF3.
        gene_parent: Parent attribute for hierarchical features in GTF/GFF3
            (e.g., "transcript_id" for exons). Default: "transcript_id" for GTF,
            "Parent" for GFF3.
        use_rust: If True, use Rust-accelerated counting (requires wasp2_rust
            extension). Falls back to Python if Rust extension not available.
            Default: True.
        precomputed_vcf_bed: Path to pre-computed VCF BED file to skip variant
            file conversion step. Useful for repeated runs on same variant file.
            Optional.
        precomputed_intersect: Path to pre-computed intersection BED file to skip
            bedtools intersect step. Useful for repeated runs. Optional.
        include_indels: If True, include insertion/deletion variants in addition
            to SNPs. Default: False (SNPs only).

    Returns:
        None. Results written to out_file.

    Raises:
        FileNotFoundError: If bam_file, variant_file, or region_file does not exist.
        ValueError: If sample ID not found in variant file, or if region_file
            format cannot be determined.
        RuntimeError: If BAM file is not sorted or indexed, or if Rust extension
            fails and use_rust=True.
        IOError: If output file cannot be written (e.g., permission denied).

    Examples:
        Basic counting:

        >>> run_count_variants(
        ...     bam_file="sample.bam",
        ...     variant_file="variants.vcf.gz",
        ...     out_file="counts.tsv"
        ... )

        RNA-seq with gene annotation:

        >>> run_count_variants(
        ...     bam_file="rnaseq.bam",
        ...     variant_file="genotypes.pgen",
        ...     region_file="genes.gtf",
        ...     samples="NA12878",
        ...     gene_feature="exon",
        ...     gene_attribute="gene_id",
        ...     out_file="gene_counts.tsv"
        ... )

        ATAC-seq with peak annotation:

        >>> run_count_variants(
        ...     bam_file="atac.bam",
        ...     variant_file="variants.bcf",
        ...     region_file="peaks.narrowPeak",
        ...     samples="NA12878",
        ...     out_file="peak_counts.tsv"
        ... )

    Notes:
        Performance Tips:
        - Use PGEN format for large variant files (>10M variants, ~25x speedup)
        - Install cyvcf2 for faster VCF parsing: pip install wasp2[cyvcf2]
        - Process chromosomes separately for very large datasets
        - Use precomputed_vcf_bed and precomputed_intersect for repeated runs

        Memory Usage:
        - Typical: 2-8 GB for whole-genome data
        - Use PGEN format to reduce memory footprint
        - Process by chromosome if encountering memory issues

        Reference Genome Compatibility:
        - BAM and variant file must use same reference genome build
        - Chromosome naming must match (chr10 vs 10)
        - Use samtools view and bcftools view to verify

    See Also:
        run_ai_analysis: Detect allelic imbalance from count data
        run_make_remap_reads: Generate reads for WASP mapping

    References:
        van de Geijn et al. (2015). WASP: allele-specific software for robust
        molecular quantitative trait locus discovery. Nature Methods 12:1061-1063.
        https://doi.org/10.1038/nmeth.3582
    """
    # Implementation...
```

### 4.2 Type Hints for Documentation

#### Current State
WASP2 has type hints in function signatures. Sphinx autodoc_typehints is enabled.

#### Best Practices

```python
from typing import Optional, Union, List, Tuple, Dict, Any
from pathlib import Path
from dataclasses import dataclass

# Use Path for file paths
def count_variants(
    bam_file: Union[str, Path],
    variant_file: Union[str, Path],
    *,  # Force keyword arguments
    region_file: Optional[Union[str, Path]] = None,
    samples: Optional[Union[str, List[str]]] = None,
    out_file: Optional[Union[str, Path]] = None,
) -> None:
    """Count alleles with type-safe interface."""
    pass

# Use dataclasses for structured returns
@dataclass
class CountResult:
    """Results from allele counting.

    Attributes:
        n_variants: Total variants processed
        n_het_snps: Heterozygous SNPs counted
        n_regions: Genomic regions overlapped
        output_file: Path to output file
        warnings: List of warning messages
    """
    n_variants: int
    n_het_snps: int
    n_regions: int
    output_file: Path
    warnings: List[str]

def count_variants_typed(...) -> CountResult:
    """Count alleles with structured return."""
    # ...
    return CountResult(
        n_variants=1000,
        n_het_snps=500,
        n_regions=200,
        output_file=Path("counts.tsv"),
        warnings=[]
    )

# Use TypedDict for dictionary returns
from typing import TypedDict

class VariantDict(TypedDict):
    """Variant information dictionary.

    Keys:
        chrom: Chromosome name
        pos: Position (1-based)
        ref: Reference allele
        alt: Alternate allele
        genotype: Sample genotype (0/1, 1/0, etc.)
    """
    chrom: str
    pos: int
    ref: str
    alt: str
    genotype: str

def get_variant(vcf_file: str, index: int) -> VariantDict:
    """Get variant by index with typed return."""
    pass
```

### 4.3 Sphinx Documentation Structure

#### Recommended Structure

```
docs/
├── source/
│   ├── index.rst                    # Landing page
│   ├── installation.rst             # Installation guide
│   ├── quickstart.rst               # 5-min tutorial
│   ├── concepts.rst                 # Background concepts (NEW)
│   │
│   ├── tutorials/                   # Tutorial documentation (NEW)
│   │   ├── index.rst
│   │   ├── basic_workflow.rst
│   │   ├── rnaseq_ase.rst
│   │   ├── atacseq_ase.rst
│   │   ├── single_cell.rst
│   │   └── troubleshooting.rst
│   │
│   ├── user_guide/                  # Existing user guides
│   │   ├── counting.rst
│   │   ├── mapping.rst
│   │   └── analysis.rst
│   │
│   ├── how_to/                      # Task-oriented guides (NEW)
│   │   ├── index.rst
│   │   ├── process_multiple_samples.rst
│   │   ├── optimize_performance.rst
│   │   ├── integrate_with_pipelines.rst
│   │   └── interpret_results.rst
│   │
│   ├── api/                         # API reference
│   │   ├── index.rst
│   │   ├── counting.rst
│   │   ├── mapping.rst
│   │   ├── analysis.rst
│   │   └── io.rst                   # I/O modules (NEW)
│   │
│   ├── cli/                         # CLI reference (NEW)
│   │   ├── index.rst
│   │   ├── wasp2_count.rst
│   │   ├── wasp2_map.rst
│   │   └── wasp2_analyze.rst
│   │
│   ├── explanations/                # Background/theory (NEW)
│   │   ├── index.rst
│   │   ├── allelic_imbalance.rst
│   │   ├── reference_bias.rst
│   │   ├── wasp_algorithm.rst
│   │   └── statistical_models.rst
│   │
│   ├── data_formats/                # Format specifications (NEW)
│   │   ├── index.rst
│   │   ├── input_formats.rst
│   │   ├── output_formats.rst
│   │   └── variant_formats.rst
│   │
│   ├── changelog.rst                # Version history
│   ├── development.rst              # Developer guide
│   ├── faq.rst                      # FAQ (NEW)
│   └── citation.rst                 # How to cite (NEW)
│
├── VCF_PERFORMANCE.md               # Existing performance doc
├── PLINK2_INTEGRATION_DESIGN.md     # Existing design doc
└── examples/                        # Code examples (NEW)
    └── notebooks/
        ├── basic_analysis.ipynb
        ├── rnaseq_workflow.ipynb
        └── visualization.ipynb
```

#### Example: CLI Reference Page (cli/wasp2_count.rst)

```rst
wasp2-count
===========

Command-line interface for the WASP2 counting module.

.. contents:: Commands
   :local:
   :depth: 2

Overview
--------

The ``wasp2-count`` command quantifies allele-specific read counts at
heterozygous SNP positions. It provides two subcommands:

* ``count-variants`` - Count alleles in bulk sequencing data
* ``count-variants-sc`` - Count alleles in single-cell data

Global Options
--------------

.. option:: --help

   Show help message and exit

.. option:: --version

   Show version number and exit

count-variants
--------------

Count allele-specific reads at heterozygous SNPs in bulk data.

Synopsis
~~~~~~~~

.. code-block:: bash

   wasp2-count count-variants [OPTIONS] BAM VARIANTS

Arguments
~~~~~~~~~

.. option:: BAM

   Path to aligned reads (BAM format). Must be sorted and indexed.

.. option:: VARIANTS

   Path to variants (VCF, BCF, or PGEN format).

Options
~~~~~~~

Input Filtering
^^^^^^^^^^^^^^^

.. option:: -s <SAMPLE>, --samples <SAMPLE>

   Sample ID(s) to filter heterozygous SNPs.

   Accepts:
      - Comma-separated list: ``-s sample1,sample2``
      - File with one sample per line: ``-s samples.txt``

   If not provided, all variants are used.

.. option:: -r <PATH>, --region <PATH>

   Filter SNPs overlapping genomic regions.

   Accepts:
      - BED format (``.bed``)
      - GTF format (``.gtf``)
      - GFF3 format (``.gff``, ``.gff3``)
      - narrowPeak format (``.narrowPeak``)

Output
^^^^^^

.. option:: -o <PATH>, --out_file <PATH>

   Output file path. Default: ``counts.tsv``

.. option:: --temp_loc <DIR>

   Directory for intermediate files. If not specified, uses system
   temporary directory and removes files after completion.

Region Annotation
^^^^^^^^^^^^^^^^^

.. option:: --gene_feature <TYPE>

   Feature type from GTF/GFF3 to count overlapping SNPs.
   Default: ``exon``

   Examples: ``exon``, ``CDS``, ``five_prime_UTR``

.. option:: --gene_attribute <NAME>

   Attribute name for feature identifier.
   Default: ``gene_id`` (GTF), ``ID`` (GFF3)

.. option:: --gene_parent <NAME>

   Parent attribute for hierarchical features.
   Default: ``transcript_id`` (GTF), ``Parent`` (GFF3)

.. option:: --use_region_names

   Use region names (4th BED column) instead of coordinates in output.

Performance
^^^^^^^^^^^

.. option:: --use-rust / --no-rust

   Enable or disable Rust acceleration. Default: ``--use-rust``

.. option:: --include-indels

   Include indels in addition to SNPs. Default: SNPs only

Advanced
^^^^^^^^

.. option:: --vcf-bed <PATH>

   Pre-computed VCF BED file (skip variant conversion)

.. option:: --intersect-bed <PATH>

   Pre-computed intersect BED file (skip intersection)

Examples
--------

Basic Counting
~~~~~~~~~~~~~~

Count alleles at all variants:

.. code-block:: bash

   wasp2-count count-variants sample.bam variants.vcf.gz

Filter by Sample
~~~~~~~~~~~~~~~~

Count only heterozygous SNPs for specific sample:

.. code-block:: bash

   wasp2-count count-variants sample.bam variants.vcf.gz \
     --samples NA12878 \
     --out_file counts.tsv

RNA-seq with Genes
~~~~~~~~~~~~~~~~~~

Annotate counts with gene information:

.. code-block:: bash

   wasp2-count count-variants rnaseq.bam genotypes.pgen \
     --samples NA12878 \
     --region genes.gtf \
     --gene_feature exon \
     --gene_attribute gene_id \
     --out_file gene_counts.tsv

ATAC-seq with Peaks
~~~~~~~~~~~~~~~~~~~

Annotate counts with ATAC-seq peaks:

.. code-block:: bash

   wasp2-count count-variants atac.bam variants.bcf \
     --samples NA12878 \
     --region peaks.narrowPeak \
     --out_file peak_counts.tsv

Output Format
-------------

Tab-separated file with the following columns:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Column
     - Description
   * - ``chr``
     - Chromosome name
   * - ``pos``
     - SNP position (1-based)
   * - ``ref``
     - Reference allele
   * - ``alt``
     - Alternate allele
   * - ``ref_count``
     - Reads supporting reference allele
   * - ``alt_count``
     - Reads supporting alternate allele
   * - ``other_count``
     - Reads with other alleles
   * - ``total_count``
     - Total overlapping reads
   * - ``region``
     - Overlapping region (if ``--region`` used)
   * - ``gene_id``
     - Gene ID (if GTF/GFF3 used)

Example output:

.. code-block:: text

   chr     pos       ref  alt  ref_count  alt_count  other_count  gene_id
   chr10   1000000   A    G    12         15         0            ENSG00000123456
   chr10   1001000   C    T    20         18         1            ENSG00000123456
   chr10   1050000   G    A    8          10         0            ENSG00000789012

Performance Tips
----------------

Use High-Performance Formats
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For large variant files (>10M variants):

1. **PGEN format** (fastest, ~25x speedup):

   .. code-block:: bash

      plink2 --vcf variants.vcf.gz --make-pgen --out variants
      wasp2-count count-variants sample.bam variants.pgen

2. **cyvcf2 backend** (7x speedup for VCF):

   .. code-block:: bash

      pip install wasp2[cyvcf2]
      wasp2-count count-variants sample.bam variants.vcf.gz

3. **BCF format** (5-8x speedup):

   .. code-block:: bash

      bcftools view -O b variants.vcf.gz > variants.bcf
      wasp2-count count-variants sample.bam variants.bcf

Process by Chromosome
~~~~~~~~~~~~~~~~~~~~~

For very large files, process chromosomes separately:

.. code-block:: bash

   for chr in {1..22} X Y; do
     wasp2-count count-variants sample.bam variants.pgen \
       --region chr${chr}.bed \
       --out_file counts_chr${chr}.tsv
   done

   # Combine results
   head -1 counts_chr1.tsv > all_counts.tsv
   tail -n +2 -q counts_chr*.tsv >> all_counts.tsv

Troubleshooting
---------------

No Output SNPs
~~~~~~~~~~~~~~

**Problem**: Output file is empty or has only header

**Diagnostic**:

.. code-block:: bash

   # Check for heterozygous SNPs
   bcftools view -s sample1 -g ^0/0,^1/1 variants.vcf.gz | head

   # Check BAM coverage
   samtools depth sample.bam | head

**Solutions**:

1. Verify sample name: ``bcftools query -l variants.vcf.gz``
2. Check chromosome naming (chr10 vs 10)
3. Ensure same reference genome for BAM and VCF

Low Count Numbers
~~~~~~~~~~~~~~~~~

**Problem**: Counts are unexpectedly low

**Diagnostic**:

.. code-block:: bash

   # Check read depth
   samtools depth sample.bam | awk '{sum+=$3; count++} END {print sum/count}'

   # Check mapping quality
   samtools flagstat sample.bam

**Solutions**:

1. Check sequencing depth (need >10x for reliable counts)
2. Verify BAM quality (remove duplicates, low-quality reads)
3. Ensure variants overlap sequenced regions

See Also
--------

* :doc:`/api/counting` - Python API documentation
* :doc:`/tutorials/rnaseq_ase` - RNA-seq tutorial
* :doc:`/tutorials/atacseq_ase` - ATAC-seq tutorial
* :doc:`wasp2_analyze` - Analyze allelic imbalance
```

### 4.4 Interactive Examples in Docstrings

Use doctest format for runnable examples:

```python
def parse_genotype(gt_string: str) -> Tuple[int, int]:
    """Parse VCF genotype string to allele indices.

    Args:
        gt_string: VCF format genotype (e.g., "0/1", "1|0", "./.")

    Returns:
        Tuple of (allele1, allele2) indices. Returns (-1, -1) for missing.

    Examples:
        >>> parse_genotype("0/1")
        (0, 1)

        >>> parse_genotype("1|0")
        (1, 0)

        >>> parse_genotype("./.")
        (-1, -1)

        >>> parse_genotype("1/1")
        (1, 1)

    Note:
        Phased (|) and unphased (/) genotypes are treated identically
        for allele extraction. Use separate functions if phasing matters.
    """
    if gt_string == "./." or gt_string == ".|.":
        return (-1, -1)

    separator = "|" if "|" in gt_string else "/"
    alleles = gt_string.split(separator)
    return (int(alleles[0]), int(alleles[1]))
```

---

## 5. Comparison with Successful Bioinformatics Tools

### 5.1 What WASP2 Can Learn From

#### STAR (RNA-seq aligner)
**Strengths**:
- Comprehensive manual (40+ pages PDF)
- Detailed parameter descriptions with biological context
- Performance benchmarks prominently displayed
- Example commands for every use case

**Apply to WASP2**:
- Create comprehensive PDF manual (in addition to web docs)
- Add biological context to parameter descriptions
- Expand benchmark section

#### salmon (RNA-seq quantification)
**Strengths**:
- Clear "Getting Started" tutorial
- Extensive FAQ section
- Algorithm explanation with diagrams
- Output format documentation with example data

**Apply to WASP2**:
- Add FAQ section (see 5.2 below)
- Create algorithm diagrams for WASP
- Expand output format documentation with examples

#### cellranger (10x Genomics single-cell)
**Strengths**:
- Use-case driven documentation structure
- Clear system requirements
- Troubleshooting decision trees
- Runtime and resource estimates

**Apply to WASP2**:
- Add runtime estimates for different data sizes
- Create troubleshooting decision trees
- Document system requirements more clearly

#### bcftools (Variant manipulation)
**Strengths**:
- Excellent man pages
- One-liner examples for common tasks
- Clear cheat sheets
- Integration examples with other tools

**Apply to WASP2**:
- Create man pages (section 3.2)
- Develop one-liner cheat sheet
- Add pipeline integration examples

### 5.2 FAQ Section Template

Create `docs/source/faq.rst`:

```rst
Frequently Asked Questions
==========================

General
-------

What is allelic imbalance?
~~~~~~~~~~~~~~~~~~~~~~~~~~

Allelic imbalance (AI) occurs when one allele of a heterozygous variant
is preferentially expressed or accessible compared to the other allele.
This can indicate:

* **cis-regulatory variants**: SNPs affecting gene regulation
* **Imprinting**: Parent-of-origin specific expression
* **X-inactivation**: Random silencing of one X chromosome
* **Technical artifacts**: Mapping bias, PCR bias

When should I use WASP2 vs GATK ASEReadCounter?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use **WASP2** if:

* You need reference bias correction (WASP mapping)
* Analyzing single-cell data
* Want statistical testing for allelic imbalance
* Need high performance (Rust acceleration)

Use **GATK ASEReadCounter** if:

* You only need raw allele counts
* Already using GATK workflows
* Don't need statistical analysis

Do I need to run WASP mapping before counting?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**It depends on your aligner and reference genome**:

* **Yes, use WASP** if you used standard aligners (STAR, BWA, bowtie2)
  and have divergent haplotypes
* **Maybe not needed** if you used allele-aware aligners or references
  (WASP-corrected STAR, diploid reference genome)

Rule of thumb: If in doubt, run WASP mapping. It's conservative and won't
hurt accuracy.

Installation
------------

Installation fails with "Rust compiler not found"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Install Rust using rustup
   curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
   source $HOME/.cargo/env

   # Retry WASP2 installation
   pip install wasp2

Can I install WASP2 without Rust?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Yes, but you'll miss significant performance benefits. WASP2 includes
Python fallbacks for all Rust-accelerated functions.

To disable Rust requirement:

.. code-block:: bash

   # Install without building Rust extension
   pip install wasp2 --no-build-isolation

   # Or set environment variable
   export WASP2_DISABLE_RUST=1

Performance will be 10-25x slower for counting and mapping operations.

Data Formats
------------

What variant formats does WASP2 support?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1

   * - Format
     - Extensions
     - Speed
     - Use Case
   * - VCF (pysam)
     - .vcf, .vcf.gz
     - Baseline (1x)
     - Default, compatibility
   * - VCF (cyvcf2)
     - .vcf, .vcf.gz
     - 7x faster
     - Production (install cyvcf2)
   * - BCF
     - .bcf
     - 5-8x faster
     - Binary VCF
   * - PGEN
     - .pgen
     - 25x faster
     - Large cohorts (install Pgenlib)

How do I convert VCF to PGEN?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Install plink2
   wget https://s3.amazonaws.com/plink2-assets/alpha3/plink2_linux_x86_64.zip
   unzip plink2_linux_x86_64.zip

   # Convert VCF to PGEN
   ./plink2 --vcf variants.vcf.gz --make-pgen --out variants

   # Use in WASP2
   wasp2-count count-variants sample.bam variants.pgen

Do BAM and VCF need to use the same reference genome?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Yes, absolutely**. Mismatched reference genomes will cause:

* Missing SNPs (different coordinates)
* Incorrect counts (different alleles)
* Chromosome naming issues (chr10 vs 10)

Verify your references:

.. code-block:: bash

   # Check BAM header
   samtools view -H sample.bam | grep "@SQ"

   # Check VCF header
   bcftools view -h variants.vcf.gz | grep "##contig"

   # Should match reference genome (e.g., both GRCh38)

Analysis
--------

How many reads do I need for allelic imbalance analysis?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Minimum recommendations**:

* **Per SNP**: ≥10 reads total (5 per allele)
* **Per gene/peak**: ≥20 reads total across all SNPs
* **For single-cell**: ≥100 cells per cell type

More reads = higher statistical power to detect imbalance.

What does "FDR < 0.05" mean in results?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

False Discovery Rate (FDR) is the expected proportion of false positives
among significant results.

* **FDR < 0.05**: Expect <5% of "significant" genes to be false positives
* **FDR < 0.01**: More stringent, <1% false positives

Use FDR instead of raw p-values when testing many genes/peaks.

Why are some genes significant with weak allelic imbalance?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

High coverage genes can show statistical significance even with small
allelic ratios (e.g., 55:45 instead of 50:50).

**Interpretation**:

* **Statistical significance** (FDR < 0.05): Effect is real, not random
* **Biological significance**: Depends on effect size and context

Filter by effect size for biologically relevant results:

.. code-block:: bash

   # Genes with strong imbalance (ratio >2:1)
   awk 'NR==1 || ($8 < 0.05 && ($5/$6 > 2 || $6/$5 > 2))' results.tsv

Single-Cell
-----------

How should I handle low coverage in single cells?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Strategies**:

1. **Aggregate by cell type**: Combine cells before analysis
2. **Lower threshold**: Use ``--min 5`` instead of default 10
3. **Filter features**: Only analyze high-coverage peaks/genes
4. **Pseudobulk**: Sum counts across cells of same type

Example aggregation:

.. code-block:: python

   import anndata as ad

   adata = ad.read_h5ad('sc_counts.h5ad')

   # Sum counts by cell type
   adata_bulk = adata.obs.groupby('celltype').sum()

Can I analyze multiple samples in single-cell data?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**It's complicated**. Single-cell barcodes are sample-specific, so
analyzing multiple samples requires:

1. **Demultiplexing**: Assign cells to samples (e.g., using genotypes)
2. **Sample-specific counting**: Run ``count-variants-sc`` per sample
3. **Combined analysis**: Merge h5ad objects with sample labels

For now, **analyze one sample at a time** and combine results downstream.

Troubleshooting
---------------

"Sample not found in VCF" error
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # List samples in VCF
   bcftools query -l variants.vcf.gz

   # Use exact sample name
   wasp2-count count-variants sample.bam variants.vcf.gz \
     --samples "SAMPLE_NAME_FROM_VCF"

"No space left on device" error
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

WASP2 creates temporary files during processing.

**Solutions**:

.. code-block:: bash

   # Use different temp directory
   wasp2-count count-variants sample.bam variants.vcf.gz \
     --temp_loc /scratch/large_disk/

   # Or clean up space
   df -h  # Check disk usage
   rm -rf /tmp/*  # Clear system temp (carefully!)

"TypeError: 'NoneType' object is not subscriptable"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This usually means a required file is missing or empty.

**Diagnostic**:

.. code-block:: bash

   # Check all files exist and are non-empty
   ls -lh sample.bam sample.bam.bai variants.vcf.gz variants.vcf.gz.tbi

   # Check VCF has data
   bcftools view variants.vcf.gz | head

**Common causes**:

* Missing BAM index (.bai)
* Missing VCF index (.tbi)
* Empty VCF file
* Corrupt BAM file

See :doc:`tutorials/troubleshooting` for more debugging tips.
```

---

## 6. Implementation Priority

### Phase 1: Quick Wins (1-2 weeks)
1. Enhanced README with badges, quick start, citation, comparison table
2. Basic FAQ section
3. Shell completion scripts
4. Example commands directory

### Phase 2: Core Documentation (2-3 weeks)
1. Tutorial series (concepts through troubleshooting)
2. Enhanced --help output (better examples and descriptions)
3. CLI reference documentation in Sphinx
4. Performance tuning guide

### Phase 3: Advanced Documentation (2-3 weeks)
1. Man pages for all commands
2. Comprehensive API docstrings (Google style)
3. Jupyter notebook examples
4. Integration guides (Nextflow, Snakemake, CWL)

### Phase 4: Polish (1 week)
1. Diagrams and illustrations
2. Video tutorials (optional)
3. Interactive documentation features
4. Translation (optional, if international audience)

---

## 7. Maintenance and Versioning

### Documentation Versioning
Use Read the Docs or GitHub Pages with version switcher:

```yaml
# .readthedocs.yml
version: 2

build:
  os: ubuntu-22.04
  tools:
    python: "3.10"

sphinx:
  configuration: docs/source/conf.py

python:
  install:
    - requirements: docs/requirements.txt

versions:
  - latest
  - stable
  - v1.2
  - v1.1
```

### Documentation Testing
```bash
# Test docstrings
python -m doctest counting/run_counting.py

# Test Sphinx build
cd docs && make clean && make html

# Check for broken links
sphinx-build -b linkcheck source build/linkcheck

# Spell check
sphinx-build -b spelling source build/spelling
```

### Documentation Metrics
Track documentation quality:
- Coverage: % of functions with docstrings
- Broken links: Regular link checking
- User feedback: GitHub issues tagged "documentation"
- Search analytics: Most searched terms (add Google Analytics)

---

## 8. Resources and References

### Style Guides
- **Google Python Style Guide**: https://google.github.io/styleguide/pyguide.html
- **NumPy Docstring Guide**: https://numpydoc.readthedocs.io/
- **Divio Documentation System**: https://documentation.divio.com/

### Tools
- **Sphinx**: https://www.sphinx-doc.org/
- **Read the Docs**: https://readthedocs.org/
- **MkDocs**: https://www.mkdocs.org/ (alternative to Sphinx)
- **Typer**: https://typer.tiangolo.com/

### Examples of Excellent Bioinformatics Documentation
- **STAR**: https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
- **salmon**: https://salmon.readthedocs.io/
- **scanpy**: https://scanpy.readthedocs.io/
- **snakemake**: https://snakemake.readthedocs.io/
- **bcftools**: http://samtools.github.io/bcftools/

---

## Summary

This plan provides a comprehensive roadmap for elevating WASP2's documentation to production-grade standards. Key recommendations:

1. **README**: Add badges, quick start, citation, comparison table, and learning paths
2. **Tutorials**: Create progressive tutorial series from 5-min quickstart to advanced workflows
3. **CLI**: Enhance --help output, create man pages, provide shell completion
4. **API**: Use Google-style docstrings with comprehensive examples and type hints
5. **Structure**: Organize docs using Divio framework (tutorials, how-to, reference, explanation)

The documentation should serve users at all levels, from newcomers exploring allele-specific analysis to power users optimizing large-scale pipelines.

Implementation can be phased over 6-8 weeks, with quick wins (README, FAQ, examples) delivering immediate value while larger efforts (full tutorial series, man pages) provide long-term benefits.
