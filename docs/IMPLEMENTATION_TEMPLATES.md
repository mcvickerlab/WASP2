# WASP2 Documentation Implementation Templates

Quick reference templates for implementing the documentation plan.

## Table of Contents
1. [README Templates](#readme-templates)
2. [Tutorial Templates](#tutorial-templates)
3. [Docstring Templates](#docstring-templates)
4. [CLI Help Templates](#cli-help-templates)
5. [Sphinx Configuration](#sphinx-configuration)

---

## README Templates

### Badge Section (Enhanced)

```markdown
<p align="center">
  <!-- Build & Testing -->
  <a href="https://github.com/Jaureguy760/WASP2-exp/actions/workflows/ci.yml">
    <img src="https://github.com/Jaureguy760/WASP2-exp/actions/workflows/ci.yml/badge.svg" alt="CI">
  </a>
  <a href="https://codecov.io/gh/Jaureguy760/WASP2-exp">
    <img src="https://codecov.io/gh/Jaureguy760/WASP2-exp/branch/main/graph/badge.svg" alt="Coverage">
  </a>

  <!-- Documentation -->
  <a href="https://jaureguy760.github.io/WASP2-exp/">
    <img src="https://img.shields.io/badge/docs-latest-blue" alt="Documentation">
  </a>
  <a href="https://github.com/Jaureguy760/WASP2-exp/actions/workflows/docs.yml">
    <img src="https://github.com/Jaureguy760/WASP2-exp/actions/workflows/docs.yml/badge.svg" alt="Docs Build">
  </a>

  <!-- Package Distribution -->
  <a href="https://pypi.org/project/wasp2/">
    <img src="https://img.shields.io/pypi/v/wasp2" alt="PyPI">
  </a>
  <a href="https://anaconda.org/bioconda/wasp2">
    <img src="https://img.shields.io/conda/vn/bioconda/wasp2" alt="Bioconda">
  </a>
  <img src="https://img.shields.io/pypi/dm/wasp2" alt="Downloads">

  <!-- Language & License -->
  <a href="https://github.com/Jaureguy760/WASP2-exp/blob/master/LICENSE">
    <img src="https://img.shields.io/badge/license-MIT-green" alt="License">
  </a>
  <img src="https://img.shields.io/badge/python-3.10+-blue" alt="Python">
  <img src="https://img.shields.io/badge/rust-1.70+-orange" alt="Rust">

  <!-- Community -->
  <img src="https://img.shields.io/github/stars/Jaureguy760/WASP2-exp?style=social" alt="Stars">
  <a href="https://github.com/Jaureguy760/WASP2-exp/issues">
    <img src="https://img.shields.io/github/issues/Jaureguy760/WASP2-exp" alt="Issues">
  </a>
</p>
```

### Quick Start Section

```markdown
## Quick Start

Get started with WASP2 in under 5 minutes:

```bash
# 1. Install WASP2
pip install wasp2

# 2. Count allele-specific reads
wasp2-count count-variants \
  sample.bam \
  variants.vcf.gz \
  --samples NA12878 \
  --out_file counts.tsv

# 3. Detect allelic imbalance
wasp2-analyze find-imbalance \
  counts.tsv \
  --out_file results.tsv

# 4. View significant results (FDR < 0.05)
awk 'NR==1 || $8 < 0.05' results.tsv | column -t | head -20
```

**What you get**: Statistical tests showing which genes/regions have significant allelic imbalance.

**Next steps**:
- [Full Tutorial](docs/tutorials/basic_workflow.md) - 30-minute walkthrough
- [RNA-seq Guide](docs/tutorials/rnaseq_ase.md) - RNA-seq specific workflow
- [Documentation](https://jaureguy760.github.io/WASP2-exp/) - Complete reference
```

### Installation Options Matrix

```markdown
## Installation

Choose the installation method that fits your needs:

| Method | Use Case | Installation Time | Command |
|--------|----------|------------------|---------|
| **PyPI** | Most users | ~1 minute | `pip install wasp2` |
| **PyPI + Performance** | Production | ~2 minutes | `pip install wasp2[cyvcf2,plink]` |
| **Conda** | Conda users | ~5 minutes | `conda install -c bioconda wasp2` |
| **From Source** | Developers | ~10 minutes | See below |
| **GitHub Codespaces** | Try without installing | ~3 minutes | Click "Code" → "Codespaces" |

### Standard Installation

```bash
pip install wasp2
```

### With Performance Enhancements

```bash
# Install with cyvcf2 (7x faster VCF parsing)
pip install wasp2[cyvcf2]

# Install with PLINK2 support (25x faster variant I/O)
pip install wasp2[plink]

# Install everything
pip install wasp2[all]
```

### Developer Installation

```bash
git clone https://github.com/Jaureguy760/WASP2-exp.git
cd WASP2-exp

# Create environment
conda env create -f environment.yml
conda activate WASP2

# Build Rust extension
export LIBCLANG_PATH=$CONDA_PREFIX/lib
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
maturin develop --release -m rust/Cargo.toml

# Install in development mode
pip install -e ".[dev,docs]"
```
```

### Citation Section

```markdown
## Citation

If you use WASP2 in published research, please cite:

**WASP2 paper** (when available):
```bibtex
@article{wasp2_2025,
  title={WASP2: High-performance allele-specific analysis with Rust acceleration},
  author={Ho, Aaron and Jaureguy, Jeff and McVicker, Graham},
  journal={Bioinformatics},
  year={2025},
  note={In preparation}
}
```

**Original WASP algorithm**:
```bibtex
@article{vandegeijn2015wasp,
  title={{WASP}: allele-specific software for robust molecular quantitative trait locus discovery},
  author={van de Geijn, Bryce and McVicker, Graham and Gilad, Yoav and Pritchard, Jonathan K},
  journal={Nature Methods},
  volume={12},
  number={11},
  pages={1061--1063},
  year={2015},
  publisher={Nature Publishing Group},
  doi={10.1038/nmeth.3582}
}
```

### Key Publications

WASP2 builds on and extends these methods:

- **Reference bias correction**: van de Geijn et al. (2015) *Nature Methods*
- **Beta-binomial testing**: Skelly et al. (2011) *Genome Research*
- **Single-cell ASE**: Larsson et al. (2019) *Nature Communications*
```

---

## Tutorial Templates

### Tutorial Front Matter Template

```markdown
# [Tutorial Title]

**Estimated Time**: XX minutes
**Difficulty**: [Beginner | Intermediate | Advanced]
**Prerequisites**:
- Prerequisite 1
- Prerequisite 2

**Dataset**:
- Description of dataset
- Download link or instructions

**Learning Objectives**:

By completing this tutorial, you will learn how to:
- [ ] Learning objective 1
- [ ] Learning objective 2
- [ ] Learning objective 3

---

## Table of Contents

1. [Background](#background)
2. [Setup](#setup)
3. [Step 1: ...](#step-1-...)
4. [Step 2: ...](#step-2-...)
5. [Interpreting Results](#interpreting-results)
6. [Troubleshooting](#troubleshooting)
7. [Next Steps](#next-steps)

---

## Background

[2-3 paragraphs explaining the biological/technical context]

---

## Setup

### Download Data

```bash
# Download example dataset
wget https://example.com/tutorial_data.tar.gz
tar -xzf tutorial_data.tar.gz
cd tutorial_data/

# Verify contents
ls -lh
```

### Expected Files

```
tutorial_data/
├── sample.bam          # Aligned reads (500 MB)
├── sample.bam.bai      # BAM index
├── variants.vcf.gz     # Genotypes (100 MB)
├── variants.vcf.gz.tbi # VCF index
└── regions.bed         # Genomic regions (1 MB)
```

---

## Step 1: [Action Verb - e.g., "Count Alleles"]

### Goal

[What you'll accomplish in this step]

### Command

```bash
wasp2-count count-variants \
  sample.bam \
  variants.vcf.gz \
  --samples NA12878 \
  --region regions.bed \
  --out_file counts.tsv
```

### Explanation

- `sample.bam` - Input aligned reads
- `variants.vcf.gz` - Genotype information for NA12878
- `--samples NA12878` - Filter to heterozygous SNPs in this sample
- `--region regions.bed` - Only count SNPs in these regions
- `--out_file counts.tsv` - Save results here

### Expected Output

```
Processing variants...
Found 10,523 heterozygous SNPs for NA12878
Overlapping 2,341 genomic regions
Counting alleles...
Processed 1,000,000 reads
Output written to counts.tsv
```

### Verification

```bash
# Check output file
head -5 counts.tsv

# Count total SNPs
wc -l counts.tsv  # Should be ~2,342 (header + 2,341 SNPs)

# Check for reasonable coverage
awk 'NR>1 {print $5+$6}' counts.tsv | \
  awk '{sum+=$1; count++} END {print "Average coverage:", sum/count}'
```

### Expected Results

- File: `counts.tsv` (approximately XXX KB)
- Total SNPs: ~2,341
- Average coverage: ~30-50 reads per SNP

---

[Repeat for each step...]

---

## Interpreting Results

### Output Format

The `results.tsv` file contains:

| Column | Description | Example Value |
|--------|-------------|---------------|
| `region` | Genomic region | chr10:1000000-1001000 |
| `n_snps` | Number of SNPs | 3 |
| `ref_total` | Total reference reads | 45 |
| `alt_total` | Total alternate reads | 55 |
| `p_value` | Statistical p-value | 0.023 |
| `fdr` | FDR-adjusted p-value | 0.045 |
| `log2_ratio` | log2(alt/ref) | 0.29 |

### What to Look For

**Significant allelic imbalance** (FDR < 0.05):
- These regions show non-random allele expression
- May indicate cis-regulatory variants
- Requires follow-up validation

**High log2_ratio** (|ratio| > 1):
- One allele >2x more expressed than other
- Strong biological effect
- Prime candidates for functional studies

**Low p-value but high FDR**:
- Not statistically significant after multiple testing correction
- May be interesting but require larger sample size

### Quality Control

```bash
# Distribution of p-values (should be uniform under null hypothesis)
awk 'NR>1 {print $5}' results.tsv | \
  sort -n | \
  awk '{print int($1*10)/10}' | \
  uniq -c

# Coverage distribution
awk 'NR>1 {print $3+$4}' results.tsv | \
  awk '{if($1<10) low++; else if($1<50) med++; else high++}
       END {print "Low (<10):", low, "Medium (10-50):", med, "High (>50):", high}'
```

---

## Troubleshooting

### Problem: No output file generated

**Diagnostic**:
```bash
# Check for error messages
echo $?  # Should be 0 for success

# Check disk space
df -h .
```

**Possible Causes**:
1. Insufficient disk space
2. Permission error
3. Invalid input files

**Solutions**:
```bash
# Free up space or change output location
wasp2-count count-variants sample.bam variants.vcf.gz \
  --temp_loc /scratch/temp/ \
  --out_file /scratch/results/counts.tsv

# Check file permissions
ls -l sample.bam variants.vcf.gz
```

---

### Problem: Very few SNPs in output

**Diagnostic**:
```bash
# Check number of het SNPs for sample
bcftools view -s NA12878 -g het variants.vcf.gz | grep -v "^#" | wc -l

# Check BAM coverage
samtools depth sample.bam | awk '{sum+=$3; n++} END {print "Mean depth:", sum/n}'
```

**Possible Causes**:
1. Wrong sample name
2. Low sequencing coverage
3. Chromosome naming mismatch (chr10 vs 10)

**Solutions**:
```bash
# List available samples
bcftools query -l variants.vcf.gz

# Check chromosome naming
samtools view -H sample.bam | grep "^@SQ" | head -3
bcftools view -h variants.vcf.gz | grep "^##contig" | head -3

# Fix if needed (rename chromosomes in VCF)
bcftools annotate --rename-chrs chr_name_conv.txt variants.vcf.gz -Oz -o fixed.vcf.gz
```

---

## Next Steps

Now that you've completed this tutorial:

1. **Try with your own data**: Adapt these commands to your dataset
2. **Explore other workflows**:
   - [ATAC-seq Analysis](atac_ase.md)
   - [Single-Cell Workflow](single_cell.md)
3. **Learn advanced features**:
   - [Performance Tuning](../how_to/optimize_performance.md)
   - [Pipeline Integration](../how_to/integrate_with_pipelines.md)
4. **Understand the methods**:
   - [WASP Algorithm](../explanations/wasp_algorithm.md)
   - [Statistical Models](../explanations/statistical_models.md)

---

## Further Reading

- Original WASP paper: van de Geijn et al. (2015) *Nature Methods*
- Beta-binomial models: Skelly et al. (2011) *Genome Research*
- WASP2 API documentation: [Counting Module](../../api/counting.rst)

---

## Feedback

Found an issue with this tutorial? Please [open an issue](https://github.com/Jaureguy760/WASP2-exp/issues/new) or suggest improvements.
```

---

## Docstring Templates

### Function Docstring (Google Style)

```python
def run_count_variants(
    bam_file: Union[str, Path],
    variant_file: Union[str, Path],
    region_file: Optional[Union[str, Path]] = None,
    samples: Optional[str] = None,
    out_file: Optional[Union[str, Path]] = None,
    min_mapping_quality: int = 10,
    min_base_quality: int = 20,
    use_rust: bool = True,
    threads: int = 1,
) -> None:
    """Count allele-specific reads at heterozygous SNP positions.

    Quantifies reads supporting reference vs. alternate alleles at heterozygous
    single nucleotide polymorphisms (SNPs). This is the first step in allelic
    imbalance analysis, producing per-SNP allele counts for downstream statistical
    testing.

    The function processes aligned reads from a BAM file and variant calls from
    a VCF/BCF/PGEN file. It can filter variants by sample genotype and annotate
    counts with genomic regions (genes, ATAC-seq peaks, etc.).

    Args:
        bam_file: Path to aligned reads in BAM format. Must be coordinate-sorted
            and indexed (.bai file required in same directory).
        variant_file: Path to variant calls. Supports VCF (.vcf, .vcf.gz),
            BCF (.bcf), and PLINK2 PGEN (.pgen) formats. For VCF/BCF, index
            files (.tbi or .csi) are recommended for faster processing.
        region_file: Path to genomic regions for SNP filtering. Accepts BED,
            GTF, GFF3, or narrowPeak formats. If provided, only SNPs overlapping
            these regions are counted. Default: None (use all SNPs).
        samples: Sample ID(s) to filter heterozygous SNPs. Accepts comma-separated
            IDs (e.g., "sample1,sample2") or path to file with one ID per line.
            If None, all variants are used regardless of genotype. Default: None.
        out_file: Output file path for allele counts (TSV format). If None,
            defaults to "counts.tsv" in current directory. Default: None.
        min_mapping_quality: Minimum mapping quality (MAPQ) for reads to be
            counted. Reads with MAPQ below this threshold are ignored. Typical
            values: 10 (permissive), 20 (moderate), 30 (strict). Default: 10.
        min_base_quality: Minimum base quality (Phred score) at SNP position
            for read to be counted. Bases below this quality are ignored.
            Typical values: 20 (moderate), 30 (strict). Default: 20.
        use_rust: If True, use Rust-accelerated counting implementation (requires
            wasp2_rust extension). Falls back to Python if extension unavailable.
            Rust implementation is ~10-25x faster. Default: True.
        threads: Number of threads for BAM I/O operations. Currently only
            supported by Rust implementation. Default: 1.

    Returns:
        None. Results are written to out_file.

    Raises:
        FileNotFoundError: If bam_file, variant_file, or region_file does not exist.
        ValueError: If sample ID not found in variant file, or if variant_file
            format cannot be determined from extension.
        RuntimeError: If BAM file is not sorted or indexed, or if Rust extension
            fails unexpectedly.
        IOError: If output file cannot be written (permission denied, disk full).
        MemoryError: If system runs out of memory (try processing by chromosome).

    Examples:
        Basic counting at all variants:

        >>> run_count_variants(
        ...     bam_file="sample.bam",
        ...     variant_file="variants.vcf.gz",
        ...     out_file="counts.tsv"
        ... )

        Filter by sample and annotate with genes:

        >>> run_count_variants(
        ...     bam_file="rnaseq.bam",
        ...     variant_file="genotypes.pgen",
        ...     region_file="genes.gtf",
        ...     samples="NA12878",
        ...     out_file="gene_counts.tsv"
        ... )

        ATAC-seq with peak annotation:

        >>> run_count_variants(
        ...     bam_file="atac.bam",
        ...     variant_file="variants.bcf",
        ...     region_file="peaks.narrowPeak",
        ...     samples="NA12878",
        ...     min_mapping_quality=30,
        ...     out_file="peak_counts.tsv"
        ... )

        Process multiple samples:

        >>> run_count_variants(
        ...     bam_file="multi_sample.bam",
        ...     variant_file="1000G.vcf.gz",
        ...     samples="NA12878,NA12891,NA12892",
        ...     out_file="multi_counts.tsv"
        ... )

    Notes:
        **Output Format:**
        Tab-separated file with columns:

        - chr: Chromosome name
        - pos: SNP position (1-based)
        - ref: Reference allele
        - alt: Alternate allele
        - ref_count: Reads supporting reference allele
        - alt_count: Reads supporting alternate allele
        - other_count: Reads with other alleles
        - total_count: Total overlapping reads
        - region: Overlapping region (if region_file provided)

        **Performance Tips:**

        - Use PGEN format for large variant files (>10M variants, ~25x speedup)
        - Install cyvcf2 for faster VCF parsing: ``pip install wasp2[cyvcf2]``
        - Process chromosomes separately for very large datasets
        - Use ``threads > 1`` with Rust implementation for faster I/O

        **Memory Considerations:**

        - Typical memory usage: 2-8 GB for whole-genome data
        - PGEN format uses less memory than VCF
        - Process by chromosome if encountering memory issues

        **Quality Control:**

        - Check BAM alignment rate: ``samtools flagstat sample.bam``
        - Verify sample names: ``bcftools query -l variants.vcf.gz``
        - Ensure matching reference genomes (BAM and VCF)
        - Check chromosome naming consistency (chr10 vs 10)

    See Also:
        run_ai_analysis: Detect allelic imbalance from count data.
        run_make_remap_reads: Generate reads for WASP mapping.
        count_variants_sc: Count alleles in single-cell data.

    References:
        van de Geijn, B., McVicker, G., Gilad, Y., & Pritchard, J. K. (2015).
        WASP: allele-specific software for robust molecular quantitative trait
        locus discovery. Nature Methods, 12(11), 1061-1063.
        https://doi.org/10.1038/nmeth.3582

    Version History:
        - v1.0.0: Initial Python implementation
        - v1.1.0: Added PGEN format support
        - v1.2.0: Rust acceleration, cyvcf2 support
        - v1.2.1: Multi-threading support in Rust
    """
    # Implementation
    pass
```

### Class Docstring Template

```python
@dataclass
class WaspCountFiles:
    """Container for WASP counting workflow files and metadata.

    Manages file paths and temporary directories for the counting workflow.
    Handles cleanup of temporary files on context exit.

    This class is typically used as a context manager to ensure proper cleanup
    of temporary files, even if an exception occurs during processing.

    Attributes:
        bam_file: Path to input BAM file
        variant_file: Path to variant file (VCF/BCF/PGEN)
        region_file: Path to region file (BED/GTF), or None
        out_file: Path to output counts file
        temp_dir: Temporary directory for intermediate files
        vcf_bed: Path to converted VCF BED file
        intersect_bed: Path to intersected BED file
        keep_temp: If True, preserve temporary files after completion

    Examples:
        Basic usage with automatic cleanup:

        >>> with WaspCountFiles(
        ...     bam_file="sample.bam",
        ...     variant_file="variants.vcf.gz",
        ...     out_file="counts.tsv"
        ... ) as files:
        ...     # Process files
        ...     process_counts(files)
        ...     # Temp files automatically cleaned up here

        Preserve temporary files for debugging:

        >>> files = WaspCountFiles(
        ...     bam_file="sample.bam",
        ...     variant_file="variants.vcf.gz",
        ...     temp_loc="/scratch/debug/",
        ...     keep_temp=True
        ... )
        >>> # Temp files preserved in /scratch/debug/

    Notes:
        - Temporary directory is created lazily on first access
        - Context manager ensures cleanup even on exceptions
        - Set keep_temp=True or specify temp_loc to preserve intermediates
        - Intermediate files can be large (similar size to input VCF)

    See Also:
        run_count_variants: Main counting workflow using this class
    """
    bam_file: Path
    variant_file: Path
    region_file: Optional[Path] = None
    out_file: Path = Path("counts.tsv")
    temp_dir: Optional[Path] = None
    vcf_bed: Optional[Path] = None
    intersect_bed: Optional[Path] = None
    keep_temp: bool = False

    def __enter__(self) -> "WaspCountFiles":
        """Set up temporary directory on context entry."""
        if self.temp_dir is None:
            self.temp_dir = Path(tempfile.mkdtemp(prefix="wasp2_"))
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Clean up temporary files on context exit."""
        if not self.keep_temp and self.temp_dir:
            shutil.rmtree(self.temp_dir, ignore_errors=True)
```

### Module Docstring Template

```python
"""Allele-specific read counting module.

This module provides functions to count reads supporting reference vs. alternate
alleles at heterozygous SNP positions. It is the first step in allelic imbalance
analysis.

The main entry point is :func:`run_count_variants`, which orchestrates the
workflow:

1. Convert variant file to BED format (:func:`vcf_to_bed`)
2. Intersect variants with genomic regions (:func:`intersect_vcf_region`)
3. Count alleles at each SNP (:func:`make_count_df`)
4. Write results to output file

Typical Usage
-------------

Basic counting::

    from counting.run_counting import run_count_variants

    run_count_variants(
        bam_file="sample.bam",
        variant_file="variants.vcf.gz",
        samples="NA12878",
        out_file="counts.tsv"
    )

With region annotation::

    run_count_variants(
        bam_file="rnaseq.bam",
        variant_file="genotypes.pgen",
        region_file="genes.gtf",
        samples="NA12878",
        out_file="gene_counts.tsv"
    )

Performance Optimization
------------------------

For large datasets:

1. **Use PGEN format** for 25x faster variant I/O::

       plink2 --vcf variants.vcf.gz --make-pgen --out variants
       run_count_variants(bam_file="sample.bam", variant_file="variants.pgen")

2. **Install cyvcf2** for 7x faster VCF parsing::

       pip install wasp2[cyvcf2]

3. **Process by chromosome** for very large files::

       for chrom in ['chr1', 'chr2', ...]:
           run_count_variants(
               bam_file="sample.bam",
               variant_file="variants.pgen",
               region_file=f"{chrom}.bed",
               out_file=f"counts_{chrom}.tsv"
           )

Module Contents
---------------

Main Functions
~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   run_count_variants
   run_count_variants_sc

Workflow Functions
~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   vcf_to_bed
   intersect_vcf_region
   make_count_df

Data Classes
~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   WaspCountFiles

See Also
--------
analysis.run_analysis : Statistical testing for allelic imbalance
mapping.run_mapping : WASP reference bias correction

References
----------
.. [1] van de Geijn et al. (2015). WASP: allele-specific software for robust
       molecular quantitative trait locus discovery. Nature Methods 12:1061-1063.

Examples
--------
Complete RNA-seq workflow:

>>> # Step 1: Count alleles
>>> from counting.run_counting import run_count_variants
>>> run_count_variants(
...     bam_file="rnaseq.bam",
...     variant_file="genotypes.pgen",
...     region_file="genes.gtf",
...     samples="NA12878",
...     out_file="gene_counts.tsv"
... )

>>> # Step 2: Analyze for allelic imbalance
>>> from analysis.run_analysis import run_ai_analysis
>>> run_ai_analysis(
...     count_file="gene_counts.tsv",
...     min_count=10,
...     out_file="gene_imbalance.tsv"
... )
"""

from .run_counting import run_count_variants
from .run_counting_sc import run_count_variants_sc
from .filter_variant_data import vcf_to_bed, intersect_vcf_region
from .count_alleles import make_count_df

__all__ = [
    "run_count_variants",
    "run_count_variants_sc",
    "vcf_to_bed",
    "intersect_vcf_region",
    "make_count_df",
]
```

---

## CLI Help Templates

### Enhanced Command Help (Typer)

```python
@app.command(
    help="""
    Count allele-specific reads at heterozygous SNP positions.

    Quantifies reads supporting reference vs. alternate alleles at heterozygous
    SNPs. This is the first step in allelic imbalance analysis.

    \b
    Quick Examples:
      # Basic counting
      wasp2-count count-variants sample.bam variants.vcf.gz

      # With sample filtering
      wasp2-count count-variants sample.bam variants.vcf.gz \\
        --samples NA12878 --out_file counts.tsv

      # RNA-seq with gene annotation
      wasp2-count count-variants rnaseq.bam genotypes.pgen \\
        --samples NA12878 --region genes.gtf --out_file gene_counts.tsv

    \b
    Output Format:
      Tab-separated file with columns:
        chr, pos, ref, alt - Variant information
        ref_count, alt_count - Reads per allele
        other_count - Reads with other alleles
        region - Overlapping region (if --region used)

    \b
    Performance Tips:
      - Use PGEN format for 25x faster I/O on large files
      - Install cyvcf2: pip install wasp2[cyvcf2] (7x VCF speedup)
      - Process by chromosome for very large datasets

    See full documentation at:
    https://jaureguy760.github.io/WASP2-exp/cli/wasp2_count.html
    """
)
def count_variants(
    bam: Annotated[
        str,
        typer.Argument(
            help="Aligned reads (BAM format, sorted and indexed)",
            metavar="BAM",
            show_default=False
        )
    ],
    variants: Annotated[
        str,
        typer.Argument(
            help="Variant calls (VCF, BCF, or PGEN format)",
            metavar="VARIANTS",
            show_default=False
        )
    ],
    samples: Annotated[
        Optional[List[str]],
        typer.Option(
            "--samples", "-s",
            help="Sample ID(s) for filtering heterozygous SNPs. "
                 "Comma-separated or file with one per line.",
            metavar="SAMPLE",
            show_default="all variants"
        )
    ] = None,
    region: Annotated[
        Optional[str],
        typer.Option(
            "--region", "-r",
            help="Genomic regions (BED, GTF, GFF3, narrowPeak). "
                 "Only count SNPs overlapping these regions.",
            metavar="PATH",
            show_default="all SNPs"
        )
    ] = None,
    out_file: Annotated[
        Optional[str],
        typer.Option(
            "--out_file", "-o",
            help="Output file path (TSV format)",
            metavar="PATH",
            show_default="counts.tsv"
        )
    ] = None,
    min_mapq: Annotated[
        int,
        typer.Option(
            "--min-mapq",
            help="Minimum mapping quality (MAPQ) for reads",
            metavar="INT",
            min=0,
            max=60,
            show_default=True
        )
    ] = 10,
    min_baseq: Annotated[
        int,
        typer.Option(
            "--min-baseq",
            help="Minimum base quality at SNP position",
            metavar="INT",
            min=0,
            max=60,
            show_default=True
        )
    ] = 20,
    use_rust: Annotated[
        bool,
        typer.Option(
            "--use-rust/--no-rust",
            help="Use Rust acceleration (10-25x faster)",
            show_default="--use-rust"
        )
    ] = True,
) -> None:
    """Count alleles at heterozygous SNPs."""

    # Parse samples
    sample_str = samples[0] if samples and len(samples) > 0 else None

    # Run counting
    run_count_variants(
        bam_file=bam,
        variant_file=variants,
        region_file=region,
        samples=sample_str,
        out_file=out_file,
        min_mapping_quality=min_mapq,
        min_base_quality=min_baseq,
        use_rust=use_rust,
    )
```

### Command Group Help

```python
app = typer.Typer(
    name="wasp2-count",
    help="""
    WASP2 Counting Module - Quantify allele-specific reads.

    This module counts reads supporting reference vs. alternate alleles at
    heterozygous SNP positions. It provides two commands:

      count-variants     Count alleles in bulk sequencing data
      count-variants-sc  Count alleles in single-cell data

    \b
    Quick Start:
      wasp2-count count-variants sample.bam variants.vcf.gz

    \b
    Common Workflows:
      RNA-seq ASE:
        wasp2-count count-variants rnaseq.bam genotypes.pgen \\
          --samples NA12878 --region genes.gtf --out_file gene_counts.tsv

      ATAC-seq:
        wasp2-count count-variants atac.bam variants.bcf \\
          --samples NA12878 --region peaks.narrowPeak --out_file peak_counts.tsv

      Single-cell:
        wasp2-count count-variants-sc sc.bam variants.pgen barcodes.txt \\
          --samples donor1 --out_file sc_counts.h5ad

    For detailed help on each command:
      wasp2-count count-variants --help
      wasp2-count count-variants-sc --help

    Full documentation: https://jaureguy760.github.io/WASP2-exp/
    """,
    no_args_is_help=True,
    add_completion=True,
)
```

---

## Sphinx Configuration

### Enhanced conf.py Additions

```python
# -- Project information (update version dynamically) -------------------------

import sys
from pathlib import Path

# Get version from pyproject.toml
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

try:
    from importlib.metadata import version
    release = version("wasp2")
except Exception:
    release = "1.2.1"  # Fallback

version = ".".join(release.split(".")[:2])  # Short version (1.2)

# -- General configuration (enhanced) ------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.autosummary",
    "sphinx.ext.coverage",
    "sphinx.ext.todo",
    "sphinx.ext.mathjax",  # For equations
    "sphinx.ext.graphviz",  # For diagrams
    "sphinx_copybutton",  # Copy code blocks
    "sphinx_tabs.tabs",  # Tabbed content
    "sphinx_design",  # Cards, grids, etc.
    "myst_parser",  # Markdown support
]

# MyST (Markdown) configuration
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "fieldlist",
    "html_admonition",
    "html_image",
    "linkify",
    "replacements",
    "smartquotes",
    "substitution",
    "tasklist",
]

# Autodoc configuration (enhanced)
autodoc_default_options = {
    "members": True,
    "member-order": "bysource",
    "special-members": "__init__,__call__",
    "undoc-members": True,
    "exclude-members": "__weakref__,__dict__,__module__",
    "show-inheritance": True,
    "inherited-members": False,
}

# Autosummary configuration
autosummary_generate = True
autosummary_imported_members = False

# Napoleon configuration (enhanced for better formatting)
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = True
napoleon_use_admonition_for_notes = True
napoleon_use_admonition_for_references = True
napoleon_use_ivar = True
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_use_keyword = True
napoleon_custom_sections = [
    ("Performance", "params_style"),
    ("Version History", "notes_style"),
]

# Intersphinx mapping (extended)
intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
    "scanpy": ("https://scanpy.readthedocs.io/en/stable/", None),
    "anndata": ("https://anndata.readthedocs.io/en/latest/", None),
}

# -- Options for HTML output (pydata theme enhanced) ---------------------------

html_theme = "pydata_sphinx_theme"

html_theme_options = {
    "github_url": "https://github.com/Jaureguy760/WASP2-exp",
    "use_edit_page_button": True,
    "show_toc_level": 2,
    "navbar_align": "left",
    "navbar_end": ["search-field", "navbar-icon-links"],
    "footer_items": ["copyright", "sphinx-version"],

    # Navigation
    "navigation_depth": 4,
    "collapse_navigation": False,
    "show_nav_level": 2,

    # Icons
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/Jaureguy760/WASP2-exp",
            "icon": "fa-brands fa-github",
            "type": "fontawesome",
        },
        {
            "name": "PyPI",
            "url": "https://pypi.org/project/wasp2/",
            "icon": "fa-solid fa-box",
            "type": "fontawesome",
        },
    ],

    # Announcement banner
    "announcement": "WASP2 v1.2.1 with Rust acceleration now available! 🚀",

    # External links
    "external_links": [
        {"name": "Tutorials", "url": "https://jaureguy760.github.io/WASP2-exp/tutorials/"},
        {"name": "Examples", "url": "https://github.com/Jaureguy760/WASP2-exp/tree/main/examples"},
    ],
}

html_context = {
    "github_user": "Jaureguy760",
    "github_repo": "WASP2-exp",
    "github_version": "main",
    "doc_path": "docs/source",
}

# Sidebars
html_sidebars = {
    "**": ["search-field", "sidebar-nav-bs", "sidebar-ethical-ads"],
}

# -- Copy button configuration --------------------------------------------------

copybutton_prompt_text = r">>> |\.\.\. |\$ "
copybutton_prompt_is_regexp = True
copybutton_only_copy_prompt_lines = True
copybutton_remove_prompts = True

# -- Code highlighting ----------------------------------------------------------

pygments_style = "sphinx"
pygments_dark_style = "monokai"

# -- LaTeX configuration (for PDF generation) -----------------------------------

latex_elements = {
    "papersize": "letterpaper",
    "pointsize": "11pt",
    "preamble": r"""
        \usepackage{amsmath}
        \usepackage{amssymb}
    """,
}

latex_documents = [
    (
        "index",
        "wasp2.tex",
        "WASP2 Documentation",
        "Aaron Ho, Jeff Jaureguy",
        "manual",
    ),
]
```

### index.rst Template (Landing Page)

```rst
WASP2 Documentation
===================

.. image:: _static/wasp2_logo.png
   :align: center
   :width: 400px
   :alt: WASP2 Logo

.. raw:: html

   <p style="text-align: center; font-size: 1.2em; margin: 20px 0;">
   High-performance allele-specific analysis of next-generation sequencing data
   </p>

----

.. grid:: 3
   :gutter: 3

   .. grid-item-card:: 🚀 Quick Start
      :link: quickstart
      :link-type: doc

      Get started with WASP2 in 5 minutes

   .. grid-item-card:: 📖 Tutorials
      :link: tutorials/index
      :link-type: doc

      Step-by-step guides for common workflows

   .. grid-item-card:: 📚 API Reference
      :link: api/index
      :link-type: doc

      Detailed API documentation

----

What is WASP2?
--------------

WASP2 is a comprehensive suite of tools for **allele-specific analysis** of
next-generation sequencing data. It addresses reference bias in read mapping
and provides statistical methods for detecting allelic imbalance.

Key Features
~~~~~~~~~~~~

.. grid:: 2
   :gutter: 2

   .. grid-item-card:: Unbiased Mapping
      :class-card: sd-border-1

      WASP algorithm corrects reference bias in RNA-seq, ATAC-seq, and ChIP-seq

   .. grid-item-card:: Statistical Testing
      :class-card: sd-border-1

      Beta-binomial models for rigorous allelic imbalance detection

   .. grid-item-card:: High Performance
      :class-card: sd-border-1

      Rust acceleration provides 10-25x speedup over pure Python

   .. grid-item-card:: Multi-Format Support
      :class-card: sd-border-1

      VCF, BCF, PGEN formats with up to 25x faster I/O

Applications
~~~~~~~~~~~~

- **RNA-seq**: Allele-specific expression (ASE) analysis
- **ATAC-seq**: Allele-specific chromatin accessibility
- **ChIP-seq**: Allele-specific transcription factor binding
- **Single-cell**: Cell-type-specific allelic imbalance

Quick Example
-------------

.. code-block:: bash

   # Install
   pip install wasp2

   # Count alleles
   wasp2-count count-variants sample.bam variants.vcf.gz \
     --samples NA12878 --out_file counts.tsv

   # Detect imbalance
   wasp2-analyze find-imbalance counts.tsv --out_file results.tsv

.. toctree::
   :maxdepth: 2
   :caption: Getting Started
   :hidden:

   installation
   quickstart
   concepts

.. toctree::
   :maxdepth: 2
   :caption: Tutorials
   :hidden:

   tutorials/index
   tutorials/basic_workflow
   tutorials/rnaseq_ase
   tutorials/atacseq_ase
   tutorials/single_cell
   tutorials/troubleshooting

.. toctree::
   :maxdepth: 2
   :caption: User Guide
   :hidden:

   user_guide/counting
   user_guide/mapping
   user_guide/analysis

.. toctree::
   :maxdepth: 2
   :caption: How-To Guides
   :hidden:

   how_to/index
   how_to/optimize_performance
   how_to/integrate_with_pipelines
   how_to/interpret_results

.. toctree::
   :maxdepth: 2
   :caption: API Reference
   :hidden:

   api/index
   api/counting
   api/mapping
   api/analysis
   api/io

.. toctree::
   :maxdepth: 2
   :caption: CLI Reference
   :hidden:

   cli/index
   cli/wasp2_count
   cli/wasp2_map
   cli/wasp2_analyze

.. toctree::
   :maxdepth: 2
   :caption: Background
   :hidden:

   explanations/index
   explanations/allelic_imbalance
   explanations/reference_bias
   explanations/wasp_algorithm
   explanations/statistical_models

.. toctree::
   :maxdepth: 1
   :caption: Reference
   :hidden:

   data_formats/index
   faq
   changelog
   citation
   development

Indices and Tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
```

---

## Quick Reference Card Template

Create `docs/CHEATSHEET.md`:

```markdown
# WASP2 Quick Reference

## Installation

```bash
pip install wasp2                    # Standard
pip install wasp2[cyvcf2,plink]     # With performance enhancements
```

## Common Commands

### Counting
```bash
# Basic
wasp2-count count-variants SAMPLE.bam VARIANTS.vcf.gz

# With sample filtering
wasp2-count count-variants SAMPLE.bam VARIANTS.vcf.gz -s SAMPLE_ID -o counts.tsv

# RNA-seq (with genes)
wasp2-count count-variants RNA.bam VARIANTS.pgen -s SAMPLE_ID -r genes.gtf -o gene_counts.tsv

# ATAC-seq (with peaks)
wasp2-count count-variants ATAC.bam VARIANTS.bcf -s SAMPLE_ID -r peaks.narrowPeak -o peak_counts.tsv
```

### Analysis
```bash
# Basic analysis
wasp2-analyze find-imbalance counts.tsv -o results.tsv

# With custom threshold
wasp2-analyze find-imbalance counts.tsv --min 20 -o results.tsv

# Gene-level analysis
wasp2-analyze find-imbalance gene_counts.tsv --groupby gene_id -o gene_results.tsv
```

### Mapping (WASP)
```bash
# Step 1: Generate reads for remapping
wasp2-map make-reads ORIGINAL.bam VARIANTS.vcf.gz -s SAMPLE_ID

# Step 2: Remap with your aligner (example with BWA)
bwa mem genome.fa *_swapped_alleles_r*.fq | samtools view -Sb - > remapped.bam

# Step 3: Filter remapped reads
wasp2-map filter-remapped remapped.bam to_remap.bam keep.bam -o wasp_filtered.bam
```

## Format Conversion

```bash
# VCF to BCF (5-8x faster)
bcftools view -O b variants.vcf.gz > variants.bcf

# VCF to PGEN (25x faster)
plink2 --vcf variants.vcf.gz --make-pgen --out variants
```

## Quick Diagnostics

```bash
# Check sample names in VCF
bcftools query -l variants.vcf.gz

# Count heterozygous SNPs
bcftools view -s SAMPLE -g het variants.vcf.gz | grep -v "^#" | wc -l

# Check BAM statistics
samtools flagstat sample.bam

# Check chromosome naming
samtools view -H sample.bam | grep "^@SQ" | head -3
bcftools view -h variants.vcf.gz | grep "^##contig" | head -3
```

## Common Patterns

```bash
# Process multiple samples
for sample in sample1 sample2 sample3; do
  wasp2-count count-variants ${sample}.bam variants.pgen -s ${sample} -o ${sample}_counts.tsv
done

# Process by chromosome
for chr in {1..22} X Y; do
  wasp2-count count-variants sample.bam variants.pgen --region chr${chr}.bed -o counts_chr${chr}.tsv
done

# Extract significant results (FDR < 0.05)
awk 'NR==1 || $8 < 0.05' results.tsv > significant.tsv
```

## Output Formats

### Counts (TSV)
```
chr     pos       ref  alt  ref_count  alt_count  other_count
chr10   1000000   A    G    12         15         0
```

### Analysis Results (TSV)
```
region           n_snps  ref_total  alt_total  p_value   fdr        log2_ratio
chr10:1M-1.5M    3       45         55         0.023     0.045      0.29
```

## Performance Tips

- Use PGEN for large files: `plink2 --vcf X.vcf.gz --make-pgen`
- Install cyvcf2: `pip install wasp2[cyvcf2]`
- Process by chromosome for very large datasets
- Use `--use-rust` (default) for 10-25x speedup

## Getting Help

```bash
wasp2-count --help
wasp2-count count-variants --help
wasp2-map --help
wasp2-analyze --help
```

Full documentation: https://jaureguy760.github.io/WASP2-exp/
```

---

This implementation guide provides copy-paste ready templates for all major documentation components. Use these as starting points and customize for specific WASP2 features and workflows.
