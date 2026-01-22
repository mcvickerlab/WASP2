Frequently Asked Questions
==========================

Common questions about WASP2 and allele-specific analysis.

.. contents:: Categories
   :local:
   :depth: 2

General Questions
-----------------

What is allelic imbalance?
~~~~~~~~~~~~~~~~~~~~~~~~~~

Allelic imbalance (AI) occurs when one allele of a heterozygous variant is
preferentially expressed or accessible compared to the other allele. In a
balanced state, you'd expect a 50:50 ratio of reads from each allele. AI
can indicate:

- **Cis-regulatory variants**: SNPs affecting gene regulation
- **Imprinting**: Parent-of-origin specific expression
- **X-inactivation**: Random silencing of one X chromosome
- **Technical artifacts**: Mapping bias (which WASP corrects)

When should I use WASP2 vs GATK ASEReadCounter?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use **WASP2** if you:

- Need reference bias correction (WASP mapping)
- Are analyzing single-cell data
- Want statistical testing for allelic imbalance
- Need high performance (Rust acceleration)
- Want support for multiple variant formats (VCF, BCF, PGEN)

Use **GATK ASEReadCounter** if you:

- Only need raw allele counts without statistical analysis
- Are already using GATK workflows
- Don't need mapping bias correction

Do I need to run WASP mapping before counting?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**It depends on your use case:**

- **Yes, use WASP** if you used standard aligners (STAR, BWA, bowtie2)
  and want accurate allelic ratios
- **Maybe not needed** if you used allele-aware aligners or are doing
  a quick exploratory analysis

**Rule of thumb**: If in doubt, run WASP mapping. It's conservative and
won't hurt accuracy - it just removes some reads.

What's the difference between WASP and WASP2?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**WASP** (original):
- Python-only implementation
- VCF-only variant support
- Separate tools for each step

**WASP2**:
- Rust-accelerated core (10-25x faster)
- Multiple variant formats (VCF, BCF, PGEN)
- Unified CLI interface
- Single-cell support
- Better statistical models
- Active development

Installation Questions
----------------------

Installation fails with "Rust compiler not found"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Install the Rust compiler using rustup:

.. code-block:: bash

   curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
   source $HOME/.cargo/env

   # Retry WASP2 installation
   pip install wasp2

Can I install WASP2 without Rust?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Yes, but you'll miss significant performance benefits. WASP2 includes
Python fallbacks for all Rust-accelerated functions.

.. code-block:: bash

   # Disable Rust extension
   export WASP2_DISABLE_RUST=1
   pip install wasp2

Performance will be 10-25x slower for counting operations.

How do I install optional dependencies?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Faster VCF parsing
   pip install wasp2[cyvcf2]

   # PLINK2 format support
   pip install wasp2[plink]

   # All optional dependencies
   pip install wasp2[all]

Data Format Questions
---------------------

What variant formats does WASP2 support?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 20 30 30

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
     - Large cohorts

How do I convert VCF to PGEN?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Install plink2
   # Download from https://www.cog-genomics.org/plink/2.0/

   # Convert VCF to PGEN
   plink2 --vcf variants.vcf.gz --make-pgen --out variants

   # Use in WASP2
   wasp2-count count-variants sample.bam variants.pgen --samples sample1

Do BAM and VCF need to use the same reference genome?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Yes, absolutely**. Mismatched reference genomes will cause:

- Missing SNPs (different coordinates)
- Incorrect counts (different alleles)
- Chromosome naming issues (chr10 vs 10)

Verify your references:

.. code-block:: bash

   # Check BAM header
   samtools view -H sample.bam | grep "@SQ" | head -3

   # Check VCF header
   bcftools view -h variants.vcf.gz | grep "##contig" | head -3

What region file formats are supported?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- **BED**: Tab-delimited (chr, start, end, name)
- **narrowPeak**: MACS2 peak format
- **GTF**: Gene annotations (GENCODE, Ensembl)
- **GFF3**: General feature format

Analysis Questions
------------------

How many reads do I need for allelic imbalance analysis?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Minimum recommendations**:

- **Per SNP**: ≥10 total reads (5 per allele)
- **Per gene/peak**: ≥20 reads across all SNPs
- **For single-cell**: ≥100 cells per cell type

More reads = higher statistical power to detect imbalance.

What does "FDR < 0.05" mean?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

False Discovery Rate (FDR) is the expected proportion of false positives
among significant results.

- **FDR < 0.05**: Expect <5% of "significant" results to be false positives
- **FDR < 0.01**: More stringent, <1% false positives

Use FDR instead of raw p-values when testing many genes/peaks.

Why are some genes significant but with weak imbalance?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

High coverage genes can show statistical significance even with small
allelic ratios (e.g., 55:45 instead of 50:50).

**Interpretation**:

- **Statistical significance** (FDR < 0.05): Effect is real, not random
- **Biological significance**: Depends on effect size and context

Filter by effect size for biologically relevant results:

.. code-block:: bash

   # Genes with strong imbalance (ratio <0.35 or >0.65)
   awk 'NR==1 || ($7 < 0.05 && ($5 < 0.35 || $5 > 0.65))' results.tsv

What statistical model does WASP2 use?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

WASP2 uses the **beta-binomial distribution** which accounts for:

- **Overdispersion**: More variation than simple binomial predicts
- **Technical noise**: PCR bias, sequencing errors
- **Biological variation**: Stochastic gene expression

This provides more accurate p-values than simple binomial tests.

Single-Cell Questions
---------------------

How should I handle low coverage in single cells?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Strategies**:

1. **Aggregate by cell type**: Combine cells before analysis
2. **Lower threshold**: Use ``--min 5`` instead of default 10
3. **Filter features**: Only analyze high-coverage peaks/genes
4. **Pseudobulk**: Sum counts across cells of same type

Can I analyze multiple donors in single-cell data?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Single-cell barcodes are sample-specific, so analyzing multiple donors requires:

1. **Demultiplexing**: Assign cells to samples
2. **Sample-specific counting**: Run ``count-variants-sc`` per sample
3. **Combined analysis**: Merge results with sample labels

For now, **analyze one donor at a time** and combine results downstream.

Performance Questions
---------------------

WASP2 is running slowly. How can I speed it up?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Install cyvcf2** (7x faster VCF parsing):

   .. code-block:: bash

      pip install wasp2[cyvcf2]

2. **Use PGEN format** (25x faster):

   .. code-block:: bash

      plink2 --vcf variants.vcf.gz --make-pgen --out variants
      wasp2-count count-variants sample.bam variants.pgen ...

3. **Process by chromosome** for very large files:

   .. code-block:: bash

      for chr in {1..22}; do
          wasp2-count count-variants sample.bam variants.pgen \
              --region chr${chr}.bed --out_file counts_chr${chr}.tsv
      done

I'm running out of memory. What can I do?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Use PGEN format** (lower memory footprint)
2. **Process chromosomes separately**
3. **Reduce threads**: ``--threads 1``
4. **Pre-filter VCF** to heterozygous SNPs only:

   .. code-block:: bash

      bcftools view -s sample1 -g het variants.vcf.gz -O z -o het_only.vcf.gz

Troubleshooting
---------------

See :doc:`troubleshooting` for detailed error messages and solutions.

"Sample not found in VCF" error
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # List samples in VCF
   bcftools query -l variants.vcf.gz

   # Use exact sample name
   wasp2-count count-variants sample.bam variants.vcf.gz \
       --samples "EXACT_SAMPLE_NAME"

No output SNPs
~~~~~~~~~~~~~~

**Diagnostic**:

.. code-block:: bash

   # Check for heterozygous SNPs
   bcftools view -s sample1 -g het variants.vcf.gz | bcftools view -H | wc -l

   # Check chromosome naming
   samtools view -H sample.bam | grep "@SQ" | head -1
   bcftools view -h variants.vcf.gz | grep "##contig" | head -1

**Solutions**:

- Verify sample name spelling
- Check chromosome naming matches (chr1 vs 1)
- Ensure VCF contains heterozygous genotypes

Getting Help
------------

If your question isn't answered here:

1. Check the :doc:`troubleshooting` guide
2. Search `GitHub Issues <https://github.com/Jaureguy760/WASP2-exp/issues>`_
3. Open a new issue with:

   - WASP2 version (``wasp2-count --version``)
   - Operating system
   - Complete error message
   - Minimal example to reproduce
