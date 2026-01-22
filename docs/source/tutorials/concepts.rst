Understanding Allele-Specific Analysis
======================================

**Time**: 15 minutes reading
**Level**: Beginner
**Prerequisites**: Basic understanding of genetics and NGS

This tutorial introduces the key concepts behind allele-specific analysis,
explains why reference bias matters, and describes how WASP2 addresses these
challenges.

.. contents:: Topics
   :local:
   :depth: 2

What is Allelic Imbalance?
--------------------------

In diploid organisms like humans, each individual carries two copies (alleles) of
most genes - one inherited from each parent. Under normal circumstances, both
alleles are expressed equally, producing a 50:50 ratio of transcripts.

**Allelic imbalance (AI)** occurs when one allele is preferentially expressed,
resulting in an unequal ratio (e.g., 70:30). This can indicate important
biological phenomena:

Causes of Allelic Imbalance
~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Cis-regulatory variants**
   Single nucleotide polymorphisms (SNPs) that affect transcription factor
   binding sites can alter expression of the nearby allele.

2. **Imprinting**
   Some genes show parent-of-origin-specific expression, where only the maternal
   or paternal allele is expressed.

3. **X-chromosome inactivation**
   In females, one X chromosome is randomly inactivated in each cell, leading
   to monoallelic expression of X-linked genes.

4. **Allele-specific methylation**
   Epigenetic modifications can silence specific alleles.

5. **Chromatin accessibility differences**
   Allele-specific chromatin states can affect transcription factor access.

Why Does Reference Bias Matter?
-------------------------------

When sequencing reads are aligned to a reference genome, there's an inherent
bias toward reads that match the reference allele:

The Problem
~~~~~~~~~~~

Consider a heterozygous SNP where the individual has both reference (A) and
alternate (G) alleles:

.. code-block:: text

   Reference genome:    ...ACGT[A]CGTA...
   Read with ref (A):   ...ACGT[A]CGTA...  -> Perfect match, maps well
   Read with alt (G):   ...ACGT[G]CGTA...  -> 1 mismatch, may fail to map

This creates **artificial allelic imbalance** because:

- Reads with alternate alleles have more mismatches
- More mismatches can lower alignment scores
- Lower scores may cause reads to be filtered or unmapped
- Result: Artificially inflated reference allele counts

The Impact
~~~~~~~~~~

Reference bias can cause:

- **False positives**: Genes appear to have allelic imbalance when they don't
- **False negatives**: Real imbalance is masked by the artificial bias
- **Quantitative errors**: Effect sizes are over- or under-estimated

Studies have shown reference bias can cause 5-15% over-representation of
reference alleles, depending on the aligner and parameters used.

The WASP Solution
-----------------

WASP (WASP Allele-Specific Pipeline) corrects reference bias through a
clever read-swapping approach:

Algorithm Overview
~~~~~~~~~~~~~~~~~~

For each read overlapping a heterozygous SNP:

1. **Identify overlap**: Find reads that cover variant positions
2. **Swap alleles**: Create a modified read with the alternate allele
3. **Re-map**: Align the swapped read to the reference genome
4. **Filter**: Keep only reads that map to the same location

.. code-block:: text

   Original read:   ...ACGT[A]CGTA...  -> maps to chr1:1000
   Swapped read:    ...ACGT[G]CGTA...  -> maps to chr1:1000? Keep!
                                       -> maps elsewhere? Discard!

Why This Works
~~~~~~~~~~~~~~

Reads that map identically regardless of which allele they carry are
**unbiased** - their mapping is not affected by the variant. By keeping
only these reads, WASP creates an unbiased set for allele counting.

**Trade-off**: WASP discards some reads, reducing coverage. However, the
remaining reads provide accurate allelic ratios.

When to Use Each WASP2 Module
-----------------------------

WASP2 provides three main modules for different stages of analysis:

Decision Tree
~~~~~~~~~~~~~

.. code-block:: text

   Do you have aligned reads (BAM)?
   ├── YES: Have they been WASP-filtered?
   │   ├── YES: Use wasp2-count directly
   │   └── NO: Run wasp2-map first, then wasp2-count
   └── NO: Align reads first (STAR, BWA, etc.), then decide above

   Do you have allele counts?
   ├── YES: Use wasp2-analyze
   └── NO: Run wasp2-count first

Module Summary
~~~~~~~~~~~~~~

**wasp2-count** (Counting Module)
   - Input: BAM file + VCF file
   - Output: Allele counts per SNP
   - Use when: You have aligned reads and want allele counts
   - Note: Works best with WASP-filtered BAM for unbiased results

**wasp2-map** (Mapping Module)
   - Input: BAM file + VCF file
   - Output: WASP-filtered BAM file
   - Use when: You need to correct reference bias
   - Note: Requires re-mapping step with your aligner

**wasp2-analyze** (Analysis Module)
   - Input: Allele counts (from wasp2-count)
   - Output: Statistical test results
   - Use when: You want to identify significant allelic imbalance
   - Note: Uses beta-binomial model for robust statistics

Statistical Models
------------------

WASP2 uses the **beta-binomial distribution** for testing allelic imbalance.

Why Beta-Binomial?
~~~~~~~~~~~~~~~~~~

Simple binomial testing assumes each read is independent, but sequencing
data has **overdispersion** - more variation than the binomial predicts due to:

- PCR amplification artifacts
- Batch effects
- Biological variation between cells

The beta-binomial model accounts for this extra variation, providing more
accurate p-values and fewer false positives.

Interpretation
~~~~~~~~~~~~~~

WASP2 reports:

- **Allelic ratio**: Proportion of reference reads (0.5 = balanced)
- **P-value**: Probability of observing this imbalance by chance
- **FDR (q-value)**: False discovery rate adjusted p-value

General guidelines:

- **FDR < 0.05**: Significant allelic imbalance
- **Ratio > 0.6 or < 0.4**: Moderate effect size
- **Ratio > 0.7 or < 0.3**: Strong effect size

Key Terminology
---------------

.. glossary::

   Allelic Imbalance (AI)
      Unequal expression or accessibility of two alleles at a heterozygous site.

   Reference Bias
      Systematic over-representation of reference alleles in sequencing data
      due to alignment algorithms favoring reference-matching reads.

   Heterozygous SNP
      A position where an individual carries two different alleles.

   WASP Filtering
      The process of removing reads whose mapping is affected by variant alleles.

   Beta-binomial Distribution
      A statistical distribution that models count data with overdispersion,
      more appropriate than binomial for sequencing data.

   Cis-regulatory Variant
      A genetic variant that affects expression of a nearby gene on the same
      chromosome.

   eQTL (expression QTL)
      A genetic locus that affects gene expression levels.

   caQTL (chromatin accessibility QTL)
      A genetic locus that affects chromatin accessibility.

Further Reading
---------------

**Original WASP Paper**
   van de Geijn B, et al. (2015). WASP: allele-specific software for robust
   molecular quantitative trait locus discovery. *Nature Methods* 12:1061-1063.
   `doi:10.1038/nmeth.3582 <https://doi.org/10.1038/nmeth.3582>`_

**Reference Bias Studies**
   - Degner JF, et al. (2009). Effect of read-mapping biases on detecting
     allele-specific expression from RNA-sequencing data. *Bioinformatics*.
   - Castel SE, et al. (2015). Tools and best practices for data processing
     in allelic expression analysis. *Genome Biology*.

**Beta-Binomial Model**
   - Skelly DA, et al. (2011). A powerful and flexible statistical framework
     for testing hypotheses of allele-specific gene expression from RNA-seq data.
     *Genome Research*.

Next Steps
----------

Now that you understand the concepts, proceed to:

- :doc:`/quickstart` - Try WASP2 with example data
- :doc:`basic_workflow` - Learn the complete pipeline
- :doc:`rnaseq_ase` - RNA-seq allele-specific expression
- :doc:`atacseq_ase` - ATAC-seq allelic accessibility
