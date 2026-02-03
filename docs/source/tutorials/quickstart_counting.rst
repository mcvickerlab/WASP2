Quickstart: Count Alleles in 5 Minutes
======================================

This tutorial demonstrates the basic WASP2 allele counting workflow using a minimal test dataset.

**What you'll learn:**

- How to count allele-specific reads from a BAM file
- Basic WASP2 command-line usage
- Understanding the output format

**Prerequisites:**

- WASP2 installed (``pip install wasp2``)
- Basic familiarity with BAM and VCF file formats

Setup
-----

First, verify WASP2 is installed:

.. code-block:: bash

   wasp2-count --version

Test Data
---------

We'll use the minimal test data included in the WASP2 repository:

- **BAM file**: Synthetic paired-end reads overlapping heterozygous variants
- **VCF file**: 6 variants with genotypes for two samples
- **GTF file**: Gene annotations for 3 genes

The test data is located in ``pipelines/nf-modules/tests/data/``.

**VCF contents:**

.. code-block:: text

   #CHROM  POS  ID   REF  ALT  QUAL  FILTER  INFO    FORMAT  sample1  sample2
   chr1    100  rs1  A    G    30    PASS    DP=50   GT      0/1      0/0
   chr1    200  rs2  C    T    30    PASS    DP=45   GT      1/1      0/1
   chr1    300  rs3  G    A    30    PASS    DP=60   GT      0/0      1/1
   chr1    400  rs4  T    C    30    PASS    DP=55   GT      0/1      0/1
   chr2    100  rs5  A    T    30    PASS    DP=40   GT      0/1      0/0
   chr2    200  rs6  G    C    30    PASS    DP=35   GT      ./.      0/1

The ``GT`` field shows genotypes:

- ``0/1``: Heterozygous (has both reference and alternate alleles)
- ``0/0``: Homozygous reference
- ``1/1``: Homozygous alternate

For allele-specific analysis, we focus on **heterozygous sites** (0/1).

Step 1: Basic Allele Counting
-----------------------------

The simplest way to count alleles is to provide a BAM file and VCF file:

.. code-block:: bash

   wasp2-count count-variants \
     pipelines/nf-modules/tests/data/minimal.bam \
     pipelines/nf-modules/tests/data/sample.vcf.gz \
     --out_file counts_basic.tsv

**Output:**

.. code-block:: text

   chr   pos  ref  alt  ref_count  alt_count  other_count
   chr1  100  A    G    1          0          0
   chr1  400  T    C    1          0          0
   chr2  100  A    T    1          0          0

Output Columns
~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Column
     - Description
   * - ``chr``
     - Chromosome
   * - ``pos``
     - Variant position (1-based)
   * - ``ref``
     - Reference allele
   * - ``alt``
     - Alternate allele
   * - ``ref_count``
     - Reads supporting reference allele
   * - ``alt_count``
     - Reads supporting alternate allele
   * - ``other_count``
     - Reads with other alleles (errors, indels)

Step 2: Filter by Sample
------------------------

When your VCF contains multiple samples, use ``--samples`` to filter for heterozygous sites in a specific sample:

.. code-block:: bash

   wasp2-count count-variants \
     pipelines/nf-modules/tests/data/minimal.bam \
     pipelines/nf-modules/tests/data/sample.vcf.gz \
     --samples sample1 \
     --out_file counts_sample1.tsv

This returns only the 3 sites where sample1 is heterozygous:

- chr1:100 (rs1)
- chr1:400 (rs4)
- chr2:100 (rs5)

Step 3: Annotate with Gene Regions
----------------------------------

Use ``--region`` to annotate variants with overlapping genomic features (genes, peaks, etc.):

.. code-block:: bash

   wasp2-count count-variants \
     pipelines/nf-modules/tests/data/minimal.bam \
     pipelines/nf-modules/tests/data/sample.vcf.gz \
     --samples sample1 \
     --region pipelines/nf-modules/tests/data/sample.gtf \
     --out_file counts_annotated.tsv

The output now includes gene annotations from the GTF file, allowing you to aggregate counts per gene for downstream analysis.

Next Steps
----------

Now that you have allele counts, you can:

1. **Analyze allelic imbalance** using ``wasp2-analyze find-imbalance``
2. **Compare between conditions** using ``wasp2-analyze compare-imbalance``
3. **Correct mapping bias** using ``wasp2-map`` (for WASP-filtered BAMs)

See Also
--------

* :doc:`/user_guide/counting` - Detailed counting options
* :doc:`/tutorials/scrna_seq` - Single-cell RNA-seq tutorial
* :doc:`/tutorials/comparative_imbalance` - Differential imbalance analysis
