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

   Path to variants. Supports VCF (.vcf, .vcf.gz), BCF (.bcf), and
   PLINK2 PGEN (.pgen) formats. VCF/BCF files should be indexed.

Options
~~~~~~~

**Input Filtering**

.. option:: -s <SAMPLE>, --samples <SAMPLE>

   Sample ID(s) to filter heterozygous SNPs.

   Accepts:
      - Comma-separated list: ``-s sample1,sample2``
      - File with one sample per line: ``-s samples.txt``

   If not provided, all variants are used regardless of genotype.

.. option:: -r <PATH>, --region <PATH>

   Filter SNPs overlapping genomic regions.

   Accepts:
      - BED format (``.bed``)
      - GTF format (``.gtf``)
      - GFF3 format (``.gff``, ``.gff3``)
      - narrowPeak format (``.narrowPeak``)

**Output**

.. option:: -o <PATH>, --out_file <PATH>

   Output file path. Default: ``counts.tsv``

.. option:: --temp_loc <DIR>

   Directory for intermediate files. If not specified, uses system
   temporary directory and removes files after completion.

**Region Annotation (GTF/GFF3)**

.. option:: --gene_feature <TYPE>

   Feature type from GTF/GFF3 to count. Default: ``exon``

   Examples: ``exon``, ``CDS``, ``five_prime_UTR``

.. option:: --gene_attribute <NAME>

   Attribute name for feature identifier.
   Default: ``gene_id`` (GTF), ``ID`` (GFF3)

.. option:: --gene_parent <NAME>

   Parent attribute for hierarchical features.
   Default: ``transcript_id`` (GTF), ``Parent`` (GFF3)

.. option:: --use_region_names

   Use region names (4th BED column) instead of coordinates in output.

**Performance**

.. option:: --use-rust / --no-rust

   Enable or disable Rust acceleration. Default: ``--use-rust``

.. option:: --include-indels

   Include indels in addition to SNPs. Default: SNPs only

**Advanced**

.. option:: --vcf-bed <PATH>

   Pre-computed VCF BED file (skip variant conversion)

.. option:: --intersect-bed <PATH>

   Pre-computed intersect BED file (skip intersection)

Examples
~~~~~~~~

**Basic counting**:

.. code-block:: bash

   wasp2-count count-variants sample.bam variants.vcf.gz

**Filter by sample**:

.. code-block:: bash

   wasp2-count count-variants sample.bam variants.vcf.gz \
       --samples NA12878 \
       --out_file counts.tsv

**RNA-seq with gene annotation**:

.. code-block:: bash

   wasp2-count count-variants sample.bam variants.vcf.gz \
       --samples NA12878 \
       --region genes.gtf \
       --gene_feature exon \
       --gene_attribute gene_id \
       --out_file gene_counts.tsv

**ATAC-seq with peaks**:

.. code-block:: bash

   wasp2-count count-variants sample.bam variants.vcf.gz \
       --samples NA12878 \
       --region peaks.narrowPeak \
       --use_region_names \
       --out_file peak_counts.tsv

**High-performance with PGEN**:

.. code-block:: bash

   wasp2-count count-variants sample.bam variants.pgen \
       --samples NA12878 \
       --out_file counts.tsv

Output Format
~~~~~~~~~~~~~

Tab-separated file with columns:

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
   * - ``region``
     - Overlapping region (if ``--region`` used)

Example output:

.. code-block:: text

   chr     pos       ref  alt  ref_count  alt_count  other_count  region
   chr10   1000000   A    G    12         15         0            ENSG00000123456
   chr10   1001000   C    T    20         18         1            ENSG00000123456

count-variants-sc
-----------------

Count allele-specific reads in single-cell data.

Synopsis
~~~~~~~~

.. code-block:: bash

   wasp2-count count-variants-sc [OPTIONS] BAM VARIANTS BARCODES

Arguments
~~~~~~~~~

.. option:: BAM

   Path to aligned reads (BAM format, typically from 10x Genomics).

.. option:: VARIANTS

   Path to variant file (VCF, BCF, or PGEN).

.. option:: BARCODES

   File with valid cell barcodes, one per line.

Options
~~~~~~~

.. option:: -s <SAMPLE>, --samples <SAMPLE>

   Sample ID for heterozygous SNP filtering. Recommended for single donor.

.. option:: -f <PATH>, --feature <PATH>

   Feature file (BED format) for region annotation.

.. option:: -o <PATH>, --out_file <PATH>

   Output file path. Default: ``allele_counts.h5ad``

.. option:: --temp_loc <DIR>

   Directory for intermediate files.

Example
~~~~~~~

.. code-block:: bash

   wasp2-count count-variants-sc \
       possorted_genome_bam.bam \
       variants.vcf.gz \
       filtered_barcodes.txt \
       --samples donor1 \
       --feature peaks.bed \
       --out_file sc_allele_counts.h5ad

Output Format
~~~~~~~~~~~~~

AnnData (h5ad) file containing:

* ``.X``: Cell x SNP count matrix (sparse)
* ``.var``: SNP annotations (chr, pos, ref, alt)
* ``.obs``: Cell annotations (barcode)

Performance Tips
----------------

Use High-Performance Formats
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For large variant files (>10M variants):

1. **PGEN format** (fastest, ~25x speedup):

   .. code-block:: bash

      plink2 --vcf variants.vcf.gz --make-pgen --out variants
      wasp2-count count-variants sample.bam variants.pgen

2. **cyvcf2 backend** (7x speedup for VCF):

   .. code-block:: bash

      pip install wasp2[cyvcf2]

3. **BCF format** (5-8x speedup):

   .. code-block:: bash

      bcftools view -O b variants.vcf.gz > variants.bcf

Process by Chromosome
~~~~~~~~~~~~~~~~~~~~~

For very large files:

.. code-block:: bash

   for chr in {1..22}; do
       wasp2-count count-variants sample.bam variants.vcf.gz \
           --region chr${chr}.bed \
           --out_file counts_chr${chr}.tsv
   done

   # Combine
   head -1 counts_chr1.tsv > all_counts.tsv
   tail -n +2 -q counts_chr*.tsv >> all_counts.tsv

See Also
--------

* :doc:`/api/counting` - Python API documentation
* :doc:`/tutorials/basic_workflow` - Complete workflow tutorial
* :doc:`wasp2_analyze` - Analyze allelic imbalance
