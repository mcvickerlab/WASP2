Counting Module
===============

Overview
--------

The counting module quantifies allele-specific read counts at heterozygous SNP positions. It's the first step in allelic imbalance analysis.

Purpose
~~~~~~~

* Count reads supporting reference vs alternate alleles
* Filter by sample genotype (heterozygous sites)
* Annotate with genomic regions (genes, peaks)
* Support single-cell RNA-seq

When to Use
~~~~~~~~~~~

Use counting when you have:
* Aligned reads (BAM file)
* Variant calls (VCF file)
* Want to quantify allele-specific expression

CLI Usage
---------

Basic Command
~~~~~~~~~~~~~

.. code-block:: bash

   wasp2-count count-variants BAM_FILE VCF_FILE

Full Options
~~~~~~~~~~~~

.. code-block:: bash

   wasp2-count count-variants \
     input.bam \
     variants.vcf \
     --samples sample1,sample2 \
     --region genes.gtf \
     --out_file counts.tsv

Input Requirements
------------------

BAM File
~~~~~~~~

* Aligned reads (single-end or paired-end)
* Indexed (.bai file in same directory)
* Sorted by coordinate

VCF File
~~~~~~~~

* Variant calls with genotype information
* Heterozygous SNPs (GT=0|1 or 1|0)
* Can include sample-specific genotypes

Optional: Region File
~~~~~~~~~~~~~~~~~~~~~

Annotate SNPs overlapping genes/peaks:

* GTF/GFF3 format (genes)
* BED format (peaks, regions)
* narrowPeak format (ATAC-seq, ChIP-seq)

Parameters
----------

``--samples`` / ``-s``
~~~~~~~~~~~~~~~~~~~~~~

Filter SNPs heterozygous in specified samples:

.. code-block:: bash

   --samples sample1,sample2,sample3
   # or
   --samples samples.txt  # one per line

``--region`` / ``-r``
~~~~~~~~~~~~~~~~~~~~~

Annotate SNPs with overlapping regions:

.. code-block:: bash

   --region genes.gtf      # Gene annotations
   --region peaks.bed      # ATAC-seq peaks
   --region regions.gff3   # Custom regions

``--out_file`` / ``-o``
~~~~~~~~~~~~~~~~~~~~~~~

Output file path (default: counts.tsv):

.. code-block:: bash

   --out_file my_counts.tsv

Output Format
-------------

Tab-separated file with columns:

Basic Columns
~~~~~~~~~~~~~

* ``chr``: Chromosome
* ``pos``: SNP position (1-based)
* ``ref``: Reference allele
* ``alt``: Alternate allele
* ``ref_count``: Reads supporting reference
* ``alt_count``: Reads supporting alternate
* ``other_count``: Reads supporting other alleles

Optional Columns (with --region)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* ``gene_id``: Overlapping gene
* ``gene_name``: Gene symbol
* ``feature``: Feature type (exon, intron, etc.)

Example Workflow
----------------

1. Basic Counting
~~~~~~~~~~~~~~~~~

.. code-block:: bash

   wasp2-count count-variants sample.bam variants.vcf

2. Filter by Sample
~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   wasp2-count count-variants \
     sample.bam \
     variants.vcf \
     --samples NA12878

3. Annotate with Genes
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   wasp2-count count-variants \
     sample.bam \
     variants.vcf \
     --samples NA12878 \
     --region genes.gtf \
     --out_file counts_annotated.tsv

Single-Cell Counting
--------------------

For single-cell RNA-seq:

.. code-block:: bash

   wasp2-count count-variants-sc \
     sc_rnaseq.bam \
     variants.vcf \
     --barcode_map barcodes.tsv

Output includes cell-type-specific counts.

Common Issues
-------------

Low Count Numbers
~~~~~~~~~~~~~~~~~

* Check BAM file coverage (``samtools depth``)
* Verify VCF contains heterozygous SNPs
* Ensure BAM and VCF use same reference genome

No Output SNPs
~~~~~~~~~~~~~~

* Check if --samples filter is too restrictive
* Verify VCF has genotype information (GT field)
* Ensure BAM file is indexed

Next Steps
----------

After counting:
* :doc:`analysis` - Detect allelic imbalance
* :doc:`mapping` - Correct reference bias with WASP
