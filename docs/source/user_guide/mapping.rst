Mapping Module
==============

Overview
--------

``wasp2-map`` implements the WASP remap-and-filter workflow for removing
reference mapping bias before allele counting.

The public CLI has two commands:

1. ``make-reads``: find reads overlapping sample variants and generate
   allele-swapped FASTQ files for remapping
2. ``filter-remapped``: keep only remapped reads that return to the same locus

There is no separate ``find-intersecting-snps`` command in WASP2. That overlap
step is part of ``make-reads``.

Typical Workflow
----------------

Step 1: Generate swapped reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   wasp2-map make-reads \
     sample.bam \
     variants.vcf.gz \
     --samples SAMPLE1 \
     --out_dir wasp_output

This writes:

* ``sample_to_remap.bam``: original reads that must be remapped
* ``sample_keep.bam``: reads that never overlapped eligible variants
* ``sample_swapped_alleles_r1.fq`` and ``sample_swapped_alleles_r2.fq``:
  swapped FASTQ reads to realign
* ``sample_wasp_data_files.json``: metadata for ``filter-remapped``

Step 2: Realign swapped reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use the same aligner and alignment settings used for the original BAM.

.. code-block:: bash

   bwa mem -M -t 8 genome.fa \
     wasp_output/sample_swapped_alleles_r1.fq \
     wasp_output/sample_swapped_alleles_r2.fq | \
     samtools sort -o wasp_output/sample_remapped.bam -

   samtools index wasp_output/sample_remapped.bam

Step 3: Filter remapped reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   wasp2-map filter-remapped \
     wasp_output/sample_remapped.bam \
     --wasp_data_json wasp_output/sample_wasp_data_files.json \
     --out_bam wasp_output/sample_wasp_filtered.bam

You can also provide ``to_remap_bam`` and ``keep_bam`` positionally instead of
``--wasp_data_json``.

Command Reference
-----------------

``make-reads``
~~~~~~~~~~~~~~

.. code-block:: bash

   wasp2-map make-reads [OPTIONS] BAM VARIANTS

Important options:

* ``--samples`` / ``-s``: sample name(s) used to select het variants
* ``--out_dir`` / ``-o``: output directory
* ``--out_json`` / ``-j``: explicit metadata JSON path
* ``--indels``: include indels as well as SNPs
* ``--threads``: BAM I/O threads

Notes:

* paired-end input is required
* phased genotypes are strongly recommended
* supported variant formats are VCF, VCF.GZ, BCF, and PGEN

``filter-remapped``
~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   wasp2-map filter-remapped [OPTIONS] REMAPPED_BAM [TO_REMAP_BAM] [KEEP_BAM]

Important options:

* ``--wasp_data_json`` / ``-j``: load ``to_remap_bam`` and ``keep_bam`` from
  ``make-reads`` metadata
* ``--out_bam`` / ``-o``: output BAM path
* ``--remap_keep_bam``: optional BAM of remapped reads that passed filtering
* ``--remap_keep_file``: optional text file of kept read names
* ``--same-locus-slop``: positional tolerance for same-locus matching
* ``--threads``: BAM I/O threads

Interpreting Outputs
--------------------

Common outcomes after ``filter-remapped``:

* reads kept because they remap to the same locus
* reads dropped because they remap elsewhere
* reads dropped because they fail to remap cleanly

The final WASP-corrected BAM is the output of ``filter-remapped`` merged with
the ``*_keep.bam`` reads that never required remapping.

Example
-------

.. code-block:: bash

   wasp2-map make-reads \
     sample.bam \
     variants.vcf.gz \
     --samples SAMPLE1 \
     --out_dir wasp_output

   bwa mem -M -t 8 genome.fa \
     wasp_output/sample_swapped_alleles_r1.fq \
     wasp_output/sample_swapped_alleles_r2.fq | \
     samtools sort -o wasp_output/sample_remapped.bam -

   samtools index wasp_output/sample_remapped.bam

   wasp2-map filter-remapped \
     wasp_output/sample_remapped.bam \
     --wasp_data_json wasp_output/sample_wasp_data_files.json \
     --out_bam wasp_output/sample_wasp_filtered.bam

Next Steps
----------

* :doc:`counting` to count alleles from the WASP-filtered BAM
* :doc:`analysis` to test for allelic imbalance
