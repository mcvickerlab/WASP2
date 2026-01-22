CLI Reference
=============

WASP2 provides three command-line tools for allele-specific analysis.

.. contents:: Tools
   :local:
   :depth: 1

Overview
--------

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - Command
     - Purpose
     - Primary Use
   * - ``wasp2-count``
     - Count allele-specific reads
     - Quantify ref/alt reads at SNPs
   * - ``wasp2-map``
     - WASP mapping filter
     - Remove reference-biased reads
   * - ``wasp2-analyze``
     - Statistical analysis
     - Detect allelic imbalance

Getting Help
------------

All WASP2 commands support ``--help``:

.. code-block:: bash

   wasp2-count --help
   wasp2-count count-variants --help
   wasp2-map --help
   wasp2-analyze --help

Command Reference
-----------------

.. toctree::
   :maxdepth: 2

   wasp2_count
   wasp2_map
   wasp2_analyze

Quick Examples
--------------

Basic allele counting
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   wasp2-count count-variants sample.bam variants.vcf.gz \
       --samples NA12878 \
       --out_file counts.tsv

WASP mapping filter
~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Step 1: Generate swapped reads
   wasp2-map make-reads sample.bam variants.vcf.gz \
       --samples NA12878 \
       --out_dir wasp_output/

   # Step 2: Remap (use your aligner)
   bwa mem reference.fa wasp_output/sample_swapped_r1.fq wasp_output/sample_swapped_r2.fq \
       | samtools sort -o wasp_output/remapped.bam

   # Step 3: Filter
   wasp2-map filter-remapped wasp_output/remapped.bam \
       --json sample_wasp_data_files.json \
       --out_bam sample_wasp_filtered.bam

Statistical analysis
~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   wasp2-analyze find-imbalance counts.tsv \
       --min 10 \
       --out_file results.tsv

Single-cell analysis
~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Count per cell
   wasp2-count count-variants-sc sample.bam variants.vcf.gz barcodes.txt \
       --samples donor1 \
       --out_file sc_counts.h5ad

   # Analyze per cell type
   wasp2-analyze find-imbalance-sc sc_counts.h5ad celltype_map.tsv \
       --min 20 \
       --out_file ai_results

Environment Variables
---------------------

WASP2 behavior can be modified with environment variables:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Variable
     - Description
   * - ``WASP2_DISABLE_RUST``
     - Set to ``1`` to disable Rust acceleration
   * - ``TMPDIR``
     - Directory for temporary files (default: system temp)

Exit Codes
----------

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Code
     - Meaning
   * - 0
     - Success
   * - 1
     - General error (missing files, invalid arguments)
   * - 2
     - Data processing error (empty output, incompatible formats)
