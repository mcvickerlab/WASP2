wasp2-map
=========

Command-line interface for the WASP2 mapping module.

.. contents:: Commands
   :local:
   :depth: 2

Overview
--------

The ``wasp2-map`` command implements the WASP algorithm to remove
reference-biased reads. It provides two subcommands:

* ``make-reads`` - Generate reads with swapped alleles for remapping
* ``filter-remapped`` - Filter reads that don't remap consistently

Workflow
~~~~~~~~

.. code-block:: text

   Original BAM
       |
       v
   [make-reads]  ->  Swapped FASTQ + metadata
       |
       v
   [Your aligner]  ->  Remapped BAM
       |
       v
   [filter-remapped]  ->  Unbiased BAM

Global Options
--------------

.. option:: --help

   Show help message and exit

.. option:: --version

   Show version number and exit

make-reads
----------

Generate reads with swapped alleles for remapping.

Synopsis
~~~~~~~~

.. code-block:: bash

   wasp2-map make-reads [OPTIONS] BAM VARIANTS

Arguments
~~~~~~~~~

.. option:: BAM

   Path to aligned reads (BAM format). Must be sorted and indexed.

.. option:: VARIANTS

   Path to variant file (VCF, BCF, or PGEN).

Options
~~~~~~~

.. option:: --threads <N>

   Number of threads. Default: 1

.. option:: -s <SAMPLE>, --samples <SAMPLE>

   Sample ID(s) for filtering heterozygous SNPs.

.. option:: -o <DIR>, --out_dir <DIR>

   Output directory. Default: current directory

.. option:: --temp_loc <DIR>

   Directory for intermediate files.

.. option:: -j <PATH>, --out_json <PATH>

   Output JSON file with metadata.
   Default: ``[BAM_PREFIX]_wasp_data_files.json``

Output Files
~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - File
     - Description
   * - ``*_swapped_alleles_r1.fq``
     - Reads with swapped alleles (read 1)
   * - ``*_swapped_alleles_r2.fq``
     - Reads with swapped alleles (read 2, paired-end)
   * - ``*_to_remap.bam``
     - Original reads that need remapping
   * - ``*_keep.bam``
     - Reads not overlapping SNPs (keep unchanged)
   * - ``*_wasp_data_files.json``
     - Metadata for filter-remapped step

Example
~~~~~~~

.. code-block:: bash

   wasp2-map make-reads sample.bam variants.vcf.gz \
       --samples NA12878 \
       --out_dir wasp_output/ \
       --threads 4

filter-remapped
---------------

Filter reads that don't map consistently after allele swapping.

Synopsis
~~~~~~~~

.. code-block:: bash

   wasp2-map filter-remapped [OPTIONS] REMAPPED_BAM

   # OR

   wasp2-map filter-remapped [OPTIONS] REMAPPED_BAM TO_REMAP_BAM KEEP_BAM

Arguments
~~~~~~~~~

.. option:: REMAPPED_BAM

   BAM file with remapped reads (from your aligner).

.. option:: TO_REMAP_BAM (optional)

   Original reads that were remapped.
   Not needed if using ``--json``.

.. option:: KEEP_BAM (optional)

   Reads that didn't need remapping.
   Not needed if using ``--json``.

Options
~~~~~~~

.. option:: --threads <N>

   Number of threads. Default: 1

.. option:: -j <PATH>, --json <PATH>

   JSON metadata file from make-reads step.
   Default: ``[BAM_PREFIX]_wasp_data_files.json``

.. option:: -o <PATH>, --out_bam <PATH>

   Output filtered BAM file.
   Default: ``[BAM_PREFIX]_wasp_filt.bam``

.. option:: --remap_keep_bam <PATH>

   Output BAM with kept remapped reads (optional).

.. option:: --remap_keep_file <PATH>

   Output text file with kept read names (optional).

Examples
~~~~~~~~

**Using JSON metadata (recommended)**:

.. code-block:: bash

   wasp2-map filter-remapped wasp_output/remapped.bam \
       --json sample_wasp_data_files.json \
       --out_bam sample_wasp_filtered.bam

**Using explicit files**:

.. code-block:: bash

   wasp2-map filter-remapped \
       wasp_output/remapped.bam \
       wasp_output/sample_to_remap.bam \
       wasp_output/sample_keep.bam \
       --out_bam sample_wasp_filtered.bam

Complete Workflow Example
-------------------------

Step 1: Generate swapped reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   wasp2-map make-reads sample.bam variants.vcf.gz \
       --samples NA12878 \
       --out_dir wasp_output/ \
       --threads 4

Step 2: Remap with your aligner
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**BWA-MEM**:

.. code-block:: bash

   bwa mem -M reference.fa \
       wasp_output/sample_swapped_alleles_r1.fq \
       wasp_output/sample_swapped_alleles_r2.fq \
   | samtools view -b -F 4 - \
   | samtools sort -o wasp_output/remapped.bam -

   samtools index wasp_output/remapped.bam

**STAR (RNA-seq)**:

.. code-block:: bash

   STAR --runThreadN 8 \
       --genomeDir /path/to/star_index \
       --readFilesIn wasp_output/sample_swapped_alleles_r1.fq \
                     wasp_output/sample_swapped_alleles_r2.fq \
       --outSAMtype BAM SortedByCoordinate \
       --outFileNamePrefix wasp_output/remapped_

   # Rename for consistency
   mv wasp_output/remapped_Aligned.sortedByCoord.out.bam wasp_output/remapped.bam
   samtools index wasp_output/remapped.bam

**Bowtie2**:

.. code-block:: bash

   bowtie2 -x /path/to/bowtie2_index \
       -1 wasp_output/sample_swapped_alleles_r1.fq \
       -2 wasp_output/sample_swapped_alleles_r2.fq \
   | samtools view -b -F 4 - \
   | samtools sort -o wasp_output/remapped.bam -

   samtools index wasp_output/remapped.bam

Step 3: Filter inconsistent mappings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   wasp2-map filter-remapped wasp_output/remapped.bam \
       --json sample_wasp_data_files.json \
       --out_bam sample_wasp_filtered.bam \
       --threads 4

Important Notes
---------------

Aligner Settings
~~~~~~~~~~~~~~~~

**Critical**: Use the same aligner and parameters for remapping as you
used for the original alignment. Different settings may cause all reads
to be filtered.

Read Names
~~~~~~~~~~

Do not modify read names during remapping. WASP matching requires
original read names to identify corresponding reads.

Single-End vs Paired-End
~~~~~~~~~~~~~~~~~~~~~~~~

WASP2 automatically detects paired-end data. Both reads in a pair
must remap consistently to be kept.

Performance
~~~~~~~~~~~

- The make-reads step is I/O bound
- The filter-remapped step benefits from multiple threads
- Consider using fast local storage for temporary files

See Also
--------

* :doc:`/tutorials/concepts` - WASP algorithm explanation
* :doc:`/tutorials/basic_workflow` - Complete workflow tutorial
* :doc:`wasp2_count` - Count alleles from filtered BAM
