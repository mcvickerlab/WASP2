Mapping Module (WASP)
=====================

Overview
--------

The WASP (Weighted Allele-Specific Mapping) algorithm corrects reference bias by remapping reads with all possible alleles.

What is Reference Bias?
~~~~~~~~~~~~~~~~~~~~~~~~

Reference bias occurs when reads containing alternate alleles align worse than reads with reference alleles, leading to false allelic imbalance signals.

WASP Solution
~~~~~~~~~~~~~

1. Identify reads overlapping heterozygous SNPs
2. Generate alternative reads (swap alleles)
3. Remap both original and swapped reads
4. Keep only reads that map to the same location

Purpose
-------

* Correct reference bias in RNA-seq, ATAC-seq
* Improve accuracy of allelic imbalance detection
* Required before allele counting

When to Use
~~~~~~~~~~~

Use WASP when:
* Reads will be used for allelic analysis
* Reference genome differs from sample genotype
* High-confidence bias correction needed

Workflow
--------

Complete WASP workflow has 3 steps:

Step 1: Find Intersecting SNPs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Identify reads overlapping heterozygous SNPs:

.. code-block:: bash

   wasp2-map find-intersecting-snps \
     input.bam \
     variants.vcf \
     --output intersecting.bam

Output: BAM file with reads overlapping SNPs.

Step 2: Generate Remapping Reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create reads with swapped alleles:

.. code-block:: bash

   wasp2-map make-reads \
     intersecting.bam \
     variants.vcf \
     --samples sample1 \
     --output remap_reads.fastq

Output: FASTQ file(s) with alternative allele sequences.

Step 3: Remap and Filter
~~~~~~~~~~~~~~~~~~~~~~~~~

User remaps with their aligner (BWA, STAR, etc.):

.. code-block:: bash

   # Example with BWA
   bwa mem -t 8 reference.fa remap_reads.fastq | \
     samtools sort -o remapped.bam -

Then filter to consistent mappings:

.. code-block:: bash

   wasp2-map filt-remapped-reads \
     intersecting.bam \
     remapped.bam \
     --output filtered.bam

Output: BAM file with bias-corrected reads.

CLI Reference
-------------

find-intersecting-snps
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   wasp2-map find-intersecting-snps [OPTIONS] BAM VCF

Options:
* ``--samples``: Filter by sample genotype
* ``--output``: Output BAM file

make-reads
~~~~~~~~~~

.. code-block:: bash

   wasp2-map make-reads [OPTIONS] BAM VCF

Options:
* ``--samples``: Sample name(s)
* ``--output``: Output FASTQ prefix
* ``--paired``: Paired-end mode

filt-remapped-reads
~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   wasp2-map filt-remapped-reads [OPTIONS] ORIGINAL REMAPPED

Options:
* ``--output``: Filtered BAM file
* ``--keep_read_file``: Save kept read IDs

Input Requirements
------------------

* **Original BAM**: Aligned reads from initial mapping
* **VCF File**: Phased heterozygous SNPs (recommended)
* **Reference Genome**: Same as used for original alignment

Output Interpretation
---------------------

WASP Filter Rate
~~~~~~~~~~~~~~~~

Typical filter rates:
* **Good**: 95-99% reads kept
* **Acceptable**: 90-95% reads kept
* **Concerning**: <90% reads kept (check data quality)

Low filter rate may indicate:
* Poor mapping quality
* High SNP density
* Problematic reference genome

Complete Example
----------------

Full WASP workflow:

.. code-block:: bash

   # Step 1: Find SNP-overlapping reads
   wasp2-map find-intersecting-snps \
     original.bam \
     phased_variants.vcf \
     --samples NA12878 \
     --output intersecting.bam

   # Step 2: Generate remapping reads
   wasp2-map make-reads \
     intersecting.bam \
     phased_variants.vcf \
     --samples NA12878 \
     --paired \
     --output remap

   # Step 3: Remap (user's aligner)
   bwa mem -t 16 reference.fa \
     remap_R1.fastq remap_R2.fastq | \
     samtools sort -o remapped.bam -
   samtools index remapped.bam

   # Step 4: Filter
   wasp2-map filt-remapped-reads \
     intersecting.bam \
     remapped.bam \
     --output filtered_wasp.bam

   # Step 5: Count alleles (use filtered BAM)
   wasp2-count count-variants \
     filtered_wasp.bam \
     phased_variants.vcf \
     --samples NA12878

Performance Tips
----------------

* Use multi-threading for remapping step
* Filter VCF to high-quality SNPs only
* Use phased genotypes when available

Common Issues
-------------

Many Reads Filtered
~~~~~~~~~~~~~~~~~~~~

* Check remapping quality (MAPQ scores)
* Verify same reference genome used
* Consider relaxing mapping parameters

Slow Remapping
~~~~~~~~~~~~~~

* Use multi-threading (``-t`` flag)
* Process chromosomes in parallel
* Consider downsampling for testing

Next Steps
----------

* :doc:`counting` - Count alleles from WASP-filtered BAM
* :doc:`analysis` - Analyze allelic imbalance
