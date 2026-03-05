Quick Start
===========

This 5-minute tutorial demonstrates basic WASP2 usage.

Prerequisites
-------------

You will need:

* A coordinate-sorted, indexed BAM file (``sample.bam`` + ``sample.bam.bai``)
* A phased VCF file with heterozygous variants (``variants.vcf.gz`` + ``.tbi``)

These are typically produced by your alignment pipeline (BWA-MEM, STAR, etc.)
followed by variant calling and phasing (GATK, WhatsHap, ShapeIt).

Count Alleles
-------------

Count allele-specific reads from a BAM file:

.. code-block:: bash

   wasp2-count count-variants \
     sample.bam \
     variants.vcf.gz \
     -s SAMPLE_ID \
     --out_file counts.tsv

Output: ``counts.tsv`` with columns:

* chr, pos, ref, alt
* ref_count, alt_count, other_count

Analyze Allelic Imbalance
--------------------------

Detect significant allelic imbalance:

.. code-block:: bash

   wasp2-analyze find-imbalance \
     counts.tsv \
     --output results.tsv

Output: ``results.tsv`` with columns:

* region, ref_count, alt_count
* p-value, FDR-corrected p-value
* Statistical metrics

Interpret Results
-----------------

Significant imbalance (FDR < 0.05) indicates:

* Preferential expression of one allele
* Potential cis-regulatory variation
* Technical artifacts (check coverage)

Next Steps
----------

* :doc:`user_guide/counting` - Detailed counting options
* :doc:`user_guide/mapping` - WASP remapping workflow
* :doc:`user_guide/analysis` - Statistical models
