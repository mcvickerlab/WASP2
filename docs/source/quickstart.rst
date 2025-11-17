Quick Start
===========

This 5-minute tutorial demonstrates basic WASP2 usage.

Example Data
------------

Use the included test data:

.. code-block:: bash

   cd WASP2-exp
   ls test_data/

Count Alleles
-------------

Count allele-specific reads from a BAM file:

.. code-block:: bash

   wasp2-count count-variants \
     test_data/CD4_ATACseq_Day1_merged_filtered.sort.bam \
     test_data/filter_chr10.vcf \
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
