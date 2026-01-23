Analysis Module
===============

Overview
--------

The analysis module detects statistically significant allelic imbalance using beta-binomial models.

Purpose
-------

* Detect allelic imbalance at genomic regions
* Control for biological and technical variation
* Support single-cell and bulk RNA-seq
* Compare imbalance between groups/conditions

Statistical Models
------------------

Beta-Binomial Model
~~~~~~~~~~~~~~~~~~~

WASP2 uses beta-binomial distribution to model:
* Overdispersion (variation beyond binomial)
* Biological variability between regions
* Technical noise in sequencing

The model:
* Null hypothesis: Equal expression from both alleles (p=0.5)
* Alternative: Allelic imbalance (p ≠ 0.5)
* FDR correction for multiple testing

Dispersion Parameter
~~~~~~~~~~~~~~~~~~~~

Two models:
1. **Single**: One dispersion parameter for all regions
2. **Linear**: Dispersion varies with read depth

CLI Usage
---------

Basic Analysis
~~~~~~~~~~~~~~

.. code-block:: bash

   wasp2-analyze find-imbalance counts.tsv

Options
~~~~~~~

.. code-block:: bash

   wasp2-analyze find-imbalance \
     counts.tsv \
     --min-count 10 \
     --pseudocount 1 \
     --model single \
     --output results.tsv

Parameters
----------

``--min-count``
~~~~~~~~~~~~~~~

Minimum total read count per region (default: 10):

.. code-block:: bash

   --min-count 20  # More stringent

``--pseudocount``
~~~~~~~~~~~~~~~~~

Pseudocount added to avoid zero counts (default: 1):

.. code-block:: bash

   --pseudocount 0  # No pseudocount

``--model``
~~~~~~~~~~~

Dispersion model (default: single):

.. code-block:: bash

   --model linear  # Depth-dependent dispersion

``--phased``
~~~~~~~~~~~~

Use phased genotype information:

.. code-block:: bash

   --phased  # Requires phased VCF

Output Format
-------------

Tab-separated file with columns:

Statistical Columns
~~~~~~~~~~~~~~~~~~~

* ``region``: Genomic region identifier
* ``ref_count``: Total reference allele counts
* ``alt_count``: Total alternate allele counts
* ``p_value``: Likelihood ratio test p-value
* ``fdr_pval``: FDR-corrected p-value
* ``effect_size``: Log2 fold-change (ref/alt)

Model Parameters
~~~~~~~~~~~~~~~~

* ``dispersion``: Beta-binomial dispersion parameter
* ``log_likelihood_null``: Null model log-likelihood
* ``log_likelihood_alt``: Alternative model log-likelihood

Interpreting Results
--------------------

Significant Imbalance
~~~~~~~~~~~~~~~~~~~~~

FDR < 0.05 indicates significant imbalance:

* **Biological**: cis-regulatory variation, ASE
* **Technical**: mapping bias (check WASP), PCR artifacts

Effect Size
~~~~~~~~~~~

* log2FC > 1: Strong imbalance (2-fold difference)
* log2FC > 2: Very strong imbalance (4-fold difference)

Single-Cell Analysis
--------------------

For single-cell data:

.. code-block:: bash

   wasp2-analyze find-imbalance-sc \
     adata.h5ad \
     --sample donor1 \
     --groups cell_type \
     --min-count 5

Output: Cell-type-specific imbalance results.

Group Comparison
----------------

Compare imbalance between conditions:

.. code-block:: bash

   wasp2-analyze compare-imbalance \
     adata.h5ad \
     --groups "control,treatment"

Output: Differential imbalance between groups.

Example Workflow
----------------

.. code-block:: bash

   # 1. Count alleles
   wasp2-count count-variants \
     wasp_filtered.bam \
     variants.vcf \
     --region genes.gtf \
     --samples NA12878 \
     --output counts.tsv

   # 2. Analyze imbalance
   wasp2-analyze find-imbalance \
     counts.tsv \
     --min-count 20 \
     --model single \
     --output imbalance.tsv

   # 3. Filter significant results
   awk '$5 < 0.05' imbalance.tsv > significant.tsv

Best Practices
--------------

Read Depth
~~~~~~~~~~

* Minimum 10 reads per region (use ``--min-count``)
* Higher depth = more power
* Consider downsampling very deep regions

Quality Control
~~~~~~~~~~~~~~~

* Use WASP-filtered reads
* Remove low-complexity regions
* Filter low-quality SNPs

Multiple Testing
~~~~~~~~~~~~~~~~

* FDR correction is automatic
* Consider Bonferroni for very important regions
* Validate top hits experimentally

Common Issues
-------------

No Significant Results
~~~~~~~~~~~~~~~~~~~~~~

* Increase sample size
* Check read depth (use deeper sequencing)
* Verify heterozygous SNPs present

Many Significant Results
~~~~~~~~~~~~~~~~~~~~~~~~

* Check for batch effects
* Verify WASP filtering was applied
* Consider stricter FDR threshold

Next Steps
----------

* Validate results with qPCR or DNA-seq
* Integrate with eQTL data
* Perform pathway enrichment analysis
