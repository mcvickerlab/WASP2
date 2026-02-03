RNA-seq Allelic Imbalance Tutorial
===================================

This tutorial demonstrates a complete workflow for detecting allele-specific expression (ASE)
in bulk RNA-seq data using WASP2.

**Estimated time:** ~30 minutes

Overview
--------

The tutorial covers the complete RNA-seq allelic imbalance analysis pipeline:

1. **Data Loading** - BAM, VCF, and gene annotations (GTF)
2. **Allele Counting** - Count reads at heterozygous SNPs within genes/exons
3. **Statistical Testing** - Beta-binomial model for allelic imbalance detection
4. **ASE Visualization** - Volcano plots and allele ratio distributions
5. **Imprinting Detection** - Identify monoallelic expression patterns
6. **eQTL Integration** - Connect ASE signals to regulatory variants

Prerequisites
-------------

**Software:**

* WASP2 (``pip install wasp2``)
* Python packages: pandas, numpy, matplotlib, seaborn, scipy

**Data:**

* Aligned BAM file (coordinate-sorted, indexed)
* Phased VCF file with heterozygous variants
* Gene annotation file (GTF format, e.g., GENCODE)

Workflow Summary
----------------

Step 1: Count Alleles at Genes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use ``wasp2-count count-variants`` to count allele-specific reads at heterozygous SNPs:

.. code-block:: bash

   wasp2-count count-variants \
       sample.bam \
       variants.vcf.gz \
       --samples SAMPLE_ID \
       --region genes.gtf \
       --out_file allele_counts.tsv

This produces a TSV file with columns:

* ``chr``, ``pos``: SNP location
* ``ref``, ``alt``: Alleles
* ``ref_count``, ``alt_count``: Read counts per allele
* ``other_count``: Reads supporting other alleles (non-ref, non-alt)
* ``gene_id``, ``gene_name``: Overlapping gene annotation
* ``feature``: Feature type (exon, intron, etc.) when using GTF

Step 2: Test for Allelic Imbalance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use ``wasp2-analyze find-imbalance`` to detect significant ASE:

.. code-block:: bash

   wasp2-analyze find-imbalance \
       allele_counts.tsv \
       --min 10 \
       --pseudocount 1 \
       --phased \
       --out_file ai_results.tsv

The beta-binomial model tests for deviation from 50:50 allele ratios, accounting for
biological overdispersion.

**Output columns:**

* ``region``: Gene identifier
* ``ref_count``, ``alt_count``: Aggregated counts
* ``p_value``: Likelihood ratio test p-value
* ``fdr_pval``: FDR-corrected p-value
* ``effect_size``: Log2 fold change (ref/alt)

Step 3: Identify Significant ASE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Filter results by significance and effect size:

.. code-block:: bash

   # Significant ASE (FDR < 0.05, |log2FC| > 1)
   # Column indices: $5 = fdr_pval, $6 = effect_size
   awk -F'\t' 'NR==1 || ($5 < 0.05 && ($6 > 1 || $6 < -1))' ai_results.tsv > significant_ase.tsv

Step 4: Detect Imprinting Patterns
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Monoallelic expression (>90:10 allele ratio) may indicate genomic imprinting:

.. code-block:: python

   import pandas as pd

   results = pd.read_csv('ai_results.tsv', sep='\t')
   total = results['ref_count'] + results['alt_count']
   results['ref_ratio'] = results['ref_count'] / total

   # Monoallelic genes
   monoallelic = results[
       (results['ref_ratio'] > 0.9) | (results['ref_ratio'] < 0.1)
   ]

Step 5: Integrate with eQTL Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Cross-reference ASE results with eQTL databases (e.g., GTEx):

.. code-block:: python

   import pandas as pd

   ase = pd.read_csv('ai_results.tsv', sep='\t')
   eqtl = pd.read_csv('gtex_eqtl.tsv', sep='\t')

   # Merge on gene ID
   integrated = ase.merge(eqtl, left_on='region', right_on='gene_id')

   # Check direction concordance between ASE and eQTL
   # ASE effect_size > 0 means REFERENCE allele is more expressed
   # eQTL slope > 0 means ALTERNATE allele INCREASES expression
   # Therefore, concordance = OPPOSITE signs (ref high in ASE = alt low in eQTL)
   integrated['concordant'] = (
       (integrated['effect_size'] > 0) != (integrated['slope'] > 0)
   )

Key Concepts
------------

Beta-Binomial Model
~~~~~~~~~~~~~~~~~~~

WASP2 uses a beta-binomial distribution to model allele counts:

* Accounts for **overdispersion** (biological variation beyond binomial sampling)
* Models **technical noise** from PCR amplification and sequencing
* Aggregates information across **multiple SNPs** per gene

The null hypothesis is equal expression from both alleles (p = 0.5).

Effect Size Interpretation
~~~~~~~~~~~~~~~~~~~~~~~~~~

* **log2FC > 1**: Reference allele 2x more expressed (strong ASE)
* **log2FC > 2**: Reference allele 4x more expressed (very strong ASE)
* **log2FC near 0**: Balanced expression (no ASE)

Significance Thresholds
~~~~~~~~~~~~~~~~~~~~~~~

* **FDR < 0.05**: Standard significance threshold
* **FDR < 0.01**: Stringent threshold for high-confidence hits
* Combine with effect size filters to focus on biologically meaningful results

Troubleshooting
---------------

Low SNP Counts
~~~~~~~~~~~~~~

If few heterozygous SNPs are detected:

* Verify VCF contains heterozygous genotypes:

  - **Phased format**: GT = ``0|1`` or ``1|0`` (pipe separator)
  - **Unphased format**: GT = ``0/1`` or ``1/0`` (slash separator)
  - Use ``--phased`` flag only with phased genotypes

* Check sample ID matches VCF sample column
* Ensure BAM and VCF use the same reference genome

No Significant Results
~~~~~~~~~~~~~~~~~~~~~~

* Increase sequencing depth (more reads = more power)
* Lower ``--min`` threshold (but interpret with caution)
* Check for batch effects or technical artifacts

Too Many Significant Results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Verify WASP mapping bias correction was applied
* Check for copy number variation (CNV) artifacts
* Use stricter FDR threshold (e.g., 0.01)

See Also
--------

* :doc:`/user_guide/counting` - Detailed counting options
* :doc:`/user_guide/analysis` - Statistical methods and parameters
* :doc:`/tutorials/scrna_seq` - Single-cell RNA-seq tutorial
* :doc:`/tutorials/comparative_imbalance` - Comparing ASE between groups
