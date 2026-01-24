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

For single-cell data, WASP2 detects allelic imbalance within specific cell populations
using aggregated counts across cells of the same type.

.. code-block:: bash

   wasp2-analyze find-imbalance-sc \
     adata.h5ad \
     barcode_map.tsv \
     --sample donor1 \
     --min-count 10

Output: Per-celltype TSV files (``ai_results_[CELLTYPE].tsv``).

Single-Cell Comparative Imbalance
---------------------------------

Overview
~~~~~~~~

The comparative imbalance analysis detects **differential allelic imbalance** between
cell types, conditions, or biological groups. This is useful for identifying:

* Cell-type-specific regulatory variation
* Sex differences in chromatin accessibility
* Condition-dependent allelic effects (e.g., treatment vs control)
* Developmental stage-specific imbalance

Statistical Model
~~~~~~~~~~~~~~~~~

The analysis uses a **likelihood ratio test (LRT)** comparing two hypotheses:

* **Null (H0)**: Both groups share the same allelic imbalance (μ_combined)
* **Alternative (H1)**: Groups have different imbalance (μ₁ ≠ μ₂)

The test statistic follows a chi-squared distribution with 1 degree of freedom:

.. code-block:: text

   LRT = -2 × (log L_null - log L_alt)
   p-value = P(χ²(df=1) > LRT)

Input Format: Count Matrix (.h5ad)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The count matrix must be an AnnData object with the following structure:

.. code-block:: text

   AnnData object (n_obs × n_vars)
   ├── .obs                 # SNP metadata (rows)
   │   ├── index            # SNP identifiers (0, 1, 2, ...)
   │   └── [sample_name]    # Genotypes: '0|1', '1|0', '0/1', '1/0'
   │
   ├── .var                 # Cell metadata (columns)
   │   └── group            # Cell type/group assignment
   │
   ├── .layers
   │   ├── "ref"            # Reference allele counts (sparse matrix)
   │   └── "alt"            # Alternate allele counts (sparse matrix)
   │
   └── .uns
       ├── feature          # DataFrame: SNP → region mapping
       └── samples          # List of sample names

**Example count matrix creation:**

.. code-block:: bash

   # Generate counts from BAM + variants + barcodes
   wasp2-count count-variants-sc \
     aligned.bam \
     variants.vcf.gz \
     barcodes.txt \
     --samples NA12878 \
     --feature peaks.bed \
     --out_file allele_counts.h5ad

Input Format: Barcode Map (TSV)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A two-column TSV file (no header) mapping cell barcodes to groups:

.. code-block:: text

   AAACGAACAGTCAGTT-1    excitatory_neurons
   AAACGAAGTCGCTCTA-1    inhibitory_neurons
   AAACGAAGTGAACCTA-1    excitatory_neurons
   AAAGGATCATCGATGT-1    astrocytes
   AAAGGATGTGCAACGA-1    microglia

**Requirements:**

* Tab-separated (``\t``)
* No header row
* Barcodes must match those in the count matrix
* Groups can be cell types, conditions, sex, or any categorical variable

Basic Usage
~~~~~~~~~~~

Compare imbalance between all groups:

.. code-block:: bash

   wasp2-analyze compare-imbalance \
     allele_counts.h5ad \
     barcode_map.tsv

Compare specific groups only:

.. code-block:: bash

   wasp2-analyze compare-imbalance \
     allele_counts.h5ad \
     barcode_map.tsv \
     --groups "excitatory_neurons,inhibitory_neurons"

Output Format
~~~~~~~~~~~~~

Results are written to ``ai_results_[GROUP1]_[GROUP2].tsv``:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Column
     - Description
   * - ``region``
     - Genomic region identifier (peak or gene)
   * - ``num_snps``
     - Number of shared heterozygous SNPs in region
   * - ``combined_mu``
     - Reference allele frequency under null hypothesis (shared)
   * - ``mu1``
     - Reference allele frequency in group 1
   * - ``mu2``
     - Reference allele frequency in group 2
   * - ``null_ll``
     - Log-likelihood under null (shared μ)
   * - ``alt_ll``
     - Log-likelihood under alternative (separate μ values)
   * - ``pval``
     - Likelihood ratio test p-value
   * - ``fdr_pval``
     - FDR-corrected p-value (Benjamini-Hochberg)

**Interpreting results:**

* ``fdr_pval < 0.05``: Significant differential imbalance
* ``|mu1 - mu2| > 0.1``: Meaningful effect size (~20% difference)
* ``mu < 0.5``: Alternate allele favored; ``mu > 0.5``: Reference allele favored

Parameters
~~~~~~~~~~

``--groups``
   Comma-separated list of groups to compare. If omitted, compares all pairwise
   combinations found in the barcode map.

``--min``
   Minimum total allele count per region per group (default: 10). Higher values
   increase statistical power but reduce the number of testable regions.

``--pseudocount``
   Pseudocount added to avoid zero counts (default: 1). Affects dispersion estimation.

``--sample``
   Sample name for heterozygous SNP filtering. Required if multiple samples are
   present in the count matrix.

``--phased``
   Use phased genotype information from VCF. Requires genotypes in ``0|1`` or ``1|0``
   format. Improves power when haplotype phase is known.

``-z/--z_cutoff``
   Remove SNPs with counts exceeding this z-score threshold. Useful for removing
   outliers caused by mapping artifacts or copy number variation.

Example: Sex Differences Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Identify chromatin accessibility regions with sex-biased allelic imbalance:

**Step 1: Prepare barcode map with sex labels**

.. code-block:: text

   # barcode_sex_map.tsv
   AAACGAACAGTCAGTT-1    male
   AAACGAAGTCGCTCTA-1    female
   AAACGAAGTGAACCTA-1    male
   AAAGGATCATCGATGT-1    female

**Step 2: Run comparative analysis**

.. code-block:: bash

   wasp2-analyze compare-imbalance \
     allele_counts.h5ad \
     barcode_sex_map.tsv \
     --groups "male,female" \
     --min 20 \
     --out_file ai_results_sex_comparison.tsv

**Step 3: Filter significant results**

.. code-block:: bash

   # Extract regions with significant sex differences
   awk -F'\t' 'NR==1 || $9 < 0.05' ai_results_male_female.tsv > significant_sex_diff.tsv

   # Find regions with large effect size
   awk -F'\t' 'NR==1 || ($4 - $5 > 0.15 || $5 - $4 > 0.15)' significant_sex_diff.tsv

Example: snATAC-seq Cell Type Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Complete workflow for analyzing cell-type-specific chromatin accessibility imbalance:

**Step 1: Count alleles from snATAC-seq BAM**

.. code-block:: bash

   # Extract valid barcodes from Cell Ranger output
   zcat filtered_barcodes.tsv.gz > barcodes.txt

   # Count alleles at heterozygous SNPs overlapping peaks
   wasp2-count count-variants-sc \
     possorted_bam.bam \
     phased_variants.vcf.gz \
     barcodes.txt \
     --samples sample1 \
     --feature atac_peaks.bed \
     --out_file snATAC_counts.h5ad

**Step 2: Create barcode-to-celltype mapping**

Export cell type annotations from your clustering analysis (e.g., Seurat, ArchR):

.. code-block:: r

   # R/Seurat example
   write.table(
     data.frame(barcode = Cells(seurat_obj),
                celltype = Idents(seurat_obj)),
     "barcode_celltype_map.tsv",
     sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
   )

**Step 3: Run single-cell imbalance analysis**

.. code-block:: bash

   # Per-celltype analysis
   wasp2-analyze find-imbalance-sc \
     snATAC_counts.h5ad \
     barcode_celltype_map.tsv \
     --sample sample1 \
     --phased \
     --min 10 \
     -z 3.0

**Step 4: Compare imbalance between cell types**

.. code-block:: bash

   # Compare specific cell types
   wasp2-analyze compare-imbalance \
     snATAC_counts.h5ad \
     barcode_celltype_map.tsv \
     --sample sample1 \
     --groups "excitatory,inhibitory,astrocyte" \
     --phased \
     --min 15

   # This produces:
   # - ai_results_excitatory_inhibitory.tsv
   # - ai_results_excitatory_astrocyte.tsv
   # - ai_results_inhibitory_astrocyte.tsv

**Step 5: Identify cell-type-specific regulatory regions**

.. code-block:: bash

   # Find peaks with differential imbalance between excitatory and inhibitory neurons
   awk -F'\t' '$9 < 0.01 && ($4 > 0.6 || $4 < 0.4)' \
     ai_results_excitatory_inhibitory.tsv > neuron_subtype_specific_peaks.tsv

Best Practices for Single-Cell Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Data Quality:**

* Use WASP-filtered BAM files to remove mapping bias
* Require ≥10 total counts per region per group (``--min 10``)
* Apply z-score filtering to remove outliers (``-z 3.0``)

**Statistical Power:**

* Merge similar cell types if individual populations have low coverage
* Use phased genotypes when available (``--phased``)
* Focus on regions with multiple SNPs for better estimates

**Interpretation:**

* Consider biological replication across samples
* Validate top hits with orthogonal methods (allele-specific CRISPR, etc.)
* Integrate with eQTL data to identify causal variants

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
