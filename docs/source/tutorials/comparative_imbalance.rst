Comparative Imbalance Analysis
==============================

Compare allelic imbalance between biological groups — cell types, conditions,
sex, treatment status — using ``wasp2-analyze compare-imbalance``.

What this tests
---------------

Standard allelic-imbalance analysis asks whether a region shows a preference
for one allele. Comparative analysis asks whether the *degree* of imbalance
differs between groups:

.. code-block:: text

   H0:  Both groups share the same allelic imbalance (mu_combined)
   H1:  Groups have different imbalance (mu_1 != mu_2)

   Test statistic:  LRT = -2 * (log L_null - log L_alt)
   P-value:         P(chi^2(df=1) > LRT)

See :doc:`/methods/statistical_models` for the full model and
:doc:`/user_guide/single_cell` for input-data formats.

Inputs
------

Two files are required:

1. **AnnData count matrix** (``.h5ad``) with ``ref`` and ``alt`` layers, a
   genotype column in ``.obs`` for your sample, and a ``group`` column in
   ``.var`` identifying the per-cell group assignment. Produced by
   ``wasp2-count count-variants-sc``.

2. **Barcode-to-group TSV** — two columns, no header, tab-separated.

   .. code-block:: text

      AAACGAACAGTCAGTT-1    excitatory_neurons
      AAACGAAGTCGCTCTA-1    inhibitory_neurons
      AAAGGATCATCGATGT-1    astrocytes

   Barcodes must match those in the count matrix exactly (including any
   ``-1`` suffix). For exporting barcode maps from Seurat or Scanpy, see
   :doc:`/user_guide/single_cell`.

Example: comparing cell types
-----------------------------

With counts in ``allele_counts.h5ad`` and cell-type assignments in
``barcode_celltype_map.tsv``, run per-group imbalance first, then compare:

.. code-block:: bash

   # Per-group imbalance within each cell type
   wasp2-analyze find-imbalance-sc \
     allele_counts.h5ad \
     barcode_celltype_map.tsv \
     --sample SAMPLE_ID --phased --min 10 -z 3

   # Pairwise comparison between two named groups
   wasp2-analyze compare-imbalance \
     allele_counts.h5ad \
     barcode_celltype_map.tsv \
     --sample SAMPLE_ID \
     --groups "excitatory_neurons,inhibitory_neurons" \
     --phased --min 15

   # Omit --groups to compare all cell types pairwise
   wasp2-analyze compare-imbalance \
     allele_counts.h5ad \
     barcode_celltype_map.tsv \
     --sample SAMPLE_ID --phased --min 15

Pairwise output files are named ``ai_results_<group1>_<group2>.tsv``.

Other common comparisons — sex, treatment, developmental stage — use the
same command with a different ``barcode_map.tsv``. The biology varies; the
command does not.

Output columns
--------------

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - Column
     - Description
   * - ``region``
     - Genomic region identifier
   * - ``num_snps``
     - Shared heterozygous SNPs used for the comparison
   * - ``combined_mu``
     - Reference-allele frequency under H0 (pooled across groups)
   * - ``mu1``, ``mu2``
     - Per-group reference-allele frequencies
   * - ``null_ll``, ``alt_ll``
     - Log-likelihoods under H0 / H1
   * - ``pval``
     - Likelihood-ratio-test p-value
   * - ``fdr_pval``
     - Benjamini–Hochberg adjusted p-value

Interpretation: ``|mu1 - mu2| > 0.1`` is a meaningful difference in allele
preference; ``fdr_pval < 0.05`` declares significance after multiple-testing
correction.

Visualization
-------------

A minimal volcano plot:

.. code-block:: python

   import pandas as pd
   import numpy as np
   import matplotlib.pyplot as plt

   results = pd.read_csv('ai_results_excitatory_inhibitory.tsv', sep='\t')
   results['effect'] = results['mu1'] - results['mu2']
   results['neglog10p'] = -np.log10(results['pval'].clip(lower=1e-300))

   fig, ax = plt.subplots(figsize=(8, 6))
   sig = results['fdr_pval'] < 0.05
   ax.scatter(results.loc[~sig, 'effect'], results.loc[~sig, 'neglog10p'],
              c='lightgrey', s=8, label='FDR >= 0.05')
   ax.scatter(results.loc[sig, 'effect'], results.loc[sig, 'neglog10p'],
              c='crimson', s=14, label='FDR < 0.05')
   ax.axhline(-np.log10(0.05), ls='--', c='black', alpha=0.4)
   ax.axvline(0, ls='-', c='black', alpha=0.3)
   ax.set(xlabel='mu1 - mu2', ylabel='-log10(p)',
          title='Excitatory vs. inhibitory neurons')
   ax.legend()
   plt.tight_layout()

For heatmaps across many cell types or overlap with eQTL/cis-QTL tables, see
the :doc:`/user_guide/analysis` guide.

Command-line reference
----------------------

.. code-block:: text

   Usage: wasp2-analyze compare-imbalance [OPTIONS] ADATA BARCODE_MAP

     ADATA        AnnData file with allele counts (.h5ad)
     BARCODE_MAP  TSV file mapping barcodes to groups

     --groups TEXT        Comma-separated groups (default: all pairwise)
     --min INT            Minimum allele count per region per group [10]
     --pseudocount INT    Pseudocount for zero counts [1]
     --sample TEXT        Sample name for genotype filtering
     --phased             Use phased genotype information
     -z, --z_cutoff INT   Z-score outlier filter
     --out_file TEXT      Output path

Good practices
--------------

- Run on WASP-filtered BAMs to remove mapping-bias artifacts
  (:doc:`/methods/mapping_filter`).
- Use ``--min 15`` or higher for robust estimates when per-group coverage
  is moderate.
- Apply ``-z 3`` to drop SNP outliers from CNVs or mapping artifacts.
- Check that groups have comparable sequencing depth; very uneven depth can
  produce false-positive differences.
- Validate top hits in an independent cohort or with an orthogonal assay.

Common issues
-------------

**Few significant results.** Increase coverage, merge similar populations to
increase per-group counts, or use phased genotypes if available.

**Too many significant results.** Check for batch effects between groups,
confirm WASP filtering was applied, and use a stricter FDR threshold.

**Memory.** For large cohorts, process chromosomes separately and concatenate
results.

See Also
--------

- :doc:`/user_guide/analysis` — analysis-CLI reference and parameters
- :doc:`/user_guide/single_cell` — input data formats, barcode exports
- :doc:`/tutorials/scrna_seq` — basic single-cell workflow
- :doc:`/methods/statistical_models` — the LRT underlying this test
