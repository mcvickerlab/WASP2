Analysis Module
===============

Overview
--------

``wasp2-analyze`` runs allelic imbalance analysis on bulk or single-cell count
outputs.

Commands:

* ``find-imbalance``: bulk allelic imbalance from TSV counts
* ``find-imbalance-sc``: per-group single-cell imbalance from ``.h5ad`` counts
* ``compare-imbalance``: differential imbalance between single-cell groups

Bulk Analysis
-------------

.. code-block:: bash

   wasp2-analyze find-imbalance \
     counts.tsv \
     --min 10 \
     --pseudocount 1 \
     --output ai_results.tsv

Useful options:

* ``--min`` / ``--min_count``: minimum total count threshold
* ``--pseudocount``: pseudocount added before modeling
* ``--model``: dispersion model (currently ``single`` or ``linear`` input)
* ``--output`` / ``--out_file`` / ``-o``: output TSV path
* ``--region_col``: explicit region column name if auto-detection is not desired
* ``--groupby``: group on an alternate annotation column, such as a parent gene column

Single-Cell Analysis
--------------------

.. code-block:: bash

   wasp2-analyze find-imbalance-sc \
     allele_counts.h5ad \
     barcode_groups.tsv \
     --sample SAMPLE1 \
     --min 10 \
     --out_file ai_results.tsv

``barcode_groups.tsv`` is a two-column TSV:

.. code-block:: text

   BARCODE<TAB>GROUP

The command writes one output file per group using the requested output prefix.

Comparative Single-Cell Analysis
--------------------------------

.. code-block:: bash

   wasp2-analyze compare-imbalance \
     allele_counts.h5ad \
     barcode_groups.tsv \
     --sample SAMPLE1 \
     --groups B_cell T_cell \
     --out_file compare_ai.tsv

This compares allelic imbalance between the requested groups and writes one TSV
per comparison. Each output has per-group columns ``mu1`` and ``mu2`` alongside
the pooled-null ``combined_mu``, plus ``null_ll`` / ``alt_ll`` / ``pval`` /
``fdr_pval``.

Notes
-----

* If your count file contains genotype columns for multiple samples, you must
  provide ``--sample`` for single-cell analysis.
* For bulk analysis, region columns are auto-detected when present in the count
  TSV. Use ``--region_col`` only when you need to override that behavior.

Defaults and contracts
----------------------

- ``pseudocount`` defaults to **1** (added to both ``ref_count`` and
  ``alt_count`` — Laplace-style shrinkage toward :math:`\mu = 0.5`,
  conservative under :math:`H_0`).
- ``min_count`` defaults to **10**; regions with :math:`N < 10` are dropped.
- Beta-binomial LRT: :math:`\rho` is held at its **null-model MLE** while
  maximizing the alternative likelihood over :math:`\mu` (profile likelihood,
  df = 1). :math:`\rho` is not jointly re-estimated under :math:`H_1`.
- Dispersion optimizer bounds :math:`\rho \in (10^{-6},\, 1-10^{-6})`; the
  linear-dispersion model clips the logit at :math:`\pm 10` for numerical
  stability on extreme :math:`N`.
- FDR correction uses :func:`scipy.stats.false_discovery_control` with
  ``method='bh'`` (Benjamini–Hochberg). ``scipy`` raises on NaN p-values,
  but a hand-written BH loop using ``np.minimum.accumulate`` silently
  propagates NaN through the cumulative minimum — always drop or impute
  NaN p-values before BH correction.

See the upstream mapping-filter docs (:doc:`/methods/mapping_filter`) for the
WASP2 canonical filter contract that applies to both remapping and counting
steps.

Outputs
-------

Typical bulk outputs include:

* region or feature identifier
* aggregated ``ref_count`` and ``alt_count``
* p-values and FDR-adjusted p-values

Typical single-cell outputs include the same statistics stratified by barcode
group.

Next Steps
----------

* :doc:`counting` to generate bulk or single-cell counts
