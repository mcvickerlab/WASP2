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
per comparison.

Notes
-----

* If your count file contains genotype columns for multiple samples, you must
  provide ``--sample`` for single-cell analysis.
* For bulk analysis, region columns are auto-detected when present in the count
  TSV. Use ``--region_col`` only when you need to override that behavior.

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
* :doc:`/tutorials/comparative_imbalance` for group-comparison workflows
