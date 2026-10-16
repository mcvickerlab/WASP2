Analysis Module
===============

Overview
--------

``wasp2-analyze`` runs allelic imbalance analysis on bulk or single-cell count
outputs. The Python CLI calls Rust for bulk inference. Single-cell imbalance
and differential inference currently run in Python/SciPy; single-cell counting
uses Rust.

Numerical compatibility retains the archived statistical procedures with wide
counts and consistent precision. Untouched legacy file readers and explicitly
widened references are evaluated separately; the latter do not claim to reproduce
historical UInt16 overflow or mixed scalar/array precision.

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
* ``--groupby`` (aliases ``--group``, ``--parent_col``, ``--parent``): group on
  an alternate annotation column, such as a parent gene column. The Rust backend
  groups SNP observations by this column before deduplication and fitting;
  this is equivalent to ``--region_col <parent_column>``.

Donor-Local Bulk Analysis
-------------------------

Locked cohort count bundles use ``donor_id`` as the donor identity. Analyze the
bundle directory, its ``count_manifest.json``, or its locked ``counts.tsv.gz``:

.. code-block:: bash

   wasp2-analyze find-imbalance cohort_counts \
     --scope per-donor \
     --unit snv \
     --model single \
     --dispersion-scope per-donor \
     --min-donor-observations 50 \
     --output donor_snv.tsv

The supported units are ``snv`` and ``feature``; ``peak`` is an alias for
``feature``. SNV analysis is unphased. Feature analysis may use ``--phased`` and
uses the locked ``region`` column, or an explicit ``--region-col`` for a
standalone count table. Overlapping feature memberships are retained for feature
tests and dispersion fitting, matching archived Python's retained SNP-feature
observations. Deduplication follows Python's chromosome, position, counts, and
feature key, including GT only when phasing is effective. Different count
observations at the same SNP are retained. SNV output aggregates those retained
observations into one result per position, as archived Python does. Allele and
SNV-ID consistency checks remain enforced.

``--dispersion-scope per-donor`` fits ``single`` or ``linear`` nuisance
parameters independently for each included donor. ``--dispersion-scope global``
fits the combined eligible donor-SNP-feature observations and reuses the returned fit ID and
parameters unchanged in each donor run. In both modes, effect inference and
Benjamini-Hochberg correction are independent within each donor.

A feature-level ``--phased`` request follows archived Python's fallback checks
within each donor after depth filtering. Missing GT, a single distinct GT, or
no numeric ``0|1``/``1|0`` genotype selects unphased inference. Allele-string GT
such as ``A|G`` is retained rather than converted to numeric phase. Provenance
``donors.phasing`` records each included donor's ``requested_phased``,
``effective_phased``, and ``fallback_reason`` before phase-aware deduplication;
``analysis.phased`` remains the
original request. A shared dispersion fit does not share phasing decisions.

QC ``eligible_snv_observations`` and dispersion-sidecar
``donor_eligible_snv_observations`` count unique eligible SNPs per donor, so
multiple memberships do not inflate the donor eligibility threshold.
``fit_observations`` counts actual retained fitting rows: all included donors
for a global fit, or the individual donor for a per-donor fit.

The ``linear`` model carries four inseparable nuisance values:
``linear_d1``, ``linear_d2``, ``linear_depth_center``, and
``linear_depth_scale``.

Linear fitting follows archived Python's balanced-null objective on
pseudocounted depth, with Nelder-Mead starting at ``(0, 0)`` and SciPy's default
simplex, tolerances, and evaluation limits. The fitted coefficients are then
converted to the standardized raw-depth representation above. Equal-depth
inputs use a unit scale, preserving their fitted predictions. The ``[-10, 10]``
predictor guard applies during fitting only; testing uses the unclipped sigmoid
prediction with the existing numerical rho boundary guard.

The fit ID includes all four, and fixed-nuisance donor
runs must echo all four exactly. The dispersion TSV and provenance JSON retain
these values. They retain optimizer or fit status fields when Rust supplies
them; otherwise ``fit_status=parameters_returned`` records only that a usable
parameter payload was returned and does not claim optimizer convergence.

The donor-local route defaults to pseudocount 1 and creates result,
``.dispersion.tsv``, ``.qc.tsv``, and ``.provenance.json`` artifacts without
overwriting existing files. ``--expected-manifest-sha256`` can pin the locked
bundle manifest to an external digest. Standalone legacy tables may use
``sample`` instead of ``donor_id``; a locked bundle may not, and a table with
both columns is rejected. If the count table has ``snv_id``, its exact text is
retained after enforcing a one-to-one mapping with donor, chromosome, position,
reference, and alternate allele. Only tables without the column receive a
synthesized ``chrom:pos:ref:alt`` identifier, so canonical zero-based input IDs
are not rewritten.

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
     --groups B_cell --groups T_cell \
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

Both routes retain the archived chi-square LRT p-values (one degree of freedom).
Unphased feature likelihoods anchor the first SNP to its reference allele and
average the two orientations only for subsequent SNPs. Input order therefore
sets the reference orientation; the first SNP is not symmetrized.

- ``pseudocount`` defaults to **1** on both the legacy route and with
  ``--scope per-donor``. It is added to both ``ref_count`` and
  ``alt_count`` — Laplace-style shrinkage toward :math:`\mu = 0.5`,
  before likelihood fitting; raw output counts exclude these pseudocounts.
- ``min_count`` defaults to **10**; regions with :math:`N < 10` are dropped.
- Beta-binomial LRT: :math:`\rho` is held at its **null-model MLE** while
  maximizing the alternative likelihood over :math:`\mu` (profile likelihood,
  df = 1). :math:`\rho` is not jointly re-estimated under :math:`H_1`.
- Scalar dispersion is fit on :math:`(0,1)`; exact endpoints receive numerical
  regularization. The linear fit uses pseudocounted depth and a :math:`\pm 10`
  predictor guard during fitting. Stored coefficients preserve that predictor
  using :math:`\rho(N)=\operatorname{expit}(d_1+d_2(N-c)/s)` with raw depth
  :math:`N`, center :math:`c=-2\,\mathrm{pseudocount}` and scale :math:`s=1`.
- Both backends use Benjamini–Hochberg FDR correction. Bulk uses the Rust
  implementation; single-cell analysis uses
  :func:`scipy.stats.false_discovery_control` with ``method='bh'``.

See the upstream mapping-filter docs (:doc:`/methods/mapping_filter`) for the
WASP2 canonical filter contract that applies to both remapping and counting
steps.

Outputs
-------

Typical bulk outputs include:

* region or feature identifier
* aggregated ``ref_count`` and ``alt_count``
* p-values and FDR-adjusted p-values

Single-cell outputs contain ``region``, ``num_snps``, ``mu``, ``null_ll``,
``alt_ll``, ``pval``, and ``fdr_pval`` in one file per group. They do not
include aggregated count columns. Differential outputs replace ``mu`` with
``combined_mu``, ``mu1``, and ``mu2``.

Next Steps
----------

* :doc:`counting` to generate bulk or single-cell counts
