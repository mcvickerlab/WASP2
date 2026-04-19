False Discovery Rate Correction
================================

WASP2 corrects for multiple testing using the Benjamini-Hochberg (BH) procedure
[Benjamini1995]_, applied via :func:`scipy.stats.false_discovery_control`.

.. code-block:: python

   from scipy.stats import false_discovery_control

   fdr_pvals = false_discovery_control(pvals, method='bh')

Each region receives both a raw p-value (``pval``) from the likelihood ratio test
and a BH-adjusted p-value (``fdr_pval``). Significance is typically declared at
``fdr_pval < 0.05``.

BH controls FDR under independence or positive regression dependence (PRDS).
The chi-squared p-values produced by WASP2's LRT satisfy these assumptions when
regions are analyzed independently. If regions share SNPs or otherwise induce
negative dependence, the more conservative Benjamini-Yekutieli procedure
(``method='by'``) can be substituted.

.. warning::

   ``scipy.stats.false_discovery_control`` raises on NaN p-values, but a manual
   BH implementation based on ``np.minimum.accumulate`` silently propagates NaN
   through the cumulative minimum. Always drop or impute NaN p-values before BH
   correction.

For studies requiring higher-power FDR control with an estimated proportion of
true nulls, :math:`\pi_0`, see Storey's q-value [Storey2003]_. WASP2 does not
apply this by default because BH is sufficient for typical ASE workloads and
requires no tuning.

Reporting
---------

When reporting results, state the correction method, threshold, and number
of tests:

    "BH correction was applied to N = ... chi-squared p-values; regions with
    ``fdr_pval < 0.05`` were declared significant."

Output columns
--------------

.. table::
   :widths: 20 20 60

   ========= ======== ====================================================
   Column    Type     Description
   ========= ======== ====================================================
   pval      float    Raw p-value from likelihood ratio test
   fdr_pval  float    BH-adjusted p-value
   ========= ======== ====================================================

References
----------

.. [Benjamini1995] Benjamini Y, Hochberg Y (1995). Controlling the false
   discovery rate: A practical and powerful approach to multiple testing.
   *Journal of the Royal Statistical Society B* 57:289-300.

.. [Storey2003] Storey JD, Tibshirani R (2003). Statistical significance for
   genomewide studies. *Proceedings of the National Academy of Sciences*
   100:9440-9445.
