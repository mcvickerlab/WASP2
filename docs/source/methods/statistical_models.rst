Statistical Model
=================

WASP2 detects allelic imbalance with a **beta-binomial likelihood-ratio
test** applied per region (gene or peak). This page documents the model,
the dispersion-estimation choices, the multiple-testing correction, and
the WASP2-specific defaults and contracts that are not covered by the
papers cited below.

For the WASP mapping filter upstream of counting, see
:doc:`mapping_filter`.

Model
-----

With :math:`X` the reference count out of :math:`N` total reads at a
site:

.. math::

   p \sim \text{Beta}(\alpha, \beta), \qquad X \mid p \sim \text{Binomial}(N, p)

WASP2 uses the mean-dispersion parameterization
:math:`\mu = \alpha/(\alpha+\beta)` and :math:`\rho = 1/(\alpha+\beta+1)`:

.. math::

   \alpha = \mu\cdot\frac{1-\rho}{\rho}, \qquad \beta = (1-\mu)\cdot\frac{1-\rho}{\rho}, \qquad \text{Var}(X) = N\mu(1-\mu)[1 + (N-1)\rho]

Binomial variance is recovered as :math:`\rho \to 0`. The framework
follows the RASQUAL / MBASED lineage ([Kumasaka2016]_, [Mayba2014]_).

Hypothesis test
---------------

For each region, WASP2 tests

.. math::

   H_0: \mu = 0.5 \qquad \text{vs.} \qquad H_1: \mu \neq 0.5

with the likelihood-ratio statistic

.. math::

   \Lambda = -2[\log\mathcal{L}_0 - \log\mathcal{L}_1] \sim \chi^2_1 \text{ under } H_0.

.. important::

   **Profile-likelihood convention.** :math:`\rho` is held at its
   null-model MLE when maximizing :math:`\mathcal{L}_1` over :math:`\mu`
   (profile likelihood in :math:`\mu`, df = 1). WASP2 does **not** jointly
   re-estimate :math:`\rho` under :math:`H_1`.

MLE for :math:`\mu` under :math:`H_1`:

.. code-block:: python

   from scipy.optimize import minimize_scalar
   from scipy.stats import betabinom
   import numpy as np

   def neg_ll_mu(mu, ref, n, rho):
       alpha = mu * (1 - rho) / rho
       beta = (1 - mu) * (1 - rho) / rho
       return -np.sum(betabinom.logpmf(ref, n, alpha, beta))

   mu_mle = minimize_scalar(neg_ll_mu, args=(ref, n, rho),
                            method='bounded', bounds=(0, 1)).x

Dispersion
----------

Estimation of :math:`\rho` is by MLE under :math:`H_0` (:math:`\mu=0.5`).
WASP2 ships two models:

- ``single`` (default) — one pooled :math:`\rho` across all regions.
  Appropriate for small-to-moderate datasets or uniform coverage.
- ``linear`` — :math:`\text{logit}(\rho_i) = \beta_0 + \beta_1 N_i`.
  Appropriate for large cohorts with wide coverage ranges, where
  low-coverage regions show different apparent dispersion from
  high-coverage ones. Compared to ``single`` via AIC/BIC.

**Implementation contract.** :math:`\rho` is bounded to the open interval
:math:`(10^{-6},\, 1-10^{-6})` during MLE to avoid the boundary
singularities of the beta-binomial. In the linear model, the logit is
clipped to :math:`[-10,\, 10]` before the inverse-logit, preserving
numerical stability on extreme :math:`N`.

Phased vs. unphased data
-------------------------

**Phased** genotypes reduce to a direct beta-binomial per SNP, with
:math:`\mu` or :math:`1-\mu` depending on which haplotype carries the
reference allele.

**Unphased** analysis marginalizes over the :math:`2^{n-1}` phase
configurations for a region of :math:`n` SNPs using a uniform prior
(dynamic programming over SNPs; approach of [Mayba2014]_):

.. math::

   P(\mathbf{X}) = \sum_{\phi} P(\mathbf{X}\mid\phi)\,P(\phi),
   \qquad P(\phi) = 0.5^{n-1}.

Multiple-testing correction
---------------------------

WASP2 applies the Benjamini–Hochberg procedure [Benjamini1995]_ via
:func:`scipy.stats.false_discovery_control`:

.. code-block:: python

   from scipy.stats import false_discovery_control
   fdr_pvals = false_discovery_control(pvals, method='bh')

BH controls FDR under independence or positive regression dependence
(PRDS); :math:`\chi^2_1` p-values from the LRT satisfy these assumptions
when regions are analyzed independently. For dependent tests across
overlapping regions, Benjamini–Yekutieli (``method='by'``) is the
conservative fallback.

.. warning::

   ``scipy.stats.false_discovery_control`` raises on NaN p-values, but a
   hand-written BH loop based on ``np.minimum.accumulate`` **silently
   propagates NaN** through the cumulative minimum — every corrected
   p-value at or above the first NaN becomes NaN. Drop or impute NaN
   p-values before BH correction.

Defaults and output
-------------------

WASP2 default parameter values that affect the results:

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Parameter
     - Default
     - Meaning
   * - ``pseudocount``
     - 1
     - Added to both ``ref_count`` and ``alt_count``. Laplace-style
       shrinkage toward :math:`\mu = 0.5` — conservative under :math:`H_0`.
   * - ``min_count``
     - 10
     - Regions with :math:`N < 10` are dropped; below this, power is too
       low and the MLE is unstable.

Output columns per region:

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Column
     - Type
     - Interpretation
   * - ``null_ll``
     - float
     - Log-likelihood under :math:`H_0` (:math:`\mu=0.5`).
   * - ``alt_ll``
     - float
     - Log-likelihood under :math:`H_1` (:math:`\mu=\hat\mu_{\text{MLE}}`).
   * - ``mu``
     - float
     - MLE of the reference-allele proportion.
   * - ``lrt``
     - float
     - Likelihood-ratio statistic.
   * - ``pval``
     - float
     - :math:`\chi^2_1` p-value.
   * - ``fdr_pval``
     - float
     - BH-adjusted p-value.

:math:`\mu > 0.5` indicates a reference-allele preference;
:math:`|\mu - 0.5|` is the effect size.

See Also
--------

- :doc:`mapping_filter` — canonical WASP filter contract, upstream of counts

References
----------

.. [Kumasaka2016] Kumasaka N, Knights AJ, Gaffney DJ (2016). Fine-mapping
   cellular QTLs with RASQUAL and ATAC-seq. *Nature Genetics* 48:206-213.

.. [Mayba2014] Mayba O, Gilbert HN, Liu J, et al. (2014). MBASED:
   allele-specific expression detection in cancer tissues and cell lines.
   *Genome Biology* 15:405.

.. [Benjamini1995] Benjamini Y, Hochberg Y (1995). Controlling the false
   discovery rate: A practical and powerful approach to multiple testing.
   *Journal of the Royal Statistical Society B* 57:289-300.
