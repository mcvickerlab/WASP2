Dispersion Parameter Estimation
================================

The beta-binomial overdispersion parameter :math:`\rho \in (0, 1)` quantifies
variance inflation beyond the binomial expectation:

.. math::

   \text{Var}(X) = N\mu(1-\mu)\left[1 + (N-1)\rho\right]

Dispersion arises from biological variation, PCR amplification, and
aggregation of counts across SNPs within a region. Estimation approaches
follow the beta-binomial ASE literature [Kumasaka2016]_; analogous methods
for negative-binomial dispersion in RNA-seq counts are reviewed in
[Robinson2010]_ and [Yu2013]_.

Single Dispersion Model (default)
---------------------------------

WASP2's default ``single`` model estimates one pooled :math:`\rho` under the
null hypothesis :math:`\mu = 0.5` by maximum likelihood:

.. math::

   \hat{\rho}_{\text{MLE}} = \arg\max_{\rho} \sum_{i=1}^{n}
   \log \text{BetaBinom}(X_i; N_i, \alpha(\rho), \beta(\rho))

with :math:`\alpha = \beta = 0.5(1-\rho)/\rho`.

.. code-block:: python

   from scipy.optimize import minimize_scalar
   from scipy.stats import betabinom
   import numpy as np

   def neg_log_likelihood(rho, ref_counts, n_counts):
       alpha = 0.5 * (1 - rho) / rho
       beta = 0.5 * (1 - rho) / rho
       return -np.sum(betabinom.logpmf(ref_counts, n_counts, alpha, beta))

   rho_mle = minimize_scalar(
       neg_log_likelihood,
       args=(ref_array, n_array),
       method='bounded',
       bounds=(1e-6, 1 - 1e-6),
   ).x

:math:`\rho` is held at this null-model MLE when computing the likelihood-ratio
test statistic (see :doc:`statistical_models`).

Linear Dispersion Model
-----------------------

For large datasets with wide coverage ranges, WASP2 provides a ``linear`` model
in which dispersion varies with total count on the logit scale:

.. math::

   \text{logit}(\rho_i) = \beta_0 + \beta_1 \cdot N_i

Low-coverage regions typically show greater apparent dispersion from sampling
noise; high-coverage regions reveal more of the true biological variance. The
linear model captures this trend with one additional parameter.

.. code-block:: python

   from scipy.optimize import minimize
   from scipy.special import expit

   def neg_ll_linear(params, ref_counts, n_counts):
       beta0, beta1 = params
       rho = expit(np.clip(beta0 + beta1 * n_counts, -10, 10))
       alpha = 0.5 * (1 - rho) / rho
       beta = 0.5 * (1 - rho) / rho
       return -np.sum(betabinom.logpmf(ref_counts, n_counts, alpha, beta))

   beta0, beta1 = minimize(neg_ll_linear, x0=(0, 0),
                           method='Nelder-Mead',
                           args=(ref_array, n_array)).x

Choosing a model
----------------

.. table::
   :widths: 30 70

   ================================ =============================================================
   Context                          Recommended model
   ================================ =============================================================
   Small dataset (< 1,000 regions)  ``single`` — linear model is underpowered
   Uniform coverage                 ``single``
   Large dataset, variable coverage ``linear``
   Model-selection check            Compare AIC/BIC; linear has one extra parameter
   ================================ =============================================================

Model comparison via AIC (:math:`\text{AIC} = 2k - 2\ell(\hat\theta)`) can be
used to check whether the extra parameter is justified on a given dataset.

Convergence
-----------

The MLE optimizer may have trouble when :math:`\rho` is very close to 0 or 1,
or when extreme outliers dominate the pooled likelihood. WASP2 bounds
:math:`\rho \in (10^{-6}, 1-10^{-6})` and clips logit values in the linear
model. Persistent convergence failures usually indicate data issues (CNVs,
mapping artifacts) rather than optimizer pathology.

See Also
--------

- :doc:`statistical_models` — how :math:`\rho` enters the LRT
- :doc:`fdr_correction` — multiple-testing correction downstream

References
----------

.. [Kumasaka2016] Kumasaka N, Knights AJ, Gaffney DJ (2016). Fine-mapping
   cellular QTLs with RASQUAL and ATAC-seq. *Nature Genetics* 48:206-213.

.. [Robinson2010] Robinson MD, Smyth GK (2010). Small-sample estimation of
   negative binomial dispersion, with applications to SAGE data.
   *Biostatistics* 9:321-332. *Analogous NB-dispersion literature.*

.. [Yu2013] Yu D, Huber W, Vitek O (2013). Shrinkage estimation of dispersion
   in Negative Binomial models for RNA-seq experiments with small sample size.
   *Bioinformatics* 29:1275-1282. *Analogous NB-dispersion literature.*
