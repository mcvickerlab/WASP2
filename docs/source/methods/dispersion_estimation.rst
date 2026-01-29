Dispersion Parameter Estimation
================================

This document describes how WASP2 estimates the overdispersion parameter
:math:`\rho` in the beta-binomial model.

.. contents:: Contents
   :local:
   :depth: 2

The Dispersion Parameter
------------------------

The dispersion parameter :math:`\rho` quantifies the excess variance beyond
the binomial expectation. In the beta-binomial model:

.. math::

   \text{Var}(X) = N\mu(1-\mu)[1 + (N-1)\rho]

where:

- :math:`\rho = 0`: Binomial variance (no overdispersion)
- :math:`\rho > 0`: Variance inflation due to correlated sampling

Estimation Methods
------------------

WASP2 supports two approaches for dispersion estimation:

1. **Maximum Likelihood Estimation (MLE)**: Optimize the likelihood function
2. **Method of Moments (MoM)**: Solve equations based on sample moments

Both methods have trade-offs in terms of efficiency, bias, and computational cost.

Maximum Likelihood Estimation
-----------------------------

MLE finds the value of :math:`\rho` that maximizes the likelihood of the
observed data.

Single Dispersion Model
^^^^^^^^^^^^^^^^^^^^^^^

WASP2's default approach estimates a single :math:`\rho` across all observations:

.. math::

   \hat{\rho}_{\text{MLE}} = \arg\max_{\rho} \sum_{i=1}^{n} \log P(X_i | N_i, \mu=0.5, \rho)

Under the null hypothesis of no imbalance (:math:`\mu = 0.5`), the log-likelihood is:

.. math::

   \ell(\rho) = \sum_{i=1}^{n} \log \text{BetaBinom}(X_i; N_i, \alpha(\rho), \beta(\rho))

where :math:`\alpha = \beta = 0.5(1-\rho)/\rho`.

**Implementation**:

.. code-block:: python

   from scipy.optimize import minimize_scalar
   from scipy.stats import betabinom
   import numpy as np

   def neg_log_likelihood(rho, ref_counts, n_counts):
       alpha = 0.5 * (1 - rho) / rho
       beta = 0.5 * (1 - rho) / rho
       return -np.sum(betabinom.logpmf(ref_counts, n_counts, alpha, beta))

   result = minimize_scalar(
       neg_log_likelihood,
       args=(ref_array, n_array),
       method='bounded',
       bounds=(1e-6, 1 - 1e-6)
   )
   rho_mle = result.x

**Properties of MLE**:

- Asymptotically unbiased as :math:`n \to \infty`
- Achieves the Cramér-Rao lower bound (efficient)
- Computationally requires numerical optimization
- May be sensitive to outliers

Linear Dispersion Model
^^^^^^^^^^^^^^^^^^^^^^^

For large datasets with variable coverage, WASP2 offers a model where
dispersion varies linearly with total count on the logit scale:

.. math::

   \text{logit}(\rho_i) = \beta_0 + \beta_1 \cdot N_i

The logit link ensures :math:`\rho_i \in (0, 1)`:

.. math::

   \rho_i = \frac{\exp(\beta_0 + \beta_1 N_i)}{1 + \exp(\beta_0 + \beta_1 N_i)}

**Motivation**:

Empirically, regions with different coverage levels may exhibit different
dispersion characteristics:

- **Low coverage**: Greater sampling noise, potentially higher apparent dispersion
- **High coverage**: More stable estimates, may reveal true biological variance

**Implementation**:

.. code-block:: python

   from scipy.optimize import minimize
   from scipy.special import expit

   def neg_ll_linear(params, ref_counts, n_counts):
       beta0, beta1 = params
       logit_rho = beta0 + beta1 * n_counts
       # Clip to avoid numerical issues
       logit_rho = np.clip(logit_rho, -10, 10)
       rho = expit(logit_rho)
       alpha = 0.5 * (1 - rho) / rho
       beta = 0.5 * (1 - rho) / rho
       return -np.sum(betabinom.logpmf(ref_counts, n_counts, alpha, beta))

   result = minimize(
       neg_ll_linear,
       x0=(0, 0),
       method='Nelder-Mead',
       args=(ref_array, n_array)
   )
   beta0, beta1 = result.x

Method of Moments
-----------------

MoM estimates :math:`\rho` by equating theoretical and sample moments.

Variance-Based Estimator
^^^^^^^^^^^^^^^^^^^^^^^^

For a beta-binomial with :math:`\mu = 0.5`, the variance is:

.. math::

   \text{Var}(X) = \frac{N}{4}[1 + (N-1)\rho]

Solving for :math:`\rho`:

.. math::

   \hat{\rho}_{\text{MoM}} = \frac{4S^2/N - 1}{N - 1}

where :math:`S^2` is the sample variance of :math:`X/N` (the allelic ratio).

**Pooled Estimator**:

For observations with varying :math:`N`:

.. math::

   \hat{\rho}_{\text{MoM}} = \frac{\sum_i (X_i - N_i/2)^2 - \sum_i N_i/4}{\sum_i N_i(N_i-1)/4}

**Properties of MoM**:

- Closed-form solution (fast computation)
- May produce negative estimates (truncate to 0)
- Less efficient than MLE, especially for small samples
- More robust to model misspecification

Comparison: MLE vs MoM
----------------------

.. table:: MLE vs MoM for Dispersion Estimation
   :widths: 25 37 38

   =================== ============================== ==============================
   Property            MLE                            MoM
   =================== ============================== ==============================
   Computation         Iterative optimization         Closed-form
   Efficiency          Optimal (achieves CRLB)        Suboptimal
   Bias                Asymptotically unbiased        May be biased for small n
   Robustness          Sensitive to outliers          More robust
   Boundary behavior   Always in (0,1)                May give ρ < 0
   WASP2 default       Yes                            No
   =================== ============================== ==============================

WASP2 uses MLE because:

1. The beta-binomial likelihood is well-behaved
2. Modern optimization is fast enough for typical datasets
3. MLE provides consistent estimates across sample sizes

Practical Considerations
------------------------

Convergence Issues
^^^^^^^^^^^^^^^^^^

The MLE optimizer may fail to converge if:

- :math:`\rho` is very close to 0 (nearly binomial data)
- :math:`\rho` is very close to 1 (extreme overdispersion)
- The data contains extreme outliers

WASP2 handles these by:

- Bounding :math:`\rho` away from 0 and 1
- Clipping logit values to avoid overflow
- Using robust optimization methods (bounded, Nelder-Mead)

Sample Size Requirements
^^^^^^^^^^^^^^^^^^^^^^^^

MLE performance depends on sample size:

.. table:: Dispersion Estimate Quality by Sample Size
   :widths: 25 25 50

   ============ ================ =====================================
   n (regions)  CV of estimate   Recommendation
   ============ ================ =====================================
   < 50         > 50%            Use pooled estimate or prior
   50-200       20-50%           MLE reasonable but uncertain
   200-1000     10-20%           MLE reliable
   > 1000       < 10%            MLE highly accurate
   ============ ================ =====================================

Model Selection
^^^^^^^^^^^^^^^

Choosing between single and linear dispersion models:

**Use Single Dispersion When**:

- Dataset is small (< 1000 regions)
- Coverage is relatively uniform
- Quick analysis is needed

**Use Linear Dispersion When**:

- Large dataset (> 10,000 regions)
- Wide range of coverage values
- Systematic coverage-dispersion relationship suspected

Model comparison can be done via AIC/BIC:

.. math::

   \text{AIC} = 2k - 2\ell(\hat{\theta})

where :math:`k` is the number of parameters (1 for single, 2 for linear).

See Also
--------

- :doc:`fdr_correction` - Multiple testing correction after estimation

References
----------

.. [Robinson2010] Robinson MD, Smyth GK (2010). Small-sample estimation of
   negative binomial dispersion, with applications to SAGE data.
   *Biostatistics* 9:321-332.

.. [Yu2013] Yu D, Huber W, Vitek O (2013). Shrinkage estimation of dispersion
   in Negative Binomial models for RNA-seq experiments with small sample size.
   *Bioinformatics* 29:1275-1282.
