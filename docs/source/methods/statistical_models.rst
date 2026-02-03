Beta-Binomial Model for Allelic Imbalance
==========================================

This document describes the statistical framework WASP2 uses to detect
allelic imbalance from allele count data.

.. contents:: Contents
   :local:
   :depth: 2

Motivation
----------

Why Not Use the Binomial?
^^^^^^^^^^^^^^^^^^^^^^^^^

The simplest model for allele counts is the binomial distribution. If reads
are sampled independently from two alleles with equal probability:

.. math::

   X \sim \text{Binomial}(N, 0.5)

where :math:`X` is the reference count and :math:`N` is the total count.

However, real allele count data exhibits **overdispersion**: the variance
exceeds the binomial expectation. Sources of overdispersion include:

- **Biological variation**: True allelic imbalance varies across cells
- **Technical noise**: PCR amplification introduces correlated errors
- **Aggregation effects**: Combining counts across SNPs within a region
- **Sampling from a population**: Different individuals have different AI

The Beta-Binomial Model
-----------------------

The beta-binomial distribution naturally accommodates overdispersion by
modeling the success probability as a random variable:

.. math::

   p &\sim \text{Beta}(\alpha, \beta) \\
   X | p &\sim \text{Binomial}(N, p)

Marginalizing over :math:`p` gives the beta-binomial:

.. math::

   P(X = k) = \binom{N}{k} \frac{B(k + \alpha, N - k + \beta)}{B(\alpha, \beta)}

where :math:`B(\cdot, \cdot)` is the beta function.

Parameterization
^^^^^^^^^^^^^^^^

WASP2 uses the **mean-dispersion parameterization**:

.. math::

   \mu &= \frac{\alpha}{\alpha + \beta} \\
   \rho &= \frac{1}{\alpha + \beta + 1}

The dispersion parameter :math:`\rho \in (0, 1)` controls overdispersion:

- :math:`\rho \to 0`: Approaches binomial (no overdispersion)
- :math:`\rho \to 1`: Maximum overdispersion (all probability at 0 or N)

The inverse transformation:

.. math::

   \alpha &= \mu \cdot \frac{1 - \rho}{\rho} \\
   \beta &= (1 - \mu) \cdot \frac{1 - \rho}{\rho}

Variance Structure
^^^^^^^^^^^^^^^^^^

The beta-binomial variance is:

.. math::

   \text{Var}(X) = N\mu(1-\mu) \left[ 1 + (N-1)\rho \right]

The factor :math:`[1 + (N-1)\rho]` is the **variance inflation factor**.
For typical values (:math:`\rho \approx 0.01`, :math:`N \approx 100`), this
gives roughly 2x the binomial variance.

Hypothesis Testing
------------------

WASP2 tests for allelic imbalance using a likelihood ratio test:

**Null Hypothesis** :math:`H_0`: No imbalance (:math:`\mu = 0.5`)

**Alternative Hypothesis** :math:`H_1`: Imbalance present (:math:`\mu \neq 0.5`)

Likelihood Functions
^^^^^^^^^^^^^^^^^^^^

Under the null hypothesis:

.. math::

   \mathcal{L}_0 = \prod_{i=1}^{n} P(X_i | N_i, \mu=0.5, \rho)

Under the alternative:

.. math::

   \mathcal{L}_1 = \prod_{i=1}^{n} P(X_i | N_i, \hat{\mu}_{\text{MLE}}, \rho)

where the product is over SNPs within a region (peak, gene).

Likelihood Ratio Test
^^^^^^^^^^^^^^^^^^^^^

The test statistic:

.. math::

   \Lambda = -2 \left[ \log \mathcal{L}_0 - \log \mathcal{L}_1 \right]

Under the null hypothesis, :math:`\Lambda` follows a chi-squared distribution
with 1 degree of freedom:

.. math::

   \Lambda \sim \chi^2_1

The p-value is:

.. math::

   p = P(\chi^2_1 > \Lambda)

MLE Estimation
^^^^^^^^^^^^^^

The MLE for :math:`\mu` under the alternative is found by numerical optimization:

.. code-block:: python

   from scipy.optimize import minimize_scalar
   from scipy.stats import betabinom

   def neg_log_likelihood(mu, ref_counts, n_counts, rho):
       alpha = mu * (1 - rho) / rho
       beta = (1 - mu) * (1 - rho) / rho
       return -np.sum(betabinom.logpmf(ref_counts, n_counts, alpha, beta))

   result = minimize_scalar(
       neg_log_likelihood,
       args=(ref_counts, n_counts, rho),
       method='bounded',
       bounds=(0, 1)
   )
   mu_mle = result.x

Implementation in WASP2
-----------------------

Single Dispersion Model
^^^^^^^^^^^^^^^^^^^^^^^

The default model estimates a single :math:`\rho` for all data:

.. code-block:: python

   from analysis.as_analysis import single_model

   results = single_model(df, region_col="region", phased=False)

This assumes homogeneous overdispersion across regions, which is often
reasonable for moderately-sized datasets.

Linear Dispersion Model
^^^^^^^^^^^^^^^^^^^^^^^

For large datasets, WASP2 offers a linear model where dispersion varies
with total count:

.. math::

   \text{logit}(\rho) = \beta_0 + \beta_1 \cdot N

This captures the observation that regions with higher coverage often show
different dispersion characteristics:

.. code-block:: python

   from analysis.as_analysis import linear_model

   results = linear_model(df, region_col="region", phased=False)

Phased vs Unphased Analysis
---------------------------

WASP2 supports both phased and unphased genotype data.

Unphased Model
^^^^^^^^^^^^^^

When genotype phase is unknown, WASP2 uses a mixture model that marginalizes
over possible phase configurations:

For a region with multiple SNPs, if we don't know which haplotype each
ref allele belongs to, we sum over phase assignments using dynamic programming:

.. math::

   P(\mathbf{X}) = \sum_{\phi \in \text{phases}} P(\mathbf{X} | \phi) \cdot P(\phi)

where :math:`\phi` indexes phase configurations with prior :math:`P(\phi) = 0.5^{n-1}`.

Phased Model
^^^^^^^^^^^^

With phased genotypes (e.g., from read-backed phasing), the model is simpler:

.. math::

   P(X_i | \text{hap}_i) = \text{BetaBinom}(X_i; N_i, \mu_{\text{hap}_i}, \rho)

where :math:`\mu_{\text{hap}_i}` is :math:`\mu` or :math:`1-\mu` depending on
which haplotype the reference allele belongs to.

Output Interpretation
---------------------

WASP2 returns the following statistics for each region:

.. table:: Beta-Binomial Analysis Output
   :widths: 20 15 65

   ============ ======== ======================================================
   Column       Type     Interpretation
   ============ ======== ======================================================
   null_ll      float    Log-likelihood under null (μ=0.5)
   alt_ll       float    Log-likelihood under alternative (μ=MLE)
   mu           float    MLE of imbalance proportion
   lrt          float    Likelihood ratio test statistic
   pval         float    p-value from χ² distribution
   fdr_pval     float    FDR-corrected p-value (BH method)
   ============ ======== ======================================================

**Interpreting μ (mu)**:

- :math:`\mu = 0.5`: No imbalance (equal allele expression)
- :math:`\mu > 0.5`: Reference allele preference
- :math:`\mu < 0.5`: Alternate allele preference
- :math:`|\mu - 0.5|`: Effect size (magnitude of imbalance)

Practical Considerations
------------------------

Pseudocounts
^^^^^^^^^^^^

WASP2 adds pseudocounts to avoid log(0) issues:

.. math::

   X' = X + c, \quad N' = N + 2c

Default :math:`c = 1` (Laplace smoothing). This slightly shrinks estimates
toward 0.5, providing conservative inference.

Minimum Count Threshold
^^^^^^^^^^^^^^^^^^^^^^^

Regions with low total counts have poor statistical power. WASP2 filters
regions with :math:`N < \text{min\_count}` (default: 10).

Power depends on coverage and effect size:

.. table:: Approximate Power (α=0.05)
   :widths: 20 20 20 20 20

   ========= ======== ======== ======== ========
   Total N   μ=0.55   μ=0.60   μ=0.65   μ=0.70
   ========= ======== ======== ======== ========
   20        5%       10%      20%      35%
   50        8%       25%      50%      75%
   100       15%      45%      80%      95%
   200       30%      75%      95%      99%
   ========= ======== ======== ======== ========

Region Aggregation
^^^^^^^^^^^^^^^^^^

Analyzing at the region level (genes, peaks) rather than individual SNPs:

- **Increases power**: More counts per test
- **Reduces multiple testing burden**: Fewer tests to correct
- **Captures regulatory effects**: ASE often affects entire genes

See Also
--------

- :doc:`dispersion_estimation` - Estimating the dispersion parameter
- :doc:`fdr_correction` - Multiple testing correction methods

References
----------

.. [Mayba2014] Mayba O, Gilbert HN, Liu J, et al. (2014). MBASED: allele-specific
   expression detection in cancer tissues and cell lines. *Genome Biology* 15:405.
