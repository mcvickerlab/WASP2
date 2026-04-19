Beta-Binomial Model for Allelic Imbalance
==========================================

WASP2 detects allelic imbalance via a beta-binomial likelihood-ratio test.
The beta-binomial accommodates overdispersion from biological variation,
PCR amplification, and aggregation of counts across SNPs within a region —
sources of variance that a simple binomial model cannot absorb.

Model
-----

With :math:`X` the reference count out of :math:`N` total at a site:

.. math::

   p &\sim \text{Beta}(\alpha, \beta) \\
   X \mid p &\sim \text{Binomial}(N, p)

WASP2 uses the mean-dispersion parameterization, :math:`\mu = \alpha/(\alpha+\beta)`
and :math:`\rho = 1/(\alpha+\beta+1)`:

.. math::

   \alpha = \mu \cdot \frac{1 - \rho}{\rho}, \qquad
   \beta  = (1 - \mu) \cdot \frac{1 - \rho}{\rho}

The variance,

.. math::

   \text{Var}(X) = N\mu(1-\mu)\left[1 + (N-1)\rho\right],

recovers the binomial as :math:`\rho \to 0`. Estimation of :math:`\rho` is
described in :doc:`dispersion_estimation`.

Hypothesis Test
---------------

For each region, WASP2 tests

.. math::

   H_0 : \mu = 0.5 \qquad \text{vs.} \qquad H_1 : \mu \neq 0.5

with the likelihood ratio statistic

.. math::

   \Lambda = -2\left[\log \mathcal{L}_0 - \log \mathcal{L}_1\right], \qquad
   \Lambda \sim \chi^2_1 \text{ under } H_0,

where :math:`\mathcal{L}_1` is maximized over :math:`\mu` with :math:`\rho`
held at its null-model MLE (profile likelihood in :math:`\mu`). The p-value
is :math:`P(\chi^2_1 > \Lambda)`.

MLE for :math:`\mu` under :math:`H_1`:

.. code-block:: python

   from scipy.optimize import minimize_scalar
   from scipy.stats import betabinom
   import numpy as np

   def neg_log_likelihood(mu, ref_counts, n_counts, rho):
       alpha = mu * (1 - rho) / rho
       beta = (1 - mu) * (1 - rho) / rho
       return -np.sum(betabinom.logpmf(ref_counts, n_counts, alpha, beta))

   mu_mle = minimize_scalar(
       neg_log_likelihood,
       args=(ref_counts, n_counts, rho),
       method='bounded',
       bounds=(0, 1),
   ).x

Dispersion models
-----------------

WASP2 provides a single pooled :math:`\rho` (default) and a coverage-dependent
``linear`` model (logit-linear in :math:`N`). See :doc:`dispersion_estimation`
for the choice criteria and code.

.. code-block:: python

   from analysis.as_analysis import single_model, linear_model

   results = single_model(df, region_col="region", phased=False)
   results = linear_model(df, region_col="region", phased=False)

Phased vs. unphased data
------------------------

**Phased** genotypes reduce to a direct beta-binomial per SNP, with
:math:`\mu` or :math:`1-\mu` depending on which haplotype carries the
reference allele:

.. math::

   P(X_i \mid \text{hap}_i) = \text{BetaBinom}(X_i; N_i, \mu_{\text{hap}_i}, \rho)

**Unphased** analysis marginalizes over phase configurations with a uniform
prior across the :math:`2^{n-1}` phase assignments for a region of :math:`n`
SNPs, using dynamic programming over SNPs (approach introduced in
[Mayba2014]_):

.. math::

   P(\mathbf{X}) = \sum_{\phi} P(\mathbf{X} \mid \phi) \, P(\phi),
   \qquad P(\phi) = 0.5^{n-1}.

Output
------

.. table:: Columns emitted by the analysis pipeline
   :widths: 20 15 65

   ============ ======== ======================================================
   Column       Type     Interpretation
   ============ ======== ======================================================
   null_ll      float    Log-likelihood under :math:`H_0` (:math:`\mu=0.5`)
   alt_ll       float    Log-likelihood under :math:`H_1` (:math:`\mu=\hat\mu_{\text{MLE}}`)
   mu           float    MLE of the reference-allele proportion
   lrt          float    Likelihood ratio statistic
   pval         float    :math:`\chi^2_1` p-value
   fdr_pval     float    BH-adjusted p-value (see :doc:`fdr_correction`)
   ============ ======== ======================================================

:math:`\mu > 0.5` indicates a reference-allele preference; :math:`\mu < 0.5`
indicates an alternate preference; :math:`|\mu - 0.5|` is the effect size.

Practical notes
---------------

**Pseudocounts.** WASP2 adds :math:`c = 1` to both allele counts by default,
shrinking estimates slightly toward 0.5. This gives conservative inference and
avoids log-zero issues at extreme counts.

**Minimum counts.** Regions with :math:`N < \text{min\_count}` (default 10) are
dropped — low coverage yields poor power and unstable MLEs.

**Region aggregation.** Analyzing at the region level (peaks, genes) rather
than individual SNPs increases power, reduces the multiple-testing burden, and
better captures regulatory effects that act across a region.

**Power.** Power depends jointly on :math:`N`, :math:`\mu`, and :math:`\rho`.
Because WASP2 does not compute analytic power estimates at runtime, we do not
tabulate values here — power simulations for a given dataset should be run at
the dataset's own estimated :math:`\hat\rho` rather than relying on generic
tables.

See Also
--------

- :doc:`dispersion_estimation` — estimating :math:`\rho`
- :doc:`fdr_correction` — multiple-testing correction
- :doc:`counting_algorithm` — how counts are produced
- :doc:`mapping_filter` — WASP remap-and-filter producing the BAM

References
----------

.. [Mayba2014] Mayba O, Gilbert HN, Liu J, et al. (2014). MBASED:
   allele-specific expression detection in cancer tissues and cell lines.
   *Genome Biology* 15:405. *Introduces the phase-marginalized mixture model.*

.. [Kumasaka2016] Kumasaka N, Knights AJ, Gaffney DJ (2016). Fine-mapping
   cellular QTLs with RASQUAL and ATAC-seq. *Nature Genetics* 48:206-213.
   *Beta-binomial LRT framework in which WASP2's analysis sits.*

.. [vandeGeijn2015] van de Geijn B, McVicker G, Gilad Y, Pritchard JK (2015).
   WASP: allele-specific software for robust molecular quantitative trait
   locus discovery. *Nature Methods* 12:1061-1063. *The mapping-bias
   correction upstream of the count step.*
