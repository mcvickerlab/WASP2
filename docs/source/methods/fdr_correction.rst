False Discovery Rate Correction
================================

This document describes the multiple testing correction methods used in WASP2
to control false positive rates when testing many genomic regions.

.. contents:: Contents
   :local:
   :depth: 2

The Multiple Testing Problem
----------------------------

When testing thousands of genomic regions for allelic imbalance, even a small
per-test false positive rate leads to many false discoveries:

.. math::

   E[\text{false positives}] = m \cdot \alpha

For :math:`m = 10{,}000` tests at :math:`\alpha = 0.05`, we expect 500 false
positives by chance alone.

**Example**: Testing 20,000 genes for ASE

- At α = 0.05: ~1,000 expected false positives
- At α = 0.01: ~200 expected false positives
- Even stringent thresholds yield many false discoveries

Error Rate Definitions
----------------------

Family-Wise Error Rate (FWER)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The probability of making **any** false discovery:

.. math::

   \text{FWER} = P(V \geq 1)

where :math:`V` is the number of false positives.

FWER control (e.g., Bonferroni) is very conservative for genomic studies.

False Discovery Rate (FDR)
^^^^^^^^^^^^^^^^^^^^^^^^^^

The expected **proportion** of false discoveries among rejections:

.. math::

   \text{FDR} = E\left[\frac{V}{R}\right]

where :math:`V` is false positives and :math:`R` is total rejections.

FDR is more appropriate for discovery-oriented genomic studies where:

- Some false positives are acceptable
- The goal is to prioritize candidates for follow-up
- The number of tests is very large

Benjamini-Hochberg Procedure
----------------------------

WASP2 uses the Benjamini-Hochberg (BH) procedure [Benjamini1995]_ for FDR control.

Algorithm
^^^^^^^^^

Given :math:`m` p-values :math:`p_1, p_2, \ldots, p_m`:

1. Sort p-values: :math:`p_{(1)} \leq p_{(2)} \leq \cdots \leq p_{(m)}`
2. Find the largest :math:`k` such that :math:`p_{(k)} \leq \frac{k}{m} \cdot q`
3. Reject all hypotheses with :math:`p_{(i)} \leq p_{(k)}`

where :math:`q` is the target FDR level (typically 0.05 or 0.1).

**Adjusted P-Values (q-values)**:

WASP2 reports BH-adjusted p-values:

.. math::

   q_i = \min_{j \geq i} \left\{ \frac{m \cdot p_{(j)}}{j} \right\}

These can be interpreted as the minimum FDR at which the hypothesis would be
rejected.

Implementation
^^^^^^^^^^^^^^

.. code-block:: python

   from scipy.stats import false_discovery_control

   # WASP2 uses the BH method
   fdr_pvals = false_discovery_control(pvals, method='bh')

   # Equivalent manual implementation
   def benjamini_hochberg(pvals):
       n = len(pvals)
       ranked = np.argsort(pvals)
       adjusted = np.empty(n)
       cummin = 1.0
       for i in range(n - 1, -1, -1):
           idx = ranked[i]
           adjusted[idx] = min(cummin, pvals[idx] * n / (i + 1))
           cummin = adjusted[idx]
       return adjusted

Properties
^^^^^^^^^^

**Advantages**:

- Controls FDR at level :math:`q` under independence
- More powerful than FWER methods (fewer false negatives)
- Simple to compute and interpret
- Works well with continuous p-values

**Assumptions**:

- P-values under the null are uniformly distributed
- Independence or positive regression dependence (PRDS)

The chi-squared p-values from WASP2's likelihood ratio test satisfy these
assumptions when regions are independent.

Alternative Methods
-------------------

While WASP2 uses BH by default, researchers may consider alternatives
for specific scenarios.

Benjamini-Yekutieli (BY)
^^^^^^^^^^^^^^^^^^^^^^^^

For arbitrary dependence between tests:

.. math::

   p_{(k)} \leq \frac{k}{m \cdot c(m)} \cdot q

where :math:`c(m) = \sum_{i=1}^{m} 1/i \approx \ln(m) + 0.577`.

BY is more conservative but valid under any dependence structure.

.. code-block:: python

   # Available in scipy
   fdr_pvals = false_discovery_control(pvals, method='by')

Storey's q-value
^^^^^^^^^^^^^^^^

Estimates the proportion of true nulls (:math:`\pi_0`) for improved power:

.. math::

   \hat{\pi}_0 = \frac{\#\{p_i > \lambda\}}{m(1 - \lambda)}

The q-value procedure is more powerful when many tests are true discoveries.

.. code-block:: python

   # Requires qvalue package
   # pip install qvalue
   from qvalue import qvalue
   q = qvalue(pvals)

Alternative Correction Methods
------------------------------

**Discrete FDR and Mid-P Adjustments**

For exact tests with discrete p-values (Fisher's exact, exact binomial),
specialized methods like Gilbert's procedure [Gilbert2005]_ or mid-p
adjustments can reduce conservativeness.

However, WASP2's likelihood ratio test produces continuous p-values from
the chi-squared distribution, so **standard BH is appropriate** and these
discrete methods are not needed.

Practical Guidelines
--------------------

Choosing an FDR Threshold
^^^^^^^^^^^^^^^^^^^^^^^^^

.. table:: FDR Threshold Guidelines
   :widths: 20 40 40

   ======== ================================= ================================
   FDR      Use Case                          Interpretation
   ======== ================================= ================================
   0.01     High-confidence discoveries       ~1% false among significant
   0.05     Standard exploratory analysis     ~5% false among significant
   0.10     Liberal discovery                 ~10% false, maximize sensitivity
   0.20     Hypothesis generation             For follow-up validation
   ======== ================================= ================================

When to Use Stricter Control
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Consider stricter FDR or FWER control when:

- Results will directly inform clinical decisions
- Follow-up validation is expensive
- The number of true positives is expected to be small
- Independence assumptions may be violated

Reporting Results
^^^^^^^^^^^^^^^^^

When reporting FDR-corrected results:

1. **Report the method**: "FDR correction was performed using the
   Benjamini-Hochberg procedure"
2. **Report the threshold**: "Significance was declared at FDR < 0.05"
3. **Report both p-values**: Include raw and adjusted p-values in supplements
4. **Report the number of tests**: "Among 15,234 tested regions..."

Output Format
-------------

WASP2 reports both raw and adjusted p-values:

.. table:: P-value Columns in Output
   :widths: 20 20 60

   ========= ======== ====================================================
   Column    Type     Description
   ========= ======== ====================================================
   pval      float    Raw p-value from likelihood ratio test
   fdr_pval  float    BH-adjusted p-value (q-value)
   ========= ======== ====================================================

**Interpretation**:

- ``pval < 0.05``: Nominally significant (not corrected)
- ``fdr_pval < 0.05``: Significant after multiple testing correction

References
----------

.. [Benjamini1995] Benjamini Y, Hochberg Y (1995). Controlling the false
   discovery rate: A practical and powerful approach to multiple testing.
   *Journal of the Royal Statistical Society B* 57:289-300.

.. [Storey2003] Storey JD, Tibshirani R (2003). Statistical significance for
   genomewide studies. *Proceedings of the National Academy of Sciences*
   100:9440-9445.

.. [Gilbert2005] Gilbert PB (2005). A modified false discovery rate
   multiple-comparisons procedure for discrete data, applied to human
   immunodeficiency virus genetics. *Journal of the Royal Statistical
   Society C* 54:143-158.
