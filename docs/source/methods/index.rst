Statistical Methods
===================

This section provides detailed documentation of the statistical methods,
algorithms, and biological rationale underlying WASP2's allele-specific
analysis pipeline.

.. toctree::
   :maxdepth: 2
   :caption: Contents

   counting_algorithm
   mapping_filter
   statistical_models
   dispersion_estimation
   fdr_correction

Overview
--------

WASP2 implements a complete pipeline for allele-specific analysis:

1. **Allele Counting**: Reads are assigned to reference or alternate alleles
   at heterozygous variant sites using base-level alignment information.

2. **Mapping Bias Correction**: The WASP algorithm removes reads that exhibit
   mapping bias by testing whether allele-swapped reads map to the same location.

3. **Statistical Testing**: Beta-binomial models account for overdispersion
   in allele count data when testing for allelic imbalance.

4. **Multiple Testing Correction**: False discovery rate control ensures
   reliable detection of true imbalanced regions.

References
----------

.. [vandeGeijn2015] van de Geijn B, McVicker G, Gilad Y, Pritchard JK (2015).
   WASP: allele-specific software for robust molecular quantitative trait
   locus discovery. *Nature Methods* 12:1061-1063.

.. [Castel2015] Castel SE, Levy-Moonshine A, Mohammadi P, Banks E, Lappalainen T (2015).
   Tools and best practices for data processing in allelic expression analysis.
   *Genome Biology* 16:195.

.. [Skelly2011] Skelly DA, Johansson M, Madeoy J, Wakefield J, Akey JM (2011).
   A powerful and flexible statistical framework for testing hypotheses of
   allele-specific gene expression from RNA-seq data. *Genome Research* 21:1728-1737.
