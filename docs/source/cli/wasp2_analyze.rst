wasp2-analyze
=============

Command-line interface for the WASP2 analysis module.

.. contents:: Commands
   :local:
   :depth: 2

Overview
--------

The ``wasp2-analyze`` command performs statistical testing for allelic
imbalance. It provides three subcommands:

* ``find-imbalance`` - Test for AI in bulk data
* ``find-imbalance-sc`` - Test for AI in single-cell data (per cell type)
* ``compare-imbalance`` - Compare AI between cell types

Global Options
--------------

.. option:: --help

   Show help message and exit

.. option:: --version

   Show version number and exit

find-imbalance
--------------

Test for allelic imbalance in bulk sequencing data.

Synopsis
~~~~~~~~

.. code-block:: bash

   wasp2-analyze find-imbalance [OPTIONS] COUNTS

Arguments
~~~~~~~~~

.. option:: COUNTS

   Allele count file from ``wasp2-count`` (TSV format).

Options
~~~~~~~

.. option:: -o <PATH>, --out_file <PATH>

   Output file path. Default: ``ai_results.tsv``

.. option:: --min <N>

   Minimum total allele count required. Default: 10

.. option:: -p <N>, --pseudocount <N>

   Pseudocount added to avoid division by zero. Default: 1

.. option:: --phased

   Use phased haplotype model (requires phased genotypes).

.. option:: --region_col <NAME>

   Column name for region identifier. Auto-detected if not specified.

.. option:: --groupby <NAME>

   Aggregate counts by this column before testing.
   Use for gene-level analysis: ``--groupby gene_id``

Examples
~~~~~~~~

**Basic analysis**:

.. code-block:: bash

   wasp2-analyze find-imbalance counts.tsv \
       --out_file results.tsv

**With minimum count filter**:

.. code-block:: bash

   wasp2-analyze find-imbalance counts.tsv \
       --min 20 \
       --out_file results.tsv

**Gene-level aggregation (RNA-seq)**:

.. code-block:: bash

   wasp2-analyze find-imbalance gene_counts.tsv \
       --min 20 \
       --groupby gene_id \
       --out_file gene_results.tsv

Output Format
~~~~~~~~~~~~~

Tab-separated file with columns:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Column
     - Description
   * - ``region``
     - Region/gene/peak identifier
   * - ``ref_count``
     - Total reference allele reads
   * - ``alt_count``
     - Total alternate allele reads
   * - ``total``
     - Total reads (ref + alt)
   * - ``ratio``
     - Reference allele fraction (ref/total)
   * - ``pvalue``
     - Beta-binomial test p-value
   * - ``fdr``
     - Benjamini-Hochberg adjusted p-value

Example output:

.. code-block:: text

   region           ref_count  alt_count  total  ratio   pvalue      fdr
   ENSG00000123456  45         38         83     0.542   0.4521      0.8234
   ENSG00000789012  120        60         180    0.667   0.00012     0.0034

find-imbalance-sc
-----------------

Test for allelic imbalance in single-cell data, per cell type.

Synopsis
~~~~~~~~

.. code-block:: bash

   wasp2-analyze find-imbalance-sc [OPTIONS] COUNTS BARCODE_MAP

Arguments
~~~~~~~~~

.. option:: COUNTS

   Single-cell allele counts (h5ad format from ``count-variants-sc``).

.. option:: BARCODE_MAP

   Two-column TSV mapping barcodes to cell types.

   Format:

   .. code-block:: text

      AAACCTGAGAAACCAT-1	CD4_T
      AAACCTGAGAAACCGC-1	CD4_T
      AAACCTGAGAAACCTA-1	CD8_T

Options
~~~~~~~

.. option:: -o <PATH>, --out_file <PATH>

   Output file prefix. Creates one file per cell type.
   Default: ``ai_results``

.. option:: --min <N>

   Minimum total count per region (aggregated across cells). Default: 10

.. option:: -p <N>, --pseudocount <N>

   Pseudocount. Default: 1

.. option:: -s <SAMPLE>, --sample <SAMPLE>

   Sample ID for heterozygous SNP filtering.
   Required if multiple samples in data.

.. option:: --phased

   Use phased haplotype model.

.. option:: --unphased

   Explicitly use unphased model (default).

.. option:: -z <N>, --z_cutoff <N>

   Remove SNPs with counts exceeding this Z-score (QC filter).

Example
~~~~~~~

.. code-block:: bash

   wasp2-analyze find-imbalance-sc sc_counts.h5ad barcode_celltype.tsv \
       --min 20 \
       --sample donor1 \
       --out_file ai_results

Output: ``ai_results_CD4_T.tsv``, ``ai_results_CD8_T.tsv``, etc.

compare-imbalance
-----------------

Compare allelic imbalance between cell types/groups.

Synopsis
~~~~~~~~

.. code-block:: bash

   wasp2-analyze compare-imbalance [OPTIONS] COUNTS BARCODE_MAP

Arguments
~~~~~~~~~

.. option:: COUNTS

   Single-cell allele counts (h5ad format).

.. option:: BARCODE_MAP

   Two-column TSV mapping barcodes to cell types.

Options
~~~~~~~

.. option:: -o <PATH>, --out_file <PATH>

   Output file prefix. Creates one file per comparison.
   Default: ``ai_results``

.. option:: --groups <G1,G2,...>, --celltypes <G1,G2,...>

   Specific groups to compare (comma-separated).
   If not specified, compares all pairs.

.. option:: --min <N>

   Minimum total count per region in each group. Default: 10

.. option:: -p <N>, --pseudocount <N>

   Pseudocount. Default: 1

.. option:: -s <SAMPLE>, --sample <SAMPLE>

   Sample ID for heterozygous SNP filtering.

.. option:: --phased

   Use phased haplotype model.

.. option:: -z <N>, --z_cutoff <N>

   Z-score cutoff for QC filtering.

Examples
~~~~~~~~

**Compare two specific cell types**:

.. code-block:: bash

   wasp2-analyze compare-imbalance sc_counts.h5ad barcode_celltype.tsv \
       --groups CD4_T,CD8_T \
       --min 10 \
       --out_file diff_ai

Output: ``diff_ai_CD4_T_CD8_T.tsv``

**Compare all cell type pairs**:

.. code-block:: bash

   wasp2-analyze compare-imbalance sc_counts.h5ad barcode_celltype.tsv \
       --min 10 \
       --out_file all_comparisons

Output: ``all_comparisons_CD4_T_CD8_T.tsv``, ``all_comparisons_CD4_T_B_cell.tsv``, etc.

Output Format
~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Column
     - Description
   * - ``region``
     - Region identifier
   * - ``ratio_group1``
     - Allelic ratio in first group
   * - ``ratio_group2``
     - Allelic ratio in second group
   * - ``ratio_diff``
     - Difference (group1 - group2)
   * - ``pvalue``
     - Statistical test p-value
   * - ``fdr``
     - Adjusted p-value

Statistical Methods
-------------------

Beta-Binomial Model
~~~~~~~~~~~~~~~~~~~

WASP2 uses the beta-binomial distribution which accounts for overdispersion
in sequencing data. This provides more accurate p-values than simple
binomial tests.

The model assumes:

.. math::

   p \sim \text{Beta}(\alpha, \beta)

   k \mid p \sim \text{Binomial}(n, p)

where:
- :math:`k` = reference allele count
- :math:`n` = total count
- :math:`p` = underlying allelic ratio

Multiple Testing Correction
~~~~~~~~~~~~~~~~~~~~~~~~~~~

FDR (False Discovery Rate) values are computed using the
Benjamini-Hochberg procedure. Use FDR < 0.05 as the significance threshold.

Interpretation Guide
--------------------

Significant Results
~~~~~~~~~~~~~~~~~~~

- **FDR < 0.05**: Statistically significant allelic imbalance
- **ratio < 0.35 or > 0.65**: Strong effect size
- **ratio = 0.5**: Balanced (no imbalance)

Filtering Recommendations
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Significant with strong effect
   awk 'NR==1 || ($7 < 0.05 && ($5 < 0.35 || $5 > 0.65))' results.tsv

   # Top 20 most imbalanced
   sort -t$'\t' -k5,5n results.tsv | head -21

See Also
--------

* :doc:`/tutorials/concepts` - Statistical background
* :doc:`/tutorials/rnaseq_ase` - RNA-seq analysis tutorial
* :doc:`/tutorials/single_cell` - Single-cell analysis tutorial
