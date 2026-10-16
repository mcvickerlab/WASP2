Single-Cell ATAC-seq Workflow
=================================

Count alleles in a barcode-tagged scATAC-seq BAM and test imbalance after
aggregating cells into annotated groups. The counter does not perform
UMI-aware molecule counting for scRNA-seq.

Inputs
----------

* Coordinate-sorted, indexed BAM with cell barcodes in the ``CB`` tag.
* Donor VCF with heterozygous SNPs and matching reference coordinates.
* BED peak file.
* One-column barcode whitelist for counting, without a header.
* Two-column barcode-to-group TSV for analysis, without a header:

.. code-block:: text

   AAACCCAAGAAACACT-1	B_cell
   AAACCCAAGAAAGCGA-1	T_cell

Derive the counting whitelist from this mapping:

.. code-block:: bash

   cut -f1 barcode_groups.tsv > barcodes.txt

Use a BAM that has undergone appropriate mapping-bias correction; see
:doc:`/user_guide/mapping`. Fragment-coordinate files alone do not contain
base-level evidence for REF and ALT allele counts.

Step 1 — Count alleles per cell
-----------------------------------

.. code-block:: bash

   wasp2-count count-variants-sc \
     sample.bam variants.vcf.gz barcodes.txt \
     --feature peaks.bed \
     --samples SAMPLE_ID \
     --out_file allele_counts.h5ad

AnnData has cells in ``obs``, SNPs and genotypes in ``var``, and ``ref``,
``alt``, and ``other`` layers of shape ``cells × SNPs``. Peak membership is
stored in ``uns['feature']``.

Step 2 — Test imbalance within each group
---------------------------------------------

.. code-block:: bash

   wasp2-analyze find-imbalance-sc \
     allele_counts.h5ad barcode_groups.tsv \
     --sample SAMPLE_ID \
     --phased --min 10 \
     --out_file imbalance.tsv

This writes one file per group, such as ``imbalance_B_cell.tsv``. Columns
are ``region``, ``num_snps``, ``mu``, ``null_ll``, ``alt_ll``, ``pval``, and
``fdr_pval``. The filename identifies the group; there is no ``cell_type``
column or reported ``effect_size``.

Phasing is optional. Omit ``--phased`` to use automatic genotype checks,
or specify ``--unphased``. A phased request falls back to the unphased
model when the selected genotypes do not satisfy the phase checks.

Step 3 — Compare groups
---------------------------

.. code-block:: bash

   wasp2-analyze compare-imbalance \
     allele_counts.h5ad barcode_groups.tsv \
     --sample SAMPLE_ID \
     --groups B_cell --groups T_cell \
     --phased \
     --out_file comparison.tsv

Omit ``--groups`` to compare all available groups pairwise. One file is
written per pair, such as ``comparison_B_cell_T_cell.tsv``. Columns include
``region``, ``num_snps``, ``combined_mu``, ``mu1``, ``mu2``, ``null_ll``,
``alt_ll``, ``pval``, and ``fdr_pval``.

Interpretation
------------------

.. code-block:: python

   import pandas as pd

   results = pd.read_csv('imbalance_B_cell.tsv', sep='\t')
   significant = results.loc[results['fdr_pval'] < 0.05]
   print(significant[['region', 'mu', 'fdr_pval']])

``mu`` is a fitted allele/haplotype proportion in the model's regional
orientation. It is not a gene-wide REF expression fraction or a log-fold
change. Comparing direction across regions or datasets requires consistent
SNP and haplotype orientation.

Troubleshooting
-------------------

**No matching barcodes.** Match whitelist and annotation barcodes to BAM
``CB`` tags exactly, including suffixes such as ``-1``.

**Sparse counts.** Aggregate cells into biologically meaningful groups and
inspect coverage before choosing a minimum-count threshold. The CLI tests
group-level aggregate counts; it has no separate per-cell model.

**Coverage outliers.** The optional ``-z`` filter removes SNPs and associated
regions with extreme coverage. Apply a consistent filtering policy when
comparing groups.

See Also
------------

* :doc:`/user_guide/single_cell` — barcode files and AnnData schema
* :doc:`/user_guide/analysis` — models and command options
* :doc:`/methods/mapping_filter` — mapping-bias filter contract
