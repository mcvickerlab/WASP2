10X scRNA-seq Tutorial
======================

End-to-end allele-specific expression (ASE) workflow for 10X Chromium
scRNA-seq. Assumes a Cell Ranger ``possorted_genome_bam.bam`` with cell
barcodes in the ``CB`` tag, a phased VCF for the donor, and a cell-type
annotation from Seurat or Scanpy.

Inputs
------

- Cell Ranger BAM + index
- Phased VCF/BCF/PGEN for the sample
- A barcode-to-group TSV (see :doc:`/user_guide/single_cell` for Seurat /
  Scanpy export code and exact format)
- GTF gene annotation

Step 1 — Count alleles per cell per gene
-----------------------------------------

.. code-block:: bash

   wasp2-count count-variants-sc \
     cellranger_output/outs/possorted_genome_bam.bam \
     phased_variants.vcf.gz \
     barcodes_celltype.tsv \
     --feature genes.gtf \
     --samples SAMPLE_ID \
     --out_file allele_counts.h5ad

Output: an AnnData ``.h5ad`` with ``ref`` / ``alt`` / ``other`` layers,
genotype columns in ``.obs``, and cell-type assignments in ``.var``. See
:doc:`/user_guide/single_cell` for the full schema.

Step 2 — Per-cell-type imbalance
--------------------------------

.. code-block:: bash

   wasp2-analyze find-imbalance-sc \
     allele_counts.h5ad \
     barcodes_celltype.tsv \
     --sample SAMPLE_ID \
     --out_file imbalance_by_celltype.tsv

Output columns: ``region``, ``cell_type``, aggregated ``ref_count`` /
``alt_count``, ``pval``, ``fdr_pval``, ``effect_size`` (log₂ ref/alt).

Step 3 — Compare cell types (optional)
--------------------------------------

.. code-block:: bash

   wasp2-analyze compare-imbalance \
     allele_counts.h5ad \
     barcodes_celltype.tsv \
     --groups "CD4_T_cell,CD8_T_cell" \
     --out_file differential_imbalance.tsv

For more on comparative analysis (multiple groups, all-pairs, volcano
plots), see :doc:`comparative_imbalance`.

Interpreting results
--------------------

.. code-block:: python

   import pandas as pd

   results = pd.read_csv('imbalance_by_celltype.tsv', sep='\t')
   sig = results[results['fdr_pval'] < 0.05]

   top = (sig.groupby('cell_type')
             .apply(lambda x: x.nsmallest(10, 'fdr_pval'))
             .reset_index(drop=True))

   print(top[['region', 'cell_type', 'effect_size', 'fdr_pval']])

Troubleshooting
---------------

**Zero barcodes matched.** Confirm barcode format in the BAM vs. the TSV —
the ``CB:Z:...`` tag often has a ``-1`` suffix that your export must match:

.. code-block:: bash

   samtools view your.bam | head -10000 | grep -o 'CB:Z:[^[:space:]]*' \
     | cut -d: -f3 | sort -u > bam_bc.txt
   cut -f1 barcodes.tsv | sort -u > file_bc.txt
   comm -12 bam_bc.txt file_bc.txt | wc -l   # should be > 0

Fix a missing suffix with:

.. code-block:: bash

   awk -F'\t' '{print $1"-1\t"$2}' barcodes_no_suffix.tsv > barcodes.tsv

**Sparse counts.** Single-cell data is sparse. Consider pseudobulk
aggregation by cell type, lower ``--min`` / ``--min_count``, or focus on
highly expressed genes.

**Memory.** For large cohorts, split the region file by chromosome and
process chunks, then concatenate results.

See Also
--------

- :doc:`/user_guide/single_cell` — barcode file format, Seurat/Scanpy export
- :doc:`/user_guide/analysis` — analysis CLI reference
- :doc:`/methods/statistical_models` — beta-binomial LRT
- :doc:`scatac_workflow` — sibling tutorial for scATAC-seq
- :doc:`comparative_imbalance` — comparing groups
