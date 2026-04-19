Single-Cell Workflow (scRNA-seq / scATAC-seq)
==============================================

End-to-end allele-specific workflow for single-cell data — 10X Chromium
scRNA-seq and 10X scATAC-seq. Pipeline is the same in both cases; the
data-type difference shows up as GTF (for scRNA-seq genes) vs. BED (for
scATAC-seq peaks) in the ``--feature`` argument.

Inputs
------

- Cell Ranger BAM with cell barcodes in the ``CB:Z:...`` tag + index
- Phased VCF/BCF/PGEN for the donor
- Barcode-to-group TSV (cell type or other assignment — see
  :doc:`/user_guide/single_cell` for Seurat/Scanpy export code and format)
- **scRNA-seq**: GTF gene annotation
- **scATAC-seq**: BED peak file (usually from Cell Ranger
  ``filtered_peak_bc_matrix`` or a consensus peak set)

Step 1 — Count alleles per cell
--------------------------------

**scRNA-seq (genes):**

.. code-block:: bash

   wasp2-count count-variants-sc \
     cellranger_output/outs/possorted_genome_bam.bam \
     phased_variants.vcf.gz \
     barcodes_celltype.tsv \
     --feature genes.gtf \
     --samples SAMPLE_ID \
     --out_file allele_counts.h5ad

**scATAC-seq (peaks):**

.. code-block:: bash

   wasp2-count count-variants-sc \
     cellranger_output/outs/possorted_bam.bam \
     phased_variants.vcf.gz \
     barcodes_celltype.tsv \
     --feature peaks.bed \
     --samples SAMPLE_ID \
     --out_file allele_counts.h5ad

Output: an AnnData ``.h5ad`` with ``ref`` / ``alt`` / ``other`` layers,
genotype columns in ``.obs``, and cell-type assignments in ``.var``. See
:doc:`/user_guide/single_cell` for the full schema.

Step 2 — Per-group imbalance
----------------------------

.. code-block:: bash

   wasp2-analyze find-imbalance-sc \
     allele_counts.h5ad \
     barcodes_celltype.tsv \
     --sample SAMPLE_ID \
     --phased --min 10 -z 3 \
     --out_file imbalance_by_celltype.tsv

Output columns: ``region``, ``cell_type``, aggregated ``ref_count`` /
``alt_count``, ``pval``, ``fdr_pval``, ``effect_size`` (log₂ ref/alt).

Step 3 — Compare groups (optional)
-----------------------------------

.. code-block:: bash

   wasp2-analyze compare-imbalance \
     allele_counts.h5ad \
     barcodes_celltype.tsv \
     --groups "CD4_T_cell,CD8_T_cell" \
     --phased \
     --out_file differential_imbalance.tsv

Omit ``--groups`` to compare all available groups pairwise. See
:doc:`/user_guide/analysis` for the full CLI reference and output columns.

Per-cell vs. pseudo-bulk
------------------------

Single-cell ATAC data is especially sparse — most cells contribute zero
reads to most peaks. Two analysis modes are common:

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - Aspect
     - Per-cell
     - Pseudo-bulk (per-cell-type)
   * - Resolution
     - Single cell
     - Cell population
   * - Power
     - Low (sparse)
     - High (aggregated)
   * - Use case
     - Outlier cells
     - Population-level imbalance

Pseudo-bulk (the default, via the barcode-to-group TSV) is the right
starting point for most scATAC experiments. Per-cell analysis is useful
when investigating rare subpopulations or outlier effects.

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

Fix a missing suffix:

.. code-block:: bash

   awk -F'\t' '{print $1"-1\t"$2}' barcodes_no_suffix.tsv > barcodes.tsv

**Sparse counts / low power.** Aggregate to pseudo-bulk by cell type,
lower ``--min`` / ``--min_count``, or focus on highly expressed genes
(scRNA-seq) / high-coverage peaks (scATAC-seq).

**Memory.** For large cohorts, split the feature file by chromosome and
process chunks:

.. code-block:: bash

   for chr in chr{1..22}; do
     grep "^${chr}\s" peaks.bed > peaks_${chr}.bed
     wasp2-count count-variants-sc sample.bam variants.vcf.gz barcodes.tsv \
       --feature peaks_${chr}.bed --out_file counts_${chr}.h5ad
   done

See Also
--------

- :doc:`/user_guide/single_cell` — barcode format, Seurat/Scanpy export
- :doc:`/user_guide/analysis` — analysis CLI reference
- :doc:`/methods/statistical_models` — beta-binomial LRT
- :doc:`bulk_workflow` — sibling tutorial for bulk RNA-seq / ATAC-seq
