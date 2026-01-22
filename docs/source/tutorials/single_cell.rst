Single-Cell Allele-Specific Analysis
=====================================

**Time**: 60 minutes
**Level**: Advanced
**Prerequisites**: Completed :doc:`basic_workflow`, single-cell data

This tutorial covers allele-specific analysis for single-cell RNA-seq and
ATAC-seq data using WASP2's specialized single-cell tools.

.. contents:: Topics
   :local:
   :depth: 2

Overview
--------

Single-cell allele-specific analysis enables:

- **Cell-type-specific ASE**: Identify genes with different allelic ratios in different cell types
- **Differential AI**: Compare allelic imbalance between conditions or cell types
- **Single-cell resolution**: Track allelic patterns at individual cell level

Challenges in Single-Cell Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Single-cell data presents unique challenges:

1. **Sparse coverage**: Most genes have few reads per cell
2. **Technical noise**: Dropout, amplification bias
3. **Cell-type heterogeneity**: Different cell types may show different patterns
4. **Large scale**: Thousands of cells, millions of data points

WASP2's single-cell tools address these through:

- Per-cell barcode tracking
- Cell-type aggregation
- Specialized statistical models
- Efficient data formats (AnnData/h5ad)

Single-Cell Workflow Overview
-----------------------------

.. code-block:: text

   10x Genomics BAM
         |
         v
   +---------------------------+
   | count-variants-sc         |  <- Per-cell allele counting
   +---------------------------+
         |
         v
   AnnData (h5ad) - Cell x SNP matrix
         |
         v
   +---------------------------+
   | find-imbalance-sc         |  <- Cell-type-specific AI
   +---------------------------+
         |
         v
   Cell-type AI results
         |
         v
   +---------------------------+
   | compare-imbalance         |  <- Differential AI (optional)
   +---------------------------+
         |
         v
   Differential AI results

Step 1: Prepare Input Files
---------------------------

Cell barcodes
~~~~~~~~~~~~~

Extract valid cell barcodes from your filtered data:

.. code-block:: bash

   # 10x Genomics format
   zcat filtered_feature_bc_matrix/barcodes.tsv.gz > cell_barcodes.txt

   # Check barcode count
   wc -l cell_barcodes.txt

**Expected**: 1,000-20,000 cells depending on experiment.

Cell-type annotations
~~~~~~~~~~~~~~~~~~~~~

Create a two-column TSV mapping barcodes to cell types:

.. code-block:: text

   # barcode_celltype_map.tsv
   AAACCTGAGAAACCAT-1	CD4_T
   AAACCTGAGAAACCGC-1	CD4_T
   AAACCTGAGAAACCTA-1	CD8_T
   AAACCTGAGAAACGTC-1	B_cell
   AAACCTGAGAAAGTGG-1	Monocyte
   ...

This can come from clustering tools like Seurat, Scanpy, or Signac.

.. code-block:: bash

   # Verify format
   head barcode_celltype_map.tsv

   # Count cells per type
   cut -f2 barcode_celltype_map.tsv | sort | uniq -c

Region annotations (ATAC-seq)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For scATAC-seq, prepare peak annotations:

.. code-block:: bash

   # Use peaks from clustering (e.g., ArchR, Signac)
   ls -l peaks.bed

Step 2: Single-Cell Allele Counting
-----------------------------------

Count alleles per cell per SNP:

.. code-block:: bash

   wasp2-count count-variants-sc \
       possorted_genome_bam.bam \
       variants.vcf.gz \
       cell_barcodes.txt \
       --samples donor1 \
       --feature peaks.bed \
       --out_file sc_allele_counts.h5ad

**Parameters**:

- ``cell_barcodes.txt``: Valid cell barcodes
- ``--samples donor1``: Sample ID for het SNP filtering
- ``--feature peaks.bed``: Region annotations (optional for RNA-seq)
- ``--out_file``: Output in AnnData format

**Output**: ``sc_allele_counts.h5ad`` containing:

- ``.X``: Cell x SNP count matrix
- ``.var``: SNP annotations (chr, pos, ref, alt)
- ``.obs``: Cell annotations (barcode)

Inspect the output
~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import anndata as ad

   adata = ad.read_h5ad('sc_allele_counts.h5ad')
   print(f"Cells: {adata.n_obs}")
   print(f"SNPs: {adata.n_vars}")
   print(f"Total counts: {adata.X.sum()}")
   print(f"Mean counts per cell: {adata.X.sum(axis=1).mean():.1f}")

Step 3: Cell-Type-Specific Analysis
-----------------------------------

Analyze allelic imbalance within each cell type:

.. code-block:: bash

   wasp2-analyze find-imbalance-sc \
       sc_allele_counts.h5ad \
       barcode_celltype_map.tsv \
       --min 20 \
       --out_file ai_results

**Parameters**:

- ``--min 20``: Minimum reads per region (aggregated across cells)
- Output: ``ai_results_CD4_T.tsv``, ``ai_results_CD8_T.tsv``, etc.

Multi-sample handling
~~~~~~~~~~~~~~~~~~~~~

If your data contains multiple donors:

.. code-block:: bash

   # Specify which sample's genotypes to use
   wasp2-analyze find-imbalance-sc \
       sc_allele_counts.h5ad \
       barcode_celltype_map.tsv \
       --sample donor1 \
       --min 20 \
       --out_file donor1_ai_results

Step 4: Compare Allelic Imbalance Between Cell Types
----------------------------------------------------

Identify regions with differential AI between cell types:

.. code-block:: bash

   wasp2-analyze compare-imbalance \
       sc_allele_counts.h5ad \
       barcode_celltype_map.tsv \
       --groups CD4_T,CD8_T \
       --min 10 \
       --out_file CD4_vs_CD8_diff_ai.tsv

**Output columns**:

- ``region``: Peak/gene ID
- ``ratio_group1``, ``ratio_group2``: Allelic ratios in each group
- ``ratio_diff``: Difference in ratios
- ``pvalue``, ``fdr``: Statistical significance

Find all pairwise comparisons
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Compare all cell type pairs
   wasp2-analyze compare-imbalance \
       sc_allele_counts.h5ad \
       barcode_celltype_map.tsv \
       --min 10 \
       --out_file all_pairwise_diff_ai

This generates files for each pair: ``all_pairwise_diff_ai_CD4_T_CD8_T.tsv``, etc.

Quality Control
---------------

Z-score filtering
~~~~~~~~~~~~~~~~~

Remove outlier SNPs with unusually high counts (potential artifacts):

.. code-block:: bash

   wasp2-analyze find-imbalance-sc \
       sc_allele_counts.h5ad \
       barcode_celltype_map.tsv \
       --min 20 \
       --z_cutoff 3.0 \
       --out_file ai_results_filtered

Coverage distribution
~~~~~~~~~~~~~~~~~~~~~

Check coverage across cell types:

.. code-block:: python

   import anndata as ad
   import pandas as pd

   # Load data
   adata = ad.read_h5ad('sc_allele_counts.h5ad')
   celltype_map = pd.read_csv('barcode_celltype_map.tsv', sep='\t',
                               header=None, names=['barcode', 'celltype'])

   # Add cell type annotations
   adata.obs = adata.obs.join(celltype_map.set_index('barcode'))

   # Coverage per cell type
   for ct in adata.obs['celltype'].unique():
       cells = adata.obs['celltype'] == ct
       total = adata[cells].X.sum()
       n_cells = cells.sum()
       print(f"{ct}: {n_cells} cells, {total:.0f} total reads")

Visualization
-------------

Python visualization with Scanpy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import anndata as ad
   import scanpy as sc
   import pandas as pd
   import matplotlib.pyplot as plt

   # Load results
   adata = ad.read_h5ad('sc_allele_counts.h5ad')
   results = pd.read_csv('ai_results_CD4_T.tsv', sep='\t')

   # Add cell type annotations
   celltype_map = pd.read_csv('barcode_celltype_map.tsv', sep='\t',
                               header=None, names=['barcode', 'celltype'])
   adata.obs = adata.obs.join(celltype_map.set_index('barcode'))

   # Plot allelic ratio for a specific region
   region_id = 'chr1:1000000-1001000'
   if region_id in adata.var_names:
       # Get ref and alt counts for this region
       ref_counts = adata[:, region_id].X.toarray().flatten()
       adata.obs['allelic_ratio'] = ref_counts / (ref_counts + 0.1)  # Add pseudocount

       sc.pl.violin(adata, 'allelic_ratio', groupby='celltype')

Compare cell types
~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import pandas as pd
   import matplotlib.pyplot as plt

   # Load differential AI results
   diff_ai = pd.read_csv('CD4_vs_CD8_diff_ai.tsv', sep='\t')

   # Volcano plot
   plt.figure(figsize=(10, 8))
   colors = ['red' if fdr < 0.05 else 'gray' for fdr in diff_ai['fdr']]
   plt.scatter(diff_ai['ratio_diff'], -np.log10(diff_ai['pvalue']),
               c=colors, alpha=0.5)
   plt.xlabel('Difference in Allelic Ratio (CD4_T - CD8_T)')
   plt.ylabel('-log10(p-value)')
   plt.axhline(-np.log10(0.05), linestyle='--', color='gray')
   plt.axvline(0, linestyle='--', color='gray')
   plt.title('Differential Allelic Imbalance: CD4 vs CD8 T cells')
   plt.savefig('diff_ai_volcano.png', dpi=150)

Troubleshooting
---------------

Very sparse matrix
~~~~~~~~~~~~~~~~~~

**Problem**: Most cells have 0 counts for most SNPs

**Solutions**:

- This is normal for single-cell data
- Aggregate by cell type for reliable statistics
- Focus on regions with higher coverage
- Consider pseudobulk analysis

.. code-block:: python

   # Check sparsity
   import numpy as np
   sparsity = (adata.X == 0).sum() / adata.X.size
   print(f"Matrix sparsity: {sparsity:.1%}")

Memory issues
~~~~~~~~~~~~~

**Problem**: Out of memory with large datasets

**Solutions**:

- Process chromosomes separately
- Subsample cells for initial analysis
- Use sparse matrix operations

.. code-block:: bash

   # Process by chromosome
   for chr in {1..22}; do
       bcftools view -r chr${chr} variants.vcf.gz -O z -o chr${chr}.vcf.gz
       wasp2-count count-variants-sc sample.bam chr${chr}.vcf.gz \
           cell_barcodes.txt --out_file counts_chr${chr}.h5ad
   done

No significant results
~~~~~~~~~~~~~~~~~~~~~~

**Problem**: No regions pass significance threshold

**Possible causes**:

- Insufficient coverage (need more cells or deeper sequencing)
- Overly stringent threshold (try FDR < 0.1)
- Biological: cell types may have similar allelic patterns

**Solutions**:

- Lower ``--min`` threshold (but interpret cautiously)
- Combine similar cell types
- Focus on highly accessible regions (ATAC) or expressed genes (RNA)

Advanced Topics
---------------

Phased analysis
~~~~~~~~~~~~~~~

If you have phased genotypes (from family data or phasing software):

.. code-block:: bash

   wasp2-analyze find-imbalance-sc \
       sc_allele_counts.h5ad \
       barcode_celltype_map.tsv \
       --phased \
       --out_file ai_results_phased

Integration with other tools
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Export to Seurat (R)**:

.. code-block:: python

   import anndata as ad

   adata = ad.read_h5ad('sc_allele_counts.h5ad')

   # Export counts matrix
   import scipy.io
   scipy.io.mmwrite('counts.mtx', adata.X.T)

   # Export metadata
   adata.var.to_csv('features.tsv', sep='\t')
   adata.obs.to_csv('barcodes.tsv', sep='\t')

**Integration with ArchR/Signac** (scATAC-seq):

Results can be added as metadata to your ArchR or Signac objects
for joint analysis with accessibility and motif data.

Next Steps
----------

- :doc:`/user_guide/analysis` - Statistical model details
- :doc:`/faq` - Common questions
- :doc:`/troubleshooting` - Detailed troubleshooting guide
