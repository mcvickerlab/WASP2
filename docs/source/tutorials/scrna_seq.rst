10X scRNA-seq Tutorial
======================

This tutorial walks through a complete WASP2 workflow for detecting allele-specific expression in 10X Genomics single-cell RNA-seq data.

Overview
--------

**Goal:** Identify genes with allele-specific expression (ASE) in different cell types from 10X Chromium scRNA-seq data.

**Input Data:**

* Cell Ranger output (BAM + filtered barcodes)
* Phased VCF file with heterozygous variants
* Cell type annotations from Seurat/Scanpy

Prerequisites
-------------

**Software:**

* WASP2 (``pip install wasp2``)
* Cell Ranger output (v3+)
* R with Seurat or Python with Scanpy

**Data:**

* Aligned BAM file with cell barcodes (CB tag)
* Phased VCF for your sample
* Completed cell type annotation

Step 1: Prepare Input Data
--------------------------

Cell Ranger Output Structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After running Cell Ranger, your output directory contains:

.. code-block:: text

   cellranger_output/
   └── outs/
       ├── possorted_genome_bam.bam
       ├── possorted_genome_bam.bam.bai
       └── filtered_feature_bc_matrix/
           ├── barcodes.tsv.gz
           ├── features.tsv.gz
           └── matrix.mtx.gz

The BAM file contains cell barcodes in the ``CB`` tag:

.. code-block:: bash

   # View CB tags in BAM
   samtools view possorted_genome_bam.bam | head -1 | tr '\t' '\n' | grep CB
   # Output: CB:Z:AAACCCAAGAAACACT-1

Step 2: Generate Barcode File
-----------------------------

From Seurat Analysis
~~~~~~~~~~~~~~~~~~~~

After running Seurat clustering and annotation:

.. code-block:: r

   library(Seurat)

   # Load your analyzed Seurat object
   seurat_obj <- readRDS("seurat_analyzed.rds")

   # Check available metadata columns
   head(seurat_obj@meta.data)

   # Extract barcodes and cell types
   barcode_df <- data.frame(
     barcode = colnames(seurat_obj),
     cell_type = seurat_obj$celltype_annotation  # Your annotation column
   )

   # Preview the data
   head(barcode_df)
   #>                 barcode    cell_type
   #> 1 AAACCCAAGAAACACT-1       B_cell
   #> 2 AAACCCAAGAAACTGT-1       B_cell
   #> 3 AAACCCAAGAAAGCGA-1   CD4_T_cell

   # Write barcode file (no header, tab-separated)
   write.table(
     barcode_df,
     file = "barcodes_celltype.tsv",
     sep = "\t",
     quote = FALSE,
     row.names = FALSE,
     col.names = FALSE
   )

From Scanpy Analysis
~~~~~~~~~~~~~~~~~~~~

After running Scanpy clustering and annotation:

.. code-block:: python

   import scanpy as sc
   import pandas as pd

   # Load your analyzed AnnData object
   adata = sc.read_h5ad("scanpy_analyzed.h5ad")

   # Check available annotations
   print(adata.obs.columns)

   # Extract barcodes and cell types
   barcode_df = pd.DataFrame({
       'barcode': adata.obs_names,
       'cell_type': adata.obs['leiden_annotation']  # Your annotation column
   })

   # Preview the data
   print(barcode_df.head())
   #                   barcode    cell_type
   # 0  AAACCCAAGAAACACT-1       B_cell
   # 1  AAACCCAAGAAACTGT-1       B_cell
   # 2  AAACCCAAGAAAGCGA-1   CD4_T_cell

   # Write barcode file (no header, tab-separated)
   barcode_df.to_csv(
       'barcodes_celltype.tsv',
       sep='\t',
       header=False,
       index=False
   )

Verify Barcode Format
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Check barcode file format
   head barcodes_celltype.tsv
   # AAACCCAAGAAACACT-1	B_cell
   # AAACCCAAGAAACTGT-1	B_cell

   # Count cells per type
   cut -f2 barcodes_celltype.tsv | sort | uniq -c | sort -rn
   #   1500 CD4_T_cell
   #   1200 B_cell
   #    800 Monocyte
   #    ...

Step 3: Count Allele-Specific Reads
-----------------------------------

Run the single-cell allele counting:

.. code-block:: bash

   wasp2-count count-variants-sc \
     cellranger_output/outs/possorted_genome_bam.bam \
     phased_variants.vcf.gz \
     --barcodes barcodes_celltype.tsv \
     --region genes.gtf \
     --samples SAMPLE_ID \
     --output allele_counts.h5ad

**Parameters:**

* ``--barcodes``: Your barcode file with cell type annotations
* ``--region``: Gene annotation file (GTF/GFF) or peak file (BED)
* ``--samples``: Sample ID matching VCF sample column
* ``--output``: Output AnnData file

Step 4: Analyze Allelic Imbalance
---------------------------------

Cell-Type-Specific Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Analyze imbalance within each cell type:

.. code-block:: bash

   wasp2-analyze find-imbalance-sc \
     allele_counts.h5ad \
     --sample SAMPLE_ID \
     --groups cell_type \
     --min-count 5 \
     --output imbalance_by_celltype.tsv

**Output columns:**

* ``region``: Gene or genomic region
* ``cell_type``: Cell type from barcode file
* ``ref_count``: Total reference allele counts
* ``alt_count``: Total alternate allele counts
* ``p_value``: Statistical significance
* ``fdr_pval``: FDR-corrected p-value
* ``effect_size``: Log2 fold change

Compare Between Cell Types
~~~~~~~~~~~~~~~~~~~~~~~~~~

Find differential allelic imbalance between cell types:

.. code-block:: bash

   wasp2-analyze compare-imbalance \
     allele_counts.h5ad \
     --groups "CD4_T_cell,CD8_T_cell" \
     --output differential_imbalance.tsv

Step 5: Interpret Results
-------------------------

Load and explore results in Python:

.. code-block:: python

   import pandas as pd
   import matplotlib.pyplot as plt

   # Load results
   results = pd.read_csv('imbalance_by_celltype.tsv', sep='\t')

   # Filter significant results (FDR < 0.05)
   significant = results[results['fdr_pval'] < 0.05]
   print(f"Found {len(significant)} significant ASE events")

   # Top genes per cell type
   top_genes = (significant
       .groupby('cell_type')
       .apply(lambda x: x.nsmallest(10, 'fdr_pval'))
       .reset_index(drop=True))

   print(top_genes[['region', 'cell_type', 'effect_size', 'fdr_pval']])

   # Visualize effect sizes
   fig, ax = plt.subplots(figsize=(10, 6))
   significant.boxplot(column='effect_size', by='cell_type', ax=ax)
   plt.title('Allelic Imbalance by Cell Type')
   plt.ylabel('Log2 Fold Change (Ref/Alt)')
   plt.savefig('ase_by_celltype.png')

Example Output
--------------

.. code-block:: text

   region          cell_type    ref_count  alt_count  fdr_pval   effect_size
   ENSG00000123456 B_cell       245        89         0.001      1.46
   ENSG00000234567 CD4_T_cell   156        312        0.003     -1.00
   ENSG00000345678 Monocyte     423        198        0.012      1.09

Troubleshooting
---------------

No Cells Matched
~~~~~~~~~~~~~~~~

If you see "0 barcodes matched":

.. code-block:: bash

   # Check BAM barcode format
   samtools view your.bam | head -1000 | grep -o 'CB:Z:[^\t]*' | head

   # Compare with your barcode file
   head barcodes.tsv

   # Common issues:
   # - Missing -1 suffix in barcode file
   # - Barcode file has header (should not)
   # - Different barcode versions (v2 vs v3)

Low Read Counts
~~~~~~~~~~~~~~~

Single-cell data is sparse. Consider:

* Using pseudo-bulk aggregation by cell type
* Lowering ``--min-count`` threshold
* Focusing on highly expressed genes

Memory Issues
~~~~~~~~~~~~~

For large datasets:

.. code-block:: bash

   # Process chromosomes separately
   for chr in chr{1..22}; do
     wasp2-count count-variants-sc \
       sample.bam \
       variants.vcf.gz \
       --barcodes barcodes.tsv \
       --region genes.gtf \
       --chrom $chr \
       --output counts_${chr}.h5ad
   done

Next Steps
----------

* Integrate with eQTL databases (GTEx, eQTLGen)
* Correlate ASE with gene expression levels
* Validate top hits with allele-specific primers
* Compare across conditions or timepoints

See Also
--------

* :doc:`/user_guide/single_cell` - Barcode file format reference
* :doc:`/user_guide/analysis` - Statistical methods
* `Seurat <https://satijalab.org/seurat/>`_ - R toolkit for scRNA-seq
* `Scanpy <https://scanpy.readthedocs.io/>`_ - Python toolkit for scRNA-seq
