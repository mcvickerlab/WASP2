Single-Cell Analysis
====================

Overview
--------

WASP2 provides specialized tools for allele-specific analysis in single-cell RNA-seq (scRNA-seq) data. This guide covers the barcode file format requirements and single-cell-specific workflows.

Barcode File Format
-------------------

WASP2 uses a two-column TSV (tab-separated values) format for barcode files. This format maps cell barcodes to cell type annotations.

Format Specification
~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

   BARCODE<TAB>CELLTYPE

**Requirements:**

* No header row
* Tab-separated (``\t``) delimiter
* Column 1: Cell barcode (string)
* Column 2: Cell type annotation (string)

**Example:**

.. code-block:: text

   CACCCAAGTGAGTTGG-1	Oligodendrocytes
   GCTTAAGCCGCGGCAT-1	Oligodendrocytes
   GTCACGGGTGGCCTAG-1	Endothelial
   AACCATGGTCACCTAA-1	Microglia
   TGAGCCGAGAAACGCC-1	Astrocytes

10X Genomics Barcodes
~~~~~~~~~~~~~~~~~~~~~

10X Chromium barcodes follow a specific format:

* 16 nucleotides followed by ``-1`` suffix (e.g., ``CACCCAAGTGAGTTGG-1``)
* The ``-1`` suffix indicates the GEM well
* Barcodes are from the 10X whitelist (~737,000 valid barcodes for v3 chemistry)

Cell Ranger Output
------------------

When using Cell Ranger output, barcodes can be found in:

.. code-block:: text

   cellranger_output/
   └── outs/
       └── filtered_feature_bc_matrix/
           └── barcodes.tsv.gz

This file contains only the barcode column. To create a WASP2-compatible barcode file, you need to add cell type annotations from your downstream analysis.

Generating Barcode Files
------------------------

From Seurat (R)
~~~~~~~~~~~~~~~

After clustering and cell type annotation in Seurat:

.. code-block:: r

   # Assuming 'seurat_obj' has cell type labels in metadata
   library(Seurat)

   # Extract barcodes and cell types
   barcode_df <- data.frame(
     barcode = colnames(seurat_obj),
     cell_type = seurat_obj$cell_type  # Your annotation column
   )

   # Write TSV without header
   write.table(
     barcode_df,
     file = "barcodes.tsv",
     sep = "\t",
     quote = FALSE,
     row.names = FALSE,
     col.names = FALSE
   )

From Scanpy (Python)
~~~~~~~~~~~~~~~~~~~~

After clustering and cell type annotation in Scanpy:

.. code-block:: python

   import pandas as pd

   # Assuming 'adata' has cell type labels in obs
   barcode_df = pd.DataFrame({
       'barcode': adata.obs_names,
       'cell_type': adata.obs['cell_type']  # Your annotation column
   })

   # Write TSV without header
   barcode_df.to_csv(
       'barcodes.tsv',
       sep='\t',
       header=False,
       index=False
   )

Simple Barcode List
~~~~~~~~~~~~~~~~~~~

If you only need to filter by barcodes without cell type annotation, you can use a single-column file:

.. code-block:: text

   CACCCAAGTGAGTTGG-1
   GCTTAAGCCGCGGCAT-1
   GTCACGGGTGGCCTAG-1

Single-Cell CLI Usage
---------------------

Count Alleles
~~~~~~~~~~~~~

.. code-block:: bash

   wasp2-count count-variants-sc \
     sample.bam \
     variants.vcf.gz \
     --barcodes barcodes.tsv \
     --region peaks.bed \
     --samples NA12878 \
     --output allele_counts.h5ad

Analyze Imbalance
~~~~~~~~~~~~~~~~~

.. code-block:: bash

   wasp2-analyze find-imbalance-sc \
     allele_counts.h5ad \
     --sample NA12878 \
     --groups cell_type \
     --min-count 5 \
     --output imbalance_results.tsv

Output Format
-------------

The single-cell counting module outputs an AnnData (``.h5ad``) file containing:

**Layers:**

* ``X``: Total allele counts (ref + alt + other)
* ``ref``: Reference allele counts
* ``alt``: Alternate allele counts
* ``other``: Other allele counts

**Observations (obs):**

* SNP information (chrom, pos, ref, alt)
* Aggregate counts per SNP

**Variables (var):**

* Cell barcodes

**Unstructured (uns):**

* Sample information
* Count statistics
* Feature-SNP mapping (if regions provided)

Best Practices
--------------

Quality Filtering
~~~~~~~~~~~~~~~~~

* Filter low-quality cells before generating barcode file
* Remove doublets and dead cells
* Use cells with sufficient UMI counts (>1000 for most protocols)

Cell Type Annotation
~~~~~~~~~~~~~~~~~~~~

* Use consistent cell type naming (no spaces, special characters)
* Consider hierarchical annotations (e.g., ``T_cell``, ``CD4_T_cell``)
* Document your annotation sources and markers

Barcode Matching
~~~~~~~~~~~~~~~~

* Ensure barcodes match exactly (including ``-1`` suffix)
* Verify barcode format matches BAM file CB tags
* Check for barcode format differences between tools

See Also
--------

* :doc:`/tutorials/scrna_seq` - Complete 10X scRNA-seq tutorial
* :doc:`analysis` - Statistical analysis methods
* :doc:`counting` - General allele counting
