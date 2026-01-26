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

* 16 nucleotides followed by ``-N`` suffix (e.g., ``CACCCAAGTGAGTTGG-1``)
* The suffix indicates the GEM well (``-1`` for single sample, ``-2``, ``-3``, etc. for aggregated samples)
* Barcodes are from the 10X whitelist (~3 million for v3 chemistry, ~737,000 for v2)

**Chemistry Versions:**

.. list-table::
   :header-rows: 1
   :widths: 20 25 55

   * - Chemistry
     - Barcode Length
     - Notes
   * - 10X v2
     - 16 bp
     - ~737,000 valid barcodes, older whitelist
   * - 10X v3/v3.1
     - 16 bp
     - ~3.5 million valid barcodes, improved capture
   * - 10X Multiome
     - 16 bp
     - Same as v3, paired ATAC+GEX

**PBMC Example (10X v3):**

.. code-block:: text

   AAACCCAAGAAACACT-1	B_cell
   AAACCCAAGAAAGCGA-1	CD4_T_cell
   AAACCCAAGAACAACT-1	CD8_T_cell
   AAACCCAAGAACCAAG-1	Monocyte
   AAACCCAAGAACGATA-1	NK_cell

**Multi-Sample Aggregated Example:**

When using Cell Ranger ``aggr`` to combine multiple samples, barcodes are distinguished by suffix:

.. code-block:: text

   AAACCCAAGAAACACT-1	B_cell	sample1
   AAACCCAAGAAACTGT-1	B_cell	sample1
   AAACCCAAGAAACACT-2	B_cell	sample2
   AAACCCAAGAAACTGT-2	B_cell	sample2
   AAACCCAAGAAACACT-3	CD4_T_cell	sample3

.. note::

   For multi-sample experiments, WASP2 uses only the first two columns (barcode, cell_type).
   The third column (sample origin) is optional metadata for your reference.

Barcode Format Validation
~~~~~~~~~~~~~~~~~~~~~~~~~

Before running WASP2, validate your barcode file format:

.. code-block:: bash

   # Check file structure (should show TAB separator)
   head -5 barcodes.tsv | cat -A
   # Expected output (^I = TAB):
   # AAACCCAAGAAACACT-1^IB_cell$

   # Verify barcode format matches 10X pattern
   head -1 barcodes.tsv | cut -f1 | grep -E '^[ACGT]{16}-[0-9]+$'
   # Should return the barcode if valid

   # Count barcodes per cell type
   cut -f2 barcodes.tsv | sort | uniq -c | sort -rn

   # Check for common issues
   # 1. No header row (first line should be a barcode, not "barcode")
   head -1 barcodes.tsv

   # 2. Correct delimiter (TAB not space/comma)
   file barcodes.tsv  # Should mention "ASCII text"

**Python Validation Script:**

.. code-block:: python

   import re

   def validate_10x_barcode_file(filepath):
       """Validate 10X scRNA-seq barcode file format."""
       pattern = re.compile(r'^[ACGT]{16}-\d+$')
       errors = []
       i = 0

       with open(filepath) as f:
           for i, line in enumerate(f, 1):
               parts = line.rstrip('\n').split('\t')

               # Check column count
               if len(parts) < 1:
                   errors.append(f"Line {i}: Empty line")
                   continue

               barcode = parts[0]

               # Check barcode format
               if not pattern.match(barcode):
                   errors.append(f"Line {i}: Invalid barcode format '{barcode}'")

               # Check for header (common mistake)
               if i == 1 and barcode.lower() in ('barcode', 'cell_barcode', 'cb'):
                   errors.append(f"Line 1: Appears to be a header row, remove it")

       if errors:
           print(f"Found {len(errors)} errors:")
           for err in errors[:10]:  # Show first 10
               print(f"  {err}")
           return False
       else:
           print(f"Validation passed: {i} barcodes")
           return True

   # Usage
   validate_10x_barcode_file('barcodes.tsv')

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

Common Format Variations
~~~~~~~~~~~~~~~~~~~~~~~~

**Cell Ranger Raw Barcodes:**

.. code-block:: bash

   # Extract filtered barcodes (single-column, add cell types later)
   zcat cellranger_output/outs/filtered_feature_bc_matrix/barcodes.tsv.gz > barcodes_raw.txt

**Barcode Suffix Handling:**

Some tools strip the ``-1`` suffix. Ensure BAM and barcode file match:

.. code-block:: bash

   # Compare formats
   samtools view sample.bam | head -1000 | grep -o 'CB:Z:[^\t]*' | cut -d: -f3 | head
   cut -f1 barcodes.tsv | head

   # Add suffix if missing
   awk -F'\t' '{print $1"-1\t"$2}' barcodes_no_suffix.tsv > barcodes.tsv

Single-Cell CLI Usage
---------------------

Count Alleles
~~~~~~~~~~~~~

.. code-block:: bash

   wasp2-count count-variants-sc \
     sample.bam \
     variants.vcf.gz \
     barcodes.tsv \
     --region peaks.bed \
     --samples NA12878 \
     --out_file allele_counts.h5ad

Analyze Imbalance
~~~~~~~~~~~~~~~~~

.. code-block:: bash

   wasp2-analyze find-imbalance-sc \
     allele_counts.h5ad \
     barcodes.tsv \
     --sample NA12878 \
     --out_file imbalance_results.tsv

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

Example Files
-------------

WASP2 includes example barcode files in the ``tests/data/`` directory:

* ``barcodes_10x_scrna.tsv`` - Standard PBMC cell types (B_cell, CD4_T_cell, etc.)
* ``barcodes_example.tsv`` - Brain tissue cell types (Neurons, Astrocytes, etc.)
* ``barcodes_10x_multi_sample.tsv`` - Multi-sample aggregated experiment with -1, -2, -3 suffixes
* ``barcodes_10x_hierarchical.tsv`` - Hierarchical cell type naming (T_cell.CD4.Naive, etc.)

These files can be used as templates or for testing your WASP2 installation.

Comparative Analysis
--------------------

After detecting allelic imbalance within individual cell populations, you can compare
imbalance **between** groups to identify cell-type-specific or condition-dependent
regulatory variation.

**Quick example:**

.. code-block:: bash

   # Compare imbalance between two cell types
   wasp2-analyze compare-imbalance \
     allele_counts.h5ad \
     barcode_celltype_map.tsv \
     --groups "excitatory_neurons,inhibitory_neurons" \
     --sample SAMPLE_ID \
     --phased

This identifies genomic regions where allelic imbalance differs significantly between
the specified groups, using a likelihood ratio test with FDR correction.

For comprehensive coverage of comparative analysis, see:

* :doc:`/tutorials/comparative_imbalance` - Detailed comparative analysis tutorial
* :doc:`analysis` - Statistical methods for comparative imbalance

See Also

--------

* :doc:`/tutorials/scrna_seq` - Complete 10X scRNA-seq tutorial
* :doc:`analysis` - Statistical analysis methods
* :doc:`counting` - General allele counting
