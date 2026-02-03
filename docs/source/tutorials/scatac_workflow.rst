Single-Cell ATAC-seq Workflow
=============================

This tutorial provides a workflow for detecting allelic imbalance in single-cell ATAC-seq data from 10x Genomics.

.. note::

   **Estimated Time**: ~30 minutes

Overview
--------

**Goal**: Identify genomic regions with allelic imbalance in chromatin accessibility at single-cell resolution.

**Input Data**:

* 10x Cell Ranger ATAC output (fragments/BAM + barcodes)
* Phased VCF with heterozygous variants
* Cell type annotations

Tutorial Sections
-----------------

1. Loading 10x scATAC Data
~~~~~~~~~~~~~~~~~~~~~~~~~~

Cell Ranger ATAC outputs needed:

.. code-block:: text

   cellranger_output/outs/
   ├── fragments.tsv.gz          # Fragment overlap counting
   ├── possorted_bam.bam         # Allele-specific counting
   ├── peaks.bed                 # Region restriction
   └── filtered_peak_bc_matrix/
       └── barcodes.tsv.gz       # Filtered barcodes

2. Cell Barcode Handling
~~~~~~~~~~~~~~~~~~~~~~~~

10x barcode format: 16 nucleotides + ``-N`` suffix (e.g., ``AAACGAACAGTCAGTT-1``)

.. code-block:: bash

   # Verify BAM and barcode file match
   samtools view your.bam | head -1000 | grep -o 'CB:Z:[^\t]*' | head
   head barcodes.tsv

3. Counting Strategies
~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - Aspect
     - Per-Cell
     - Pseudo-Bulk
   * - Resolution
     - Single-cell
     - Cell population
   * - Power
     - Low (sparse)
     - High (aggregated)
   * - Use case
     - Outlier cells
     - Population imbalance

**Recommendation**: Use pseudo-bulk for most scATAC experiments.

.. code-block:: bash

   # Count alleles at heterozygous variants
   wasp2-count count-variants-sc \
     possorted_bam.bam \
     variants.vcf.gz \
     barcodes_celltype.tsv \
     --region peaks.bed \
     --samples SAMPLE_ID \
     --out_file allele_counts.h5ad

**Output**: ``allele_counts.h5ad`` - AnnData with layers: ``X``, ``ref``, ``alt``, ``other``

4. Statistical Considerations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

WASP2 handles sparse data through:

* **Dispersion model**: Accounts for overdispersion in allele counts
* **Minimum count filters**: ``--min 10`` ensures sufficient data
* **FDR correction**: Benjamini-Hochberg for multiple testing
* **Outlier removal**: ``-z 3`` filters CNV/mapping artifacts

**Key parameters**:

* ``--phased``: Use phased genotype information (requires ``0|1`` or ``1|0`` format in VCF)

5. Visualization
~~~~~~~~~~~~~~~~

The notebook includes functions for:

* Allelic ratio heatmaps
* Volcano plots
* Cell type comparison heatmaps

6. Cell-Type-Specific Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Step 1: Find imbalance within cell types
   wasp2-analyze find-imbalance-sc \
     allele_counts.h5ad \
     barcodes_celltype.tsv \
     --sample SAMPLE_ID --phased --min 10 -z 3
   # Output: ai_results_<celltype>.tsv per cell type

   # Step 2: Compare between cell types
   wasp2-analyze compare-imbalance \
     allele_counts.h5ad \
     barcodes_celltype.tsv \
     --sample SAMPLE_ID --groups "CellTypeA,CellTypeB" --phased
   # Output: ai_results_<celltype1>_<celltype2>.tsv

**Output columns**: region, ref_count, alt_count, p_value, fdr_pval, effect_size

Troubleshooting
---------------

No Barcodes Matched
~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Add -1 suffix if missing
   awk -F'\t' '{print $1"-1\t"$2}' barcodes_no_suffix.tsv > barcodes.tsv

Memory Issues
~~~~~~~~~~~~~

Process chromosomes separately with ``--region peaks_chr1.bed``.

Low Power
~~~~~~~~~

* Merge similar cell types
* Use pseudo-bulk aggregation
* Ensure phased genotypes

See Also
--------

* :doc:`/tutorials/scrna_seq` - 10X scRNA-seq tutorial
* :doc:`/tutorials/comparative_imbalance` - Comparative analysis
* :doc:`/user_guide/single_cell` - Data format reference
