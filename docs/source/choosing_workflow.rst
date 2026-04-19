Choosing the Right Workflow
===========================

WASP2 supports four major data types. Use this guide to find your workflow.

.. list-table::
   :header-rows: 1
   :widths: 25 25 25 25

   * - Data Type
     - Input
     - Goal
     - Start Here
   * - **Bulk RNA-seq**
     - BAM + phased VCF
     - Allele-specific expression (ASE)
     - :doc:`tutorials/bulk_workflow`
   * - **Bulk ATAC-seq**
     - BAM + phased VCF
     - Allele-specific chromatin accessibility
     - :doc:`tutorials/bulk_workflow`
   * - **scRNA-seq (10x)**
     - Cell Ranger BAM + VCF + barcodes
     - Per-cell or per-cell-type ASE
     - :doc:`tutorials/single_cell_workflow`
   * - **scATAC-seq (10x)**
     - Fragments/BAM + VCF + barcodes
     - Single-cell allelic imbalance in ATAC peaks
     - :doc:`tutorials/single_cell_workflow`

Decision Flowchart
------------------

**Step 1: What sequencing assay did you run?**

* RNA-seq → go to Step 2
* ATAC-seq → go to Step 3

**Step 2: Bulk or single-cell RNA-seq?**

* Bulk RNA-seq → :doc:`tutorials/bulk_workflow`
* 10x Chromium scRNA-seq → :doc:`tutorials/single_cell_workflow`
* Other single-cell protocol → see :doc:`user_guide/single_cell`

**Step 3: Bulk or single-cell ATAC-seq?**

* Bulk ATAC-seq → :doc:`tutorials/bulk_workflow` (use BED peak file as ``--region``)
* 10x scATAC-seq → :doc:`tutorials/single_cell_workflow`

Do I Need to Run the WASP Remapping Step?
------------------------------------------

The remapping step (``wasp2-map``) corrects **reference mapping bias** — reads
carrying the alternative allele are harder to map than reference-allele reads,
causing false-positive imbalance signals.

**You need remapping if:**

* Your BAM was aligned with a standard aligner (BWA-MEM, STAR, HISAT2, bowtie2)
* You want the most rigorous allele-specific analysis
* You are studying regions near known variants (high variant density)

**You can skip remapping if:**

* Your BAM was already produced by an unbiased pipeline
* You are doing a quick exploratory analysis
* You are using simulated or controlled data

See :doc:`user_guide/mapping` for the full remapping workflow.

What VCF Do I Need?
--------------------

WASP2 requires a **phased VCF** with heterozygous variants for the sample(s)
you are analyzing. Supported formats:

* VCF/BCF (bgzip + tabix indexed)
* PLINK2 PGEN files (with ``.pvar`` + ``.psam``)

See :doc:`user_guide/counting` for VCF format requirements and examples using
bcftools to subset and phase.

Nextflow Pipelines
------------------

If you prefer a managed workflow with automatic parallelization, containerization,
and output publishing, use the bundled Nextflow pipelines instead of the CLI:

* **nf-rnaseq** — bulk RNA-seq allele-specific expression
* **nf-atacseq** — bulk ATAC-seq allele-specific chromatin accessibility
* **nf-scatac** — single-cell ATAC-seq allelic imbalance
* **nf-outrider** — outlier expression detection with allele-aware correction

See the pipeline-specific documentation for samplesheet format and parameter
reference.
