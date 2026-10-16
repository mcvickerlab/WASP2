Choosing the Right Workflow
===========================

WASP2 supports bulk RNA-seq, bulk ATAC-seq, and barcode-based scATAC-seq.
Use this guide to find your workflow.

.. list-table::
   :header-rows: 1
   :widths: 25 25 25 25

   * - Data Type
     - Input
     - Goal
     - Start Here
   * - **Bulk RNA-seq**
     - BAM + donor VCF
     - Allele-specific expression (ASE)
     - :doc:`tutorials/bulk_workflow`
   * - **Bulk ATAC-seq**
     - BAM + donor VCF
     - Allele-specific chromatin accessibility
     - :doc:`tutorials/bulk_workflow`
   * - **scATAC-seq (10x)**
     - CB-tagged BAM + VCF + barcodes
     - Single-cell allelic imbalance in ATAC peaks
     - :doc:`tutorials/single_cell_workflow`

Decision Flowchart
------------------

**Step 1: What sequencing assay did you run?**

* RNA-seq → go to Step 2
* ATAC-seq → go to Step 3

**Step 2: Bulk or single-cell RNA-seq?**

* Bulk RNA-seq → :doc:`tutorials/bulk_workflow`
* Single-cell RNA-seq → the counter does not implement UMI-aware molecule
  counting; the scATAC tutorial is not an scRNA-seq workflow.

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

WASP2 uses heterozygous variants for the selected samples. Phasing is
optional: use the phased likelihood only with suitable phased genotypes.
Supported formats:

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
* **nf-scatac** — scATAC coverage workflow; see its documentation for the
  current BAM/inference limitations
* **nf-outrider** — outlier expression detection with allele-aware correction

See the pipeline-specific documentation for samplesheet format and parameter
reference.
