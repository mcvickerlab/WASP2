Frequently Asked Questions
==========================

Installation
------------

**Which install method should I use?**

For most users: ``mamba install -c conda-forge -c bioconda wasp2`` (Bioconda).
This installs WASP2 and all dependencies (samtools, bcftools, bedtools) in one step.

Use ``pip install wasp2`` if you are on a system without conda or want a specific
Python environment. You will need to install samtools, bcftools, and bedtools separately.

Use Docker or Singularity on HPC clusters or when you need full reproducibility.

**What Python versions are supported?**

Python 3.10, 3.11, 3.12, and 3.13. Pre-built wheels are available for all four
on Linux (x86_64, aarch64) and macOS (Intel, Apple Silicon).

**Why do I get an error about missing samtools/bcftools/bedtools?**

The PyPI wheel bundles the Rust extension and htslib but not the system binaries.
Install them via conda (``mamba install -c bioconda samtools bcftools bedtools``)
or your system package manager.

Input Data
----------

**Do I need phased genotypes?**

Yes. WASP2 assigns reads to haplotypes using phased heterozygous variants. Without
phase information, WASP2 cannot distinguish which allele a read came from. Use
WhatsHap, SHAPEIT4, or Eagle2 to phase your VCF before running WASP2.

**What VCF formats does WASP2 support?**

* VCF or BCF (bgzip-compressed + tabix-indexed: ``.vcf.gz`` + ``.tbi``)
* PLINK2 PGEN format (``.pgen`` + ``.pvar`` + ``.psam``)

Multi-sample VCFs are supported; use ``-s SAMPLE_ID`` to specify the target sample.

**Can I use an unphased VCF?**

The counting step (``wasp2-count``) will still run but the allele assignments will
be arbitrary. The statistical results will have reduced power and increased false
positives. Always use phased genotypes when possible.

**My BAM doesn't have read groups. Will WASP2 work?**

Yes, for counting. Read groups are not required for allele counting. For the
remapping step (``wasp2-map``), the sample ID is needed to look up variants in
a multi-sample VCF — pass it explicitly with ``-s SAMPLE_ID``.

Running WASP2
-------------

**How long does each step take?**

Typical runtimes on a single core for a 30× whole-genome BAM (~100M reads):

* ``wasp2-map make-reads``: 2–4 hours
* Re-alignment (external): depends on aligner
* ``wasp2-map filter-remapped``: 1–2 hours
* ``wasp2-count count-variants``: 30–60 minutes
* ``wasp2-analyze find-imbalance``: < 5 minutes

Use the Nextflow pipelines for automatic parallelization across chromosomes/samples.

**Can I run WASP2 on multiple samples at once?**

Yes. WASP2 CLI processes one sample at a time; run multiple samples in parallel
with a job scheduler (SLURM, PBS) or use the Nextflow pipelines which handle
parallelization automatically.

**What is the ``--region`` flag for?**

Restrict counting to a specific genomic region (e.g., ``chr1:1000000-2000000``).
Useful for testing on a subset of data or for chromosome-level parallelization.

Single-Cell
-----------

**What single-cell chemistries are supported?**

All 10x Genomics Chromium chemistries (scRNA v1/v2/v3, scATAC v1/v2) and any
other protocol with a cell barcode tag in the BAM (CB tag by default). See
:doc:`user_guide/single_cell` for barcode format details.

**Do I need Cell Ranger output?**

No, but it is the most common input. WASP2 needs:

* A BAM with cell barcodes in a BAM tag (default: ``CB``)
* A whitelist of valid barcodes (optional but recommended)
* A phased VCF

Any aligner that produces CB-tagged BAMs will work (STARsolo, Alevin-fry, etc.).

**How do I get per-cell-type results?**

Run WASP2 on the full BAM to get per-cell allele counts, then use the output
with your cell type annotations in Python (AnnData/Scanpy) to aggregate by
cell type. See :doc:`tutorials/single_cell_workflow` for an example.

Output and Results
------------------

**What does the p-value in the output represent?**

The p-value comes from a likelihood ratio test comparing the beta-binomial model
under allelic imbalance vs. the null model of balanced expression. The test is
calibrated for the overdispersion typical of RNA-seq count data.

**What FDR threshold should I use?**

The standard threshold is FDR < 0.05. For discovery analyses you may want
FDR < 0.1. For validation or follow-up experiments, consider FDR < 0.01.
See :doc:`user_guide/analysis` for the BH procedure and the NaN-propagation warning.

**My output has very few significant sites. What's wrong?**

Common causes:

* Low coverage at heterozygous sites (increase ``--min_count``)
* Too few heterozygous variants in the VCF
* VCF and BAM use different chromosome naming conventions (``chr1`` vs ``1``)
* VCF is not phased

**My output has too many significant sites (inflated FDR).**

This typically means mapping bias is driving the signal. Run the WASP remapping
step (``wasp2-map``) before counting. See :doc:`user_guide/mapping`.

**For ATAC-seq, do I need to use WASP-remapped BAMs?**

Yes. WASP2 counting applies only the unmapped filter (see :doc:`methods/mapping_filter`
"Canonical Filter Contract"); it does **not** correct reference mapping bias
on its own. You must run ``wasp2-map make-reads`` + re-alignment +
``wasp2-map filter-remapped`` first, then pass the resulting
``*_wasp_filt_rmdup.bam`` to ``wasp2-count``. Counting on raw BWA output
leaves reference bias uncorrected — reads carrying the alt allele are
systematically under-represented.

The same requirement applies to RNA-seq and scATAC-seq. The only difference
is the aligner used in the re-alignment step (STAR for RNA, BWA for ATAC).

Troubleshooting
---------------

**I get "chromosome not found" errors.**

VCF and BAM must use the same chromosome naming convention. If your VCF uses
``chr1`` and your BAM uses ``1`` (or vice versa), use ``bcftools annotate --rename-chrs``
to harmonize the VCF.

**The Rust extension fails to load.**

This happens if the wheel was built for a different platform or Python version.
Try reinstalling: ``pip install --force-reinstall wasp2``. If building from source,
run ``pixi run verify`` to rebuild.

**WASP2 runs but produces an empty counts file.**

Check that:

* The BAM is coordinate-sorted and indexed (``.bai`` file present)
* The VCF overlaps the regions in your BAM
* The sample name passed with ``-s`` matches a sample in the VCF

Use ``bcftools query -l variants.vcf.gz`` to list VCF sample names.
