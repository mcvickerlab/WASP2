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

The development branch tests Python 3.10, 3.11, 3.12, and 3.13. Available
published wheels depend on the release and platform; see the release files.

**Why do I get an error about missing samtools/bcftools/bedtools?**

The PyPI wheel bundles the Rust extension and htslib but not the system binaries.
Install them via conda (``mamba install -c bioconda samtools bcftools bedtools``)
or your system package manager.

Input Data
----------

**Do I need phased genotypes?**

No. Counting uses bases at heterozygous variant positions and does not require
phasing. Bulk analysis supports phased and unphased likelihoods; select
``--phased`` only when using suitable phased genotypes. Single-cell analysis
checks genotypes automatically unless explicitly set to ``--unphased``.

**What VCF formats does WASP2 support?**

* VCF or BCF (bgzip-compressed + tabix-indexed: ``.vcf.gz`` + ``.tbi``)
* PLINK2 PGEN format (``.pgen`` + ``.pvar`` + ``.psam``)

Multi-sample VCFs are supported; use ``-s SAMPLE_ID`` to specify the target sample.

**Can I use an unphased VCF?**

Yes. REF and ALT counts remain defined by the variant alleles. Use unphased
inference when haplotype phase is unavailable; this marginalizes possible
orientations rather than treating arbitrary assignments as known phase.

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

Supply a region **file**, such as ``peaks.bed`` or ``genes.gtf``, to restrict
and annotate overlapping variants. ``--region`` is not a ``chr:start-end``
interval string. Single-cell region input uses BED peaks.

Single-Cell
-----------

**What single-cell chemistries are supported?**

The validated workflow is scATAC-seq with a barcode-tagged BAM. The current
counter does not implement UMI-aware molecule counting for scRNA-seq.
See :doc:`user_guide/single_cell` for barcode files and supported outputs.

**Do I need Cell Ranger output?**

No. Counting requires a sorted, indexed BAM with cell barcodes in a tag
(default ``CB``), donor variants, and a one-column barcode whitelist.
Analysis uses a separate two-column barcode-to-group mapping. Fragment
coordinates alone cannot provide REF/ALT base evidence.

**How do I get per-cell-type results?**

Run WASP2 on the full BAM to get per-cell allele counts, then use the output
with your cell type annotations in Python (AnnData/Scanpy) to aggregate by
cell type. See :doc:`tutorials/single_cell_workflow` for an example.

Output and Results
------------------

**What does the p-value in the output represent?**

The p-value comes from a likelihood ratio test comparing the beta-binomial model
under allelic imbalance vs. the null model of balanced expression. The model
estimates overdispersion from the supplied counts.

**What FDR threshold should I use?**

The standard threshold is FDR < 0.05. For discovery analyses you may want
FDR < 0.1. For validation or follow-up experiments, consider FDR < 0.01.
See :doc:`user_guide/analysis` for the BH procedure and the NaN-propagation warning.

**My output has very few significant sites. What's wrong?**

Common causes:

* Low coverage at heterozygous sites (inspect depth and the minimum-count threshold)
* Too few heterozygous variants in the VCF
* VCF and BAM use different chromosome naming conventions (``chr1`` vs ``1``)

**My output has too many significant sites (inflated FDR).**

This typically means mapping bias is driving the signal. Run the WASP remapping
step (``wasp2-map``) before counting. See :doc:`user_guide/mapping`.

**For ATAC-seq, do I need to use WASP-remapped BAMs?**

Yes. WASP2 counting applies only the unmapped filter (see :doc:`methods/mapping_filter`
"Canonical Filter Contract"); it does **not** correct reference mapping bias
on its own. You must run ``wasp2-map make-reads`` + re-alignment +
``wasp2-map filter-remapped`` first, then pass the resulting
final filtered BAM to ``wasp2-count``. Counting on raw BWA output
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
