Counting Module
===============

Overview
--------

``wasp2-count`` counts reads supporting reference and alternate alleles at
variant positions in BAM files.

It provides three commands:

* ``count-variants`` for bulk data
* ``count-cohort`` for a locked bulk cohort with one final BAM per donor
* ``count-variants-sc`` for single-cell data with ``CB``-tagged barcodes

Bulk Counting
-------------

Basic usage:

.. code-block:: bash

   wasp2-count count-variants sample.bam variants.vcf.gz --out_file counts.tsv

With sample filtering and region annotation:

.. code-block:: bash

   wasp2-count count-variants \
     sample.bam \
     variants.vcf.gz \
     --samples SAMPLE1 \
     --region genes.gtf \
     --out_file counts.tsv

Supported region files:

* BED
* MACS2 ``narrowPeak`` / ``broadPeak``
* GTF
* GFF3

For GTF/GFF3 inputs, WASP2 derives interval annotations from feature rows and
defaults to ``gene`` features when present.

Useful options:

* ``--samples`` / ``-s``: select het sites for one or more samples
* ``--region`` / ``-r``: restrict/annotate variants by overlapping regions
* ``--gene_feature``: choose the GTF/GFF3 feature type
* ``--gene_attribute``: choose the GTF/GFF3 attribute used as the feature ID
* ``--gene_parent``: choose the parent/grouping attribute for gene annotations
* ``--use_region_names``: prefer region names instead of coordinate strings
* ``--include-indels``: count indels in addition to SNPs

Output columns always include:

* ``chrom``
* ``pos`` or ``pos0`` / ``pos`` depending on input path
* ``ref``
* ``alt``
* ``ref_count``
* ``alt_count``
* ``other_count``

When sample filtering is active, genotype columns are included. When region
annotation is active, region or gene columns are included as well.

Locked Cohort Counting
----------------------

``count-cohort`` makes the donor-to-VCF-to-BAM relationship explicit. Its
manifest must contain exactly three tab-separated columns:

.. code-block:: text

   donor_id	vcf_sample	bam
   donor1	VCF_sample_1	/path/to/donor1.bam
   donor2	VCF_sample_2	/path/to/donor2.bam

Each row must identify one final merged, WASP-filtered, deduplicated, indexed
BAM. The cohort VCF must also be indexed and contain every ``vcf_sample``.
The locked cohort path applies WASP2's standard biallelic-SNV selection:
multi-allelic records and indels are excluded and that policy is recorded in
``count_manifest.json``.

Independent SNV counts restricted to ATAC peaks:

.. code-block:: bash

   wasp2-count count-cohort \
     donors.tsv variants.vcf.gz cohort_snv_counts \
     --unit snv --regions peaks.bed

Peak-feature counts:

.. code-block:: bash

   wasp2-count count-cohort \
     donors.tsv variants.vcf.gz cohort_feature_counts \
     --unit feature --regions peaks.bed

The command refuses to overwrite an existing output directory. It claims the
directory atomically and publishes ``count_manifest.json`` last as the commit
marker after hashing and rechecking every input. The count manifest is the
authoritative entry point for downstream donor-aware analysis.

Single-Cell ATAC Counting
-------------------------

Single-cell counting is designed for **scATAC-seq** data. It requires a BAM
with ``CB`` tags and a positional barcode file containing one barcode per line.

.. code-block:: bash

   wasp2-count count-variants-sc \
     sc_atac.bam \
     variants.vcf.gz \
     barcodes.tsv \
     --samples sample1 \
     --feature peaks.bed \
     --out_file allele_counts.h5ad

Important points:

* ``barcodes.tsv`` is a positional argument, not ``--barcode_map``
* ``--feature`` and ``--region`` are aliases on the single-cell command
* Accepts BED and MACS2 peak files (GTF/GFF3 are supported only by the bulk ``count-variants`` command)

The output is an AnnData ``.h5ad`` file with:

* sparse count layers for ``ref``, ``alt``, and ``other``
* variant metadata in ``adata.obs``
* barcode names in ``adata.var_names``
* feature-to-variant mapping in ``adata.uns["feature"]`` when annotations are used

Examples
--------

Count variants without regional annotation:

.. code-block:: bash

   wasp2-count count-variants \
     filtered.bam \
     variants.vcf.gz \
     --samples SAMPLE1 \
     --out_file counts.tsv

Count variants inside peaks:

.. code-block:: bash

   wasp2-count count-variants \
     filtered.bam \
     variants.vcf.gz \
     --samples SAMPLE1 \
     --region peaks.bed \
     --out_file counts_peaks.tsv

Count variants inside genes:

.. code-block:: bash

   wasp2-count count-variants \
     filtered.bam \
     variants.vcf.gz \
     --samples SAMPLE1 \
     --region genes.gtf \
     --gene_feature gene \
     --gene_attribute gene_id \
     --out_file counts_genes.tsv

Next Steps
----------

* :doc:`analysis` for statistical testing of allelic imbalance
* :doc:`/user_guide/single_cell` for barcode grouping and single-cell workflows
