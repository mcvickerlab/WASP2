Bulk Workflow (RNA-seq / ATAC-seq)
===================================

End-to-end walkthrough of WASP2 on bulk sequencing data: WASP mapping filter,
allele counting at regions, and the beta-binomial LRT for allelic imbalance.
The same pipeline works for RNA-seq at genes (GTF) and ATAC-seq at peaks
(BED); data-type differences are called out where they matter.

Inputs:

- Aligned BAM (coordinate-sorted, indexed)
- VCF with the donor's heterozygous variants (phased for RNA-seq haplotypes)
- Region file: GTF/GFF for RNA-seq genes, BED for ATAC-seq peaks

Step 1 — WASP mapping filter
----------------------------

Reference alleles map more easily than alternate alleles, biasing ref counts
upward. The WASP remap-and-filter step removes reads whose mapping position
depends on which allele they carry. The test is "would this read map to the
same location if we swapped ref↔alt at the SNP?" — only reads that do are
kept. After filtering,

.. math::

   P(\text{map} \mid \text{ref}) \approx P(\text{map} \mid \text{alt})

under the usual assumption that the aligner's scoring is deterministic and
position-dependent (see :doc:`/methods/mapping_filter` and [vandeGeijn2015]_).

.. code-block:: bash

   # Step 1a: produce allele-swapped FASTQ + BAMs
   wasp2-map make-reads sample.bam variants.vcf.gz \
     --samples SAMPLE1 \
     --out_dir wasp_output

   # Step 1b: remap with the SAME aligner and parameters used originally
   bwa mem -M -t 8 genome.fa \
     wasp_output/sample_swapped_alleles_r1.fq \
     wasp_output/sample_swapped_alleles_r2.fq \
     | samtools sort -o wasp_output/sample_remapped.bam -
   samtools index wasp_output/sample_remapped.bam

   # Step 1c: drop reads whose swapped version mapped elsewhere
   wasp2-map filter-remapped \
     wasp_output/sample_remapped.bam \
     --wasp_data_json wasp_output/sample_wasp_data_files.json \
     --out_bam wasp_output/sample_wasp_filtered.bam

   # Step 1d: merge with reads that never overlapped a SNP
   samtools merge -f wasp_output/sample_final.bam \
     wasp_output/sample_wasp_filtered.bam \
     wasp_output/sample_keep.bam
   samtools index wasp_output/sample_final.bam

.. important::

   The remap step **must** use the same aligner, version, and parameters as
   the original alignment. A different scoring function invalidates the
   filter contract.

Filter rates are **approximately** 1–10% of SNP-overlapping reads dropped
(developer experience; varies by sequencing platform, aligner, and variant
density). Outlier rates should be investigated — they typically signal
CNVs, aligner-version mismatch, or variant-call errors rather than a WASP
failure.

Step 2 — Count alleles
----------------------

**RNA-seq at genes** (GTF):

.. code-block:: bash

   wasp2-count count-variants \
     wasp_output/sample_final.bam \
     variants.vcf.gz \
     --samples SAMPLE1 \
     --region genes.gtf \
     --out_file allele_counts.tsv

**ATAC-seq at peaks** (BED):

.. code-block:: bash

   wasp2-count count-variants \
     wasp_output/sample_final.bam \
     variants.vcf.gz \
     --samples SAMPLE1 \
     --region peaks.bed \
     --out_file allele_counts.tsv

Output columns: ``chrom``, ``pos``, ``ref``, ``alt``, ``ref_count``,
``alt_count``, ``other_count``, plus ``region`` (peak name) or
``gene_id``/``gene_name`` (GTF). See :doc:`/user_guide/counting` for region
filtering and multi-sample options.

Step 3 — Test for allelic imbalance
-----------------------------------

.. code-block:: bash

   wasp2-analyze find-imbalance \
     allele_counts.tsv \
     --min 10 \
     --pseudocount 1 \
     --phased \
     --out_file ai_results.tsv

The beta-binomial LRT tests :math:`\mu = 0.5` against :math:`\mu \neq 0.5`
per region, correcting for overdispersion and multiple testing via
Benjamini–Hochberg (:doc:`/methods/statistical_models`).

``--phased`` applies only when the VCF GT field uses ``0|1`` / ``1|0``; pass
it off for unphased data.

Output columns: ``region``, ``num_snps``, aggregated ``ref_count`` /
``alt_count``, ``pval``, ``fdr_pval``, ``mu``, ``effect_size``
(:math:`\log_2(\text{ref/alt})`).

Filtering and interpretation
----------------------------

.. code-block:: bash

   # FDR < 0.05 and |log2FC| > 1
   awk -F'\t' 'NR==1 || ($5 < 0.05 && ($6 > 1 || $6 < -1))' \
     ai_results.tsv > significant_ase.tsv

Monoallelic expression (ref fraction > 0.9 or < 0.1) is a common hallmark of
genomic imprinting on RNA-seq and of strong cis-regulatory variants on ATAC.

eQTL integration (RNA-seq only)
-------------------------------

.. code-block:: python

   import pandas as pd

   ase = pd.read_csv('ai_results.tsv', sep='\t')
   eqtl = pd.read_csv('gtex_eqtl.tsv', sep='\t')

   merged = ase.merge(eqtl, left_on='region', right_on='gene_id')

   # ASE effect_size > 0 => REF allele more expressed.
   # eQTL slope > 0      => ALT allele increases expression.
   # Concordance = opposite signs.
   merged['concordant'] = (merged['effect_size'] > 0) != (merged['slope'] > 0)

Troubleshooting
---------------

**Few SNPs found.** Confirm the VCF sample column matches ``--samples``, the
VCF has heterozygous genotypes for this donor, and BAM+VCF use the same
reference build.

**Low power / few significant results.** Increase depth, lower ``--min``
(with caution), or aggregate at a coarser region level.

**Too many significant results.** Confirm the WASP filter was applied, check
for batch effects, CNVs, or aligner-version drift, and tighten FDR.

See Also
--------

- :doc:`/user_guide/mapping` — full mapping CLI reference
- :doc:`/user_guide/counting` — counting CLI reference
- :doc:`/user_guide/analysis` — analysis CLI reference
- :doc:`/methods/mapping_filter` — canonical WASP filter contract
- :doc:`/methods/statistical_models` — the LRT and beta-binomial model
- :doc:`/tutorials/single_cell_workflow` — single-cell RNA-seq / ATAC-seq
- :doc:`/tutorials/comparative_imbalance` — comparing groups
