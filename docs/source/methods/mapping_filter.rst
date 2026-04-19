WASP Mapping Bias Correction
============================

This document describes the WASP algorithm for correcting reference mapping bias
in allele-specific analysis.

.. contents:: Contents
   :local:
   :depth: 2

The Problem: Reference Mapping Bias
-----------------------------------

Standard read aligners map sequencing reads against a reference genome. When a
read originates from the alternate allele at a heterozygous site, it carries a
mismatch relative to the reference:

.. code-block:: text

   Reference:  ...ACGT[A]CGTA...
   Read (ref): ...ACGT[A]CGTA...  → Maps perfectly (0 mismatches)
   Read (alt): ...ACGT[G]CGTA...  → Maps with 1 mismatch

This asymmetry causes **reference mapping bias**: reads carrying the reference
allele are more likely to map successfully and with higher quality, leading to
inflated reference allele counts.

The effect is particularly pronounced when:

- The variant is near other polymorphisms (haplotype effects)
- The read has low overall quality
- The aligner uses strict mismatch penalties
- The region has repetitive sequence

Uncorrected mapping bias inflates reference allele counts, causing false
positive ASE signals and biased effect sizes in QTL mapping.

The WASP Algorithm
------------------

WASP (WASP Allele-Specific Pipeline) [vandeGeijn2015]_ corrects mapping bias
through a **remap-and-filter** strategy:

Algorithm Overview
^^^^^^^^^^^^^^^^^^

1. **Identify overlapping reads**: Find reads that overlap heterozygous SNPs
2. **Swap alleles**: For each overlapping read, create a version with the
   alternate allele swapped to the reference (and vice versa)
3. **Remap**: Align the swapped reads to the reference genome
4. **Filter**: Keep only reads that map to the **same location** after swapping

Mathematical Justification
^^^^^^^^^^^^^^^^^^^^^^^^^^

Let :math:`M(r)` be the mapping location of read :math:`r`, and let :math:`r'`
be the allele-swapped version of :math:`r`.

A read passes the WASP filter if and only if:

.. math::

   M(r) = M(r')

This criterion ensures that the read would have mapped identically regardless
of which allele it carried, eliminating differential mappability.

Under the assumption that the aligner's scoring is deterministic and
depends only on the read sequence and reference coordinate — the standard
case for seed-and-extend aligners at typical parameter settings — the
mapping probability becomes approximately equal for either allele after
filtering:

.. math::

   P(\text{map} \mid \text{ref allele}) \approx P(\text{map} \mid \text{alt allele})

See [vandeGeijn2015]_ §Methods for the original argument. The equality is
exact under deterministic mapping and approximate for aligners with
stochastic tie-breaking or heuristic multi-mapper resolution.

Implementation Details
----------------------

Step 1: Create Reads for Remapping
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

WASP2 identifies reads overlapping variants using interval trees (coitrees)
for efficient coordinate queries:

.. code-block:: bash

   wasp2-map make-reads sample.bam variants.vcf --samples SAMPLE1

For each read overlapping a heterozygous site:

1. Extract the original read sequence
2. Identify all variant positions within the read
3. Generate haplotype combinations with swapped alleles
4. Write swapped reads to FASTQ for remapping

**Haplotype Generation**

When a read overlaps multiple heterozygous sites, all combinations must be
tested. For :math:`n` het sites, there are :math:`2^n` haplotypes:

.. code-block:: text

   Read overlaps 2 het sites: A/G and C/T

   Original:    ...A...C...
   Haplotype 1: ...G...C...  (swap first)
   Haplotype 2: ...A...T...  (swap second)
   Haplotype 3: ...G...T...  (swap both)

WASP2 caps the number of haplotypes per read (default: 64) to prevent
combinatorial explosion with highly polymorphic regions.

Step 2: Remap with Original Aligner
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The swapped reads must be remapped using the **same aligner and parameters**
as the original mapping:

.. code-block:: bash

   bwa mem -M genome.fa swapped_r1.fq swapped_r2.fq | \
     samtools view -bS - > remapped.bam
   samtools sort -o remapped.sorted.bam remapped.bam
   samtools index remapped.sorted.bam

**Critical**: Using different alignment parameters will invalidate the
WASP correction.

Step 3: Filter Remapped Reads
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The WASP filter compares original and remapped positions:

.. code-block:: bash

   wasp2-map filter-remapped \
     remapped.bam \
     --wasp_data_json sample_wasp_data_files.json \
     --out_bam output.bam

A read passes if:

1. The remapped read exists (didn't fail to map)
2. The mapping position matches the original within tolerance
3. For paired-end reads, both mates satisfy the above

Canonical Filter Contract (v1.4.1+)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

WASP2 filter-remapped and count-variants apply **only the unmapped filter**
(SAM flag ``0x4``). Secondary, supplementary, duplicate, QC-fail, and
not-proper-pair reads are **not** filtered by WASP2 itself.

This is deliberate, following the canonical WASP contract from
van de Geijn et al. 2015:

    *"This program does not perform filtering of reads based on mappability.
    It is assumed that the input BAM files are filtered appropriately prior
    to calling this script."* — ``bmvdgeijn/WASP`` ``CHT/bam2h5.py``

If your pipeline needs those defensive filters, apply them upstream:

.. code-block:: bash

   samtools view -F 0x904 in.bam -o filtered.bam
   # 0x904 = secondary (0x100) | supplementary (0x800) | duplicate (0x400)

Earlier WASP2 releases (v1.2.0–v1.4.0) applied six SAM flag filters
internally. See ``CHANGELOG.md`` for details.

**Counting step.** The same canonical contract applies to the allele
counting step (``wasp2-count count-variants``): only the unmapped filter
(``0x4``) is applied. All other SAM-flag decisions are the caller's
responsibility — pre-filter the BAM upstream if defensive filtering is
needed.

**Same-Locus Test**

For SNPs, exact position matching is required:

.. math::

   |M_{\text{original}} - M_{\text{remapped}}| = 0

For indels, a small tolerance (slop) accommodates alignment ambiguity:

.. math::

   |M_{\text{original}} - M_{\text{remapped}}| \leq \text{slop}

Paired-End Considerations
^^^^^^^^^^^^^^^^^^^^^^^^^

For paired-end reads, WASP2 requires both mates to pass:

- Both mates must remap successfully
- Both mates must map to the same location
- Insert size must remain consistent

This is more stringent than single-end filtering but ensures the read pair
as a unit is unbiased.

Rust Acceleration
-----------------

WASP2 implements the filtering step in Rust for performance:

.. code-block:: python

   from wasp2_rust import filter_bam_wasp

   kept, filtered, total = filter_bam_wasp(
       to_remap_bam="original.bam",
       remapped_bam="remapped.bam",
       remap_keep_bam="output.bam",
       threads=4,
       same_locus_slop=0,  # Exact matching for SNPs
   )

The Rust implementation provides:

- **Parallel BAM I/O**: Multi-threaded reading and writing
- **Streaming comparison**: Memory-efficient position matching
- **htslib integration**: Native BAM format support

Performance
-----------

The Rust implementation is substantially faster than a pure-Python filter
path (roughly an order of magnitude on coordinate-sorted BAMs of 1–100M
reads in internal testing). Exact throughput depends on storage type,
thread count, and read length; run on your own hardware for accurate
numbers.

Expected Filter Rates
^^^^^^^^^^^^^^^^^^^^^

Filter rates vary with data type, variant density, aligner, and tolerance
settings. As a rough developer-experience guide only:

- RNA-seq: on the order of a few percent to ~15%; higher near splice
  junctions and in indel-dense regions.
- ATAC-seq / ChIP-seq: on the order of a few percent.
- WGS: typically lower than RNA-seq at comparable variant density.

Outlier filter rates (well above these ranges) usually indicate a mismatch
between the original and remap aligner versions/parameters, CNV-rich
regions, or errors in the variant call set — not a WASP failure.

Limitations and Considerations
------------------------------

Indel Handling
^^^^^^^^^^^^^^

The original WASP algorithm was designed for SNPs. Indels present challenges:

- Alignment ambiguity at indel boundaries
- Multiple valid alignments for the same read
- Gap penalties interact with variant detection

WASP2 supports indel mode with configurable parameters:

.. code-block:: bash

   wasp2-map make-reads sample.bam variants.vcf \
     --include-indels --max-indel-len 10

Structural Variants
^^^^^^^^^^^^^^^^^^^

WASP is not designed for structural variants (large deletions, inversions,
translocations). These require specialized methods.

Reference Panel Quality
^^^^^^^^^^^^^^^^^^^^^^^

The effectiveness of WASP depends on having accurate variant calls:

- Missing variants leave mapping bias uncorrected
- False positive variants cause unnecessary read filtering
- Imputation errors can introduce systematic biases

Use high-quality variant calls from the same sample or a well-matched
reference panel.

Pipeline Integration
--------------------

The typical WASP2 workflow:

.. code-block:: bash

   # Step 1: Initial mapping
   bwa mem -M genome.fa reads_r1.fq reads_r2.fq | \
     samtools sort -o sample.bam -

   # Step 2: Create swapped reads
   wasp2-map make-reads sample.bam variants.vcf \
     --samples SAMPLE1 --out_dir wasp_temp

   # Step 3: Remap swapped reads (SAME ALIGNER!)
   bwa mem -M genome.fa wasp_temp/swapped_r1.fq wasp_temp/swapped_r2.fq | \
     samtools sort -o remapped.bam -

   # Step 4: Filter biased reads
   wasp2-map filter-remapped \
     wasp_temp/to_remap.bam remapped.bam wasp_filtered.bam

   # Step 5: Count alleles on filtered BAM
   wasp2-count count-variants wasp_filtered.bam variants.vcf \
     --samples SAMPLE1 --regions peaks.bed

See Also
--------

- :doc:`/user_guide/analysis` — beta-binomial LRT + FDR CLI and defaults

References
----------

.. [vandeGeijn2015] van de Geijn B, McVicker G, Gilad Y, Pritchard JK (2015).
   WASP: allele-specific software for robust molecular quantitative trait
   locus discovery. *Nature Methods* 12:1061-1063.
   https://doi.org/10.1038/nmeth.3582
