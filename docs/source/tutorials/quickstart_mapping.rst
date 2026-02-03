Quickstart: WASP Mapping Filter
================================

Learn WASP2's mapping bias correction in 5 minutes.

.. contents:: Contents
   :local:
   :depth: 2

Overview
--------

**Goal:** Understand and apply the WASP mapping filter to remove reference bias from your alignment data.

**Time:** ~5 minutes to read, ~30 minutes to run on typical data

**Prerequisites:**

* WASP2 installed (``pip install wasp2``)
* Aligned BAM file (coordinate-sorted)
* VCF file with heterozygous variants

The Problem: Reference Mapping Bias
-----------------------------------

When reads are aligned to a reference genome, there's an inherent asymmetry:

.. code-block:: text

   Reference:  ...ACGT[A]CGTA...  (reference allele: A)
   Read (ref): ...ACGT[A]CGTA...  → Perfect match (0 mismatches)
   Read (alt): ...ACGT[G]CGTA...  → 1 mismatch penalty

**Result**: Reads carrying the alternate allele are more likely to:

- Fail to map entirely
- Map with lower quality scores
- Map to incorrect locations

This causes **inflated reference allele counts**, leading to false positive ASE signals.

The Solution: WASP Remap-and-Filter
-----------------------------------

WASP corrects this by testing whether each read would map identically
regardless of which allele it carries:

1. **Identify**: Find reads overlapping heterozygous SNPs
2. **Swap**: Create versions with alleles swapped (ref→alt, alt→ref)
3. **Remap**: Align swapped reads with the same aligner
4. **Filter**: Keep only reads that map to the **same location** after swapping

After filtering, the probability of mapping is equal for both alleles:

.. math::

   P(\text{map} | \text{ref allele}) = P(\text{map} | \text{alt allele})

Quick Workflow
--------------

Step 1: Create Swapped Reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Identify reads overlapping heterozygous SNPs and generate allele-swapped versions:

.. code-block:: bash

   wasp2-map make-reads sample.bam variants.vcf.gz \
     --samples SAMPLE1 \
     --out-dir wasp_output/

This produces (where ``sample`` is your BAM file prefix):

* ``wasp_output/sample_to_remap.bam``: Original reads needing remapping
* ``wasp_output/sample_keep.bam``: Reads not overlapping variants (kept as-is)
* ``wasp_output/sample_swapped_alleles_r1.fq``: Allele-swapped read 1
* ``wasp_output/sample_swapped_alleles_r2.fq``: Allele-swapped read 2
* ``wasp_output/sample_wasp_data_files.json``: Metadata for filter step

Step 2: Remap Swapped Reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Critical**: Use the **same aligner and parameters** as your original mapping!

.. code-block:: bash

   # Example with BWA (replace 'sample' with your BAM file prefix)
   bwa mem -M -t 8 genome.fa \
     wasp_output/sample_swapped_alleles_r1.fq \
     wasp_output/sample_swapped_alleles_r2.fq | \
     samtools sort -o wasp_output/sample_remapped.bam -

   samtools index wasp_output/sample_remapped.bam

Using different alignment parameters will invalidate the WASP correction.

Step 3: Filter Remapped Reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The WASP filter compares original and remapped positions:

.. code-block:: bash

   wasp2-map filter-remapped \
     wasp_output/sample_to_remap.bam \
     wasp_output/sample_remapped.bam \
     wasp_output/sample_wasp_filtered.bam

Understanding Filter Statistics
-------------------------------

The WASP filter reports three key metrics:

.. table:: WASP Filter Metrics
   :widths: 20 50 30

   ================= ============================================ ==============
   Metric            Description                                  Typical Value
   ================= ============================================ ==============
   **Kept reads**    Reads that passed the filter                 90-99%
   **Removed (moved)** Reads that mapped to different locations   1-8%
   **Removed (missing)** Reads that failed to remap               <1%
   ================= ============================================ ==============

Interpreting Filter Rates
~~~~~~~~~~~~~~~~~~~~~~~~~

* **95-99% kept**: Good - typical for most data types
* **90-95% kept**: Acceptable - may indicate difficult regions
* **<90% kept**: Investigate - check data quality or variant calls

Before/After Example
--------------------

At a site with mapping bias:

.. table:: Example: Before and After WASP
   :widths: 30 35 35

   =============== =============== ===============
   Metric          Before WASP     After WASP
   =============== =============== ===============
   Reference reads 150             95
   Alternate reads 80              85
   Ref fraction    0.65            0.53
   =============== =============== ===============

The biased site (0.65 ref fraction) is corrected to near-balanced (0.53).

Complete Workflow Script
------------------------

.. code-block:: bash

   #!/bin/bash
   set -e

   # Input files
   BAM="sample.bam"
   VCF="variants.vcf.gz"
   SAMPLE="SAMPLE1"
   GENOME="genome.fa"
   OUTDIR="wasp_output"

   # Extract BAM prefix (filename without .bam extension)
   PREFIX=$(basename $BAM .bam)

   mkdir -p $OUTDIR

   # Step 1: Create allele-swapped reads
   echo "Step 1: Creating swapped reads..."
   wasp2-map make-reads $BAM $VCF \
     --samples $SAMPLE \
     --out-dir $OUTDIR/

   # Step 2: Remap with same aligner
   echo "Step 2: Remapping swapped reads..."
   bwa mem -M -t 8 $GENOME \
     $OUTDIR/${PREFIX}_swapped_alleles_r1.fq \
     $OUTDIR/${PREFIX}_swapped_alleles_r2.fq | \
     samtools sort -o $OUTDIR/${PREFIX}_remapped.bam -
   samtools index $OUTDIR/${PREFIX}_remapped.bam

   # Step 3: Filter biased reads
   echo "Step 3: Filtering biased reads..."
   wasp2-map filter-remapped \
     $OUTDIR/${PREFIX}_to_remap.bam \
     $OUTDIR/${PREFIX}_remapped.bam \
     $OUTDIR/${PREFIX}_wasp_filtered.bam

   # Step 4: Merge with non-overlapping reads
   echo "Step 4: Merging final BAM..."
   samtools merge -f $OUTDIR/${PREFIX}_final.bam \
     $OUTDIR/${PREFIX}_wasp_filtered.bam \
     $OUTDIR/${PREFIX}_keep.bam
   samtools index $OUTDIR/${PREFIX}_final.bam

   echo "Done! WASP-filtered BAM: $OUTDIR/${PREFIX}_final.bam"

Rust Acceleration
-----------------

WASP2 includes a high-performance Rust backend that accelerates the filter step:

.. table:: Performance Comparison
   :widths: 30 35 35

   ============= =============== ===============
   Dataset Size  Python          Rust
   ============= =============== ===============
   1M reads      ~5 minutes      ~30 seconds
   10M reads     ~50 minutes     ~5 minutes
   100M reads    ~8 hours        ~50 minutes
   ============= =============== ===============

The Rust backend is used automatically when available.

Next Steps
----------

After WASP filtering:

1. **Count alleles** on the filtered BAM:

   .. code-block:: bash

      wasp2-count count-variants wasp_filtered.bam variants.vcf

2. **Analyze allelic imbalance**:

   .. code-block:: bash

      wasp2-analyze find-imbalance counts.tsv

See Also
--------

* :doc:`/user_guide/mapping` - Detailed mapping module documentation
* :doc:`/methods/mapping_filter` - Algorithm details and mathematics
* :doc:`/tutorials/scrna_seq` - Single-cell RNA-seq workflow

Summary
-------

.. table:: Key Takeaways
   :widths: 25 75

   ============ ===============================================
   Concept      Key Point
   ============ ===============================================
   **Problem**  Reference bias inflates ref allele counts
   **Solution** WASP remap-and-filter removes biased reads
   **Workflow** make-reads → remap → filter-remapped
   **Expected** 90-99% reads pass filter
   **Result**   Unbiased allele counts for ASE analysis
   ============ ===============================================
