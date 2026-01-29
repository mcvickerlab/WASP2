Allele Counting Algorithm
=========================

This document describes how WASP2 assigns reads to reference and alternate
alleles at heterozygous variant sites.

.. contents:: Contents
   :local:
   :depth: 2

Overview
--------

The allele counting algorithm forms the foundation of allele-specific analysis.
For each heterozygous SNP, WASP2 examines aligned reads and counts how many
support the reference allele, alternate allele, or neither.

Biological Rationale
--------------------

In a diploid organism with a heterozygous site (genotype A/G), reads originating
from each chromosome should carry the corresponding allele. Under the null
hypothesis of no allelic imbalance, we expect equal representation of both
alleles:

.. math::

   E[\text{ref\_count}] = E[\text{alt\_count}] = \frac{N}{2}

where :math:`N` is the total number of reads covering the variant.

Deviations from this expectation may indicate:

- **Allele-specific expression (ASE)**: Differential transcription between alleles
- **Allele-specific chromatin accessibility**: Differential regulatory activity
- **Allele-specific binding**: Differential protein-DNA interactions
- **Technical artifacts**: Mapping bias, amplification bias

Algorithm Details
-----------------

Position-Based Alignment
^^^^^^^^^^^^^^^^^^^^^^^^

For each variant position, WASP2 queries the BAM file to retrieve all reads
overlapping that genomic coordinate:

1. **Coordinate lookup**: Use BAM index to efficiently retrieve reads at position
2. **CIGAR parsing**: Walk through the CIGAR string to find the query position
   corresponding to the reference position
3. **Base extraction**: Extract the nucleotide at the query position

.. code-block:: python

   # Simplified pseudocode
   for read in bam.fetch(chrom, pos, pos + 1):
       query_pos = find_aligned_position(read, pos)
       if query_pos is not None:
           base = read.query_sequence[query_pos]
           if base == ref_allele:
               ref_count += 1
           elif base == alt_allele:
               alt_count += 1
           else:
               other_count += 1

CIGAR-Aware Position Finding
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The alignment between query (read) and reference positions is not always 1:1
due to insertions, deletions, and soft clipping. WASP2 uses the aligned pairs
from the CIGAR string:

.. table:: CIGAR Operations and Position Mapping
   :widths: 20 40 40

   ==================== =========================== =======================
   CIGAR Operation      Consumes Reference          Consumes Query
   ==================== =========================== =======================
   M (match/mismatch)   Yes                         Yes
   I (insertion)        No                          Yes
   D (deletion)         Yes                         No
   S (soft clip)        No                          Yes
   H (hard clip)        No                          No
   N (skip)             Yes                         No
   ==================== =========================== =======================

For a read with CIGAR ``50M2D30M``:

- Positions 0-49 in the reference align to positions 0-49 in the query
- Positions 50-51 in the reference have no query alignment (deletion)
- Positions 52-81 in the reference align to positions 50-79 in the query

Handling Edge Cases
^^^^^^^^^^^^^^^^^^^

**Deletions spanning the variant**
  If the variant position falls within a deletion, the read cannot be assigned
  to either allele and is counted as "other".

**Insertions adjacent to the variant**
  Insertions can shift the query position. The algorithm correctly handles
  this by using aligned pairs rather than simple arithmetic.

**Soft-clipped reads**
  Soft-clipped bases at read ends do not consume reference positions but
  are present in the query sequence. The algorithm accounts for this.

Rust Acceleration
-----------------

WASP2 uses a Rust-accelerated BAM counter for performance:

.. code-block:: python

   from wasp2_rust import BamCounter

   counter = BamCounter("sample.bam")
   regions = [("chr1", 12345, "A", "G"), ("chr1", 12400, "C", "T")]
   counts = counter.count_alleles(regions, min_qual=0, threads=4)

The Rust implementation provides:

- **Parallel processing**: Count multiple regions simultaneously
- **Memory efficiency**: Stream through BAM without loading all reads
- **htslib integration**: Direct access to BAM index for efficient queries

Performance Characteristics
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. table:: Counting Performance
   :widths: 30 35 35

   =================== ======================== ========================
   Dataset Size        Python (pysam)           Rust (wasp2_rust)
   =================== ======================== ========================
   10,000 SNPs         ~45 seconds              ~5 seconds
   100,000 SNPs        ~7 minutes               ~40 seconds
   1,000,000 SNPs      ~70 minutes              ~6 minutes
   =================== ======================== ========================

*Benchmarks on typical ATAC-seq data with 100M reads, single thread.*

Quality Filtering
-----------------

Optional base quality filtering can exclude low-confidence base calls:

.. math::

   Q = -10 \log_{10}(P_{\text{error}})

By default, WASP2 uses ``min_qual=0`` (no filtering) to match legacy WASP
behavior. For stringent analysis, ``min_qual=20`` (1% error rate) is recommended.

Output Format
-------------

The counting algorithm produces a table with the following columns:

.. table:: Count Output Columns
   :widths: 20 20 60

   ============= ======== ============================================
   Column        Type     Description
   ============= ======== ============================================
   chrom         string   Chromosome name
   pos           int      1-based genomic position
   ref           string   Reference allele
   alt           string   Alternate allele
   ref_count     int      Reads supporting reference allele
   alt_count     int      Reads supporting alternate allele
   other_count   int      Reads with neither allele (errors, indels)
   region        string   Associated region (peak, gene) if provided
   ============= ======== ============================================

Region Assignment
-----------------

When a BED file of regions (peaks, genes) is provided, variants are associated
with overlapping regions:

1. **Intersection**: Use bedtools or coitrees to find variant-region overlaps
2. **Aggregation**: Sum counts across all variants within each region
3. **Statistical testing**: Perform imbalance analysis at the region level

This aggregation increases statistical power for detecting allelic imbalance,
as individual SNPs often have low coverage.

See Also
--------

- :doc:`statistical_models` - Statistical testing of allele counts
