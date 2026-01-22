Data Formats
============

WASP2 supports multiple input and output formats for flexibility and performance.

.. contents:: Format Categories
   :local:
   :depth: 1

.. toctree::
   :maxdepth: 2

   input_formats
   output_formats

Quick Reference
---------------

Input Formats
~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 15 25 30 30

   * - Type
     - Format
     - Extension
     - Notes
   * - Reads
     - BAM
     - .bam
     - Must be sorted and indexed
   * - Variants
     - VCF
     - .vcf, .vcf.gz
     - Standard, indexed
   * - Variants
     - BCF
     - .bcf
     - Binary VCF, 5-8x faster
   * - Variants
     - PGEN
     - .pgen
     - PLINK2 format, 25x faster
   * - Regions
     - BED
     - .bed
     - Tab-delimited coordinates
   * - Regions
     - GTF
     - .gtf
     - Gene annotations
   * - Regions
     - narrowPeak
     - .narrowPeak
     - MACS2 peak format

Output Formats
~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 25 55

   * - Module
     - Format
     - Description
   * - count
     - TSV
     - Allele counts per SNP
   * - count-sc
     - h5ad
     - AnnData with cell x SNP matrix
   * - map
     - BAM
     - WASP-filtered aligned reads
   * - analyze
     - TSV
     - Statistical test results

Format Selection Guide
----------------------

Variant Format Selection
~~~~~~~~~~~~~~~~~~~~~~~~

Choose based on your needs:

**VCF (default)**
   - Universal compatibility
   - Best for small datasets or testing
   - Install cyvcf2 for 7x speedup

**BCF**
   - Binary VCF, no information loss
   - 5-8x faster than VCF
   - Good balance of speed and compatibility

**PGEN**
   - Fastest format (25x speedup)
   - Lowest memory usage
   - Best for large cohorts (>10M variants)
   - Requires separate conversion step

Region Format Selection
~~~~~~~~~~~~~~~~~~~~~~~

**BED**
   - Simple coordinates (peaks, windows)
   - ATAC-seq peaks

**GTF/GFF3**
   - Gene annotations with hierarchy
   - RNA-seq analysis

**narrowPeak**
   - MACS2 peak format
   - Direct from peak caller
