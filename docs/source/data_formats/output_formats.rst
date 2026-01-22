Output Formats
==============

WASP2 produces various output files depending on the module used. This page documents the structure and content of each output format.

Counting Output
---------------

Allele Counts TSV
~~~~~~~~~~~~~~~~~

The main output from ``wasp2-count count-variants`` is a tab-separated file containing allele counts per variant.

**Columns:**

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Column
     - Type
     - Description
   * - chrom
     - str
     - Chromosome name
   * - pos
     - int
     - 1-based genomic position
   * - ref
     - str
     - Reference allele
   * - alt
     - str
     - Alternate allele
   * - ref_count
     - int
     - Number of reads supporting reference allele
   * - alt_count
     - int
     - Number of reads supporting alternate allele
   * - total_count
     - int
     - Total read coverage (ref_count + alt_count)
   * - sample
     - str
     - Sample identifier (if multiple samples)

**Example:**

.. code-block:: text

   chrom	pos	ref	alt	ref_count	alt_count	total_count	sample
   chr1	12345	A	G	15	12	27	NA12878
   chr1	23456	C	T	8	10	18	NA12878

Single-Cell Counts
~~~~~~~~~~~~~~~~~~

For single-cell data (``wasp2-count count-variants-sc``), output includes cell barcode information.

**Additional columns:**

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Column
     - Type
     - Description
   * - barcode
     - str
     - Cell barcode sequence
   * - umi_count
     - int
     - Number of unique UMIs (if UMI-aware counting enabled)

Analysis Output
---------------

Allelic Imbalance Results
~~~~~~~~~~~~~~~~~~~~~~~~~

Output from ``wasp2-analyze find-imbalance`` contains statistical test results.

**Columns:**

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Column
     - Type
     - Description
   * - chrom
     - str
     - Chromosome name
   * - pos
     - int
     - 1-based genomic position
   * - ref
     - str
     - Reference allele
   * - alt
     - str
     - Alternate allele
   * - ref_count
     - int
     - Reference allele count
   * - alt_count
     - int
     - Alternate allele count
   * - total_count
     - int
     - Total coverage
   * - ref_ratio
     - float
     - Reference allele ratio (ref_count / total_count)
   * - pvalue
     - float
     - P-value from beta-binomial test
   * - fdr
     - float
     - Benjamini-Hochberg adjusted p-value
   * - significant
     - bool
     - Whether FDR < threshold (default 0.05)

**Example:**

.. code-block:: text

   chrom	pos	ref	alt	ref_count	alt_count	total_count	ref_ratio	pvalue	fdr	significant
   chr1	12345	A	G	45	12	57	0.789	0.0001	0.0012	True
   chr1	23456	C	T	25	27	52	0.481	0.8234	0.9123	False

Group Comparison Results
~~~~~~~~~~~~~~~~~~~~~~~~

Output from ``wasp2-analyze compare-groups`` for differential allelic imbalance.

**Additional columns:**

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Column
     - Type
     - Description
   * - group1_ratio
     - float
     - Mean reference ratio in group 1
   * - group2_ratio
     - float
     - Mean reference ratio in group 2
   * - delta_ratio
     - float
     - Difference between group ratios
   * - comparison_pvalue
     - float
     - P-value for group comparison
   * - comparison_fdr
     - float
     - FDR-adjusted comparison p-value

Mapping Output
--------------

WASP-Filtered BAM
~~~~~~~~~~~~~~~~~

The mapping module outputs BAM files with reads that pass WASP filtering.

**BAM Tags:**

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - Tag
     - Description
   * - vW
     - WASP filter status: 1=pass, 2=fail (maps to different location), 3=fail (maps to multiple locations)
   * - vA
     - Number of variants overlapping the read
   * - vG
     - Variant genotypes for overlapping variants

Remap Statistics
~~~~~~~~~~~~~~~~

Summary statistics from the remapping process.

**Contents:**

- Total reads processed
- Reads passing WASP filter
- Reads failing due to mapping differences
- Reads failing due to multi-mapping
- Variant overlap statistics

File Compression
----------------

All TSV outputs can be gzip-compressed by specifying ``.gz`` extension:

.. code-block:: bash

   wasp2-count count-variants sample.bam variants.vcf.gz \
       --out_file counts.tsv.gz

Index Files
-----------

For large outputs, WASP2 can generate index files for efficient random access:

- ``.tbi`` - Tabix index for TSV files
- ``.csi`` - CSI index for large chromosome coordinates

Generate indexes with:

.. code-block:: bash

   wasp2-count count-variants sample.bam variants.vcf.gz \
       --out_file counts.tsv.gz \
       --create_index
