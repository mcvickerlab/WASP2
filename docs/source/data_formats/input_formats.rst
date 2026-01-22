Input Formats
=============

Detailed specifications for WASP2 input files.

.. contents:: Format Types
   :local:
   :depth: 2

BAM Files
---------

Aligned sequencing reads in Binary Alignment Map format.

Requirements
~~~~~~~~~~~~

- **Sorted**: Coordinate-sorted (not name-sorted)
- **Indexed**: Accompanying .bai index file
- **Aligned**: Reads must be aligned to reference genome

Verification
~~~~~~~~~~~~

.. code-block:: bash

   # Check if sorted
   samtools view -H sample.bam | grep "@HD" | grep "SO:coordinate"

   # Check index exists
   ls -l sample.bam.bai

   # Create index if missing
   samtools index sample.bam

Single-Cell BAM
~~~~~~~~~~~~~~~

For 10x Genomics data, additional tags are used:

- ``CB``: Cell barcode
- ``UB``: UMI barcode

WASP2 automatically extracts cell barcodes from the CB tag.

Variant Formats
---------------

VCF (Variant Call Format)
~~~~~~~~~~~~~~~~~~~~~~~~~

Standard format for variant data.

**Requirements**:

- Bgzip compressed (.vcf.gz) with tabix index (.tbi)
- Contains genotype information (GT field)
- Sample IDs match those used in analysis

**Example**:

.. code-block:: text

   ##fileformat=VCFv4.2
   ##contig=<ID=chr1,length=248956422>
   ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
   #CHROM  POS     ID      REF  ALT  QUAL  FILTER  INFO    FORMAT  NA12878
   chr1    100000  rs123   A    G    .     PASS    .       GT      0/1
   chr1    200000  rs456   C    T    .     PASS    .       GT      1/1

**Index creation**:

.. code-block:: bash

   # Compress if not already
   bgzip variants.vcf

   # Create index
   tabix -p vcf variants.vcf.gz
   # or
   bcftools index -t variants.vcf.gz

**Performance enhancement**:

.. code-block:: bash

   # Install cyvcf2 for 7x faster parsing
   pip install wasp2[cyvcf2]

BCF (Binary VCF)
~~~~~~~~~~~~~~~~

Binary version of VCF format.

**Advantages**:

- 5-8x faster than VCF.gz
- No information loss
- Smaller file size

**Conversion**:

.. code-block:: bash

   # VCF to BCF
   bcftools view -O b variants.vcf.gz -o variants.bcf

   # Create index
   bcftools index variants.bcf

**Usage**:

.. code-block:: bash

   wasp2-count count-variants sample.bam variants.bcf --samples NA12878

PGEN (PLINK2 Binary)
~~~~~~~~~~~~~~~~~~~~

High-performance genotype format from PLINK2.

**Advantages**:

- 25x faster than VCF
- Lower memory usage
- Best for large cohorts

**Requirements**:

- Three files: .pgen, .pvar, .psam
- Install pgenlib: ``pip install wasp2[plink]``

**Conversion from VCF**:

.. code-block:: bash

   # Download plink2 if needed
   # https://www.cog-genomics.org/plink/2.0/

   # Convert VCF to PGEN
   plink2 --vcf variants.vcf.gz \
       --make-pgen \
       --out variants

   # Files created: variants.pgen, variants.pvar, variants.psam

**Usage**:

.. code-block:: bash

   wasp2-count count-variants sample.bam variants.pgen --samples NA12878

**Note**: PGEN stores only genotypes, not full VCF annotations.

Format Comparison
~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 20 20 20 20

   * - Format
     - Read Speed
     - Memory
     - File Size
     - Compatibility
   * - VCF.gz
     - 1x (baseline)
     - Medium
     - Large
     - Universal
   * - VCF.gz + cyvcf2
     - 7x
     - Medium
     - Large
     - Universal
   * - BCF
     - 5-8x
     - Medium
     - Medium
     - High
   * - PGEN
     - 25x
     - Low
     - Small
     - Medium

Region Formats
--------------

BED Format
~~~~~~~~~~

Simple tab-delimited format for genomic intervals.

**Columns** (minimum 3, up to 12):

1. Chromosome
2. Start position (0-based)
3. End position (1-based, exclusive)
4. Name (optional, used with ``--use_region_names``)
5. Score (optional)
6. Strand (optional)

**Example**:

.. code-block:: text

   chr1    100000  101000  peak_1  500 .
   chr1    200000  201500  peak_2  750 .
   chr2    50000   52000   peak_3  300 .

**Usage**:

.. code-block:: bash

   wasp2-count count-variants sample.bam variants.vcf.gz \
       --region peaks.bed \
       --use_region_names

GTF Format
~~~~~~~~~~

Gene annotation format used for RNA-seq analysis.

**Key fields**:

- Feature type (column 3): gene, transcript, exon, CDS
- Attributes (column 9): gene_id, transcript_id, gene_name

**Example**:

.. code-block:: text

   chr1  HAVANA  gene        11869   14409   .  +  .  gene_id "ENSG00000223972"; gene_name "DDX11L1";
   chr1  HAVANA  transcript  11869   14409   .  +  .  gene_id "ENSG00000223972"; transcript_id "ENST00000456328";
   chr1  HAVANA  exon        11869   12227   .  +  .  gene_id "ENSG00000223972"; transcript_id "ENST00000456328";

**Usage**:

.. code-block:: bash

   wasp2-count count-variants sample.bam variants.vcf.gz \
       --region genes.gtf \
       --gene_feature exon \
       --gene_attribute gene_id

**Recommended sources**:

- GENCODE: https://www.gencodegenes.org/
- Ensembl: https://www.ensembl.org/

GFF3 Format
~~~~~~~~~~~

General Feature Format version 3, similar to GTF.

**Differences from GTF**:

- Attributes use ``=`` instead of space
- Uses ``ID`` and ``Parent`` attributes

**Usage**:

.. code-block:: bash

   wasp2-count count-variants sample.bam variants.vcf.gz \
       --region genes.gff3 \
       --gene_feature exon \
       --gene_attribute ID

narrowPeak Format
~~~~~~~~~~~~~~~~~

MACS2 peak format, extension of BED.

**Columns**:

1-3: Standard BED (chr, start, end)
4: Name
5: Score
6: Strand
7: Signal value
8: pValue (-log10)
9: qValue (-log10)
10: Peak summit position

**Example**:

.. code-block:: text

   chr1    100000  101000  peak_1  500  .  15.5  8.2  6.1  450
   chr1    200000  201500  peak_2  750  .  22.3  12.4 10.8 600

**Direct usage** (no conversion needed):

.. code-block:: bash

   wasp2-count count-variants sample.bam variants.vcf.gz \
       --region peaks.narrowPeak

Cell Barcode Files
------------------

For single-cell analysis.

Format
~~~~~~

Plain text file with one cell barcode per line.

**Example**:

.. code-block:: text

   AAACCTGAGAAACCAT-1
   AAACCTGAGAAACCGC-1
   AAACCTGAGAAACCTA-1
   AAACCTGAGAAACGTC-1

**Source** (10x Genomics):

.. code-block:: bash

   # From filtered feature-barcode matrix
   zcat filtered_feature_bc_matrix/barcodes.tsv.gz > cell_barcodes.txt

Cell Type Mapping Files
-----------------------

For single-cell analysis by cell type.

Format
~~~~~~

Two-column TSV: barcode, cell_type

**Example**:

.. code-block:: text

   AAACCTGAGAAACCAT-1	CD4_T
   AAACCTGAGAAACCGC-1	CD4_T
   AAACCTGAGAAACCTA-1	CD8_T
   AAACCTGAGAAACGTC-1	B_cell
   AAACCTGAGAAAGTGG-1	Monocyte

**Creation from Seurat** (R):

.. code-block:: r

   write.table(
     data.frame(barcode = colnames(seurat_obj),
                celltype = seurat_obj$celltype),
     "barcode_celltype_map.tsv",
     sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
   )

**Creation from Scanpy** (Python):

.. code-block:: python

   import pandas as pd

   df = adata.obs[['celltype']].reset_index()
   df.columns = ['barcode', 'celltype']
   df.to_csv('barcode_celltype_map.tsv', sep='\t', index=False, header=False)
