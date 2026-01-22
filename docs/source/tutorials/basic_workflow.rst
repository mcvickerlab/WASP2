Basic Workflow Tutorial
=======================

**Time**: 30 minutes
**Level**: Beginner to Intermediate
**Prerequisites**: WASP2 installed, basic command line skills

This tutorial walks through a complete WASP2 workflow from aligned reads to
statistical analysis of allelic imbalance.

.. contents:: Steps
   :local:
   :depth: 2

Overview
--------

The complete WASP2 pipeline consists of:

1. **Quality Control**: Verify input files
2. **WASP Mapping** (optional): Correct reference bias
3. **Allele Counting**: Quantify ref/alt reads per SNP
4. **Statistical Analysis**: Detect significant allelic imbalance
5. **Interpretation**: Understand and visualize results

.. code-block:: text

   Raw Reads (FASTQ)
       |
       v
   Standard Alignment (STAR/BWA/bowtie2)
       |
       v
   +-------------------+
   |   WASP Mapping    |  <- Optional but recommended
   |   (wasp2-map)     |
   +-------------------+
       |
       v
   Unbiased BAM
       |
       v
   +-------------------+
   |  Allele Counting  |
   |  (wasp2-count)    |
   +-------------------+
       |
       v
   Allele Counts (TSV)
       |
       v
   +-------------------+
   | Statistical Test  |
   | (wasp2-analyze)   |
   +-------------------+
       |
       v
   Imbalance Results

Data Requirements
-----------------

Before starting, ensure you have:

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - File
     - Description
     - Requirements
   * - BAM
     - Aligned sequencing reads
     - Sorted and indexed (.bai file)
   * - VCF/BCF/PGEN
     - Variant genotypes
     - Contains sample of interest, indexed
   * - Region file (optional)
     - Gene/peak annotations
     - GTF, GFF3, BED, or narrowPeak format

Step 1: Quality Control
-----------------------

Verify your input files before processing.

Check BAM file
~~~~~~~~~~~~~~

.. code-block:: bash

   # Verify BAM is sorted and indexed
   samtools view -H sample.bam | head -20

   # Check alignment statistics
   samtools flagstat sample.bam

   # Verify index exists
   ls -l sample.bam.bai

**Expected output**: Sorted header lines, reasonable alignment rates (>70%),
and existing .bai index file.

Check VCF file
~~~~~~~~~~~~~~

.. code-block:: bash

   # List samples in VCF
   bcftools query -l variants.vcf.gz

   # Count total variants
   bcftools view -H variants.vcf.gz | wc -l

   # Count heterozygous SNPs for your sample
   bcftools view -s NA12878 -g het variants.vcf.gz | bcftools view -H | wc -l

   # Verify index exists
   ls -l variants.vcf.gz.tbi

**Expected output**: Your sample ID in the list, reasonable variant counts,
heterozygous SNPs present, and existing index file.

Verify coordinate compatibility
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Check chromosome naming in BAM
   samtools view -H sample.bam | grep "@SQ" | head -3

   # Check chromosome naming in VCF
   bcftools view -h variants.vcf.gz | grep "##contig" | head -3

**Important**: Chromosome names must match (e.g., both use "chr1" or both use "1").

Step 2: WASP Mapping (Optional)
-------------------------------

If your BAM hasn't been WASP-filtered, run this step to remove reference bias.

.. note::

   Skip this step if you:

   - Already have WASP-filtered reads
   - Used an allele-aware aligner
   - Are doing a quick initial analysis (can add later)

Step 2a: Generate swapped reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   wasp2-map make-reads sample.bam variants.vcf.gz \
       --samples NA12878 \
       --out_dir wasp_output/ \
       --threads 4

**Output files**:

- ``wasp_output/sample_swapped_alleles_r1.fq`` - Reads with swapped alleles
- ``wasp_output/sample_swapped_alleles_r2.fq`` - (paired-end only)
- ``wasp_output/sample_to_remap.bam`` - Original reads that need remapping
- ``wasp_output/sample_keep.bam`` - Reads that don't overlap SNPs
- ``sample_wasp_data_files.json`` - Metadata file

Step 2b: Remap swapped reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use your original aligner with the same parameters:

.. code-block:: bash

   # Example with BWA
   bwa mem -M reference.fa \
       wasp_output/sample_swapped_alleles_r1.fq \
       wasp_output/sample_swapped_alleles_r2.fq \
   | samtools view -b -F 4 - \
   | samtools sort -o wasp_output/sample_remapped.bam -

   samtools index wasp_output/sample_remapped.bam

Step 2c: Filter inconsistent mappings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   wasp2-map filter-remapped wasp_output/sample_remapped.bam \
       --json sample_wasp_data_files.json \
       --out_bam sample_wasp_filtered.bam \
       --threads 4

**Output**: ``sample_wasp_filtered.bam`` - Unbiased reads ready for counting

Step 3: Allele Counting
-----------------------

Count allele-specific reads at heterozygous SNPs.

Basic counting
~~~~~~~~~~~~~~

.. code-block:: bash

   wasp2-count count-variants sample_wasp_filtered.bam variants.vcf.gz \
       --samples NA12878 \
       --out_file counts.tsv

With region annotation (recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For RNA-seq with gene annotations:

.. code-block:: bash

   wasp2-count count-variants sample_wasp_filtered.bam variants.vcf.gz \
       --samples NA12878 \
       --region genes.gtf \
       --gene_feature exon \
       --gene_attribute gene_id \
       --out_file gene_counts.tsv

For ATAC-seq with peak annotations:

.. code-block:: bash

   wasp2-count count-variants sample_wasp_filtered.bam variants.vcf.gz \
       --samples NA12878 \
       --region peaks.narrowPeak \
       --out_file peak_counts.tsv

**Output format** (counts.tsv):

.. code-block:: text

   chr     pos       ref  alt  ref_count  alt_count  other_count  region
   chr1    1000000   A    G    15         12         0            ENSG00000123456
   chr1    1001000   C    T    8          10         1            ENSG00000123456
   chr1    2000000   G    A    20         18         0            ENSG00000789012

Step 4: Statistical Analysis
----------------------------

Test for significant allelic imbalance.

Basic analysis
~~~~~~~~~~~~~~

.. code-block:: bash

   wasp2-analyze find-imbalance counts.tsv \
       --out_file results.tsv

With filtering options
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   wasp2-analyze find-imbalance counts.tsv \
       --min 10 \
       --out_file results.tsv

**Parameters**:

- ``--min 10``: Require at least 10 total reads per site
- ``--pseudocount 1``: Add pseudocount to avoid division by zero

Gene-level aggregation (RNA-seq)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   wasp2-analyze find-imbalance gene_counts.tsv \
       --min 20 \
       --groupby gene_id \
       --out_file gene_results.tsv

**Output format** (results.tsv):

.. code-block:: text

   region           ref_count  alt_count  total  ratio   pvalue      fdr
   ENSG00000123456  45         38         83     0.542   0.4521      0.8234
   ENSG00000789012  120        60         180    0.667   0.00012     0.0034
   ENSG00000345678  25         75         100    0.250   0.00001     0.0005

Step 5: Interpretation
----------------------

Identify significant results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Genes with significant allelic imbalance (FDR < 0.05)
   awk 'NR==1 || $7 < 0.05' results.tsv

   # Count significant genes
   awk 'NR>1 && $7 < 0.05' results.tsv | wc -l

   # Sort by effect size (most imbalanced)
   sort -t$'\t' -k5,5n results.tsv | head -20

Understanding the output
~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 50 30

   * - Column
     - Description
     - Interpretation
   * - ref_count
     - Reads with reference allele
     -
   * - alt_count
     - Reads with alternate allele
     -
   * - total
     - Total reads (ref + alt)
     - Higher = more confident
   * - ratio
     - ref_count / total
     - 0.5 = balanced
   * - pvalue
     - Beta-binomial test p-value
     - Lower = more significant
   * - fdr
     - Benjamini-Hochberg adjusted p-value
     - < 0.05 typically significant

Quality checks
~~~~~~~~~~~~~~

1. **Coverage distribution**: Ensure reasonable read depths

   .. code-block:: bash

      # Check distribution of total counts
      cut -f4 results.tsv | tail -n +2 | sort -n | uniq -c | head -20

2. **Ratio distribution**: Should be centered near 0.5

   .. code-block:: bash

      # Check ratio distribution
      awk 'NR>1 {print int($5*10)/10}' results.tsv | sort | uniq -c

3. **Known controls**: Validate against expected imbalanced genes (e.g., imprinted genes)

Troubleshooting
---------------

No output SNPs
~~~~~~~~~~~~~~

**Symptoms**: Empty output file or only header

**Diagnostic**:

.. code-block:: bash

   # Check sample is in VCF
   bcftools query -l variants.vcf.gz | grep -w "NA12878"

   # Check for heterozygous SNPs
   bcftools view -s NA12878 -g het variants.vcf.gz | head

**Solutions**:

- Verify sample name spelling
- Check chromosome naming matches between BAM and VCF
- Ensure VCF contains heterozygous genotypes

Low counts everywhere
~~~~~~~~~~~~~~~~~~~~~

**Symptoms**: Most sites have < 10 reads

**Solutions**:

- Check sequencing depth: ``samtools depth sample.bam | awk '{sum+=$3; n++} END {print sum/n}'``
- Lower ``--min`` threshold if appropriate
- Focus on highly expressed genes/accessible peaks

Slow performance
~~~~~~~~~~~~~~~~

**Symptoms**: Processing takes hours

**Solutions**:

1. Install cyvcf2: ``pip install wasp2[cyvcf2]`` (7x faster)
2. Convert to PGEN format (25x faster):

   .. code-block:: bash

      plink2 --vcf variants.vcf.gz --make-pgen --out variants
      wasp2-count count-variants sample.bam variants.pgen ...

Next Steps
----------

Continue with data-type-specific tutorials:

- :doc:`rnaseq_ase` - RNA-seq allele-specific expression analysis
- :doc:`atacseq_ase` - ATAC-seq allelic chromatin accessibility
- :doc:`single_cell` - Single-cell allele-specific analysis

Or explore advanced topics:

- :doc:`/user_guide/counting` - Detailed counting options
- :doc:`/user_guide/analysis` - Statistical model details
- :doc:`/cli/wasp2_count` - Complete CLI reference
