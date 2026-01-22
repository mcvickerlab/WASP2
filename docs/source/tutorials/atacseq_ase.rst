ATAC-seq Allelic Chromatin Accessibility
========================================

**Time**: 45 minutes
**Level**: Intermediate
**Prerequisites**: Completed :doc:`basic_workflow`, ATAC-seq data

This tutorial covers allele-specific chromatin accessibility analysis with
ATAC-seq data, helping you identify regulatory variants affecting chromatin state.

.. contents:: Topics
   :local:
   :depth: 2

Overview
--------

ATAC-seq measures chromatin accessibility - regions where DNA is open and
accessible to transcription factors. Allele-specific ATAC-seq analysis
identifies genomic regions where accessibility differs between alleles.

Biological Applications
~~~~~~~~~~~~~~~~~~~~~~~

1. **Identify caQTLs**: Chromatin accessibility quantitative trait loci
2. **Map regulatory variants**: SNPs that affect TF binding
3. **Understand disease variants**: How GWAS variants affect chromatin
4. **Study allelic TF binding**: Differential transcription factor access

Key Differences from RNA-seq
----------------------------

.. list-table::
   :header-rows: 1
   :widths: 25 35 40

   * - Aspect
     - RNA-seq
     - ATAC-seq
   * - **Features**
     - Genes/Transcripts
     - Peaks/Open regions
   * - **Annotation**
     - GTF/GFF3
     - BED/narrowPeak
   * - **Coverage**
     - Exons
     - Open chromatin
   * - **Expected AI**
     - Imprinting, eQTLs
     - caQTLs, TF binding
   * - **Read length**
     - Typically longer
     - Often shorter (paired-end recommended)
   * - **Fragment sizes**
     - Variable
     - Characteristic nucleosome pattern

ATAC-seq Preprocessing
----------------------

Peak calling
~~~~~~~~~~~~

Before allele counting, you need to define accessible regions (peaks):

.. code-block:: bash

   # Call peaks with MACS2
   macs2 callpeak \
       -t sample.bam \
       -f BAMPE \
       -g hs \
       -n sample \
       --outdir peaks/ \
       -q 0.01 \
       --nomodel \
       --shift -100 \
       --extsize 200

**Output**: ``peaks/sample_peaks.narrowPeak``

Quality check peaks
~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Count total peaks
   wc -l peaks/sample_peaks.narrowPeak

   # Check peak size distribution
   awk '{print $3-$2}' peaks/sample_peaks.narrowPeak | sort -n | \
       awk 'BEGIN{min=999999;max=0} {sum+=$1; if($1<min)min=$1; if($1>max)max=$1} \
            END{print "Min:", min, "Max:", max, "Mean:", sum/NR}'

   # View top peaks
   sort -k5,5rn peaks/sample_peaks.narrowPeak | head

**Expected**: 50,000-150,000 peaks for a good ATAC-seq library.

Step-by-Step Workflow
---------------------

Step 1: WASP mapping
~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Generate swapped allele reads
   wasp2-map make-reads sample.bam variants.vcf.gz \
       --samples NA12878 \
       --out_dir wasp_output/ \
       --threads 4

   # Remap with BWA or bowtie2
   bwa mem -M reference.fa \
       wasp_output/sample_swapped_alleles_r1.fq \
       wasp_output/sample_swapped_alleles_r2.fq \
   | samtools view -b -F 4 - \
   | samtools sort -o wasp_output/sample_remapped.bam -

   samtools index wasp_output/sample_remapped.bam

   # Filter remapped reads
   wasp2-map filter-remapped wasp_output/sample_remapped.bam \
       --json sample_wasp_data_files.json \
       --out_bam sample_wasp_filtered.bam

Step 2: Peak-level allele counting
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   wasp2-count count-variants sample_wasp_filtered.bam variants.vcf.gz \
       --samples NA12878 \
       --region peaks/sample_peaks.narrowPeak \
       --use_region_names \
       --out_file peak_counts.tsv

**Key parameters**:

- ``--region peaks/sample_peaks.narrowPeak``: Use ATAC-seq peaks
- ``--use_region_names``: Use peak IDs (4th column) in output

Step 3: Statistical analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   wasp2-analyze find-imbalance peak_counts.tsv \
       --min 10 \
       --out_file peak_ase_results.tsv

Step 4: Identify significant peaks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Peaks with significant allelic imbalance
   awk 'NR==1 || $7 < 0.05' peak_ase_results.tsv > significant_peaks.tsv

   # Count significant peaks
   awk 'NR>1' significant_peaks.tsv | wc -l

   # Extract coordinates for downstream analysis
   awk 'NR>1 && $7 < 0.05 {print $1}' peak_ase_results.tsv > significant_peak_ids.txt

Expected Results
----------------

For a typical ATAC-seq dataset:

- **Peaks with informative SNPs**: 5,000-20,000
- **Peaks with sufficient coverage**: 2,000-10,000
- **Peaks with significant AI**: 100-1,000 (FDR < 0.05)

Downstream Analysis
-------------------

Motif enrichment analysis
~~~~~~~~~~~~~~~~~~~~~~~~~

Identify transcription factor motifs disrupted by variants:

.. code-block:: bash

   # Extract sequences around imbalanced SNPs
   awk 'NR>1 && $7 < 0.05 {print $1"\t"$2-1"\t"$2}' peak_ase_results.tsv \
       > imbalanced_snps.bed

   # Run HOMER motif analysis
   findMotifsGenome.pl imbalanced_snps.bed hg38 motif_results/ -size 200

   # Or use MEME Suite
   bedtools getfasta -fi genome.fa -bed imbalanced_snps.bed > imbalanced_seqs.fa
   meme imbalanced_seqs.fa -dna -oc meme_results/ -nostatus -mod zoops

Annotate with genomic features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Annotate peaks with HOMER
   annotatePeaks.pl significant_peaks.bed hg38 > annotated_peaks.txt

   # Or use bedtools with GENCODE
   bedtools intersect -a significant_peaks.bed -b gencode.v38.annotation.gtf -wa -wb \
       > peaks_gene_overlap.bed

Link to gene expression
~~~~~~~~~~~~~~~~~~~~~~~

Compare allelic chromatin accessibility with allelic expression:

.. code-block:: bash

   # Find peaks near genes with ASE
   bedtools window -a significant_peaks.bed -b significant_ase_genes.bed -w 100000 \
       > peaks_near_ase_genes.bed

Integration with GWAS
~~~~~~~~~~~~~~~~~~~~~

Check if imbalanced peaks overlap disease-associated variants:

.. code-block:: bash

   # Download GWAS catalog
   wget https://www.ebi.ac.uk/gwas/api/search/downloads/full

   # Find overlaps
   bedtools intersect -a significant_peaks.bed -b gwas_catalog.bed -wa -wb \
       > peaks_gwas_overlap.bed

Troubleshooting
---------------

Very few peaks with SNPs
~~~~~~~~~~~~~~~~~~~~~~~~

**Problem**: Most peaks have no informative heterozygous SNPs

**Solutions**:

- Use imputed genotypes (higher SNP density)
- Verify VCF and BAM use same reference genome
- Check peak calling worked (should have 50K+ peaks)

Low coverage in peaks
~~~~~~~~~~~~~~~~~~~~~

**Problem**: Most peaks have < 10 reads

**Possible causes**:

- Low sequencing depth (need 30-50M reads)
- Poor library quality (check insert size distribution)
- Over-transposition (too much Tn5)

**Diagnostic**:

.. code-block:: bash

   # Check coverage in peaks
   bedtools coverage -a peaks/sample_peaks.narrowPeak -b sample.bam \
       | awk '{print $NF}' | sort -n | uniq -c | head -20

Many false positives
~~~~~~~~~~~~~~~~~~~~

**Problem**: Results include many suspicious peaks

**Solutions**:

- Apply stricter FDR threshold (0.01 instead of 0.05)
- Require minimum effect size (ratio < 0.35 or > 0.65)
- Filter peaks in problematic regions (repeat elements, centromeres)
- Use WASP filtering to remove reference bias

.. code-block:: bash

   # Stricter filtering
   awk 'NR==1 || ($7 < 0.01 && ($5 < 0.35 || $5 > 0.65))' peak_ase_results.tsv \
       > high_confidence_peaks.tsv

Visualization
-------------

Python example
~~~~~~~~~~~~~~

.. code-block:: python

   import pandas as pd
   import matplotlib.pyplot as plt
   import numpy as np

   # Load results
   df = pd.read_csv('peak_ase_results.tsv', sep='\t')

   # Histogram of allelic ratios
   fig, axes = plt.subplots(1, 2, figsize=(12, 5))

   # All peaks
   axes[0].hist(df['ratio'], bins=50, edgecolor='black')
   axes[0].axvline(0.5, color='red', linestyle='--')
   axes[0].set_xlabel('Allelic Ratio')
   axes[0].set_ylabel('Count')
   axes[0].set_title('All Peaks')

   # Significant peaks only
   sig = df[df['fdr'] < 0.05]
   axes[1].hist(sig['ratio'], bins=30, edgecolor='black', color='orange')
   axes[1].axvline(0.5, color='red', linestyle='--')
   axes[1].set_xlabel('Allelic Ratio')
   axes[1].set_title(f'Significant Peaks (n={len(sig)})')

   plt.tight_layout()
   plt.savefig('atac_ase_histogram.png', dpi=150)

Genome browser visualization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create a bigWig track showing allelic imbalance:

.. code-block:: bash

   # Create bedGraph with imbalance scores
   awk 'NR>1 {print $1"\t"$2-1"\t"$2"\t"$5-0.5}' peak_ase_results.tsv \
       > imbalance.bedGraph

   # Convert to bigWig
   bedGraphToBigWig imbalance.bedGraph chrom.sizes imbalance.bw

Next Steps
----------

- :doc:`single_cell` - Single-cell ATAC-seq analysis
- :doc:`/user_guide/counting` - Advanced counting options
- :doc:`/faq` - Common questions
