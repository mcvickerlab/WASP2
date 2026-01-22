RNA-seq Allele-Specific Expression
===================================

**Time**: 45 minutes
**Level**: Intermediate
**Prerequisites**: Completed :doc:`basic_workflow`, RNA-seq data

This tutorial covers allele-specific expression (ASE) analysis with RNA-seq data,
including gene-level aggregation and interpretation of results.

.. contents:: Topics
   :local:
   :depth: 2

Overview
--------

RNA-seq ASE analysis identifies genes where one allele is preferentially expressed.
This can reveal:

- **Cis-regulatory variants**: SNPs affecting transcription factor binding
- **Imprinted genes**: Parent-of-origin expression patterns
- **X-inactivation effects**: Monoallelic X-linked gene expression
- **eQTL effects**: Expression quantitative trait loci

Biological Questions
~~~~~~~~~~~~~~~~~~~~

1. Which genes show preferential expression of one allele?
2. Are there parent-of-origin effects (imprinting)?
3. Do allelic ratios differ between tissues or conditions?

RNA-seq Specific Considerations
-------------------------------

Key differences from general workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Aspect
     - RNA-seq Specifics
   * - Annotation
     - Use GTF/GFF3 with gene features (exon, CDS)
   * - Feature level
     - Typically analyze at gene level, not SNP level
   * - Coverage
     - Varies by expression; housekeeping genes have more reads
   * - Splicing
     - Consider transcript-level analysis for splice variants

Alignment recommendations
~~~~~~~~~~~~~~~~~~~~~~~~~

For RNA-seq data, use splice-aware aligners:

- **STAR** (recommended): Fast, accurate, widely used
- **HISAT2**: Good alternative, lower memory
- **TopHat2**: Legacy option

.. code-block:: bash

   # STAR alignment example
   STAR --runThreadN 8 \
       --genomeDir /path/to/star_index \
       --readFilesIn sample_R1.fastq.gz sample_R2.fastq.gz \
       --readFilesCommand zcat \
       --outSAMtype BAM SortedByCoordinate \
       --outFileNamePrefix sample_

Step-by-Step Workflow
---------------------

Step 1: Prepare GTF annotation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

WASP2 works best with standard GTF files (e.g., GENCODE, Ensembl):

.. code-block:: bash

   # Download GENCODE annotations (example)
   wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
   gunzip gencode.v38.annotation.gtf.gz

   # Verify format
   head -50 gencode.v38.annotation.gtf

Step 2: WASP mapping (recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For unbiased ASE analysis:

.. code-block:: bash

   # Generate swapped allele reads
   wasp2-map make-reads sample_Aligned.sortedByCoord.out.bam variants.vcf.gz \
       --samples NA12878 \
       --out_dir wasp_output/ \
       --threads 8

   # Remap with STAR
   STAR --runThreadN 8 \
       --genomeDir /path/to/star_index \
       --readFilesIn wasp_output/sample_swapped_alleles_r1.fq wasp_output/sample_swapped_alleles_r2.fq \
       --outSAMtype BAM SortedByCoordinate \
       --outFileNamePrefix wasp_output/sample_remapped_

   # Filter remapped reads
   wasp2-map filter-remapped wasp_output/sample_remapped_Aligned.sortedByCoord.out.bam \
       --json sample_wasp_data_files.json \
       --out_bam sample_wasp_filtered.bam

Step 3: Gene-level allele counting
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   wasp2-count count-variants sample_wasp_filtered.bam variants.vcf.gz \
       --samples NA12878 \
       --region gencode.v38.annotation.gtf \
       --gene_feature exon \
       --gene_attribute gene_id \
       --gene_parent transcript_id \
       --out_file gene_snp_counts.tsv

**Key parameters for RNA-seq**:

- ``--gene_feature exon``: Count SNPs in exons (coding regions)
- ``--gene_attribute gene_id``: Use Ensembl gene IDs
- ``--gene_parent transcript_id``: Track transcript association

Step 4: Gene-level statistical analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Aggregate SNPs to gene level:

.. code-block:: bash

   wasp2-analyze find-imbalance gene_snp_counts.tsv \
       --min 20 \
       --groupby gene_id \
       --out_file gene_ase_results.tsv

**Parameters explained**:

- ``--min 20``: Require at least 20 reads per gene (more stringent for genes)
- ``--groupby gene_id``: Sum counts across all SNPs in each gene

Step 5: Interpret results
~~~~~~~~~~~~~~~~~~~~~~~~~

Extract significant genes:

.. code-block:: bash

   # Genes with significant ASE (FDR < 0.05)
   awk 'NR==1 || $7 < 0.05' gene_ase_results.tsv > significant_ase_genes.tsv

   # Count significant genes
   wc -l significant_ase_genes.tsv

   # Top 20 most imbalanced genes
   sort -t$'\t' -k5,5n gene_ase_results.tsv | head -21

Expected Results
----------------

For a typical RNA-seq dataset (30-50M reads):

- **Genes with sufficient coverage**: 10,000-15,000
- **Genes with significant ASE**: 500-2,000 (FDR < 0.05)
- **Strong ASE (ratio < 0.3 or > 0.7)**: 100-500

Validation with known imprinted genes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Check for expected patterns in known imprinted genes:

.. code-block:: bash

   # Known imprinted genes (examples)
   grep -E "SNRPN|H19|IGF2|MEG3|PEG3|MEST" gene_ase_results.tsv

These genes should show strong allelic imbalance if your sample
expresses them and has informative heterozygous SNPs.

Common Gene Categories
~~~~~~~~~~~~~~~~~~~~~~

Genes with ASE often fall into these categories:

1. **Imprinted genes**: H19, IGF2, SNRPN, KCNQ1OT1, PEG3, MEG3
2. **X-linked genes** (in females): Random X-inactivation effects
3. **Genes with eQTLs**: Known expression QTLs from GTEx/eQTLGen
4. **Highly polymorphic genes**: HLA genes, olfactory receptors

Troubleshooting
---------------

Low gene coverage
~~~~~~~~~~~~~~~~~

**Problem**: Many genes have < 20 reads

**Solutions**:

- Lower ``--min`` threshold to 10 (but interpret cautiously)
- Ensure BAM file has good alignment rate (> 70%)
- Check if VCF sample has heterozygous SNPs in expressed genes
- Consider deeper sequencing for future experiments

No heterozygous SNPs in genes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Problem**: Most genes have 0 informative SNPs

**Diagnostic**:

.. code-block:: bash

   # Count SNPs per gene
   cut -f8 gene_snp_counts.tsv | sort | uniq -c | sort -rn | head

**Solutions**:

- Use a more complete VCF (e.g., imputed genotypes)
- Verify VCF sample matches RNA-seq sample
- Check species and genome build match

Unexpected patterns
~~~~~~~~~~~~~~~~~~~

**Problem**: Results don't match expectations (e.g., known imprinted genes show 50:50)

**Possible causes**:

- Tissue-specific imprinting (not all tissues show imprinting)
- Loss of imprinting in cell lines
- Sample is homozygous at informative SNPs
- Technical issues with genotyping

Advanced Analysis
-----------------

Transcript-level analysis
~~~~~~~~~~~~~~~~~~~~~~~~~

For splice-variant specific ASE:

.. code-block:: bash

   wasp2-count count-variants sample_wasp_filtered.bam variants.vcf.gz \
       --samples NA12878 \
       --region gencode.v38.annotation.gtf \
       --gene_feature exon \
       --gene_attribute transcript_id \
       --out_file transcript_snp_counts.tsv

   wasp2-analyze find-imbalance transcript_snp_counts.tsv \
       --min 10 \
       --groupby transcript_id \
       --out_file transcript_ase_results.tsv

Comparing conditions
~~~~~~~~~~~~~~~~~~~~

To compare ASE between conditions (e.g., treatment vs control):

1. Run ASE analysis separately for each condition
2. Compare allelic ratios for the same genes
3. Look for genes where ratio changes significantly

.. code-block:: bash

   # Example: join results from two conditions
   join -t$'\t' -j1 \
       <(sort -k1 control_ase_results.tsv) \
       <(sort -k1 treatment_ase_results.tsv) \
   > combined_ase.tsv

Visualization (Python)
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import pandas as pd
   import matplotlib.pyplot as plt

   # Load results
   df = pd.read_csv('gene_ase_results.tsv', sep='\t')

   # Volcano plot
   plt.figure(figsize=(10, 8))
   colors = ['red' if fdr < 0.05 else 'gray' for fdr in df['fdr']]
   plt.scatter(df['ratio'] - 0.5, -np.log10(df['pvalue']), c=colors, alpha=0.5)
   plt.xlabel('Allelic Imbalance (ratio - 0.5)')
   plt.ylabel('-log10(p-value)')
   plt.axhline(-np.log10(0.05), linestyle='--', color='gray')
   plt.title('RNA-seq Allele-Specific Expression')
   plt.savefig('ase_volcano.png', dpi=150)

Next Steps
----------

- :doc:`single_cell` - Single-cell RNA-seq ASE analysis
- :doc:`/user_guide/analysis` - Statistical model details
- :doc:`/faq` - Common questions about ASE analysis
