Comparative Imbalance Analysis Tutorial
=======================================

This tutorial provides a comprehensive guide to detecting **differential allelic imbalance**
between cell types, conditions, or biological groups using WASP2's comparative analysis module.

Overview
--------

**What is Comparative Imbalance Analysis?**

While standard allelic imbalance (AI) analysis detects whether a genomic region shows
preferential expression of one allele, comparative imbalance analysis asks a different
question: **Does the degree of imbalance differ between groups?**

This is powerful for identifying:

* **Cell-type-specific regulatory variation** - Regions where different cell types show
  distinct allelic preferences
* **Condition-dependent effects** - Treatment-induced changes in allelic regulation
* **Sex differences** - Chromatin regions with sex-biased allelic accessibility
* **Developmental dynamics** - Stage-specific changes in allelic regulation

**Statistical Approach**

WASP2 uses a **likelihood ratio test (LRT)** to compare two hypotheses:

.. code-block:: text

   Null Hypothesis (H0):     Both groups share the same allelic imbalance (μ_combined)
   Alternative Hypothesis (H1): Groups have different imbalance (μ₁ ≠ μ₂)

   Test Statistic:  LRT = -2 × (log L_null - log L_alt)
   P-value:         P(χ²(df=1) > LRT)

The test accounts for overdispersion using beta-binomial modeling and applies
Benjamini-Hochberg FDR correction for multiple testing.

Prerequisites
-------------

**Software:**

* WASP2 (``pip install wasp2``)
* Python with scanpy/pandas for visualization (optional)
* R with Seurat for cell type annotation (optional)

**Data Requirements:**

* AnnData count matrix (``.h5ad``) with allele counts per cell per SNP
* Barcode-to-group mapping file (TSV)
* Groups can be: cell types, conditions, sex, treatment status, etc.

Input Data Format
-----------------

Count Matrix (.h5ad)
~~~~~~~~~~~~~~~~~~~~

Your AnnData object should have this structure:

.. code-block:: text

   AnnData object (n_snps × n_cells)
   ├── .obs                 # SNP metadata (rows)
   │   ├── index            # SNP identifiers
   │   └── [sample_name]    # Genotypes: '0|1', '1|0', '0/1', '1/0'
   │
   ├── .var                 # Cell metadata (columns)
   │   └── group            # Cell type/group assignment
   │
   ├── .layers
   │   ├── "ref"            # Reference allele counts (sparse matrix)
   │   └── "alt"            # Alternate allele counts (sparse matrix)
   │
   └── .uns
       ├── feature          # DataFrame: SNP → region mapping
       └── samples          # List of sample names

**Create counts from BAM + VCF:**

.. code-block:: bash

   wasp2-count count-variants-sc \
     aligned.bam \
     phased_variants.vcf.gz \
     barcodes.txt \
     --samples SAMPLE_ID \
     --feature peaks.bed \
     --out_file allele_counts.h5ad

Barcode Map (TSV)
~~~~~~~~~~~~~~~~~

A two-column tab-separated file (no header) mapping cell barcodes to groups:

.. code-block:: text

   AAACGAACAGTCAGTT-1    excitatory_neurons
   AAACGAAGTCGCTCTA-1    inhibitory_neurons
   AAACGAAGTGAACCTA-1    excitatory_neurons
   AAAGGATCATCGATGT-1    astrocytes
   AAAGGATGTGCAACGA-1    microglia

**Important:** Barcodes must exactly match those in the count matrix (including any ``-1`` suffix).

Tutorial 1: Cell Type Comparison
--------------------------------

This tutorial demonstrates comparing allelic imbalance between neuronal subtypes.

Step 1: Prepare Input Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Export cell type annotations from Seurat:**

.. code-block:: r

   library(Seurat)

   # Load your analyzed Seurat object
   seurat_obj <- readRDS("brain_snATAC.rds")

   # Create barcode-to-celltype mapping
   barcode_df <- data.frame(
     barcode = colnames(seurat_obj),
     celltype = Idents(seurat_obj)
   )

   # Write without header
   write.table(
     barcode_df,
     "barcode_celltype_map.tsv",
     sep = "\t", quote = FALSE,
     row.names = FALSE, col.names = FALSE
   )

**Verify the barcode file:**

.. code-block:: bash

   # Check format
   head barcode_celltype_map.tsv
   # AAACGAACAGTCAGTT-1	excitatory_neurons
   # AAACGAAGTCGCTCTA-1	inhibitory_neurons

   # Count cells per type
   cut -f2 barcode_celltype_map.tsv | sort | uniq -c | sort -rn
   #   2500 excitatory_neurons
   #   1800 inhibitory_neurons
   #   1200 astrocytes
   #    800 oligodendrocytes

Step 2: Run Per-Group Imbalance Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, analyze imbalance within each cell type:

.. code-block:: bash

   wasp2-analyze find-imbalance-sc \
     allele_counts.h5ad \
     barcode_celltype_map.tsv \
     --sample SAMPLE_ID \
     --phased \
     --min 10 \
     -z 3

This produces per-celltype result files (e.g., ``ai_results_excitatory_neurons.tsv``).

Step 3: Run Comparative Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Compare imbalance between specific cell types:

.. code-block:: bash

   # Compare excitatory vs inhibitory neurons
   wasp2-analyze compare-imbalance \
     allele_counts.h5ad \
     barcode_celltype_map.tsv \
     --sample SAMPLE_ID \
     --groups "excitatory_neurons,inhibitory_neurons" \
     --phased \
     --min 15

**Compare all pairwise combinations:**

.. code-block:: bash

   # Omit --groups to compare all cell types
   wasp2-analyze compare-imbalance \
     allele_counts.h5ad \
     barcode_celltype_map.tsv \
     --sample SAMPLE_ID \
     --phased \
     --min 15

This produces output files for each pairwise comparison:

* ``ai_results_excitatory_neurons_inhibitory_neurons.tsv``
* ``ai_results_excitatory_neurons_astrocytes.tsv``
* ``ai_results_inhibitory_neurons_astrocytes.tsv``
* ...

Step 4: Interpret Results
~~~~~~~~~~~~~~~~~~~~~~~~~

**Output columns explained:**

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - Column
     - Description
   * - ``region``
     - Genomic region (peak or gene) identifier
   * - ``num_snps``
     - Number of shared heterozygous SNPs used for comparison
   * - ``combined_mu``
     - Reference allele frequency under null hypothesis (shared between groups)
   * - ``mu1``
     - Reference allele frequency in group 1 (e.g., excitatory neurons)
   * - ``mu2``
     - Reference allele frequency in group 2 (e.g., inhibitory neurons)
   * - ``null_ll``
     - Log-likelihood under null hypothesis (shared μ)
   * - ``alt_ll``
     - Log-likelihood under alternative hypothesis (separate μ values)
   * - ``pval``
     - Likelihood ratio test p-value
   * - ``fdr_pval``
     - FDR-corrected p-value (Benjamini-Hochberg)

**Filtering significant results:**

.. code-block:: bash

   # Significant differential imbalance (FDR < 0.05)
   awk -F'\t' 'NR==1 || $9 < 0.05' ai_results_excitatory_neurons_inhibitory_neurons.tsv \
     > significant_differential_AI.tsv

   # Large effect size (>15% difference in allele frequency)
   awk -F'\t' 'NR==1 || ($4 - $5 > 0.15 || $5 - $4 > 0.15)' significant_differential_AI.tsv \
     > large_effect_differential_AI.tsv

**Interpret μ values:**

* ``mu < 0.5``: Alternate allele favored
* ``mu > 0.5``: Reference allele favored
* ``|mu1 - mu2| > 0.1``: Meaningful difference (~20% shift in allele preference)

Tutorial 2: Sex Differences Analysis
------------------------------------

Identify regions with sex-biased allelic imbalance in chromatin accessibility.

Step 1: Create Sex-Labeled Barcode Map
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import pandas as pd
   import scanpy as sc

   # Load your annotated data
   adata = sc.read_h5ad("processed_snATAC.h5ad")

   # Create barcode-to-sex mapping
   barcode_df = pd.DataFrame({
       'barcode': adata.obs_names,
       'sex': adata.obs['donor_sex']  # 'male' or 'female'
   })

   # Write without header
   barcode_df.to_csv('barcode_sex_map.tsv', sep='\t', header=False, index=False)

Step 2: Run Comparative Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   wasp2-analyze compare-imbalance \
     allele_counts.h5ad \
     barcode_sex_map.tsv \
     --sample SAMPLE_ID \
     --groups "male,female" \
     --phased \
     --min 20 \
     --out_file ai_results_sex_comparison.tsv

Step 3: Identify Sex-Biased Regions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Extract significant sex differences
   awk -F'\t' 'NR==1 || $9 < 0.01' ai_results_sex_comparison.tsv > sex_biased_regions.tsv

   # Count by chromosome (expect enrichment on X)
   cut -f1 sex_biased_regions.tsv | grep -E "^chr" | cut -d: -f1 | sort | uniq -c

Tutorial 3: Treatment vs Control
--------------------------------

Compare allelic imbalance before and after drug treatment.

Step 1: Prepare Condition Labels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import pandas as pd

   # Load metadata with treatment status
   metadata = pd.read_csv("sample_metadata.csv")

   # Create barcode-to-condition mapping
   barcode_df = pd.DataFrame({
       'barcode': metadata['cell_barcode'],
       'condition': metadata['treatment_status']  # 'treated' or 'control'
   })

   barcode_df.to_csv('barcode_treatment_map.tsv', sep='\t', header=False, index=False)

Step 2: Run Analysis
~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   wasp2-analyze compare-imbalance \
     allele_counts.h5ad \
     barcode_treatment_map.tsv \
     --sample SAMPLE_ID \
     --groups "treated,control" \
     --min 15 \
     --out_file ai_results_treatment.tsv

Step 3: Identify Treatment-Responsive Regions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import pandas as pd

   # Load results
   results = pd.read_csv('ai_results_treated_control.tsv', sep='\t')

   # Significant treatment effects
   significant = results[results['fdr_pval'] < 0.05]
   print(f"Found {len(significant)} treatment-responsive regions")

   # Direction of change
   treatment_gain = significant[significant['mu1'] > significant['mu2'] + 0.1]
   treatment_loss = significant[significant['mu2'] > significant['mu1'] + 0.1]

   print(f"Regions with increased ref allele in treatment: {len(treatment_gain)}")
   print(f"Regions with decreased ref allele in treatment: {len(treatment_loss)}")

Visualization Examples
----------------------

Volcano Plot
~~~~~~~~~~~~

.. code-block:: python

   import pandas as pd
   import matplotlib.pyplot as plt
   import numpy as np

   # Load results
   results = pd.read_csv('ai_results_excitatory_neurons_inhibitory_neurons.tsv', sep='\t')

   # Calculate effect size (difference in mu)
   results['effect_size'] = results['mu1'] - results['mu2']
   results['-log10_pval'] = -np.log10(results['pval'] + 1e-300)

   # Create volcano plot
   fig, ax = plt.subplots(figsize=(10, 8))

   # Non-significant points
   ns = results['fdr_pval'] >= 0.05
   ax.scatter(results.loc[ns, 'effect_size'],
              results.loc[ns, '-log10_pval'],
              c='gray', alpha=0.5, s=10, label='Not significant')

   # Significant points
   sig = results['fdr_pval'] < 0.05
   ax.scatter(results.loc[sig, 'effect_size'],
              results.loc[sig, '-log10_pval'],
              c='red', alpha=0.7, s=20, label='FDR < 0.05')

   ax.axhline(-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
   ax.axvline(0, color='black', linestyle='-', alpha=0.3)

   ax.set_xlabel('Effect Size (μ₁ - μ₂)', fontsize=12)
   ax.set_ylabel('-log₁₀(p-value)', fontsize=12)
   ax.set_title('Differential Allelic Imbalance:\nExcitatory vs Inhibitory Neurons', fontsize=14)
   ax.legend()

   plt.tight_layout()
   plt.savefig('differential_AI_volcano.png', dpi=150)

Heatmap of Top Hits
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import pandas as pd
   import seaborn as sns
   import matplotlib.pyplot as plt

   # Load results from multiple comparisons
   comparisons = [
       ('excitatory', 'inhibitory'),
       ('excitatory', 'astrocyte'),
       ('inhibitory', 'astrocyte'),
   ]

   # Collect mu values for top regions
   all_results = {}
   for g1, g2 in comparisons:
       df = pd.read_csv(f'ai_results_{g1}_{g2}.tsv', sep='\t')
       all_results[(g1, g2)] = df.set_index('region')

   # Find regions significant in any comparison
   sig_regions = set()
   for df in all_results.values():
       sig_regions.update(df[df['fdr_pval'] < 0.05].index[:20])  # Top 20 each

   # Build heatmap matrix (mu values per cell type)
   celltypes = ['excitatory', 'inhibitory', 'astrocyte']
   heatmap_data = pd.DataFrame(index=list(sig_regions), columns=celltypes)

   for region in sig_regions:
       for g1, g2 in comparisons:
           if region in all_results[(g1, g2)].index:
               row = all_results[(g1, g2)].loc[region]
               heatmap_data.loc[region, g1] = row['mu1']
               heatmap_data.loc[region, g2] = row['mu2']

   # Plot heatmap
   fig, ax = plt.subplots(figsize=(8, 12))
   sns.heatmap(heatmap_data.astype(float), cmap='RdBu_r', center=0.5,
               vmin=0, vmax=1, ax=ax, cbar_kws={'label': 'Ref Allele Frequency (μ)'})
   ax.set_title('Cell-Type-Specific Allelic Imbalance', fontsize=14)
   plt.tight_layout()
   plt.savefig('differential_AI_heatmap.png', dpi=150)

Command-Line Reference
----------------------

Full Parameter List
~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   wasp2-analyze compare-imbalance --help

   Usage: wasp2-analyze compare-imbalance [OPTIONS] ADATA BARCODE_MAP

   Arguments:
     ADATA        AnnData file with allele counts (.h5ad)
     BARCODE_MAP  TSV file mapping barcodes to groups

   Options:
     --groups TEXT        Comma-separated groups to compare (default: all)
     --min INTEGER        Minimum allele count per region per group (default: 10)
     --pseudocount INT    Pseudocount for zero counts (default: 1)
     --sample TEXT        Sample name for genotype filtering
     --phased             Use phased genotype information
     -z, --z_cutoff INT   Remove outlier SNPs above z-score threshold
     --out_file TEXT      Output file path

Best Practices
--------------

Data Quality
~~~~~~~~~~~~

* **Use WASP-filtered BAMs** to remove mapping bias artifacts
* **Require sufficient counts** (``--min 15`` or higher for robust estimates)
* **Apply z-score filtering** (``-z 3``) to remove outliers from CNVs or mapping artifacts

Statistical Power
~~~~~~~~~~~~~~~~~

* **Merge similar groups** if individual populations have low cell counts
* **Use phased genotypes** when available for improved power
* **Focus on regions with multiple SNPs** for more reliable estimates

Interpretation
~~~~~~~~~~~~~~

* **Biological replication** - Validate across independent samples
* **Effect size matters** - Consider the absolute difference between μ₁ and μ₂ alongside p-values
* **Integrate with eQTL data** - Connect to known regulatory variants
* **Orthogonal validation** - Confirm top hits with targeted methods

Common Issues
-------------

Low Power / Few Significant Results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Increase sequencing depth
* Merge similar cell types to increase counts per group
* Lower ``--min`` threshold (with caution)
* Use phased genotypes if available

Too Many Significant Results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Check for batch effects between groups
* Verify WASP filtering was applied
* Use stricter FDR threshold (e.g., 0.01)
* Check that groups have similar sequencing depth

Memory Issues
~~~~~~~~~~~~~

Process chromosomes separately:

.. code-block:: bash

   for chr in chr{1..22}; do
     wasp2-count count-variants-sc \
       sample.bam variants.vcf.gz barcodes.tsv \
       --region peaks_${chr}.bed \
       --out_file counts_${chr}.h5ad

     wasp2-analyze compare-imbalance \
       counts_${chr}.h5ad \
       barcode_celltype_map.tsv \
       --out_file results_${chr}.tsv
   done

See Also
--------

* :doc:`/user_guide/analysis` - Statistical methods and parameters
* :doc:`/user_guide/single_cell` - Single-cell data formats
* :doc:`/tutorials/scrna_seq` - Basic scRNA-seq tutorial
