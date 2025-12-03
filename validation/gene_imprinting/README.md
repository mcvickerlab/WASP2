# Gene Imprinting Validation

This directory contains the gene imprinting validation framework for WASP2 indel support.
Imprinted genes serve as biological positive controls - they should show extreme allelic
ratios (~100:0) due to parent-of-origin silencing.

## Known Imprinted Genes Tested

- **H19** (chr11) - Maternally expressed
- **IGF2** (chr11) - Paternally expressed
- **SNRPN** (chr15) - Paternally expressed
- **PEG3** (chr19) - Paternally expressed
- **PLAGL1** (chr6) - Paternally expressed
- **CDKN1C** (chr11) - Maternally expressed
- **MAGEL2** (chr15) - Paternally expressed
- **MEST** (chr7) - Paternally expressed
- **UBE3A** (chr15) - Maternally expressed in brain
- **DLK1** (chr14) - Paternally expressed
- **GNAS** (chr20) - Complex imprinting
- **TP73** (chr1) - Maternally expressed

## Data Locations

### Aaron's Original Validation Data
- **BAM**: `/iblm/netapp/data3/aho/alignment/GM12878_rna_v2/GM12878_merged.sorted.bam`
- **VCF**: `/iblm/netapp/data1/aho/variants/NA12878.vcf.gz` (573,836 indels)
- **GTF**: `/iblm/netapp/data3/aho/project_data/wasp2/imprinted_rna/data/geneimprint.gtf`
- **Reference**: `/iblm/netapp/data1/aho/ref_genomes/index/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/genome.fa`

### Aaron's Imprinting Notebook
- **Notebook**: `/iblm/netapp/home/aho/projects/wasp/testing/performance/test_imprinted.ipynb`
- **ATAC-seq results**: `/iblm/netapp/home/aho/projects/wasp/testing/performance/data/GM12878_ATACseq_50k_merged`
- **Plotting notebooks**: `/iblm/netapp/home/aho/projects/wasp/testing/performance/plot_ai_results.ipynb`

## Validation Approach

1. **Intersect variants with imprinted gene regions**
   - Extract heterozygous SNPs and indels in imprinted genes

2. **Run WASP2 allelic imbalance analysis**
   - Process with full WASP2 pipeline (remap + filter)
   - Count alleles at each variant position

3. **Compare imprinted vs control genes**
   - Imprinted genes: Expected ~100:0 ratios (extreme AI)
   - Control genes: Expected ~50:50 ratios (no AI)

4. **Validate indel handling**
   - Compare SNP-only vs SNP+indel results
   - Ensure indels in imprinted genes show same extreme ratios

## Quick Start

```bash
# Run the comparison script
python gm12878_benchmark/scripts/compare_to_aaron.py

# Run the full pipeline test (3 genes)
bash gm12878_benchmark/scripts/run_full_pipeline_tier1.sh
```

## Expected Results

From Aaron's validation:
- Imprinted genes show reference allele fractions of ~0.0 or ~1.0
- Control genes show fractions centered around 0.5
- p-values for imprinted genes are highly significant

## Related Files

- `../gm12878_benchmark/` - Full benchmark suite with indel validation
- `../../simulate_indel_ase_v2.py` - Synthetic indel ASE simulation
- `../../BIOLOGICAL_VALIDATION.md` - Documentation of validation strategy
