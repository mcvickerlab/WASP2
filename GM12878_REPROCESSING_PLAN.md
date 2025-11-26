# GM12878 Imprinted Gene Reprocessing with Indels

**Status**: ✅ **ALL DATA FOUND - READY TO PROCESS!**

**Date**: 2025-11-25

---

## Executive Summary

Found Aaron's complete GM12878 dataset used for October 2024 imprinted gene analysis:

✅ **RNA-seq BAM**: 3.0 GB, 52.3M reads, GRCh38
✅ **Phased VCF**: 36 MB, **573,836 indels** + 3.6M SNPs, Illumina Platinum Genomes
✅ **Gene annotation**: 1.3 GB GTF, 249 imprinted genes
✅ **Processing notebook**: Exact commands used by Aaron

**Aaron's results**: SNP-only (118,241 variants, 0 indels)
**Our goal**: Reprocess with indels using `--include_indels` flag

---

## Data Files (All Found!)

### 1. RNA-seq Data ✅

**BAM File**:
```bash
/iblm/netapp/data3/aho/alignment/GM12878_rna_v2/GM12878_merged.sorted.bam
```

**Details**:
- Size: 3.0 GB
- Reads: 52,299,874 mapped (100% mapping rate)
- Reference: GRCh38 (NCBI)
- Quality: MAPQ ≥30, duplicates removed, no chrM
- Date: Feb 22, 2022

**BAM Index**:
```bash
/iblm/netapp/data3/aho/alignment/GM12878_rna_v2/GM12878_merged.sorted.bam.bai
```

**Source**: 11 SRA runs merged (SRR306998-SRR307012)
**GEO Accession**: Likely GSM754335

---

### 2. Phased VCF ✅

**VCF File**:
```bash
/iblm/netapp/data1/aho/variants/NA12878.vcf.gz
```

**Details**:
- Size: 36 MB compressed
- Reference: GRCh38 (hg38)
- Source: **Illumina Platinum Genomes v2016-1.0**
- Total variants: 4,167,900
  - **SNPs**: 3,596,231
  - **Indels**: 573,836 ✅
- Phasing: YES (family-based, trio)
- Quality: Multi-method consensus (bwa_freebayes, bwa_gatk, bwa_platypus, isaac_strelka)

**VCF Index**:
```bash
/iblm/netapp/data1/aho/variants/NA12878.vcf.gz.tbi
```

**Pedigree**: NA12878 (daughter) = NA12892 (mother) + NA12891 (father)

**Indel Examples** (from VCF):
```
chr1:788418   CAG→C        (2bp deletion)
chr1:789568   TATGGA→T     (5bp deletion)
chr1:789680   A→AGGAATGGAATGGAAT  (15bp insertion)
chr1:790933   CGAATGGAATG→C      (10bp deletion)
chr1:822427   36bp→1bp     (35bp deletion)
```

---

### 3. Gene Annotation ✅

**GTF File**:
```bash
/iblm/netapp/data3/aho/project_data/wasp2/imprinted_rna/data/geneimprint.gtf
```

**Details**:
- Size: 1.3 GB
- Genes: 249 imprinted genes
- Source: geneimprint.org
- Features: exons and transcripts

**Gene List**:
```bash
/iblm/netapp/data3/aho/project_data/wasp2/imprinted_rna/data/geneimprint.txt
```

**Classic genes included**: H19, IGF2, SNRPN, XIST, TSIX, KCNQ1OT1, MEG3, DLK1, GNAS, etc.

---

## Aaron's Processing Pipeline (SNP-only)

From his notebook: `/iblm/netapp/home/aho/projects/wasp/testing/performance/test_imprinted.ipynb`

### Step 1: Count Variants (SNPs only)

```bash
test_vcf="/iblm/netapp/data1/aho/variants/NA12878.vcf.gz"
test_bam="/iblm/netapp/data3/aho/alignment/GM12878_rna_v2/GM12878_merged.sorted.bam"
filt_gtf="/iblm/netapp/data3/aho/project_data/wasp2/imprinted_rna/data/geneimprint.gtf"

test_feature="transcript"  # or "exon"
parent_attribute="gene_name"
test_attribute="transcript_id"
test_sample="NA12878"

outdir="/iblm/netapp/data3/aho/project_data/wasp2/imprinted_rna/outputs"
out_file="${outdir}/GM12878_geneimprint_${test_feature}_gene_counts.tsv"

python /iblm/netapp/home/aho/dev/WASP2/src/counting count-variants \
    ${test_bam} ${test_vcf} \
    -s ${test_sample} \
    -r ${filt_gtf} \
    -o ${out_file} \
    --gene_feature ${test_feature} \
    --gene_attribute ${test_attribute} \
    --gene_parent ${parent_attribute}
    # NOTE: NO --include_indels FLAG!
```

**Result**: 118,241 variants counted (SNPs only, 0 indels)

### Step 2: Find Allelic Imbalance (SNPs only)

```bash
ai_file="${outdir}/GM12878_geneimprint_${test_feature}_gene_ai_results.tsv"

python /iblm/netapp/home/aho/dev/WASP2/src/analysis find-imbalance \
    --phased \
    --out ${ai_file} \
    --group ${parent_attribute} \
    ${out_file}
```

**Result**: Statistical analysis, LRT, p-values, FDR correction
**Key findings**: XIST (65 REF vs 1,479 ALT), H19, IGF2, SNRPN show expected imbalance

---

## Our Reprocessing Pipeline (WITH INDELS)

### Output Directory

```bash
outdir="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/gm12878_with_indels"
mkdir -p ${outdir}
```

### Step 1: Count Variants (SNPs + Indels) ⭐

```bash
# Input files (same as Aaron's)
test_vcf="/iblm/netapp/data1/aho/variants/NA12878.vcf.gz"
test_bam="/iblm/netapp/data3/aho/alignment/GM12878_rna_v2/GM12878_merged.sorted.bam"
filt_gtf="/iblm/netapp/data3/aho/project_data/wasp2/imprinted_rna/data/geneimprint.gtf"

# Parameters
test_feature="transcript"  # or "exon"
parent_attribute="gene_name"
test_attribute="transcript_id"
test_sample="NA12878"

# Output
out_file="${outdir}/GM12878_geneimprint_transcript_gene_counts_WITH_INDELS.tsv"

# Run counting WITH INDELS
python src/counting count-variants \
    ${test_bam} ${test_vcf} \
    -s ${test_sample} \
    -r ${filt_gtf} \
    -o ${out_file} \
    --gene_feature ${test_feature} \
    --gene_attribute ${test_attribute} \
    --gene_parent ${parent_attribute} \
    --include_indels  # ⭐ KEY DIFFERENCE!
```

**Expected output**:
- SNPs: ~118K (same as Aaron)
- Indels: ???K (new!)
- Total: ~118K + ???K variants

### Step 2: Find Allelic Imbalance (SNPs + Indels)

```bash
ai_file="${outdir}/GM12878_geneimprint_transcript_gene_ai_results_WITH_INDELS.tsv"

python src/analysis find-imbalance \
    --phased \
    --out ${ai_file} \
    --group ${parent_attribute} \
    ${out_file}
```

### Step 3: Comparison Analysis

```python
import pandas as pd

# Load Aaron's SNP-only results
snp_only = pd.read_csv(
    "/iblm/netapp/data3/aho/project_data/wasp2/imprinted_rna/outputs/GM12878_geneimprint_transcript_gene_ai_results.tsv",
    sep="\t"
)

# Load our SNP+indel results
snp_indel = pd.read_csv(
    "${outdir}/GM12878_geneimprint_transcript_gene_ai_results_WITH_INDELS.tsv",
    sep="\t"
)

# Load counts to check indel contribution
counts_snp_only = pd.read_csv(
    "/iblm/netapp/data3/aho/project_data/wasp2/imprinted_rna/outputs/GM12878_geneimprint_transcript_gene_counts.tsv",
    sep="\t"
)

counts_snp_indel = pd.read_csv(
    "${outdir}/GM12878_geneimprint_transcript_gene_counts_WITH_INDELS.tsv",
    sep="\t"
)

# Count indels
# (check variant position format - may need CIGAR parsing)

# Compare results
print(f"SNP-only variants: {len(counts_snp_only)}")
print(f"SNP+indel variants: {len(counts_snp_indel)}")
print(f"New indel-based variants: {len(counts_snp_indel) - len(counts_snp_only)}")

# Compare gene-level AI results
merged = snp_only.merge(snp_indel, on="gene_name", suffixes=("_snp", "_indel"))
print(f"\nGenes with AI (SNP-only): {len(snp_only)}")
print(f"Genes with AI (SNP+indel): {len(snp_indel)}")

# Novel discoveries
novel = snp_indel[~snp_indel["gene_name"].isin(snp_only["gene_name"])]
print(f"Novel AI genes (indel-driven): {len(novel)}")
```

---

## Validation Questions

### Can We Identify Indel-Specific ASE?

**Approach**:
1. Compare variant counts per gene (SNP-only vs SNP+indel)
2. Find genes where indels contribute >20% of ASE signal
3. Check if novel AI genes have mostly indel variants

**Key genes to check**:
- H19: Known imprinted, check if indels add ASE signal
- IGF2: Paternally imprinted, may have regulatory indels
- XIST: X-inactivation, extreme maternal bias
- SNRPN: Prader-Willi locus

### Performance Comparison

**Metrics**:
- Total variants detected: SNP-only vs SNP+indel
- Genes with significant AI: SNP-only vs SNP+indel
- Concordance: Do both methods agree on classic genes (H19, IGF2, etc.)?
- Novel discoveries: Genes showing AI only when indels included

---

## Timeline

### Quick Run (1-2 hours):
1. Run count-variants with `--include_indels` (~30-60 min)
2. Run find-imbalance analysis (~10-20 min)
3. Quick comparison to Aaron's results (~10 min)

### Full Analysis (1 day):
4. Detailed comparison script (2-3 hours)
5. Statistical validation (2-3 hours)
6. Generate comparison figures (1-2 hours)
7. Write up results for manuscript (2-3 hours)

---

## Expected Manuscript Impact

### Figure: SNP-only vs SNP+indel Comparison

**Panel A**: Variant counts by gene
- X-axis: Gene (sorted by total variants)
- Y-axis: Variant count
- Two bars: SNP-only (blue), SNP+indel (red)
- Shows indel contribution per gene

**Panel B**: ASE concordance
- Scatter plot: SNP-only ratio vs SNP+indel ratio
- Points should cluster near diagonal (concordance)
- Highlight genes where indels change result

**Panel C**: Novel discoveries
- List of genes showing AI only with indels
- Biological relevance (if known)

### Text:

> "To validate indel support, we reprocessed 249 imprinted genes in GM12878 lymphoblastoid cells (52.3M RNA-seq reads) using Illumina Platinum Genomes phased variants (573K indels + 3.6M SNPs). Including indels increased variant detection by X% and identified Y additional genes with allelic imbalance. Classic imprinted genes (H19, IGF2, SNRPN, XIST) showed consistent ASE patterns with both SNP-only and SNP+indel analysis, validating indel support does not introduce artifacts. Indels contributed >20% of ASE signal in Z genes, demonstrating the importance of including structural variants in allelic imbalance analysis."

---

## Processing Script (Complete)

Save as: `reprocess_gm12878_with_indels.sh`

```bash
#!/bin/bash
set -e

# Input files
BAM="/iblm/netapp/data3/aho/alignment/GM12878_rna_v2/GM12878_merged.sorted.bam"
VCF="/iblm/netapp/data1/aho/variants/NA12878.vcf.gz"
GTF="/iblm/netapp/data3/aho/project_data/wasp2/imprinted_rna/data/geneimprint.gtf"

# Output directory
OUTDIR="gm12878_with_indels"
mkdir -p ${OUTDIR}

# Parameters
SAMPLE="NA12878"
FEATURE="transcript"
ATTRIBUTE="transcript_id"
PARENT="gene_name"

# Step 1: Count variants (WITH INDELS!)
echo "Step 1: Counting variants (SNPs + indels)..."
COUNTS="${OUTDIR}/GM12878_geneimprint_transcript_gene_counts_WITH_INDELS.tsv"

python src/counting count-variants \
    ${BAM} ${VCF} \
    -s ${SAMPLE} \
    -r ${GTF} \
    -o ${COUNTS} \
    --gene_feature ${FEATURE} \
    --gene_attribute ${ATTRIBUTE} \
    --gene_parent ${PARENT} \
    --include_indels

echo "Done! Variants counted: $(wc -l < ${COUNTS})"

# Step 2: Find allelic imbalance
echo "Step 2: Finding allelic imbalance..."
AI_RESULTS="${OUTDIR}/GM12878_geneimprint_transcript_gene_ai_results_WITH_INDELS.tsv"

python src/analysis find-imbalance \
    --phased \
    --out ${AI_RESULTS} \
    --group ${PARENT} \
    ${COUNTS}

echo "Done! AI results: ${AI_RESULTS}"

# Step 3: Compare to SNP-only
echo "Step 3: Comparing to SNP-only results..."
SNP_ONLY="/iblm/netapp/data3/aho/project_data/wasp2/imprinted_rna/outputs/GM12878_geneimprint_transcript_gene_counts.tsv"

echo "SNP-only variants: $(wc -l < ${SNP_ONLY})"
echo "SNP+indel variants: $(wc -l < ${COUNTS})"
echo "Difference: $(($(wc -l < ${COUNTS}) - $(wc -l < ${SNP_ONLY})))"

echo "All done! Results in ${OUTDIR}/"
```

**Run it**:
```bash
chmod +x reprocess_gm12878_with_indels.sh
./reprocess_gm12878_with_indels.sh
```

---

## Verification Checklist

Before running:
- [ ] Verify BAM exists and is readable: `samtools view -H ${BAM} | head`
- [ ] Verify VCF exists and has indels: `zcat ${VCF} | grep -v "^#" | awk 'length($4)>1 || length($5)>1' | head`
- [ ] Verify GTF exists: `head ${GTF}`
- [ ] Check WASP2 has `--include_indels` flag: `python src/counting count-variants --help | grep indel`

After running:
- [ ] Check variant counts increased (should be >118K)
- [ ] Check classic genes still show AI (H19, IGF2, SNRPN, XIST)
- [ ] Identify genes with indel-driven ASE
- [ ] Generate comparison figures
- [ ] Write manuscript text

---

## Bottom Line

✅ **ALL DATA FOUND AND READY**
✅ **VCF HAS 573,836 INDELS**
✅ **EXACT COMMANDS DOCUMENTED**
✅ **CAN REPROCESS IN 1-2 HOURS**

**Next step**: Run `./reprocess_gm12878_with_indels.sh` to validate indel support on real biological data!

This will provide strong validation for the manuscript:
- Simulation (gold standard) ← Already built
- Real biological data (GM12878 imprinted genes) ← About to run
- Comparison (SNP-only vs SNP+indel) ← Will demonstrate value

**Publication-ready validation in ~1 day!** 🎉
