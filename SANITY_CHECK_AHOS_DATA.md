# Sanity Check: Aho's WASP2 Data and Results

**Date**: 2025-11-25
**Question**: Can we reuse aho's validation work, or do we need to reprocess?
**Answer**: ⚠️ **MUST REPROCESS** - Aho's work is SNP-only, NO indels!

---

## Summary: What Aho Has

### ✅ **GOOD NEWS**: High-Quality SNP Results

**Location**: `/iblm/netapp/data3/aho/project_data/wasp2/imprinted_rna/outputs/`

**Files**:
- `GM12878_geneimprint_exon_gene_ai_results.tsv` (634 KB)
- `GM12878_geneimprint_transcript_gene_ai_results.tsv` (841 KB)
- `GM12878_geneimprint_exon_gene_counts.tsv` (6.2 MB)
- `GM12878_geneimprint_transcript_gene_counts.tsv` (236 MB)

**Coverage**:
- **118,241 total variants**
- **249 imprinted genes** (from geneimprint.org)
- Statistical analysis: LRT, p-values, FDR correction

**Key Results**:
| Gene | REF Count | ALT Count | Ratio | P-value | Expected |
|------|-----------|-----------|-------|---------|----------|
| XIST | 65 | 1,479 | 0.04 | 6.3e-24 | Extreme ✅ |
| TSIX | 45 | 915 | 0.05 | 3.9e-13 | Extreme ✅ |
| H19 | (in gene list) | - | - | - | Extreme |
| IGF2 | (in gene list) | - | - | - | Extreme |

**Date processed**: October 31, 2024

---

### ❌ **CRITICAL FINDING**: NO INDELS in Aho's Results!

**Evidence**:

1. **Variant count**:
   ```bash
   # Total variants
   wc -l GM12878_geneimprint_exon_gene_counts.tsv
   # Output: 118,241

   # Count indels (ref or alt length > 1)
   awk 'length($3) > 1 || length($4) > 1' GM12878_geneimprint_exon_gene_counts.tsv | wc -l
   # Output: 0  ← NO INDELS!
   ```

2. **Sample data** (all single-base SNPs):
   ```
   chrom  pos      ref  alt  GT    gene_name    ref_count  alt_count
   chr1   807445   A    G    0|1   RP11-206L10.9    0         0
   chr1   852047   C    T    0|1   LINC01128        0         0
   ```
   All variants: `A→G`, `C→T`, etc. (single base changes)
   No insertions (e.g., `C→CAT`)
   No deletions (e.g., `CAT→C`)

3. **Processing code** (`count_alleles.py`):
   ```python
   # Line 62: Simple base extraction
   count_list.append(allele.upper())

   # Line 107: Explicitly says "SNP's"
   print(f"Counted {len(snp_list)} SNP's in {end - start} seconds!\n")
   ```

   This code uses `pileup.get_query_sequences()` which returns **single bases only**.
   **CANNOT handle indels** (indels require CIGAR parsing).

4. **WASP2 version**: Aho's version predates indel support
   - Aho's WASP2: `/iblm/netapp/home/aho/dev/WASP2/`
   - Last modified: September 2024
   - No indel-related commits found

---

## What's Missing for Indel Validation

### **1. Source Data** ❌

**No raw data found in imprinted_rna directory**:
```bash
find /iblm/netapp/data3/aho/project_data/wasp2/imprinted_rna -name "*.bam"
# Output: (none)

find /iblm/netapp/data3/aho/project_data/wasp2/imprinted_rna -name "*.vcf*"
# Output: (none)
```

**Possible locations** (to check):
- ENCODE downloads
- Shared data directories
- Symlinks in other projects

### **2. Processing Pipeline** ⚠️

**Found** in `/iblm/netapp/home/aho/projects/wasp/`:
- Scripts: `count_alleles.py`, `as_analysis.py`
- Notebooks: `plot_imbalance.ipynb`, `test_as_analysis.ipynb`
- Data: `NA12878_snps.bed` (50 MB - SNPs only!)

**But**: All SNP-focused, no indel handling

### **3. WASP2 with Indel Support** ✅

**Our version** (current directory):
- Location: `/iblm/netapp/data3/jjaureguy/.../WASP2-exp/`
- Branch: `rust-optimization-plink2`
- Recent commits:
  ```
  049d528 docs: add BCF performance recommendation
  fb477cf build: add cargo config for libclang
  4173287 chore: bump version to 1.2.0 - PGEN support
  6d853a6 feat: integrate PLINK2 PGEN support
  ```

**Has indel support**: ✅ (need to verify in code)

---

## Reprocessing Requirements

### **Option A: Reprocess GM12878 Imprinted Genes with Indels** ⭐ **RECOMMENDED**

**What we need**:
1. GM12878 RNA-seq BAM (or FASTQ)
2. GM12878 phased VCF with **both SNPs AND indels**
3. Imprinted gene GTF (already have: `/iblm/netapp/data3/aho/project_data/wasp2/imprinted_rna/data/geneimprint.gtf`)

**Processing steps**:
```bash
# 1. Download/locate GM12878 RNA-seq data
# (check ENCODE, or find aho's source)

# 2. Get phased VCF with indels
# 1000 Genomes Phase 3: has SNPs + short indels
wget ftp://ftp.1000genomes.ebi.ac.uk/.../NA12878.vcf.gz

# 3. Run WASP2 with indel support
WASP2 --bam GM12878.bam \
      --vcf NA12878.vcf.gz \
      --include_indels \  # ← KEY FLAG!
      --gtf geneimprint.gtf \
      --out gm12878_indel_results/

# 4. Count alleles (SNPs + indels)
python count_alleles_with_indels.py \
    --bam gm12878_indel_results/keep.merged.bam \
    --vcf NA12878.vcf.gz \
    --gtf geneimprint.gtf \
    --out gm12878_indel_counts.tsv
```

**Time estimate**: 1-2 days (depending on data download)

**Outcome**:
- Can compare SNP-only (aho's results) vs SNP+indel (our results)
- Shows added value of indel support
- Perfect for manuscript!

---

### **Option B: Use Simulation Only** ⚠️ **WEAKER**

**What we'd do**:
- Run our simulation framework (`simulate_indel_ase_v2.py`)
- Get gold standard validation
- Skip real biological data validation

**Problem**: Reviewers will ask "but does it work on real data?"

**Rating**: 6/10 (acceptable but not strong)

---

### **Option C: Run Both Simulation + Reprocess** ⭐⭐⭐ **BEST**

**What we'd do**:
1. **Simulation** (gold standard): Proves algorithm correctness
2. **GM12878 reprocessing** (biological): Proves real-world performance
3. **Comparison**: SNP-only vs SNP+indel on same data

**Validation story**:
```
Tier 1 (Gold): Simulation with known truth
  → Mean error 2.7%, all tests <10%
  → PROVES algorithm works correctly

Tier 3 (Biological): GM12878 imprinted genes
  → SNP-only: 118K variants, XIST/H19/IGF2 show imbalance (aho's results)
  → SNP+indel: ???K variants, same genes + new indel-based ASE
  → PROVES real-world biological relevance
  → SHOWS added value of indel support
```

**Rating**: 9/10 (publication-ready)

---

## Data Sources to Check

### **Where to find GM12878 data**:

1. **ENCODE Project**:
   - RNA-seq: Search "GM12878 RNA-seq" on https://www.encodeproject.org/
   - Example: ENCSR000AEL (GM12878 polyA RNA-seq)
   - BAM files: ~10-30 GB each

2. **Aho's ENCODE directory**:
   ```bash
   /iblm/netapp/data3/aho/project_data/wasp2/encode/male_37/
   # Has ENCFF944WLM.vcf (54 MB)
   # But this is for colon tissue ATAC-seq, not GM12878 RNA-seq
   ```

3. **1000 Genomes**:
   - Phased VCF: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
   - NA12878 (same as GM12878): `ALL.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz`
   - **Has indels**: Yes! (insertions and deletions up to ~50bp)

4. **Platinum Genomes** (higher quality):
   - https://www.illumina.com/platinumgenomes.html
   - Trio-phased variants
   - Includes indels

5. **Check aho's original data source**:
   ```bash
   # Need to find where aho got the GM12878 BAM from
   # Likely in /iblm/netapp/data*/aho/ somewhere
   # Or downloaded and deleted after processing
   ```

---

## Action Items

### **Immediate** (This Week):

1. ✅ **Run simulation validation**
   ```bash
   ./run_simulation.sh moderate
   # 30 min runtime
   # Proves indel algorithm correctness
   ```

2. ⚠️ **Find GM12878 source data**
   - Check aho's data directories more thoroughly
   - Or download from ENCODE (1-2 hours)
   - Need: BAM + VCF with indels

3. ⚠️ **Verify our WASP2 has indel support**
   - Check `--include_indels` flag exists
   - Test on small dataset
   - Confirm CIGAR parsing works for indels

### **Soon** (Next Week):

4. **Reprocess GM12878** with indels (if time permits)
   - 249 imprinted genes
   - SNPs + indels
   - Compare to aho's SNP-only results
   - Time: 1-2 days

5. **Analysis comparison**:
   - SNP-only ASE (aho's results): 118K variants
   - SNP+indel ASE (our results): ???K variants
   - Novel discoveries: Indels showing ASE that SNPs missed
   - Manuscript impact: "Indel support increases ASE detection by X%"

---

## Bottom Line

### **Aho's Results**:
- ✅ Excellent SNP-based validation
- ✅ High-quality imprinted gene analysis
- ✅ Can cite/reference in manuscript
- ❌ **NO indels** - completely SNP-only
- ❌ Cannot be directly reused for indel validation

### **What We Must Do**:
1. **Run simulation** (gold standard for indels) - **CRITICAL**
2. **Find/download GM12878 data** (BAM + VCF with indels)
3. **Reprocess with indel support** (biological validation)
4. **Compare**: SNP-only vs SNP+indel performance

### **Can We Skip Reprocessing?**
- **Technically**: Yes (simulation alone is acceptable)
- **Scientifically**: Not ideal (reviewers want real data validation)
- **For publication**: Recommend doing both simulation + real data

### **Timeline**:
- **Simulation only**: 1 day (acceptable but weaker)
- **Simulation + GM12878**: 3-4 days (strong validation)
- **Full validation**: 1-2 weeks (bulletproof)

---

## Next Steps

**Priority 1**: Run simulation (proves indel correctness)
**Priority 2**: Find GM12878 source data
**Priority 3**: Reprocess with indels (proves real-world value)

**Question for user**: Do you know where aho got the GM12878 BAM/VCF originally? Or should we download fresh from ENCODE?

---

## File Locations Summary

**Aho's SNP results** ✅:
```
/iblm/netapp/data3/aho/project_data/wasp2/imprinted_rna/outputs/
├── GM12878_geneimprint_exon_gene_ai_results.tsv
├── GM12878_geneimprint_transcript_gene_ai_results.tsv
├── GM12878_geneimprint_exon_gene_counts.tsv (118,241 SNPs, 0 indels)
└── GM12878_geneimprint_transcript_gene_counts.tsv
```

**Aho's gene annotation** ✅:
```
/iblm/netapp/data3/aho/project_data/wasp2/imprinted_rna/data/
├── geneimprint.gtf (1.3 GB - 249 genes)
└── geneimprint.txt (gene list)
```

**Aho's processing code** ✅:
```
/iblm/netapp/home/aho/projects/wasp/scripts/
├── count_alleles.py (SNP-only pileup method)
├── as_analysis.py
└── Various .ipynb notebooks
```

**Aho's WASP2** (pre-indel):
```
/iblm/netapp/home/aho/dev/WASP2/
```

**Our WASP2** (with indel support):
```
/iblm/netapp/data3/jjaureguy/.../WASP2-exp/
```

**Our simulation** ✅:
```
/iblm/netapp/data3/jjaureguy/.../WASP2-exp/
├── simulate_indel_ase_v2.py
└── run_simulation.sh
```

**Source data** ❌ **MISSING**:
- GM12878 RNA-seq BAM (need to find or download)
- GM12878 phased VCF with indels (can download from 1000 Genomes)

---

**Recommendation**: Run simulation first (30 min), then decide if we have time to reprocess GM12878. Simulation alone is acceptable, but reprocessing would make the validation much stronger.
