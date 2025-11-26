# Existing ASE Simulators and Validation Approaches

**Question**: Are there existing simulators we can use instead of building our own?

**TL;DR**: Yes, but they're either **overkill** or **don't handle indels well**. A simple custom simulator is actually the best approach for WASP2.

---

## Existing RNA-seq/ASE Simulators

### 1. **CAMPAREE** ⭐ **Most Relevant**

**Source**: [CAMPAREE: BMC Genomics 2021](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-021-07934-2)

**What it does**:
- Configurable And Modular Program Allowing RNA Expression Emulation
- Simulates diploid RNA samples at the molecule level
- **Supports sample-dependent features**: indels, SNVs consistent across reads
- **Realistic allele-specific expression**
- Can generate reads from both parental alleles

**Pros**:
- ✅ Explicitly supports ASE
- ✅ Handles indels
- ✅ Uses empirical data for realism
- ✅ Can work with BEERS2

**Cons**:
- ⚠️ Complex to set up
- ⚠️ Requires learning new tool
- ⚠️ May be overkill for simple validation
- ⚠️ Not Python-native (harder to integrate)

**Verdict**: **Could use, but probably overkill**

---

### 2. **BEERS2** (RNA-Seq Benchmark for Expression)

**Source**: [BEERS2: Briefings in Bioinformatics 2024](https://academic.oup.com/bib/article/25/3/bbae164/7644138)

**What it does**:
- High fidelity in silico RNA-seq modeling
- When combined with CAMPAREE: supports genomic variants and diploid expression
- Simulates sequencing errors, biases

**Pros**:
- ✅ Very realistic
- ✅ With CAMPAREE: supports ASE

**Cons**:
- ⚠️ Complex workflow
- ⚠️ Requires CAMPAREE for ASE
- ⚠️ Heavy tool for simple task

**Verdict**: **Too complex for our needs**

---

### 3. **Polyester** (R Package)

**Source**: [Polyester: Bioinformatics 2015](https://academic.oup.com/bioinformatics/article/31/17/2778/183245)
**GitHub**: [alyssafrazee/polyester](https://github.com/alyssafrazee/polyester)

**What it does**:
- Simulates RNA-seq reads with differential expression
- Easy to use (one R function call)
- Models technical biases

**Pros**:
- ✅ Easy to use
- ✅ Well-documented
- ✅ Widely used

**Cons**:
- ❌ **Does NOT natively support diploid/ASE**
- ❌ **Does NOT handle indels well**
- ⚠️ Requires R (not Python)
- ⚠️ Would need CAMPAREE wrapper for ASE

**Verdict**: **Not suitable without modification**

---

### 4. **SymSim** (Single-Cell)

**Source**: [SymSim: Nature Communications 2019](https://www.nature.com/articles/s41467-019-10500-w)

**What it does**:
- Single-cell RNA-seq simulator
- Models allelic variation

**Verdict**: **Not relevant (single-cell specific)**

---

### 5. **SimSeq**

**Source**: [SimSeq: Bioinformatics 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4481850/)

**What it does**:
- Nonparametric RNA-seq simulator
- Uses real data distributions

**Verdict**: **No ASE support**

---

## How Original WASP Validated (2015)

**Source**: [WASP: Nature Methods 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4626402/)

### **Their Simulation Approach**:

**Quote from paper**:
> "We simulated all possible overlapping reads from both haplotypes at heterozygous SNPs in GM12878 (a completely genotyped and phased lymphoblastoid cell line), allowing reads to contain mismatches at a predefined sequencing error rate."

**What they did**:
1. **Used real genotypes** (GM12878 - fully phased)
2. **Generated synthetic reads** at each heterozygous SNP
3. **Created reads from both haplotypes** with controlled ratios
4. **Tested mapping bias**:
   - N-masked genome: **biased** (false positives)
   - Personalized genome (AlleleSeq): **biased**
   - WASP: **almost perfectly balanced** ✅

**Statistical validation**:
- Simulated 100bp reads under null (OR=1) and alternative (OR>1) models
- 90% null sites, 10% alternative sites
- Tested if WASP correctly identifies imbalanced sites

### **Key Insight**: WASP's original validation was **simple read generation**, not complex simulation

---

## ASE Tools and Their Validation

### **ASElux**
**Source**: [ASElux: Bioinformatics 2018](https://academic.oup.com/bioinformatics/article/34/8/1313/4653700)

- Counts allelic reads directly from unmapped data
- **Validation**: Compared to other aligners on simulated data
- **No custom simulator** - used existing tools

### **MBASED**
**Source**: [MBASED: Genome Biology 2014](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0405-3)

- Meta-analysis based ASE detection
- **Validation**: Real data + known imprinted genes
- **No simulation** - biological controls

### **ASEReadCounter (GATK)**
- Standard tool for counting alleles
- **Validation**: Part of GATK best practices
- Uses real data validation

---

## Indel-Specific Challenges

**From**: [Tools and Best Practices: Genome Biology 2015](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0762-6)

**Key findings**:
> "Small indels are another mechanism of allelic expression, but they tend to cause an alignment error leading to bias in ASE estimation."

> "Simulated data indicates that **12% of SNPs** and **46% of indels** had allele ratio bias >5%"

**Common practice**:
- Most tools **filter out indels** near SNPs
- Tools like VADT: "filters out indels, SNPs within a certain distance of an indel"
- **Why?** Indels cause alignment errors

**This is exactly why WASP2 indel support is novel and needs validation!**

---

## Recommended Approach for WASP2

### **Option 1: Simple Custom Simulation** ⭐⭐⭐ **RECOMMENDED**

**What**: Your `simulate_indel_ase.py` script

**Pros**:
- ✅ Simple (< 300 lines)
- ✅ Python-native (easy integration)
- ✅ Tailored to WASP2 needs
- ✅ Tests exactly what you need (indel position mapping, quality)
- ✅ Fast to implement (already done!)
- ✅ Easy for reviewers to understand
- ✅ **Follows WASP's original approach** (simple read generation)

**Cons**:
- Not as "realistic" as CAMPAREE (but doesn't matter for validation)

**Why this is fine**:
- WASP original paper used simple simulation ✅
- You're testing bias removal, not modeling biology
- Ground truth is what matters, not realism
- Reviewers will understand it easily

---

### **Option 2: Use CAMPAREE** ⚠️

**When to use**: If reviewers demand "industry standard" simulator

**How**:
1. Install CAMPAREE + BEERS2
2. Configure diploid simulation with indels
3. Generate reads with known allelic ratios
4. Run WASP2 pipeline
5. Validate

**Time investment**: 2-3 days setup + learning

**Benefit**: Can cite established simulator

**Downside**: Overkill, harder to explain

---

### **Option 3: Hybrid Approach** ⭐⭐ **GOOD COMPROMISE**

**Use both**:
1. **Simple custom simulation** (main validation)
   - Quick
   - Tailored
   - Easy to understand

2. **CAMPAREE** (supplementary validation)
   - If time permits
   - Addresses "is this realistic?" concern
   - Can add to supplement

**Manuscript**:
> "We validated using custom simulation (main text) and confirmed results with CAMPAREE (Fig. S1D)."

---

## What Published Papers Actually Do

### **Survey of ASE validation approaches**:

| Tool | Validation Approach | Simulator Used |
|------|-------------------|----------------|
| WASP (2015) | Simple read generation | **Custom** |
| ASElux (2018) | Simulated data | Not specified |
| MBASED (2014) | Imprinted genes | **None** (real data) |
| MixALime (2024) | Count generation | **Custom** (distributions) |
| ASEReadCounter | GATK benchmarks | Various |

**Pattern**: Most tools use **custom validation**, not off-the-shelf simulators!

---

## Final Recommendation

### **Use Your Custom Simulator** (`simulate_indel_ase.py`)

**Reasons**:
1. ✅ **Precedent**: WASP original paper did this
2. ✅ **Simplicity**: Easy to understand and explain
3. ✅ **Appropriate**: Tests bias removal, not biology
4. ✅ **Fast**: Already implemented
5. ✅ **Targeted**: Specifically tests indel handling

**For manuscript**:
> "Following the approach of van de Geijn et al. (WASP, 2015), we generated synthetic paired-end reads (150bp) with known allelic ratios at heterozygous SNPs, insertions, and deletions. We simulated reads from both haplotypes with controlled ratios (1:1, 2:1, 4:1) and processed them through the WASP2 pipeline to measure recovery accuracy."

**If reviewers ask**:
> "We followed the validation strategy established by the original WASP publication (van de Geijn et al., 2015), which used targeted read simulation rather than complex genome-wide simulators. This approach directly validates bias removal at variant sites, which is the core function of WASP2."

---

## If You Want to Use CAMPAREE (Optional)

**Installation**:
```bash
# Clone CAMPAREE
git clone https://github.com/kaitlinchaung/camparee
pip install -e camparee

# Will also need BEERS2
git clone https://github.com/itmat/BEERS2
```

**Usage** (simplified):
```python
# Configure diploid genome with indels
camparee configure \
    --ref genome.fa \
    --variants your_indels.vcf \
    --output config.yaml

# Generate reads with ASE
camparee simulate \
    --config config.yaml \
    --ase-ratio 2.0 \
    --output simulated_reads.fq
```

**Time**: 1-2 days to learn + set up

**Benefit**: Can cite established tool in supplement

---

## Bottom Line

### **Don't reinvent the wheel, but...**

✅ **Your custom simulator IS appropriate** because:
- WASP original paper did the same
- You're testing bias removal, not modeling biology
- Targeted validation is better than generic simulation
- Easier to understand and reproduce

❌ **Don't use CAMPAREE/BEERS unless**:
- You have extra time (2-3 days)
- Reviewers specifically request "established simulator"
- You want supplementary validation

### **Recommended Strategy**:

**Main validation** (2 days):
1. ✅ Custom simulation (`simulate_indel_ase.py`)
2. ✅ Imprinted genes (real data)
3. ✅ SNP/indel consistency

**If time/reviewer demands** (additional 2-3 days):
4. CAMPAREE validation (supplement)

**This gives you**:
- Fast path to publication (option 1-3)
- Fallback if reviewers push (option 4)
- Established precedent (WASP 2015)
- Clear, understandable validation

---

## Sources:

- [WASP: Nature Methods 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4626402/) - Original WASP validation
- [CAMPAREE: BMC Genomics 2021](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-021-07934-2) - ASE-capable simulator
- [BEERS2: Briefings in Bioinformatics 2024](https://academic.oup.com/bib/article/25/3/bbae164/7644138) - RNA-seq simulator
- [Polyester: Bioinformatics 2015](https://academic.oup.com/bioinformatics/article/31/17/2778/183245) - RNA-seq simulator (no ASE)
- [Tools and Best Practices: Genome Biology 2015](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0762-6) - Indel challenges
- [ASElux: Bioinformatics 2018](https://academic.oup.com/bioinformatics/article/34/8/1313/4653700) - ASE detection tool
- [MBASED: Genome Biology 2014](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0405-3) - ASE meta-analysis

**Your approach is sound and has precedent!** ✅
