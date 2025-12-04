# WASP2 Simulation Benchmark Implementation Plan

**Created:** 2025-12-03
**Goal:** Publication-ready simulation benchmarks for Nature Methods / Bioinformatics
**Current State:** 30 tests (minimum tier), single-end, synthetic genome

---

## Executive Summary

| Gap | Priority | Effort | Impact |
|-----|----------|--------|--------|
| Paired-end simulation | P0 | 2-3 hrs | Critical - matches real data |
| Run comprehensive tier (810 tests) | P0 | 2 hrs runtime | Required baseline |
| GATK ASEReadCounter comparison | P0 | 2-3 hrs | Required for publication |
| Publication metrics (correlation, bias) | P1 | 2 hrs | Strengthens claims |
| Real genome (chr22) | P1 | 3-4 hrs | Addresses reviewer concern |
| Multi-variant reads | P2 | 2-3 hrs | Edge case coverage |
| Complex indel scenarios | P2 | 2-3 hrs | Validates indel claims |

---

## Work Streams (Parallelizable)

### Stream A: Core Simulation Improvements
**Can run independently**

#### A1. Paired-End Simulation Module
**File:** `simulate_paired_end_ase.py`

```python
# Key changes from v2:
# 1. Generate R1 + R2 FASTQ files
# 2. Proper insert size distribution (mean=300, std=50)
# 3. Variant can be in R1, R2, or both
# 4. Output: sample_R1.fq.gz, sample_R2.fq.gz
```

**Test matrix:**
- Insert sizes: 200, 300, 400 bp
- Variant positions: R1-only, R2-only, spanning both
- Read lengths: 75, 100, 150 bp

#### A2. Comprehensive Tier Execution
**File:** `run_simulation.sh`

```bash
# Just run existing framework at higher tier
python simulate_indel_ase_v2.py --tier comprehensive --keep
# Runtime: ~2 hours
# Output: 810 test results
```

#### A3. Real Genome Integration
**File:** `simulate_real_genome_ase.py`

```python
# Use chr22 from GRCh38
# Extract real het sites from 1000 Genomes
# Simulate reads at those positions
# More realistic: repeats, GC content, real variant context
```

---

### Stream B: Competitor Comparison
**Can run independently after Stream A produces data**

#### B1. GATK ASEReadCounter Integration
**File:** `benchmark_vs_gatk.py`

```python
def run_gatk_ase(bam_file, vcf_file, ref_fasta, output_table):
    """Run GATK ASEReadCounter on same data."""
    cmd = [
        "gatk", "ASEReadCounter",
        "-R", ref_fasta,
        "-I", bam_file,
        "-V", vcf_file,
        "-O", output_table,
        "--min-mapping-quality", "10",
        "--min-base-quality", "20"
    ]
    subprocess.run(cmd, check=True)

def compare_counts(wasp2_counts, gatk_counts, ground_truth):
    """Compare WASP2 vs GATK vs ground truth."""
    # Metrics: correlation, RMSE, bias
    pass
```

#### B2. Comparison Metrics
**Output:** `benchmark_comparison.csv`

| Metric | WASP2 | GATK | Description |
|--------|-------|------|-------------|
| Pearson r | ? | ? | Correlation with ground truth |
| RMSE | ? | ? | Root mean squared error |
| REF bias | ? | ? | (REF-ALT)/(REF+ALT) deviation from 0 |
| Runtime (s) | ? | ? | Processing time |

---

### Stream C: Publication Metrics & Figures
**Can run after A and B complete**

#### C1. Enhanced Metrics Module
**File:** `simulation_metrics.py`

```python
def compute_publication_metrics(results_df, ground_truth_df):
    """Compute all metrics needed for publication."""

    metrics = {
        # Accuracy metrics
        "pearson_r": pearsonr(observed, expected),
        "spearman_rho": spearmanr(observed, expected),
        "rmse": np.sqrt(mean_squared_error(observed, expected)),
        "mae": mean_absolute_error(observed, expected),

        # Bias metrics
        "ref_bias": np.mean(ref_counts / (ref_counts + alt_counts)) - 0.5,
        "ref_bias_std": np.std(ref_counts / (ref_counts + alt_counts)),

        # Detection metrics (for ASE calling)
        "sensitivity": tp / (tp + fn),
        "specificity": tn / (tn + fp),
        "precision": tp / (tp + fp),
        "f1_score": 2 * (precision * recall) / (precision + recall),

        # By variant type
        "snp_accuracy": ...,
        "ins_accuracy": ...,
        "del_accuracy": ...,
    }

    return metrics
```

#### C2. Publication Figures
**File:** `generate_figures.py`

```python
# Figure 1: Ground truth correlation
# - Scatter: observed vs expected allelic ratio
# - Stratified by variant type (SNP, INS, DEL)
# - WASP2 vs GATK comparison

# Figure 2: REF/ALT bias
# - Before WASP filtering vs after
# - By coverage level

# Figure 3: Performance by variant type
# - Box plots: error distribution
# - SNP vs INS vs DEL

# Figure 4: Edge cases
# - Large indels
# - Multi-variant reads
# - Low coverage scenarios
```

---

### Stream D: Edge Cases & Stress Tests
**Can run in parallel with other streams**

#### D1. Multi-Variant Scenarios
**File:** `simulate_multi_variant.py`

```python
# Scenarios where reads span multiple variants
multi_variant_configs = [
    # Two SNPs close together
    {"name": "snp_snp_50bp", "variants": [
        {"type": "SNP", "pos": 1000},
        {"type": "SNP", "pos": 1050}
    ]},

    # SNP + indel
    {"name": "snp_ins_75bp", "variants": [
        {"type": "SNP", "pos": 1000},
        {"type": "INS", "pos": 1075}
    ]},

    # Complex: 3 variants in one read
    {"name": "triple_variant", "variants": [
        {"type": "SNP", "pos": 1000},
        {"type": "DEL", "pos": 1050},
        {"type": "SNP", "pos": 1100}
    ]},
]
```

#### D2. Stress Test Scenarios
**File:** `simulate_edge_cases.py`

| Scenario | Description | Why Important |
|----------|-------------|---------------|
| Large DEL (50bp) | Tests CIGAR parsing | WASP2 indel claim |
| Large INS (30bp) | Tests sequence handling | WASP2 indel claim |
| Adjacent indels | DEL then INS within 10bp | Alignment edge case |
| Homopolymer indel | A→AA in AAAA region | Known aligner weakness |
| Very low coverage (5x) | Statistical edge case | Real-world scenario |
| Very high coverage (500x) | Computational stress | Scalability |
| Near read boundary | Variant at pos 1 or 149 | Edge of read |

---

## Branching Strategy

```
master
  │
  └── ropc-indels (current)
        │
        ├── sim/paired-end      # Stream A1: Paired-end simulation
        │
        ├── sim/comprehensive   # Stream A2: Run comprehensive tier
        │
        ├── sim/gatk-compare    # Stream B: GATK comparison
        │
        ├── sim/metrics         # Stream C: Publication metrics
        │
        └── sim/edge-cases      # Stream D: Edge cases
```

**Merge order:**
1. `sim/comprehensive` → `ropc-indels` (baseline results)
2. `sim/paired-end` → `ropc-indels` (core improvement)
3. `sim/gatk-compare` → `ropc-indels` (competitor comparison)
4. `sim/metrics` → `ropc-indels` (publication polish)
5. `sim/edge-cases` → `ropc-indels` (final strengthening)

---

## Parallelization Opportunities

### Independent Tasks (can run simultaneously)

| Task | Resources | Time | Dependencies |
|------|-----------|------|--------------|
| Run comprehensive tier | 1 job, 8 cores | 2 hrs | None |
| Download chr22 + prep | 1 job | 30 min | None |
| Write paired-end module | Dev time | 2-3 hrs | None |
| Write GATK comparison | Dev time | 2 hrs | None |
| Write metrics module | Dev time | 2 hrs | None |

### Sequential Tasks (must wait)

| Task | Depends On | Time |
|------|------------|------|
| Run paired-end sim | paired-end module | 2 hrs |
| Run GATK on sim data | comprehensive tier done | 1 hr |
| Generate figures | all sim data ready | 1 hr |
| Run edge cases | paired-end module | 2 hrs |

---

## Immediate Actions (Today)

### Can Start Now (Parallel)

1. **Job 1:** Run comprehensive tier on existing simulation
   ```bash
   qsub -N sim_comprehensive run_simulation.sh comprehensive
   ```

2. **Job 2:** Download chr22 reference
   ```bash
   wget -P data/ https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz
   ```

3. **Dev work:** Start writing paired-end simulation module

### Dependencies to Check

- [ ] GATK installed? `which gatk`
- [ ] BWA index for chr22 needed
- [ ] Sufficient disk space for simulation outputs

---

## Success Criteria

### Minimum for Bioinformatics (Oxford)
- [ ] 810+ tests (comprehensive tier)
- [ ] Paired-end simulation
- [ ] GATK comparison showing concordance
- [ ] Correlation r > 0.95 with ground truth
- [ ] REF bias < 2% after WASP filtering

### Target for Nature Methods
- [ ] All above, plus:
- [ ] Real genome (chr22) simulation
- [ ] Multi-variant read scenarios
- [ ] Edge case coverage
- [ ] Publication-quality figures
- [ ] Biological validation (imprinted genes)

---

## File Structure After Implementation

```
WASP2-exp/
├── simulation/
│   ├── simulate_paired_end_ase.py      # NEW: Paired-end simulation
│   ├── simulate_real_genome_ase.py     # NEW: Chr22-based simulation
│   ├── simulate_multi_variant.py       # NEW: Multi-variant scenarios
│   ├── simulate_edge_cases.py          # NEW: Stress tests
│   ├── benchmark_vs_gatk.py            # NEW: GATK comparison
│   ├── simulation_metrics.py           # NEW: Publication metrics
│   ├── generate_figures.py             # NEW: Figure generation
│   └── run_all_simulations.sh          # NEW: Master runner
├── simulate_indel_ase_v2.py            # EXISTING: Keep as baseline
├── results/
│   ├── comprehensive_YYYYMMDD/         # Full tier results
│   ├── paired_end_YYYYMMDD/            # Paired-end results
│   ├── gatk_comparison_YYYYMMDD/       # GATK comparison
│   └── figures/                        # Publication figures
└── SIMULATION_BENCHMARK_PLAN.md        # This document
```

---

## Estimated Timeline

| Day | Tasks | Deliverable |
|-----|-------|-------------|
| Day 1 | Run comprehensive, download chr22, write paired-end module | Baseline results |
| Day 2 | Run paired-end sim, write GATK comparison | Paired-end results |
| Day 3 | Run GATK comparison, write metrics | Comparison data |
| Day 4 | Edge cases, generate figures | Complete benchmark |
| Day 5 | Analysis, documentation | Publication-ready |

---

## Questions to Resolve

1. **Branching:** One branch per stream, or feature branches within `ropc-indels`?
2. **Priority:** GATK comparison vs paired-end first?
3. **Real genome:** chr22 or use existing GM12878 variants?
4. **Parallelism:** How many SGE jobs can we run simultaneously?
5. **Storage:** Where to store large simulation outputs?
