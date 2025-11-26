# Parallel Validation Plan - WASP2 Indel Support

**Date**: 2025-11-25
**Current Branch**: `rust-optimization-plink2`

---

## 🎯 **Overall Goal**

Validate WASP2 indel support through:
1. **Simulation** (gold standard) - Proves correctness
2. **GM12878 Platinum Genomes** (real benchmark) - Proves real-world performance
3. **Statistical rigor** - Publication-ready analysis

---

## 🔀 **Parallel Workstreams (Git Branches)**

### **Branch 1: `simulation-statistical-rigor`**
**Owner**: Sub-agent or you
**Time**: 8 hours
**Dependencies**: None (can start immediately)

**Tasks**:
1. Update `simulate_indel_ase_v2.py` to add:
   - Confidence intervals (bootstrapping)
   - Formal hypothesis tests (t-tests, Wilcoxon)
   - Sample size justification (power analysis)
   - Statistical summary reporting

2. Add analysis script: `analyze_simulation_results.py`
   - Reads simulation output
   - Calculates CIs, p-values
   - Generates publication-quality figures
   - Outputs statistical summary table

3. Documentation:
   - Update `SIMULATION_FRAMEWORK_EXPLAINED.md`
   - Add statistical methods section

**Output**:
- Enhanced simulation framework
- Statistical analysis script
- Ready to merge to main

---

### **Branch 2: `run-baseline-simulation`**
**Owner**: Background process (nohup)
**Time**: 30-60 min runtime
**Dependencies**: None (use current code)

**Tasks**:
1. Setup run script:
   ```bash
   #!/bin/bash
   # run_simulation_baseline.sh

   OUTDIR="simulation_results_baseline_$(date +%Y%m%d_%H%M%S)"
   mkdir -p ${OUTDIR}

   nohup ./run_simulation.sh moderate > ${OUTDIR}/simulation.log 2>&1 &
   echo $! > ${OUTDIR}/simulation.pid

   echo "Simulation running in background"
   echo "PID: $(cat ${OUTDIR}/simulation.pid)"
   echo "Log: ${OUTDIR}/simulation.log"
   echo "Monitor: tail -f ${OUTDIR}/simulation.log"
   ```

2. Run it:
   ```bash
   chmod +x run_simulation_baseline.sh
   ./run_simulation_baseline.sh
   ```

3. Monitor progress:
   ```bash
   tail -f simulation_results_baseline_*/simulation.log
   ```

**Output**:
- Baseline simulation results (270 tests)
- Can use while waiting for statistical enhancements
- Proves algorithm works

---

### **Branch 3: `gm12878-benchmark`**
**Owner**: Sub-agent or you
**Time**: 4-6 hours (processing) + overnight (if needed)
**Dependencies**: None (data already found)

**Tasks**:
1. Create benchmark directory structure:
   ```bash
   mkdir -p gm12878_benchmark/{scripts,results,logs}
   ```

2. Write processing script: `gm12878_benchmark/scripts/run_benchmark.sh`
   ```bash
   #!/bin/bash
   set -e

   # Input data (Aaron's)
   BAM="/iblm/netapp/data3/aho/alignment/GM12878_rna_v2/GM12878_merged.sorted.bam"
   VCF="/iblm/netapp/data1/aho/variants/NA12878.vcf.gz"
   GTF="/iblm/netapp/data3/aho/project_data/wasp2/imprinted_rna/data/geneimprint.gtf"

   OUTDIR="gm12878_benchmark/results"
   LOGDIR="gm12878_benchmark/logs"

   # Run with indels
   echo "Running WASP2 counting with indels..."
   python src/counting count-variants \
       ${BAM} ${VCF} \
       -s NA12878 \
       -r ${GTF} \
       -o ${OUTDIR}/gm12878_counts_WITH_INDELS.tsv \
       --gene_feature transcript \
       --gene_attribute transcript_id \
       --gene_parent gene_name \
       --include_indels \
       2>&1 | tee ${LOGDIR}/counting.log

   echo "Running AI analysis..."
   python src/analysis find-imbalance \
       --phased \
       --out ${OUTDIR}/gm12878_ai_WITH_INDELS.tsv \
       --group gene_name \
       ${OUTDIR}/gm12878_counts_WITH_INDELS.tsv \
       2>&1 | tee ${LOGDIR}/analysis.log

   echo "Done!"
   ```

3. Create comparison script: `gm12878_benchmark/scripts/compare_to_aaron.py`
   ```python
   import pandas as pd

   # Load Aaron's SNP-only results
   aaron = pd.read_csv(
       "/iblm/netapp/data3/aho/project_data/wasp2/imprinted_rna/outputs/GM12878_geneimprint_transcript_gene_counts.tsv",
       sep="\t"
   )

   # Load our SNP+indel results
   ours = pd.read_csv(
       "gm12878_benchmark/results/gm12878_counts_WITH_INDELS.tsv",
       sep="\t"
   )

   # Count variants
   print(f"Aaron (SNP-only): {len(aaron):,} variants")
   print(f"Ours (SNP+indel): {len(ours):,} variants")
   print(f"Additional indels: {len(ours) - len(aaron):,}")

   # Check indel proportion
   # TODO: Implement indel detection from variant format

   # Compare gene-level AI
   # TODO: Load AI results and compare
   ```

4. Run with nohup:
   ```bash
   nohup bash gm12878_benchmark/scripts/run_benchmark.sh > gm12878_benchmark/logs/run.log 2>&1 &
   ```

**Output**:
- GM12878 results with indels
- Comparison to Aaron's SNP-only results
- Publication-quality benchmark validation

---

### **Branch 4: `optimization-threading`**
**Owner**: Sub-agent or you (low priority)
**Time**: 2-3 hours
**Dependencies**: None

**Tasks**:
1. Add threading to simulation:
   ```python
   # In simulate_indel_ase_v2.py
   import multiprocessing

   def align_with_bwa(...):
       threads = min(multiprocessing.cpu_count(), 16)
       subprocess.run(['bwa', 'mem', '-t', str(threads), ...])
   ```

2. Test performance:
   - Benchmark: Current vs threaded
   - Verify correctness (same results)

3. Optional: Profile Rust code for bottlenecks
   - Only if we need faster processing

**Output**:
- 2-3x faster simulation
- Faster GM12878 processing

---

### **Branch 5: `rust-indel-optimization`** (Optional)
**Owner**: Sub-agent (if needed)
**Time**: 1-2 days
**Dependencies**: After benchmark runs

**Tasks**:
1. Profile indel processing code
2. Identify bottlenecks
3. Optimize if needed

**Output**:
- Faster indel handling (if needed)

---

## 🚀 **Execution Plan**

### **Phase 1: Immediate (Start Now)**

**In parallel**:
```bash
# Terminal 1: Start baseline simulation
./run_simulation_baseline.sh
# Monitor: tail -f simulation_results_baseline_*/simulation.log

# Terminal 2: Start GM12878 benchmark
cd gm12878_benchmark/scripts
nohup bash run_benchmark.sh > ../logs/run.log 2>&1 &
# Monitor: tail -f ../logs/run.log

# Terminal 3: Work on statistical enhancements
git checkout -b simulation-statistical-rigor
# Edit simulate_indel_ase_v2.py
# Add CIs, hypothesis tests, etc.
```

### **Phase 2: While Things Run (2-4 hours)**

1. **Monitor background jobs**:
   - Check simulation progress
   - Check GM12878 processing

2. **Implement statistical rigor**:
   - Add bootstrap CIs
   - Add hypothesis tests
   - Create analysis script

3. **Setup environment** (if needed):
   - Verify WASP2 dependencies
   - Test `--include_indels` flag works

### **Phase 3: Analysis (After jobs complete)**

1. **Analyze simulation results**:
   - Run statistical analysis
   - Generate figures
   - Check error rates

2. **Analyze GM12878 results**:
   - Compare to Aaron's SNP-only
   - Count indel contribution
   - Check classic genes (H19, IGF2, etc.)

3. **Create comparison report**:
   - Simulation validates correctness
   - GM12878 validates real-world performance
   - Combined = publication-ready

---

## 📊 **Git Branch Strategy**

```
main (or master)
│
├── rust-optimization-plink2 (current)
│   └── simulation-statistical-rigor (new)
│       └── [Add CIs, tests, analysis]
│
├── run-baseline-simulation (tag, not branch)
│   └── [Just run current code]
│
├── gm12878-benchmark (new)
│   └── [Setup and run GM12878 validation]
│
└── optimization-threading (new, low priority)
    └── [Add threading, profiling]
```

**Merge strategy**:
1. `simulation-statistical-rigor` → `rust-optimization-plink2` (when done)
2. `gm12878-benchmark` → `rust-optimization-plink2` (when done)
3. `optimization-threading` → main (optional)

---

## 🤖 **Sub-Agent Assignments**

### **Agent 1: Statistical Enhancement**
```
Task: Enhance simulation with statistical rigor
Branch: simulation-statistical-rigor
Time: 8 hours
Files to modify:
- simulate_indel_ase_v2.py
- New: analyze_simulation_results.py
- Update: SIMULATION_FRAMEWORK_EXPLAINED.md
```

### **Agent 2: GM12878 Benchmark Setup**
```
Task: Create and run GM12878 benchmark pipeline
Branch: gm12878-benchmark
Time: 4-6 hours + processing time
Files to create:
- gm12878_benchmark/scripts/run_benchmark.sh
- gm12878_benchmark/scripts/compare_to_aaron.py
- gm12878_benchmark/scripts/analyze_results.py
```

### **Agent 3: Threading Optimization** (Optional)
```
Task: Add threading to speed up simulation
Branch: optimization-threading
Time: 2-3 hours
Files to modify:
- simulate_indel_ase_v2.py (add threading)
- Benchmark performance
```

---

## ✅ **Success Criteria**

### **Simulation**:
- [ ] 270 tests complete
- [ ] Mean error <5%
- [ ] All tests <10% error
- [ ] Confidence intervals calculated
- [ ] Statistical significance tests pass
- [ ] Publication-quality figures generated

### **GM12878 Benchmark**:
- [ ] Processing complete (SNPs + indels)
- [ ] Variant counts: >118K (Aaron's SNP-only baseline)
- [ ] Indel count identified
- [ ] Classic genes show expected AI (H19, IGF2, SNRPN, XIST)
- [ ] Comparison report generated
- [ ] Publication text drafted

### **Combined Validation**:
- [ ] Simulation proves correctness
- [ ] GM12878 proves real-world performance
- [ ] Comparison shows value of indel support
- [ ] Manuscript-ready figures and tables

---

## 🕐 **Timeline**

### **Day 1** (Today):
- Hour 1: Start simulation (background)
- Hour 1: Start GM12878 processing (background)
- Hours 2-8: Implement statistical rigor
- Evening: Check job progress

### **Day 2**:
- Morning: Analyze simulation results
- Afternoon: Analyze GM12878 results
- Evening: Create comparison report

### **Day 3** (Optional):
- Threading optimization
- Final polishing
- Manuscript writing

**Total**: 2-3 days for complete, publication-ready validation

---

## 📝 **Next Immediate Steps**

1. **Create branches**:
   ```bash
   git checkout -b simulation-statistical-rigor
   git checkout -b gm12878-benchmark
   git checkout -b optimization-threading
   git checkout rust-optimization-plink2  # back to main work branch
   ```

2. **Setup run scripts**:
   - Create `run_simulation_baseline.sh`
   - Create `gm12878_benchmark/scripts/run_benchmark.sh`

3. **Start parallel work**:
   - Launch simulation (nohup)
   - Launch GM12878 processing (nohup)
   - Work on statistical enhancements (interactive)

---

**Ready to start? Which workstream do you want to tackle first?**

Options:
1. Start both background jobs running (simulation + GM12878)
2. Create statistical enhancement branch
3. Setup threading optimization
4. All of the above in parallel using sub-agents
