# Agent 1: Run Comprehensive Simulation (810 Tests)

## Mission
Run the existing simulation framework at comprehensive tier to establish baseline results. **No code changes required** - execution and validation only.

---

## Repository Context

**GitHub:** https://github.com/Jaureguy760/WASP2-exp.git
**Branch:** `sim/comprehensive`
**Parent Branch:** `ropc-indels`

**Working Directory:**
```
/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp
```

**Conda Environment:** `WASP2_dev2`

---

## Directory Structure (Relevant Files)

```
WASP2-exp/
├── simulate_indel_ase_v2.py      # ← MAIN SIMULATION SCRIPT
├── run_simulation.sh              # ← WRAPPER SCRIPT
├── simulation_results/            # ← OUTPUT DIRECTORY (create if missing)
├── src/
│   └── mapping/
│       └── run_mapping.py         # WASP2 pipeline entry point
└── rust/
    └── src/
        └── unified_pipeline.rs    # Rust acceleration
```

---

## Key File: simulate_indel_ase_v2.py

### What it does:
1. Creates synthetic reference genome (1Mb random sequence)
2. Generates variants (SNPs, insertions, deletions)
3. Creates synthetic reads with known allelic ratios
4. Aligns with BWA (creates realistic CIGARs)
5. Runs WASP2 pipeline
6. Counts alleles and validates against ground truth

### Tiers:
| Tier | Tests | Variants | Coverage Levels | Runtime |
|------|-------|----------|-----------------|---------|
| minimum | 90 | 9 | 1 (50x) | ~10 min |
| moderate | 270 | 9 | 3 (20x, 50x, 100x) | ~30 min |
| **comprehensive** | **810** | **27** | **3 (20x, 50x, 100x)** | **~2 hrs** |

### Variant Types Tested:
- **SNP:** A→G, T→C, G→A
- **INS:** C→CAT (3bp), A→AGGG (4bp), T→TTTAA (5bp), up to 10bp
- **DEL:** GCC→G (2bp), ATATA→A (4bp), up to 19bp

### Allelic Ratios Tested:
- 1.0 (balanced)
- 2.0 (2:1 REF:ALT)
- 4.0 (4:1 REF:ALT)

---

## Step-by-Step Execution

### Step 1: Environment Setup

```bash
# SSH to cluster
ssh jjaureguy@iblm-cluster.salk.edu

# Activate conda
source /iblm/netapp/home/jjaureguy/mambaforge/etc/profile.d/conda.sh
conda activate WASP2_dev2

# Navigate to repo
cd /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp

# Checkout branch
git checkout sim/comprehensive
git pull origin sim/comprehensive
```

### Step 2: Verify Dependencies

```bash
# Check Rust module
python -c "import wasp2_rust; print('wasp2_rust: OK')"

# Check BWA
which bwa && bwa 2>&1 | head -3

# Check samtools
which samtools && samtools --version | head -1

# Check disk space (need ~5GB for comprehensive)
df -h .
```

### Step 3: Create Output Directory

```bash
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUTDIR="simulation_results/comprehensive_${TIMESTAMP}"
mkdir -p ${OUTDIR}
mkdir -p simulation_results  # Ensure parent exists
echo "Output: ${OUTDIR}"
```

### Step 4: Create SGE Job Script

```bash
cat > sge_comprehensive_sim.sh << 'JOBSCRIPT'
#!/bin/bash
#$ -N wasp2_sim_comp
#$ -V
#$ -pe iblm 8
#$ -l h_vmem=16G
#$ -j y
#$ -o simulation_results/
#$ -cwd

set -e

source /iblm/netapp/home/jjaureguy/mambaforge/etc/profile.d/conda.sh
conda activate WASP2_dev2

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUTDIR="simulation_results/comprehensive_${TIMESTAMP}"
mkdir -p ${OUTDIR}

echo "Starting comprehensive simulation at $(date)"
echo "Output directory: ${OUTDIR}"

python simulate_indel_ase_v2.py \
    --tier comprehensive \
    --workdir ${OUTDIR} \
    --keep

echo "Completed at $(date)"
echo "Results: ${OUTDIR}/simulation_results.csv"
JOBSCRIPT
```

### Step 5: Submit Job

```bash
qsub sge_comprehensive_sim.sh
```

### Step 6: Monitor Progress

```bash
# Check job status
qstat -u jjaureguy

# Watch output log (job ID will vary)
tail -f simulation_results/wasp2_sim_comp.o*

# Expected progress output:
# - "Creating reference genome..."
# - "Generating synthetic FASTQ..."
# - "Aligning reads with BWA..."
# - "RUNNING WASP2 PIPELINE"
# - "Counting alleles..."
# - "VALIDATION RESULTS"
```

---

## Expected Outputs

### Files Created:
```
${OUTDIR}/
├── simulation_results.csv    # ← MAIN RESULTS (validate this)
├── reference.fa              # Synthetic reference
├── reference.fa.fai          # Reference index
├── variants.vcf.gz           # VCF with test variants
├── variants.vcf.gz.tbi       # VCF index
├── synthetic.fq              # Generated reads
├── aligned.sorted.bam        # BWA alignments
├── aligned.sorted.bam.bai    # BAM index
└── wasp2_output/             # WASP2 pipeline output
    ├── keep.bam
    ├── to.remap.bam
    └── remap.fq.gz
```

### simulation_results.csv Format:
```csv
chrom,pos,variant_type,coverage,replicate,true_ratio,ref_count,alt_count,total_reads,observed_ratio,error,error_pct,status
chr1,50000,SNP,50,0,1.0,25,24,49,1.042,0.042,4.2,PASS
chr1,50000,SNP,50,1,1.0,24,26,50,0.923,0.077,7.7,PASS
...
```

---

## Validation Criteria

### Required Checks:

```bash
# 1. Check row count (should be 810)
wc -l ${OUTDIR}/simulation_results.csv
# Expected: 811 (810 + header)

# 2. Check pass rate
python -c "
import pandas as pd
df = pd.read_csv('${OUTDIR}/simulation_results.csv')
print(f'Total tests: {len(df)}')
print(f'Pass rate: {(df.status == \"PASS\").mean()*100:.1f}%')
print()
print('By variant type:')
for vtype in ['SNP', 'INS', 'DEL']:
    sub = df[df.variant_type == vtype]
    print(f'  {vtype}: {len(sub)} tests, {(sub.status==\"PASS\").mean()*100:.0f}% pass')
"
```

### Success Criteria:
- [ ] 810 tests completed
- [ ] Overall pass rate > 90%
- [ ] SNP pass rate > 95%
- [ ] INS pass rate > 90%
- [ ] DEL pass rate > 85%
- [ ] No ERROR or CRASH in log file

---

## Commit and Push Results

```bash
# Only commit the results CSV (not large BAM files)
git add simulation_results/comprehensive_*/simulation_results.csv
git add sge_comprehensive_sim.sh

git commit -m "data: comprehensive simulation results (810 tests)

Tier: comprehensive
Tests: 810
Pass rate: XX%

By variant type:
- SNP: XX% pass
- INS: XX% pass
- DEL: XX% pass

Runtime: ~X hours
"

git push origin sim/comprehensive
```

---

## Troubleshooting

### Common Issues:

**1. BWA index missing:**
```bash
bwa index ${OUTDIR}/reference.fa
```

**2. Out of memory:**
- Reduce to moderate tier first: `--tier moderate`
- Request more memory: `#$ -l h_vmem=32G`

**3. WASP2 pipeline fails:**
```bash
# Check Rust module
python -c "from wasp2_rust import unified_make_reads_py; print('OK')"

# Rebuild if needed
cd rust && maturin develop --release && cd ..
```

**4. Job killed by SGE:**
- Check `qstat -j <job_id>` for reason
- Usually memory or time limit

---

## What NOT To Do

- ❌ Do NOT modify simulate_indel_ase_v2.py
- ❌ Do NOT change test parameters
- ❌ Do NOT skip validation
- ❌ Do NOT commit large BAM files to git
- ❌ Do NOT run interactively on login node (use SGE)

---

## Handoff

When complete, notify that:
1. Results CSV is committed to `sim/comprehensive`
2. Pass rate and breakdown documented in commit message
3. Ready for Agent 4 (metrics) to consume
