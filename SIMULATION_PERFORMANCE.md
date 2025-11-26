# Simulation Framework Performance Analysis

## Realistic Runtime Estimates

### **Actual Bottlenecks**

The simulation has these major steps:

| Step | Runtime | Bottleneck? |
|------|---------|-------------|
| 1. Create reference genome | ~2 sec | ❌ No |
| 2. Create VCF | ~1 sec | ❌ No |
| 3. Generate FASTQ | ~5-10 sec | ❌ No |
| 4. **BWA alignment** | **~3-8 min** | ✅ **YES** |
| 5. BAM sort/index | ~1-2 min | ⚠️ Minor |
| 6. WASP2 find_intersecting_snps | ~2-5 min | ⚠️ Minor |
| 7. **BWA re-alignment** | **~3-8 min** | ✅ **YES** |
| 8. BAM sort/index (again) | ~1-2 min | ⚠️ Minor |
| 9. WASP2 filtering | ~1-2 min | ❌ No |
| 10. Count alleles | ~30 sec | ❌ No |

**Total bottleneck time**: ~70% spent in BWA alignment (steps 4 + 7)

---

## Tier-by-Tier Performance

### **Minimum Tier** (90 tests)

**Configuration**:
- 9 unique variant positions (3 types × 3 ratios)
- 50x coverage per variant
- 10 replicates (same variants, different seeds)
- **Total unique reads**: ~450 reads
- **Total reads with replicates**: ~4,500 reads (but same variants!)

**Runtime breakdown**:
```
Reference creation:           2 sec
VCF creation:                 1 sec
FASTQ generation:            10 sec  (4,500 reads)
BWA alignment:             3-5 min   ← BOTTLENECK
BAM processing:              1 min
WASP2 processing:          1-2 min
BWA re-alignment:          2-3 min   ← BOTTLENECK
Final processing:            1 min
Counting:                   30 sec

Total: ~8-12 minutes
```

**Actual estimate**: **~10 minutes** (not the 10 min I claimed, close!)

---

### **Moderate Tier** (270 tests)

**Configuration**:
- 27 unique variant positions (3 types × 3 ratios × 3 coverage)
- Average 50x coverage
- **Total unique reads**: ~1,350 reads
- **Total reads with replicates**: ~13,500 reads

**Runtime breakdown**:
```
Reference creation:           2 sec
VCF creation:                 1 sec
FASTQ generation:            30 sec  (13,500 reads)
BWA alignment:             5-10 min  ← BOTTLENECK
BAM processing:            2-3 min
WASP2 processing:          3-5 min
BWA re-alignment:          5-10 min  ← BOTTLENECK
Final processing:          2-3 min
Counting:                    1 min

Total: ~20-35 minutes
```

**Actual estimate**: **~25-30 minutes** (pretty close to 30 min claim!)

---

### **Comprehensive Tier** (810 tests)

**Configuration**:
- 81 unique variant positions
- **Total unique reads**: ~4,000 reads
- **Total reads with replicates**: ~40,000 reads

**Runtime breakdown**:
```
Reference creation:            2 sec
VCF creation:                  1 sec
FASTQ generation:           1-2 min  (40,000 reads)
BWA alignment:            15-25 min  ← BOTTLENECK
BAM processing:            5-8 min
WASP2 processing:          8-12 min
BWA re-alignment:         15-25 min  ← BOTTLENECK
Final processing:          5-8 min
Counting:                   2-3 min

Total: ~50-85 minutes
```

**Actual estimate**: **~1-1.5 hours** (not 2 hours - I overestimated!)

---

## Why Is It This Slow?

### **1. BWA Alignment (70% of runtime)**

BWA has to:
- Build suffix array index (first time only)
- Align each read to reference
- Handle paired-end logic
- Generate SAM output

**For 13,500 reads** (moderate tier):
- BWA throughput: ~500-1000 reads/sec (single-threaded)
- With 4 threads: ~2000-4000 reads/sec
- Time: 13,500 / 3000 = **~4-5 minutes**
- **But we do this TWICE** (initial + remap) = **~10 minutes total**

### **2. BAM Sorting/Indexing (15% of runtime)**

After alignment, we must:
- Convert SAM → BAM
- Sort BAM by coordinate
- Index BAM

**For 13,500 reads**:
- Sort: ~1-2 min
- Index: ~30 sec
- **Done 3-4 times** (initial, WASP2 output, remapped) = **~5-7 min total**

### **3. WASP2 Processing (10% of runtime)**

Python code running:
- Read BAM
- Find variants
- Build position maps (CIGAR parsing)
- Generate quality scores
- Create remapped reads

**For 13,500 reads**:
- Python throughput: ~1000-2000 reads/sec
- Time: ~5-10 min

### **4. Everything Else (5% of runtime)**

- FASTQ generation: Fast (pure Python writes)
- VCF creation: Fast (small file)
- Reference creation: Fast (random sequence)
- Counting: Fast (simple BAM iteration)

---

## Optimizations to Speed It Up

### **1. Increase BWA Threads** (Easy - 2-3x speedup!)

**Current**:
```python
subprocess.run(['bwa', 'mem', '-t', '4', ...])  # 4 threads
```

**Optimized**:
```python
import multiprocessing
threads = multiprocessing.cpu_count()  # Use all CPUs
subprocess.run(['bwa', 'mem', '-t', str(threads), ...])
```

**Impact**:
- 4 threads → 16 threads on modern server
- BWA alignment: 10 min → **3 min** (3x faster)
- **Total moderate tier**: 30 min → **15 min**

---

### **2. Use minimap2 Instead of BWA** (Easy - 3-5x speedup!)

minimap2 is MUCH faster for simulated data:

**Current**:
```bash
bwa mem -t 4 reference.fa reads.fq
# Throughput: ~3,000 reads/sec
```

**Optimized**:
```bash
minimap2 -ax sr reference.fa reads.fq
# Throughput: ~10,000-15,000 reads/sec
```

**Impact**:
- BWA: 10 min
- minimap2: **2-3 min** (4x faster)
- **Total moderate tier**: 30 min → **12-15 min**

---

### **3. Skip Intermediate BAM Sorts** (Moderate - 20% speedup)

We sort BAM multiple times. We could:
- Keep SAM in memory
- Only sort final output
- Use `samtools sort -@` with threads

**Impact**:
- BAM processing: 7 min → **4 min**
- **Total savings**: ~3-5 min

---

### **4. Parallelize Across Variants** (Hard - 5-10x speedup!)

**Current approach**:
```
Generate all reads → Align all → Process all → Done
```

**Parallelized approach**:
```
Variant 1: Generate → Align → Process ┐
Variant 2: Generate → Align → Process ├─ Run in parallel
Variant 3: Generate → Align → Process ┘
```

**Implementation**:
```python
from multiprocessing import Pool

def process_variant_batch(variants):
    # Generate, align, process one batch
    ...

# Split variants into batches
batches = [variants[i:i+10] for i in range(0, len(variants), 10)]

# Process in parallel
with Pool(processes=8) as pool:
    results = pool.map(process_variant_batch, batches)
```

**Impact**:
- With 8 cores: **5-10x speedup**
- **Moderate tier**: 30 min → **3-5 min**
- **Comprehensive tier**: 1.5 hrs → **10-15 min**

---

### **5. Use In-Memory Processing** (Hard - 30% speedup)

Avoid writing intermediate files:
- Keep FASTQs in memory (StringIO)
- Use SAM pipes instead of files
- Only write final BAM

**Impact**: ~20-30% faster

---

## Recommended Optimizations

### **Quick Wins** (Easy to implement):

1. **Use all CPU threads**:
```python
# In align_with_bwa():
threads = min(multiprocessing.cpu_count(), 16)  # Cap at 16
subprocess.run(['bwa', 'mem', '-t', str(threads), ...])
```

2. **Add progress bars**:
```python
from tqdm import tqdm
for gt in tqdm(ground_truth, desc="Processing variants"):
    ...
```

**Implementation time**: 10 minutes
**Speedup**: 2-3x
**New runtimes**:
- Minimum: 10 min → **4-5 min**
- Moderate: 30 min → **12-15 min**
- Comprehensive: 90 min → **30-40 min**

---

### **If You Need It Faster** (More work):

3. **Switch to minimap2**:
```python
subprocess.run(['minimap2', '-ax', 'sr', '-t', str(threads), ref, fastq])
```

**Implementation time**: 30 minutes (test and validate)
**Speedup**: 4-5x
**New runtimes**:
- Minimum: 10 min → **2-3 min**
- Moderate: 30 min → **6-8 min**
- Comprehensive: 90 min → **15-20 min**

---

### **If You're Desperate** (Complex):

4. **Parallelize across variants**:
- Use multiprocessing
- Split into batches
- Combine results

**Implementation time**: 2-3 hours
**Speedup**: 8-10x
**New runtimes**:
- Minimum: 10 min → **1-2 min**
- Moderate: 30 min → **3-5 min**
- Comprehensive: 90 min → **8-10 min**

---

## Current vs Optimized Performance

| Tier | Current | Quick Wins | minimap2 | Parallel |
|------|---------|------------|----------|----------|
| **Minimum (90)** | 10 min | 4 min | 2 min | 1 min |
| **Moderate (270)** | 30 min | 12 min | 6 min | 3 min |
| **Comprehensive (810)** | 90 min | 35 min | 18 min | 8 min |

---

## Is Current Speed Acceptable?

### **For Development/Testing**: ✅ YES
- 10 min for minimum tier is totally fine
- Can run while getting coffee

### **For Publication**: ✅ YES
- 30 min for moderate tier is acceptable
- Run overnight if needed
- Only run once for manuscript

### **For Iteration/Debugging**: ⚠️ MAYBE
- If you need to run many times during development
- Consider quick wins (use all threads)

---

## Bottom Line

**Current performance**:
- ✅ Minimum: ~10 min (acceptable)
- ✅ Moderate: ~30 min (acceptable for publication)
- ⚠️ Comprehensive: ~90 min (bit slow, but rarely needed)

**Main bottleneck**: BWA alignment (~70% of time)

**Easy speedup**: Use all CPU threads (2-3x faster)

**If you need faster**: Switch to minimap2 (4-5x faster)

**My recommendation**:
- Start with current implementation
- If 30 min is too slow, add threading (10 min fix)
- If still too slow, try minimap2 (30 min fix)

**For publication**: Current speed is TOTALLY FINE. 30 minutes is nothing compared to weeks of waiting for reviews! 😄

---

## Want Me to Add Threading Now?

I can add the threading optimization in ~5 minutes if you want:

```python
# Just change this line in align_with_bwa():
threads = 4  # Current

# To this:
import multiprocessing
threads = min(multiprocessing.cpu_count(), 16)  # Use all CPUs (cap at 16)
```

**Would cut runtime in half**: 30 min → 15 min

Let me know!
