# Agent 2: Implement Paired-End Simulation

## Mission
Create a paired-end read simulation module. Current simulation generates single-end reads only, but WASP2 is designed for paired-end RNA-seq data. This is a **critical gap** for publication.

---

## Repository Context

**GitHub:** https://github.com/Jaureguy760/WASP2-exp.git
**Branch:** `sim/paired-end`
**Parent Branch:** `ropc-indels`

**Working Directory:**
```
/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp
```

**Conda Environment:** `WASP2_dev2`

---

## Problem Statement

**Current State:**
- `simulate_indel_ase_v2.py` generates SINGLE-END reads only
- Writes one FASTQ file
- Aligns with `bwa mem ref.fa reads.fq` (single-end mode)

**Required State:**
- Generate PAIRED-END reads (R1 and R2)
- Proper insert size distribution (mean=300bp, std=50bp)
- Variant can appear in R1, R2, or spanning both
- Align with `bwa mem ref.fa R1.fq R2.fq` (paired-end mode)
- Run through WASP2 paired-end pipeline

---

## Existing Code Reference

### File: simulate_indel_ase_v2.py (Key Functions to Adapt)

```python
# Current single-end read generation (lines 92-130):
def create_read_sequence(
    ref_seq: str,
    var_pos: int,
    allele: str,
    read_length: int = 150,
    error_rate: float = 0.01
) -> Tuple[str, np.ndarray]:
    """
    Create read sequence with specific allele at variant position.

    Returns:
        Tuple of (sequence, quality_scores)
    """
    # Random position offset so variant isn't always centered
    offset = random.randint(-30, 30)
    read_start = max(0, var_pos - read_length // 2 + offset)

    # Build read sequence
    left_seq = ref_seq[read_start:var_pos]
    right_seq = ref_seq[var_pos + len(allele):read_start + read_length]
    read_seq = left_seq + allele + right_seq

    # Truncate or pad to exact read length
    if len(read_seq) > read_length:
        read_seq = read_seq[:read_length]
    elif len(read_seq) < read_length:
        pad_len = read_length - len(read_seq)
        read_seq += ref_seq[read_start + len(read_seq):read_start + len(read_seq) + pad_len]

    # Add sequencing errors
    read_seq = add_sequencing_errors(read_seq, error_rate)

    # Generate quality scores
    qualities = generate_quality_scores(len(read_seq))

    return read_seq, qualities


# Current FASTQ writing (lines 133-139):
def write_fastq_record(fq_handle, read_id: str, sequence: str, qualities: np.ndarray):
    """Write single FASTQ record."""
    qual_string = ''.join([chr(q + 33) for q in qualities])
    fq_handle.write(f"@{read_id}\n")
    fq_handle.write(f"{sequence}\n")
    fq_handle.write(f"+\n")
    fq_handle.write(f"{qual_string}\n")


# Current BWA alignment (lines 201-253) - SINGLE-END:
def align_with_bwa(ref_fasta, fastq_file, output_bam, threads=4):
    subprocess.run([
        'bwa', 'mem',
        '-t', str(threads),
        ref_fasta,
        fastq_file  # ← SINGLE file
    ], stdout=sam, check=True)
```

### GroundTruth Dataclass (keep this):
```python
@dataclass
class GroundTruth:
    """Known truth for validation."""
    chrom: str
    pos: int
    ref_allele: str
    alt_allele: str
    true_ratio: float       # REF/ALT expression ratio
    variant_type: str       # "SNP", "INS", "DEL"
    coverage: int           # Number of read PAIRS
    replicate: int
    seed: int
```

---

## Implementation Specification

### New File: `simulation/simulate_paired_end_ase.py`

### Required Functions:

#### 1. reverse_complement()
```python
def reverse_complement(seq: str) -> str:
    """Return reverse complement of DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(b, 'N') for b in reversed(seq))
```

#### 2. generate_insert_size()
```python
def generate_insert_size(mean: int = 300, std: int = 50, min_size: int = 150) -> int:
    """Sample insert size from normal distribution.

    Args:
        mean: Mean insert size (default 300bp, typical for Illumina)
        std: Standard deviation (default 50bp)
        min_size: Minimum valid insert size (must be > read_length)

    Returns:
        Insert size in bp
    """
    size = int(np.random.normal(mean, std))
    return max(min_size, size)
```

#### 3. create_paired_reads() - CORE FUNCTION
```python
def create_paired_reads(
    ref_seq: str,
    var_pos: int,
    ref_allele: str,
    alt_allele: str,
    use_alt: bool,
    read_length: int = 150,
    insert_mean: int = 300,
    insert_std: int = 50,
    error_rate: float = 0.01
) -> Tuple[Tuple[str, str], Tuple[str, str]]:
    """
    Generate paired-end reads with variant at specified position.

    The fragment is positioned so the variant is covered. The variant
    may appear in R1, R2, or both depending on fragment position.

    Args:
        ref_seq: Reference sequence
        var_pos: Variant position (0-based)
        ref_allele: Reference allele
        alt_allele: Alternative allele
        use_alt: If True, use alt_allele; else use ref_allele
        read_length: Read length (default 150bp)
        insert_mean: Mean insert size
        insert_std: Insert size std dev
        error_rate: Sequencing error rate

    Returns:
        ((r1_seq, r1_qual), (r2_seq, r2_qual))

    Fragment layout:
        |<----- insert_size ----->|
        |----R1----|              |----R2----|
        5'---------|==VARIANT==|----------3'
                   ^var_pos
    """
    allele = alt_allele if use_alt else ref_allele
    insert_size = generate_insert_size(insert_mean, insert_std, read_length + 50)

    # Position fragment so variant is inside
    # Random offset within valid range
    max_offset = insert_size - len(allele) - 10
    var_offset = random.randint(10, max(10, max_offset))
    frag_start = var_pos - var_offset
    frag_start = max(0, frag_start)
    frag_end = frag_start + insert_size

    # Build fragment with allele inserted
    left_of_var = ref_seq[frag_start:var_pos]
    right_of_var = ref_seq[var_pos + len(ref_allele):frag_end]
    fragment = left_of_var + allele + right_of_var

    # Ensure fragment is correct length
    if len(fragment) < insert_size:
        # Pad with reference
        fragment += ref_seq[frag_start + len(fragment):frag_start + insert_size]
    fragment = fragment[:insert_size]

    # R1: first read_length bases (forward strand)
    r1_seq = fragment[:read_length]

    # R2: last read_length bases (reverse complement)
    r2_seq = reverse_complement(fragment[-read_length:])

    # Add sequencing errors
    r1_seq = add_sequencing_errors(r1_seq, error_rate)
    r2_seq = add_sequencing_errors(r2_seq, error_rate)

    # Generate quality scores
    r1_qual = generate_quality_string(len(r1_seq))
    r2_qual = generate_quality_string(len(r2_seq))

    return ((r1_seq, r1_qual), (r2_seq, r2_qual))
```

#### 4. write_paired_fastq()
```python
def write_paired_fastq(
    r1_handle,
    r2_handle,
    read_id: str,
    r1_data: Tuple[str, str],
    r2_data: Tuple[str, str]
):
    """Write paired FASTQ records with proper naming."""
    r1_seq, r1_qual = r1_data
    r2_seq, r2_qual = r2_data

    # R1 file (forward)
    r1_handle.write(f"@{read_id}/1\n")
    r1_handle.write(f"{r1_seq}\n")
    r1_handle.write("+\n")
    r1_handle.write(f"{r1_qual}\n")

    # R2 file (reverse)
    r2_handle.write(f"@{read_id}/2\n")
    r2_handle.write(f"{r2_seq}\n")
    r2_handle.write("+\n")
    r2_handle.write(f"{r2_qual}\n")
```

#### 5. align_paired_with_bwa()
```python
def align_paired_with_bwa(
    ref_fasta: str,
    r1_fastq: str,
    r2_fastq: str,
    output_bam: str,
    threads: int = 4
) -> str:
    """Align paired-end reads with BWA MEM."""
    import pysam

    # Check BWA index
    if not Path(f"{ref_fasta}.bwt").exists():
        subprocess.run(['bwa', 'index', ref_fasta], check=True)

    # Align paired-end
    sam_file = output_bam.replace('.bam', '.sam')
    with open(sam_file, 'w') as sam:
        subprocess.run([
            'bwa', 'mem',
            '-t', str(threads),
            ref_fasta,
            r1_fastq,  # ← R1
            r2_fastq   # ← R2
        ], stdout=sam, check=True)

    # Convert to sorted BAM
    pysam.view('-bS', '-o', output_bam, sam_file, catch_stdout=False)
    sorted_bam = output_bam.replace('.bam', '.sorted.bam')
    pysam.sort('-o', sorted_bam, output_bam)
    pysam.index(sorted_bam)

    Path(sam_file).unlink()
    Path(output_bam).unlink()

    return sorted_bam
```

#### 6. run_wasp2_paired_pipeline()
```python
def run_wasp2_paired_pipeline(
    bam_file: str,
    vcf_file: str,
    ref_fasta: str,
    output_dir: Path
) -> str:
    """Run WASP2 unified pipeline on paired-end data."""
    import sys
    sys.path.insert(0, 'src')
    from mapping.run_mapping import run_make_remap_reads_unified

    # Run unified pipeline (handles paired-end internally)
    stats = run_make_remap_reads_unified(
        bam_file=bam_file,
        variant_file=vcf_file,
        samples=None,  # No sample filtering for simulation
        out_dir=str(output_dir),
        threads=4,
        compression_threads=2,
        use_parallel=True
    )

    print(f"WASP2 stats: {stats}")
    return str(output_dir / "remap_r1.fq.gz")
```

---

## Test Cases

### Must Test:
| Scenario | Description | Validation |
|----------|-------------|------------|
| Variant in R1 only | Fragment positioned so var in first 150bp | Check R1 has variant, R2 is ref |
| Variant in R2 only | Fragment positioned so var in last 150bp | Check R2 has variant, R1 is ref |
| Variant in both | Short insert, var near middle | Both reads show variant |
| Proper pairing | Read names match | `samtools flagstat` shows proper pairs |
| Insert size distribution | Mean ~300bp | Extract from BAM, check distribution |

### Test Script:
```python
def test_paired_end_simulation():
    """Quick validation of paired-end simulation."""
    import pysam

    # Run minimum tier
    # ... run simulation ...

    # Check BAM for proper pairs
    bam = pysam.AlignmentFile("output/aligned.sorted.bam")
    stats = {"proper_pair": 0, "total": 0}

    for read in bam.fetch():
        stats["total"] += 1
        if read.is_proper_pair:
            stats["proper_pair"] += 1

    proper_rate = stats["proper_pair"] / stats["total"]
    print(f"Proper pair rate: {proper_rate:.1%}")
    assert proper_rate > 0.90, "Too few proper pairs!"

    # Check insert sizes
    insert_sizes = []
    for read in bam.fetch():
        if read.is_proper_pair and read.is_read1:
            insert_sizes.append(abs(read.template_length))

    mean_insert = np.mean(insert_sizes)
    std_insert = np.std(insert_sizes)
    print(f"Insert size: {mean_insert:.0f} ± {std_insert:.0f}")
    assert 250 < mean_insert < 350, "Insert size out of range!"
```

---

## Step-by-Step Implementation

### Step 1: Setup
```bash
source /iblm/netapp/home/jjaureguy/mambaforge/etc/profile.d/conda.sh
conda activate WASP2_dev2
cd /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp

git checkout sim/paired-end
git pull origin sim/paired-end

mkdir -p simulation
```

### Step 2: Create the Module
```bash
# Create simulation/simulate_paired_end_ase.py
# Use the specifications above
```

### Step 3: Test Minimum Tier
```bash
python simulation/simulate_paired_end_ase.py \
    --tier minimum \
    --workdir /tmp/pe_test \
    --keep

# Validate outputs
ls /tmp/pe_test/
# Should see: synthetic_R1.fq, synthetic_R2.fq, aligned.sorted.bam, etc.

# Check read pairing
samtools flagstat /tmp/pe_test/aligned.sorted.bam
# Should show high "properly paired" percentage
```

### Step 4: Validate Insert Sizes
```bash
python -c "
import pysam
import numpy as np

bam = pysam.AlignmentFile('/tmp/pe_test/aligned.sorted.bam')
inserts = [abs(r.template_length) for r in bam if r.is_proper_pair and r.is_read1]
print(f'Insert size: {np.mean(inserts):.0f} ± {np.std(inserts):.0f} bp')
print(f'Range: {min(inserts)} - {max(inserts)} bp')
"
```

### Step 5: Run Through WASP2
```bash
# The simulation should run the full WASP2 pipeline
# Check that unified pipeline handles paired-end correctly
```

### Step 6: Compare with Single-End
```bash
# Run same variants with single-end (existing) and paired-end (new)
# Results should be similar (both testing same ground truth)
```

---

## Commit Template

```bash
git add simulation/simulate_paired_end_ase.py
git add simulation/test_paired_end.py  # if you create tests

git commit -m "feat: add paired-end simulation module

Implements paired-end read generation for WASP2 validation:
- R1/R2 FASTQ output with proper /1 /2 suffixes
- Configurable insert size (default: 300±50bp)
- Variants can appear in R1, R2, or both mates
- BWA paired-end alignment integration
- WASP2 unified pipeline integration

Test results (minimum tier):
- Proper pair rate: XX%
- Insert size: XXX ± XX bp
- Pass rate: XX%
"

git push origin sim/paired-end
```

---

## Success Criteria

- [ ] Generates valid R1.fq and R2.fq files
- [ ] Read names match between R1 and R2 (with /1 and /2 suffixes)
- [ ] Insert size distribution: mean ~300bp, std ~50bp
- [ ] BWA reports >90% properly paired
- [ ] WASP2 unified pipeline runs successfully
- [ ] Simulation pass rate >90% (comparable to single-end)
- [ ] Code committed to sim/paired-end branch

---

## Handoff

When complete:
1. Push to `sim/paired-end`
2. Note pass rate and insert size distribution
3. Ready for Agent 4 to run and generate metrics
