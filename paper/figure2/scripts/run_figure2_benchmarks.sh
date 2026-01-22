#!/bin/bash
#$ -N figure2_bench
#$ -cwd
#$ -j y
#$ -V
#$ -o paper/figure2/logs/figure2_benchmarks.$JOB_ID.log
#$ -pe iblm 8
#$ -l h_vmem=8G
#$ -l h_rt=4:00:00

###############################################################################
# Figure 2 Benchmarks - WASP2-Rust vs GATK vs phASER
#
# Runs allele counting benchmarks for:
# - Panel A: Speed comparison (WASP2-Rust, GATK ASEReadCounter, phASER)
# - Panel B: Count comparison across tools
# - Panel C: Bias reduction (original vs remapped BAM)
#
# Datasets:
# - HG00731 RNA-seq (56M reads)
# - GM12878 ATAC-seq (159M reads)
###############################################################################

# NOTE: Avoid `set -u` here because conda's shell hooks are not nounset-safe.
set -eo pipefail

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------

REPO_ROOT="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp"
cd "$REPO_ROOT"

# Ensure the in-repo Python package (`src/wasp2`) is importable inside conda.
export PYTHONPATH="${REPO_ROOT}/src:${PYTHONPATH:-}"

# Which dataset(s) to run. Main-text Figure 2 is RNA-seq (HG00731).
# ATAC-seq support is optional/supplementary.
DATASET="${DATASET:-hg00731}"  # hg00731 | gm12878 | both

# Set to 1 to recompute "original" counts even if output files already exist.
# Useful when `paper/figure2/data/hg00731/original.bam` is regenerated.
FORCE_ORIGINAL="${FORCE_ORIGINAL:-0}"

# Optional: also produce WASP2-rust counter outputs on a WASP1-filtered BAM, if available.
# This is used for an additional sanity check / supplement comparing WASP1 vs WASP2 filtering.
INCLUDE_WASP1="${INCLUDE_WASP1:-0}"

# Create output directories
mkdir -p paper/figure2/data/hg00731
mkdir -p paper/figure2/data/gm12878
mkdir -p paper/figure2/logs

# Output files
TIMING_JSON="paper/figure2/data/timing_results.json"
LOG_FILE="paper/figure2/logs/benchmarks_$(date +%Y-%m-%d_%H-%M-%S).log"

# Redirect all output to log
exec > >(tee -a "$LOG_FILE") 2>&1

echo "=================================================="
echo "Figure 2 Benchmarks - Started at $(date)"
echo "=================================================="
echo "Threads: $NSLOTS"
echo "Host: $(hostname)"
echo "Dataset: $DATASET"
echo ""

# -----------------------------------------------------------------------------
# Activate conda environment
# -----------------------------------------------------------------------------

echo "[$(date +%H:%M:%S)] Activating WASP2_dev2 environment..."
source /iblm/netapp/home/jjaureguy/mambaforge/etc/profile.d/conda.sh
conda activate WASP2_dev2

# Verify tools
echo "[$(date +%H:%M:%S)] Verifying tools..."
which python
which gatk
echo ""

# -----------------------------------------------------------------------------
# Input Files
# -----------------------------------------------------------------------------

# Reference genome used by the STAR index for HG00731; includes a local `.dict` for GATK.
REF_GENOME="$REPO_ROOT/paper/figure2/tools/hg38.fa"

# HG00731 RNA-seq
find_latest_dir_with_file() {
    local glob_pattern="$1"
    local required_file="$2"
    local d
    for d in $(ls -d ${glob_pattern} 2>/dev/null | sort -r); do
        if [ -f "${d}/${required_file}" ]; then
            echo "${d}"
            return 0
        fi
    done
    return 1
}

HG00731_WASP2_RESULTS_DIR=""
HG00731_BAM_WASP_FILTERED=""
HG00731_BAM_WASP_FILTERED_WITH_RG="$REPO_ROOT/paper/figure2/data/hg00731/wasp_filtered.with_rg.bam"
HG00731_BAM_ORIGINAL_WITH_RG="$REPO_ROOT/paper/figure2/data/hg00731/original.with_rg.bam"
HG00731_BAM_ORIGINAL="$REPO_ROOT/paper/figure2/data/hg00731/original.bam"
if [[ "$DATASET" == "hg00731" || "$DATASET" == "both" ]]; then
    HG00731_WASP2_RESULTS_DIR="$(
        find_latest_dir_with_file \
            "$REPO_ROOT/benchmarking/star_wasp_comparison/results/wasp2rust_fair_snv_*" \
            "wasp_filtered.bam" \
        || true
    )"

    if [ -z "${HG00731_WASP2_RESULTS_DIR}" ]; then
        echo "ERROR: Could not find a WASP2-Rust FAIR results dir containing wasp_filtered.bam" >&2
        echo "Searched: $REPO_ROOT/benchmarking/star_wasp_comparison/results/wasp2rust_fair_snv_*" >&2
        exit 1
    fi

    HG00731_BAM_WASP_FILTERED="${HG00731_WASP2_RESULTS_DIR}/wasp_filtered.bam"

    # Prefer the exact STAR "original.bam" produced by the same FAIR run (if exported with KEEP_BAMS=1).
    if [ -f "${HG00731_WASP2_RESULTS_DIR}/original.bam" ]; then
        HG00731_BAM_ORIGINAL="${HG00731_WASP2_RESULTS_DIR}/original.bam"
        HG00731_BAM_ORIGINAL_WITH_RG="$REPO_ROOT/paper/figure2/data/hg00731/original_from_fair.with_rg.bam"
    fi
fi
HG00731_VCF_ORIG="$REPO_ROOT/benchmarking/star_wasp_comparison/data/HG00731_het_only_chr.vcf.gz"
# GATK 4.0.x chokes on VCFv4.3; rewrite header to VCFv4.2 (variants unchanged) and bgzip+tabix.
HG00731_VCF="$REPO_ROOT/paper/figure2/data/hg00731/HG00731_het_only_chr.vcf.v4.2.gz"
HG00731_SAMPLE="HG00731"

# GM12878 ATAC-seq
GM12878_BAM_ORIGINAL="/iblm/netapp/data1/aho/atac/Buenrostro2013/merged/GM12878_ATACseq_50k/GM12878_ATACseq_50k_merged.sorted.bam"
GM12878_BAM_REMAP="$REPO_ROOT/benchmarking/atac_gm12878/results/wasp2rust_snp_fixed_2025-12-15_03-43-24/remap_keep.bam"
GM12878_VCF="/iblm/netapp/data1/aho/variants/NA12878.vcf.gz"
GM12878_SAMPLE="NA12878"

# phASER tool
PHASER_DIR="$REPO_ROOT/benchmarking/phaser_tool/phaser"

# Verify input files
echo "[$(date +%H:%M:%S)] Verifying input files..."

if [ ! -f "$REF_GENOME" ]; then
    echo "ERROR: Missing required file: $REF_GENOME"
    exit 1
fi
echo "  OK: $REF_GENOME"

if [[ "$DATASET" == "hg00731" || "$DATASET" == "both" ]]; then
    for file in "$HG00731_BAM_WASP_FILTERED" "$HG00731_BAM_ORIGINAL" "$HG00731_VCF_ORIG"; do
        if [ ! -f "$file" ]; then
            echo "ERROR: Missing required HG00731 file: $file"
            exit 1
        fi
        echo "  OK: $file"
    done
fi

if [[ "$DATASET" == "gm12878" || "$DATASET" == "both" ]]; then
    for file in "$GM12878_BAM_ORIGINAL" "$GM12878_BAM_REMAP" "$GM12878_VCF"; do
        if [ ! -f "$file" ]; then
            echo "ERROR: Missing required GM12878 file: $file"
            exit 1
        fi
        echo "  OK: $file"
    done
fi
echo ""

# -----------------------------------------------------------------------------
# Helper Function: Time Command
# -----------------------------------------------------------------------------

time_command() {
    local name="$1"
    shift
    local cmd="$@"

    echo "[$(date +%H:%M:%S)] Running: $name" >&2
    echo "  Command: $cmd" >&2

    local start_time=$(date +%s.%N)
    # Send the command's stdout to stderr so command substitution captures only elapsed time.
    eval "$cmd" 1>&2
    local exit_code=$?
    local end_time=$(date +%s.%N)
    local elapsed=$(echo "$end_time - $start_time" | bc)

    if [ $exit_code -eq 0 ]; then
        echo "  SUCCESS: $name completed in ${elapsed}s" >&2
        echo "$elapsed"
        return 0
    else
        echo "  ERROR: $name failed with exit code $exit_code" >&2
        return $exit_code
    fi
}

# -----------------------------------------------------------------------------
# Prepare VCFs for tool compatibility
# -----------------------------------------------------------------------------

prepare_vcf_v42_bgzip() {
    local in_vcf="$1"
    local out_bgz="$2"

    echo "[$(date +%H:%M:%S)] Preparing GATK-compatible VCF: $out_bgz"
    mkdir -p "$(dirname "$out_bgz")"
    local tmp_bgz="${out_bgz}.tmp"
    local tmp_vcf="${out_bgz}.tmp.vcf"
    rm -f "$tmp_bgz" "${tmp_bgz}.tbi"
    rm -f "$tmp_vcf"

    # Rewrite header for tool compatibility:
    # - GATK 4.0.x requires VCFv4.2 (not 4.3)
    # - GATK requires ##contig lines to include length=
    python - <<PY
import pathlib

in_vcf = pathlib.Path("$in_vcf")
ref_fai = pathlib.Path("$REF_GENOME.fai")
tmp_vcf = pathlib.Path("$tmp_vcf")

lengths = {}
with ref_fai.open() as f:
    for line in f:
        fields = line.rstrip("\\n").split("\\t")
        if len(fields) >= 2:
            lengths[fields[0]] = int(fields[1])

header_lines = []
contig_ids_seen = set()
had_contig = False

def contig_line(cid: str) -> str | None:
    ln = lengths.get(cid)
    if ln is None:
        return None
    return f"##contig=<ID={cid},length={ln}>\\n"

with tmp_vcf.open("w") as out, open(in_vcf, "rb") as fin:
    # Read gzip/bgzip transparently by shell before python? No: input is gz; just read bytes and decode linewise.
    # We'll stream by invoking `zcat` outside python if needed; here use simple approach: rely on in_vcf being BGZF/gzip.
    import gzip
    f = gzip.open(fin, "rt")
    last_key = None
    for line in f:
        if not line.startswith("#"):
            # flush any missing contig lines before body
            for h in header_lines:
                out.write(h)
            # First record: also handle possible duplicate positions
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 5:
                ref = fields[3]
                alt = fields[4]
                # Keep SNPs only: avoids indels spanning SNP positions in GATK ASEReadCounter.
                if len(ref) != 1 or len(alt) != 1 or "," in alt:
                    # Skip this first record; keep scanning until we find the first SNP record.
                    continue
                last_key = (fields[0], fields[1])
                out.write(line)
            break

        if line.startswith("##fileformat=VCFv4.3"):
            header_lines.append("##fileformat=VCFv4.2\\n")
            continue

        if line.startswith("##contig=<ID="):
            had_contig = True
            # parse ID field
            inner = line.strip()[10:-1]  # remove '##contig=<' and trailing '>'
            parts = dict(p.split("=", 1) for p in inner.split(",") if "=" in p)
            cid = parts.get("ID")
            if not cid:
                continue
            contig_ids_seen.add(cid)
            if "length" in parts:
                header_lines.append(line)
            else:
                cl = contig_line(cid)
                header_lines.append(cl if cl is not None else line)
            continue

        if line.startswith("#CHROM"):
            if not had_contig:
                # Insert contig lines for all reference contigs (safe superset).
                for cid in lengths:
                    header_lines.append(contig_line(cid))
            header_lines.append(line)
            continue

        header_lines.append(line)

    # write remainder
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 5:
            continue
        ref = fields[3]
        alt = fields[4]
        if len(ref) != 1 or len(alt) != 1 or "," in alt:
            continue
        key = (fields[0], fields[1])
        # VCF is expected to be coordinate-sorted; when the same POS appears multiple
        # times, keep only the first record so GATK sees a single VariantContext.
        if key == last_key:
            continue
        last_key = key
        out.write(line)
PY

    bgzip -c "$tmp_vcf" > "$tmp_bgz"
    tabix -f -p vcf "$tmp_bgz"
    rm -f "$tmp_vcf"
    mv -f "$tmp_bgz" "$out_bgz"
    mv -f "${tmp_bgz}.tbi" "${out_bgz}.tbi"
}

prepare_bam_with_read_groups() {
    local in_bam="$1"
    local out_bam="$2"
    local sample="$3"

    if [ -f "$out_bam" ] && [ -f "${out_bam}.bai" ]; then
        return 0
    fi

    echo "[$(date +%H:%M:%S)] Adding read groups for GATK: $out_bam" >&2
    mkdir -p "$(dirname "$out_bam")"
    gatk AddOrReplaceReadGroups \
        -I "$in_bam" \
        -O "$out_bam" \
        -RGID 1 \
        -RGLB lib1 \
        -RGPL ILLUMINA \
        -RGPU unit1 \
        -RGSM "$sample"
    samtools index -@ "${NSLOTS:-1}" "$out_bam"
}

# -----------------------------------------------------------------------------
# Benchmark 1: HG00731 RNA-seq (WASP2-Rust)
# -----------------------------------------------------------------------------

echo "=================================================="
echo "BENCHMARK 1: HG00731 RNA-seq - WASP2-Rust"
echo "=================================================="

HG00731_WASP2_OUT_FILTERED="paper/figure2/data/hg00731/wasp2_counts.filtered.tsv"
HG00731_WASP2_OUT_ORIGINAL="paper/figure2/data/hg00731/wasp2_counts.original.tsv"

time_wasp2_hg00731=0
time_gatk_hg00731=0
time_phaser_hg00731=0
if [[ "$DATASET" == "hg00731" || "$DATASET" == "both" ]]; then
    prepare_vcf_v42_bgzip "$HG00731_VCF_ORIG" "$HG00731_VCF"
    prepare_bam_with_read_groups "$HG00731_BAM_WASP_FILTERED" "$HG00731_BAM_WASP_FILTERED_WITH_RG" "$HG00731_SAMPLE"
    # Ensure the STAR-only "original" BAM has read groups for GATK/phASER.
    # Regenerate if missing or if the underlying original BAM changed.
    if [ ! -f "$HG00731_BAM_ORIGINAL_WITH_RG" ] || [ "$HG00731_BAM_ORIGINAL" -nt "$HG00731_BAM_ORIGINAL_WITH_RG" ]; then
        prepare_bam_with_read_groups "$HG00731_BAM_ORIGINAL" "$HG00731_BAM_ORIGINAL_WITH_RG" "$HG00731_SAMPLE"
    fi
    time_wasp2_hg00731=$(time_command "WASP2-Rust HG00731 (filtered)" \
        "python paper/figure2/scripts/wasp2_rust_ase_counts.py \
            --bam '$HG00731_BAM_WASP_FILTERED_WITH_RG' \
            --vcf '$HG00731_VCF' \
            --sample '$HG00731_SAMPLE' \
            --out '$HG00731_WASP2_OUT_FILTERED' \
            --threads $NSLOTS \
            --min-mapq 10 \
            --min-baseq 20" \
        || echo "0")
    cp -f "$HG00731_WASP2_OUT_FILTERED" "paper/figure2/data/hg00731/wasp2_counts.tsv"

    if [ "${FORCE_ORIGINAL}" = "1" ] || [ ! -f "$HG00731_WASP2_OUT_ORIGINAL" ]; then
        time_command "WASP2-Rust HG00731 (original)" \
            "python paper/figure2/scripts/wasp2_rust_ase_counts.py \
                --bam '$HG00731_BAM_ORIGINAL_WITH_RG' \
                --vcf '$HG00731_VCF' \
                --sample '$HG00731_SAMPLE' \
                --out '$HG00731_WASP2_OUT_ORIGINAL' \
                --threads $NSLOTS \
                --min-mapq 10 \
                --min-baseq 20" \
            || true
    fi
fi

echo "WASP2-Rust HG00731 (filtered): ${time_wasp2_hg00731}s"
echo ""

# -----------------------------------------------------------------------------
# Optional: WASP1-filtered BAM counts with WASP2 counter (for filtering comparison)
# -----------------------------------------------------------------------------

if [[ "$DATASET" == "hg00731" || "$DATASET" == "both" ]] && [ "${INCLUDE_WASP1}" = "1" ]; then
    echo "=================================================="
    echo "OPTIONAL: HG00731 RNA-seq - WASP1-filtered BAM counts (WASP2 counter)"
    echo "=================================================="

    WASP1_RESULTS_DIR="$(
        find_latest_dir_with_file \
            "$REPO_ROOT/benchmarking/star_wasp_comparison/results/wasp1_*" \
            "wasp1_final_sorted.bam" \
        || true
    )"
    if [ -z "${WASP1_RESULTS_DIR}" ]; then
        echo "WARN: INCLUDE_WASP1=1 but no WASP1 results dir with wasp1_final_sorted.bam found"
    else
        HG00731_WASP1_BAM="${WASP1_RESULTS_DIR}/wasp1_final_sorted.bam"
        HG00731_WASP1_WITH_RG="$REPO_ROOT/paper/figure2/data/hg00731/wasp1_filtered.with_rg.bam"
        HG00731_WASP1_OUT_FILTERED="paper/figure2/data/hg00731/wasp1_counts.filtered.tsv"
        prepare_bam_with_read_groups "$HG00731_WASP1_BAM" "$HG00731_WASP1_WITH_RG" "$HG00731_SAMPLE"
        time_command "WASP2-Rust counter on WASP1-filtered BAM" \
            "python paper/figure2/scripts/wasp2_rust_ase_counts.py \
                --bam '$HG00731_WASP1_WITH_RG' \
                --vcf '$HG00731_VCF' \
                --sample '$HG00731_SAMPLE' \
                --out '$HG00731_WASP1_OUT_FILTERED' \
                --threads $NSLOTS \
                --min-mapq 10 \
                --min-baseq 20" \
            || true
        echo "Wrote: $HG00731_WASP1_OUT_FILTERED"
    fi
    echo ""
fi

# -----------------------------------------------------------------------------
# Benchmark 2: HG00731 RNA-seq (GATK ASEReadCounter)
# -----------------------------------------------------------------------------

echo "=================================================="
echo "BENCHMARK 2: HG00731 RNA-seq - GATK ASEReadCounter"
echo "=================================================="

HG00731_GATK_OUT="paper/figure2/data/hg00731/gatk_counts.filtered.table"
HG00731_GATK_ORIG_OUT="paper/figure2/data/hg00731/gatk_counts.original.table"

if [[ "$DATASET" == "hg00731" || "$DATASET" == "both" ]]; then
    time_gatk_hg00731=$(time_command "GATK HG00731" \
        "gatk ASEReadCounter \
            -R '$REF_GENOME' \
            -I '$HG00731_BAM_WASP_FILTERED_WITH_RG' \
            -V '$HG00731_VCF' \
            -O '$HG00731_GATK_OUT' \
            --min-mapping-quality 10 \
            --min-base-quality 20" \
        || echo "0")

    if [ "${FORCE_ORIGINAL}" = "1" ] || [ ! -f "$HG00731_GATK_ORIG_OUT" ]; then
        time_command "GATK HG00731 (original)" \
            "gatk ASEReadCounter \
                -R '$REF_GENOME' \
                -I '$HG00731_BAM_ORIGINAL_WITH_RG' \
                -V '$HG00731_VCF' \
                -O '$HG00731_GATK_ORIG_OUT' \
                --min-mapping-quality 10 \
                --min-base-quality 20" \
            || true
    fi
fi

echo "GATK HG00731: ${time_gatk_hg00731}s"
echo ""

# -----------------------------------------------------------------------------
# Benchmark 3: HG00731 RNA-seq (phASER)
# -----------------------------------------------------------------------------

echo "=================================================="
echo "BENCHMARK 3: HG00731 RNA-seq - phASER"
echo "=================================================="

HG00731_PHASER_PREFIX_FILTERED="paper/figure2/data/hg00731/phaser.filtered"
HG00731_PHASER_PREFIX_ORIGINAL="paper/figure2/data/hg00731/phaser.original"

if [[ "$DATASET" == "hg00731" || "$DATASET" == "both" ]]; then
    time_phaser_hg00731=$(time_command "phASER HG00731 (filtered)" \
        "python '$PHASER_DIR/phaser.py' \
            --vcf '$HG00731_VCF' \
            --bam '$HG00731_BAM_WASP_FILTERED_WITH_RG' \
            --sample '$HG00731_SAMPLE' \
            --mapq 10 \
            --baseq 20 \
            --paired_end 1 \
            --threads $NSLOTS \
            --o '$HG00731_PHASER_PREFIX_FILTERED'" \
        || echo "0")

    # Backwards-compatible prefix expected by older scripts: write `phaser.*` as copies of filtered outputs.
    # Use bash string ops to avoid regex escaping issues.
    for f in paper/figure2/data/hg00731/phaser.filtered.*; do
        [ -e "$f" ] || continue
        base="$(basename "$f")"
        suffix="${base#phaser.filtered}"  # e.g. ".allelic_counts.txt"
        dest="paper/figure2/data/hg00731/phaser${suffix}"
        if [ "$f" != "$dest" ]; then
            cp -f "$f" "$dest" || true
        fi
    done

    if [ "${FORCE_ORIGINAL}" = "1" ] || [ ! -f "${HG00731_PHASER_PREFIX_ORIGINAL}.allelic_counts.txt" ]; then
        time_command "phASER HG00731 (original)" \
            "python '$PHASER_DIR/phaser.py' \
                --vcf '$HG00731_VCF' \
                --bam '$HG00731_BAM_ORIGINAL_WITH_RG' \
                --sample '$HG00731_SAMPLE' \
                --mapq 10 \
                --baseq 20 \
                --paired_end 1 \
                --threads $NSLOTS \
                --o '$HG00731_PHASER_PREFIX_ORIGINAL'" \
            || true
    fi
fi

echo "phASER HG00731: ${time_phaser_hg00731}s"
echo ""

# -----------------------------------------------------------------------------
# Benchmark 4: GM12878 ATAC-seq (WASP2-Rust)
# -----------------------------------------------------------------------------

echo "=================================================="
echo "BENCHMARK 4: GM12878 ATAC-seq - WASP2-Rust"
echo "=================================================="

GM12878_WASP2_OUT="paper/figure2/data/gm12878/wasp2_counts.tsv"

time_wasp2_gm12878=0
time_gatk_gm12878=0
time_phaser_gm12878=0
if [[ "$DATASET" == "gm12878" || "$DATASET" == "both" ]]; then
    time_wasp2_gm12878=$(time_command "WASP2-Rust GM12878" \
        "python paper/figure2/scripts/wasp2_rust_ase_counts.py \
            --bam '$GM12878_BAM_REMAP' \
            --vcf '$GM12878_VCF' \
            --sample '$GM12878_SAMPLE' \
            --out '$GM12878_WASP2_OUT' \
            --threads $NSLOTS \
            --min-mapq 10 \
            --min-baseq 20" \
        || echo "0")
fi

echo "WASP2-Rust GM12878: ${time_wasp2_gm12878}s"
echo ""

# -----------------------------------------------------------------------------
# Benchmark 5: GM12878 ATAC-seq (GATK ASEReadCounter)
# -----------------------------------------------------------------------------

echo "=================================================="
echo "BENCHMARK 5: GM12878 ATAC-seq - GATK ASEReadCounter"
echo "=================================================="

GM12878_GATK_OUT="paper/figure2/data/gm12878/gatk_counts.table"

if [[ "$DATASET" == "gm12878" || "$DATASET" == "both" ]]; then
    time_gatk_gm12878=$(time_command "GATK GM12878" \
        "gatk ASEReadCounter \
            -R '$REF_GENOME' \
            -I '$GM12878_BAM_REMAP' \
            -V '$GM12878_VCF' \
            -O '$GM12878_GATK_OUT' \
            --min-mapping-quality 10 \
            --min-base-quality 20" \
        || echo "0")
fi

echo "GATK GM12878: ${time_gatk_gm12878}s"
echo ""

# -----------------------------------------------------------------------------
# Benchmark 6: GM12878 ATAC-seq (phASER)
# -----------------------------------------------------------------------------

echo "=================================================="
echo "BENCHMARK 6: GM12878 ATAC-seq - phASER"
echo "=================================================="

GM12878_PHASER_PREFIX="paper/figure2/data/gm12878/phaser"

if [[ "$DATASET" == "gm12878" || "$DATASET" == "both" ]]; then
    time_phaser_gm12878=$(time_command "phASER GM12878" \
        "python '$PHASER_DIR/phaser.py' \
            --vcf '$GM12878_VCF' \
            --bam '$GM12878_BAM_REMAP' \
            --sample '$GM12878_SAMPLE' \
            --mapq 10 \
            --baseq 20 \
            --paired_end 1 \
            --threads $NSLOTS \
            --o '$GM12878_PHASER_PREFIX'" \
        || echo "0")
fi

echo "phASER GM12878: ${time_phaser_gm12878}s"
echo ""

# -----------------------------------------------------------------------------
# Save Timing Results
# -----------------------------------------------------------------------------

echo "=================================================="
echo "Saving timing results to JSON..."
echo "=================================================="

if [[ "$DATASET" == "hg00731" ]]; then
  cat > "$TIMING_JSON" <<EOF
{
  "timestamp": "$(date -Iseconds)",
  "threads": $NSLOTS,
  "host": "$(hostname)",
  "datasets": {
    "hg00731_rnaseq": {
      "description": "HG00731 RNA-seq, 56M reads",
      "bam": "$HG00731_BAM_WASP_FILTERED",
      "vcf": "$HG00731_VCF",
      "sample": "$HG00731_SAMPLE",
      "timing": {
        "wasp2_rust_s": $time_wasp2_hg00731,
        "gatk_s": $time_gatk_hg00731,
        "phaser_s": $time_phaser_hg00731
      }
    }
  }
}
EOF
elif [[ "$DATASET" == "gm12878" ]]; then
  cat > "$TIMING_JSON" <<EOF
{
  "timestamp": "$(date -Iseconds)",
  "threads": $NSLOTS,
  "host": "$(hostname)",
  "datasets": {
    "gm12878_atacseq": {
      "description": "GM12878 ATAC-seq, 159M reads",
      "bam": "$GM12878_BAM_REMAP",
      "vcf": "$GM12878_VCF",
      "sample": "$GM12878_SAMPLE",
      "timing": {
        "wasp2_rust_s": $time_wasp2_gm12878,
        "gatk_s": $time_gatk_gm12878,
        "phaser_s": $time_phaser_gm12878
      }
    }
  }
}
EOF
else
  cat > "$TIMING_JSON" <<EOF
{
  "timestamp": "$(date -Iseconds)",
  "threads": $NSLOTS,
  "host": "$(hostname)",
  "datasets": {
    "hg00731_rnaseq": {
      "description": "HG00731 RNA-seq, 56M reads",
      "bam": "$HG00731_BAM_WASP_FILTERED",
      "vcf": "$HG00731_VCF",
      "sample": "$HG00731_SAMPLE",
      "timing": {
        "wasp2_rust_s": $time_wasp2_hg00731,
        "gatk_s": $time_gatk_hg00731,
        "phaser_s": $time_phaser_hg00731
      }
    },
    "gm12878_atacseq": {
      "description": "GM12878 ATAC-seq, 159M reads",
      "bam": "$GM12878_BAM_REMAP",
      "vcf": "$GM12878_VCF",
      "sample": "$GM12878_SAMPLE",
      "timing": {
        "wasp2_rust_s": $time_wasp2_gm12878,
        "gatk_s": $time_gatk_gm12878,
        "phaser_s": $time_phaser_gm12878
      }
    }
  }
}
EOF
fi

echo "Saved timing results to: $TIMING_JSON"
cat "$TIMING_JSON"
echo ""

# -----------------------------------------------------------------------------
# Optional: Auto-generate merged tables + plots
# -----------------------------------------------------------------------------

FINALIZE="${FINALIZE:-0}"
if [ "${FINALIZE}" = "1" ]; then
    echo "=================================================="
    echo "FINALIZE=1: Generating merged tables + plots"
    echo "=================================================="

    # Merge tool counts (Panel B input)
    if [[ "$DATASET" == "hg00731" ]]; then
        python paper/figure2/scripts/generate_count_comparison.py --dataset hg00731
        python paper/figure2/scripts/generate_before_after_counts.py --dataset hg00731
        python paper/figure2/scripts/generate_delta_counts.py --dataset hg00731
        python paper/figure2/scripts/generate_figure2.py --dataset hg00731 --panel-c-style summary_bars
    elif [[ "$DATASET" == "gm12878" ]]; then
        python paper/figure2/scripts/generate_count_comparison.py --dataset gm12878
        python paper/figure2/scripts/generate_bias_comparison.py --dataset gm12878 || true
        python paper/figure2/scripts/generate_figure2.py --dataset gm12878 --panel-c-style summary_bars
    elif [[ "$DATASET" == "both" ]]; then
        python paper/figure2/scripts/generate_count_comparison.py --dataset both
        python paper/figure2/scripts/generate_before_after_counts.py --dataset hg00731
        python paper/figure2/scripts/generate_delta_counts.py --dataset hg00731
        python paper/figure2/scripts/generate_figure2.py --dataset both --panel-c-style summary_bars
    fi

    echo "FINALIZE done."
    echo ""
fi

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------

echo "=================================================="
echo "Figure 2 Benchmarks Complete - $(date)"
echo "=================================================="
echo ""
echo "Output files generated:"
echo "  Timing:       $TIMING_JSON"
echo "  HG00731:"
echo "    WASP2 (filtered): $HG00731_WASP2_OUT_FILTERED"
echo "    WASP2 (original): $HG00731_WASP2_OUT_ORIGINAL"
echo "    GATK:       $HG00731_GATK_OUT"
echo "    phASER (filtered): ${HG00731_PHASER_PREFIX_FILTERED}.allelic_counts.txt"
echo "    phASER (original): ${HG00731_PHASER_PREFIX_ORIGINAL}.allelic_counts.txt"
echo "  GM12878:"
echo "    WASP2:      $GM12878_WASP2_OUT"
echo "    GATK:       $GM12878_GATK_OUT"
echo "    phASER:     ${GM12878_PHASER_PREFIX}.haplotypic_counts.txt"
echo ""
echo "Next steps:"
echo "  1. Run generate_count_comparison.py to merge counts"
echo "  2. Run generate_bias_comparison.py for bias analysis"
echo "  3. Run generate_figure2.py to create plots"
echo ""
echo "Log file: $LOG_FILE"
echo "=================================================="
