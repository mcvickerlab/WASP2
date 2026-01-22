#!/bin/bash
# Produce WASP2-rust counter outputs on a WASP1-filtered BAM (HG00731).
#
# This is a small helper job to generate:
#   paper/figure2/data/hg00731/wasp1_counts.filtered.tsv
#
# It assumes a completed WASP1 run exists under:
#   benchmarking/star_wasp_comparison/results/wasp1_*/wasp1_final_sorted.bam
#
#$ -N fig2_wasp1_counts
#$ -cwd
#$ -j y
#$ -V
#$ -o paper/figure2/logs/fig2_wasp1_counts.$JOB_ID.log
#$ -pe iblm 8
#$ -l h_vmem=8G
#$ -l h_rt=2:00:00

set -eo pipefail

REPO_ROOT="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp"
cd "$REPO_ROOT"

export PYTHONPATH="${REPO_ROOT}/src:${PYTHONPATH:-}"

source /iblm/netapp/home/jjaureguy/mambaforge/etc/profile.d/conda.sh
conda activate WASP2_dev2

THREADS="${NSLOTS:-8}"

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

WASP1_RESULTS_DIR="$(
    find_latest_dir_with_file \
        "$REPO_ROOT/benchmarking/star_wasp_comparison/results/wasp1_*" \
        "wasp1_final_sorted.bam" \
    || true
)"
if [ -z "${WASP1_RESULTS_DIR}" ]; then
  echo "ERROR: No WASP1 results dir with wasp1_final_sorted.bam found" >&2
  exit 1
fi

WASP1_BAM="${WASP1_RESULTS_DIR}/wasp1_final_sorted.bam"
OUT_WITH_RG="$REPO_ROOT/paper/figure2/data/hg00731/wasp1_filtered.with_rg.bam"
OUT_TSV="$REPO_ROOT/paper/figure2/data/hg00731/wasp1_counts.filtered.tsv"

REF_GENOME="$REPO_ROOT/paper/figure2/tools/hg38.fa"
VCF_ORIG="$REPO_ROOT/benchmarking/star_wasp_comparison/data/HG00731_het_only_chr.vcf.gz"
VCF_V42="$REPO_ROOT/paper/figure2/data/hg00731/HG00731_het_only_chr.vcf.v4.2.gz"
SAMPLE="HG00731"

mkdir -p "$REPO_ROOT/paper/figure2/data/hg00731" "$REPO_ROOT/paper/figure2/logs"

echo "Using WASP1 BAM: $WASP1_BAM"

echo "Ensuring GATK-compatible VCF exists: $VCF_V42"
if [ ! -f "$VCF_V42" ]; then
  python - <<PY
import pathlib
in_vcf = pathlib.Path("$VCF_ORIG")
ref_fai = pathlib.Path("$REF_GENOME.fai")
out_bgz = pathlib.Path("$VCF_V42")
tmp_vcf = pathlib.Path(str(out_bgz) + ".tmp.vcf")
lengths = {}
with ref_fai.open() as f:
    for line in f:
        fields = line.rstrip("\n").split("\t")
        if len(fields) >= 2:
            lengths[fields[0]] = int(fields[1])

lines = in_vcf.read_bytes()
import gzip
text = gzip.decompress(lines).decode("utf-8")
out_lines = []
for line in text.splitlines(True):
    if line.startswith("##fileformat=VCFv"):
        out_lines.append("##fileformat=VCFv4.2\n")
        continue
    if line.startswith("##contig=<ID=") and "length=" not in line:
        # parse ID
        try:
            inside = line[len("##contig=<") :].split(">", 1)[0]
            parts = dict(x.split("=", 1) for x in inside.split(",") if "=" in x)
            cid = parts.get("ID")
            if cid and cid in lengths:
                out_lines.append(f"##contig=<ID={cid},length={lengths[cid]}>\n")
                continue
        except Exception:
            pass
    out_lines.append(line)
tmp_vcf.write_text("".join(out_lines))
PY
  bgzip -f -c "$VCF_V42.tmp.vcf" > "$VCF_V42"
  tabix -f -p vcf "$VCF_V42"
  rm -f "$VCF_V42.tmp.vcf"
fi

echo "Adding read groups for GATK/phASER compatibility..."
gatk AddOrReplaceReadGroups \
    -I "$WASP1_BAM" \
    -O "$OUT_WITH_RG" \
    -RGID 1 -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM "$SAMPLE"
samtools index -@ "$THREADS" "$OUT_WITH_RG"

echo "Counting alleles (WASP2 counter) on WASP1-filtered BAM..."
python paper/figure2/scripts/wasp2_rust_ase_counts.py \
    --bam "$OUT_WITH_RG" \
    --vcf "$VCF_V42" \
    --sample "$SAMPLE" \
    --out "$OUT_TSV" \
    --threads "$THREADS" \
    --min-mapq 10 \
    --min-baseq 20

echo "Wrote: $OUT_TSV"
