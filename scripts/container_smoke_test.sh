#!/bin/bash
# =============================================================================
# WASP2 Container Smoke Test
# =============================================================================
# Validates WASP2 tools inside a Docker/Singularity container.
# Designed to be bundled at /opt/wasp2/scripts/container_smoke_test.sh
# with test data at /opt/wasp2/test-data/ (copied from tests/shared_data/).
#
# When run outside a container, uses tests/shared_data/ directly.
#
# Usage (container):
#   docker run wasp2:test /opt/wasp2/scripts/container_smoke_test.sh
#   singularity exec wasp2.sif /opt/wasp2/scripts/container_smoke_test.sh
#
# Usage (local):
#   bash scripts/container_smoke_test.sh
#
# Exit codes:
#   0 = all tests passed
#   1 = one or more tests failed
# =============================================================================

set -euo pipefail

# Detect data directory: container path or repo path
if [[ -n "${WASP2_TEST_DATA_DIR:-}" ]]; then
    DATA_DIR="$WASP2_TEST_DATA_DIR"
elif [[ -d "/opt/wasp2/test-data" ]]; then
    DATA_DIR="/opt/wasp2/test-data"
elif [[ -d "tests/shared_data" ]]; then
    DATA_DIR="tests/shared_data"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    DATA_DIR="$(cd "$SCRIPT_DIR/../tests/shared_data" 2>/dev/null && pwd)" || {
        echo "ERROR: Cannot find test data directory"
        exit 1
    }
fi

TMP_DIR=$(mktemp -d "${TMPDIR:-/tmp}/wasp2_container_smoke_XXXXXX")
PASS=0
FAIL=0

cleanup() {
    rm -rf "$TMP_DIR"
}
trap cleanup EXIT

echo "==================================================================="
echo " WASP2 Container Smoke Test"
echo "==================================================================="
echo "Data dir: $DATA_DIR"
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# Test 1: CLI binaries exist and print version
# ─────────────────────────────────────────────────────────────────────────────
echo "[1/7] Version checks..."

for cmd in wasp2-count wasp2-map wasp2-analyze; do
    if $cmd --version; then
        echo "  PASS: $cmd --version"
        PASS=$((PASS + 1))
    else
        echo "  FAIL: $cmd --version"
        FAIL=$((FAIL + 1))
    fi
done
python - <<'PY'
import os
from importlib.metadata import version
import wasp2
import wasp2_rust

import subprocess
import re

installed = version("wasp2")
assert wasp2.__version__ == wasp2_rust.__version__ == installed
expected = os.environ.get("EXPECTED_WASP2_VERSION", installed)
assert installed == expected, (installed, expected)
for command in ("wasp2-count", "wasp2-map", "wasp2-analyze"):
    output = subprocess.check_output([command, "--version"], stderr=subprocess.STDOUT, text=True)
    match = re.search(r"WASP2\s*│\s*([^\s│]+)", output)
    assert match and match.group(1) == expected, (command, output, expected)
print(f"  PASS: Python, Rust and distribution versions agree: {installed}")
PY
PASS=$((PASS + 1))
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# Test 2: Python imports
# ─────────────────────────────────────────────────────────────────────────────
echo "[2/7] Python imports..."

if python -c "import wasp2_rust; print('Rust extension OK')" 2>/dev/null; then
    echo "  PASS: wasp2_rust import"
    PASS=$((PASS + 1))
else
    echo "  FAIL: wasp2_rust import"
    FAIL=$((FAIL + 1))
fi

if python -c "from counting.run_counting import run_count_variants; print('counting OK')" 2>/dev/null; then
    echo "  PASS: counting module import"
    PASS=$((PASS + 1))
else
    echo "  FAIL: counting module import"
    FAIL=$((FAIL + 1))
fi
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# Test 3: Count variants with the bundled synthetic BAM/VCF
# ─────────────────────────────────────────────────────────────────────────────
echo "[3/7] Count variants with bundled BAM/VCF..."

if [[ -f "$DATA_DIR/sample1.bam" && -f "$DATA_DIR/variants.vcf.gz" ]]; then
    if wasp2-count count-variants \
        "$DATA_DIR/sample1.bam" \
        "$DATA_DIR/variants.vcf.gz" \
        --samples sample1 \
        --out "$TMP_DIR/counts.tsv" \
        --temp_loc "$TMP_DIR"; then

        python - "$TMP_DIR/counts.tsv" "$DATA_DIR/sample1.bam" "$DATA_DIR/variants.vcf.gz" <<'PYCOUNTS'
import csv
import sys
from collections import Counter
import pysam

# Independent Python reference: consume each read name at the first fetched SNP,
# with no base-quality filter, matching the original counting contract.
with pysam.AlignmentFile(sys.argv[2]) as bam, pysam.VariantFile(sys.argv[3]) as vcf:
    expected = {}
    seen = set()
    for variant in vcf:
        if variant.samples["sample1"]["GT"] not in ((0, 1), (1, 0)):
            continue
        bases = Counter()
        for read in bam.fetch(variant.contig, variant.start, variant.stop):
            identity = (variant.contig, read.query_name)
            if identity in seen:
                continue
            seen.add(identity)
            for query_pos, reference_pos in read.get_aligned_pairs(matches_only=True):
                if reference_pos == variant.start:
                    bases[read.query_sequence[query_pos]] += 1
        ref, alt = variant.ref, variant.alts[0]
        expected[(variant.contig, str(variant.pos), ref, alt)] = (bases[ref], bases[alt], sum(bases.values()) - bases[ref] - bases[alt])
with open(sys.argv[1]) as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
actual = {(r["chrom"], r["pos"], r["ref"], r["alt"]): tuple(int(r[k]) for k in ("ref_count", "alt_count", "other_count")) for r in rows}
assert len(rows) == len(actual) == len(expected) > 0
assert actual == expected, (actual, expected)
print("  PASS: exact keyed allele counts against independent BAM base counts")
PYCOUNTS
        PASS=$((PASS + 1))
    else
        echo "  FAIL: count-variants exited with error"
        FAIL=$((FAIL + 1))
    fi
else
    echo "  FAIL: Required test BAM/VCF missing"
    FAIL=$((FAIL + 1))
fi
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# Test 4: External tools (samtools, bcftools, bedtools)
# ─────────────────────────────────────────────────────────────────────────────
echo "[4/7] External tool availability..."

for tool in samtools bcftools bedtools bwa; do
    if command -v $tool > /dev/null 2>&1; then
        echo "  PASS: $tool available"
        PASS=$((PASS + 1))
    else
        echo "  FAIL: $tool not found"
        FAIL=$((FAIL + 1))
    fi
done
echo ""

# Exercise the original ATAC failure: reordered columns, nonnumeric region IDs,
# and both phase orientations. Compare phased and unphased CLI runs to catch
# silently ignored --phased as well as positional parsing and ignored grouping.
echo "[5/7] ATAC region analysis..."
python - "$TMP_DIR/regions.tsv" <<'PY'
import csv
import sys

with open(sys.argv[1], "w") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["ref", "alt", "region", "chrom", "start", "pos", "GT", "ref_count", "alt_count"])
    for pos, gt, ref, alt in [(101, "0|1", 18, 2), (201, "1|0", 3, 17),
                              (301, "0|1", 22, 18), (401, "1|0", 25, 15)]:
        region = "chr1:100-250" if pos < 300 else "chr1:300-450"
        writer.writerow(["A", "G", region, "chr1", pos - 1, pos, gt, ref, alt])
PY
wasp2-analyze find-imbalance "$TMP_DIR/regions.tsv" --phased \
    --region_col region --out_file "$TMP_DIR/phased.tsv"
wasp2-analyze find-imbalance "$TMP_DIR/regions.tsv" \
    --region_col region --out_file "$TMP_DIR/unphased.tsv"
python - "$TMP_DIR/phased.tsv" "$TMP_DIR/unphased.tsv" <<'PY'
import csv
import math
import sys
import numpy as np
from scipy.stats import chi2, false_discovery_control

def read(path):
    with open(path) as handle:
        return {row["region"]: row for row in csv.DictReader(handle, delimiter="\t")}

phased, unphased = map(read, sys.argv[1:])
assert set(phased) == set(unphased) == {"chr1:100-250", "chr1:300-450"}
for region, counts in [("chr1:100-250", (21, 19)), ("chr1:300-450", (47, 33))]:
    for results in (phased, unphased):
        row = results[region]
        assert (int(row["ref_count"]), int(row["alt_count"])) == counts
        assert int(row["snp_count"]) == 2
        assert all(math.isfinite(float(row[key])) for key in ["mu", "null_ll", "alt_ll", "pval", "fdr_pval"])
for results in (phased, unphased):
    pvalues = []
    for row in results.values():
        lrt = max(0.0, 2 * (float(row["alt_ll"]) - float(row["null_ll"])))
        assert math.isclose(float(row["lrt"]), lrt, rel_tol=1e-8, abs_tol=1e-8)
        pvalues.append(chi2.sf(lrt, 1))
    assert np.allclose([float(row["pval"]) for row in results.values()], pvalues, rtol=1e-8, atol=1e-8)
    assert np.allclose([float(row["fdr_pval"]) for row in results.values()], false_discovery_control(pvalues), rtol=1e-8, atol=1e-8)
assert any(abs(float(phased[k]["alt_ll"]) - float(unphased[k]["alt_ll"])) > 1e-6 for k in phased)
print("  PASS: region grouping, raw counts and phased likelihood through CLI")
PY
PASS=$((PASS + 1))
echo ""

# Mapping must run through the installed CLI, not merely import successfully.
echo "[6/7] Mapping and remapping..."
# The older shared BAM has no proper-pair flags. Generate a small, uniquely
# mappable proper-pair fixture so this test actually exercises allele swapping.
python - "$TMP_DIR" <<'PYFIXTURE'
import random
import sys
from pathlib import Path
import pysam
root = Path(sys.argv[1])
rng = random.Random(7341)
sequence = "".join(rng.choices("ACGT", k=3000))
(root / "reference.fa").write_text(">chr_smoke\n" + sequence + "\n")
pysam.faidx(str(root / "reference.fa"))
header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr_smoke", "LN": len(sequence)}]}
reads = []
for name, start in [("overlap", 250), ("unaffected", 800)]:
    for mate, position in enumerate((start, start + 200)):
        read = pysam.AlignedSegment()
        read.query_name = name
        read.flag = 99 if mate == 0 else 147
        read.reference_id = read.next_reference_id = 0
        read.reference_start = position
        read.next_reference_start = start + 200 if mate == 0 else start
        read.template_length = 300 if mate == 0 else -300
        read.mapping_quality = 60
        read.cigarstring = "100M"
        read.query_sequence = sequence[position:position + 100]
        read.query_qualities = pysam.qualitystring_to_array("I" * 100)
        reads.append(read)
with pysam.AlignmentFile(root / "mapping.bam", "wb", header=header) as bam:
    for read in sorted(reads, key=lambda r: r.reference_start):
        bam.write(read)
pysam.index(str(root / "mapping.bam"))
ref = sequence[300]
alt = next(base for base in "ACGT" if base != ref)
(root / "mapping.vcf").write_text(
    '##fileformat=VCFv4.2\n##contig=<ID=chr_smoke,length=3000>\n'
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n'
    f'chr_smoke\t301\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t0|1\n')
PYFIXTURE
wasp2-map make-reads "$TMP_DIR/mapping.bam" "$TMP_DIR/mapping.vcf" \
    --samples sample1 --paired --phased --out_dir "$TMP_DIR" \
    --temp_loc "$TMP_DIR" --out_json "$TMP_DIR/mapping.json"
python - "$TMP_DIR/mapping.json" <<'PYMAPPING'
import json
import sys
from collections import Counter
import pysam
with open(sys.argv[1]) as handle:
    paths = json.load(handle)
with pysam.AlignmentFile(paths["to_remap_bam"]) as bam:
    original = Counter((r.query_name, r.is_read1) for r in bam)
assert original, "fixture must exercise remapping"
root = __import__("pathlib").Path(sys.argv[1]).parent
with pysam.FastaFile(root / "reference.fa") as reference:
    sequence = reference.fetch("chr_smoke")
with pysam.VariantFile(root / "mapping.vcf") as variants:
    variant = next(variants)
ref_seq = sequence[250:350]
alt_seq = ref_seq[:50] + variant.alts[0] + ref_seq[51:]
mate_seq = sequence[450:550].translate(str.maketrans("ACGT", "TGCA"))[::-1]
expected_sequences = {"remap_fq1": Counter([alt_seq]), "remap_fq2": Counter([mate_seq])}
for key in ("remap_fq1", "remap_fq2"):
    with pysam.FastxFile(paths[key]) as fastq:
        records = list(fastq)
    assert Counter(r.sequence for r in records) == expected_sequences[key], (key, Counter(r.sequence for r in records), expected_sequences[key])
    assert all(r.quality == "I" * 100 for r in records)
    assert len({r.name for r in records}) == len(records)
print("  PASS: nonempty paired swapped FASTQs and remap BAM")
PYMAPPING
bwa index "$TMP_DIR/reference.fa"
bwa mem -I 300,10,350,250 "$TMP_DIR/reference.fa" "$TMP_DIR/mapping_swapped_alleles_r1.fq" \
    "$TMP_DIR/mapping_swapped_alleles_r2.fq" | samtools sort -o "$TMP_DIR/remapped.bam"
wasp2-map filter-remapped "$TMP_DIR/remapped.bam" --json "$TMP_DIR/mapping.json" \
    --out_bam "$TMP_DIR/filtered.bam"
python - "$TMP_DIR" <<'PYFILTER'
import sys
from pathlib import Path
from collections import Counter
import pysam
root = Path(sys.argv[1])
def identities(path):
    with pysam.AlignmentFile(path) as bam:
        return Counter((r.query_name, r.is_read1) for r in bam if not r.is_secondary and not r.is_supplementary)
original = identities(root / "mapping_keep.bam") + identities(root / "mapping_to_remap.bam")
kept = identities(root / "mapping_keep.bam")
filtered = identities(root / "filtered.bam")
assert filtered == original, (filtered, original)
assert not (kept - filtered), "unaffected reads must survive filtering"
with pysam.AlignmentFile(root / "filtered.bam") as bam:
    assert bam.has_index()
print("  PASS: indexed filtered BAM preserves unaffected read identities")
PYFILTER
PASS=$((PASS + 1))

# A moved allele-swapped pair and an incomplete pair must both be rejected.
python - "$TMP_DIR" <<'PYNEGATIVE'
import sys
from pathlib import Path
import pysam
root = Path(sys.argv[1])
with pysam.AlignmentFile(root / "remapped.bam") as bam:
    header = bam.header
    reads = list(bam)
with pysam.AlignmentFile(root / "missing.bam", "wb", header=header) as bam:
    for read in reads:
        if read.is_read2:
            bam.write(read)
with pysam.AlignmentFile(root / "moved.bam", "wb", header=header) as bam:
    for original in reads:
        read = pysam.AlignedSegment.fromstring(original.to_string(), header)
        if read.is_read1:
            read.reference_start += 1
        else:
            read.next_reference_start += 1
        bam.write(read)
PYNEGATIVE
for scenario in moved missing; do
    wasp2-map filter-remapped "$TMP_DIR/${scenario}.bam" --json "$TMP_DIR/mapping.json" \
        --out_bam "$TMP_DIR/filtered-${scenario}.bam"
done
python - "$TMP_DIR" <<'PYREJECT'
import sys
from pathlib import Path
from collections import Counter
import pysam
root = Path(sys.argv[1])
for scenario in ("moved", "missing"):
    with pysam.AlignmentFile(root / f"filtered-{scenario}.bam") as bam:
        kept = Counter((read.query_name, read.is_read1) for read in bam)
    assert kept == Counter({("unaffected", True): 1, ("unaffected", False): 1}), (scenario, kept)
print("  PASS: moved and incomplete remap pairs are rejected")
PYREJECT
PASS=$((PASS + 1))

echo "[7/7] Single-cell counting..."
python - "$TMP_DIR" <<'PYSCINPUT'
import sys
from pathlib import Path
import pysam
root = Path(sys.argv[1])
with pysam.AlignmentFile(root / "mapping.bam") as bam:
    header = bam.header
    reads = [r for r in bam if r.query_name == "overlap"]
with pysam.VariantFile(root / "mapping.vcf") as variants:
    variant = next(variants)
records = []
for cell, alt in [("cellA-1", False), ("cellB-1", True)]:
    for original in reads:
        read = pysam.AlignedSegment.fromstring(original.to_string(), header)
        read.query_name = cell
        read.set_tag("CB", cell)
        if alt and read.is_read1:
            bases = list(read.query_sequence)
            bases[50] = variant.alts[0]
            read.query_sequence = "".join(bases)
            read.query_qualities = pysam.qualitystring_to_array("I" * 100)
        records.append(read)
with pysam.AlignmentFile(root / "cells.bam", "wb", header=header) as bam:
    for read in sorted(records, key=lambda r: r.reference_start):
        bam.write(read)
pysam.index(str(root / "cells.bam"))
(root / "barcodes.tsv").write_text("cellA-1\ncellB-1\n")
(root / "features.bed").write_text("chr_smoke\t250\t550\tpeak1\n")
PYSCINPUT
wasp2-count count-variants-sc "$TMP_DIR/cells.bam" "$TMP_DIR/mapping.vcf" \
    "$TMP_DIR/barcodes.tsv" --samples sample1 --feature "$TMP_DIR/features.bed" \
    --out "$TMP_DIR/cells.h5ad" --temp_loc "$TMP_DIR"
python - "$TMP_DIR" <<'PYSC'
import os
import sys
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
root = Path(sys.argv[1])
data = ad.read_h5ad(root / "cells.h5ad")
assert data.shape == (2, 1), data.shape
assert list(data.obs_names) == ["cellA-1", "cellB-1"]
assert np.array_equal(data.layers["ref"].toarray(), [[1], [0]])
assert np.array_equal(data.layers["alt"].toarray(), [[0], [1]])
assert np.array_equal(data.X.toarray(), [[1], [1]])
assert "feature" in data.uns
if os.environ.get("WASP2_TEST_EXPORTS") == "1":
    data.write_zarr(root / "cells.zarr")
    reopened = ad.read_zarr(root / "cells.zarr")
    assert reopened.shape == data.shape
    assert np.array_equal(reopened.X.toarray(), data.X.toarray())
    frame = pd.read_csv(root / "phased.tsv", sep="\t")
    frame.to_parquet(root / "results.parquet", index=False)
    pd.testing.assert_frame_equal(pd.read_parquet(root / "results.parquet"), frame)
    print("  PASS: Zarr and Parquet export round trips")
print("  PASS: exact single-cell counts, barcodes, axes and feature metadata")
PYSC
PASS=$((PASS + 1))

# ─────────────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────────────
TOTAL=$((PASS + FAIL))
echo "==================================================================="
echo " Results: $PASS/$TOTAL passed"
echo "==================================================================="

if [[ $FAIL -gt 0 ]]; then
    echo " WARNING: $FAIL test(s) failed"
    exit 1
else
    echo " All container smoke tests passed!"
    exit 0
fi
