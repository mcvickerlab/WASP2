import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"

for p in (ROOT, SRC):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))


def _parse_total_seqs_from_name(name: str) -> int:
    # {orig}_WASP_{pos1}_{pos2}_{seq}_{total}[...]/1
    core = name[:-2] if name.endswith("/1") or name.endswith("/2") else name
    suffix = core.split("_WASP_", 1)[1]
    return int(suffix.split("_")[3])


@pytest.mark.unit
def test_quickbench_indel_trim_invariants(tmp_path: Path) -> None:
    """INDEL-mode produces N+1 trim-combos for a +2bp insertion and preserves read length."""
    wasp2_rust = pytest.importorskip("wasp2_rust")

    # Skip if benchmarking module not available (not included in release)
    try:
        from benchmarking.quickbench.fastq_utils import iter_fastq
        from benchmarking.quickbench.synthetic_dataset import (
            quickbench_indel_variants,
            write_bed,
            write_synthetic_bam_indel,
        )
    except ImportError:
        pytest.skip("benchmarking module not available (not included in release)")

    import pysam

    bam = tmp_path / "synthetic_indel.bam"
    bed = tmp_path / "variants_indel.bed"
    write_synthetic_bam_indel(bam)
    variants = quickbench_indel_variants()
    write_bed(bed, variants)

    out_dir = tmp_path / "unified"
    out_dir.mkdir()
    out_r1 = out_dir / "r1.fq"
    out_r2 = out_dir / "r2.fq"

    wasp2_rust.unified_make_reads_py(
        str(bam),
        str(bed),
        str(out_r1),
        str(out_r2),
        max_seqs=256,
        threads=1,
        compression_threads=1,
        compress_output=False,
        indel_mode=True,
        max_indel_size=50,
    )

    with pysam.AlignmentFile(str(bam), "rb") as bf:
        recs = [r for r in bf.fetch(until_eof=True) if r.query_name == "pairI"]
    r1 = next(r for r in recs if r.is_read1)
    r2 = next(r for r in recs if r.is_read2)
    r1_seq = r1.query_sequence
    r2_seq = r2.query_sequence
    read_len = len(r1_seq)

    v = variants[0]
    offset = v.start - r1.reference_start
    ref_len = len(v.ref)
    extended = r1_seq[:offset] + v.alt + r1_seq[offset + ref_len :]
    expected_trimmed = {extended[i : i + read_len] for i in range(0, 3)}

    mate1_seqs: set[str] = set()
    mate2_seqs: set[str] = set()
    mate1_totals: set[int] = set()
    mate2_totals: set[int] = set()

    for fq in (out_r1, out_r2):
        for name, seq, qual in iter_fastq(fq):
            if name.split("_WASP_", 1)[0] != "pairI":
                continue
            if len(seq) != read_len or len(qual) != read_len:
                raise AssertionError(
                    f"Length mismatch for {name}: seq={len(seq)} qual={len(qual)} expected={read_len}"
                )
            if name.endswith("/1"):
                mate1_seqs.add(seq)
                mate1_totals.add(_parse_total_seqs_from_name(name))
            else:
                mate2_seqs.add(seq)
                mate2_totals.add(_parse_total_seqs_from_name(name))

    assert mate1_seqs == expected_trimmed
    assert mate2_seqs == {r2_seq}
    assert mate1_totals == {3}
    assert mate2_totals == {3}

