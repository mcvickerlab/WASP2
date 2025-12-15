#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import shutil
import sys
import tempfile
import time
from pathlib import Path


def _add_repo_src_to_path() -> Path:
    root = Path(__file__).resolve().parents[2]
    src = root / "src"
    if str(src) not in sys.path:
        sys.path.insert(0, str(src))
    return root


def run_snv_parity(out_dir: Path, keep_tmp: bool) -> int:
    root = _add_repo_src_to_path()

    from benchmarking.quickbench.fastq_utils import counter_diff, fastq_counter
    from benchmarking.quickbench.synthetic_dataset import (
        quickbench_snv_variants,
        write_bed,
        write_synthetic_bam,
    )
    from mapping.intersect_variant_data import intersect_reads, process_bam
    from mapping.make_remap_reads import write_remap_bam

    try:
        import wasp2_rust
    except ImportError as e:
        print(f"ERROR: wasp2_rust not importable: {e}", file=sys.stderr)
        return 2

    out_dir.mkdir(parents=True, exist_ok=True)

    work_dir = Path(
        tempfile.mkdtemp(
            prefix="wasp2_quickbench_", dir=str(out_dir) if keep_tmp else None
        )
    )

    try:
        bam = work_dir / "synthetic.bam"
        bed = work_dir / "variants_snv.bed"
        write_synthetic_bam(bam)
        write_bed(bed, quickbench_snv_variants())

        # -----------------------------
        # Baseline: multi-pass pipeline
        # -----------------------------
        baseline_dir = work_dir / "baseline"
        baseline_dir.mkdir()
        to_remap_bam = baseline_dir / "to_remap.bam"
        keep_bam = baseline_dir / "keep.bam"
        remap_reads_txt = baseline_dir / "remap_reads.txt"
        intersect_bed = baseline_dir / "intersect.bed"
        baseline_r1 = baseline_dir / "baseline_r1.fq"
        baseline_r2 = baseline_dir / "baseline_r2.fq"

        t0 = time.perf_counter()
        process_bam(
            bam_file=str(bam),
            vcf_bed=str(bed),
            remap_bam=str(to_remap_bam),
            remap_reads=str(remap_reads_txt),
            keep_bam=str(keep_bam),
            is_paired=True,
            threads=1,
        )
        intersect_reads(
            remap_bam=str(to_remap_bam),
            vcf_bed=str(bed),
            out_bed=str(intersect_bed),
            num_samples=1,
        )
        write_remap_bam(
            bam_file=str(to_remap_bam),
            intersect_file=str(intersect_bed),
            r1_out=str(baseline_r1),
            r2_out=str(baseline_r2),
            samples=["SYNTH"],
            max_seqs=64,
            include_indels=False,
        )
        baseline_s = time.perf_counter() - t0

        # -----------------------
        # Candidate: unified Rust
        # -----------------------
        unified_dir = work_dir / "unified"
        unified_dir.mkdir()
        unified_r1 = unified_dir / "unified_r1.fq"
        unified_r2 = unified_dir / "unified_r2.fq"
        unified_remap_names = unified_dir / "remap_names.txt"
        unified_keep_no_flip_names = unified_dir / "keep_no_flip_names.txt"

        t1 = time.perf_counter()
        stats = wasp2_rust.unified_make_reads_py(
            str(bam),
            str(bed),
            str(unified_r1),
            str(unified_r2),
            max_seqs=64,
            threads=1,
            compression_threads=1,
            compress_output=False,
            indel_mode=False,
            keep_no_flip_names_path=str(unified_keep_no_flip_names),
            remap_names_path=str(unified_remap_names),
        )
        unified_s = time.perf_counter() - t1

        # Compare FASTQ outputs (as multisets, order-independent)
        baseline_counter = fastq_counter(baseline_r1, baseline_r2)
        unified_counter = fastq_counter(unified_r1, unified_r2)

        only_baseline, only_unified = counter_diff(baseline_counter, unified_counter)

        summary = {
            "work_dir": str(work_dir),
            "bam": str(bam),
            "bed": str(bed),
            "baseline": {
                "r1": str(baseline_r1),
                "r2": str(baseline_r2),
                "seconds": baseline_s,
                "records_total": sum(baseline_counter.values()),
            },
            "unified": {
                "r1": str(unified_r1),
                "r2": str(unified_r2),
                "seconds": unified_s,
                "records_total": sum(unified_counter.values()),
                "stats": stats,
                "remap_names_path": str(unified_remap_names),
                "keep_no_flip_names_path": str(unified_keep_no_flip_names),
            },
            "match": len(only_baseline) == 0 and len(only_unified) == 0,
        }

        out_json = out_dir / "quickbench_snv_parity.json"
        out_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

        if not summary["match"]:
            print("SNV parity mismatch:")
            print(f"- Only in baseline: {len(only_baseline)} canonical records")
            print(f"- Only in unified:  {len(only_unified)} canonical records")
            for rec, n in only_baseline[:10]:
                print(f"  baseline +{n}: {rec}")
            for rec, n in only_unified[:10]:
                print(f"  unified +{n}: {rec}")
            print(f"Wrote: {out_json}")
            return 1

        print(f"✅ SNV parity OK (wrote {out_json})")
        return 0
    finally:
        if keep_tmp:
            print(f"Keeping quickbench dir: {work_dir}")
        else:
            shutil.rmtree(work_dir, ignore_errors=True)


def _parse_total_seqs_from_name(name: str) -> int:
    if name.endswith("/1") or name.endswith("/2"):
        core = name[:-2]
    else:
        core = name
    try:
        suffix = core.split("_WASP_", 1)[1]
        parts = suffix.split("_")
        return int(parts[3])
    except Exception as e:  # pragma: no cover
        raise ValueError(f"Failed to parse total_seqs from FASTQ name: {name}") from e


def run_indel_parity(out_dir: Path, keep_tmp: bool) -> int:
    _add_repo_src_to_path()

    from benchmarking.quickbench.fastq_utils import counter_diff, fastq_counter
    from benchmarking.quickbench.synthetic_dataset import (
        quickbench_indel_variants,
        write_bed,
        write_synthetic_bam_indel,
    )
    from mapping.intersect_variant_data import intersect_reads, process_bam
    from mapping.make_remap_reads import write_remap_bam

    try:
        import wasp2_rust
    except ImportError as e:
        print(f"ERROR: wasp2_rust not importable: {e}", file=sys.stderr)
        return 2

    out_dir.mkdir(parents=True, exist_ok=True)
    work_dir = Path(
        tempfile.mkdtemp(
            prefix="wasp2_quickbench_", dir=str(out_dir) if keep_tmp else None
        )
    )

    try:
        bam = work_dir / "synthetic_indel.bam"
        bed = work_dir / "variants_indel.bed"
        write_synthetic_bam_indel(bam)
        write_bed(bed, quickbench_indel_variants())

        baseline_dir = work_dir / "baseline"
        baseline_dir.mkdir()
        to_remap_bam = baseline_dir / "to_remap.bam"
        keep_bam = baseline_dir / "keep.bam"
        remap_reads_txt = baseline_dir / "remap_reads.txt"
        intersect_bed = baseline_dir / "intersect.bed"
        baseline_r1 = baseline_dir / "baseline_r1.fq"
        baseline_r2 = baseline_dir / "baseline_r2.fq"

        t0 = time.perf_counter()
        process_bam(
            bam_file=str(bam),
            vcf_bed=str(bed),
            remap_bam=str(to_remap_bam),
            remap_reads=str(remap_reads_txt),
            keep_bam=str(keep_bam),
            is_paired=True,
            threads=1,
        )
        intersect_reads(
            remap_bam=str(to_remap_bam),
            vcf_bed=str(bed),
            out_bed=str(intersect_bed),
            num_samples=1,
        )
        write_remap_bam(
            bam_file=str(to_remap_bam),
            intersect_file=str(intersect_bed),
            r1_out=str(baseline_r1),
            r2_out=str(baseline_r2),
            samples=["SYNTH"],
            max_seqs=64,
            include_indels=True,
        )
        baseline_s = time.perf_counter() - t0

        unified_dir = work_dir / "unified"
        unified_dir.mkdir()
        unified_r1 = unified_dir / "unified_r1.fq"
        unified_r2 = unified_dir / "unified_r2.fq"

        t1 = time.perf_counter()
        stats = wasp2_rust.unified_make_reads_py(
            str(bam),
            str(bed),
            str(unified_r1),
            str(unified_r2),
            max_seqs=64,
            threads=1,
            compression_threads=1,
            compress_output=False,
            indel_mode=False,  # parity to baseline (no trim combos)
        )
        unified_s = time.perf_counter() - t1

        baseline_counter = fastq_counter(baseline_r1, baseline_r2)
        unified_counter = fastq_counter(unified_r1, unified_r2)
        only_baseline, only_unified = counter_diff(baseline_counter, unified_counter)

        summary = {
            "work_dir": str(work_dir),
            "bam": str(bam),
            "bed": str(bed),
            "baseline": {
                "r1": str(baseline_r1),
                "r2": str(baseline_r2),
                "seconds": baseline_s,
                "records_total": sum(baseline_counter.values()),
            },
            "unified": {
                "r1": str(unified_r1),
                "r2": str(unified_r2),
                "seconds": unified_s,
                "records_total": sum(unified_counter.values()),
                "stats": stats,
            },
            "match": len(only_baseline) == 0 and len(only_unified) == 0,
        }

        out_json = out_dir / "quickbench_indel_parity.json"
        out_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

        if not summary["match"]:
            print("INDEL parity mismatch:")
            print(f"- Only in baseline: {len(only_baseline)} canonical records")
            print(f"- Only in unified:  {len(only_unified)} canonical records")
            print(f"Wrote: {out_json}")
            return 1

        print(f"✅ INDEL parity OK (wrote {out_json})")
        return 0
    finally:
        if keep_tmp:
            print(f"Keeping quickbench dir: {work_dir}")
        else:
            shutil.rmtree(work_dir, ignore_errors=True)


def run_indel_trim_invariants(out_dir: Path, keep_tmp: bool) -> int:
    _add_repo_src_to_path()

    from benchmarking.quickbench.fastq_utils import iter_fastq
    from benchmarking.quickbench.synthetic_dataset import (
        quickbench_indel_variants,
        write_bed,
        write_synthetic_bam_indel,
    )

    try:
        import wasp2_rust
    except ImportError as e:
        print(f"ERROR: wasp2_rust not importable: {e}", file=sys.stderr)
        return 2

    import pysam

    out_dir.mkdir(parents=True, exist_ok=True)
    work_dir = Path(
        tempfile.mkdtemp(
            prefix="wasp2_quickbench_", dir=str(out_dir) if keep_tmp else None
        )
    )

    try:
        bam = work_dir / "synthetic_indel.bam"
        bed = work_dir / "variants_indel.bed"
        write_synthetic_bam_indel(bam)
        variants = quickbench_indel_variants()
        write_bed(bed, variants)

        unified_dir = work_dir / "unified_indel_trim"
        unified_dir.mkdir()
        unified_r1 = unified_dir / "unified_r1.fq"
        unified_r2 = unified_dir / "unified_r2.fq"

        stats = wasp2_rust.unified_make_reads_py(
            str(bam),
            str(bed),
            str(unified_r1),
            str(unified_r2),
            max_seqs=256,
            threads=1,
            compression_threads=1,
            compress_output=False,
            indel_mode=True,
            max_indel_size=50,
        )

        # Pull original sequences for pairI from BAM
        with pysam.AlignmentFile(str(bam), "rb") as bf:
            recs = [r for r in bf.fetch(until_eof=True) if r.query_name == "pairI"]
        r1 = next(r for r in recs if r.is_read1)
        r2 = next(r for r in recs if r.is_read2)
        r1_seq = r1.query_sequence
        r2_seq = r2.query_sequence

        # Compute expected trimmed sequences for the +2bp insertion haplotype.
        v = variants[0]
        offset = v.start - r1.reference_start
        assert offset >= 0
        ref_len = len(v.ref)
        extended = r1_seq[:offset] + v.alt + r1_seq[offset + ref_len :]
        expected_trimmed = {extended[i : i + len(r1_seq)] for i in range(0, 3)}

        # Validate FASTQs
        per_mate_seqs: dict[tuple[str, int], set[str]] = {}
        per_mate_totals: dict[tuple[str, int], set[int]] = {}
        for fq in (unified_r1, unified_r2):
            for name, seq, qual in iter_fastq(fq):
                if name.split("_WASP_", 1)[0] != "pairI":
                    continue
                mate = 1 if name.endswith("/1") else 2
                per_mate_seqs.setdefault(("pairI", mate), set()).add(seq)
                per_mate_totals.setdefault(("pairI", mate), set()).add(
                    _parse_total_seqs_from_name(name)
                )
                if len(seq) != len(r1_seq) or len(qual) != len(r1_seq):
                    raise AssertionError(
                        f"Length mismatch for {name}: seq={len(seq)} qual={len(qual)} expected={len(r1_seq)}"
                    )

        r1_set = per_mate_seqs.get(("pairI", 1), set())
        r2_set = per_mate_seqs.get(("pairI", 2), set())
        if r1_set != expected_trimmed:
            raise AssertionError(
                f"Unexpected trimmed R1 sequences for pairI.\nExpected: {sorted(expected_trimmed)}\nGot: {sorted(r1_set)}"
            )
        if r2_set != {r2_seq}:
            raise AssertionError(
                f"Unexpected R2 sequences for pairI. Expected original only: {r2_seq}"
            )

        if per_mate_totals.get(("pairI", 1), set()) != {3}:
            raise AssertionError(f"pairI/1 total_seqs != 3: {per_mate_totals.get(('pairI', 1))}")
        if per_mate_totals.get(("pairI", 2), set()) != {3}:
            raise AssertionError(f"pairI/2 total_seqs != 3: {per_mate_totals.get(('pairI', 2))}")

        summary = {
            "work_dir": str(work_dir),
            "bam": str(bam),
            "bed": str(bed),
            "unified": {"r1": str(unified_r1), "r2": str(unified_r2), "stats": stats},
            "expected_trimmed_count": len(expected_trimmed),
        }
        out_json = out_dir / "quickbench_indel_trim_invariants.json"
        out_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

        print(f"✅ INDEL trim invariants OK (wrote {out_json})")
        return 0
    finally:
        if keep_tmp:
            print(f"Keeping quickbench dir: {work_dir}")
        else:
            shutil.rmtree(work_dir, ignore_errors=True)


def main() -> int:
    p = argparse.ArgumentParser(description="WASP2 quickbench (fast sanity + parity checks)")
    p.add_argument(
        "--out-dir",
        type=Path,
        default=Path("benchmarking/quickbench_out"),
        help="Where to write quickbench outputs (default: benchmarking/quickbench_out).",
    )
    p.add_argument(
        "--keep-tmp",
        action="store_true",
        help="Keep the temporary working directory for inspection.",
    )
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("snv-parity", help="Compare unified output vs multi-pass output on a synthetic SNV dataset.")
    sub.add_parser(
        "indel-parity",
        help="Compare unified output vs multi-pass output on a synthetic INDEL dataset (no trim combos).",
    )
    sub.add_parser(
        "indel-trim-invariants",
        help="Validate trim-combo invariants on a synthetic +2bp insertion dataset.",
    )

    args = p.parse_args()
    if args.cmd == "snv-parity":
        return run_snv_parity(args.out_dir, args.keep_tmp)
    if args.cmd == "indel-parity":
        return run_indel_parity(args.out_dir, args.keep_tmp)
    if args.cmd == "indel-trim-invariants":
        return run_indel_trim_invariants(args.out_dir, args.keep_tmp)
    raise AssertionError(f"Unhandled command: {args.cmd}")


if __name__ == "__main__":
    raise SystemExit(main())
