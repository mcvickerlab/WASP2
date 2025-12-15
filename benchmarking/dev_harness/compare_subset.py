#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple


@dataclass(frozen=True)
class Region:
    chrom: str
    start0: Optional[int] = None  # 0-based inclusive
    end0: Optional[int] = None  # 0-based exclusive


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _ensure_import_paths() -> None:
    root = _repo_root()
    src = root / "src"
    for p in (root, src):
        if str(p) not in sys.path:
            sys.path.insert(0, str(p))


def _parse_region(region: str) -> Region:
    # Accept:
    # - chr1
    # - chr1:1000-2000 (samtools 1-based inclusive coordinates)
    if ":" not in region:
        return Region(chrom=region)
    chrom, rest = region.split(":", 1)
    if "-" not in rest:
        raise ValueError(f"Invalid region (expected chr:start-end): {region}")
    start1_s, end1_s = rest.split("-", 1)
    start1 = int(start1_s.replace(",", ""))
    end1 = int(end1_s.replace(",", ""))
    if start1 <= 0 or end1 <= 0 or end1 < start1:
        raise ValueError(f"Invalid region coordinates: {region}")
    # Convert to 0-based half-open for BED filtering.
    start0 = start1 - 1
    end0 = end1
    return Region(chrom=chrom, start0=start0, end0=end0)


def _subset_bam(in_bam: Path, out_bam: Path, region: str) -> None:
    out_bam.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_bam.with_suffix(".tmp.bam")
    subprocess.run(["samtools", "view", "-b", "-h", str(in_bam), region, "-o", str(tmp)], check=True)
    subprocess.run(["samtools", "sort", "-o", str(out_bam), str(tmp)], check=True)
    subprocess.run(["samtools", "index", str(out_bam)], check=True)
    tmp.unlink(missing_ok=True)


def _subset_bed(in_bed: Path, out_bed: Path, region: Region) -> None:
    out_bed.parent.mkdir(parents=True, exist_ok=True)
    with in_bed.open("r", encoding="utf-8") as fin, out_bed.open("w", encoding="utf-8") as fout:
        for line in fin:
            if not line or line.startswith("#") or line.strip() == "":
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            if chrom != region.chrom:
                continue
            if region.start0 is not None and region.end0 is not None:
                try:
                    start = int(parts[1])
                    stop = int(parts[2])
                except ValueError:
                    continue
                # Keep variants that overlap the requested interval.
                if stop <= region.start0 or start >= region.end0:
                    continue
            fout.write(line)


def main() -> int:
    _ensure_import_paths()

    p = argparse.ArgumentParser(
        description="Compare unified vs multi-pass remap FASTQs on a real BAM/BED subset."
    )
    p.add_argument("--bam", type=Path, required=True, help="Coordinate-sorted input BAM")
    p.add_argument("--bed", type=Path, required=True, help="Variant BED (chrom start stop ref alt GT)")
    p.add_argument(
        "--region",
        type=str,
        required=True,
        help="Subset region (samtools syntax), e.g. chr1:1,000,000-2,000,000",
    )
    p.add_argument(
        "--out-dir",
        type=Path,
        default=Path("benchmarking/dev_harness_out"),
        help="Output directory for artifacts and summary JSON",
    )
    p.add_argument(
        "--threads",
        type=int,
        default=8,
        help="Threads for unified pipeline (baseline multi-pass uses internal defaults)",
    )
    p.add_argument(
        "--indel-mode",
        action="store_true",
        help="Run unified in indel_mode=True (trim combos). Parity compare is disabled in this mode.",
    )
    p.add_argument(
        "--keep-tmp",
        action="store_true",
        help="Keep the temporary subset working directory for inspection.",
    )
    args = p.parse_args()

    import wasp2_rust

    from benchmarking.quickbench.fastq_utils import counter_diff, fastq_counter
    from mapping.intersect_variant_data import intersect_reads, process_bam
    from mapping.make_remap_reads import write_remap_bam

    region = _parse_region(args.region)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    tmp_parent: Optional[str] = str(args.out_dir) if args.keep_tmp else None
    work_dir = Path(tempfile.mkdtemp(prefix="wasp2_dev_harness_", dir=tmp_parent))
    try:
        subset_bam = work_dir / "subset.bam"
        subset_bed = work_dir / "subset.bed"

        t_subset = time.perf_counter()
        _subset_bam(args.bam, subset_bam, args.region)
        _subset_bed(args.bed, subset_bed, region)
        subset_s = time.perf_counter() - t_subset

        unified_dir = work_dir / "unified"
        unified_dir.mkdir()
        unified_r1 = unified_dir / "unified_r1.fq"
        unified_r2 = unified_dir / "unified_r2.fq"

        t_unified = time.perf_counter()
        unified_stats = wasp2_rust.unified_make_reads_parallel_py(
            str(subset_bam),
            str(subset_bed),
            str(unified_r1),
            str(unified_r2),
            max_seqs=256 if args.indel_mode else 64,
            threads=max(1, args.threads),
            compression_threads=1,
            compress_output=False,
            indel_mode=bool(args.indel_mode),
        )
        unified_s = time.perf_counter() - t_unified

        summary: dict = {
            "input": {"bam": str(args.bam), "bed": str(args.bed), "region": args.region},
            "subset": {"bam": str(subset_bam), "bed": str(subset_bed), "seconds": subset_s},
            "unified": {
                "seconds": unified_s,
                "r1": str(unified_r1),
                "r2": str(unified_r2),
                "stats": unified_stats,
                "indel_mode": bool(args.indel_mode),
            },
        }

        if args.indel_mode:
            summary["compare"] = {
                "enabled": False,
                "reason": "indel_mode=True generates trim combos; baseline multi-pass does not.",
            }
        else:
            baseline_dir = work_dir / "baseline"
            baseline_dir.mkdir()
            to_remap_bam = baseline_dir / "to_remap.bam"
            keep_bam = baseline_dir / "keep.bam"
            remap_reads_txt = baseline_dir / "remap_reads.txt"
            intersect_bed = baseline_dir / "intersect.bed"
            baseline_r1 = baseline_dir / "baseline_r1.fq"
            baseline_r2 = baseline_dir / "baseline_r2.fq"

            t_base = time.perf_counter()
            process_bam(
                bam_file=str(subset_bam),
                vcf_bed=str(subset_bed),
                remap_bam=str(to_remap_bam),
                remap_reads=str(remap_reads_txt),
                keep_bam=str(keep_bam),
                is_paired=True,
                threads=1,
            )
            intersect_reads(
                remap_bam=str(to_remap_bam),
                vcf_bed=str(subset_bed),
                out_bed=str(intersect_bed),
                num_samples=1,
            )
            write_remap_bam(
                bam_file=str(to_remap_bam),
                intersect_file=str(intersect_bed),
                r1_out=str(baseline_r1),
                r2_out=str(baseline_r2),
                samples=["SUBSET"],
                max_seqs=64,
                include_indels=True,
            )
            baseline_s = time.perf_counter() - t_base

            baseline_counter = fastq_counter(baseline_r1, baseline_r2)
            unified_counter = fastq_counter(unified_r1, unified_r2)
            only_baseline, only_unified = counter_diff(baseline_counter, unified_counter)

            summary["baseline"] = {
                "seconds": baseline_s,
                "r1": str(baseline_r1),
                "r2": str(baseline_r2),
                "records_total": sum(baseline_counter.values()),
            }
            summary["compare"] = {
                "enabled": True,
                "match": len(only_baseline) == 0 and len(only_unified) == 0,
                "only_in_baseline": len(only_baseline),
                "only_in_unified": len(only_unified),
            }

        out_json = args.out_dir / "subset_compare_summary.json"
        out_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
        print(f"Wrote: {out_json}")

        if summary.get("compare", {}).get("enabled") and not summary["compare"]["match"]:
            return 1
        return 0
    finally:
        if args.keep_tmp:
            print(f"Keeping work dir: {work_dir}")
        else:
            shutil.rmtree(work_dir, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
