#!/usr/bin/env python3
"""
Combine per-task unified SNV scaling TSVs into a single COMBINED.tsv.

Each array task in `benchmark_unified_scaling.sh` writes a TSV to:
  benchmarking/results/unified_scaling_<timestamp>/wasp2_unified_scaling.tsv

By default this script keeps only the latest run per n_reads (so older reruns
don't pollute the plot). Use --all-runs to include everything.
"""

from __future__ import annotations

import argparse
import datetime as dt
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def parse_run_time_from_dir(dir_name: str) -> Optional[dt.datetime]:
    prefix = "unified_scaling_"
    if not dir_name.startswith(prefix):
        return None
    ts = dir_name[len(prefix) :]
    for fmt in ("%Y-%m-%d_%H-%M-%S", "%Y-%m-%d_%H-%M"):
        try:
            return dt.datetime.strptime(ts, fmt)
        except ValueError:
            continue
    return None


def read_records(tsv_path: Path, run_time: dt.datetime) -> List[dict]:
    records: List[dict] = []
    with tsv_path.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("timestamp\t") or line.startswith("timestamp"):
                continue
            parts = line.split("\t")
            if len(parts) < 7:
                continue
            timestamp, n_reads, seed, total_s, unified_s, remap_s, filter_s = parts[:7]
            try:
                n_reads_i = int(float(n_reads))
                seed_i = int(float(seed))
                total_f = float(total_s)
                unified_f = float(unified_s)
                remap_f = float(remap_s)
                filter_f = float(filter_s)
            except ValueError:
                continue
            records.append(
                {
                    "timestamp": timestamp,
                    "n_reads": n_reads_i,
                    "seed": seed_i,
                    "total_s": total_f,
                    "unified_s": unified_f,
                    "remap_s": remap_f,
                    "filter_s": filter_f,
                    "run_time": run_time,
                    "path": str(tsv_path),
                }
            )
    return records


def combine(pattern: str, out_path: Path, all_runs: bool) -> Tuple[List[dict], List[Path]]:
    files = sorted(Path(".").glob(pattern))
    records: List[dict] = []
    for p in files:
        run_time = parse_run_time_from_dir(p.parent.name)
        if run_time is None:
            run_time = dt.datetime.fromtimestamp(p.stat().st_mtime)
        records.extend(read_records(p, run_time))

    if not records:
        return [], files

    if not all_runs:
        latest_by_n_reads: Dict[int, dt.datetime] = {}
        for r in records:
            n = r["n_reads"]
            t = r["run_time"]
            if n not in latest_by_n_reads or t > latest_by_n_reads[n]:
                latest_by_n_reads[n] = t
        records = [r for r in records if r["run_time"] == latest_by_n_reads[r["n_reads"]]]

    records.sort(key=lambda r: (r["n_reads"], r["seed"], r["timestamp"]))

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as out:
        out.write("timestamp\tn_reads\tseed\ttotal_s\tunified_s\tremap_s\tfilter_s\n")
        for r in records:
            out.write(
                f"{r['timestamp']}\t{r['n_reads']}\t{r['seed']}\t"
                f"{r['total_s']}\t{r['unified_s']}\t{r['remap_s']}\t{r['filter_s']}\n"
            )
    return records, files


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pattern",
        default="benchmarking/results/unified_scaling_*/wasp2_unified_scaling.tsv",
        help="Glob pattern for per-task TSVs.",
    )
    parser.add_argument(
        "--out",
        default="benchmarking/results/wasp2_unified_scaling_COMBINED.tsv",
        help="Output combined TSV path.",
    )
    parser.add_argument(
        "--all-runs",
        action="store_true",
        help="Include all runs (do not de-duplicate by n_reads).",
    )
    args = parser.parse_args()

    records, files = combine(args.pattern, Path(args.out), args.all_runs)
    if not records:
        print(f"No records found. Looked for {len(files)} files matching {args.pattern}.")
        return
    n_unique = len({r["n_reads"] for r in records})
    print(f"Wrote {len(records)} records ({n_unique} unique n_reads) to {args.out}")


if __name__ == "__main__":
    main()

