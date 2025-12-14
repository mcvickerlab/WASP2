#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


KEYS = [
    "total_reads",
    "pairs_processed",
    "pairs_with_variants",
    "pairs_with_snvs_only",
    "pairs_with_indels_only",
    "pairs_with_snvs_and_indels",
    "pairs_kept",
    "pairs_keep_no_flip",
    "pairs_skipped_unmappable",
    "pairs_haplotype_failed",
    "orphan_reads",
    "haplotypes_written",
    "tree_build_ms",
    "bam_stream_ms",
    "overlap_query_ms",
    "pair_process_ms",
    "send_ms",
    "writer_thread_ms",
    "wall_s",
]


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Log a concise one-line summary from a WASP2 unified_stats.json"
    )
    parser.add_argument("stats_json", type=Path, help="Path to unified_stats.json")
    parser.add_argument(
        "--label",
        default="",
        help="Optional label (e.g., STEP2, atac_indel, rnaseq_snv)",
    )
    parser.add_argument(
        "--prefix",
        default="UNIFIED_STATS",
        help="Line prefix for greppable logs (default: UNIFIED_STATS)",
    )
    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    if not args.stats_json.exists():
        print(f"ERROR: stats JSON not found: {args.stats_json}", file=sys.stderr)
        return 2

    try:
        data = json.loads(args.stats_json.read_text())
    except Exception as exc:  # noqa: BLE001
        print(f"ERROR: failed to parse JSON: {args.stats_json}: {exc}", file=sys.stderr)
        return 2

    parts: list[str] = []
    if args.prefix:
        parts.append(args.prefix)
    if args.label:
        parts.append(f"label={args.label}")
    parts.append(f"stats_json={args.stats_json}")

    for key in KEYS:
        if key in data:
            value = data[key]
        else:
            value = "NA" if key == "wall_s" else 0
        parts.append(f"{key}={value}")

    print(" ".join(parts))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
