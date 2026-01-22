#!/usr/bin/env python3
"""
Generate Figure 3 Panel D data (gene imprinting) from the existing validation pipeline.

Outputs:
  paper/figure3/data/gene_significance.tsv

This wraps:
  validation/gene_imprinting/run_atac_analysis.py          -> results/atac_analysis/variant_counts.tsv
  validation/gene_imprinting/run_significance_analysis.py -> results/atac_analysis/gene_significance.tsv
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
ATAC_RESULTS_DIR = REPO_ROOT / "results" / "atac_analysis"
OUT_TSV = REPO_ROOT / "paper" / "figure3" / "data" / "gene_significance.tsv"


def _run(cmd: list[str]) -> None:
    result = subprocess.run(cmd, cwd=REPO_ROOT)
    if result.returncode != 0:
        raise SystemExit(result.returncode)


def main() -> None:
    out_dir = OUT_TSV.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Generating imprinting data for Figure 3 Panel D…", flush=True)
    _run([sys.executable, "validation/gene_imprinting/run_atac_analysis.py"])
    _run([sys.executable, "validation/gene_imprinting/run_significance_analysis.py"])

    src = ATAC_RESULTS_DIR / "gene_significance.tsv"
    if not src.exists():
        raise FileNotFoundError(f"Expected file not found: {src}")

    shutil.copy2(src, OUT_TSV)
    print(f"Saved: {OUT_TSV}", flush=True)


if __name__ == "__main__":
    main()

