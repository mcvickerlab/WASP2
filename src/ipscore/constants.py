"""
iPSCORE data path constants and configuration.

All paths reference the Frazer lab data infrastructure.
Environment variables can override defaults for testing or alternate deployments:
- IPSCORE_BASE: Base path for iPSCORE data
- WASP2_BASE: Base path for WASP2 evaluation data
"""

import os
from pathlib import Path
from typing import TypedDict


class DatasetConfig(TypedDict):
    """Configuration for a single iPSCORE dataset."""

    tissue: str
    assay: str
    expected_samples: int
    counts_path: Path
    master_csv: Path | None


# Base paths for iPSCORE data (with environment variable overrides)
IPSCORE_BASE = Path(os.environ.get("IPSCORE_BASE", "/iblm/netapp/data4/Frazer_collab"))
WASP2_BASE = Path(os.environ.get("WASP2_BASE", "/iblm/netapp/data3/jjaureguy/gvl_files/wasp2"))

# QTL summary statistics location
QTL_DATA_PATH = IPSCORE_BASE / "seqdata_datasets/carter_collab/ipscore/iPSCORE_files"

# Sample master CSV locations
SAMPLE_MANIFEST_PATHS = {
    "CVPC": IPSCORE_BASE
    / "ipscs/KF_transfer_dirs/KF_transfer_material/cvpc/final_si_dfs/CVPC_master.csv",
    "PPC": IPSCORE_BASE / "ipscs/KF_transfer_dirs/KF_transfer/final_si_dfs/PPC_master.csv",
    "iPSC": IPSCORE_BASE
    / "ipscs/KF_transfer_dirs/KF_transfer/final_si_dfs/GSE203377_ipsc_master.csv",
}

# WASP allelic counts paths for each dataset
WASP_COUNTS_PATHS = {
    "CVPC_RNA": WASP2_BASE
    / "WASP2_extensive_evaluation/WASP2_current/cvpc/data_counts/snv_only_genome/rna",
    "CVPC_ATAC": WASP2_BASE
    / "WASP2_extensive_evaluation/WASP2_current/cvpc/data_counts/snv_only_genome/atac",
    "PPC_RNA": Path(
        "/iblm/netapp/data4/shared_dir/hyena_dna_collab/downstream_tasks/iPSCORE_ppc/wasp_counts"
    ),
    "PPC_ATAC": IPSCORE_BASE / "ipscs/datasets/processed/ipscs_PPC/wasp/wasp_asoc_counts",
    "iPSC_RNA": Path(
        "/iblm/netapp/data4/shared_dir/hyena_dna_collab/downstream_tasks/iPSCORE_ipscs/wasp_counts"
    ),
}

# Dataset configurations with expected sample counts
DATASETS: dict[str, DatasetConfig] = {
    "CVPC_RNA": {
        "tissue": "CVPC",
        "assay": "RNA",
        "expected_samples": 137,
        "counts_path": WASP_COUNTS_PATHS["CVPC_RNA"],
        "master_csv": SAMPLE_MANIFEST_PATHS["CVPC"],
    },
    "CVPC_ATAC": {
        "tissue": "CVPC",
        "assay": "ATAC",
        "expected_samples": 137,
        "counts_path": WASP_COUNTS_PATHS["CVPC_ATAC"],
        "master_csv": SAMPLE_MANIFEST_PATHS["CVPC"],
    },
    "PPC_RNA": {
        "tissue": "PPC",
        "assay": "RNA",
        "expected_samples": 106,
        "counts_path": WASP_COUNTS_PATHS["PPC_RNA"],
        "master_csv": SAMPLE_MANIFEST_PATHS["PPC"],
    },
    "PPC_ATAC": {
        "tissue": "PPC",
        "assay": "ATAC",
        "expected_samples": 108,
        "counts_path": WASP_COUNTS_PATHS["PPC_ATAC"],
        "master_csv": SAMPLE_MANIFEST_PATHS["PPC"],
    },
    "iPSC_RNA": {
        "tissue": "iPSC",
        "assay": "RNA",
        "expected_samples": 220,
        "counts_path": WASP_COUNTS_PATHS["iPSC_RNA"],
        "master_csv": SAMPLE_MANIFEST_PATHS["iPSC"],
    },
}

# QTL file paths
QTL_FILES = {
    "all_caqtls": QTL_DATA_PATH / "all_caqtls.txt",
    "ipsc_finemapped": QTL_DATA_PATH / "ipsc_caqtl_finemapped.txt.gz",
    "ipsc_sumstats": QTL_DATA_PATH / "ipsc_caqtl_sumstats.txt.gz",
    "cvpc_caqtls_stats": QTL_DATA_PATH / "CVPC_downstream_caqtls_stats.txt.gz",
}

# QTL counts per tissue (from Issue #40)
QTL_COUNTS = {
    "CVPC": 13249,
    "PPC": 12236,
    "iPSC": 11074,
}

# Tissue labels used in all_caqtls.txt
TISSUE_LABELS = ["CVPC", "PPC", "iPSC"]
