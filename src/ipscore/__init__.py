"""
iPSCORE Multi-Tissue Allelic Imbalance Resource.

This module provides utilities for:
- Data inventory and validation across iPSCORE tissues (CVPC, PPC, iPSC)
- Sample manifest generation and harmonization
- QTL data loading and fine-mapping integration
- ML-ready output generation for GenVarLoader

Supports analysis of:
- CVPC RNA-seq (137 samples)
- CVPC ATAC-seq (137 samples)
- PPC RNA-seq (106 samples)
- PPC ATAC-seq (108 samples)
- iPSC RNA-seq (220 samples)

Total: 463 RNA samples + 245 ATAC samples = 708 sample-assays
"""

__version__ = "0.1.0"

from .data_inventory import (
    DataInventory,
    validate_inventory,
)
from .qtl_loader import (
    QTLLoader,
    create_qtl_loader,
    load_all_caqtls,
    load_finemapped_qtls,
)
from .sample_manifest import (
    SampleManifest,
    create_unified_manifest,
)

__all__ = [
    "DataInventory",
    "validate_inventory",
    "SampleManifest",
    "create_unified_manifest",
    "QTLLoader",
    "create_qtl_loader",
    "load_all_caqtls",
    "load_finemapped_qtls",
]
