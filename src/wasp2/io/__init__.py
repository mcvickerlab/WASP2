"""
I/O module for WASP2.

Provides data structures and readers for variant files (VCF, PGEN).
"""

from .variant_source import (
    Genotype,
    Variant,
    VariantGenotype,
    VariantSource,
)

# Import format handlers to register them with factory
from . import vcf_source  # noqa: F401

# Import PGEN handler if pgenlib is available
try:
    from . import pgen_source  # noqa: F401
except ImportError:
    pass  # pgenlib not available - PGEN support disabled

__all__ = [
    "Genotype",
    "Variant",
    "VariantGenotype",
    "VariantSource",
]
