"""
WASP2: Allele-Specific Pipeline, Version 2.

A Python package for allele-specific analysis of sequencing data.
"""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("wasp2")
except PackageNotFoundError:
    __version__ = "0+unknown"
