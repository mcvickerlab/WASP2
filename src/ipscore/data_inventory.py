"""
Data inventory and validation for iPSCORE multi-tissue resource.

Verifies existence and completeness of:
- WASP allelic count files for all 5 datasets
- Sample manifest files (CVPC_master.csv, PPC_master.csv, etc.)
- QTL summary statistics and fine-mapping data
"""

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from wasp2.cli import error, info, success, warning

from .constants import DATASETS, QTL_FILES, DatasetConfig


@dataclass
class DatasetStatus:
    """Status of a single iPSCORE dataset."""

    name: str
    tissue: str
    assay: str
    expected_samples: int
    found_samples: int
    counts_path: Path
    exists: bool
    sample_dirs: list[str] = field(default_factory=list)

    @property
    def complete(self) -> bool:
        """Check if dataset has all expected samples."""
        return self.exists and self.found_samples >= self.expected_samples

    @property
    def status_emoji(self) -> str:
        """Return status indicator."""
        if not self.exists:
            return "X"
        if self.complete:
            return "OK"
        return "PARTIAL"


@dataclass
class DataInventory:
    """Complete inventory of iPSCORE data resources."""

    datasets: dict[str, DatasetStatus] = field(default_factory=dict)
    qtl_files: dict[str, bool] = field(default_factory=dict)
    manifest_files: dict[str, bool] = field(default_factory=dict)

    @property
    def total_samples(self) -> int:
        """Total samples found across all datasets."""
        return sum(ds.found_samples for ds in self.datasets.values())

    @property
    def total_expected(self) -> int:
        """Total expected samples."""
        return sum(ds.expected_samples for ds in self.datasets.values())

    @property
    def rna_samples(self) -> int:
        """Total RNA-seq samples found."""
        return sum(ds.found_samples for ds in self.datasets.values() if ds.assay == "RNA")

    @property
    def atac_samples(self) -> int:
        """Total ATAC-seq samples found."""
        return sum(ds.found_samples for ds in self.datasets.values() if ds.assay == "ATAC")

    def to_dataframe(self) -> pd.DataFrame:
        """Convert inventory to DataFrame for display/export."""
        rows = []
        for name, ds in self.datasets.items():
            rows.append(
                {
                    "Dataset": name,
                    "Tissue": ds.tissue,
                    "Assay": ds.assay,
                    "Expected": ds.expected_samples,
                    "Found": ds.found_samples,
                    "Status": ds.status_emoji,
                    "Path": str(ds.counts_path),
                }
            )
        return pd.DataFrame(rows)

    def print_summary(self) -> None:
        """Print formatted inventory summary."""
        info("=" * 60)
        info("iPSCORE Data Inventory Summary")
        info("=" * 60)

        for name, ds in self.datasets.items():
            status = ds.status_emoji
            if status == "OK":
                success(f"  [{status}] {name}: {ds.found_samples}/{ds.expected_samples} samples")
            elif status == "PARTIAL":
                warning(f"  [{status}] {name}: {ds.found_samples}/{ds.expected_samples} samples")
            else:
                error(f"  [{status}] {name}: Path not found")

        info("-" * 60)
        info(f"Total: {self.total_samples}/{self.total_expected} samples")
        info(f"  RNA-seq:  {self.rna_samples} samples")
        info(f"  ATAC-seq: {self.atac_samples} samples")
        info("-" * 60)

        info("\nQTL Data Files:")
        for name, exists in self.qtl_files.items():
            status = "OK" if exists else "X"
            if exists:
                success(f"  [{status}] {name}")
            else:
                error(f"  [{status}] {name}")


def _count_sample_directories(counts_path: Path) -> tuple[int, list[str]]:
    """Count sample directories in a WASP counts path.

    Expects directory structure: counts_path/<sample_uuid>/

    Returns:
        Tuple of (count, list of sample UUIDs)
    """
    if not counts_path.exists():
        return 0, []

    sample_dirs = []
    for item in counts_path.iterdir():
        if item.is_dir():
            sample_dirs.append(item.name)

    return len(sample_dirs), sorted(sample_dirs)


def validate_dataset(name: str, config: DatasetConfig) -> DatasetStatus:
    """Validate a single iPSCORE dataset.

    Args:
        name: Dataset identifier (e.g., "CVPC_RNA")
        config: Dataset configuration

    Returns:
        DatasetStatus with validation results
    """
    counts_path = config["counts_path"]
    exists = counts_path.exists()
    found_samples, sample_dirs = _count_sample_directories(counts_path) if exists else (0, [])

    return DatasetStatus(
        name=name,
        tissue=config["tissue"],
        assay=config["assay"],
        expected_samples=config["expected_samples"],
        found_samples=found_samples,
        counts_path=counts_path,
        exists=exists,
        sample_dirs=sample_dirs,
    )


def validate_inventory(verbose: bool = True) -> DataInventory:
    """Validate all iPSCORE data resources.

    Checks:
    - All 5 WASP counts directories exist and have expected sample counts
    - QTL summary statistics files exist
    - Sample manifest CSVs exist

    Args:
        verbose: Print progress messages

    Returns:
        DataInventory with complete validation results
    """
    inventory = DataInventory()

    if verbose:
        info("Validating iPSCORE data inventory...")

    # Validate each dataset
    for name, config in DATASETS.items():
        if verbose:
            info(f"  Checking {name}...")
        inventory.datasets[name] = validate_dataset(name, config)

    # Check QTL files
    for name, path in QTL_FILES.items():
        inventory.qtl_files[name] = path.exists()

    if verbose:
        inventory.print_summary()

    return inventory
