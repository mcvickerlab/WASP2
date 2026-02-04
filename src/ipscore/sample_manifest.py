"""
Unified sample manifest generation for iPSCORE multi-tissue resource.

Parses and harmonizes sample metadata from:
- CVPC_master.csv (137 samples)
- PPC_master.csv (106-108 samples)
- GSE203377_ipsc_master.csv (220 samples)

Creates a unified manifest linking:
- Sample UUIDs
- Tissue type (CVPC, PPC, iPSC)
- Assay type (RNA, ATAC)
- WASP counts file paths
- Genotype IDs (for VCF linkage)
"""

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from wasp2.cli import info, success, warning

from .constants import DATASETS, WASP_COUNTS_PATHS


@dataclass
class SampleManifest:
    """Unified sample manifest for iPSCORE data."""

    df: pd.DataFrame

    def __len__(self) -> int:
        return len(self.df)

    @property
    def tissues(self) -> list[str]:
        """Unique tissues in manifest."""
        return self.df["tissue"].unique().tolist()

    @property
    def assays(self) -> list[str]:
        """Unique assays in manifest."""
        return self.df["assay"].unique().tolist()

    def filter(
        self,
        tissue: str | None = None,
        assay: str | None = None,
    ) -> "SampleManifest":
        """Filter manifest by tissue and/or assay.

        Args:
            tissue: Filter to specific tissue (CVPC, PPC, iPSC)
            assay: Filter to specific assay (RNA, ATAC)

        Returns:
            Filtered SampleManifest
        """
        df = self.df.copy()
        if tissue:
            df = df[df["tissue"] == tissue]
        if assay:
            df = df[df["assay"] == assay]
        return SampleManifest(df=df)

    def to_csv(self, path: Path | str) -> None:
        """Export manifest to CSV."""
        self.df.to_csv(path, index=False)

    def to_tsv(self, path: Path | str) -> None:
        """Export manifest to TSV."""
        self.df.to_csv(path, sep="\t", index=False)

    def print_summary(self) -> None:
        """Print manifest summary statistics."""
        info("=" * 60)
        info("Unified Sample Manifest Summary")
        info("=" * 60)

        for tissue in self.tissues:
            tissue_df = self.df[self.df["tissue"] == tissue]
            info(f"\n{tissue}:")
            for assay in self.assays:
                count = len(tissue_df[tissue_df["assay"] == assay])
                if count > 0:
                    info(f"  {assay}: {count} samples")

        info("-" * 60)
        info(f"Total: {len(self.df)} sample-assay combinations")


def _load_master_csv(path: Path, tissue: str) -> pd.DataFrame:
    """Load and parse a master CSV file.

    Args:
        path: Path to master CSV
        tissue: Tissue label (CVPC, PPC, iPSC)

    Returns:
        DataFrame with standardized columns
    """
    if not path.exists():
        warning(f"Master CSV not found: {path}")
        return pd.DataFrame()

    df = pd.read_csv(path)

    # Extract the sample UUID column (first column is typically 'uuid')
    uuid_col = df.columns[0]
    samples_df = pd.DataFrame({"sample_id": df[uuid_col].astype(str)})
    samples_df["tissue"] = tissue

    # Extract genotype ID if available (typically 'best.1' or similar column)
    genotype_cols = [c for c in df.columns if "best" in c.lower()]
    if genotype_cols:
        samples_df["genotype_id"] = df[genotype_cols[0]].astype(str)
    else:
        samples_df["genotype_id"] = samples_df["sample_id"]

    return samples_df


def _find_sample_counts(
    sample_id: str,
    tissue: str,
    assay: str,
) -> Path | None:
    """Find WASP counts path for a specific sample.

    Args:
        sample_id: Sample UUID
        tissue: Tissue type
        assay: Assay type

    Returns:
        Path to counts directory if exists, None otherwise
    """
    dataset_key = f"{tissue}_{assay}"
    if dataset_key not in WASP_COUNTS_PATHS:
        return None

    counts_base = WASP_COUNTS_PATHS[dataset_key]
    sample_path = counts_base / sample_id

    if sample_path.exists():
        return sample_path
    return None


def _discover_samples_from_counts(
    tissue: str,
    assay: str,
) -> list[dict[str, str | Path]]:
    """Discover samples directly from WASP counts directories.

    When master CSVs don't have all samples, we can discover them
    from the actual count file directories.

    Args:
        tissue: Tissue type
        assay: Assay type

    Returns:
        List of sample records with paths
    """
    dataset_key = f"{tissue}_{assay}"
    if dataset_key not in WASP_COUNTS_PATHS:
        return []

    counts_base = WASP_COUNTS_PATHS[dataset_key]
    if not counts_base.exists():
        return []

    samples = []
    for sample_dir in counts_base.iterdir():
        if sample_dir.is_dir():
            samples.append(
                {
                    "sample_id": sample_dir.name,
                    "tissue": tissue,
                    "assay": assay,
                    "counts_path": str(sample_dir),
                    "genotype_id": sample_dir.name,  # Default to sample_id
                }
            )

    return samples


def create_unified_manifest(
    validate_paths: bool = True,
    discover_from_counts: bool = True,
    verbose: bool = True,
) -> SampleManifest:
    """Create unified sample manifest from all iPSCORE data sources.

    Args:
        validate_paths: Verify that counts paths exist
        discover_from_counts: Discover samples from counts dirs if not in master CSV
        verbose: Print progress messages

    Returns:
        SampleManifest with all samples across tissues and assays
    """
    all_records: list[dict[str, str | Path]] = []

    if verbose:
        info("Creating unified iPSCORE sample manifest...")

    # Process each dataset
    for dataset_name, config in DATASETS.items():
        tissue = config["tissue"]
        assay = config["assay"]

        if verbose:
            info(f"  Processing {dataset_name}...")

        if discover_from_counts:
            # Primary method: discover from actual counts directories
            samples = _discover_samples_from_counts(tissue, assay)

            if samples:
                all_records.extend(samples)
                if verbose:
                    success(f"    Found {len(samples)} samples from counts directory")
            else:
                warning(f"    No samples found for {dataset_name}")
        else:
            # Alternative: use master CSV and look up paths
            master_csv = config.get("master_csv")
            if master_csv and Path(master_csv).exists():
                master_df = _load_master_csv(Path(master_csv), tissue)

                for _, row in master_df.iterrows():
                    counts_path = _find_sample_counts(row["sample_id"], tissue, assay)
                    if counts_path or not validate_paths:
                        all_records.append(
                            {
                                "sample_id": row["sample_id"],
                                "tissue": tissue,
                                "assay": assay,
                                "counts_path": str(counts_path) if counts_path else "",
                                "genotype_id": row.get("genotype_id", row["sample_id"]),
                            }
                        )

    # Create DataFrame
    if all_records:
        df = pd.DataFrame(all_records)
        # Sort for reproducibility
        df = df.sort_values(["tissue", "assay", "sample_id"]).reset_index(drop=True)
    else:
        df = pd.DataFrame(columns=["sample_id", "tissue", "assay", "counts_path", "genotype_id"])

    manifest = SampleManifest(df=df)

    if verbose:
        manifest.print_summary()

    return manifest
