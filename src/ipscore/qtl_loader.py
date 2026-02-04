"""
QTL data loading and harmonization for iPSCORE resource.

Handles:
- all_caqtls.txt (36,559 caQTLs across 3 tissues)
- ipsc_caqtl_finemapped.txt.gz (fine-mapped with PIPs and credible sets)
- Tissue-specific summary statistics
"""

from dataclasses import dataclass

import pandas as pd

from wasp2.cli import error, info, success

from .constants import QTL_COUNTS, QTL_FILES


@dataclass
class QTLLoader:
    """Container for loaded QTL data."""

    caqtls: pd.DataFrame | None = None
    finemapped: pd.DataFrame | None = None

    @property
    def total_qtls(self) -> int:
        """Total caQTLs loaded."""
        if self.caqtls is None:
            return 0
        return len(self.caqtls)

    @property
    def tissues(self) -> list[str]:
        """Unique tissues in QTL data."""
        if self.caqtls is None:
            return []
        return self.caqtls["tissue"].unique().tolist()

    def filter_by_tissue(self, tissue: str) -> pd.DataFrame:
        """Filter caQTLs by tissue.

        Args:
            tissue: Tissue label (CVPC, PPC, iPSC)

        Returns:
            Filtered DataFrame
        """
        if self.caqtls is None:
            return pd.DataFrame()
        return self.caqtls[self.caqtls["tissue"] == tissue].copy()

    def get_qtl_counts(self) -> dict[str, int]:
        """Get QTL counts per tissue."""
        if self.caqtls is None:
            return {}
        return self.caqtls["tissue"].value_counts().to_dict()

    def print_summary(self) -> None:
        """Print QTL loading summary."""
        info("=" * 60)
        info("QTL Data Summary")
        info("=" * 60)

        if self.caqtls is not None:
            info(f"\nTotal caQTLs: {self.total_qtls:,}")
            for tissue, count in self.get_qtl_counts().items():
                expected = QTL_COUNTS.get(tissue, "?")
                info(f"  {tissue}: {count:,} (expected: {expected:,})")

        if self.finemapped is not None:
            info(f"\nFine-mapped records: {len(self.finemapped):,}")
            if "Credible Set" in self.finemapped.columns:
                credible_count = self.finemapped["Credible Set"].sum()
                info(f"  In 99% credible sets: {credible_count:,}")


def load_all_caqtls(verbose: bool = True) -> pd.DataFrame:
    """Load the combined caQTL file for all 3 tissues.

    File format (from Issue #40):
        element_id  type  element_chrom  element_start  element_end  id  tissue

    Args:
        verbose: Print progress messages

    Returns:
        DataFrame with all caQTLs and standardized columns
    """
    qtl_path = QTL_FILES["all_caqtls"]

    if not qtl_path.exists():
        error(f"QTL file not found: {qtl_path}")
        return pd.DataFrame()

    if verbose:
        info(f"Loading caQTLs from {qtl_path.name}...")

    df = pd.read_csv(
        qtl_path,
        sep="\t",
        dtype={
            "element_id": str,
            "type": "int8",
            "element_chrom": "category",
            "element_start": "int32",
            "element_end": "int32",
            "id": str,
            "tissue": "category",
        },
    )

    # Standardize column names
    df = df.rename(
        columns={
            "element_id": "peak_id",
            "element_chrom": "chrom",
            "element_start": "start",
            "element_end": "end",
            "id": "variant_id",
        }
    )

    # Parse variant_id to extract position info
    # Format: VAR_chr_pos_ref_alt
    if "variant_id" in df.columns:
        variant_parts = df["variant_id"].str.split("_", expand=True)
        if variant_parts.shape[1] >= 3:
            df["var_chrom"] = "chr" + variant_parts[1].astype(str)
            df["var_pos"] = pd.to_numeric(variant_parts[2], errors="coerce")

    if verbose:
        success(f"Loaded {len(df):,} caQTLs across {df['tissue'].nunique()} tissues")

    return df


def load_finemapped_qtls(
    verbose: bool = True,
    nrows: int | None = None,
) -> pd.DataFrame:
    """Load fine-mapped iPSC caQTL data with PIPs and credible sets.

    File: ipsc_caqtl_finemapped.txt.gz (164 MB)

    Expected columns:
        SNP ID, Position, Element ID, SNP.PP (PIP), Credible Set

    Args:
        verbose: Print progress messages
        nrows: Limit number of rows (for testing)

    Returns:
        DataFrame with fine-mapping results
    """
    finemapped_path = QTL_FILES["ipsc_finemapped"]

    if not finemapped_path.exists():
        error(f"Fine-mapped file not found: {finemapped_path}")
        return pd.DataFrame()

    if verbose:
        info(f"Loading fine-mapped data from {finemapped_path.name}...")

    df = pd.read_csv(
        finemapped_path,
        sep="\t",
        nrows=nrows,
        dtype={
            "SNP ID": str,
            "Position": "int32",
            "Element ID": str,
        },
    )

    # Standardize column names
    col_mapping = {
        "SNP ID": "variant_id",
        "Position": "position",
        "Element ID": "peak_id",
        "SNP.PP": "pip",
        "Credible Set": "in_credible_set",
    }
    df = df.rename(columns={k: v for k, v in col_mapping.items() if k in df.columns})

    # Convert credible set to boolean if it's a string
    if "in_credible_set" in df.columns:
        if df["in_credible_set"].dtype == object:
            df["in_credible_set"] = df["in_credible_set"].str.upper() == "TRUE"

    if verbose:
        success(f"Loaded {len(df):,} fine-mapped records")
        if "in_credible_set" in df.columns:
            credible_count = df["in_credible_set"].sum()
            info(f"  Variants in 99% credible sets: {credible_count:,}")

    return df


def create_qtl_loader(
    load_finemapped: bool = True,
    finemapped_nrows: int | None = None,
    verbose: bool = True,
) -> QTLLoader:
    """Create a QTLLoader with all iPSCORE QTL data.

    Args:
        load_finemapped: Whether to load fine-mapped data (large file)
        finemapped_nrows: Limit fine-mapped rows (for testing)
        verbose: Print progress messages

    Returns:
        QTLLoader with loaded data
    """
    if verbose:
        info("=" * 60)
        info("Loading iPSCORE QTL Data")
        info("=" * 60)

    loader = QTLLoader()

    # Load combined caQTLs
    loader.caqtls = load_all_caqtls(verbose=verbose)

    # Optionally load fine-mapped data
    if load_finemapped:
        loader.finemapped = load_finemapped_qtls(verbose=verbose, nrows=finemapped_nrows)

    if verbose:
        loader.print_summary()

    return loader
