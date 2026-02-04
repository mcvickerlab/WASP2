"""
CLI entry point for iPSCORE multi-tissue analysis module.

Usage:
    wasp2-ipscore inventory   # Validate data inventory
    wasp2-ipscore manifest    # Generate unified sample manifest
    wasp2-ipscore qtls        # Load and summarize QTL data
"""

from pathlib import Path
from typing import Literal

import typer

from wasp2.cli import error, info, success, warning

from .constants import TISSUE_LABELS

app = typer.Typer(
    name="wasp2-ipscore",
    help="iPSCORE Multi-Tissue Allelic Imbalance Resource utilities",
    add_completion=False,
)


@app.command()
def inventory(
    output: Path | None = typer.Option(
        None,
        "--output",
        "-o",
        help="Output path for inventory report (TSV)",
    ),
    quiet: bool = typer.Option(
        False,
        "--quiet",
        "-q",
        help="Suppress progress messages",
    ),
) -> None:
    """Validate iPSCORE data inventory.

    Checks existence and completeness of:
    - WASP allelic count files for all 5 datasets
    - QTL summary statistics files
    - Sample manifest files
    """
    from .data_inventory import validate_inventory

    inv = validate_inventory(verbose=not quiet)

    if output:
        df = inv.to_dataframe()
        df.to_csv(output, sep="\t", index=False)
        success(f"Inventory report saved to {output}")


@app.command()
def manifest(
    output: Path = typer.Option(
        ...,
        "--output",
        "-o",
        help="Output path for unified manifest (CSV or TSV)",
    ),
    output_format: Literal["csv", "tsv"] = typer.Option(
        "csv",
        "--format",
        "-f",
        help="Output format",
    ),
    validate_paths: bool = typer.Option(
        True,
        "--validate/--no-validate",
        help="Validate that counts paths exist",
    ),
    quiet: bool = typer.Option(
        False,
        "--quiet",
        "-q",
        help="Suppress progress messages",
    ),
) -> None:
    """Generate unified sample manifest.

    Creates a manifest linking:
    - Sample UUIDs across all tissues
    - Tissue and assay type
    - WASP counts file paths
    - Genotype IDs
    """
    from .sample_manifest import create_unified_manifest

    mani = create_unified_manifest(
        validate_paths=validate_paths,
        verbose=not quiet,
    )

    if output_format == "csv":
        mani.to_csv(output)
    else:
        mani.to_tsv(output)

    success(f"Manifest saved to {output} ({len(mani)} samples)")


@app.command()
def qtls(
    output: Path | None = typer.Option(
        None,
        "--output",
        "-o",
        help="Output path for QTL summary (TSV)",
    ),
    tissue: str | None = typer.Option(
        None,
        "--tissue",
        "-t",
        help="Filter to specific tissue (CVPC, PPC, iPSC)",
    ),
    include_finemapped: bool = typer.Option(
        False,
        "--finemapped",
        help="Load fine-mapped data (slow, 164MB file)",
    ),
    quiet: bool = typer.Option(
        False,
        "--quiet",
        "-q",
        help="Suppress progress messages",
    ),
) -> None:
    """Load and summarize QTL data.

    Loads:
    - all_caqtls.txt (36,559 caQTLs across 3 tissues)
    - Optionally: fine-mapped iPSC data with PIPs
    """
    from .qtl_loader import create_qtl_loader

    # Validate tissue parameter if provided
    if tissue and tissue not in TISSUE_LABELS:
        error(f"Invalid tissue: '{tissue}'. Valid options: {', '.join(TISSUE_LABELS)}")
        raise typer.Exit(code=1)

    loader = create_qtl_loader(
        load_finemapped=include_finemapped,
        verbose=not quiet,
    )

    if tissue and loader.caqtls is not None:
        filtered = loader.filter_by_tissue(tissue)
        info(f"\nFiltered to {tissue}: {len(filtered):,} QTLs")

        if output:
            filtered.to_csv(output, sep="\t", index=False)
            success(f"Filtered QTLs saved to {output}")
    elif output and loader.caqtls is not None:
        loader.caqtls.to_csv(output, sep="\t", index=False)
        success(f"All QTLs saved to {output}")


@app.command()
def validate(
    quiet: bool = typer.Option(
        False,
        "--quiet",
        "-q",
        help="Suppress progress messages",
    ),
) -> None:
    """Run full validation of all iPSCORE data resources.

    Combines inventory validation with QTL data checks.
    Reports any missing or incomplete data.
    """
    from .data_inventory import validate_inventory
    from .qtl_loader import create_qtl_loader
    from .sample_manifest import create_unified_manifest

    info("=" * 70)
    info("iPSCORE Multi-Tissue Resource Validation")
    info("=" * 70)

    # Step 1: Inventory
    info("\n[Step 1/3] Validating data inventory...")
    inv = validate_inventory(verbose=not quiet)

    # Step 2: Sample manifest
    info("\n[Step 2/3] Creating sample manifest...")
    mani = create_unified_manifest(verbose=not quiet)

    # Step 3: QTL data
    info("\n[Step 3/3] Loading QTL data...")
    qtl_loader = create_qtl_loader(
        load_finemapped=False,  # Skip large file in validation
        verbose=not quiet,
    )

    # Summary
    info("\n" + "=" * 70)
    info("VALIDATION SUMMARY")
    info("=" * 70)

    total_ok = sum(1 for ds in inv.datasets.values() if ds.complete)
    total_ds = len(inv.datasets)

    if total_ok == total_ds:
        success(f"Datasets: {total_ok}/{total_ds} complete")
    else:
        warning(f"Datasets: {total_ok}/{total_ds} complete")

    info(f"Samples in manifest: {len(mani)}")
    info(f"QTLs loaded: {qtl_loader.total_qtls:,}")

    # Check against expected totals with explicit reporting
    missing = inv.total_expected - inv.total_samples
    if missing > 0:
        warning(f"Missing samples: {missing}")

    if inv.total_samples >= inv.total_expected * 0.95:
        success("Data inventory validation PASSED")
    else:
        error(
            f"Data inventory validation FAILED: "
            f"{inv.total_samples}/{inv.total_expected} samples found"
        )
        raise typer.Exit(code=1)


if __name__ == "__main__":
    app()
