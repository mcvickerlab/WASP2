"""
WASP2 Counting CLI Module.

Provides command-line interface for counting alleles in bulk and single-cell
sequencing data across heterozygous variants.
"""

from pathlib import Path
from typing import List, Optional
from typing_extensions import Annotated

import typer
import sys
from rich.console import Console
from rich.table import Table

# Local Imports
from .run_counting import run_count_variants
from .run_counting_sc import run_count_variants_sc

# Shared CLI utilities
from wasp2.cli_utils import (
    validate_file_exists,
    handle_error,
    check_file_status,
    check_output_writable,
    create_dry_run_panel,
    create_input_validation_table,
    add_bam_to_validation_table,
)

# Version import for --version flag
try:
    from wasp2 import __version__
except ImportError:
    __version__ = "unknown"

console = Console()


def version_callback(value: bool) -> None:
    """Display version information and exit."""
    if value:
        console.print(f"wasp2-count version {__version__}")
        raise typer.Exit()


def _dry_run_report(
    bam: str,
    variants: str,
    region_file: Optional[str],
    samples: Optional[str],
    out_file: Optional[str],
    use_rust: bool,
    include_indels: bool,
) -> None:
    """Display dry-run validation report."""
    console.print(create_dry_run_panel())

    # Input validation table
    table = create_input_validation_table()
    add_bam_to_validation_table(table, bam)

    var_status, var_size = check_file_status(variants)
    table.add_row("Variants", variants, var_status, var_size)

    if region_file:
        reg_status, reg_size = check_file_status(region_file)
        table.add_row("Regions", region_file, reg_status, reg_size)

    console.print(table)
    console.print()

    # Configuration table
    config_table = Table(title="Configuration", show_header=False, box=None)
    config_table.add_column("Setting", style="cyan")
    config_table.add_column("Value", style="green")

    config_table.add_row("Sample(s)", samples or "[dim]auto-detect[/dim]")
    config_table.add_row("Output file", out_file or "counts.tsv")
    config_table.add_row("Rust acceleration", "[green]Yes[/green]" if use_rust else "[yellow]No[/yellow]")
    config_table.add_row("Include indels", "Yes" if include_indels else "No")

    out_path = Path(out_file) if out_file else Path("counts.tsv")
    out_dir = out_path.parent or Path(".")
    config_table.add_row("Output writable", check_output_writable(out_dir))

    console.print(config_table)
    console.print()

    # Rust extension check
    try:
        import wasp2_rust  # noqa: F401
        console.print("[green]✓[/green] Rust extension available")
    except ImportError:
        if use_rust:
            console.print("[yellow]![/yellow] Rust extension not found - will use Python fallback")
        else:
            console.print("[dim]Rust extension not used (--no-rust specified)[/dim]")

    console.print()
    console.print("[yellow]Dry run complete.[/yellow] Remove --dry-run to execute.")


app = typer.Typer(
    name="wasp2-count",
    help="Count alleles in sequencing data (bulk and single-cell).",
    rich_markup_mode="rich",
    no_args_is_help=True,
    pretty_exceptions_short=False,
    epilog=(
        "[dim]Examples:[/dim]\n\n"
        "  [green]Bulk counting with VCF:[/green]\n"
        "    wasp2-count count-variants aligned.bam variants.vcf.gz -s sample1 -r peaks.bed\n\n"
        "  [green]Single-cell counting:[/green]\n"
        "    wasp2-count count-variants-sc aligned.bam variants.vcf.gz barcodes.txt -f features.bed\n\n"
        "[dim]Documentation: https://github.com/mcvickerlab/WASP2[/dim]"
    ),
)


@app.callback()
def main(
    version: Optional[bool] = typer.Option(
        None,
        "--version",
        "-V",
        callback=version_callback,
        is_eager=True,
        help="Show version and exit.",
    ),
) -> None:
    """
    [bold green]WASP2 Count[/bold green]: Process allele-specific read counts.

    Count reference and alternate alleles at heterozygous SNP positions in
    sequencing data. Supports both bulk and single-cell workflows with
    optional region filtering.

    [dim]Use 'wasp2-count COMMAND --help' for command-specific options.[/dim]
    """
    pass


@app.command("count-variants")
def count_variants(
    bam: Annotated[
        str,
        typer.Argument(
            help="Input BAM file (coordinate-sorted, indexed recommended)."
        )
    ],
    variants: Annotated[
        str,
        typer.Argument(
            help="Variant file: VCF (.vcf, .vcf.gz), BCF (.bcf), or PGEN (.pgen)."
        )
    ],
    samples: Annotated[
        Optional[List[str]],
        typer.Option(
            "--samples",
            "--sample",
            "--samps",
            "-s",
            help=(
                "Sample(s) from variant file. "
                "Accepts comma-separated string or file with one sample per line."
            )
        )
    ] = None,
    region_file: Annotated[
        Optional[str],
        typer.Option(
            "--region",
            "--regions",
            "--region_file",
            "--regions_file",
            "-r",
            help=(
                "Only count variants overlapping these regions. "
                "Accepts BED, narrowPeak, broadPeak, GTF, or GFF3 files."
            )
        )
    ] = None,
    out_file: Annotated[
        Optional[str],
        typer.Option(
            "--out_file",
            "--outfile",
            "--out",
            "-o",
            help="Output file for counts. [default: counts.tsv]"
        )
    ] = None,
    temp_loc: Annotated[
        Optional[str],
        typer.Option(
            "--temp_loc",
            "--temp",
            "-t",
            help="Directory for intermediate files. [default: system temp]"
        )
    ] = None,
    use_region_names: Annotated[
        bool,
        typer.Option(
            "--use_region_names",
            help="Use region names (column 4) instead of coordinates in output."
        )
    ] = False,
    gene_feature: Annotated[
        Optional[str],
        typer.Option(
            "--gene_feature",
            "--feature",
            "--feat",
            help="Feature type in GTF/GFF3 for counting. [default: exon]"
        )
    ] = None,
    gene_attribute: Annotated[
        Optional[str],
        typer.Option(
            "--gene_attribute",
            "--attribute",
            "--attributes",
            "--attrs",
            "--attr",
            help="Attribute name from GTF/GFF3 to use as ID. [default: auto]"
        )
    ] = None,
    gene_parent: Annotated[
        Optional[str],
        typer.Option(
            "--gene_parent",
            "--parent",
            "--parent_feature",
            "--parent_attribute",
            help="Parent attribute in GTF/GFF3. [default: transcript_id or Parent]"
        )
    ] = None,
    use_rust: Annotated[
        bool,
        typer.Option(
            "--use-rust/--no-rust",
            help="Use Rust acceleration for BAM counting (faster, recommended)."
        )
    ] = True,
    vcf_bed: Annotated[
        Optional[str],
        typer.Option(
            "--vcf-bed",
            help="Precomputed VCF BED file to skip vcf_to_bed step."
        )
    ] = None,
    intersect_bed: Annotated[
        Optional[str],
        typer.Option(
            "--intersect-bed",
            help="Precomputed intersect BED file to skip bedtools intersect."
        )
    ] = None,
    include_indels: Annotated[
        bool,
        typer.Option(
            "--include-indels/--no-indels",
            help="Include indels in addition to SNPs. [default: SNPs only]"
        )
    ] = False,
    dry_run: Annotated[
        bool,
        typer.Option(
            "--dry-run",
            help="Validate inputs and show what would be done without running."
        )
    ] = False,
) -> None:
    """
    Count alleles at heterozygous variants in bulk sequencing data.

    Counts reference and alternate allele reads at each heterozygous position,
    optionally filtering to specific genomic regions (peaks, genes, etc.).
    """
    # Validate input files
    validate_file_exists(bam, "BAM file")
    validate_file_exists(variants, "Variant file")
    if region_file:
        validate_file_exists(region_file, "Region file")

    # Parse sample string
    sample_str = samples[0] if samples else None

    # Dry run mode - validate and report
    if dry_run:
        _dry_run_report(
            bam=bam,
            variants=variants,
            region_file=region_file,
            samples=sample_str,
            out_file=out_file,
            use_rust=use_rust,
            include_indels=include_indels,
        )
        raise typer.Exit(0)

    # Run counting with error handling
    try:
        console.print("[green]Starting allele counting...[/green]")
        run_count_variants(
            bam_file=bam,
            variant_file=variants,
            region_file=region_file,
            samples=sample_str,
            use_region_names=use_region_names,
            out_file=out_file,
            temp_loc=temp_loc,
            gene_feature=gene_feature,
            gene_attribute=gene_attribute,
            gene_parent=gene_parent,
            use_rust=use_rust,
            precomputed_vcf_bed=vcf_bed,
            precomputed_intersect=intersect_bed,
            include_indels=include_indels
        )
        output_path = out_file or "counts.tsv"
        console.print(f"[green]Done![/green] Output written to: [cyan]{output_path}[/cyan]")
    except Exception as e:
        # Extra hints specific to counting
        extra_hints = {
            "does not match": "Ensure the BAM file is properly indexed and variant file is valid.",
            "header": "Ensure the BAM file is properly indexed and variant file is valid.",
        }
        handle_error(e, "allele counting", extra_hints)


@app.command("count-variants-sc")
def count_variants_sc(
    bam: Annotated[
        str,
        typer.Argument(
            help="Input BAM file (coordinate-sorted, with CB tag for barcodes)."
        )
    ],
    variants: Annotated[
        str,
        typer.Argument(
            help="Variant file: VCF (.vcf, .vcf.gz), BCF (.bcf), or PGEN (.pgen)."
        )
    ],
    barcodes: Annotated[
        str,
        typer.Argument(
            help="File with one cell barcode per line (used as cell index)."
        )
    ],
    samples: Annotated[
        Optional[List[str]],
        typer.Option(
            "--samples",
            "--sample",
            "--samps",
            "-s",
            help=(
                "Sample(s) from variant file. "
                "Recommended: use one sample at a time for single-cell data."
            )
        )
    ] = None,
    feature_file: Annotated[
        Optional[str],
        typer.Option(
            "--feature",
            "--features",
            "--feat",
            "-f",
            "--region",
            "--regions",
            "-r",
            help=(
                "Feature file defining regions to count. "
                "Accepts BED or MACS2 peak files."
            )
        )
    ] = None,
    out_file: Annotated[
        Optional[str],
        typer.Option(
            "--out_file",
            "--outfile",
            "--out",
            "-o",
            help="Output AnnData file (.h5ad). [default: allele_counts.h5ad]"
        )
    ] = None,
    temp_loc: Annotated[
        Optional[str],
        typer.Option(
            "--temp_loc",
            "--temp",
            "-t",
            help="Directory for intermediate files. [default: system temp]"
        )
    ] = None,
) -> None:
    """
    Count alleles at heterozygous variants in single-cell data.

    Produces an AnnData (.h5ad) file with allele counts per cell and variant,
    suitable for single-cell allelic imbalance analysis.
    """
    # Validate input files
    validate_file_exists(bam, "BAM file")
    validate_file_exists(variants, "Variant file")
    validate_file_exists(barcodes, "Barcode file")
    if feature_file:
        validate_file_exists(feature_file, "Feature file")

    # Parse sample string
    sample_str = samples[0] if samples else None

    # Run single-cell counting with error handling
    try:
        console.print("[green]Starting single-cell allele counting...[/green]")
        run_count_variants_sc(
            bam_file=bam,
            variant_file=variants,
            barcode_file=barcodes,
            feature_file=feature_file,
            samples=sample_str,
            out_file=out_file,
            temp_loc=temp_loc
        )
        output_path = out_file or "allele_counts.h5ad"
        console.print(f"[green]Done![/green] Output written to: [cyan]{output_path}[/cyan]")
    except Exception as e:
        handle_error(e, "single-cell allele counting")


if __name__ == "__main__":
    root_dir = Path(__file__).parent
    sys.path.append(str(root_dir))
    app()
