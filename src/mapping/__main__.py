"""
WASP2 Mapping CLI Module.

Provides command-line interface for the WASP mapping pipeline, which corrects
mapping biases by re-mapping reads with swapped alleles.
"""

from pathlib import Path
from typing import List, Optional
from typing_extensions import Annotated

import typer
import sys
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

# Local Imports
from .run_mapping import run_make_remap_reads, run_wasp_filt

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
        console.print(f"wasp2-map version {__version__}")
        raise typer.Exit()


def _dry_run_mapping_report(
    bam: str,
    variants: str,
    samples: Optional[str],
    out_dir: Optional[str],
    is_paired: Optional[bool],
    is_phased: Optional[bool],
    include_indels: bool,
) -> None:
    """Display dry-run validation report for mapping."""
    console.print(create_dry_run_panel("Validating inputs without running mapping pipeline"))

    # Input validation table
    table = create_input_validation_table()
    add_bam_to_validation_table(table, bam)

    var_status, var_size = check_file_status(variants)
    table.add_row("Variants", variants, var_status, var_size)

    console.print(table)
    console.print()

    # Configuration table
    config_table = Table(title="Configuration", show_header=False, box=None)
    config_table.add_column("Setting", style="cyan")
    config_table.add_column("Value", style="green")

    config_table.add_row("Sample(s)", samples or "[dim]auto-detect[/dim]")
    config_table.add_row("Output directory", out_dir or "[dim]current directory[/dim]")

    # Format read type
    if is_paired is True:
        read_type = "paired-end"
    elif is_paired is False:
        read_type = "single-end"
    else:
        read_type = "[dim]auto-detect[/dim]"
    config_table.add_row("Read type", read_type)

    # Format phasing
    if is_phased is True:
        phased_str = "Yes"
    elif is_phased is False:
        phased_str = "No"
    else:
        phased_str = "[dim]auto-detect[/dim]"
    config_table.add_row("Phased", phased_str)

    config_table.add_row("Include indels", "Yes" if include_indels else "No")

    # Check output directory writability
    out_path = Path(out_dir) if out_dir else Path(".")
    config_table.add_row("Output writable", check_output_writable(out_path))

    console.print(config_table)
    console.print()
    console.print("[yellow]Dry run complete.[/yellow] Remove --dry-run to execute.")


app = typer.Typer(
    name="wasp2-map",
    help="WASP mapping pipeline for bias correction.",
    rich_markup_mode="rich",
    no_args_is_help=True,
    pretty_exceptions_short=False,
    epilog=(
        "[dim]Examples:[/dim]\n\n"
        "  [green]Generate remap reads:[/green]\n"
        "    wasp2-map make-reads aligned.bam variants.vcf.gz -s sample1 -o remap_dir/\n\n"
        "  [green]Filter remapped reads:[/green]\n"
        "    wasp2-map filter-remapped remapped.bam -j wasp_data.json -o filtered.bam\n\n"
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
    [bold magenta]WASP2 Map[/bold magenta]: Unbiased allele-specific read mapping.

    Corrects mapping biases by generating reads with swapped alleles,
    re-mapping them, and filtering reads that map to different locations.

    [dim]Workflow:[/dim]
    1. make-reads: Generate alternate-allele reads for remapping
    2. (External) Realign the remap FASTQs with your aligner
    3. filter-remapped: Filter reads that map to different positions

    [dim]Use 'wasp2-map COMMAND --help' for command-specific options.[/dim]
    """
    pass


@app.command("make-reads")
def make_reads(
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
    out_dir: Annotated[
        Optional[str],
        typer.Option(
            "--out_dir",
            "--outdir",
            "--out",
            "-o",
            help="Output directory for remap FASTQs and data files."
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
    out_json: Annotated[
        Optional[str],
        typer.Option(
            "--out_json",
            "--json",
            "--outjson",
            "-j",
            help="Output JSON file with WASP data paths. [default: <bam>_wasp_data_files.json]"
        )
    ] = None,
    is_paired: Annotated[
        Optional[bool],
        typer.Option(
            "--paired/--single",
            help="Specify paired-end or single-end reads. [default: auto-detect]"
        )
    ] = None,
    is_phased: Annotated[
        Optional[bool],
        typer.Option(
            "--phased/--unphased",
            help=(
                "Specify if variants are phased. Phased data strongly recommended. "
                "[default: auto-detect]"
            )
        )
    ] = None,
    include_indels: Annotated[
        bool,
        typer.Option(
            "--indels/--snps-only",
            help="Include indels in addition to SNPs. [default: SNPs only]"
        )
    ] = False,
    max_indel_len: Annotated[
        int,
        typer.Option(
            "--max-indel-len",
            help="Maximum indel length (bp) to process. Longer indels are skipped.",
            min=1
        )
    ] = 10,
    insert_qual: Annotated[
        int,
        typer.Option(
            "--insert-qual",
            help="Quality score (Phred) for inserted bases in alternate reads.",
            min=0,
            max=60
        )
    ] = 30,
    max_seqs: Annotated[
        int,
        typer.Option(
            "--max-seqs",
            help=(
                "Maximum alternate sequences per read. "
                "Reads with more variants are skipped to prevent combinatorial explosion."
            ),
            min=1
        )
    ] = 64,
    threads: Annotated[
        int,
        typer.Option(
            "--threads",
            help="Threads for BAM I/O operations.",
            min=1
        )
    ] = 1,
    dry_run: Annotated[
        bool,
        typer.Option(
            "--dry-run",
            help="Validate inputs and show what would be done without running."
        )
    ] = False,
) -> None:
    """
    Generate reads with swapped alleles for remapping.

    Creates FASTQ files containing reads with alternate alleles at heterozygous
    positions. These should be realigned with your aligner, then filtered with
    'filter-remapped' to remove mapping-biased reads.
    """
    # Validate input files
    validate_file_exists(bam, "BAM file")
    validate_file_exists(variants, "Variant file")

    # Parse sample string
    sample_str = samples[0] if samples else None

    # Dry run mode - validate and report
    if dry_run:
        _dry_run_mapping_report(
            bam=bam,
            variants=variants,
            samples=sample_str,
            out_dir=out_dir,
            is_paired=is_paired,
            is_phased=is_phased,
            include_indels=include_indels,
        )
        raise typer.Exit(0)

    # Run with error handling
    try:
        console.print("[magenta]Generating remap reads...[/magenta]")
        run_make_remap_reads(
            bam_file=bam,
            variant_file=variants,
            samples=sample_str,
            out_dir=out_dir,
            temp_loc=temp_loc,
            out_json=out_json,
            is_paired=is_paired,
            is_phased=is_phased,
            include_indels=include_indels,
            max_indel_len=max_indel_len,
            insert_qual=insert_qual,
            max_seqs=max_seqs,
            threads=threads
        )
        console.print("[green]Done![/green] Remap reads generated successfully.")
        console.print(
            "[dim]Next step: Realign the remap FASTQs, then run 'wasp2-map filter-remapped'[/dim]"
        )
    except Exception as e:
        # Extra hints specific to mapping
        extra_hints = {
            "index": "Ensure the BAM file has an index (.bai). Run: samtools index <bam>",
            "bai": "Ensure the BAM file has an index (.bai). Run: samtools index <bam>",
            "wasp2_rust": "The Rust extension is not installed. Try: pip install wasp2 --force-reinstall",
        }
        handle_error(e, "remap read generation", extra_hints)


@app.command("filter-remapped")
def filter_remapped(
    remapped_bam: Annotated[
        str,
        typer.Argument(
            help="BAM file with realigned remap reads."
        )
    ],
    to_remap_bam: Annotated[
        Optional[str],
        typer.Argument(
            help="Original to_remap BAM from make-reads (or use --json)."
        )
    ] = None,
    keep_bam: Annotated[
        Optional[str],
        typer.Argument(
            help="BAM with reads that did not need remapping (or use --json)."
        )
    ] = None,
    wasp_data_json: Annotated[
        Optional[str],
        typer.Option(
            "--wasp_data_json",
            "--json",
            "-j",
            help="JSON file from make-reads with paths to to_remap_bam and keep_bam."
        )
    ] = None,
    out_bam: Annotated[
        Optional[str],
        typer.Option(
            "--out_bam",
            "--outbam",
            "--out",
            "-o",
            help="Output filtered BAM file. [default: auto-generated name]"
        )
    ] = None,
    remap_keep_bam: Annotated[
        Optional[str],
        typer.Option(
            "--remap_keep_bam",
            help="Also output BAM with remapped reads that passed filter."
        )
    ] = None,
    remap_keep_file: Annotated[
        Optional[str],
        typer.Option(
            "--remap_keep_file",
            help="Also output text file with names of kept reads."
        )
    ] = None,
    threads: Annotated[
        int,
        typer.Option(
            "--threads",
            help="Threads for BAM I/O (Rust filter supports multithreading).",
            min=1
        )
    ] = 1,
    use_rust: Annotated[
        bool,
        typer.Option(
            "--use-rust/--no-rust",
            help="Use Rust acceleration if available (faster, recommended)."
        )
    ] = True,
    same_locus_slop: Annotated[
        int,
        typer.Option(
            "--same-locus-slop",
            help=(
                "Tolerance (bp) for 'same locus' test. "
                "Use 2-3 for indels (micro-homology), 0 for strict SNP matching."
            ),
            min=0
        )
    ] = 0,
) -> None:
    """
    Filter remapped reads using the WASP algorithm.

    Compares original and remapped read positions. Reads that map to different
    locations after allele-swapping are removed to eliminate mapping bias.
    """
    # Validate input files
    validate_file_exists(remapped_bam, "Remapped BAM file")

    # Check that we have the required input files
    if wasp_data_json is None and (to_remap_bam is None or keep_bam is None):
        console.print(
            Panel(
                "[red]Error:[/red] Missing required input files.\n\n"
                "You must provide either:\n"
                "  1. --json <wasp_data.json> (recommended), or\n"
                "  2. Both to_remap_bam and keep_bam positional arguments\n\n"
                "[dim]The JSON file is created by 'wasp2-map make-reads'.[/dim]",
                title="Missing Input",
                border_style="red",
            )
        )
        raise typer.Exit(code=1)

    if wasp_data_json:
        validate_file_exists(wasp_data_json, "WASP data JSON file")
    if to_remap_bam:
        validate_file_exists(to_remap_bam, "To-remap BAM file")
    if keep_bam:
        validate_file_exists(keep_bam, "Keep BAM file")

    # Run with error handling
    try:
        console.print("[magenta]Filtering remapped reads...[/magenta]")
        run_wasp_filt(
            remapped_bam,
            to_remap_bam=to_remap_bam,
            keep_bam=keep_bam,
            wasp_out_bam=out_bam,
            remap_keep_bam=remap_keep_bam,
            remap_keep_file=remap_keep_file,
            wasp_data_json=wasp_data_json,
            threads=threads,
            use_rust=use_rust,
            same_locus_slop=same_locus_slop,
        )
        output_path = out_bam or "wasp_filtered.bam"
        console.print(f"[green]Done![/green] Filtered BAM written to: [cyan]{output_path}[/cyan]")
    except Exception as e:
        extra_hints = {
            "index": "Ensure the BAM file has an index (.bai). Run: samtools index <bam>",
            "wasp2_rust": "The Rust extension is not installed. Try: pip install wasp2 --force-reinstall",
        }
        handle_error(e, "WASP filtering", extra_hints)


if __name__ == "__main__":
    root_dir = Path(__file__).parent
    sys.path.append(str(root_dir))
    app()
