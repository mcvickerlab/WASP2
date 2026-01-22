"""
WASP2 Unified Command Line Interface.

Provides a polished CLI experience with subcommands for all WASP2 operations:
- count: Count alleles in sequencing data
- analyze: Analyze allelic imbalance
- map: Mapping pipeline for bias correction
"""

from typing import Optional

import typer
from rich.console import Console
from rich.panel import Panel
from rich.text import Text

from wasp2 import __version__

# Create Rich console for styled output
console = Console()

# Version callback for --version flag
def version_callback(value: bool) -> None:
    """Display version information and exit."""
    if value:
        console.print(
            Panel(
                Text.from_markup(
                    f"[bold blue]WASP2[/bold blue] version [green]{__version__}[/green]\n"
                    f"[dim]Allele-Specific Pipeline for Unbiased Read Mapping[/dim]"
                ),
                title="Version Info",
                border_style="blue",
            )
        )
        raise typer.Exit()


# Main application with Rich formatting
app = typer.Typer(
    name="wasp2",
    help="WASP2: Allele-Specific Analysis Pipeline",
    rich_markup_mode="rich",
    no_args_is_help=True,
    add_completion=True,
    pretty_exceptions_enable=True,
    pretty_exceptions_show_locals=False,
)


@app.callback()
def main(
    version: Optional[bool] = typer.Option(
        None,
        "--version",
        "-V",
        callback=version_callback,
        is_eager=True,
        help="Show version information and exit.",
    ),
) -> None:
    """
    [bold blue]WASP2[/bold blue]: Toolkit for allele-specific analysis and unbiased read mapping.

    [dim]A comprehensive pipeline for:[/dim]
    - Counting alleles across heterozygous SNPs
    - Analyzing allelic imbalance in genomic regions
    - Correcting mapping biases in NGS data

    [yellow]Documentation:[/yellow] https://github.com/mcvickerlab/WASP2
    """
    pass


# Import and register subcommands
def register_subcommands() -> None:
    """Register all subcommand modules."""
    # Import subcommand apps
    from counting.__main__ import app as counting_app
    from mapping.__main__ import app as mapping_app
    from analysis.__main__ import app as analysis_app

    # Add subcommands with descriptive help
    app.add_typer(
        counting_app,
        name="count",
        help="[green]Count[/green] alleles in sequencing data (bulk and single-cell)",
    )
    app.add_typer(
        analysis_app,
        name="analyze",
        help="[cyan]Analyze[/cyan] allelic imbalance across genomic regions",
    )
    app.add_typer(
        mapping_app,
        name="map",
        help="[magenta]Map[/magenta] reads with WASP bias correction pipeline",
    )


# Register on import
register_subcommands()


def cli() -> None:
    """Entry point for the CLI."""
    app()


if __name__ == "__main__":
    cli()
