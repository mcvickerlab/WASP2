"""
WASP2 Unified Command Line Interface.

Provides a polished CLI experience with subcommands for all WASP2 operations:
- count: Count alleles in sequencing data
- analyze: Analyze allelic imbalance
- map: Mapping pipeline for bias correction
- info: Display system and environment information
- config: Manage configuration settings
"""

import platform
import sys
from pathlib import Path
from typing import Optional

import typer
from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.text import Text

from wasp2 import __version__
from wasp2.config import (
    WASP2Config,
    get_config,
    get_config_path,
    load_config,
    save_config,
    setup_logging,
)

console = Console()
_verbosity: int = 0


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
    verbose: int = typer.Option(
        0,
        "--verbose",
        "-v",
        count=True,
        help="Increase verbosity. Use -v for info, -vv for debug, -vvv for trace.",
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
    global _verbosity
    _verbosity = verbose

    # Set up logging based on verbosity
    config = get_config()
    log_file = config.log_file if hasattr(config, "log_file") else None
    setup_logging(verbose, log_file)


@app.command("info")
def info_cmd() -> None:
    """
    Display system information and WASP2 environment details.

    Shows Python version, installed dependencies, Rust extension status,
    and current configuration settings.
    """
    # System info table
    sys_table = Table(title="System Information", show_header=False, box=None)
    sys_table.add_column("Key", style="cyan")
    sys_table.add_column("Value", style="green")

    sys_table.add_row("WASP2 Version", __version__)
    sys_table.add_row("Python Version", platform.python_version())
    sys_table.add_row("Platform", platform.platform())
    sys_table.add_row("Python Executable", sys.executable)

    # Check Rust extension
    try:
        import wasp2_rust
        rust_status = f"[green]✓ Installed[/green] (v{getattr(wasp2_rust, '__version__', 'unknown')})"
    except ImportError:
        rust_status = "[yellow]✗ Not installed[/yellow] (using Python fallback)"

    sys_table.add_row("Rust Extension", rust_status)

    console.print(sys_table)
    console.print()

    # Dependencies table
    dep_table = Table(title="Key Dependencies", show_header=True)
    dep_table.add_column("Package", style="cyan")
    dep_table.add_column("Version", style="green")
    dep_table.add_column("Status", style="dim")

    deps = [
        ("numpy", "numpy"),
        ("pandas", "pandas"),
        ("pysam", "pysam"),
        ("pybedtools", "pybedtools"),
        ("typer", "typer"),
        ("rich", "rich"),
        ("anndata", "anndata"),
        ("cyvcf2", "cyvcf2"),
        ("pgenlib", "pgenlib"),
    ]

    for name, module in deps:
        try:
            mod = __import__(module)
            version = getattr(mod, "__version__", "installed")
            dep_table.add_row(name, version, "[green]✓[/green]")
        except ImportError:
            dep_table.add_row(name, "-", "[dim]not installed[/dim]")

    console.print(dep_table)
    console.print()

    # Config info
    config_path = get_config_path()
    config = load_config()

    config_table = Table(title="Configuration", show_header=False, box=None)
    config_table.add_column("Key", style="cyan")
    config_table.add_column("Value", style="green")

    config_table.add_row("Config File", str(config_path))
    config_table.add_row("  Exists", "[green]Yes[/green]" if config_path.exists() else "[dim]No (using defaults)[/dim]")
    config_table.add_row("Default Threads", str(config.threads))
    config_table.add_row("Use Rust", "[green]Yes[/green]" if config.use_rust else "[yellow]No[/yellow]")
    config_table.add_row("Log Level", config.log_level)

    console.print(config_table)


config_app = typer.Typer(
    name="config",
    help="Manage WASP2 configuration settings.",
    rich_markup_mode="rich",
    no_args_is_help=True,
)


@config_app.command("show")
def config_show() -> None:
    """Display current configuration settings."""
    config = load_config()
    config_path = get_config_path()

    console.print(f"[dim]Config file: {config_path}[/dim]\n")

    table = Table(title="Current Settings", show_header=True)
    table.add_column("Setting", style="cyan")
    table.add_column("Value", style="green")
    table.add_column("Description", style="dim")

    descriptions = {
        "threads": "Default thread count for parallel operations",
        "use_rust": "Use Rust acceleration when available",
        "temp_dir": "Default temporary directory",
        "output_format": "Default output format (tsv, json, yaml)",
        "log_level": "Logging verbosity (DEBUG, INFO, WARNING, ERROR)",
        "log_file": "Optional log file path",
        "include_indels": "Include indels by default",
        "pseudocount": "Default pseudocount for imbalance analysis",
        "min_count": "Default minimum count threshold",
    }

    for key, value in config.to_dict().items():
        desc = descriptions.get(key, "")
        display_value = str(value) if value is not None else "[dim]not set[/dim]"
        table.add_row(key, display_value, desc)

    console.print(table)


@config_app.command("set")
def config_set(
    key: str = typer.Argument(help="Configuration key to set"),
    value: str = typer.Argument(help="Value to set"),
) -> None:
    """Set a configuration value.

    Example: wasp2 config set threads 4
    """
    config = load_config()

    # Validate key exists
    if not hasattr(config, key):
        valid_keys = list(config.to_dict().keys())
        console.print(f"[red]Error:[/red] Unknown config key: {key}")
        console.print(f"[dim]Valid keys: {', '.join(valid_keys)}[/dim]")
        raise typer.Exit(code=1)

    # Type conversion based on field type
    current_value = getattr(config, key)
    try:
        if isinstance(current_value, bool):
            # Handle boolean strings
            parsed = value.lower() in ("true", "yes", "1", "on")
        elif isinstance(current_value, int):
            parsed = int(value)
        else:
            parsed = value if value.lower() != "none" else None
    except ValueError:
        console.print(f"[red]Error:[/red] Invalid value '{value}' for {key}")
        raise typer.Exit(code=1)

    setattr(config, key, parsed)
    save_config(config)

    console.print(f"[green]✓[/green] Set {key} = {parsed}")
    console.print(f"[dim]Saved to {get_config_path()}[/dim]")


@config_app.command("reset")
def config_reset(
    force: bool = typer.Option(False, "--force", "-f", help="Skip confirmation"),
) -> None:
    """Reset configuration to defaults."""
    if not force:
        confirm = typer.confirm("Reset all settings to defaults?")
        if not confirm:
            raise typer.Abort()

    config = WASP2Config()
    save_config(config)
    console.print("[green]✓[/green] Configuration reset to defaults")


@config_app.command("path")
def config_path() -> None:
    """Show the configuration file path."""
    path = get_config_path()
    console.print(f"Config file: [cyan]{path}[/cyan]")
    if path.exists():
        console.print("[green]✓[/green] File exists")
    else:
        console.print("[dim]File does not exist (using defaults)[/dim]")


app.add_typer(config_app, name="config", help="Manage configuration settings")


def register_subcommands() -> None:
    """Register all subcommand modules."""
    from counting.__main__ import app as counting_app
    from mapping.__main__ import app as mapping_app
    from analysis.__main__ import app as analysis_app

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
