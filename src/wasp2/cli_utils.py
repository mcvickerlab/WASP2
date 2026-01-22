"""
Shared CLI utilities for WASP2 modules.

Provides common functions for file validation, error handling, and dry-run
reporting used across counting, mapping, and analysis CLI modules.
"""

from pathlib import Path
from typing import Optional
import os

from rich.console import Console
from rich.panel import Panel
from rich.table import Table

console = Console()


def validate_file_exists(path: str, description: str) -> None:
    """Validate that a file exists and exit with error if not found.

    Args:
        path: Path to the file to check
        description: Human-readable description of the file (e.g., "BAM file")
    """
    import typer
    if not Path(path).is_file():
        console.print(
            Panel(
                f"[red]Error:[/red] {description} not found:\n"
                f"  [yellow]{path}[/yellow]\n\n"
                f"Please check that the file path is correct and the file exists.",
                title="File Not Found",
                border_style="red",
            )
        )
        raise typer.Exit(code=1)


def format_file_size(size_bytes: int) -> str:
    """Format file size in human-readable format."""
    if size_bytes > 1_000_000_000:
        return f"{size_bytes / 1_000_000_000:.1f} GB"
    elif size_bytes > 1_000_000:
        return f"{size_bytes / 1_000_000:.1f} MB"
    return f"{size_bytes / 1_000:.1f} KB"


def check_file_status(path: str) -> tuple[str, str]:
    """Check if file exists and return status and formatted size.

    Returns:
        Tuple of (status_string, size_string) where status is Rich-formatted
    """
    p = Path(path)
    if p.is_file():
        size_str = format_file_size(p.stat().st_size)
        return "[green]✓ Found[/green]", size_str
    return "[red]✗ Not found[/red]", "-"


def check_output_writable(out_path: Path) -> str:
    """Check if output path is writable and return Rich-formatted status."""
    if out_path.exists() and os.access(out_path, os.W_OK):
        return "[green]Yes[/green]"
    if not out_path.exists():
        parent = out_path.parent or Path(".")
        if parent.exists() and os.access(parent, os.W_OK):
            return "[green]Yes (will create)[/green]"
    return "[red]No - check permissions[/red]"


def create_dry_run_panel(subtitle: str = "Validating inputs without running analysis") -> Panel:
    """Create the dry-run mode header panel."""
    return Panel(
        f"[bold]Dry Run Mode[/bold]\n[dim]{subtitle}[/dim]",
        border_style="yellow",
    )


def create_input_validation_table() -> Table:
    """Create a standard table for input file validation."""
    table = Table(title="Input Validation", show_header=True)
    table.add_column("File", style="cyan")
    table.add_column("Path", style="dim")
    table.add_column("Status", style="green")
    table.add_column("Size", style="dim")
    return table


def add_bam_to_validation_table(table: Table, bam_path: str) -> None:
    """Add BAM and BAM index rows to validation table."""
    bam_status, bam_size = check_file_status(bam_path)
    table.add_row("BAM", bam_path, bam_status, bam_size)

    # Check BAM index (.bai or with .bam replaced)
    bam_idx = bam_path + ".bai"
    if not Path(bam_idx).exists():
        bam_idx = bam_path.replace(".bam", ".bai")
    idx_status, idx_size = check_file_status(bam_idx)
    table.add_row("BAM Index", bam_idx, idx_status, idx_size)


def handle_error(
    e: Exception,
    context: str,
    extra_hints: Optional[dict[str, str]] = None
) -> None:
    """Handle exceptions with user-friendly error messages and hints.

    Args:
        e: The exception that occurred
        context: Description of what was being done (e.g., "allele counting")
        extra_hints: Additional error patterns and hints to check
    """
    import typer
    error_msg = str(e)

    # Base hints for common errors
    hints = {
        "No module named": "Try installing missing dependencies: pip install wasp2[plink,cyvcf2]",
        "Permission denied": "Check file permissions or try a different output location.",
        "No space left": "Free up disk space or use --temp to specify a different temp directory.",
    }

    # Add extra hints if provided
    if extra_hints:
        hints.update(extra_hints)

    # Find matching hint
    hint = "Check input files and parameters. Use --help for usage information."
    for pattern, suggestion in hints.items():
        if pattern.lower() in error_msg.lower():
            hint = suggestion
            break

    console.print(
        Panel(
            f"[red]Error during {context}:[/red]\n"
            f"  {error_msg}\n\n"
            f"[dim]Hint: {hint}[/dim]",
            title="Error",
            border_style="red",
        )
    )
    raise typer.Exit(code=1)
