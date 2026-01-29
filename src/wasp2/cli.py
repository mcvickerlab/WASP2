"""WASP2 CLI utilities with Rich output formatting.

This module provides centralized CLI output functions with:
- Consistent Rich-formatted output (colors, spinners, progress bars)
- Verbosity control (verbose/normal/quiet modes)
- Progress tracking for long-running operations
- Shared CLI callbacks for version and verbosity
"""

from __future__ import annotations

import sys
from enum import IntEnum
from typing import TYPE_CHECKING, Any

from rich.console import Console

if TYPE_CHECKING:
    from collections.abc import Callable

from rich.progress import (
    BarColumn,
    Progress,
    SpinnerColumn,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
)
from rich.table import Table

# Global console instance
console = Console(stderr=True)


# Verbosity levels
class Verbosity(IntEnum):
    """Verbosity levels for CLI output."""

    QUIET = 0  # Only errors
    NORMAL = 1  # Standard output (default)
    VERBOSE = 2  # Detailed output


# Module-level verbosity setting
_verbosity: Verbosity = Verbosity.NORMAL


def set_verbosity(level: Verbosity | int) -> None:
    """Set the global verbosity level.

    Parameters
    ----------
    level : Verbosity | int
        Verbosity level (0=quiet, 1=normal, 2=verbose).
    """
    global _verbosity
    _verbosity = Verbosity(level)


def get_verbosity() -> Verbosity:
    """Get the current verbosity level."""
    return _verbosity


def is_quiet() -> bool:
    """Check if output is in quiet mode."""
    return _verbosity == Verbosity.QUIET


def is_verbose() -> bool:
    """Check if output is in verbose mode."""
    return _verbosity == Verbosity.VERBOSE


# Output functions
def info(message: str, verbose_only: bool = False) -> None:
    """Print an info message (blue text).

    Parameters
    ----------
    message : str
        Message to print.
    verbose_only : bool
        If True, only print in verbose mode.
    """
    if is_quiet() or (verbose_only and not is_verbose()):
        return
    console.print(f"[blue]{message}[/blue]")


def success(message: str) -> None:
    """Print a success message (green text with checkmark)."""
    if is_quiet():
        return
    console.print(f"[green]✓[/green] {message}")


def warning(message: str) -> None:
    """Print a warning message (yellow text)."""
    if is_quiet():
        return
    console.print(f"[yellow]⚠[/yellow] [yellow]{message}[/yellow]")


def error(message: str) -> None:
    """Print an error message (red text). Always shown."""
    console.print(f"[red]✗[/red] [red]{message}[/red]")


def status(message: str) -> None:
    """Print a status message (cyan text). Respects verbosity."""
    if is_quiet():
        return
    console.print(f"[cyan]{message}[/cyan]")


def detail(message: str) -> None:
    """Print detailed output (dim text). Only in verbose mode."""
    if not is_verbose():
        return
    console.print(f"[dim]{message}[/dim]")


def rust_status(message: str) -> None:
    """Print a Rust-related status message with crab emoji."""
    if is_quiet():
        return
    console.print(f"[orange1]🦀[/orange1] {message}")


# Progress tracking
def create_progress() -> Progress:
    """Create a Rich Progress instance for tracking operations.

    Returns
    -------
    Progress
        Configured Progress instance with spinner and time elapsed.
    """
    return Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        TimeElapsedColumn(),
        console=console,
        disable=is_quiet(),
    )


def create_spinner_progress() -> Progress:
    """Create a simple spinner Progress (no bar, for indeterminate tasks).

    Returns
    -------
    Progress
        Configured Progress instance with spinner only.
    """
    return Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        TimeElapsedColumn(),
        console=console,
        disable=is_quiet(),
    )


# Table formatting
def create_table(title: str | None = None, **kwargs: Any) -> Table:
    """Create a Rich Table for displaying results.

    Parameters
    ----------
    title : str | None
        Optional table title.
    **kwargs
        Additional arguments passed to Table constructor.

    Returns
    -------
    Table
        Configured Table instance.
    """
    return Table(title=title, show_header=True, header_style="bold cyan", **kwargs)


def print_table(table: Table) -> None:
    """Print a Rich Table, respecting verbosity."""
    if is_quiet():
        return
    console.print(table)


# Context manager for timed operations
class TimedOperation:
    """Context manager for timing operations with status output.

    Usage
    -----
    >>> with TimedOperation("Processing chromosomes"):
    ...     # do work
    ... # prints: ✓ Processing chromosomes (2.34s)
    """

    def __init__(self, description: str, show_spinner: bool = True) -> None:
        """Initialize timed operation.

        Parameters
        ----------
        description : str
            Description of the operation.
        show_spinner : bool
            Whether to show a spinner during execution.
        """
        self.description = description
        self.show_spinner = show_spinner and not is_quiet()
        self._progress: Progress | None = None
        self._task_id: int | None = None
        self._start_time: float = 0.0

    def __enter__(self) -> TimedOperation:
        """Start the timed operation."""
        import time

        self._start_time = time.perf_counter()
        if self.show_spinner:
            self._progress = create_spinner_progress()
            self._progress.__enter__()
            self._task_id = self._progress.add_task(self.description, total=None)
        return self

    def __exit__(self, exc_type: Any, exc_val: Any, exc_tb: Any) -> None:
        """End the timed operation and print result."""
        import time

        elapsed = time.perf_counter() - self._start_time
        if self._progress is not None:
            self._progress.__exit__(exc_type, exc_val, exc_tb)

        if exc_type is None:
            success(f"{self.description} ({elapsed:.2f}s)")
        else:
            error(f"{self.description} failed ({elapsed:.2f}s)")


def print_file_path(label: str, path: str) -> None:
    """Print a file path with label.

    Parameters
    ----------
    label : str
        Label describing the file.
    path : str
        Path to the file.
    """
    if is_quiet():
        return
    console.print(f"  [dim]{label}:[/dim] [cyan]{path}[/cyan]")


def print_summary(title: str, items: dict[str, Any]) -> None:
    """Print a summary with key-value pairs.

    Parameters
    ----------
    title : str
        Summary title.
    items : dict[str, Any]
        Key-value pairs to display.
    """
    if is_quiet():
        return
    console.print(f"\n[bold]{title}[/bold]")
    for key, value in items.items():
        console.print(f"  [dim]{key}:[/dim] {value}")


# Version info helper
def print_version_info(
    version: str,
    python_version: str,
    dependencies: dict[str, str],
    rust_available: bool,
) -> None:
    """Print version information in a formatted table.

    Parameters
    ----------
    version : str
        WASP2 version.
    python_version : str
        Python version.
    dependencies : dict[str, str]
        Dictionary of dependency names to versions.
    rust_available : bool
        Whether Rust backend is available.
    """
    table = create_table(title="WASP2 Version Information")
    table.add_column("Component", style="cyan")
    table.add_column("Version", style="green")

    table.add_row("WASP2", version)
    table.add_row("Python", python_version)

    for dep, ver in dependencies.items():
        table.add_row(dep, ver)

    rust_status_str = "[green]available[/green]" if rust_available else "[red]not available[/red]"
    table.add_row("Rust backend", rust_status_str)

    console.print(table)


# CLI callbacks for Typer apps
def verbosity_callback(verbose: bool, quiet: bool) -> None:
    """Set verbosity level based on CLI flags.

    Parameters
    ----------
    verbose : bool
        Enable verbose output.
    quiet : bool
        Suppress all output except errors.
    """
    if quiet:
        set_verbosity(Verbosity.QUIET)
    elif verbose:
        set_verbosity(Verbosity.VERBOSE)
    else:
        set_verbosity(Verbosity.NORMAL)


def version_callback(
    value: bool,
    extra_deps: dict[str, str] | None = None,
) -> None:
    """Show version information and exit.

    Parameters
    ----------
    value : bool
        Whether version flag was passed.
    extra_deps : dict[str, str] | None
        Additional dependencies to show (e.g., {"pysam": pysam.__version__}).
    """
    if not value:
        return

    import rich
    import typer

    from wasp2 import __version__

    deps = {
        "rich": rich.__version__,
        "typer": typer.__version__,
    }
    if extra_deps:
        deps.update(extra_deps)

    try:
        from wasp2_rust import __version__ as rust_version

        rust_available = True
        deps["wasp2_rust"] = rust_version
    except ImportError:
        rust_available = False

    print_version_info(
        version=__version__,
        python_version=f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}",
        dependencies=deps,
        rust_available=rust_available,
    )
    raise typer.Exit()


def create_version_callback(
    extra_deps_func: Callable[[], dict[str, str]] | None = None,
) -> Callable[[bool], None]:
    """Create a version callback with optional extra dependencies.

    Parameters
    ----------
    extra_deps_func : Callable[[], dict[str, str]] | None
        Function that returns extra dependencies to display.
        Called lazily only when --version is used.

    Returns
    -------
    Callable[[bool], None]
        Version callback function for Typer.
    """

    def _callback(value: bool) -> None:
        extra_deps = extra_deps_func() if extra_deps_func and value else None
        version_callback(value, extra_deps)

    return _callback
