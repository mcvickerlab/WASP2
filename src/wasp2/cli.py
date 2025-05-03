from pathlib import Path
from typing import List, Optional
from typing_extensions import Annotated

import typer
import sys

# Local Imports
from wasp2.analysis.__main__ import app as analysis_app
from wasp2.counting.__main__ import app as counting_app
from wasp2.mapping.__main__ import app as mapping_app

# Create a Typer app instance with a brief description.
cli: typer.Typer = typer.Typer(help="WASP2: Toolkit for allele-specific analysis.")

# Register subcommands (command groups) with their respective Typer app objects.
cli.add_typer(counting_app, name="counting", help="Commands for generating allele-specific read counts.")
cli.add_typer(analysis_app, name="analysis", help="Commands for analyzing allelic imbalance.")
cli.add_typer(mapping_app, name="mapping", help="Commands for unbiased read mapping and remapping.")

def main() -> None:
    """
    Entry point for the WASP2 CLI.

    This function initializes the WASP2 Typer application and dispatches commands
    to the appropriate subcommand group.

    **Usage Examples:**
      - Count variants in a BAM file:
        ```
        WASP2 counting count-variants [BAM] [VCF] [OPTIONS]
        ```
      - Analyze allelic imbalance from count data:
        ```
        WASP2 analysis find-imbalance [COUNTS] [OPTIONS]
        ```
      - Run the mapping pipeline:
        ```
        WASP2 mapping make-reads [BAM] [VCF] [OPTIONS]
        ```

    Returns:
        None
    """
    cli()

if __name__ == "__main__":
    main()
