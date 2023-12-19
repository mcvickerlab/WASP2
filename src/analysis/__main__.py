from pathlib import Path
from typing import List, Optional
from typing_extensions import Annotated

import typer
import sys

# Local Imports
from run_analysis import run_ai_analysis

# app = typer.Typer()
# app = typer.Typer(pretty_exceptions_show_locals=False)
app = typer.Typer(pretty_exceptions_short=False)


# TODO GOTTA TEST THIS

# What should i name this?
@app.command()
def find_imbalance(
    counts: Annotated[str, typer.Argument(help="Count File")],
    min: Annotated[
        Optional[int],
        typer.Option(
            "--min",
            "--min_count",
            help=("Minimum allele count for measuring imbalance."
                  " (Default: 10)"
                  )
            )
        ] = None,
    out_file: Annotated[
        Optional[str],
        typer.Option(
            "--out_file",
            "--outfile",
            "--output",
            "--out",
            "-o",
            help=(
                "Output file for analysis. "
                "Defaults to ai_results.tsv"
                ),
            )] = None,
    model: Annotated[
        Optional[str],
        typer.Option(
            "-m",
            "--model",
            help=(
                "Model used for measuring optimization parameter when finding imbalance. "
                "Choice of 'single' or 'linear'. "
                "(Default: 'single')"
                ),
            )] = None,
    region_col: Annotated[
        Optional[str],
        typer.Option(
            "--region_col",
            help=(
                "Name of region column for current data..."
                "Choice of 'region' or 'peak' for ATAC-seq. " 
                "Plans for 'genes' for RNA-seq and 'SNP' for per SNP."
                "(Default: Auto-parses if none provided)"
                ),
            )] = None,
):
    
    # Run
    run_ai_analysis(count_file=counts,
                    min_count=min,
                    model=model,
                    out_file=out_file,
                    region_col=region_col)
    
    # TODO TEST CASES FOR TYPER


if __name__ == "__main__":
    root_dir = Path(__file__).parent
    sys.path.append(str(root_dir))
    app()