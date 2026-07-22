from enum import Enum
from typing import Annotated

import typer

from wasp2.cli import create_version_callback, verbosity_callback

from .run_analysis import run_ai_analysis


class AnalysisScopeChoice(str, Enum):
    per_donor = "per-donor"


class AnalysisUnitChoice(str, Enum):
    snv = "snv"
    feature = "feature"
    peak = "peak"


class AnalysisModelChoice(str, Enum):
    single = "single"
    linear = "linear"


class DispersionScopeChoice(str, Enum):
    per_donor = "per-donor"
    global_ = "global"


def _get_analysis_deps() -> dict[str, str]:
    """Get analysis-specific dependency versions."""
    import polars
    import scipy

    return {"polars": polars.__version__, "scipy": scipy.__version__}


_version_callback = create_version_callback(_get_analysis_deps)

app = typer.Typer(
    pretty_exceptions_short=False,
    rich_markup_mode="rich",
    help="[bold]WASP2 Analysis[/bold] - Detect and compare allelic imbalance.",
    epilog="[dim]Example: wasp2-analyze find-imbalance counts.tsv -o results.tsv[/dim]",
)


@app.callback(invoke_without_command=True)
def main(
    ctx: typer.Context,
    version: Annotated[
        bool,
        typer.Option(
            "--version",
            "-V",
            callback=_version_callback,
            is_eager=True,
            help="Show version and dependency information.",
        ),
    ] = False,
    verbose: Annotated[
        bool,
        typer.Option("--verbose", "-v", help="Enable verbose output with detailed progress."),
    ] = False,
    quiet: Annotated[
        bool,
        typer.Option("--quiet", "-q", help="Suppress all output except errors."),
    ] = False,
) -> None:
    """WASP2 allelic imbalance analysis commands."""
    verbosity_callback(verbose, quiet)


@app.command()
def find_imbalance(
    counts: Annotated[
        str,
        typer.Argument(help="Count table, locked count bundle, or count_manifest.json."),
    ],
    min: Annotated[
        int | None,
        typer.Option(
            "--min",
            "--min_count",
            help="Minimum allele count for measuring imbalance. (Default: 10)",
        ),
    ] = None,
    pseudocount: Annotated[
        int | None,
        typer.Option(
            "-p",
            "--ps",
            "--pseudo",
            "--pseudocount",
            help=(
                "Pseudocount added when measuring allelic imbalance. Defaults to 0 for "
                "--scope per-donor and 1 for the legacy route."
            ),
        ),
    ] = None,
    out_file: Annotated[
        str | None,
        typer.Option(
            "--out_file",
            "--outfile",
            "--output",
            "--out",
            "-o",
            help="Output file for analysis. Defaults to ai_results.tsv",
        ),
    ] = None,
    phased: Annotated[
        bool | None,
        typer.Option(
            "--phased",
            help=(
                "Calculate allelic imbalance using the phased haplotype model. "
                "Genotype info must phased and included in allelic count data!"
                "\nBy default, calculates unphased AI assuming equal liklihood for each haplotype."
            ),
        ),
    ] = False,
    model: Annotated[
        AnalysisModelChoice | None,
        typer.Option(
            "-m",
            "--model",
            help=("Dispersion model used for imbalance inference. Defaults to 'single'."),
        ),
    ] = None,
    region_col: Annotated[
        str | None,
        typer.Option(
            "--region-col",
            "--region_col",
            help=(
                "Feature identifier column. Locked feature bundles use 'region'; standalone "
                "feature tables require this option."
            ),
        ),
    ] = None,
    groupby: Annotated[
        str | None,
        typer.Option(
            "--groupby",
            "--group",
            "--parent_col",
            "--parent",
            help=(
                "Report allelic imbalance by parent group instead of feature level in RNA-seq counts. "
                "Name of parent column. Not valid if no parent column or if using ATAC-seq peaks. "
                "(Default: Report by feature level instead of parent level)"
            ),
        ),
    ] = None,
    per_variant: Annotated[
        bool,
        typer.Option(
            "--per-variant",
            "--snv-solo",
            "--per_variant",
            help=(
                "Test each SNV independently (per-variant) instead of grouping by region/peak. "
                "Forces per-variant even when a region column is present. "
                "Cannot be combined with --region_col."
            ),
        ),
    ] = False,
    scope: Annotated[
        AnalysisScopeChoice | None,
        typer.Option(
            "--scope",
            help=(
                "Analyze each donor independently. Requires --unit and canonical donor_id; "
                "standalone legacy tables may use sample instead."
            ),
        ),
    ] = None,
    unit: Annotated[
        AnalysisUnitChoice | None,
        typer.Option(
            "--unit",
            help="Statistical unit: snv or feature. peak is an alias for feature.",
        ),
    ] = None,
    dispersion_scope: Annotated[
        DispersionScopeChoice,
        typer.Option(
            "--dispersion-scope",
            help="Fit dispersion globally across eligible donor-SNV rows or per donor.",
        ),
    ] = DispersionScopeChoice.per_donor,
    min_donor_observations: Annotated[
        int,
        typer.Option(
            "--min-donor-observations",
            min=1,
            help="Exclude donors with fewer eligible unique SNV observations.",
        ),
    ] = 50,
    expected_manifest_sha256: Annotated[
        str | None,
        typer.Option(
            "--expected-manifest-sha256",
            help="Require a locked bundle manifest with this SHA-256 digest.",
        ),
    ] = None,
) -> None:
    run_ai_analysis(
        count_file=counts,
        min_count=min,
        model=model.value if model is not None else None,
        pseudocount=pseudocount,
        phased=phased,
        out_file=out_file,
        region_col=region_col,
        groupby=groupby,
        per_variant=per_variant,
        scope=scope.value if scope is not None else None,
        unit=unit.value if unit is not None else None,
        dispersion_scope=dispersion_scope.value,
        min_donor_observations=min_donor_observations,
        expected_manifest_sha256=expected_manifest_sha256,
    )


@app.command()
def find_imbalance_sc(
    counts: Annotated[str, typer.Argument(help="Count File")],
    bc_map: Annotated[
        str,
        typer.Argument(
            help="Two Column TSV file mapping specific barcodes to some grouping/celltype. Each line following format [BARCODE]\\t[GROUP]"
        ),
    ],
    min: Annotated[
        int | None,
        typer.Option(
            "--min",
            "--min_count",
            help="Minimum allele count per region for measuring imbalance. (Default: 10)",
        ),
    ] = None,
    pseudocount: Annotated[
        int | None,
        typer.Option(
            "-p",
            "--ps",
            "--pseudo",
            "--pseudocount",
            help="Pseudocount added when measuring allelic imbalance. (Default: 1)",
        ),
    ] = None,
    sample: Annotated[
        str | None,
        typer.Option(
            "--sample",
            "--samp",
            "-s",
            help="Use heterozygous genotypes for this sample in count file. Automatically parses if data contains 0 or 1 sample. REQUIRED IF COUNT DATA CONTAINS MULTIPLE SAMPLES.",
        ),
    ] = None,
    groups: Annotated[
        list[str] | None,
        typer.Option(
            "--groups",
            "--group",
            "--celltypes",
            "--g",
            help="Specific groups in barcode file/bc_map to analyze allelic imbalance in. Uses all groups in barcode file/bc_map by default.",
        ),
    ] = None,
    phased: Annotated[
        bool | None,
        typer.Option(
            "--phased/--unphased",
            help="If genotypes are phased use phasing information to measure imbalance. Otherwise assume all haplotypes are equally likely. Autoparses genotype data by default.",
        ),
    ] = None,
    out_file: Annotated[
        str | None,
        typer.Option(
            "--out_file",
            "--outfile",
            "--output",
            "--out",
            "-o",
            help="Output file for analysis. Defaults to ai_results_[GROUP].tsv",
        ),
    ] = None,
    z_cutoff: Annotated[
        int | None,
        typer.Option(
            "-z",
            "--z_cutoff",
            "--zscore_cutoff",
            "--remove_outliers",
            "--remove_extreme",
            "--z_boundary",
            "--zcore_boundary",
            help="Remove SNPs and associated regions whose counts exceed Z-Score cutoff. (Default: None)",
        ),
    ] = None,
) -> None:
    try:
        from .run_analysis_sc import run_ai_analysis_sc
    except (ImportError, ModuleNotFoundError) as error:
        raise typer.BadParameter(
            "Single-cell analysis requires: pip install 'wasp2[single-cell]'"
        ) from error

    run_ai_analysis_sc(
        count_file=counts,
        bc_map=bc_map,
        min_count=min,
        pseudocount=pseudocount,
        phase=phased,
        sample=sample,
        groups=groups,
        out_file=out_file,
        z_cutoff=z_cutoff,
    )


@app.command()
def compare_imbalance(
    counts: Annotated[str, typer.Argument(help="Count File")],
    bc_map: Annotated[
        str,
        typer.Argument(
            help="Two Column TSV file mapping specific barcodes to some grouping/celltype. Each line following format [BARCODE]\\t[GROUP]"
        ),
    ],
    min: Annotated[
        int | None,
        typer.Option(
            "--min",
            "--min_count",
            help="Minimum allele count for measuring imbalance. (Default: 10)",
        ),
    ] = None,
    pseudocount: Annotated[
        int | None,
        typer.Option(
            "-p",
            "--ps",
            "--pseudo",
            "--pseudocount",
            help="Pseudocount added when measuring allelic imbalance. (Default: 1)",
        ),
    ] = None,
    sample: Annotated[
        str | None,
        typer.Option(
            "--sample",
            "--samp",
            "-s",
            help="Use heterozygous genotypes for this sample in count file. Automatically parses if data contains 0 or 1 sample. REQUIRED IF COUNT DATA CONTAINS MULTIPLE SAMPLES.",
        ),
    ] = None,
    groups: Annotated[
        list[str] | None,
        typer.Option(
            "--groups",
            "--group",
            "--celltypes",
            "--g",
            help="Specific groups in barcode file/bc_map to compare allelic imbalance between. If providing input, requires a minimum of 2 groups. Uses all group combinations by default.",
        ),
    ] = None,
    phased: Annotated[
        bool | None,
        typer.Option(
            "--phased/--unphased",
            help="If genotypes are phased use phasing information to measure imbalance. Otherwise assume all haplotypes are equally likely. Autoparses genotype data by default.",
        ),
    ] = None,
    out_file: Annotated[
        str | None,
        typer.Option(
            "--out_file",
            "--outfile",
            "--output",
            "--out",
            "-o",
            help="Output file for comparisons. Defaults to ai_results_[GROUP1]_[GROUP2].tsv",
        ),
    ] = None,
    z_cutoff: Annotated[
        int | None,
        typer.Option(
            "-z",
            "--z_cutoff",
            "--zscore_cutoff",
            "--remove_outliers",
            "--remove_extreme",
            "--z_boundary",
            "--zcore_boundary",
            help="Remove SNPs and associated regions whose counts exceed Z-Score cutoff. (Default: None)",
        ),
    ] = None,
) -> None:
    try:
        from .run_compare_ai import run_ai_comparison
    except (ImportError, ModuleNotFoundError) as error:
        raise typer.BadParameter(
            "Single-cell comparison requires: pip install 'wasp2[single-cell]'"
        ) from error

    run_ai_comparison(
        count_file=counts,
        bc_map=bc_map,
        min_count=min,
        pseudocount=pseudocount,
        phase=phased,
        sample=sample,
        groups=groups,
        out_file=out_file,
        z_cutoff=z_cutoff,
    )
