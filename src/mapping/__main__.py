from typing import Annotated

import typer

from wasp2.cli import create_version_callback, verbosity_callback

from .run_mapping import run_make_remap_reads, run_wasp_filt


def _get_mapping_deps() -> dict[str, str]:
    """Get mapping-specific dependency versions."""
    import polars
    import pysam

    return {"polars": polars.__version__, "pysam": pysam.__version__}


_version_callback = create_version_callback(_get_mapping_deps)

app = typer.Typer(
    pretty_exceptions_short=False,
    rich_markup_mode="rich",
    help="[bold]WASP2 Mapping[/bold] - Generate and filter remapped reads for allele-specific analysis.",
    epilog="[dim]Example: wasp2-map make-reads sample.bam variants.vcf.gz -o remap_dir/[/dim]",
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
    """WASP2 read mapping commands for allele swapping and filtering."""
    verbosity_callback(verbose, quiet)


@app.command()
def make_reads(
    bam: Annotated[str, typer.Argument(help="BAM file")],
    variants: Annotated[str, typer.Argument(help="Variant file (VCF, VCF.GZ, BCF, or PGEN)")],
    samples: Annotated[
        list[str] | None,
        typer.Option(
            "--samples",
            "--sample",
            "--samps",
            "-s",
            help=(
                "One or more samples to use in variant file. "
                "Accepts comma delimited string, or file with one sample per line"
            ),
        ),
    ] = None,
    out_dir: Annotated[
        str | None,
        typer.Option(
            "--out_dir", "--outdir", "--out", "-o", help="Output directory for data to be remapped"
        ),
    ] = None,
    temp_loc: Annotated[
        str | None,
        typer.Option(
            "--temp_loc",
            "--temp",
            "-t",
            help=(
                "Directory for keeping intermediary files."
                "Defaults to removing intermediary files using temp directory"
            ),
        ),
    ] = None,
    out_json: Annotated[
        str | None,
        typer.Option(
            "--out_json",
            "--json",
            "--outjson",
            "-j",
            help=(
                "Output json containing wasp file info to this file instead of default. "
                "Defaults to [BAM_PREFIX]_wasp_data_files.json"
            ),
        ),
    ] = None,
    is_paired: Annotated[
        bool | None,
        typer.Option(
            "--paired/--single",
            help="Reads are paired or single. Will autoparse by default (SINGLE END NOT SUPPORTED YET)",
        ),
    ] = None,
    is_phased: Annotated[
        bool | None,
        typer.Option(
            "--phased/--unphased",
            help=(
                "If variant file is phased/unphased. Will autoparse by default "
                "(PHASED STRONGLY RECOMMENDED-SINGLE END NOT SUPPORTED YET)"
            ),
        ),
    ] = None,
    include_indels: Annotated[
        bool,
        typer.Option(
            "--indels/--snps-only",
            help=(
                "Include indels in addition to SNPs. "
                "Default is SNPs only for backward compatibility. Indel support uses variable-length approach."
            ),
        ),
    ] = False,
    max_indel_len: Annotated[
        int,
        typer.Option(
            "--max-indel-len",
            help="Maximum indel length to process (bp). Indels longer than this are skipped.",
            min=1,
        ),
    ] = 10,
    insert_qual: Annotated[
        int,
        typer.Option(
            "--insert-qual",
            help="Quality score for inserted bases (Phred scale). Used when creating alternate reads.",
            min=0,
            max=60,
        ),
    ] = 30,
    max_seqs: Annotated[
        int,
        typer.Option(
            "--max-seqs",
            help="Maximum number of alternate sequences per read. Reads with more variants are skipped.",
            min=1,
        ),
    ] = 64,
    threads: Annotated[
        int, typer.Option("--threads", help="Threads for BAM I/O operations", min=1)
    ] = 1,
) -> None:
    """Generate reads with swapped alleles for remapping."""
    sample_str = samples[0] if samples else None
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
        threads=threads,
    )


@app.command()
def filter_remapped(
    remapped_bam: Annotated[str, typer.Argument(help="remapped BAM File")],
    to_remap_bam: Annotated[
        str | None, typer.Argument(help="to_remap_bam used to generate swapped alleles")
    ] = None,
    keep_bam: Annotated[
        str | None, typer.Argument(help="BAM containing reads that were not remapped")
    ] = None,
    wasp_data_json: Annotated[
        str | None,
        typer.Option(
            "--wasp_data_json",
            "--json",
            "-j",
            help="json containing wasp file info to load to_remap_bam and keep_bam",
        ),
    ] = None,
    out_bam: Annotated[
        str | None,
        typer.Option(
            "--out_bam",
            "--outbam",
            "--out",
            "-o",
            help="File to output filt bam. Will be created in default name and loc if not provided",
        ),
    ] = None,
    remap_keep_bam: Annotated[
        str | None,
        typer.Option(
            "--remap_keep_bam", help="Also output remapped bam file containing kept reads"
        ),
    ] = None,
    remap_keep_file: Annotated[
        str | None,
        typer.Option("--remap_keep_file", help="Also output txt file with kept read names"),
    ] = None,
    threads: Annotated[
        int, typer.Option("--threads", help="Threads for BAM I/O (Rust filter supports >1)", min=1)
    ] = 1,
    use_rust: Annotated[
        bool,
        typer.Option(
            "--use-rust/--no-rust",
            help="Use Rust acceleration if available (respects WASP2_DISABLE_RUST)",
        ),
    ] = True,
    same_locus_slop: Annotated[
        int,
        typer.Option(
            "--same-locus-slop",
            help=(
                "Tolerance (bp) for 'same locus' test. "
                "Allows remapped reads to differ by this many bp. "
                "Use 2-3 for indels to handle micro-homology shifts. Use 0 for strict SNP-only matching."
            ),
            min=0,
        ),
    ] = 0,
) -> None:
    """Filter remapped reads using WASP algorithm."""
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
