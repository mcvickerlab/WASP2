from pathlib import Path
from typing import List, Optional
from typing_extensions import Annotated

import typer
import sys

# Local Imports
from .run_mapping import run_make_remap_reads, run_wasp_filt


app = typer.Typer()
# app = typer.Typer(pretty_exceptions_show_locals=False)
# app = typer.Typer(pretty_exceptions_short=False)

@app.command()
def make_reads(
    bam: Annotated[str, typer.Argument(help="BAM file")],
    variants: Annotated[str, typer.Argument(help="Variant file (VCF, VCF.GZ, BCF, or PGEN)")],
    samples: Annotated[
        Optional[List[str]],
        typer.Option(
            "--samples",
            "--sample",
            "--samps",
            "-s",
            help=(
                "One or more samples to use in variant file. "
                "Accepts comma delimited string, "
                "or file with one sample per line"
            )
            )] = None,
    out_dir: Annotated[
        Optional[str],
        typer.Option(
            "--out_dir",
            "--outdir",
            "--out",
            "-o",
            help="Output directory for data to be remapped"
            )] = None,
    temp_loc: Annotated[
        Optional[str],
        typer.Option(
            "--temp_loc",
            "--temp",
            "-t",
            help=(
                "Directory for keeping intermediary files."
                "Defaults to removing intermediary files using temp directory")
            )] = None,
    # is_phased: Annotated[Optional[bool], typer.Argument()] = None,
    out_json: Annotated[
        Optional[str],
        typer.Option("--out_json",
                     "--json",
                     "--outjson",
                     "-j",
                     help=(
                         "Output json containing wasp file info to this file instead of default. "
                         "Defaults to [BAM_PREFIX]_wasp_data_files.json"
                     )
                     )] = None,
    is_paired: Annotated[
        Optional[bool],
        typer.Option("--paired/--single",
                     help=(
                         "Reads are paired or single. "
                         "Will autoparse by default "
                         "(SINGLE END NOT SUPPORTED YET)"
                         )
                     )] = None,
    is_phased: Annotated[
        Optional[bool],
        typer.Option("--phased/--unphased",
                     help=(
                         "If variant file is phased/unphased. "
                         "Will autoparse by default "
                         "(PHASED STRONGLY RECOMMENDED-SINGLE END NOT SUPPORTED YET)"
                         )
                     )] = None,
    include_indels: Annotated[
        bool,
        typer.Option("--indels/--snps-only",
                     help=(
                         "Include indels in addition to SNPs. "
                         "Default is SNPs only for backward compatibility. "
                         "Indel support uses variable-length approach."
                         )
                     )] = False,
    max_indel_len: Annotated[
        int,
        typer.Option("--max-indel-len",
                     help=(
                         "Maximum indel length to process (bp). "
                         "Indels longer than this are skipped. "
                         "Prevents excessive computational burden."
                         ),
                     min=1
                     )] = 10,
    insert_qual: Annotated[
        int,
        typer.Option("--insert-qual",
                     help=(
                         "Quality score for inserted bases (Phred scale). "
                         "Used when creating alternate reads with insertions."
                         ),
                     min=0,
                     max=60
                     )] = 30,
    max_seqs: Annotated[
        int,
        typer.Option("--max-seqs",
                     help=(
                         "Maximum number of alternate sequences per read. "
                         "Reads with more variants are skipped. "
                         "Prevents combinatorial explosion."
                         ),
                     min=1
                     )] = 64,
    threads: Annotated[
        int,
        typer.Option(
            "--threads",
            help="Threads for BAM I/O operations",
            min=1
        )
    ] = 1,
) -> None:
    """Generate reads with swapped alleles for remapping."""

    # Parse sample string
    sample_str: Optional[str]
    if samples is not None and len(samples) > 0:
        sample_str = samples[0]
    else:
        sample_str = None

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




@app.command()
def filter_remapped(
    remapped_bam: Annotated[str, typer.Argument(help="remapped BAM File")],
    to_remap_bam: Annotated[
        Optional[str], typer.Argument(
            help="to_remap_bam used to generate swapped alleles")
        ] = None,
    keep_bam: Annotated[
        Optional[str], typer.Argument(
            help="BAM containing reads that were not remapped")
        ] = None,
    wasp_data_json: Annotated[
        Optional[str],
        typer.Option("--wasp_data_json",
                     "--json",
                     "-j",
                     help=(
                         "json containing wasp file info to load to_remap_bam and keep_bam"
                         ))
        ] = None,
    out_bam: Annotated[
        Optional[str],
        typer.Option(
            "--out_bam",
            "--outbam",
            "--out",
            "-o",
            help=(
                "File to output filt bam. "
                "Will be created in default name and loc if not provided"
                )
            )] = None,
    remap_keep_bam: Annotated[
        Optional[str],
        typer.Option(
            "--remap_keep_bam",
            help=(
                "Also output remapped bam file containing kept reads"
                )
            )] = None,
    remap_keep_file: Annotated[
        Optional[str],
        typer.Option(
            "--remap_keep_file",
            help=(
                "Also output txt file with kept read names"
                )
            )] = None,
    threads: Annotated[
        int,
        typer.Option(
            "--threads",
            help="Threads for BAM I/O (Rust filter supports >1)",
            min=1
        )
    ] = 1,
    use_rust: Annotated[
        bool,
        typer.Option(
            "--use-rust/--no-rust",
            help="Use Rust acceleration if available (respects WASP2_DISABLE_RUST)",
        )
    ] = True,
    same_locus_slop: Annotated[
        int,
        typer.Option(
            "--same-locus-slop",
            help=(
                "Tolerance (bp) for 'same locus' test. "
                "Allows remapped reads to differ by this many bp. "
                "Use 2-3 for indels to handle micro-homology shifts. "
                "Use 0 for strict SNP-only matching."
            ),
            min=0
        )
    ] = 0,
) -> None:
    """Filter remapped reads using WASP algorithm."""

    # Checks
    # print(remapped_bam)
    # print(to_remap_bam)
    # print(keep_bam)
    # print(wasp_data_json)
    # print(out_bam)
    # print(remap_keep_bam)
    # print(remap_keep_file)
    
    # Run WASP Filt
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
    

if __name__ == "__main__":
    root_dir = Path(__file__).parent
    sys.path.append(str(root_dir))
    app()
