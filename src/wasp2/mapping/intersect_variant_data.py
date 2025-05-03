import timeit
import subprocess
from pathlib import Path
from typing import Optional, List

import numpy as np
import polars as pl

import pysam
from pysam.libcalignmentfile import AlignmentFile

from pybedtools import BedTool


def vcf_to_bed(vcf_file: str, out_bed: str, samples: Optional[List[str]] = None) -> str:
    """
    Convert a VCF file to BED format using bcftools.

    This function builds and executes bcftools commands to filter SNPs from a VCF file
    and output the results in BED format. The command options depend on whether sample
    filtering is requested.

    Parameters
    ----------
    vcf_file : str
        Path to the input VCF file.
    out_bed : str
        Path to the output BED file.
    samples : list of str, optional
        List of sample names to filter on. If None, genotypes are dropped.

    Returns
    -------
    str
        The path to the output BED file.
    """
    # Maybe change this later?
    # out_bed = f"{out_dir}/filt_variants.bed"
    
    # Base commands
    view_cmd = [
        "bcftools", "view", str(vcf_file),
        "-m2", "-M2", "-v", "snps", "-Ou"
    ]

    query_cmd = [
        "bcftools", "query",
        "-o", str(out_bed),
        "-f"
    ]
    
    # Parse based on number of samples
    if samples is None:
        # 0 samps, no GTs
        view_cmd.append("--drop-genotypes")
        query_cmd.append("%CHROM\t%POS0\t%END\t%REF\t%ALT\n")
        
        view_process = subprocess.run(view_cmd, stdout=subprocess.PIPE, check=True)
        
    else:
        # Samples
        samples_arg = ",".join(samples)
        num_samples = len(samples)
        
        if num_samples > 1:
            # Multisamp
            view_cmd.extend([
                "-s", samples_arg,
                "--min-ac", "1",
                "--max-ac", str((num_samples * 2) - 1)
            ])
            
            view_process = subprocess.run(view_cmd, stdout=subprocess.PIPE, check=True)
                    
        else:
            # Single Samp subset
            view_cmd.extend(["-s", samples_arg])
            subset_process = subprocess.run(view_cmd, stdout=subprocess.PIPE, check=True)
            
            # Get het genotypes
            new_view_cmd = ["bcftools", "view", "--genotype", "het", "-Ou"]
            view_process = subprocess.run(new_view_cmd, input=subset_process.stdout,
                                          stdout=subprocess.PIPE, check=True)
        
        query_cmd.append("%CHROM\t%POS0\t%END\t%REF\t%ALT[\t%TGT]\n")
    
    # Run Subprocess
    query_process = subprocess.run(query_cmd, input=view_process.stdout, check=True)
    
    return out_bed


# TODO FIX ALL OF THESE TO USE A CLASS
# Process single and pe bam
def process_bam(bam_file: str, vcf_bed: str, remap_bam: str, remap_reads: str, keep_bam: str, is_paired: bool = True) -> str:
    """
    Process a BAM file to filter reads that overlap specified regions.

    This function uses pysam.view to filter reads in the input BAM file (bam_file) that
    overlap regions specified in the BED file (vcf_bed). For paired-end data, it extracts
    unique read names and writes them to a file (remap_reads), and then uses these names to
    create a new BAM file (keep_bam) containing only the desired read pairs. The filtered BAM
    (remap_bam) is then sorted and indexed.

    Parameters
    ----------
    bam_file : str
        Path to the original BAM file.
    vcf_bed : str
        Path to the BED file used for filtering.
    remap_bam : str
        Path to the intermediate BAM file containing reads that overlap the BED regions.
    remap_reads : str
        Path to the file where unique read names (that overlap regions) will be written.
    keep_bam : str
        Path to the output BAM file containing only reads with names in the remap_reads file.
    is_paired : bool, optional
        Indicates if the BAM file contains paired-end reads. Defaults to True.

    Returns
    -------
    str
        The path to the sorted and indexed filtered BAM file (remap_bam).

    Notes
    -----
    The following commented-out code is preserved:
    
        # def process_bam(bam_file, vcf_bed, out_dir=None, is_paired=True):
        #     out_bam = str(Path(out_dir) / "to_remap.bam")
        #     ...
        #     # pysam.view("-N", read_file, "-o", out_bam,
        #     #            str(bam_file), catch_stdout=False)
    """
    # TODO set is_paired to None, and auto check paired vs single
    # print("Filtering reads that overlap regions of interest")
    pysam.view("-F", "4", "-L", str(vcf_bed), "-o",
               remap_bam, str(bam_file), catch_stdout=False)

    if is_paired:
        # Not needed...but suppresses warning
        pysam.index(str(remap_bam), catch_stdout=False)

        # Extract reads names that overlap het snps
        with AlignmentFile(remap_bam, "rb") as bam, open(remap_reads, "w") as file:
            unique_reads = np.unique(
                [read.query_name for read in bam.fetch(until_eof=True)]
            )
            file.write("\n".join(unique_reads))

        # Extract all pairs using read names
        pysam.view("-N", remap_reads, "-o", remap_bam, "-U", keep_bam,
                   str(bam_file), catch_stdout=False)
        
    pysam.sort(remap_bam, "-o", remap_bam, catch_stdout=False)
    pysam.index(remap_bam, catch_stdout=False)

    # print("BAM file filtered!")
    return remap_bam


# def process_bam(bam_file, vcf_bed, out_dir=None, is_paired=True):
#     out_bam = str(Path(out_dir) / "to_remap.bam")
#
#     # TODO set is_paired to None, and auto check paired vs single
#     # print("Filtering reads that overlap regions of interest")
#     pysam.view("-F", "4", "-L", str(vcf_bed), "-o",
#                out_bam, str(bam_file), catch_stdout=False)
#
#     if is_paired:
#         # Not needed...but suppresses warning
#         pysam.index(str(out_bam), catch_stdout=False)
#
#         # Extract reads names that overlap het snps
#         read_file = str(Path(out_dir) / "to_remap.txt")
#
#         with AlignmentFile(out_bam, "rb") as bam, open(read_file, "w") as file:
#             unique_reads = np.unique(
#                 [read.query_name for read in bam.fetch(until_eof=True)]
#             )
#             file.write("\n".join(unique_reads))
#
#         # Extract all pairs using read names
#         keep_bam = str(Path(out_dir) / "keep.bam")
#         pysam.view("-N", read_file, "-o", out_bam, "-U", keep_bam,
#                    str(bam_file), catch_stdout=False)
#
#         # pysam.view("-N", read_file, "-o", out_bam,
#         #            str(bam_file), catch_stdout=False)
#
#     pysam.sort(out_bam, "-o", out_bam, catch_stdout=False)
#     pysam.index(out_bam, catch_stdout=False)
#
#     # print("BAM file filtered!")
#     return out_bam


def intersect_reads(remap_bam: str, vcf_bed: str, out_bed: str) -> str:
    """
    Intersect a BAM file with a BED file using bedtools.

    This function uses pybedtools to intersect reads from a BAM file (remap_bam) with regions
    specified in a BED file (vcf_bed), and writes the intersection results to an output BED file.

    Parameters
    ----------
    remap_bam : str
        Path to the remapped BAM file.
    vcf_bed : str
        Path to the BED file (derived from VCF) to use for intersection.
    out_bed : str
        Path to the output BED file where the intersection results will be written.

    Returns
    -------
    str
        The path to the output intersection BED file.
    """
    # Create Intersections
    a = BedTool(remap_bam)
    b = BedTool(vcf_bed)

    # out_bed = str(Path(out_dir) / "intersect.bed")
    
    # Perform intersections
    # a.intersect(b, wb=True, bed=True, sorted=True, output=str(out_bed))
    a.intersect(b, wb=True, bed=True, sorted=False, output=str(out_bed))

    # print("Created Intersection File")
    return out_bed


def make_intersect_df(intersect_file: str, samples: List[str], is_paired: bool = True) -> pl.DataFrame:
    """
    Create a DataFrame from an intersection file with sample-specific allele data.

    This function reads an intersection file using Polars, processes the data by selecting and renaming
    columns, and splits genotype information for each sample. It also splits mate information from the
    read names, casts the mate column to UInt8, and removes duplicate entries.

    Parameters
    ----------
    intersect_file : str
        Path to the intersection file.
    samples : list of str
        List of sample names corresponding to genotype data.
    is_paired : bool, optional
        Flag indicating if the data comes from paired-end sequencing. Defaults to True.

    Returns
    -------
    polars.DataFrame
        A unique, processed DataFrame containing columns "chrom", "read", "start", "stop", and genotype
        columns for each sample, along with split mate information.
    
    Notes
    -----
    The following commented-out code is preserved:
    
        # def make_intersect_df(intersect_file, samples, is_paired=True):
        #     # Create Dataframe
        #     df = pl.scan_csv(intersect_file, separator="\t", has_header=False)
        #     ...
        #     # df = df.unique() # Remove possible dups
    """
    # Create Dataframe
    df = pl.scan_csv(intersect_file,
                     separator="\t",
                     has_header=False,
                     infer_schema_length=0)
    
    # Parse sample data
    num_samps = len(samples)
    
    subset_cols = [df.columns[i] for i in np.r_[0, 3, 1, 2, -num_samps:0]]
    new_cols = ["chrom", "read", "start", "stop", *samples]
    
    rename_cols = {old_col: new_col for old_col, new_col in zip(subset_cols, new_cols)}
    
    base_schema = [
        pl.col("chrom").cast(pl.Categorical),
        pl.col("read").cast(pl.Utf8),
        pl.col("start").cast(pl.UInt32),
        pl.col("stop").cast(pl.UInt32)
    ]
    
    sample_schema = [pl.col(samp).cast(pl.Utf8) for samp in samples]
    col_schema = [*base_schema, *sample_schema]
    
    # Make sure types are correct
    df = df.select(subset_cols).rename(rename_cols).with_columns(col_schema)

    expr_list = []
    cast_list = []
    
    for s in samples:
        a1 = f"{s}_a1"
        a2 = f"{s}_a2"
        # Add split per sample
        expr_list.append(
            pl.col(s).str.split_exact(by="|", n=1).struct.rename_fields([a1, a2])
        )
        # Cast new genotype columns
        cast_list.append(pl.col(a1).cast(pl.Categorical))
        cast_list.append(pl.col(a2).cast(pl.Categorical))
    
    # Split mate expr
    expr_list.append(
        pl.col("read").str.split_exact(by="/", n=1).struct.rename_fields(["read", "mate"])
    )
    
    cast_list.append(pl.col("mate").cast(pl.UInt8))
    
    df = df.with_columns(expr_list).unnest([*samples, "read"]).with_columns(cast_list)
    
    # Should I remove instead of keep first?
    df = df.unique(["chrom", "read", "mate", "start", "stop"], keep="first")  # Doesnt remove dup snp in pair?
    
    return df.collect()
