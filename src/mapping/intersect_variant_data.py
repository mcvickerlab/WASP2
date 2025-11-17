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
    
    # Maybe change this later?
    # out_bed = f"{out_dir}/filt_variants.bed"
    
    # Base commands
    view_cmd = ["bcftools", "view", str(vcf_file),
                "-m2", "-M2", "-v", "snps", "-Ou"
               ]

    query_cmd = ["bcftools", "query",
                 "-o", str(out_bed),
                 "-f"]
    
    # Parse based on num samps
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
            view_cmd.extend(["-s", samples_arg,
                             "--min-ac", "1",
                             "--max-ac", str((num_samples * 2) - 1)])
            
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
def process_bam(
    bam_file: str,
    vcf_bed: str,
    remap_bam: str,
    remap_reads: str,
    keep_bam: str,
    is_paired: bool = True
) -> str:

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
                [read.query_name for read in bam.fetch(until_eof=True)])
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

#     # TODO set is_paired to None, and auto check paired vs single
#     # print("Filtering reads that overlap regions of interest")
#     pysam.view("-F", "4", "-L", str(vcf_bed), "-o",
#                out_bam, str(bam_file), catch_stdout=False)

#     if is_paired:
#         # Not needed...but suppresses warning
#         pysam.index(str(out_bam), catch_stdout=False)

#         # Extract reads names that overlap het snps
#         read_file = str(Path(out_dir) / "to_remap.txt")

#         with AlignmentFile(out_bam, "rb") as bam, open(read_file, "w") as file:
#             unique_reads = np.unique(
#                 [read.query_name for read in bam.fetch(until_eof=True)])
#             file.write("\n".join(unique_reads))

#         # Extract all pairs using read names
#         keep_bam = str(Path(out_dir) / "keep.bam")
#         pysam.view("-N", read_file, "-o", out_bam, "-U", keep_bam,
#                    str(bam_file), catch_stdout=False)
        
#         # pysam.view("-N", read_file, "-o", out_bam,
#         #            str(bam_file), catch_stdout=False)
        

#     pysam.sort(out_bam, "-o", out_bam, catch_stdout=False)
#     pysam.index(out_bam, catch_stdout=False)

#     # print("BAM file filtered!")
#     return out_bam


def intersect_reads(remap_bam: str, vcf_bed: str, out_bed: str) -> str:
    # Create Intersections
    a = BedTool(remap_bam)
    b = BedTool(vcf_bed)

    # out_bed = str(Path(out_dir) / "intersect.bed")

    # Perform intersections
    # a.intersect(b, wb=True, bed=True, sorted=True, output=str(out_bed))
    a.intersect(b, wb=True, bed=True, sorted=False, output=str(out_bed))

    # print("Created Intersection File")

    return out_bed


# Probs should move this to a method
# def filter_intersect_data(bam_file, vcf_file, out_dir, samples=None, is_paired=True):

#     # Get het snps
#     het_start = timeit.default_timer()

#     het_bed_file = vcf_to_bed(vcf_file, samples, out_dir)
#     # het_bed_file = vcf_to_bed(vcf_file, out_dir)
#     print(f"Finished in {timeit.default_timer() - het_start:.2f} seconds!\n")

#     # Filter bam reads intersecting snps
#     bam_start = timeit.default_timer()

#     het_bam_file = process_bam(
#         bam_file, het_bed_file, out_dir, is_paired=is_paired)
#     print(f"Finished in {timeit.default_timer() - bam_start:.2f} seconds!\n")

#     # Get reads overlapping snps
#     snp_start = timeit.default_timer()

#     read_intersect_file = intersect_reads(
#         het_bam_file, het_bed_file, out_dir)
#     print(f"Finished in {timeit.default_timer() - snp_start:.2f} seconds!\n")

#     return het_bam_file, read_intersect_file


# Should this be here?
# def make_intersect_df(intersect_file, samples, is_paired=True):
    
#     # Create Dataframe
#     df = pl.scan_csv(intersect_file, separator="\t", has_header=False)
    
#     # Parse sample data
#     num_samps = len(samples)
    
#     subset_cols = [df.columns[i] for i in np.r_[0, 3, 1, 2, -num_samps:0]]
#     new_cols = ["chrom", "read", "start", "stop", *samples]
#     rename_cols = {old_col: new_col for old_col, new_col in zip(subset_cols, new_cols)}
    
#     # Make sure types are correct
#     df = df.select(subset_cols).rename(rename_cols).with_columns(
#         [
#             pl.col(col).cast(pl.UInt32) if (col == "start") or (col == "stop")
#             else pl.col(col).cast(pl.Utf8) for col in new_cols
#         ]
#     )
    
#     # TODO CHANGE THESE TO BE A BIT CATEGORICAL
#     # df = df.select(subset_cols).rename(
#     #     rename_cols).with_columns(
#     #         [
#     #             pl.col("chrom").cast(pl.Categorical),
#     #             pl.col("pos").cast(pl.UInt32),
#     #             pl.col("ref").cast(pl.Categorical),
#     #             pl.col("alt").cast(pl.Categorical)
#     #             ]
#     #         )
    
#     # Split sample alleles expr
#     # Maybe don't do this for multi
#     expr_list = [
#         pl.col(s).str.split_exact(
#             by="|", n=1).struct.rename_fields([f"{s}_a1", f"{s}_a2"])
#         for s in df.columns[4:]
#     ]

#     # Split mate expr
#     expr_list.append(
#         pl.col("read").str.split_exact(
#             by="/", n=1).struct.rename_fields(["read", "mate"])
#     )


#     df = df.with_columns(expr_list).unnest(
#         [*df.columns[4:], "read"]).with_columns(
#         pl.col("mate").cast(pl.UInt8))

#     # df = df.unique() # Remove possible dups
#     # should i remove instead of keep first?
#     # df = df.unique(["chrom", "read", "start", "stop"], keep="first") # Remove dup snps
#     df = df.unique(["chrom", "read", "mate", "start", "stop"], keep="first") # Doesnt remove dup snp in pair?
#     df = df.collect()
    
#     return df


def make_intersect_df(intersect_file: str, samples: List[str], is_paired: bool = True) -> pl.DataFrame:

    # Create Dataframe
    df = pl.scan_csv(intersect_file,
                     separator="\t",
                     has_header=False,
                     infer_schema_length=0
                    )
    
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
            pl.col(s).str.split_exact(
                by="|", n=1).struct.rename_fields([a1, a2])
        )
        
        # cast new gt cols
        cast_list.append(pl.col(a1).cast(pl.Categorical))
        cast_list.append(pl.col(a2).cast(pl.Categorical))

    # Split mate expr
    expr_list.append(
        pl.col("read").str.split_exact(
            by="/", n=1).struct.rename_fields(["read", "mate"])
    )
    
    cast_list.append(pl.col("mate").cast(pl.UInt8))
    
    df = df.with_columns(expr_list).unnest(
        [*samples, "read"]).with_columns(
        cast_list
    )


    # should i remove instead of keep first?
    df = df.unique(["chrom", "read", "mate", "start", "stop"], keep="first") # Doesnt remove dup snp in pair?
    
    return df.collect()