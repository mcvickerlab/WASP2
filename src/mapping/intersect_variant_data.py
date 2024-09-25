import timeit
import subprocess
from pathlib import Path

import numpy as np
import polars as pl

import pysam
from pysam.libcalignmentfile import AlignmentFile

from pybedtools import BedTool

def vcf_to_bed(vcf_file, out_bed, samples=None):
    
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


def process_bam(bam_file, vcf_bed, remap_bam, remap_reads, keep_bam, threads=0, is_paired=True):

    # TODO set is_paired to None, and auto check paired vs single
    # print("Filtering reads that overlap regions of interest")

    
    subprocess.run(
        [
            "samtools", "view", "-@", str(threads),
            "-F", "4", "-L", str(vcf_bed), 
            "-o", str(remap_bam), str(bam_file)
        ],
        check=True)

    if is_paired:
        
        # Not needed...but suppresses warning
        subprocess.run(
            ["samtools", "index", "-@", str(threads), str(remap_bam)],
            check=True)
        
        # Extract reads names that overlap het snps
        with AlignmentFile(remap_bam, "rb") as bam, open(remap_reads, "w") as file:
            
            # unique_reads = np.unique(
            #     [read.query_name for read in bam.fetch(until_eof=True)])
            
            # np unique actually kinda slow, use a set instead
            unique_reads = {read.query_name for read in bam.fetch(until_eof=True)}
            
            
            file.write("\n".join(unique_reads))


        # Extract all pairs using read names
        # TODO MIGHT HAVE TO USE INDEXED READ NAMES AND STUFF
        subprocess.run(
            [
                "samtools", "view", "-@", str(threads),
                "-N", remap_reads, "-o", remap_bam, "-U", keep_bam,
                str(bam_file)
                ],
            check=True)

    subprocess.run(
        ["samtools", "sort", "-@", str(threads), "-o", remap_bam, remap_bam],
        check=True)

    subprocess.run(
        ["samtools", "index", "-@", str(threads), str(remap_bam)],
        check=True)

    return remap_bam


# def process_bam(bam_file, vcf_bed, remap_bam, remap_reads, keep_bam, is_paired=True):

#     # TODO set is_paired to None, and auto check paired vs single
#     # print("Filtering reads that overlap regions of interest")
#     pysam.view("-F", "4", "-L", str(vcf_bed), "-o",
#                remap_bam, str(bam_file), catch_stdout=False)

#     if is_paired:
#         # Not needed...but suppresses warning
#         pysam.index(str(remap_bam), catch_stdout=False)

#         # Extract reads names that overlap het snps

#         with AlignmentFile(remap_bam, "rb") as bam, open(remap_reads, "w") as file:
#             unique_reads = np.unique(
#                 [read.query_name for read in bam.fetch(until_eof=True)])
#             file.write("\n".join(unique_reads))

#         # Extract all pairs using read names
#         pysam.view("-N", remap_reads, "-o", remap_bam, "-U", keep_bam,
#                    str(bam_file), catch_stdout=False)
        

#     pysam.sort(remap_bam, "-o", remap_bam, catch_stdout=False)
#     pysam.index(remap_bam, catch_stdout=False)

#     # print("BAM file filtered!")
#     return remap_bam


def intersect_reads(remap_bam, vcf_bed, out_bed):
    # Create Intersections
    a = BedTool(remap_bam)
    b = BedTool(vcf_bed)

    # out_bed = str(Path(out_dir) / "intersect.bed")

    # Perform intersections
    
    # sorted might be necessary...
    a.intersect(b, wb=True, bed=True, sorted=True, output=str(out_bed))
    # a.intersect(b, wb=True, bed=True, sorted=False, output=str(out_bed))

    # print("Created Intersection File")

    return out_bed


def make_intersect_df(intersect_file, samples, is_paired=True):
    
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