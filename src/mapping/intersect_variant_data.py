import timeit
import subprocess
from pathlib import Path
from typing import Optional, List, Union

import numpy as np
import polars as pl

import pysam
from pysam.libcalignmentfile import AlignmentFile

from pybedtools import BedTool

# Import from new wasp2.io module for multi-format support
from wasp2.io import variants_to_bed as _variants_to_bed


def vcf_to_bed(
    vcf_file: Union[str, Path],
    out_bed: Union[str, Path],
    samples: Optional[List[str]] = None
) -> str:
    """Convert variant file to BED format.

    Supports VCF, VCF.GZ, BCF, and PGEN formats via the VariantSource API.

    Note: Parameter name 'vcf_file' is kept for backward compatibility,
    but accepts any supported variant format (VCF, BCF, PGEN).

    Args:
        vcf_file: Path to variant file (VCF, VCF.GZ, BCF, or PGEN)
        out_bed: Output BED file path
        samples: Optional list of sample IDs. If provided, filters to het sites.

    Returns:
        Path to output BED file as string
    """
    # Use new unified interface
    # include_gt=True for mapping (needs genotypes for allele assignment)
    result = _variants_to_bed(
        variant_file=vcf_file,
        out_bed=out_bed,
        samples=samples,
        include_gt=True,
        het_only=True if samples else False,
    )
    return str(result)

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