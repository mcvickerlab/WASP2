import sys
import timeit
import subprocess
import warnings


from pathlib import Path
from typing import Optional, List, Union

import numpy as np
import polars as pl

# Import from new wasp2.io module for multi-format support
from wasp2.io import variants_to_bed as _variants_to_bed


def vcf_to_bed(
    vcf_file: Union[str, Path],
    out_bed: Union[str, Path],
    samples: Optional[List[str]] = None,
    include_gt: bool = True,
    include_indels: bool = False,
    max_indel_len: int = 10
) -> str:
    """Convert variant file to BED format.

    Supports VCF, VCF.GZ, BCF, and PGEN formats via the VariantSource API.
    This is the unified version that replaces the duplicate implementation.

    Note: Parameter name 'vcf_file' is kept for backward compatibility,
    but accepts any supported variant format (VCF, BCF, PGEN).

    Args:
        vcf_file: Path to variant file (VCF, VCF.GZ, BCF, or PGEN)
        out_bed: Output BED file path
        samples: Optional list of sample IDs. If provided, filters to het sites.
        include_gt: Include genotype column in output (default True)
        include_indels: Include indels in addition to SNPs (default False)
        max_indel_len: Maximum indel length in bp (default 10)

    Returns:
        Path to output BED file as string
    """
    # Use new unified interface
    result = _variants_to_bed(
        variant_file=vcf_file,
        out_bed=out_bed,
        samples=samples,
        include_gt=include_gt,
        het_only=True if samples else False,
        include_indels=include_indels,
        max_indel_len=max_indel_len,
    )
    return str(result)


def gtf_to_bed(gtf_file, out_bed, feature, attribute):
    
    # Use gtf col names
    gtf_cols = [
        "seqname", "source", "feature",
        "start", "end", "score",
        "strand", "frame", "attribute"]
    
    
    # Cant use lazyframe in case of compressed
    df = pl.read_csv(gtf_file, separator="\t",
                     comment_prefix="#",
                     has_header=False,
                     new_columns=gtf_cols)
    
    # Extract from attribute col
    attr_regex = fr'{attribute}[=\s]\"?\'?(.*?)\"?\'?;' # works for gtf/gff3
    
    # Extract feature only and attributes
    df = df.filter(pl.col("feature") == feature
                  ).with_columns(
        pl.col("start").sub(1),
        pl.col("attribute").str.extract(attr_regex).alias(attribute)
    ).select(["seqname", "start", "end", attribute])
    
    # TODO Extra validation and may want to return some data?
    
    # Write to BED
    df.write_csv(out_bed, separator="\t", include_header=False)
    
    return out_bed


# Perform intersection
def intersect_vcf_region(vcf_file, region_file, out_file):
    
    # Parse region file before or after???
    intersect_cmd = ["bedtools", "intersect" , "-a", str(vcf_file),
                     "-b", str(region_file), "-wb"]
    

    # write intersect out
    with open(out_file, "w") as file:
        intersect_process = subprocess.run(intersect_cmd, stdout=file, check=True)


# TODO, update old software to use this new version
# Convert Intersect file to df
def parse_intersect_region_new(intersect_file, samples=None, use_region_names=False, region_col=None):
    
    if region_col is None:
        region_col = "region"


    # Default number of columns
    vcf_cols = ["chrom", "pos0", "pos", "ref", "alt"] # Default columns for vcf
    vcf_schema = [pl.Categorical, pl.UInt32, pl.UInt32,
                  pl.Categorical, pl.Categorical]


    if samples is not None:
        vcf_cols.extend(samples)
        vcf_schema.extend([pl.Categorical] * len(samples))


    vcf_ncols = len(vcf_cols)

    # Process with gt
    df = pl.scan_csv(intersect_file, separator="\t",
                     has_header=False, infer_schema_length=0,
                     new_columns=vcf_cols, dtypes=vcf_schema
                    )


    # Check how many region columns
    subset_cols = [vcf_cols[0], *vcf_cols[2:]] # skip pos0
    schema = df.collect_schema()
    intersect_ncols = len(schema.names())


    # Intersected with peak, check if region col needs to be made
    if intersect_ncols > vcf_ncols:

        subset_cols.append(region_col)

        # Contains a fourth column to be used as regions
        if use_region_names and (intersect_ncols - vcf_ncols) > 3:

            df = df.rename({df.columns[vcf_ncols+3]: region_col})

        else:
            df = df.with_columns(
                pl.concat_str(
                    [
                        pl.col(i) for i in schema.names()[vcf_ncols:vcf_ncols+3]
                    ],
                    separator="_"
                ).alias(region_col)
            )

        # Retrieve region col
        df = df.select(subset_cols)
        
    return df.unique(maintain_order=True).collect()


def parse_intersect_region(intersect_file, use_region_names=False, region_col=None):
    
    df = pl.scan_csv(intersect_file, separator="\t",
                     has_header=False, infer_schema_length=0)
    
    # If we need to use coords as name
    use_coords = False
    
    if region_col is None:
        region_col = "region"

    # No regions, only variants
    if len(df.columns) <= 5:
        # No regions, only variants
        subset_cols = [df.columns[i] for i in [0, 2, 3, 4]]
        new_cols = ["chrom", "pos", "ref", "alt"]
    
    elif use_region_names and len(df.columns) >= 9:
        # Use included names in region file
        subset_cols = [df.columns[i] for i in [0, 2, 3, 4, 8]]
        new_cols = ["chrom", "pos", "ref", "alt", region_col]
        
    elif len(df.columns) >= 8:
        # Either no names included or use coords instead
        subset_cols = [df.columns[i] for i in [0, 2, 3, 4, 5, 6, 7]]
        new_cols = ["chrom", "pos", "ref", "alt",
                    "region_chrom", "region_start", "region_end"]
        use_coords = True

    else:
        raise ValueError(f"Could not recognize BED format. Expected 3-6 columns, got {len(df.columns)} columns")

    # Parse dataframe columns
    rename_cols = {old_col: new_col for old_col, new_col in zip(subset_cols, new_cols)}
    df = df.select(subset_cols).rename(
        rename_cols).with_columns(
            [
                pl.col("chrom").cast(pl.Categorical),
                pl.col("pos").cast(pl.UInt32),
                pl.col("ref").cast(pl.Categorical),
                pl.col("alt").cast(pl.Categorical)
                ]
            )
    
    # Create coords
    if use_coords:
        df = df.with_columns(
            pl.concat_str(
                [pl.col(i) for i in new_cols[-3::]],
                separator="_"
            ).alias("region")
        ).select(
            ["chrom", "pos", "ref", "alt", "region"])
    
    return df.unique(maintain_order=True).collect()

