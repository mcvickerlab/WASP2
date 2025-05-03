import sys
import timeit
import subprocess
import warnings

from pathlib import Path

import numpy as np
import polars as pl
from typing import Optional, List  # or any other needed types

# same as in mapping...should create unified utils
def vcf_to_bed(vcf_file: str, out_bed: str, samples: list[str] = None, include_gt: bool = True) -> str:
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
    samples : list[str], optional
        List of sample names to filter on. If None, genotypes are dropped.
    include_gt : bool, optional
        If True, include genotype information in the output BED file. Defaults to True.

    Returns
    -------
    str
        The path to the output BED file.
    
    Notes
    -----
    The following commented-out code is preserved:
    
        # Maybe change this later?
        # out_bed = f"{out_dir}/filt_variants.bed"
    """
    # Maybe change this later?
    # out_bed = f"{out_dir}/filt_variants.bed"
    
    # Base commands
    view_cmd = ["bcftools", "view", str(vcf_file),
                "-m2", "-M2", "-v", "snps", "-Ou"]
    
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
        
        # If we include GT
        if include_gt:
            # Changed %TGT to GT, ref/alt -> 0/1
            query_cmd.append("%CHROM\t%POS0\t%END\t%REF\t%ALT[\t%GT]\n")
        else:
            query_cmd.append("%CHROM\t%POS0\t%END\t%REF\t%ALT\n")
    
    # Run Subprocess
    query_process = subprocess.run(query_cmd, input=view_process.stdout, check=True)
    
    return out_bed


def gtf_to_bed(gtf_file: str, out_bed: str, feature: str, attribute: str) -> str:
    """
    Convert a GTF file to BED format by filtering on a specific feature and attribute.

    This function reads a GTF file into a pandas DataFrame, filters the rows based on the
    provided feature(s) (if any), and writes the resulting DataFrame to a new GTF file in the
    specified output directory.
    
    Parameters
    ----------
    gtf_file : str
        Path to the input GTF file.
    out_bed : str
        Path where the output BED file will be written.
    feature : str
        Feature type to filter on (e.g., "exon", "CDS").
    attribute : str
        The attribute key to extract from the GTF attribute column.
    
    Returns
    -------
    str
        The path to the output BED file.
    
    Notes
    -----
    The following commented-out code is preserved:
    
        # Use gtf col names
    """
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
    attr_regex = fr'{attribute}[=\s]\"?\'?(.*?)\"?\'?;'  # works for gtf/gff3
    
    # Extract feature only and attributes
    df = df.filter(pl.col("feature") == feature).with_columns(
        pl.col("start").sub(1),
        pl.col("attribute").str.extract(attr_regex).alias(attribute)
    ).select(["seqname", "start", "end", attribute])
    
    # TODO Extra validation and may want to return some data?
    
    # Write to BED
    df.write_csv(out_bed, separator="\t", include_header=False)
    
    return out_bed


def intersect_vcf_region(vcf_file: str, region_file: str, out_file: str) -> None:
    """
    Perform an intersection between a VCF file and a region file using bedtools.

    This function runs the bedtools intersect command to identify variants in the VCF file
    that overlap regions defined in the region file, and writes the output to a specified file.
    
    Parameters
    ----------
    vcf_file : str
        Path to the VCF file.
    region_file : str
        Path to the region file (e.g., BED, Peaks, GTF).
    out_file : str
        Path where the intersection output will be written.
    
    Returns
    -------
    None
    """
    # Parse region file before or after???
    intersect_cmd = [
        "bedtools", "intersect", "-a", str(vcf_file),
        "-b", str(region_file), "-wb"
    ]
    
    # write intersect out
    with open(out_file, "w") as file:
        intersect_process = subprocess.run(intersect_cmd, stdout=file, check=True)


def parse_intersect_region_new(
    intersect_file: str,
    samples: Optional[List[str]] = None,
    use_region_names: bool = False,
    region_col: Optional[str] = None,
) -> pl.DataFrame:
    """
    Parse an intersection file into a Polars DataFrame with region information.

    This function reads an intersection file (e.g., produced by intersect_vcf_region) using Polars,
    defines column names and types based on the number of samples, and creates a new region column
    if necessary.
    
    Parameters
    ----------
    intersect_file : str
        Path to the intersection file.
    samples : list[str], optional
        List of sample names to include.
    use_region_names : bool, optional
        If True and sufficient extra columns are present, use the provided region names.
        Otherwise, construct region names from coordinate columns.
    region_col : str, optional
        Name to assign to the region column. Defaults to "region" if not provided.
    
    Returns
    -------
    polars.DataFrame
        A unique, ordered Polars DataFrame containing the intersection data.
    """
    if region_col is None:
        region_col = "region"

    # Default number of columns
    vcf_cols: List[str] = ["chrom", "pos0", "pos", "ref", "alt"]  # Default columns for vcf
    vcf_schema = [pl.Categorical, pl.UInt32, pl.UInt32, pl.Categorical, pl.Categorical]

    if samples is not None:
        vcf_cols.extend(samples)
        vcf_schema.extend([pl.Categorical] * len(samples))

    vcf_ncols = len(vcf_cols)

    # Process with gt
    df = pl.scan_csv(
        intersect_file,
        separator="\t",
        has_header=False,
        infer_schema_length=0,
        new_columns=vcf_cols,
        dtypes=vcf_schema,
    )

    # Check how many region columns
    subset_cols = [vcf_cols[0], *vcf_cols[2:]]  # skip pos0
    intersect_ncols = len(df.columns)

    # Intersected with peak, check if region col needs to be made
    if intersect_ncols > vcf_ncols:
        subset_cols.append(region_col)

        # Contains a fourth column to be used as regions
        if use_region_names and (intersect_ncols - vcf_ncols) > 3:
            df = df.rename({df.columns[vcf_ncols+3]: region_col})
        else:
            df = df.with_columns(
                pl.concat_str([pl.col(i) for i in df.columns[vcf_ncols:vcf_ncols+3]], separator="_").alias(region_col)
            )

        # Retrieve region col
        df = df.select(subset_cols)
        
    return df.unique(maintain_order=True).collect()


def parse_intersect_region(
    intersect_file: str,
    use_region_names: bool = False,
    region_col: Optional[str] = None,
) -> pl.DataFrame:
    """
    Parse an intersection file into a Polars DataFrame, determining region names based on column count.

    This function reads an intersection file with Polars and renames or constructs a region column
    depending on the number of columns present and whether region names are provided. If coordinates
    are used, a region name is constructed by concatenating coordinate columns.
    
    Parameters
    ----------
    intersect_file : str
        Path to the intersection file.
    use_region_names : bool, optional
        If True and there are enough columns, use the provided region names.
        Otherwise, use coordinates to construct a region name.
    region_col : str, optional
        The name to use for the region column. Defaults to "region" if not provided.
    
    Returns
    -------
    polars.DataFrame
        A unique, ordered Polars DataFrame with columns "chrom", "pos", "ref", "alt", and "region".
    """
    df = pl.scan_csv(intersect_file, separator="\t", has_header=False, infer_schema_length=0)
    
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
        # CHANGE TO RAISE ERROR
        print("COULD NOT RECOGNIZE FORMAT OR WRONG NUMBER OF COLS")
        return
    
    # Parse dataframe columns
    rename_cols = {old_col: new_col for old_col, new_col in zip(subset_cols, new_cols)}
    df = df.select(subset_cols).rename(rename_cols).with_columns(
        [
            pl.col("chrom").cast(pl.Categorical),
            pl.col("pos").cast(pl.UInt32),
            pl.col("ref").cast(pl.Categorical),
            pl.col("alt").cast(pl.Categorical),
        ]
    )
    
    # Create coords
    if use_coords:
        df = df.with_columns(
            pl.concat_str([pl.col(i) for i in new_cols[-3::]], separator="_").alias("region")
        ).select(["chrom", "pos", "ref", "alt", "region"])
    
    return df.unique(maintain_order=True).collect()
