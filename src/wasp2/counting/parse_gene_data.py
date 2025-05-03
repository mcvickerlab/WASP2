import sys
import warnings
from pathlib import Path
from typing import Any, Optional

import polars as pl
from collections import namedtuple


class WaspGeneData:
    """
    Class to hold gene-related data for analysis.

    This class stores the path to the gene file and optional parameters for filtering,
    such as feature, attribute, and parent attribute. It may later be extended to create
    a DataFrame directly.
    """

    def __init__(
        self,
        gene_file: str,
        feature: Optional[str] = None,
        attribute: Optional[str] = None,
        parent_attribute: Optional[str] = None,
    ) -> None:
        """
        Initialize a WaspGeneData object.

        Parameters
        ----------
        gene_file : str
            Path to the gene file (e.g., GTF file).
        feature : str, optional
            Feature type to filter on (e.g., "exon", "transcript", "gene").
        attribute : str, optional
            Attribute to extract from the GTF attribute column.
        parent_attribute : str, optional
            Parent attribute to extract for grouping (e.g., "Parent", "gene_id").
        """
        self.gene_file: str = gene_file
        self.feature: Optional[str] = feature
        self.attribute: Optional[str] = attribute
        self.parent_attribute: Optional[str] = parent_attribute

        # Maybe should create dataframe in here...

    def update_data(self, data: Any) -> None:
        """
        Update attributes of the instance using values from a namedtuple.

        This method updates instance attributes with matching keys found in the provided
        namedtuple 'data'.

        Parameters
        ----------
        data : namedtuple
            A namedtuple with fields corresponding to attributes of the instance.

        Returns
        -------
        None
        """
        # Update attributes with namedtuple after parsing
        # Only updates matching keys
        for key in data._fields:
            if hasattr(self, key):
                setattr(self, key, getattr(data, key))


def parse_gene_file(
    gene_file: str,
    feature: Optional[str] = None,
    attribute: Optional[str] = None,
    parent_attribute: Optional[str] = None,
) -> Any:
    """
    Parse a GTF gene file and extract relevant information.

    This function reads a gene file (GTF) using Polars, applies a feature filter, and
    extracts specified attributes using regular expressions. It also attempts to determine
    default values for feature, attribute, and parent attribute if not provided.

    Parameters
    ----------
    gene_file : str
        Path to the gene file (GTF).
    feature : str, optional
        Feature to filter on (e.g., "exon", "transcript", "gene"). If not provided, a default
        is chosen based on the unique features present in the file.
    attribute : str, optional
        Attribute to extract from the GTF attribute column. Defaults to a feature-specific ID
        (e.g., "exon_id") or falls back to "ID" or "Name" if available.
    parent_attribute : str, optional
        Parent attribute to extract. Defaults to "Parent", "transcript_id", or "gene_id" based
        on what is found in the file. If none found, falls back to the attribute.

    Returns
    -------
    parsed_data : namedtuple
        A namedtuple with fields:
            - gene_df: The parsed Polars DataFrame containing gene data.
            - feature: The feature used for filtering.
            - attribute: The attribute extracted.
            - parent_attribute: The name of the parent attribute column (or its derived alias).
    """
    # Use gtf col names
    gtf_cols = [
        "seqname", "source", "feature",
        "start", "end", "score",
        "strand", "frame", "attribute"
    ]
    
    # Cant use lazyframe in case of compressed
    df = pl.read_csv(
        gene_file,
        separator="\t",
        comment_prefix="#",
        has_header=False,
        new_columns=gtf_cols,
    )
    
    # Check if we need to preparse feature
    if feature is None:
        feature_list = df.select(pl.col("feature").unique()).to_series()
        if "exon" in feature_list:
            feature = "exon"
        elif "transcript" in feature_list:
            feature = "transcript"
        elif "gene" in feature_list:
            feature = "gene"
        else:
            print(f"exon, gene or transcript not found in feature list: \n{feature_list}")
    
    # feature filter
    df = df.filter(pl.col("feature") == feature)
    
    # Parse attributes
    if attribute is None:
        if df.get_column("attribute").str.contains(f"{feature}_id").all() is True:
            attribute = f"{feature}_id"
        elif df.get_column("attribute").str.contains("ID").all() is True:
            attribute = "ID"
        elif df.get_column("attribute").str.contains("Name").all() is True:
            attribute = "Name"
        else:
            # TODO return an error
            # TODO maybe just use region or coords as a feature
            print(f"No 'ID', '{feature}_id' or 'Name' attribute found. Please include ID")
    
    # TODO: Figure out best way to handle parent attribute
    # Parse parent attributes
    # Figure out best way to parse and handle this
    if parent_attribute is None:
        # Defaults to gene (possibly transcript???)
        if df.get_column("attribute").str.contains("Parent").all() is True:
            parent_attribute = "Parent"
        elif df.get_column("attribute").str.contains("transcript_id").all() is True:
            parent_attribute = "transcript_id"
        elif df.get_column("attribute").str.contains("gene_id").all() is True:
            parent_attribute = "gene_id"
        else:
            parent_attribute = attribute
    
    # TODO: Allow for count output without parent column
    if parent_attribute == attribute:
        parent_col = f"groupby_{attribute}"
    else:
        parent_col = parent_attribute
        
    # Extract relevant attributes
    attr_regex = fr'{attribute}[=\s]\"?\'?(.*?)\"?\'?;'
    parent_regex = fr'{parent_attribute}[=\s]\"?\'?(.*?)\"?\'?;'
    
    df = df.with_columns(
        pl.col("start").sub(1),
        pl.col("attribute").str.extract(attr_regex).alias(attribute),
        pl.col("attribute").str.extract(parent_regex).alias(parent_col)
    ).select(["seqname", "start", "end", attribute, parent_col])
    
    # metadata...maybe should create a method
    parsed_data = namedtuple(
        "parsed_data",
        ["gene_df", "feature", "attribute", "parent_attribute"]
    )
    
    return parsed_data(df, feature, attribute, parent_col)


def make_gene_data(
    gene_file: str,
    out_bed: str,
    feature: Optional[str] = None,
    attribute: Optional[str] = None,
    parent_attribute: Optional[str] = None,
) -> WaspGeneData:
    """
    Parse a gene file and write a BED file for gene intersections.

    This function creates a WaspGeneData object from the provided gene file and parameters.
    It then parses the gene file using :func:`parse_gene_file`, updates the gene data object,
    writes the resulting DataFrame as a BED file, and returns the updated gene data object.

    Parameters
    ----------
    gene_file : str
        Path to the gene file (e.g., GTF).
    out_bed : str
        Path where the BED file will be written.
    feature : str, optional
        Feature type to filter on.
    attribute : str, optional
        Attribute to extract.
    parent_attribute : str, optional
        Parent attribute to extract.

    Returns
    -------
    WaspGeneData
        The updated WaspGeneData object containing parsed gene information.
    """
    gene_data = WaspGeneData(
        gene_file=gene_file,
        feature=feature,
        attribute=attribute,
        parent_attribute=parent_attribute,
    )
    
    parsed_file_data = parse_gene_file(
        gene_data.gene_file,
        feature=gene_data.feature,
        attribute=gene_data.attribute,
        parent_attribute=gene_data.parent_attribute,
    )
    
    # Update gene_data obj
    gene_data.update_data(parsed_file_data)
    
    # Write out_bed
    parsed_file_data.gene_df.write_csv(
        out_bed, separator="\t", include_header=False
    )
    
    return gene_data


def parse_intersect_genes(
    intersect_file: str,
    attribute: Optional[str] = None,
    parent_attribute: Optional[str] = None,
) -> pl.DataFrame:
    """
    Parse an intersection file for genes into a Polars DataFrame.

    This function reads an intersection file produced after performing gtf_to_bed and intersection,
    selects and renames relevant columns, and casts the "pos" column to UInt32. It returns a unique
    and ordered DataFrame.

    Parameters
    ----------
    intersect_file : str
        Path to the intersection file.
    attribute : str, optional
        Attribute column name to use. Defaults to "ID" if not provided.
    parent_attribute : str, optional
        Parent attribute column name to use. Defaults to "Parent" if not provided.

    Returns
    -------
    polars.DataFrame
        A unique, ordered DataFrame with columns: "chrom", "pos", "ref", "alt", attribute, parent_attribute.
    """
    if attribute is None:
        attribute = "ID"
    
    if parent_attribute is None:
        parent_attribute = "Parent"

    # AFTER performing gtf_to_bed and intersecting!
    df = pl.scan_csv(
        intersect_file,
        separator="\t",
        has_header=False,
        infer_schema_length=0,
    )
    
    # Might want to do checks for wrong number of columns
    subset_cols = [df.columns[i] for i in [0, 2, 3, 4, -2, -1]]
    new_cols = ["chrom", "pos", "ref", "alt", attribute, parent_attribute]
    
    # Parse dataframe columns
    rename_cols = {old_col: new_col for old_col, new_col in zip(subset_cols, new_cols)}
    df = df.select(subset_cols).rename(rename_cols).with_columns(
        pl.col("pos").cast(pl.UInt32)
    )
    
    return df.unique(maintain_order=True).collect()


def parse_intersect_genes_new(
    intersect_file: str,
    attribute: Optional[str] = None,
    parent_attribute: Optional[str] = None,
) -> pl.DataFrame:
    """
    Parse an intersection file for genes (new format) into a Polars DataFrame.

    This function reads an intersection file and expects at least 10 columns (or 11 if GT is included).
    It selects, renames, and casts columns appropriately based on the presence of genotype (GT) data.
    It returns a unique, ordered DataFrame with standardized column names.

    Parameters
    ----------
    intersect_file : str
        Path to the intersection file.
    attribute : str, optional
        Attribute column name to use. Defaults to "ID" if not provided.
    parent_attribute : str, optional
        Parent attribute column name to use. Defaults to "Parent" if not provided.

    Returns
    -------
    polars.DataFrame
        A unique, ordered DataFrame with columns including "chrom", "pos", "ref", "alt", and optionally "GT",
        followed by attribute and parent_attribute.
    """
    if attribute is None:
        attribute = "ID"
    
    if parent_attribute is None:
        parent_attribute = "Parent"

    # AFTER performing gtf_to_bed and intersecting!
    df = pl.scan_csv(
        intersect_file,
        separator="\t",
        has_header=False,
        infer_schema_length=0,
    )

    vcf_schema = [
        pl.col("chrom").cast(pl.Categorical),
        pl.col("pos").cast(pl.UInt32),
        pl.col("ref").cast(pl.Categorical),
        pl.col("alt").cast(pl.Categorical)
    ]

    # Expect at least 10 columns, 11 if GT included
    if len(df.columns) > 10:
        subset_cols = [df.columns[i] for i in [0, 2, 3, 4, 5, -2, -1]]
        new_cols = ["chrom", "pos", "ref", "alt", "GT", attribute, parent_attribute]
        vcf_schema.append(pl.col("GT").cast(pl.Categorical))
    else:
        subset_cols = [df.columns[i] for i in [0, 2, 3, 4, -2, -1]]
        new_cols = ["chrom", "pos", "ref", "alt", attribute, parent_attribute]
    
    # Parse dataframe columns
    rename_cols = {old_col: new_col for old_col, new_col in zip(subset_cols, new_cols)}
    df = df.select(subset_cols).rename(rename_cols).with_columns(vcf_schema)
    
    return df.unique(maintain_order=True).collect()
