import sys
import warnings

from pathlib import Path

import polars as pl

from collections import namedtuple

# Hold relevant gene data
class WaspGeneData:
    
    def __init__(self, gene_file,
                 feature=None,
                 attribute=None,
                 parent_attribute=None):
        
        self.gene_file = gene_file
        self.feature = feature
        self.attribute = attribute
        self.parent_attribute = parent_attribute
        
        # Maybe should create dataframe in here...
    

    def update_data(self, data):
        
        # Update attributes with namedtuple after parsing
        # Only updates matching keys
        for key in data._fields:
            if hasattr(self, key):
                setattr(self, key,
                        getattr(data, key)
                       )


def parse_gene_file(gene_file, feature=None, attribute=None, parent_attribute=None):
    
    # Use gtf col names
    gtf_cols = [
        "seqname", "source", "feature",
        "start", "end", "score",
        "strand", "frame", "attribute"]
    
    
    # Cant use lazyframe in case of compressed
    df = pl.read_csv(gene_file, separator="\t",
                     comment_prefix="#",
                     has_header=False,
                     new_columns=gtf_cols)
    
    # Check if we need to preparse feature
    if feature is None:
        
        feature_list = df.select(pl.col("feature").unique()).to_series()

        if "exon" in feature_list:
            feature = "exon"
        elif "gene" in feature_list:
            feature = "gene"
        else:
            # TODO return an error
            print("Exon and Gene not found in feature list")
    
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
        
        # Defaults to gene(possibly transcript???)

        if df.get_column("attribute").str.contains("Parent").all() is True:
            parent_attribute = "Parent"
        elif df.get_column("attribute").str.contains("gene_id").all() is True:
            parent_attribute = "gene_id"
        elif df.get_column("attribute").str.contains("transcript_id").all() is True:
            parent_attribute = "transcript_id"
        else:
            parent_attribute = attribute
    
    # Extract relevant attributes
    attr_regex = fr'{attribute}[=\s]\"?\'?(.*?)\"?\'?;' # OG attribute extract
    parent_regex = fr'{parent_attribute}[=\s]\"?\'?(.*?)\"?\'?;' # OG attribute extract
    
    df = df.with_columns(
        pl.col("start").sub(1),
        pl.col("attribute").str.extract(attr_regex).alias(attribute),
        pl.col("attribute").str.extract(parent_regex).alias(parent_attribute)
    ).select(["seqname", "start", "end", attribute, parent_attribute])
    
    # metadata...maybe should create a method
    parsed_data = namedtuple(
        "parsed_data",
        ["gene_df", "feature", "attribute", "parent_attribute"]
    )
    
    return parsed_data(df, feature, attribute, parent_attribute)


def parse_intersect_genes(intersect_file, attribute=None, parent_attribute=None):
    
    if attribute is None:
        attribute = "ID"
    
    if parent_attribute is None:
        parent_attribute = "Parent"


    df = pl.scan_csv(intersect_file, separator="\t",
                     has_header=False, infer_schema_length=0)
    
    # Should i poossibly consider diff number of cols?
    
    # Might want to do checks for wrong number of columns
    subset_cols = [df.columns[i] for i in [0, 2, 3, 4, -2, -1]]
    new_cols = ["chrom", "pos", "ref", "alt", attribute, parent_attribute]
    
    # Parse dataframe columns
    rename_cols = {old_col: new_col for old_col, new_col in zip(subset_cols, new_cols)}
    df = df.select(subset_cols).rename(
        rename_cols).with_columns(pl.col("pos").cast(pl.UInt32))
    
    return df.unique(maintain_order=True).collect()


