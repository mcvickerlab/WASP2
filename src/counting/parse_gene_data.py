from collections import namedtuple

import polars as pl


# Hold relevant gene data
class WaspGeneData:
    def __init__(self, gene_file, feature=None, attribute=None, parent_attribute=None):
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
                setattr(self, key, getattr(data, key))


def parse_gene_file(gene_file, feature=None, attribute=None, parent_attribute=None):
    # Use gtf col names
    gtf_cols = [
        "seqname",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attribute",
    ]

    # Cant use lazyframe in case of compressed
    df = pl.read_csv(
        gene_file, separator="\t", comment_prefix="#", has_header=False, new_columns=gtf_cols
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
        # Defaults to gene(possibly transcript???)

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
    attr_regex = rf"{attribute}[=\s]\"?\'?(.*?)\"?\'?;"
    parent_regex = rf"{parent_attribute}[=\s]\"?\'?(.*?)\"?\'?;"

    df = df.with_columns(
        pl.col("start").sub(1),
        pl.col("attribute").str.extract(attr_regex).alias(attribute),
        pl.col("attribute").str.extract(parent_regex).alias(parent_col),
    ).select(["seqname", "start", "end", attribute, parent_col])

    # metadata...maybe should create a method
    parsed_data = namedtuple("parsed_data", ["gene_df", "feature", "attribute", "parent_attribute"])

    return parsed_data(df, feature, attribute, parent_col)


# Parse and create gtf_bed for intersection
# and return parsed WaspGeneData obj
def make_gene_data(
    gene_file,
    out_bed,
    feature=None,
    attribute=None,
    parent_attribute=None,
):
    gene_data = WaspGeneData(
        gene_file=gene_file, feature=feature, attribute=attribute, parent_attribute=parent_attribute
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
    parsed_file_data.gene_df.write_csv(out_bed, separator="\t", include_header=False)

    return gene_data


def parse_intersect_genes(intersect_file, attribute=None, parent_attribute=None):
    if attribute is None:
        attribute = "ID"

    if parent_attribute is None:
        parent_attribute = "Parent"

    # AFTER performing gtf_to_bed and intersecting!
    df = pl.scan_csv(intersect_file, separator="\t", has_header=False, infer_schema_length=0)

    # Should i poossibly consider diff number of cols?

    # Might want to do checks for wrong number of columns
    subset_cols = [df.columns[i] for i in [0, 2, 3, 4, -2, -1]]
    new_cols = ["chrom", "pos", "ref", "alt", attribute, parent_attribute]

    # Parse dataframe columns
    rename_cols = dict(zip(subset_cols, new_cols))
    df = df.select(subset_cols).rename(rename_cols).with_columns(pl.col("pos").cast(pl.UInt32))

    return df.unique(maintain_order=True).collect()


def parse_intersect_genes_new(intersect_file, attribute=None, parent_attribute=None):
    if attribute is None:
        attribute = "ID"

    if parent_attribute is None:
        parent_attribute = "Parent"

    # AFTER performing gtf_to_bed and intersecting!
    df = pl.scan_csv(intersect_file, separator="\t", has_header=False, infer_schema_length=0)

    vcf_schema = [
        pl.col("chrom").cast(pl.Categorical),
        pl.col("pos").cast(pl.UInt32),
        pl.col("ref").cast(pl.Categorical),
        pl.col("alt").cast(pl.Categorical),
    ]

    # Expect at min 10 cols, 11 if GT included
    if len(df.columns) > 10:
        subset_cols = [df.columns[i] for i in [0, 2, 3, 4, 5, -2, -1]]
        new_cols = ["chrom", "pos", "ref", "alt", "GT", attribute, parent_attribute]
        vcf_schema.append(pl.col("GT").cast(pl.Categorical))
    else:
        subset_cols = [df.columns[i] for i in [0, 2, 3, 4, -2, -1]]
        new_cols = ["chrom", "pos", "ref", "alt", attribute, parent_attribute]

    # Parse dataframe columns
    rename_cols = dict(zip(subset_cols, new_cols))
    df = df.select(subset_cols).rename(rename_cols).with_columns(vcf_schema)

    return df.unique(maintain_order=True).collect()
