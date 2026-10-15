"""Gene annotation parsing and data management."""

from __future__ import annotations

import re
from pathlib import Path
from typing import NamedTuple

import polars as pl


class ParsedGeneData(NamedTuple):
    """Parsed gene data from GTF/GFF3 file."""

    gene_df: pl.DataFrame
    feature: str
    attribute: str
    parent_attribute: str


class WaspGeneData:
    """Container for gene annotation file paths and configuration.

    Attributes
    ----------
    gene_file : str | Path
        Path to the gene annotation file.
    feature : str | None
        Feature type to extract.
    attribute : str | None
        Attribute for region names.
    parent_attribute : str | None
        Parent attribute for hierarchical features.
    """

    def __init__(
        self,
        gene_file: str | Path,
        feature: str | None = None,
        attribute: str | None = None,
        parent_attribute: str | None = None,
    ) -> None:
        self.gene_file = gene_file
        self.feature = feature
        self.attribute = attribute
        self.parent_attribute = parent_attribute

    def update_data(self, data: ParsedGeneData) -> None:
        """Update attributes with parsed data.

        Parameters
        ----------
        data : ParsedGeneData
            Parsed gene data to update from.
        """
        for key in data._fields:
            if hasattr(self, key):
                setattr(self, key, getattr(data, key))


def _attribute_values(attribute: str) -> pl.Expr:
    """Extract an exact GTF/GFF3 attribute key, including the final attribute."""
    pattern = rf"(?:^|;)\s*{re.escape(attribute)}(?:=|\s+)\s*([^;]*)"
    return pl.col("attribute").str.extract(pattern).str.strip_chars(" \t\"'")


def parse_gene_file(
    gene_file: str | Path,
    feature: str | None = None,
    attribute: str | None = None,
    parent_attribute: str | None = None,
) -> ParsedGeneData:
    """Parse GTF/GFF3 gene annotation file.

    Parameters
    ----------
    gene_file : str | Path
        Path to GTF/GFF3 file.
    feature : str | None, optional
        Feature type to extract (exon, then transcript, then gene if None).
    attribute : str | None, optional
        Attribute for region names (auto-detected if None).
    parent_attribute : str | None, optional
        Parent attribute for hierarchical features (auto-detected if None).

    Returns
    -------
    ParsedGeneData
        Named tuple with (gene_df, feature, attribute, parent_attribute).
    """
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

    if feature is None:
        features = df.get_column("feature")
        feature = next((name for name in ("exon", "transcript", "gene") if name in features), None)
    if feature is None:
        raise ValueError("Could not determine feature type: expected exon, transcript or gene")

    df = df.filter(pl.col("feature") == feature)
    if df.is_empty():
        raise ValueError(f"No {feature!r} features found in {gene_file}")

    def has_attribute(name: str) -> bool:
        return bool(df.select(_attribute_values(name).is_not_null().any()).item())

    if attribute is None:
        attribute = next(
            (name for name in (f"{feature}_id", "ID", "Name") if has_attribute(name)), None
        )
    if attribute is None:
        raise ValueError(
            f"No feature identifier found for {feature!r} in {gene_file}; "
            f"expected {feature}_id, ID or Name. Specify --gene_attribute."
        )

    if parent_attribute is None:
        parent_attribute = next(
            (name for name in ("Parent", "transcript_id", "gene_id") if has_attribute(name)),
            attribute,
        )
    parent_col = f"groupby_{attribute}" if parent_attribute == attribute else parent_attribute

    # Validate before writing a BED so missing IDs cannot reach interval intersection.
    for name in dict.fromkeys((attribute, parent_attribute)):
        invalid = df.filter(_attribute_values(name).is_null() | (_attribute_values(name) == ""))
        if not invalid.is_empty():
            examples = ", ".join(
                f"{row['seqname']}:{row['start']}-{row['end']}"
                for row in invalid.head(3).iter_rows(named=True)
            )
            raise ValueError(
                f"Missing or empty annotation attribute {name!r} on {invalid.height} "
                f"{feature!r} record(s) in {gene_file}: {examples}. "
                "Provide complete identifiers or select --gene_attribute/--gene_parent."
            )

    df = df.with_columns(
        pl.col("start").sub(1),
        _attribute_values(attribute).alias(attribute),
        _attribute_values(parent_attribute).alias(parent_col),
    ).select(["seqname", "start", "end", attribute, parent_col])

    return ParsedGeneData(df, feature, attribute, parent_col)


def make_gene_data(
    gene_file: str | Path,
    out_bed: str | Path,
    feature: str | None = None,
    attribute: str | None = None,
    parent_attribute: str | None = None,
) -> WaspGeneData:
    """Parse gene file and create BED for intersection.

    Parameters
    ----------
    gene_file : str | Path
        Path to GTF/GFF3 file.
    out_bed : str | Path
        Output BED file path.
    feature : str | None, optional
        Feature type to extract.
    attribute : str | None, optional
        Attribute for region names.
    parent_attribute : str | None, optional
        Parent attribute for hierarchical features.

    Returns
    -------
    WaspGeneData
        Container with parsed gene data and configuration.
    """
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


def parse_intersect_genes(
    intersect_file: str | Path,
    attribute: str | None = None,
    parent_attribute: str | None = None,
) -> pl.DataFrame:
    """Parse gene intersection file (legacy version).

    Parameters
    ----------
    intersect_file : str | Path
        Path to bedtools intersection output.
    attribute : str | None, optional
        Attribute column name, by default 'ID'.
    parent_attribute : str | None, optional
        Parent attribute column name, by default 'Parent'.

    Returns
    -------
    pl.DataFrame
        Parsed intersection data.
    """
    if attribute is None:
        attribute = "ID"

    if parent_attribute is None:
        parent_attribute = "Parent"

    # Guard against empty intersection file (0 variants in region)
    intersect_path = Path(intersect_file)
    if not intersect_path.exists() or intersect_path.stat().st_size == 0:
        return pl.DataFrame(
            schema={
                "chrom": pl.Categorical,
                "pos": pl.UInt32,
                "ref": pl.Categorical,
                "alt": pl.Categorical,
                attribute: pl.Utf8,
                parent_attribute: pl.Utf8,
            }
        )

    df = pl.scan_csv(intersect_file, separator="\t", has_header=False, infer_schema_length=0)

    subset_cols = [df.columns[i] for i in [0, 2, 3, 4, -2, -1]]
    new_cols = ["chrom", "pos", "ref", "alt", attribute, parent_attribute]

    # Parse dataframe columns
    rename_cols = dict(zip(subset_cols, new_cols))
    df = df.select(subset_cols).rename(rename_cols).with_columns(pl.col("pos").cast(pl.UInt32))

    return df.unique(maintain_order=True).collect()


def parse_intersect_genes_new(
    intersect_file: str | Path,
    attribute: str | None = None,
    parent_attribute: str | None = None,
) -> pl.DataFrame:
    """Parse gene intersection file with typed columns.

    Parameters
    ----------
    intersect_file : str | Path
        Path to bedtools intersection output.
    attribute : str | None, optional
        Attribute column name, by default 'ID'.
    parent_attribute : str | None, optional
        Parent attribute column name, by default 'Parent'.

    Returns
    -------
    pl.DataFrame
        Parsed intersection data with typed columns.
    """
    if attribute is None:
        attribute = "ID"

    if parent_attribute is None:
        parent_attribute = "Parent"

    # Guard against empty intersection file (0 variants in region)
    intersect_path = Path(intersect_file)
    if not intersect_path.exists() or intersect_path.stat().st_size == 0:
        return pl.DataFrame(
            schema={
                "chrom": pl.Categorical,
                "pos": pl.UInt32,
                "ref": pl.Categorical,
                "alt": pl.Categorical,
                attribute: pl.Utf8,
                parent_attribute: pl.Utf8,
            }
        )

    # AFTER performing gtf_to_bed and intersecting!
    df = pl.scan_csv(intersect_file, separator="\t", has_header=False, infer_schema_length=0)
    schema_names = list(df.collect_schema().names())

    vcf_schema = [
        pl.col("chrom").cast(pl.Categorical),
        pl.col("pos").cast(pl.UInt32),
        pl.col("ref").cast(pl.Categorical),
        pl.col("alt").cast(pl.Categorical),
    ]

    # Expect at min 10 cols, 11 if GT included
    if len(schema_names) > 10:
        subset_cols = [schema_names[i] for i in [0, 2, 3, 4, 5, -2, -1]]
        new_cols = ["chrom", "pos", "ref", "alt", "GT", attribute, parent_attribute]
        vcf_schema.append(pl.col("GT").cast(pl.Categorical))
    else:
        subset_cols = [schema_names[i] for i in [0, 2, 3, 4, -2, -1]]
        new_cols = ["chrom", "pos", "ref", "alt", attribute, parent_attribute]

    # Parse dataframe columns
    rename_cols = dict(zip(subset_cols, new_cols))
    df = df.select(subset_cols).rename(rename_cols).with_columns(vcf_schema)

    return df.unique(maintain_order=True).collect()
