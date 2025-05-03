"""
Author: Aaron Ho
Python Version: 3.9
"""

# Default Python package Imports
from pathlib import Path
from csv import DictReader, reader
from typing import Optional

# External package imports
import pandas as pd

# Local script imports
from wasp2.analysis.as_analysis import get_imbalance


# TODO GOTTA IMPLEMENT THIS

class WaspAnalysisData:
    """
    Data structure for configuring and running allele imbalance analysis.

    This class stores analysis parameters and input file details. It is used to configure the allele
    imbalance analysis on count data. Parameters include thresholds, pseudocounts, phasing information,
    and the dispersion model to use. The class also attempts to automatically infer the region and grouping
    columns from the count file header.
    
    Attributes
    ----------
    count_file : str
        Path to the count file.
    min_count : int
        Minimum count threshold (default: 10).
    pseudocount : int
        Pseudocount added to avoid division by zero (default: 1).
    phased : bool
        Flag indicating if the data is phased (default: False).
    model : str
        Dispersion model to use ("single" or "linear"; default: "single").
    out_file : str
        Path to the output file for analysis results.
    region_col : str
        Column name that contains region information.
    groupby : Optional[str]
        Column name used for grouping the data. If None, grouping is performed by feature.
    """

    def __init__(
        self,
        count_file: str,
        min_count: Optional[int] = None,
        pseudocount: Optional[int] = None,
        phased: Optional[bool] = None,
        model: Optional[str] = None,
        out_file: Optional[str] = None,
        region_col: Optional[str] = None,
        groupby: Optional[str] = None,
    ) -> None:
        """
        Initialize a WaspAnalysisData instance with provided analysis parameters.

        Parameters
        ----------
        count_file : str
            Path to the count file.
        min_count : int, optional
            Minimum count threshold. Defaults to 10 if not provided.
        pseudocount : int, optional
            Pseudocount to add to avoid zero counts. Defaults to 1 if not provided.
        phased : bool, optional
            Flag indicating whether the data is phased. Defaults to False.
        model : str, optional
            Dispersion model to use ("single" or "linear"). Defaults to "single" if not provided.
        out_file : str, optional
            Path to the output file. Defaults to "ai_results.tsv" in the current working directory.
        region_col : str, optional
            Column name containing region information. If not provided, it is inferred from the header.
        groupby : str, optional
            Column name for grouping data. If provided, it overrides the default grouping by feature.
        
        Notes
        -----
        - The count file header is read to determine the number of columns and to infer the region and
          groupby columns.
        - If "GT" is present in the header, the minimum number of columns is set to 8; otherwise 7.
        - The implementation assumes a specific column ordering in the count file.
        """
        # User input data
        self.count_file: str = count_file
        self.min_count: Optional[int] = min_count
        self.pseudocount: Optional[int] = pseudocount
        self.phased: Optional[bool] = phased
        self.model: Optional[str] = model
        self.out_file: Optional[str] = out_file

        # Group by feature by default
        self.region_col: Optional[str] = region_col
        self.groupby: Optional[str] = groupby  # group by region or parent?

        # TODO parse VCF for phased instead of default unphased
        if not self.phased:
            self.phased = False

        # Default to single dispersion model
        if (self.model is None) or (self.model not in {"single", "linear"}):
            self.model = "single"

        # Default min count of 10 
        if self.min_count is None:
            self.min_count = 10

        if self.pseudocount is None:
            # self.pseudocount = 0 # either 0 or 1 for default
            self.pseudocount = 1

        # Read header only for validation
        with open(self.count_file) as f:
            count_cols = next(reader(f, delimiter="\t"))

        # 7 columns at minimum, 10 at maximum
        # 3 required: chr, pos, ref, alt
        # 3 optional: <GT>, <region>, <parent>
        # 3 required: ref_count, alt_count, other_count 
        # [chr, pos, ref, alt, <GT>, <region>, <parent>, ref_c, alt_c, other_c]
        if "GT" in count_cols:
            min_cols = 8
            region_idx = 5
        else:
            min_cols = 7
            region_idx = 4

        # Check regions
        if self.region_col is None:
            if len(count_cols) > min_cols:
                self.region_col = count_cols[region_idx]

        # By default group by feature rather than parent?
        if self.groupby is not None:
            # If denoting to group by feature
            if (self.region_col is None) or (self.groupby == self.region_col):
                self.groupby = None
            elif ((len(count_cols) > (min_cols + 1)) and
                  self.groupby in {count_cols[region_idx + 1], "Parent", "parent"}):
                self.groupby = count_cols[region_idx + 1]  # Set group
            else:
                # Maybe throw error instead
                print(f"{self.groupby} not found in columns \n{count_cols}")
                self.groupby = None

        # Create default outfile 
        if self.out_file is None:
            self.out_file = str(Path.cwd() / "ai_results.tsv")


# The following commented-out code represents an alternative initialization approach.
# class WaspAnalysisData:
#
#     def __init__(self, count_file,
#                  min_count=None,
#                  model=None,
#                  phased=None,
#                  out_dir=None,
#                  out_file=None,
#                  region_col=None,
#                  features=None):
#         
#         # User input data
#         self.count_file = count_file
#         self.min_count = min_count
#         self.model = model
#         self.phased = phased # TODO
#         self.out_file = out_file
#         self.out_dir = out_dir  # should i replace this with out file???
#         self.region_col = region_col
#         self.features = features # TODO and also add rna-seq support back
#         
#         # I need to also add other things for single cell back
#         
#
#         # Default to single dispersion model
#         if self.model is None:
#             self.model = "single"
#         
#         # Default min count of 10 
#         if self.min_count is None:
#             self.min_count = 10
#         
#         
#         # Automatically parse region col
#         # Should i do this after the df is created?
#         if self.region_col is None:
#             
#             # Read header only
#             with open(self.count_file) as f:
#                 count_cols = next(reader(f, delimiter = "\t"))
#             
#             # Check region_col from file
#             if "region" in count_cols:
#                 self.region_col = "region" # default atac naming
#             elif "peak" in count_cols:
#                 self.region_col = "peak" # from previous implementation
#             elif "genes" in count_cols:
#                 self.region_col = "genes"
#             else:
#                 # SNPs only
#                 # df["region"] = df["chrom"] + "_" + df["pos"].astype(str)
#                 self.region_col = "region" # should i name as snp?
#
#         # Create default outfile 
#         if self.out_file is None:
#             self.out_file = str(Path.cwd() / "ai_results.tsv") # do this after


def run_ai_analysis(
    count_file: str,
    min_count: Optional[int] = None,
    pseudocount: Optional[int] = None,
    phased: Optional[bool] = None,
    model: Optional[str] = None,
    out_file: Optional[str] = None,
    region_col: Optional[str] = None,
    groupby: Optional[str] = None,
) -> None:
    """
    Run the allele imbalance analysis pipeline.

    This function creates a WaspAnalysisData instance with the provided parameters, then runs the allele
    imbalance analysis using the `get_imbalance` function. The resulting DataFrame is sorted by FDR-corrected
    p-value and written to a TSV file.

    Parameters
    ----------
    count_file : str
        Path to the count file.
    min_count : int, optional
        Minimum count threshold. Defaults to 10.
    pseudocount : int, optional
        Pseudocount added to avoid zero counts. Defaults to 1.
    phased : bool, optional
        Flag indicating whether the data is phased. Defaults to False.
    model : str, optional
        Dispersion model to use ("single" or "linear"). Defaults to "single" if not provided.
    out_file : str, optional
        Path to the output file for analysis results. Defaults to "ai_results.tsv" in the current working directory.
    region_col : str, optional
        Column name that contains region information.
    groupby : str, optional
        Column name for grouping the data, if applicable.

    Returns
    -------
    None

    Examples
    --------
    >>> run_ai_analysis("counts.tsv", min_count=10, pseudocount=1, phased=False, model="single")
    """
    # Store analysis data and parameters
    ai_files = WaspAnalysisData(
        count_file,
        min_count=min_count,
        pseudocount=pseudocount,
        phased=phased,
        model=model,
        out_file=out_file,
        region_col=region_col,
        groupby=groupby,
    )
    
    # Run analysis pipeline
    ai_df = get_imbalance(
        ai_files.count_file,
        min_count=ai_files.min_count,
        pseudocount=ai_files.pseudocount,
        method=ai_files.model,
        phased=ai_files.phased,
        region_col=ai_files.region_col,
        groupby=ai_files.groupby,
    )
    
    # Optionally sort results by FDR p-value
    ai_df = ai_df.sort_values(by="fdr_pval", ascending=True)
    
    # Write results to file
    ai_df.to_csv(ai_files.out_file, sep="\t", header=True, index=False)
