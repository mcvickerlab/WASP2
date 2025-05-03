"""
Author: Aaron Ho
Python Version: 3.8
"""

# Default Python package Imports
import sys
import warnings
from collections import namedtuple
from pathlib import Path
from typing import Any, List, Optional, Union

import numpy as np
import pandas as pd
import anndata as ad

# local imports
from wasp2.analysis.as_analysis_sc import get_imbalance_sc, adata_count_qc


class WaspAnalysisSC:
    """
    Class to store and manage single-cell allele imbalance analysis parameters and input data.

    This class encapsulates configuration parameters for running allele imbalance analysis on single-cell
    data using an AnnData object. It stores input file paths, threshold values, grouping information, and
    other analysis settings. It also provides a method to update its attributes using a namedtuple.

    Attributes
    ----------
    adata_file : str
        Path to the AnnData file containing count data.
    bc_map : str
        Path to the barcode mapping file.
    min_count : int
        Minimum allele count threshold. Default is 10.
    pseudocount : int
        Pseudocount added to avoid zero counts. Default is 1.
    sample : str
        Sample identifier for the analysis.
    groups : Union[str, List[str]]
        List of groups to analyze.
    model : str
        Dispersion model used for analysis. Default is "single".
    out_file : str
        Output file path for analysis results.
    phased : bool
        Flag indicating whether genotype data is phased.
    z_cutoff : Optional[Union[int, float]]
        Z-score cutoff used in prefiltering the AnnData object.
    out_dir : Path
        Directory derived from out_file.
    prefix : str
        Filename prefix derived from out_file.
    """

    def __init__(
        self,
        adata_file: str,
        bc_map: str,
        min_count: Optional[int] = None,
        pseudocount: Optional[int] = None,
        phased: Optional[bool] = None,
        sample: Optional[str] = None,
        groups: Optional[Union[str, List[str]]] = None,
        model: Optional[str] = None,
        out_file: Optional[str] = None,
        z_cutoff: Optional[Union[int, float]] = None,
    ) -> None:
        """
        Initialize a WaspAnalysisSC instance with the provided parameters.

        Parameters
        ----------
        adata_file : str
            Path to the AnnData file containing count data.
        bc_map : str
            Path to the barcode mapping file.
        min_count : int, optional
            Minimum allele count threshold. Defaults to 10 if not provided.
        pseudocount : int, optional
            Pseudocount value added to counts to avoid division by zero. Defaults to 1.
        phased : bool, optional
            Flag indicating whether genotype data is phased. Defaults to False.
        sample : str, optional
            Sample identifier.
        groups : str or list, optional
            Groups to analyze. Can be a comma-delimited string or a list.
        model : str, optional
            Dispersion model to use. Only "single" is currently supported.
        out_file : str, optional
            Output file path for analysis results. Defaults to "ai_results.tsv" in the current directory.
        z_cutoff : int or float, optional
            Z-score cutoff for prefiltering the AnnData object (e.g., 4 or 5).

        Raises
        ------
        ValueError
            If min_count or pseudocount are not non-negative integers.
        """
        # User input data
        self.adata_file: str = adata_file
        self.bc_map: str = bc_map
        self.min_count: Optional[int] = min_count
        self.pseudocount: Optional[int] = pseudocount
        self.sample: Optional[str] = sample
        self.groups: Optional[Union[str, List[str]]] = groups
        self.model: Optional[str] = model
        self.out_file: Optional[str] = out_file
        self.phased: Optional[bool] = phased
        self.z_cutoff: Optional[Union[int, float]] = z_cutoff  # Should i default to something like 4 or 5?

        # Default to single dispersion model
        # TODO ADD GROUP DISP and other model types
        if (self.model is None) or (self.model not in {"single"}):
            self.model = "single"

        # Default min count of 10
        if self.min_count is None:
            self.min_count = 10

        if self.pseudocount is None:
            # self.pseudocount = 0 # either 0 or 1 for default
            self.pseudocount = 1

        # Make sure min_count and pseudocount are valid
        if not all([(i >= 0) and isinstance(i, int) for i in (self.min_count, self.pseudocount)]):
            raise ValueError("min_count and pseudocount must be non-negative integers")

        # Handle group inputs as strings and convert to list if necessary
        if isinstance(self.groups, str):
            # Check if groups is a file or a comma-delimited string
            if Path(self.groups).is_file():
                with open(self.groups) as group_file:
                    self.groups = [line.strip() for line in group_file]
            else:
                self.groups = [s.strip() for s in self.groups.split(",")]

        # Create default outfile if not provided
        if self.out_file is None:
            self.out_file = str(Path.cwd() / "ai_results.tsv")

        # Process output names for groups
        self.out_dir: Path = Path(self.out_file).parent
        self.prefix: str = Path(self.out_file).stem

    def update_data(self, data: Any) -> None:
        """
        Update instance attributes using values from a namedtuple.

        This method updates the attributes of the instance with matching keys from the provided
        namedtuple 'data'.

        Parameters
        ----------
        data : namedtuple
            A namedtuple containing fields that match some of the attributes of the instance.

        Returns
        -------
        None
        """
        # Update attributes with namedtuple after parsing
        # Only updates matching keys
        for key in data._fields:
            if hasattr(self, key):
                setattr(self, key, getattr(data, key))


def process_adata_inputs(
    adata: ad.AnnData,
    ai_files: Optional[WaspAnalysisSC] = None,
    bc_map: Optional[str] = None,
    sample: Optional[str] = None,
    groups: Optional[Union[str, List[str]]] = None,
    phased: Optional[bool] = None,
) -> Any:
    """
    Process and validate AnnData inputs for single-cell analysis.

    This function processes an AnnData object along with optional parameters. It checks for sample
    and genotype information, subsets the AnnData object to include only heterozygous genotypes, and
    updates region indices if applicable. Barcode mapping is applied if provided. Finally, it returns a
    namedtuple containing the processed AnnData object, sample, groups, and phased flag.

    Parameters
    ----------
    adata : anndata.AnnData
        AnnData object containing count data.
    ai_files : WaspAnalysisSC, optional
        Instance of WaspAnalysisSC containing analysis parameters. If provided, its attributes are used.
    bc_map : str, optional
        Path to the barcode mapping file.
    sample : str, optional
        Sample identifier.
    groups : list or str, optional
        Groups to analyze.
    phased : bool, optional
        Flag indicating whether genotype data is phased.

    Returns
    -------
    namedtuple
        A namedtuple with fields ["adata", "sample", "groups", "phased"] representing the processed inputs.

    Raises
    ------
    KeyError
        If the provided sample is not found in the AnnData object.
    ValueError
        If multiple samples exist in the AnnData object and no sample is provided.
    """
    if ai_files is not None:
        bc_map = ai_files.bc_map
        sample = ai_files.sample
        groups = ai_files.groups
        phased = ai_files.phased

    # Check genotype and phasing input 
    if "samples" not in adata.uns_keys():
        if sample is not None:
            raise KeyError(f"Sample '{sample}' provided, but no samples found in count data")
        phased = False
    elif sample is None:
        sample_list = adata.uns["samples"]
        if len(sample_list) != 1:
            raise ValueError(
                "Genotype Ambiguous: Count data contains multiple samples, but none provided. "
                f"Please input a sample from list: {sample_list}"
            )
        else:
            sample = sample_list[0]
    else:
        sample_list = adata.uns["samples"]
        if sample not in sample_list:
            raise KeyError(
                f"Sample: '{sample}' not found in dataset. "
                f"Please input a sample from list: {sample_list}"
            )
        else:
            # Subset to include only heterozygous genotypes
            if not any(i in ['1|0', '0|1', '1/0', '0/1'] for i in adata.obs[sample].unique()):
                raise ValueError(f"Heterozygous genotypes not found for sample: {sample}.")
            adata = adata[adata.obs[sample].isin(['1|0', '0|1', '1/0', '0/1'])].copy()
            adata.obs = adata.obs.reset_index(drop=True)  # Reset index after subsetting

            # Reindex regions after filtering genotype calls
            if "feature" in adata.uns_keys():
                # idx_df = adata.obs[["index"]].reset_index(drop=True).copy().reset_index(names="filt_index")
                adata.uns["feature"] = adata.uns["feature"].merge(
                    adata.obs[["index"]].reset_index(names="filt_index"),
                    on="index")[["region", "filt_index"]].rename(columns={"filt_index": "index"})
                # Update adata.obs index column as well
                adata.obs["index"] = adata.obs.index

    # Check phasing if True or None
    if phased is not False:
        if {'0|1', '1|0'} == set(adata.obs[sample].unique()):
            phased = True
        else:
            phased = False
            warning_msg = (
                f"Phased model selected for unphased genotypes ({adata.obs[sample].unique()}). "
                "Switching to unphased model"
            )
            warnings.warn(warning_msg)

    # Add groups if barcode mapping provided
    if bc_map is not None:
        map_df = pd.read_csv(bc_map, sep="\t", header=None, names=["group"], index_col=0, dtype="category")
        adata.var = adata.var.join(map_df, how="left")

    # No existing groups or mapping provided
    if "group" not in adata.var_keys():
        raise KeyError("groups not found in dataset, please provide a barcode mapping")
    elif groups is not None:
        valid_groups = list(adata.var["group"].dropna().unique())
        new_groups = [i for i in groups if i in valid_groups]
        if len(new_groups) == 0:
            raise KeyError(f"Provided groups {groups} not found.")
        elif len(new_groups) < len(groups):
            skipped_groups = [i for i in groups if i not in new_groups]
            print(f"Skipping missing groups: {skipped_groups}")
            groups = new_groups
        else:
            groups = new_groups
    else:
        groups = list(adata.var["group"].dropna().unique())

    # Return processed inputs as a namedtuple
    adata_inputs = namedtuple("adata_inputs", ["adata", "sample", "groups", "phased"])
    return adata_inputs(adata, sample, groups, phased)


def run_ai_analysis_sc(
    count_file: str,
    bc_map: str,
    min_count: Optional[int] = None,
    pseudocount: Optional[int] = None,
    phase: Optional[bool] = None,
    sample: Optional[str] = None,
    groups: Optional[Union[str, List[str]]] = None,
    out_file: Optional[str] = None,
    z_cutoff: Optional[Union[int, float]] = None,
) -> None:
    """
    Run the single-cell allele imbalance analysis pipeline.

    This function creates a WaspAnalysisSC instance with the provided parameters, processes the AnnData
    input, applies quality control filtering, and computes allele imbalance statistics for specified groups.
    The resulting dataframes are written as TSV files to the output directory.

    Parameters
    ----------
    count_file : str
        Path to the AnnData file (e.g., an .h5ad file) containing count data.
    bc_map : str
        Path to the barcode mapping file.
    min_count : int, optional
        Minimum allele count threshold. Defaults to 10.
    pseudocount : int, optional
        Pseudocount added to counts. Defaults to 1.
    phase : bool, optional
        Flag indicating whether the data is phased.
    sample : str, optional
        Sample identifier.
    groups : str or list, optional
        Groups to analyze. Can be provided as a comma-delimited string or a list.
    out_file : str, optional
        Output file path for analysis results. Defaults to "ai_results.tsv" in the current directory.
    z_cutoff : int or float, optional
        Z-score cutoff for prefiltering the AnnData object.

    Returns
    -------
    None

    Examples
    --------
    >>> run_ai_analysis_sc("data.h5ad", "bc_map.tsv", min_count=10, pseudocount=1, phase=False, sample="Sample1", groups="Group1,Group2")
    """
    # Create data class that holds input data
    ai_files = WaspAnalysisSC(
        adata_file=count_file,
        bc_map=bc_map,
        min_count=min_count,
        pseudocount=pseudocount,
        phased=phase,
        sample=sample,
        groups=groups,
        model="single",
        out_file=out_file,
        z_cutoff=z_cutoff,
    )
    
    adata_inputs = process_adata_inputs(ad.read_h5ad(ai_files.adata_file), ai_files=ai_files)
    
    # Debugging prints (commented out)
    # print(*vars(ai_files).items(), sep="\n")
    # print(adata_inputs)
    
    # Update class attributes based on processed inputs
    ai_files.update_data(adata_inputs)
    
    # Prefilter AnnData count data using quality control
    adata = adata_count_qc(adata_inputs.adata,
                           z_cutoff=ai_files.z_cutoff,
                           gt_error=None)
    
    # Compute allele imbalance statistics for each group
    df_dict = get_imbalance_sc(
        adata,
        min_count=ai_files.min_count,
        pseudocount=ai_files.pseudocount,
        phased=ai_files.phased,
        sample=ai_files.sample,
        groups=ai_files.groups,
    )
    
    # Create output directory and write results for each group
    out_path = Path(ai_files.out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    for key, value in df_dict.items():
        group_out_file = out_path / f"{ai_files.prefix}_{str(key).replace('/', '-')}.tsv"
        value.sort_values(by="pval", ascending=True).to_csv(
            group_out_file, sep="\t", header=True, index=False)
        
    print(
        f"Allelic Imbalance measured for {len(df_dict)} groups!\n"
        f"Results written to: {out_path}/{ai_files.prefix}_[GROUP].tsv"
    )
