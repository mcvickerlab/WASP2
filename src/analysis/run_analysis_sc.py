# Default Python package Imports
import sys
import warnings

from pathlib import Path
from typing import Optional, List, Dict, Union, Any, NamedTuple

import numpy as np
import pandas as pd
import anndata as ad
from anndata import AnnData

# local imports
from .as_analysis_sc import get_imbalance_sc, adata_count_qc

# Class that stores relevant data
class WaspAnalysisSC:

    def __init__(
        self,
        adata_file: Union[str, Path],
        bc_map: Union[str, Path],
        min_count: Optional[int] = None,
        pseudocount: Optional[int] = None,
        phased: Optional[bool] = None,
        sample: Optional[str] = None,
        groups: Optional[Union[str, List[str]]] = None,
        model: Optional[str] = None,
        out_file: Optional[Union[str, Path]] = None,
        z_cutoff: Optional[float] = None
    ) -> None:
        
        # User input data
        self.adata_file = adata_file
        self.bc_map = bc_map
        self.min_count = min_count
        self.pseudocount = pseudocount
        self.sample = sample
        self.groups = groups
        self.model = model
        self.out_file = out_file
        self.phased = phased
        self.z_cutoff = z_cutoff # Should i default to something like 4 or 5?

        # Default to single dispersion model
        # TODO ADD GROUP DISP and other model types
        if ((self.model is None) or
            (self.model not in {"single"})):
            
            self.model = "single"
        
        # Default min count of 10
        if self.min_count is None:
            self.min_count = 10


        if self.pseudocount is None:
            # self.pseudocount = 0 # either 0 or 1 for default
            self.pseudocount = 1
        
        
        # Make sure min and pseudocounts are valid
        if not all([(i >= 0) and isinstance(i, int) 
                    for i in (self.min_count, self.pseudocount)]):
            raise ValueError("min_count and pseudocount must be non-negative integers")
        

        # Handle group inputs as strings to list
        if isinstance(self.groups, str):
            
            # Check if group file or comma delim string
            if Path(self.groups).is_file():
                
                with open(self.groups) as group_file:
                    self.groups = [line.strip() for line in group_file]
            
            else:
                self.groups = [s.strip() for s in self.groups.split(",")]


        # Create default outfile 
        if self.out_file is None:
            self.out_file = str(Path.cwd() / "ai_results.tsv") # do this after

        
        # Process output names for groups
        self.out_dir = Path(self.out_file).parent
        self.prefix = Path(self.out_file).stem


    def update_data(self, data: NamedTuple) -> None:

        # Update attributes with namedtuple after parsing
        # Only updates matching keys
        for key in data._fields:
            if hasattr(self, key):
                setattr(self, key,
                        getattr(data, key)
                       )


# Define namedtuple for adata inputs
class AdataInputs(NamedTuple):
    adata: AnnData
    sample: str
    groups: List[str]
    phased: bool


# Process adata inputs
def process_adata_inputs(
    adata: AnnData,
    ai_files: Optional[WaspAnalysisSC] = None,
    bc_map: Optional[Union[str, Path]] = None,
    sample: Optional[str] = None,
    groups: Optional[List[str]] = None,
    phased: Optional[bool] = None
) -> AdataInputs:

    if ai_files is not None:
        bc_map = ai_files.bc_map
        sample = ai_files.sample
        # ai_files.groups is already converted to List[str] in __init__ if it was a string
        groups = ai_files.groups if isinstance(ai_files.groups, list) else None
        phased = ai_files.phased
    
    # Check genotype and phasing input 
    if "samples" not in adata.uns_keys():
        
        if sample is not None:
            raise KeyError(
                f"Sample '{sample}' provided, but no samples found in count data")
        
        phased = False

    elif sample is None:
        sample_list = adata.uns["samples"]
        
        if len(sample_list) != 1:
            raise ValueError(
                "Genotype Ambiguous: Count data contains mutiple samples, but none provided. "
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
            # We gotta subset to include het genotypes now
            if not any(i in ['1|0', '0|1', '1/0', '0/1'] for i in adata.obs[sample].unique()):
                raise ValueError(f"Heterozygous genotypes not found for sample: {sample}.")

            # adata = adata[adata.obs[sample].isin(['1|0', '0|1', '1/0', '0/1'])]
            
            # Using copy instead of view stops implicit mod warning, need to check memory usage
            adata = adata[adata.obs[sample].isin(['1|0', '0|1', '1/0', '0/1'])].copy()
            adata.obs = adata.obs.reset_index(drop=True) # Have to reset index every time i subset adata

            # Have to reindex the regions after filtering GT's
            if "feature" in adata.uns_keys():
                
                # idx_df = adata.obs[["index"]].reset_index(
                #     drop=True).copy().reset_index(names="filt_index")
                
                adata.uns["feature"] = adata.uns["feature"].merge(
                    adata.obs[["index"]].reset_index(names="filt_index"),
                    on="index")[["region", "filt_index"]].rename(
                    columns={"filt_index": "index"})
                    
                # Need to update adata.obs index col as well
                adata.obs["index"] = adata.obs.index

    # Check phasing if True or None
    if phased is not False:
        
        if {'0|1', '1|0'} == set(adata.obs[sample].unique()):
            phased = True
        else:
            phased = False
            warning_msg = (
                f"Phased model selected for unphased genotypes ({adata.obs[sample].unique()}). "
                "Switching to unphased model")
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

    # Ensure all required values are set (type narrowing for mypy)
    assert sample is not None, "sample must be set by this point"
    assert groups is not None, "groups must be set by this point"
    assert phased is not None, "phased must be set by this point"

    # Return properly typed namedtuple
    return AdataInputs(adata, sample, groups, phased)


# Parse user inputs and run entire pipeline
def run_ai_analysis_sc(
    count_file: Union[str, Path],
    bc_map: Union[str, Path],
    min_count: Optional[int] = None,
    pseudocount: Optional[int] = None,
    phase: Optional[bool] = None,
    sample: Optional[str] = None,
    groups: Optional[Union[str, List[str]]] = None,
    out_file: Optional[Union[str, Path]] = None,
    z_cutoff: Optional[float] = None
) -> None:
    
    # Create data class that holds input data
    ai_files = WaspAnalysisSC(adata_file=count_file,
                              bc_map=bc_map,
                              min_count=min_count,
                              pseudocount=pseudocount,
                              phased=phase,
                              sample=sample,
                              groups=groups,
                              model="single",
                              out_file=out_file,
                              z_cutoff=z_cutoff
                              )
    
    adata_inputs = process_adata_inputs(ad.read_h5ad(ai_files.adata_file), ai_files=ai_files)
    
    
    # print(*vars(ai_files).items(), sep="\n") # For debugging
    # print(adata_inputs) # For debugging
    
    # Update class attributes
    ai_files.update_data(adata_inputs)
    
    # adata = adata_inputs.adata # Hold parsed adata file obj in memory
    
    # Prefilter and hold adata data in memory
    adata = adata_count_qc(adata_inputs.adata,
                           z_cutoff=ai_files.z_cutoff,
                           gt_error=None
                           )

    # Type narrowing: after update_data, these values should be properly set
    assert ai_files.min_count is not None, "min_count should be set in __init__"
    assert ai_files.pseudocount is not None, "pseudocount should be set in __init__"
    assert ai_files.phased is not None, "phased should be set by process_adata_inputs"
    assert isinstance(ai_files.groups, list), "groups should be a list after update_data"

    # Create dictionary of resulting dataframes
    df_dict = get_imbalance_sc(adata,
                               min_count=ai_files.min_count,
                               pseudocount=ai_files.pseudocount,
                               phased=ai_files.phased,
                               sample=ai_files.sample,
                               groups=ai_files.groups)
    
    # Write outputs
    out_path = Path(ai_files.out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    for key, value in df_dict.items():
        group_out_file = out_path / f"{ai_files.prefix}_{key.replace('/', '-')}.tsv"
        
        value.sort_values(by="pval", ascending=True).to_csv(
            group_out_file, sep="\t", header=True, index=False)
        
    print(
        (f"Allelic Imbalance measured for {len(df_dict)} groups!\n"
        f"Results written to: {out_path}/{ai_files.prefix}_[GROUP].tsv")
    )
    