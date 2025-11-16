"""
Author: Aaron Ho
Python Version: 3.9
"""

# Default Python package Imports
from pathlib import Path
from csv import DictReader, reader

# External package imports
import pandas as pd

# Local script imports
from as_analysis import get_imbalance



# TODO GOTTA IMPLEMENT THIS

class WaspAnalysisData:

    def __init__(self, count_file,
                 min_count=None,
                 pseudocount=None,
                 phased=None,
                 model=None,
                 out_file=None,
                 region_col=None,
                 groupby=None,
                ):
        
        # User input data
        self.count_file = count_file
        self.min_count = min_count
        self.pseudocount = pseudocount
        self.phased = phased
        self.model = model
        self.out_file = out_file

        # Group by feature by default
        self.region_col = region_col
        self.groupby = groupby # group by region or parent?
        
        # TODO parse vcf for phased instead of default unphased
        if not self.phased:
            self.phased = False


        # Default to single dispersion model
        if ((self.model is None) or 
            (self.model not in {"single", "linear"})):
            
            self.model = "single"
        
        # Default min count of 10 
        if self.min_count is None:
            self.min_count = 10

        if self.pseudocount is None:
            # self.pseudocount = 0 # either 0 or 1 for default
            self.pseudocount = 1
        
        # Read header only for validation
        with open(self.count_file) as f:
            count_cols = next(reader(f, delimiter = "\t"))
        
        # 7 columns at minimum, 10 at maximum
        # 3required : chr, pos, ref, alt
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

            elif ((len(count_cols) > (min_cols+1)) and
                  self.groupby in {count_cols[region_idx+1], "Parent", "parent"}):

                self.groupby = count_cols[region_idx+1] # Set group
            else:
                # Maybe throw error instead
                print(f"{self.groupby} not found in columns \n{count_cols}")
                self.groupby = None


        # Create default outfile 
        if self.out_file is None:
            self.out_file = str(Path.cwd() / "ai_results.tsv") # do this after


def run_ai_analysis(count_file,
                    min_count=None,
                    pseudocount=None,
                    phased=None,
                    model=None,
                    out_file=None,
                    region_col=None,
                    groupby=None):
    
    # Store analysis data and params
    ai_files = WaspAnalysisData(count_file,
                                min_count=min_count,
                                pseudocount=pseudocount,
                                phased=phased,
                                model=model,
                                out_file=out_file,
                                region_col=region_col,
                                groupby=groupby
                                )
    
    # Run analysis pipeline
    ai_df = get_imbalance(ai_files.count_file,
                          min_count=ai_files.min_count,
                          pseudocount=ai_files.pseudocount,
                          method=ai_files.model,
                          phased=ai_files.phased,
                          region_col=ai_files.region_col,
                          groupby=ai_files.groupby
                          )
    
    # Maybe give option to sort or not sort by pval
    ai_df = ai_df.sort_values(by="fdr_pval", ascending=True)
    
    # Write results
    ai_df.to_csv(ai_files.out_file, sep="\t", header=True, index=False)
