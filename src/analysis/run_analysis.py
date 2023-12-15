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
                 model=None,
                 phased=None,
                 out_dir=None,
                 out_file=None,
                 region_col=None,
                 features=None):
        
        # User input data
        self.count_file = count_file
        self.min_count = min_count
        self.model = model
        self.phased = phased # TODO
        self.out_file = out_file
        self.out_dir = out_dir  # should i replace this with out file???
        self.region_col = region_col
        self.features = features # TODO and also add rna-seq support back
        
        # I need to also add other things for single cell back
        

        # Default to single dispersion model
        if self.model is None:
            self.model = "single"
        
        # Default min count of 10 
        if self.min_count is None:
            self.min_count = 10
        
        
        # Automatically parse region col
        # Should i do this after the df is created?
        if self.region_col is None:
            
            # Read header only
            with open(self.count_file) as f:
                count_cols = next(reader(f, delimiter = "\t"))
            
            # Check region_col from file
            if "region" in count_cols:
                self.region_col = "region" # default atac naming
            elif "peak" in count_cols:
                self.region_col = "peak" # from previous implementation
            elif "genes" in count_cols:
                self.region_col = "genes"
            else:
                # SNPs only
                # df["region"] = df["chrom"] + "_" + df["pos"].astype(str)
                self.region_col = "region" # should i name as snp?

                
        # Create default outfile 
        if self.out_file is None:
            self.out_file = str(Path.cwd() / "ai_results.tsv") # do this after


def run_ai_analysis(count_file,
                    min_count=None,
                    model=None,
                    out_file=None,
                    region_col=None):
    
    # Store analysis data and params
    ai_files = WaspAnalysisData(count_file,
                                min_count=min_count,
                                model=model,
                                phased=None,
                                out_dir=None,
                                out_file=out_file,
                                region_col=region_col,
                                features=None)
    
    # Run analysis pipeline
    ai_df = get_imbalance(ai_files.count_file,
                          min_count=ai_files.min_count,
                          method=ai_files.model,
                          region_col=ai_files.region_col,
                          is_gene=False)

    
    # Write results
    ai_df.to_csv(ai_files.out_file, sep="\t", header=True, index=False)
