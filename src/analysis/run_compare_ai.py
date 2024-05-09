from pathlib import Path

import anndata as ad
import pandas as pd

from run_analysis_sc import WaspAnalysisSC, process_adata_inputs
from compare_ai import get_compared_imbalance

def run_ai_comparison(count_file,
                       bc_map,
                       min_count=None,
                       pseudocount=None,
                       phase=None,
                       sample=None,
                       groups=None,
                       out_file=None):
    
    
    # Might be smart to change some of the defaults in the class
    # Create data class that holds input data
    ai_files = WaspAnalysisSC(adata_file=count_file,
                              bc_map=bc_map,
                              min_count=min_count,
                              pseudocount=pseudocount,
                              phased=phase,
                              sample=sample,
                              groups=groups,
                              model="single",
                              out_file=out_file
                              )
    
    adata_inputs = process_adata_inputs(ad.read_h5ad(ai_files.adata_file), ai_files=ai_files)
    
    
    print(*vars(ai_files).items(), sep="\n") # For debugging
    print(adata_inputs) # For debugging
    
    # Update class attributes
    ai_files.update_data(adata_inputs)
    
    adata = adata_inputs.adata # Hold parsed adata file obj in memory
    
    df_dict = get_compared_imbalance(adata,
                                     min_count=ai_files.min_count,
                                     pseudocount=ai_files.pseudocount,
                                     phased=ai_files.phased,
                                     sample=ai_files.sample,
                                     groups=ai_files.groups)
    
    # Write outputs
    out_path = Path(ai_files.out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    compared_set = set()

    for key, value in df_dict.items():
        compared_set.update(key)
        
        compare_out_file = out_path / f"{ai_files.prefix}_{'_'.join(key).replace('/', '-')}.tsv"

        value.sort_values(by="pval", ascending=True).to_csv(
            compare_out_file, sep="\t", header=True, index=False)

    print(
        (f"Performed {len(df_dict)} allelic imbalance comparisons between {len(compared_set)} groups!\n"
         f"Results written to: {out_path}/{ai_files.prefix}_[GROUP1]_[GROUP2].tsv")
        )
