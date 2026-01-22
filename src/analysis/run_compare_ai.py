from pathlib import Path
from typing import Optional, Union, List

import anndata as ad
from anndata import AnnData
import pandas as pd

from .as_analysis_sc import adata_count_qc
from .run_analysis_sc import WaspAnalysisSC, process_adata_inputs, AdataInputs
from .compare_ai import get_compared_imbalance

def run_ai_comparison(
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
    
    
    # Might be smart to change some of the defaults in the class
    # Create data class that holds input data
    ai_files: WaspAnalysisSC = WaspAnalysisSC(
        adata_file=count_file,
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

    adata_inputs: AdataInputs = process_adata_inputs(ad.read_h5ad(ai_files.adata_file), ai_files=ai_files)

    # Update class attributes
    ai_files.update_data(adata_inputs)

    # adata = adata_inputs.adata # Hold parsed adata file obj in memory

    # Prefilter and hold adata data in memory
    adata: AnnData = adata_count_qc(
        adata_inputs.adata,
        z_cutoff=ai_files.z_cutoff,
        gt_error=None
    )

    # After __init__ and update_data, these attributes are guaranteed to be non-None
    assert ai_files.min_count is not None
    assert ai_files.pseudocount is not None
    assert ai_files.phased is not None
    assert isinstance(ai_files.groups, list)

    df_dict: dict[tuple[str, str], pd.DataFrame] = get_compared_imbalance(
        adata,
        min_count=ai_files.min_count,
        pseudocount=ai_files.pseudocount,
        phased=ai_files.phased,
        sample=ai_files.sample,
        groups=ai_files.groups
    )
    
    # Write outputs
    out_path: Path = Path(ai_files.out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    compared_set: set[str] = set()

    for key, value in df_dict.items():
        compared_set.update(key)

        compare_out_file: Path = out_path / f"{ai_files.prefix}_{'_'.join(key).replace('/', '-')}.tsv"

        value.sort_values(by="pval", ascending=True).to_csv(
            compare_out_file, sep="\t", header=True, index=False)

    print(
        (f"Performed {len(df_dict)} allelic imbalance comparisons between {len(compared_set)} groups!\n"
         f"Results written to: {out_path}/{ai_files.prefix}_[GROUP1]_[GROUP2].tsv")
        )
