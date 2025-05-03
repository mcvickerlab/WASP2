from pathlib import Path
from typing import List, Optional, Union

import anndata as ad
import pandas as pd

from wasp2.analysis.as_analysis_sc import adata_count_qc
from wasp2.analysis.run_analysis_sc import WaspAnalysisSC, process_adata_inputs
from wasp2.analysis.compare_ai import get_compared_imbalance


def run_ai_comparison(
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
    Run allelic imbalance comparisons between groups using single-cell data.

    This function creates a WaspAnalysisSC instance from the provided parameters, processes the
    AnnData input, applies quality control filtering, and computes allelic imbalance comparisons
    between groups. The resulting comparison DataFrames are written as TSV files in the designated
    output directory.

    Parameters
    ----------
    count_file : str
        Path to the AnnData file (e.g., a .h5ad file) containing count data.
    bc_map : str
        Path to the barcode mapping file.
    min_count : int, optional
        Minimum allele count threshold. Defaults to the value specified within WaspAnalysisSC if None.
    pseudocount : int, optional
        Pseudocount added to avoid division by zero. Defaults to the value specified within WaspAnalysisSC if None.
    phase : bool, optional
        Flag indicating whether genotype data is phased. Defaults to the value specified within WaspAnalysisSC if None.
    sample : str, optional
        Sample identifier.
    groups : str or list, optional
        Groups to analyze. Can be a comma-delimited string or a list of group names.
    out_file : str, optional
        Output file path for analysis results. If not provided, defaults to "ai_results.tsv" in the current directory.
    z_cutoff : int or float, optional
        Z-score cutoff for quality control filtering of the AnnData object.

    Returns
    -------
    None

    Notes
    -----
    The function performs the following steps:
      1. Creates a WaspAnalysisSC instance with the specified parameters.
      2. Processes the AnnData input using `process_adata_inputs` and updates the instance attributes.
      3. Applies quality control filtering using `adata_count_qc`.
      4. Computes allelic imbalance comparisons using `get_compared_imbalance`.
      5. Writes the comparison results as TSV files in the designated output directory.

    Examples
    --------
    >>> run_ai_comparison("data.h5ad", "barcode_map.tsv", min_count=10, pseudocount=1,
    ...                   phase=False, sample="Sample1", groups="Group1,Group2",
    ...                   out_file="results.tsv", z_cutoff=4)
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

    print(*vars(ai_files).items(), sep="\n")  # For debugging
    print(adata_inputs)  # For debugging

    # Update class attributes
    ai_files.update_data(adata_inputs)

    # Prefilter and hold AnnData data in memory
    adata = adata_count_qc(adata_inputs.adata, z_cutoff=ai_files.z_cutoff, gt_error=None)

    # Compute allelic imbalance comparisons between groups
    df_dict = get_compared_imbalance(
        adata,
        min_count=ai_files.min_count,
        pseudocount=ai_files.pseudocount,
        phased=ai_files.phased,
        sample=ai_files.sample,
        groups=ai_files.groups,
    )

    # Write outputs
    out_path = Path(ai_files.out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    compared_set = set()

    for key, value in df_dict.items():
        compared_set.update(key)
        compare_out_file = out_path / f"{ai_files.prefix}_{'_'.join(key).replace('/', '-')}.tsv"
        value.sort_values(by="pval", ascending=True).to_csv(compare_out_file, sep="\t", header=True, index=False)

    print(
        f"Performed {len(df_dict)} allelic imbalance comparisons between {len(compared_set)} groups!\n"
        f"Results written to: {out_path}/{ai_files.prefix}_[GROUP1]_[GROUP2].tsv"
    )
