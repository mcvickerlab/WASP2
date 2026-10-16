# nf-scatac output

Fragment inputs provide overlap coverage. BAM inputs provide observed reference/alternate allele counts and optional haplotype layers. **This workflow produces no allelic imbalance results:** fragment coverage contains no observed alleles, and BAM-to-inference aggregation is not yet implemented. The `imbalance` output channel is empty.

## Fragment coverage

`fragment_coverage/{sample}_fragment_coverage.tsv` contains:

| Column | Meaning |
|--------|---------|
| `barcode` | Cell barcode |
| `chrom`, `pos` | SNP chromosome and 1-based position |
| `ref`, `alt` | Variant annotation from the VCF; these do not identify an observed allele |
| `overlap_count` | Fragment multiplicity summed over overlaps at this SNP and cell |

`pseudobulk/{sample}_pseudobulk_coverage.tsv` aggregates the above rows across cells and contains `chrom`, `pos`, `ref`, `alt`, `overlap_count`, and `n_cells`. SNPs must meet `--min_cells_per_snp`. These files are **not** allele-count input for `find-imbalance`.

`count_stats/{sample}_count_stats.tsv` and `pseudobulk/{sample}_aggregation_stats.tsv` contain counting and filtering summaries.

## AnnData

Fragment inputs produce `{sample}_coverage.h5ad` (and optional `_coverage.zarr`). `X` contains overlap coverage with cells in `obs` and SNPs in `var`; allele-specific layers are absent. Cell QC is published under `cell_qc/`.

BAM inputs produce allele-count AnnData with `ref` and `alt` layers. The haplotype-layer module uses VCF phase information for `hap1`/`hap2`. These counts can be passed separately to the WASP2 single-cell analysis CLI after supplying validated cell-group metadata.

## Other outputs

`variants/` contains the VCF-derived SNP BED. `pipeline_info/` contains execution reports and traces. Consult the trace for the processes actually run for each input path.

```python
import pandas as pd

coverage = pd.read_csv("results/fragment_coverage/sample_fragment_coverage.tsv", sep="\t")
cell_by_snp = coverage.pivot_table(
    index="barcode", columns=["chrom", "pos"], values="overlap_count", fill_value=0
)
```

An empty coverage table can indicate barcode mismatch, chromosome naming differences, no SNP overlaps, or filtering thresholds that remove all observations. Coverage alone must not be interpreted as reference-allele support.
