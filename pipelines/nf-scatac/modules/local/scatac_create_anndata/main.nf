/*
 * SCATAC_CREATE_ANNDATA - Convert allele counts to AnnData H5AD format
 *
 * Creates AnnData object with per-cell allele counts at heterozygous SNPs.
 * Outputs:
 *   - H5AD file with sparse matrices in layers (X=total, ref, alt)
 *   - Cell metadata (obs) with QC metrics
 *   - Variant metadata (var) with SNP annotations
 *   - Optional Zarr format for GenVarLoader integration
 */

process SCATAC_CREATE_ANNDATA {
    tag "$meta.id"
    label 'process_medium'

    conda "${projectDir}/../nf-modules/modules/wasp2/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://jaureguy760/wasp2:latest' :
        'jaureguy760/wasp2:latest' }"

    input:
    tuple val(meta), path(counts)
    val(create_zarr)

    output:
    tuple val(meta), path("*.h5ad")           , emit: anndata
    tuple val(meta), path("*.zarr")           , emit: zarr, optional: true
    tuple val(meta), path("*_cell_qc.tsv")    , emit: cell_qc
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def zarr_flag = create_zarr ? "True" : "False"
    """
    python3 << 'PYEOF'
import pandas as pd
import numpy as np
from scipy import sparse
import anndata as ad
import warnings
warnings.filterwarnings('ignore')

# Read allele counts TSV
# Expected columns: barcode, chrom, pos, ref, alt, overlap_count
df = pd.read_csv("${counts}", sep='\\t')

if df.empty:
    # Create minimal empty AnnData for edge case
    adata = ad.AnnData(
        X=sparse.csr_matrix((0, 0)),
        obs=pd.DataFrame(index=[]),
        var=pd.DataFrame(index=[])
    )
    adata.write_h5ad("${prefix}_allelic.h5ad")
    pd.DataFrame(columns=['barcode', 'n_snps', 'total_counts']).to_csv(
        "${prefix}_cell_qc.tsv", sep='\\t', index=False
    )
else:
    # Create unique SNP identifiers
    df['snp_id'] = df['chrom'].astype(str) + ':' + df['pos'].astype(str) + ':' + df['ref'] + '>' + df['alt']

    # Get unique barcodes and SNPs
    barcodes = df['barcode'].unique()
    snp_ids = df['snp_id'].unique()
    n_cells = len(barcodes)
    n_snps = len(snp_ids)

    # Create index mappings
    barcode_to_idx = {bc: i for i, bc in enumerate(barcodes)}
    snp_to_idx = {snp: i for i, snp in enumerate(snp_ids)}

    # Build sparse matrix (cells x SNPs) with overlap counts
    row_idx = df['barcode'].map(barcode_to_idx).values
    col_idx = df['snp_id'].map(snp_to_idx).values
    data = df['overlap_count'].values

    X = sparse.csr_matrix((data, (row_idx, col_idx)), shape=(n_cells, n_snps))

    # Build cell metadata (obs)
    cell_stats = df.groupby('barcode').agg({
        'snp_id': 'nunique',
        'overlap_count': 'sum'
    }).rename(columns={'snp_id': 'n_snps', 'overlap_count': 'total_counts'})
    cell_stats = cell_stats.reindex(barcodes)

    obs = pd.DataFrame({
        'n_snps': cell_stats['n_snps'].values,
        'total_counts': cell_stats['total_counts'].values,
        'chemistry': '${meta.chemistry ?: "10x-atac-v2"}',
        'sample_id': '${meta.id}'
    }, index=barcodes)
    obs.index.name = 'barcode'

    # Build variant metadata (var)
    snp_info = df.drop_duplicates('snp_id').set_index('snp_id')[['chrom', 'pos', 'ref', 'alt']]
    snp_info = snp_info.reindex(snp_ids)

    var = pd.DataFrame({
        'chrom': snp_info['chrom'].values,
        'pos': snp_info['pos'].values,
        'ref': snp_info['ref'].values,
        'alt': snp_info['alt'].values
    }, index=snp_ids)
    var.index.name = 'snp_id'

    # Create AnnData object
    adata = ad.AnnData(
        X=X,
        obs=obs,
        var=var
    )

    # Add metadata
    adata.uns['sample_id'] = '${meta.id}'
    adata.uns['pipeline'] = 'nf-scatac'
    adata.uns['data_type'] = 'scATAC_allelic_counts'

    # Write H5AD
    adata.write_h5ad("${prefix}_allelic.h5ad")

    # Write cell QC metrics
    obs.reset_index().to_csv("${prefix}_cell_qc.tsv", sep='\\t', index=False)

    # Optionally write Zarr for GenVarLoader
    if ${zarr_flag}:
        adata.write_zarr("${prefix}_allelic.zarr")

print(f"Created AnnData: {n_cells} cells x {n_snps} SNPs")
PYEOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        anndata: \$(python -c "import anndata; print(anndata.__version__)")
        scipy: \$(python -c "import scipy; print(scipy.__version__)")
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Create minimal stub H5AD
    python3 << 'PYEOF'
import pandas as pd
import numpy as np
from scipy import sparse
import anndata as ad

# Minimal test data
adata = ad.AnnData(
    X=sparse.csr_matrix([[5, 3], [2, 8]]),
    obs=pd.DataFrame({'n_snps': [2, 2], 'total_counts': [8, 10]}, index=['AAACGAAC-1', 'AAACGAAT-1']),
    var=pd.DataFrame({'chrom': ['chr1', 'chr1'], 'pos': [100, 200]}, index=['chr1:100:A>G', 'chr1:200:C>T'])
)
adata.write_h5ad("${prefix}_allelic.h5ad")

pd.DataFrame({
    'barcode': ['AAACGAAC-1', 'AAACGAAT-1'],
    'n_snps': [2, 2],
    'total_counts': [8, 10]
}).to_csv("${prefix}_cell_qc.tsv", sep='\\t', index=False)
PYEOF

    if [ "${create_zarr}" == "true" ]; then
        mkdir -p ${prefix}_allelic.zarr
        touch ${prefix}_allelic.zarr/.zgroup
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
        anndata: stub
        scipy: stub
        pandas: stub
    END_VERSIONS
    """
}
