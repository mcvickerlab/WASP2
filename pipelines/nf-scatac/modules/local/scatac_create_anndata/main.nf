/*
 * SCATAC_CREATE_ANNDATA - Convert allele counts to AnnData H5AD format
 *
 * Creates AnnData object with per-cell allele counts at heterozygous SNPs.
 * Outputs:
 *   - H5AD file with sparse matrix (X=total overlap counts per cell/SNP)
 *   - Cell metadata (obs) with QC metrics
 *   - Variant metadata (var) with SNP annotations
 *   - Optional Zarr format for GenVarLoader integration
 *
 * Note: Fragment files contain coordinates only, not sequences, so we count
 * total overlaps rather than allele-specific counts. Handles empty input
 * gracefully by creating a minimal valid H5AD with 0x0 matrix.
 */

process SCATAC_CREATE_ANNDATA {
    tag "$meta.id"
    label 'process_medium'

    conda "${projectDir}/../nf-modules/modules/wasp2/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://ghcr.io/jaureguy760/wasp2:latest' :
        'ghcr.io/jaureguy760/wasp2:latest' }"

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
import sys
import pandas as pd
from scipy import sparse
import anndata as ad

# Read allele counts TSV
df = pd.read_csv("${counts}", sep='\\t')

# Validate required columns
required_cols = ['barcode', 'chrom', 'pos', 'ref', 'alt', 'overlap_count']
missing = set(required_cols) - set(df.columns)
if missing:
    print(f"ERROR: Input TSV missing required columns: {missing}", file=sys.stderr)
    print(f"Found columns: {list(df.columns)}", file=sys.stderr)
    sys.exit(1)

# Handle empty input (early return pattern)
if df.empty:
    ad.AnnData(X=sparse.csr_matrix((0, 0))).write_h5ad("${prefix}_allelic.h5ad")
    pd.DataFrame(columns=['barcode', 'n_snps', 'total_counts']).to_csv(
        "${prefix}_cell_qc.tsv", sep='\\t', index=False)
    print("Created empty AnnData (no data)")
else:
    # Build SNP identifiers and get unique values
    df['snp_id'] = df['chrom'].astype(str) + ':' + df['pos'].astype(str) + ':' + df['ref'] + '>' + df['alt']
    barcodes, snp_ids = df['barcode'].unique(), df['snp_id'].unique()

    # Build sparse matrix using categorical codes for efficiency
    row_idx = pd.Categorical(df['barcode'], categories=barcodes).codes
    col_idx = pd.Categorical(df['snp_id'], categories=snp_ids).codes
    X = sparse.csr_matrix((df['overlap_count'].values, (row_idx, col_idx)),
                          shape=(len(barcodes), len(snp_ids)))

    # Build cell metadata directly from groupby
    cell_stats = df.groupby('barcode', sort=False).agg(
        n_snps=('snp_id', 'nunique'),
        total_counts=('overlap_count', 'sum')
    ).reindex(barcodes)
    cell_stats['chemistry'] = '${meta.chemistry ?: "10x-atac-v2"}'
    cell_stats['sample_id'] = '${meta.id}'
    cell_stats.index.name = 'barcode'

    # Build variant metadata
    var = df.drop_duplicates('snp_id').set_index('snp_id')[['chrom', 'pos', 'ref', 'alt']].reindex(snp_ids)
    var.index.name = 'snp_id'

    # Create and write AnnData
    adata = ad.AnnData(X=X, obs=cell_stats, var=var,
                       uns={'sample_id': '${meta.id}', 'pipeline': 'nf-scatac',
                            'data_type': 'scATAC_allelic_counts'})
    adata.write_h5ad("${prefix}_allelic.h5ad")
    cell_stats.reset_index().to_csv("${prefix}_cell_qc.tsv", sep='\\t', index=False)

    if ${zarr_flag}:
        adata.write_zarr("${prefix}_allelic.zarr")

    print(f"Created AnnData: {len(barcodes)} cells x {len(snp_ids)} SNPs")
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
    python3 << 'PYEOF'
import pandas as pd
from scipy import sparse
import anndata as ad

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
