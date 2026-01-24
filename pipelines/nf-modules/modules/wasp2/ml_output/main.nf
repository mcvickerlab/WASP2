process WASP2_ML_OUTPUT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/../environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://jaureguy760/wasp2:latest' :
        'jaureguy760/wasp2:latest' }"

    input:
    tuple val(meta), path(counts)
    val(output_format)  // comma-separated: "zarr,parquet,anndata"

    output:
    tuple val(meta), path("*.zarr", type: 'dir'), emit: zarr, optional: true
    tuple val(meta), path("*.parquet")          , emit: parquet, optional: true
    tuple val(meta), path("*.h5ad")             , emit: anndata, optional: true
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python3
    import pandas as pd
    import numpy as np
    import sys

    # Configuration
    counts_file = "${counts}"
    prefix = "${prefix}"
    sample_id = "${meta.id}"
    formats = "${output_format}".lower().split(',')

    # Read WASP2 counts TSV
    df = pd.read_csv(counts_file, sep='\\t')

    # Ensure required columns exist
    for col in ['chrom', 'pos', 'ref', 'alt', 'ref_count', 'alt_count']:
        if col not in df.columns:
            print(f"Warning: Missing column {col}, creating empty", file=sys.stderr)
            df[col] = 0 if 'count' in col else ''

    # Compute derived columns
    df['total_count'] = df['ref_count'] + df['alt_count']
    df['ref_ratio'] = np.where(df['total_count'] > 0, df['ref_count'] / df['total_count'], 0.5)
    df['hap1_count'] = df['ref_count']
    df['hap2_count'] = df['alt_count']

    n_variants = len(df)

    # Zarr output (GenVarLoader compatible)
    if 'zarr' in formats:
        import zarr
        z = zarr.open(f"{prefix}.zarr", mode='w')
        for col, dtype in [('chrom', str), ('pos', 'i8'), ('ref', str), ('alt', str),
                           ('ref_count', 'i4'), ('alt_count', 'i4'), ('hap1_count', 'i4'),
                           ('hap2_count', 'i4'), ('total_count', 'i4'), ('ref_ratio', 'f4')]:
            data = df[col].values.astype(dtype) if dtype == str else df[col].values
            z.create_dataset(col, data=data, chunks=True, dtype=dtype if dtype != str else None)
        z.attrs.update({'sample_id': sample_id, 'format': 'wasp2_genvarloader', 'version': '1.0'})
        print(f"Created {prefix}.zarr with {n_variants} variants", file=sys.stderr)

    # Parquet output (Polars/DuckDB compatible)
    if 'parquet' in formats:
        df.to_parquet(f"{prefix}.parquet", index=False, compression='snappy')
        print(f"Created {prefix}.parquet with {n_variants} variants", file=sys.stderr)

    # AnnData output (scverse compatible)
    if 'anndata' in formats or 'h5ad' in formats:
        import anndata as ad
        import scipy.sparse as sp

        X = sp.csr_matrix(df['total_count'].values.reshape(1, -1))
        obs = pd.DataFrame({'sample_id': [sample_id]}, index=[sample_id])
        var = pd.DataFrame({
            'chrom': df['chrom'].values, 'pos': df['pos'].values,
            'ref': df['ref'].values, 'alt': df['alt'].values,
            'region': df.get('region', pd.Series([''] * n_variants)).values,
        }, index=[f"{r.chrom}_{r.pos}_{r.ref}_{r.alt}" for r in df.itertuples()])

        adata = ad.AnnData(X=X, obs=obs, var=var)
        for layer in ['ref_count', 'alt_count', 'hap1_count', 'hap2_count']:
            adata.layers[layer.replace('_count', '')] = sp.csr_matrix(df[layer].values.reshape(1, -1))
        adata.obsm['ref_ratio'] = df['ref_ratio'].values.reshape(1, -1)
        adata.write_h5ad(f"{prefix}.h5ad")
        print(f"Created {prefix}.h5ad with 1 sample x {n_variants} variants", file=sys.stderr)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wasp2: \$(python -c "import wasp2; print(wasp2.__version__)" 2>/dev/null || echo "dev")
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        zarr: \$(python -c "import zarr; print(zarr.__version__)" 2>/dev/null || echo "N/A")
        anndata: \$(python -c "import anndata; print(anndata.__version__)" 2>/dev/null || echo "N/A")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def formats = output_format.toLowerCase().split(',')
    """
    ${formats.contains('zarr') ? "mkdir -p ${prefix}.zarr && touch ${prefix}.zarr/.zarray" : ''}
    ${formats.contains('parquet') ? "touch ${prefix}.parquet" : ''}
    ${formats.contains('anndata') || formats.contains('h5ad') ? "touch ${prefix}.h5ad" : ''}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wasp2: dev
        pandas: stub
        zarr: stub
        anndata: stub
    END_VERSIONS
    """
}
