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
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def formats = output_format.split(',').collect { it.trim().toLowerCase() }
    def zarr_cmd = formats.contains('zarr') ? "convert_to_zarr ${counts} ${prefix}.zarr" : ''
    def parquet_cmd = formats.contains('parquet') ? "convert_to_parquet ${counts} ${prefix}.parquet" : ''
    def anndata_cmd = formats.contains('anndata') || formats.contains('h5ad') ? "convert_to_anndata ${counts} ${prefix}.h5ad" : ''
    """
    #!/usr/bin/env python3
    import pandas as pd
    import numpy as np
    import sys
    from pathlib import Path

    # Read WASP2 counts TSV
    counts_file = "${counts}"
    prefix = "${prefix}"
    formats = "${output_format}".lower().split(',')

    df = pd.read_csv(counts_file, sep='\\t')

    # Ensure required columns exist
    required_cols = ['chrom', 'pos', 'ref', 'alt', 'ref_count', 'alt_count']
    for col in required_cols:
        if col not in df.columns:
            print(f"Warning: Missing column {col}, creating empty", file=sys.stderr)
            df[col] = 0 if 'count' in col else ''

    # Add computed columns
    df['total_count'] = df['ref_count'] + df['alt_count']
    df['ref_ratio'] = np.where(df['total_count'] > 0, df['ref_count'] / df['total_count'], 0.5)
    df['hap1_count'] = df['ref_count']  # Default: ref = hap1
    df['hap2_count'] = df['alt_count']  # Default: alt = hap2

    # Convert to Zarr (GenVarLoader compatible)
    if 'zarr' in formats:
        import zarr
        z = zarr.open(f"{prefix}.zarr", mode='w')
        # Store as structured array
        z.create_dataset('chrom', data=df['chrom'].values.astype(str), chunks=True)
        z.create_dataset('pos', data=df['pos'].values, chunks=True, dtype='i8')
        z.create_dataset('ref', data=df['ref'].values.astype(str), chunks=True)
        z.create_dataset('alt', data=df['alt'].values.astype(str), chunks=True)
        z.create_dataset('ref_count', data=df['ref_count'].values, chunks=True, dtype='i4')
        z.create_dataset('alt_count', data=df['alt_count'].values, chunks=True, dtype='i4')
        z.create_dataset('hap1_count', data=df['hap1_count'].values, chunks=True, dtype='i4')
        z.create_dataset('hap2_count', data=df['hap2_count'].values, chunks=True, dtype='i4')
        z.create_dataset('total_count', data=df['total_count'].values, chunks=True, dtype='i4')
        z.create_dataset('ref_ratio', data=df['ref_ratio'].values, chunks=True, dtype='f4')
        # Add metadata
        z.attrs['sample_id'] = "${meta.id}"
        z.attrs['format'] = 'wasp2_genvarloader'
        z.attrs['version'] = '1.0'
        print(f"Created {prefix}.zarr with {len(df)} variants", file=sys.stderr)

    # Convert to Parquet (Polars/DuckDB compatible)
    if 'parquet' in formats:
        df.to_parquet(f"{prefix}.parquet", index=False, compression='snappy')
        print(f"Created {prefix}.parquet with {len(df)} variants", file=sys.stderr)

    # Convert to AnnData/H5AD (scverse compatible)
    if 'anndata' in formats or 'h5ad' in formats:
        import anndata as ad
        import scipy.sparse as sp

        # For single-sample: create 1-row AnnData
        # X = total counts, layers = ref/alt/hap1/hap2
        n_vars = len(df)
        X = sp.csr_matrix(df['total_count'].values.reshape(1, -1))

        # Create obs (sample metadata)
        obs = pd.DataFrame({'sample_id': ["${meta.id}"]}, index=["${meta.id}"])

        # Create var (variant metadata)
        var = pd.DataFrame({
            'chrom': df['chrom'].values,
            'pos': df['pos'].values,
            'ref': df['ref'].values,
            'alt': df['alt'].values,
            'region': df.get('region', pd.Series([''] * n_vars)).values,
        }, index=[f"{r.chrom}_{r.pos}_{r.ref}_{r.alt}" for r in df.itertuples()])

        adata = ad.AnnData(X=X, obs=obs, var=var)

        # Add layers for allele counts
        adata.layers['ref'] = sp.csr_matrix(df['ref_count'].values.reshape(1, -1))
        adata.layers['alt'] = sp.csr_matrix(df['alt_count'].values.reshape(1, -1))
        adata.layers['hap1'] = sp.csr_matrix(df['hap1_count'].values.reshape(1, -1))
        adata.layers['hap2'] = sp.csr_matrix(df['hap2_count'].values.reshape(1, -1))

        # Store ref_ratio in obsm for convenience
        adata.obsm['ref_ratio'] = df['ref_ratio'].values.reshape(1, -1)

        adata.write_h5ad(f"{prefix}.h5ad")
        print(f"Created {prefix}.h5ad with 1 sample x {n_vars} variants", file=sys.stderr)

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
    def formats = output_format.split(',').collect { it.trim().toLowerCase() }
    def zarr_cmd = formats.contains('zarr') ? "mkdir -p ${prefix}.zarr && touch ${prefix}.zarr/.zarray" : ''
    def parquet_cmd = formats.contains('parquet') ? "touch ${prefix}.parquet" : ''
    def anndata_cmd = formats.contains('anndata') || formats.contains('h5ad') ? "touch ${prefix}.h5ad" : ''
    """
    ${zarr_cmd}
    ${parquet_cmd}
    ${anndata_cmd}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wasp2: dev
        pandas: stub
        zarr: stub
        anndata: stub
    END_VERSIONS
    """
}
