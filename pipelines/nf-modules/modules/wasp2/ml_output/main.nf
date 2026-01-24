process WASP2_ML_OUTPUT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://ghcr.io/jaureguy760/wasp2:latest' :
        'ghcr.io/jaureguy760/wasp2:latest' }"

    input:
    tuple val(meta), path(counts)
    val(output_format)  // comma-separated: "zarr,parquet,anndata" (or "h5ad")

    output:
    tuple val(meta), path("*.zarr", type: 'dir'), emit: zarr, optional: true
    tuple val(meta), path("*.parquet")          , emit: parquet, optional: true
    tuple val(meta), path("*.h5ad")             , emit: anndata, optional: true
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = (task.ext.prefix ?: "${meta.id}").replaceAll(/[^a-zA-Z0-9._-]/, '_')
    """
    #!/usr/bin/env python3
    import pandas as pd
    import numpy as np
    import sys

    # Configuration
    counts_file = "${counts}"
    prefix = "${prefix}"
    sample_id = "${meta.id}"
    VALID_FORMATS = {'zarr', 'parquet', 'anndata', 'h5ad'}

    # Validate and parse output formats
    formats = [f.strip() for f in "${output_format}".lower().split(',') if f.strip()]
    if not formats:
        print("ERROR: No output formats specified", file=sys.stderr)
        sys.exit(1)

    unknown_formats = set(formats) - VALID_FORMATS
    if unknown_formats:
        print(f"ERROR: Unknown output formats: {unknown_formats}", file=sys.stderr)
        print(f"Valid formats: {VALID_FORMATS}", file=sys.stderr)
        sys.exit(1)

    # Validate library availability for requested formats
    if 'zarr' in formats:
        try:
            import zarr
        except ImportError as e:
            print(f"ERROR: zarr format requested but library unavailable: {e}", file=sys.stderr)
            sys.exit(1)

    if 'anndata' in formats or 'h5ad' in formats:
        try:
            import anndata as ad
            import scipy.sparse as sp
        except ImportError as e:
            print(f"ERROR: anndata format requested but library unavailable: {e}", file=sys.stderr)
            sys.exit(1)

    # Read WASP2 counts TSV with error handling
    try:
        df = pd.read_csv(counts_file, sep='\\t')
    except FileNotFoundError:
        print(f"ERROR: Input file not found: {counts_file}", file=sys.stderr)
        sys.exit(1)
    except pd.errors.EmptyDataError:
        print(f"ERROR: Input file is empty: {counts_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Failed to read input file '{counts_file}': {e}", file=sys.stderr)
        sys.exit(1)

    # Validate required columns exist - fail fast on malformed input
    required_cols = ['chrom', 'pos', 'ref', 'alt', 'ref_count', 'alt_count']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"ERROR: Input file missing required columns: {missing_cols}", file=sys.stderr)
        print(f"Found columns: {list(df.columns)}", file=sys.stderr)
        print("This may indicate upstream WASP2 process failure.", file=sys.stderr)
        sys.exit(1)

    # Validate data content
    if len(df) == 0:
        print(f"ERROR: Input file contains no data rows: {counts_file}", file=sys.stderr)
        sys.exit(1)

    # Validate numeric columns
    for col in ['ref_count', 'alt_count']:
        if not pd.api.types.is_numeric_dtype(df[col]):
            print(f"ERROR: Column '{col}' contains non-numeric data", file=sys.stderr)
            sys.exit(1)
        if (df[col] < 0).any():
            print(f"ERROR: Column '{col}' contains negative values", file=sys.stderr)
            sys.exit(1)

    # Compute derived columns
    df['total_count'] = df['ref_count'] + df['alt_count']

    # Handle zero-count variants with logging
    zero_count_mask = df['total_count'] == 0
    n_zero = zero_count_mask.sum()
    if n_zero > 0:
        print(f"Warning: {n_zero} variants have zero total count, ref_ratio set to NaN", file=sys.stderr)
    df['ref_ratio'] = np.where(df['total_count'] > 0, df['ref_count'] / df['total_count'], np.nan)

    df['hap1_count'] = df['ref_count']
    df['hap2_count'] = df['alt_count']
    n_variants = len(df)

    # Zarr output (GenVarLoader compatible)
    if 'zarr' in formats:
        try:
            z = zarr.open(f"{prefix}.zarr", mode='w')
            for col, dtype in [('chrom', str), ('pos', 'i8'), ('ref', str), ('alt', str),
                               ('ref_count', 'i4'), ('alt_count', 'i4'), ('hap1_count', 'i4'),
                               ('hap2_count', 'i4'), ('total_count', 'i4'), ('ref_ratio', 'f4')]:
                data = df[col].values.astype(dtype) if dtype == str else df[col].values
                z.create_dataset(col, data=data, chunks=True, dtype=dtype if dtype != str else None)
            z.attrs.update({'sample_id': sample_id, 'format': 'wasp2_genvarloader', 'version': '1.0'})
            print(f"Created {prefix}.zarr with {n_variants} variants", file=sys.stderr)
        except Exception as e:
            print(f"ERROR: Failed to create Zarr output: {e}", file=sys.stderr)
            sys.exit(1)

    # Parquet output (Polars/DuckDB compatible)
    if 'parquet' in formats:
        try:
            df.to_parquet(f"{prefix}.parquet", index=False, compression='snappy')
            print(f"Created {prefix}.parquet with {n_variants} variants", file=sys.stderr)
        except Exception as e:
            print(f"ERROR: Failed to create Parquet output: {e}", file=sys.stderr)
            sys.exit(1)

    # AnnData output (scverse compatible)
    if 'anndata' in formats or 'h5ad' in formats:
        try:
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
        except Exception as e:
            print(f"ERROR: Failed to create AnnData output: {e}", file=sys.stderr)
            sys.exit(1)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wasp2: \$(python -c "import wasp2; print(wasp2.__version__)" 2>/dev/null || echo "dev")
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        zarr: \$(python -c "import zarr; print(zarr.__version__)" 2>/dev/null || echo "N/A")
        anndata: \$(python -c "import anndata; print(anndata.__version__)" 2>/dev/null || echo "N/A")
    END_VERSIONS
    """

    stub:
    def prefix = (task.ext.prefix ?: "${meta.id}").replaceAll(/[^a-zA-Z0-9._-]/, '_')
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
