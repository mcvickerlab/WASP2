/*
 * WASP2_COUNT_SC - Single-cell allele-specific variant counting
 *
 * Counts allele-specific reads at heterozygous SNPs for single-cell data.
 * Uses cell barcodes from BAM tags to assign counts to individual cells.
 * Outputs H5AD with per-cell ref/alt allele counts.
 */

process WASP2_COUNT_SC {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/../environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wasp2:1.2.1--pyhdfd78af_0' :
        'biocontainers/wasp2:1.2.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(vcf), path(vcf_index)
    path(barcodes)  // Cell barcodes file (one barcode per line)
    path(features)  // Optional: BED file with regions to restrict analysis

    output:
    tuple val(meta), path("*.h5ad")        , emit: counts
    tuple val(meta), path("*_stats.tsv")   , emit: stats, optional: true
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = (task.ext.prefix ?: "${meta.id}").replaceAll(/[^a-zA-Z0-9._-]/, '_')
    def sample = meta.sample?.toString()?.replaceAll(/[^a-zA-Z0-9._-]/, '_') ?: ''
    def sample_arg = sample ? "-s ${sample}" : ''
    def feature_arg = features.name != 'NO_FILE' ? "-f ${features}" : ''
    """
    set -euo pipefail

    wasp2-count \\
        count-variants-sc \\
        ${bam} \\
        ${vcf} \\
        ${barcodes} \\
        ${sample_arg} \\
        ${feature_arg} \\
        -o ${prefix}_allele_counts.h5ad \\
        ${args}

    # Validate output was created
    if [ ! -f "${prefix}_allele_counts.h5ad" ]; then
        echo "ERROR: Output file ${prefix}_allele_counts.h5ad was not created" >&2
        exit 1
    fi

    # Generate counting statistics
    python3 << 'PYEOF'
import anndata as ad
import pandas as pd

adata = ad.read_h5ad("${prefix}_allele_counts.h5ad")
stats = {
    'metric': ['total_cells', 'total_snps', 'total_ref_counts', 'total_alt_counts'],
    'value': [
        adata.n_obs,
        adata.n_vars,
        int(adata.layers.get('ref', adata.X).sum()) if 'ref' in adata.layers else 0,
        int(adata.layers.get('alt', adata.X).sum()) if 'alt' in adata.layers else 0
    ]
}
pd.DataFrame(stats).to_csv("${prefix}_stats.tsv", sep='\\t', index=False)
PYEOF

    # Get version
    WASP2_VERSION=\$(wasp2-count --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+' || echo "unknown")

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wasp2: \${WASP2_VERSION}
        python: \$(python --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+')
        anndata: \$(python -c "import anndata; print(anndata.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = (task.ext.prefix ?: "${meta.id}").replaceAll(/[^a-zA-Z0-9._-]/, '_')
    """
    python3 << 'PYEOF'
import numpy as np
import pandas as pd
from scipy import sparse
import anndata as ad

# Create stub AnnData with allele-specific layers
n_cells, n_snps = 10, 50
X = sparse.random(n_cells, n_snps, density=0.3, format='csr')
ref_counts = sparse.random(n_cells, n_snps, density=0.3, format='csr')
alt_counts = sparse.random(n_cells, n_snps, density=0.3, format='csr')

obs = pd.DataFrame({
    'n_snps': np.random.randint(10, 50, n_cells),
    'total_ref': np.random.randint(100, 1000, n_cells),
    'total_alt': np.random.randint(100, 1000, n_cells)
}, index=[f'AAACGAAC-{i}' for i in range(n_cells)])

var = pd.DataFrame({
    'chrom': ['chr1'] * n_snps,
    'pos': range(100000, 100000 + n_snps * 1000, 1000),
    'ref': ['A'] * n_snps,
    'alt': ['G'] * n_snps
}, index=[f'chr1:{100000 + i*1000}:A>G' for i in range(n_snps)])

adata = ad.AnnData(X=X, obs=obs, var=var)
adata.layers['ref'] = ref_counts
adata.layers['alt'] = alt_counts
adata.write_h5ad("${prefix}_allele_counts.h5ad")
PYEOF

    echo -e "metric\\tvalue" > ${prefix}_stats.tsv
    echo -e "total_cells\\t10" >> ${prefix}_stats.tsv
    echo -e "total_snps\\t50" >> ${prefix}_stats.tsv
    echo -e "total_ref_counts\\t500" >> ${prefix}_stats.tsv
    echo -e "total_alt_counts\\t480" >> ${prefix}_stats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wasp2: stub
        python: stub
        anndata: stub
    END_VERSIONS
    """
}
