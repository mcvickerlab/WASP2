/*
 * SCATAC_ADD_HAPLOTYPE_LAYERS - Add haplotype layers to AnnData
 *
 * Takes AnnData with ref/alt layers and phased VCF, adds hap1/hap2 layers.
 * Haplotype assignment is based on VCF phasing: 0|1 means ref=hap1, alt=hap2;
 * 1|0 means alt=hap1, ref=hap2.
 *
 * For unphased SNPs (missing from VCF or without '|' delimiter), default
 * assignment is used: hap1=ref, hap2=alt. This is tracked in:
 *   - adata.var['is_phased']: per-SNP phasing status
 *   - adata.uns['unphased_snps']: count of unphased SNPs
 *   - adata.uns['unphased_default']: documents the default behavior
 *
 * Outputs AnnData with layers: X (total), ref, alt, hap1, hap2
 */

process SCATAC_ADD_HAPLOTYPE_LAYERS {
    tag "$meta.id"
    label 'process_medium'

    conda "${projectDir}/../nf-modules/modules/wasp2/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wasp2:1.2.1--pyhdfd78af_0' :
        'biocontainers/wasp2:1.2.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(anndata)
    tuple val(meta2), path(vcf), path(vcf_index)
    val(create_zarr)

    output:
    tuple val(meta), path("*_with_haplotypes.h5ad"), emit: anndata
    tuple val(meta), path("*.zarr")                , emit: zarr, optional: true
    tuple val(meta), path("*_cell_qc.tsv")         , emit: cell_qc
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sample = meta.sample ?: meta.id
    def zarr_flag = create_zarr ? "True" : "False"
    """
    python3 << 'PYEOF'
import sys
import numpy as np
import pandas as pd
from scipy import sparse
import anndata as ad

try:
    import pysam
except ImportError:
    print("ERROR: pysam is required for haplotype layer creation but is not installed", file=sys.stderr)
    print("Install with: pip install pysam", file=sys.stderr)
    sys.exit(1)

# Read input AnnData
adata = ad.read_h5ad("${anndata}")

# Validate input has ref/alt layers
if 'ref' not in adata.layers or 'alt' not in adata.layers:
    print("ERROR: Input AnnData must have 'ref' and 'alt' layers", file=sys.stderr)
    print(f"Found layers: {list(adata.layers.keys())}", file=sys.stderr)
    sys.exit(1)

# Parse phasing from VCF
sample_name = "${sample}"
phasing = {}  # snp_id -> (hap1_allele, hap2_allele)

vcf = pysam.VariantFile("${vcf}")

# Find sample index - fail if sample not found
try:
    sample_idx = list(vcf.header.samples).index(sample_name)
except ValueError:
    available = list(vcf.header.samples)
    print(f"ERROR: Sample '{sample_name}' not found in VCF.", file=sys.stderr)
    print(f"Available samples: {available}", file=sys.stderr)
    print("Specify the correct sample name in your samplesheet.", file=sys.stderr)
    sys.exit(1)

for rec in vcf:
    gt = rec.samples[sample_idx].get('GT', None)
    if gt is None:
        continue
    # Check if phased (tuple with phased=True or string with '|')
    phased = rec.samples[sample_idx].phased if hasattr(rec.samples[sample_idx], 'phased') else False
    if not phased:
        continue
    # Build SNP ID to match AnnData var index
    snp_id = f"{rec.chrom}:{rec.pos}:{rec.ref}>{rec.alts[0]}"
    # gt[0] is first haplotype, gt[1] is second haplotype
    # 0 = ref, 1 = alt
    hap1_allele = gt[0]  # 0=ref, 1=alt
    hap2_allele = gt[1]
    phasing[snp_id] = (hap1_allele, hap2_allele)
vcf.close()

# Create hap1/hap2 layers based on phasing
# IMPORTANT: Unphased SNPs default to ref=hap1, alt=hap2 (arbitrary assignment).
# This is tracked in adata.uns['unphased_snps'] and adata.var['is_phased'].
ref_matrix = adata.layers['ref']
alt_matrix = adata.layers['alt']

# Convert to lil_matrix for efficient row/column assignment
hap1 = sparse.lil_matrix(ref_matrix.shape, dtype=ref_matrix.dtype)
hap2 = sparse.lil_matrix(ref_matrix.shape, dtype=ref_matrix.dtype)

n_phased = 0
is_phased = []  # Track per-SNP phasing status
for i, snp_id in enumerate(adata.var_names):
    if snp_id in phasing:
        hap1_allele, hap2_allele = phasing[snp_id]
        n_phased += 1
        is_phased.append(True)
    else:
        # Unphased SNPs: arbitrary default of ref=hap1, alt=hap2
        hap1_allele, hap2_allele = 0, 1
        is_phased.append(False)
    # Assign counts: 0=ref, 1=alt for each haplotype
    hap1[:, i] = ref_matrix[:, i] if hap1_allele == 0 else alt_matrix[:, i]
    hap2[:, i] = ref_matrix[:, i] if hap2_allele == 0 else alt_matrix[:, i]

# Convert to csr for storage efficiency
adata.layers['hap1'] = sparse.csr_matrix(hap1)
adata.layers['hap2'] = sparse.csr_matrix(hap2)

# Track phasing status per SNP
adata.var['is_phased'] = is_phased

# Update uns metadata
n_unphased = adata.n_vars - n_phased
adata.uns['phased_snps'] = n_phased
adata.uns['unphased_snps'] = n_unphased
adata.uns['total_snps'] = adata.n_vars
adata.uns['phasing_rate'] = n_phased / adata.n_vars if adata.n_vars > 0 else 0
adata.uns['unphased_default'] = 'ref=hap1, alt=hap2'  # Document the default behavior
adata.uns['pipeline'] = 'nf-scatac'
adata.uns['data_type'] = 'scATAC_allelic_counts_phased'

# Update X to be total (ref + alt)
if sparse.issparse(ref_matrix):
    adata.X = ref_matrix + alt_matrix
else:
    adata.X = sparse.csr_matrix(ref_matrix + alt_matrix)

# Write output
adata.write_h5ad("${prefix}_with_haplotypes.h5ad")

# Generate cell QC metrics
def to_array(mat):
    """Convert sparse matrix sum result to 1D array."""
    return mat.A1 if sparse.issparse(mat) else np.asarray(mat).flatten()

cell_qc = pd.DataFrame({
    'barcode': adata.obs_names,
    'n_snps': to_array((adata.X > 0).sum(axis=1)),
    'total_counts': to_array(adata.X.sum(axis=1)),
    'ref_counts': to_array(adata.layers['ref'].sum(axis=1)),
    'alt_counts': to_array(adata.layers['alt'].sum(axis=1)),
    'hap1_counts': to_array(adata.layers['hap1'].sum(axis=1)),
    'hap2_counts': to_array(adata.layers['hap2'].sum(axis=1)),
})
cell_qc.to_csv("${prefix}_cell_qc.tsv", sep='\\t', index=False)

if ${zarr_flag}:
    adata.write_zarr("${prefix}_with_haplotypes.zarr")

print(f"Created AnnData with haplotype layers: {adata.n_obs} cells x {adata.n_vars} SNPs")
if adata.n_vars > 0:
    phasing_pct = 100 * n_phased / adata.n_vars
    print(f"Phased {n_phased}/{adata.n_vars} SNPs ({phasing_pct:.1f}%)")
    if n_unphased > 0:
        print(f"Note: {n_unphased} SNPs were unphased (defaulted to hap1=ref, hap2=alt)")
    if n_phased == 0:
        print("=" * 70, file=sys.stderr)
        print("WARNING: No SNPs were phased!", file=sys.stderr)
        print("Hap1/hap2 layers contain ARBITRARY assignments (hap1=ref, hap2=alt).", file=sys.stderr)
        print("Check that your VCF contains phased genotypes (with '|' delimiter).", file=sys.stderr)
        print("See adata.var['is_phased'] and adata.uns['unphased_default'] for details.", file=sys.stderr)
        print("=" * 70, file=sys.stderr)
else:
    print("WARNING: No SNPs in input AnnData.", file=sys.stderr)
print(f"Layers: {list(adata.layers.keys())}")
PYEOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        anndata: \$(python -c "import anndata; print(anndata.__version__)")
        scipy: \$(python -c "import scipy; print(scipy.__version__)")
        pysam: \$(python -c "import pysam; print(pysam.__version__)" 2>/dev/null || echo "not installed")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 << 'PYEOF'
import numpy as np
import pandas as pd
from scipy import sparse
import anndata as ad

np.random.seed(42)  # Deterministic stub for reproducible snapshots

# Create stub AnnData with all layers
n_cells, n_snps = 10, 50
X = sparse.random(n_cells, n_snps, density=0.3, format='csr', random_state=42)
ref = sparse.random(n_cells, n_snps, density=0.3, format='csr', random_state=43)
alt = sparse.random(n_cells, n_snps, density=0.3, format='csr', random_state=44)
hap1 = sparse.random(n_cells, n_snps, density=0.3, format='csr', random_state=45)
hap2 = sparse.random(n_cells, n_snps, density=0.3, format='csr', random_state=46)

obs = pd.DataFrame({
    'n_snps': np.random.randint(10, 50, n_cells),
    'total_counts': np.random.randint(100, 1000, n_cells)
}, index=[f'AAACGAAC-{i}' for i in range(n_cells)])

var = pd.DataFrame({
    'chrom': ['chr1'] * n_snps,
    'pos': range(100000, 100000 + n_snps * 1000, 1000),
    'ref': ['A'] * n_snps,
    'alt': ['G'] * n_snps
}, index=[f'chr1:{100000 + i*1000}:A>G' for i in range(n_snps)])

adata = ad.AnnData(X=X, obs=obs, var=var)
adata.layers['ref'] = ref
adata.layers['alt'] = alt
adata.layers['hap1'] = hap1
adata.layers['hap2'] = hap2
adata.uns['phased_snps'] = 40
adata.uns['total_snps'] = n_snps
adata.uns['phasing_rate'] = 0.8
adata.write_h5ad("${prefix}_with_haplotypes.h5ad")

pd.DataFrame({
    'barcode': obs.index,
    'n_snps': obs['n_snps'],
    'total_counts': obs['total_counts'],
    'ref_counts': np.random.randint(50, 500, n_cells),
    'alt_counts': np.random.randint(50, 500, n_cells),
    'hap1_counts': np.random.randint(50, 500, n_cells),
    'hap2_counts': np.random.randint(50, 500, n_cells)
}).to_csv("${prefix}_cell_qc.tsv", sep='\\t', index=False)
PYEOF

    if [ "${create_zarr}" == "true" ]; then
        mkdir -p ${prefix}_with_haplotypes.zarr
        touch ${prefix}_with_haplotypes.zarr/.zgroup
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
        anndata: stub
        scipy: stub
        pysam: stub
    END_VERSIONS
    """
}
