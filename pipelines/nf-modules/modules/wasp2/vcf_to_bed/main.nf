process WASP2_VCF_TO_BED {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/../environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wasp2:1.2.1--pyhdfd78af_0' :
        'biocontainers/wasp2:1.2.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    val(samples)

    output:
    tuple val(meta), path("*.variants.bed"), emit: bed
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sample_arg = samples ? "--samples ${samples}" : ''
    """
    # Extract heterozygous SNPs to BED format using WASP2 v1.4.0 API
    python3 << 'EOF'
from wasp2.io.variant_source import VariantSource
from pathlib import Path

source = VariantSource.open("${vcf}")
samples_list = "${samples}".split(",") if "${samples}" else None

# Use to_bed() method (v1.4.0+) — iter_variants() returns VariantGenotype
# objects without chrom/pos attributes in v1.4.0
out_path = source.to_bed(
    Path("${prefix}.variants.bed"),
    samples=samples_list,
    het_only=True,
    include_genotypes=False
)

# Deduplicate BED entries (same SNP het in multiple samples)
with open(out_path) as f:
    lines = sorted(set(f.readlines()))
with open(out_path, 'w') as f:
    f.writelines(lines)
EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wasp2: \$(python -c "import wasp2; print(wasp2.__version__)" 2>/dev/null || echo "dev")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.variants.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wasp2: \$(python -c "import wasp2; print(wasp2.__version__)" 2>/dev/null || echo "dev")
    END_VERSIONS
    """
}
