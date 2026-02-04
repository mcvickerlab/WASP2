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
    # Extract heterozygous SNPs to BED format
    python3 << 'EOF'
from wasp2.io.variant_source import VariantSource

source = VariantSource.open("${vcf}")
samples_list = "${samples}".split(",") if "${samples}" else None

with open("${prefix}.variants.bed", "w") as f:
    for var in source.iter_variants(samples=samples_list, het_only=True):
        # BED format: chrom, start (0-based), end (1-based), ref, alt
        f.write(f"{var.chrom}\\t{var.pos - 1}\\t{var.pos}\\t{var.ref}\\t{var.alt}\\n")
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
