/*
 * SCATAC_COUNT_ALLELES - Count fragment overlaps at SNP positions
 *
 * Counts fragment overlaps from 10x fragments.tsv.gz at heterozygous SNP positions.
 * Note: Fragment files contain only coordinates, not sequences, so we count total
 * overlaps per barcode/SNP. The overlap_count represents coverage at SNP sites;
 * allele-specific counts require downstream phasing or sequence-level analysis.
 */

process SCATAC_COUNT_ALLELES {
    tag "$meta.id"
    label 'process_medium'

    conda "${projectDir}/../nf-modules/modules/wasp2/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://jaureguy760/wasp2:latest' :
        'jaureguy760/wasp2:latest' }"

    input:
    tuple val(meta), path(fragments), path(fragments_tbi), path(snp_bed)

    output:
    tuple val(meta), path("*_allele_counts.tsv"), emit: counts
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    set -o pipefail

    bedtools intersect -a ${fragments} -b ${snp_bed} -wa -wb | \\
    awk -v OFS='\\t' 'BEGIN {
        print "barcode", "chrom", "pos", "ref", "alt", "overlap_count"
    }
    {
        # Fragment: chrom,start,end,barcode,read_support | SNP BED: chrom,start,end,ref,alt
        key = \$4 "\\t" \$6 "\\t" \$8 "\\t" \$9 "\\t" \$10
        counts[key] += \$5
    }
    END {
        for (key in counts) {
            print key, counts[key]
        }
    }' > ${prefix}_allele_counts.tsv

    # Validate output file was created and has header
    if [[ ! -s ${prefix}_allele_counts.tsv ]]; then
        echo "ERROR: Allele counts file is empty. Check bedtools intersect output." >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed 's/bedtools v//')
        wasp2: \$(python -c "import wasp2; print(wasp2.__version__)" || echo "unavailable")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -e "barcode\\tchrom\\tpos\\tref\\talt\\toverlap_count" > ${prefix}_allele_counts.tsv
    echo -e "AAACGAACAGTCAGTT-1\\tchr1\\t100000\\tA\\tG\\t8" >> ${prefix}_allele_counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: stub
        wasp2: stub
    END_VERSIONS
    """
}
