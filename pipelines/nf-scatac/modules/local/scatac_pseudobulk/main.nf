/*
 * SCATAC_PSEUDOBULK - Aggregate per-cell counts to pseudo-bulk for statistical analysis
 *
 * Aggregates per-cell allele counts to sample-level pseudo-bulk counts.
 * This increases statistical power for allelic imbalance testing by combining
 * sparse per-cell data into denser aggregate counts.
 *
 * Output format has columns matching WASP2_COUNT_ALLELES (ref_count, alt_count).
 * Note: Since fragment files lack sequence data, total overlaps are placed in
 * ref_count while alt_count is zero. This is a limitation of fragment-based
 * scATAC analysis - true allele-specific counting requires BAM-level data.
 */

process SCATAC_PSEUDOBULK {
    tag "$meta.id"
    label 'process_low'

    conda "${projectDir}/../nf-modules/modules/wasp2/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://jaureguy760/wasp2:latest' :
        'jaureguy760/wasp2:latest' }"

    input:
    tuple val(meta), path(cell_counts)

    output:
    tuple val(meta), path("*_pseudobulk_counts.tsv"), emit: counts
    tuple val(meta), path("*_aggregation_stats.tsv"), emit: stats
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def min_cells = params.min_cells_per_snp ?: 3
    """
    set -euo pipefail

    # Aggregate per-cell counts to pseudo-bulk and generate stats in one pass
    awk -v OFS='\\t' -v min_cells="${min_cells}" -v prefix="${prefix}" '
    BEGIN {
        print "chrom", "pos", "ref", "alt", "ref_count", "alt_count" > prefix "_pseudobulk_counts.tsv"
    }
    NR > 1 {
        key = \$2 OFS \$3 OFS \$4 OFS \$5
        total[key] += \$6
        cells_per_snp[key]++
        input_cells[\$1]++
        input_snps[\$2 ":" \$3]++
    }
    END {
        # Write filtered pseudo-bulk counts and count how many pass filter
        filtered_count = 0
        for (key in total) {
            if (cells_per_snp[key] >= min_cells) {
                print key, total[key], 0 >> prefix "_pseudobulk_counts.tsv"
                filtered_count++
            }
        }

        # Write aggregation stats
        print "metric", "value" > prefix "_aggregation_stats.tsv"
        print "total_cells_input", length(input_cells) >> prefix "_aggregation_stats.tsv"
        print "total_snps_input", length(input_snps) >> prefix "_aggregation_stats.tsv"
        print "snps_after_filtering", filtered_count >> prefix "_aggregation_stats.tsv"
    }' ${cell_counts}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -1 | sed 's/GNU Awk //' | cut -d',' -f1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -e "chrom\\tpos\\tref\\talt\\tref_count\\talt_count" > ${prefix}_pseudobulk_counts.tsv
    echo -e "chr1\\t100000\\tA\\tG\\t45\\t0" >> ${prefix}_pseudobulk_counts.tsv
    echo -e "chr1\\t200000\\tC\\tT\\t32\\t0" >> ${prefix}_pseudobulk_counts.tsv

    echo -e "metric\\tvalue" > ${prefix}_aggregation_stats.tsv
    echo -e "total_cells_input\\t10" >> ${prefix}_aggregation_stats.tsv
    echo -e "total_snps_input\\t2" >> ${prefix}_aggregation_stats.tsv
    echo -e "snps_after_filtering\\t2" >> ${prefix}_aggregation_stats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: stub
    END_VERSIONS
    """
}
