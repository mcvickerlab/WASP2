/*
 * Aggregate fragment overlap coverage across cells at each SNP position.
 * Fragment coordinates do not identify observed alleles. Outputs are coverage
 * summaries and must not be submitted to allelic imbalance inference.
 */

process SCATAC_PSEUDOBULK {
    tag "$meta.id"
    label 'process_low'

    conda "${projectDir}/../nf-modules/modules/wasp2/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://ghcr.io/mcvickerlab/wasp2:1.5.0' :
        'ghcr.io/mcvickerlab/wasp2:1.5.0' }"

    input:
    tuple val(meta), path(cell_counts)

    output:
    tuple val(meta), path("*_pseudobulk_coverage.tsv"), emit: counts
    tuple val(meta), path("*_aggregation_stats.tsv"), emit: stats
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def min_cells = params.min_cells_per_snp != null ? params.min_cells_per_snp : 3
    """
    set -euo pipefail

    # Aggregate per-cell counts to pseudo-bulk and generate stats in one pass
    awk -v OFS='\\t' -v min_cells="${min_cells}" -v prefix="${prefix}" '
    BEGIN {
        print "chrom", "pos", "ref", "alt", "overlap_count", "n_cells" > (prefix "_pseudobulk_coverage.tsv")
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
                print key, total[key], cells_per_snp[key] >> (prefix "_pseudobulk_coverage.tsv")
                filtered_count++
            }
        }

        # Write aggregation stats
        print "metric", "value" > (prefix "_aggregation_stats.tsv")
        print "total_cells_input", length(input_cells) >> (prefix "_aggregation_stats.tsv")
        print "total_snps_input", length(input_snps) >> (prefix "_aggregation_stats.tsv")
        print "snps_after_filtering", filtered_count >> (prefix "_aggregation_stats.tsv")
    }' ${cell_counts}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$( (awk --version 2>/dev/null || awk -version) | sed -n '1p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -e "chrom\\tpos\\tref\\talt\\toverlap_count\\tn_cells" > ${prefix}_pseudobulk_coverage.tsv
    echo -e "chr1\\t100000\\tA\\tG\\t45\\t3" >> ${prefix}_pseudobulk_coverage.tsv
    echo -e "chr1\\t200000\\tC\\tT\\t32\\t3" >> ${prefix}_pseudobulk_coverage.tsv

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
