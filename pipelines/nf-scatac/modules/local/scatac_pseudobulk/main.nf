/*
 * SCATAC_PSEUDOBULK - Aggregate per-cell counts to pseudo-bulk for statistical analysis
 *
 * Aggregates per-cell allele counts to sample-level pseudo-bulk counts.
 * This increases statistical power for allelic imbalance testing by combining
 * sparse per-cell data into denser aggregate counts.
 *
 * Output format matches WASP2_COUNT_ALLELES for compatibility with WASP2_ANALYZE_IMBALANCE.
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
    # Aggregate per-cell counts to pseudo-bulk (sample-level)
    # Input: barcode, chrom, pos, ref, alt, overlap_count
    # Output: chrom, pos, ref, alt, total_count (WASP2-compatible format)

    awk -v OFS='\\t' -v min_cells="${min_cells}" '
    BEGIN {
        # WASP2-compatible header for count_alleles output
        print "chrom", "pos", "ref", "alt", "ref_count", "alt_count"
    }
    NR > 1 {
        # Create SNP key
        key = \$2 "\\t" \$3 "\\t" \$4 "\\t" \$5

        # Aggregate counts across all cells
        # For fragments (no allele info), we store total as ref_count
        # TODO: When phased VCF available, split into hap1/hap2
        total_counts[key] += \$6
        cell_counts[key]++
    }
    END {
        for (key in total_counts) {
            # Filter SNPs with insufficient cell coverage
            if (cell_counts[key] >= min_cells) {
                # Output total as ref_count, 0 as alt_count (no allele resolution in fragments)
                print key, total_counts[key], 0
            }
        }
    }' ${cell_counts} > ${prefix}_pseudobulk_counts.tsv

    # Generate aggregation statistics
    awk -v OFS='\\t' '
    BEGIN {
        print "metric", "value"
    }
    NR > 1 {
        snps++
        total_frags += \$5
    }
    END {
        print "snps_after_filtering", snps
        print "total_fragment_overlaps", total_frags
    }' ${prefix}_pseudobulk_counts.tsv > ${prefix}_aggregation_stats.tsv

    # Also count cells and SNPs from input
    awk -v OFS='\\t' '
    NR > 1 {
        cells[\$1]++
        snps[\$2 ":" \$3]++
    }
    END {
        print "total_cells_input", length(cells) >> "${prefix}_aggregation_stats.tsv"
        print "total_snps_input", length(snps) >> "${prefix}_aggregation_stats.tsv"
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
    echo -e "snps_after_filtering\\t2" >> ${prefix}_aggregation_stats.tsv
    echo -e "total_fragment_overlaps\\t77" >> ${prefix}_aggregation_stats.tsv
    echo -e "total_cells_input\\t10" >> ${prefix}_aggregation_stats.tsv
    echo -e "total_snps_input\\t2" >> ${prefix}_aggregation_stats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: stub
    END_VERSIONS
    """
}
