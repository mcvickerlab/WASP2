/*
 * SCATAC_COUNT_ALLELES - Count fragment overlaps at SNP positions
 *
 * Counts fragment overlaps from 10x fragments.tsv.gz at heterozygous SNP positions.
 * Supports optional cell barcode filtering and peak region filtering.
 *
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
    path(barcodes)  // Optional: file with valid barcodes (one per line)
    path(peaks)     // Optional: BED file with peak regions to restrict analysis

    output:
    tuple val(meta), path("*_allele_counts.tsv"), emit: counts
    tuple val(meta), path("*_count_stats.tsv")  , emit: stats
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def min_frags = params.min_fragments_per_cell ?: 1000
    def filter_barcodes = barcodes.name != 'NO_FILE' ? "true" : "false"
    def filter_peaks = peaks.name != 'NO_FILE' ? "true" : "false"
    """
    set -o pipefail

    # Optionally filter SNPs to peak regions
    if [ "${filter_peaks}" == "true" ]; then
        bedtools intersect -a ${snp_bed} -b ${peaks} -u > snps_in_peaks.bed
        SNP_BED="snps_in_peaks.bed"
    else
        SNP_BED="${snp_bed}"
    fi

    # Set barcode file (empty string disables filtering)
    BARCODE_FILE=\$( [ "${filter_barcodes}" == "true" ] && echo "${barcodes}" || echo "" )

    # Count fragment overlaps at SNP positions
    bedtools intersect -a ${fragments} -b \$SNP_BED -wa -wb | \\
    awk -v OFS='\\t' -v bc_file="\$BARCODE_FILE" -v min_frags="${min_frags}" '
    BEGIN {
        if (bc_file != "") { while ((getline bc < bc_file) > 0) valid[bc]=1; close(bc_file) }
    }
    {
        if (bc_file != "" && !(\$4 in valid)) next
        key = \$4 OFS \$6 OFS \$8 OFS \$9 OFS \$10
        counts[key] += \$5; bc_total[\$4] += \$5
    }
    END {
        print "barcode", "chrom", "pos", "ref", "alt", "overlap_count"
        for (k in counts) { split(k,p,OFS); if (bc_total[p[1]] >= min_frags) print k, counts[k] }
    }' > ${prefix}_allele_counts.tsv

    # Generate counting statistics
    awk -v OFS='\\t' 'BEGIN { print "metric", "value" }
    NR > 1 { bc[\$1]++; snp[\$2 ":" \$3]++; tot += \$6 }
    END {
        print "total_barcodes", length(bc)
        print "total_snps", length(snp)
        print "total_fragment_overlaps", tot
        print "mean_snps_per_cell", length(bc) > 0 ? length(snp)/length(bc) : 0
    }' ${prefix}_allele_counts.tsv > ${prefix}_count_stats.tsv

    [ \$(wc -l < ${prefix}_allele_counts.tsv) -lt 2 ] && echo "WARNING: No overlaps found" >&2 || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed 's/bedtools v//')
        awk: \$(awk --version | head -1 | sed 's/GNU Awk //' | cut -d',' -f1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -e "barcode\\tchrom\\tpos\\tref\\talt\\toverlap_count" > ${prefix}_allele_counts.tsv
    echo -e "AAACGAACAGTCAGTT-1\\tchr1\\t100000\\tA\\tG\\t8" >> ${prefix}_allele_counts.tsv
    echo -e "AAACGAACAGTCAGTT-1\\tchr1\\t200000\\tC\\tT\\t5" >> ${prefix}_allele_counts.tsv
    echo -e "AAACGAATCTGCGGCA-1\\tchr1\\t100000\\tA\\tG\\t12" >> ${prefix}_allele_counts.tsv

    echo -e "metric\\tvalue" > ${prefix}_count_stats.tsv
    echo -e "total_barcodes\\t2" >> ${prefix}_count_stats.tsv
    echo -e "total_snps\\t2" >> ${prefix}_count_stats.tsv
    echo -e "total_fragment_overlaps\\t25" >> ${prefix}_count_stats.tsv
    echo -e "mean_snps_per_cell\\t1.0" >> ${prefix}_count_stats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: stub
        awk: stub
    END_VERSIONS
    """
}
