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

    # Prepare SNP regions (optionally intersect with peaks)
    if [ "${filter_peaks}" == "true" ]; then
        # Only count SNPs that fall within peak regions
        bedtools intersect -a ${snp_bed} -b ${peaks} -u > snps_in_peaks.bed
        SNP_BED="snps_in_peaks.bed"
    else
        SNP_BED="${snp_bed}"
    fi

    # Load barcode whitelist if provided
    if [ "${filter_barcodes}" == "true" ]; then
        BARCODE_FILE="${barcodes}"
    else
        BARCODE_FILE=""
    fi

    # Count fragment overlaps at SNP positions with optional barcode filtering
    bedtools intersect -a ${fragments} -b \$SNP_BED -wa -wb | \\
    awk -v OFS='\\t' -v barcode_file="\$BARCODE_FILE" -v min_frags="${min_frags}" '
    BEGIN {
        # Load barcode whitelist if provided
        if (barcode_file != "") {
            while ((getline bc < barcode_file) > 0) {
                valid_bc[bc] = 1
            }
            close(barcode_file)
            filter_bc = 1
        } else {
            filter_bc = 0
        }
    }
    {
        # Fragment: chrom,start,end,barcode,read_support | SNP BED: chrom,start,end,ref,alt
        barcode = \$4
        read_support = \$5
        snp_chrom = \$6
        snp_pos = \$8
        ref = \$9
        alt = \$10

        # Skip if barcode filtering is enabled and barcode not in whitelist
        if (filter_bc && !(barcode in valid_bc)) next

        # Accumulate counts per barcode-SNP pair
        key = barcode "\\t" snp_chrom "\\t" snp_pos "\\t" ref "\\t" alt
        counts[key] += read_support
        barcode_total[barcode] += read_support
    }
    END {
        # Output header
        print "barcode", "chrom", "pos", "ref", "alt", "overlap_count"

        # Output counts, optionally filtering by minimum fragments per cell
        for (key in counts) {
            split(key, parts, "\\t")
            bc = parts[1]
            if (barcode_total[bc] >= min_frags) {
                print key, counts[key]
            }
        }
    }' > ${prefix}_allele_counts.tsv

    # Generate counting statistics
    awk -v OFS='\\t' 'BEGIN {
        print "metric", "value"
    }
    NR > 1 {
        barcodes[\$1]++
        snps[\$2 ":" \$3]++
        total_counts += \$6
    }
    END {
        print "total_barcodes", length(barcodes)
        print "total_snps", length(snps)
        print "total_fragment_overlaps", total_counts
        print "mean_snps_per_cell", (length(barcodes) > 0 ? length(snps)/length(barcodes) : 0)
    }' ${prefix}_allele_counts.tsv > ${prefix}_count_stats.tsv

    # Validate output file was created and has content
    line_count=\$(wc -l < ${prefix}_allele_counts.tsv)
    if [[ \$line_count -lt 2 ]]; then
        echo "WARNING: No overlaps found after filtering. Output contains only header." >&2
        # Don't fail - empty results are valid for some samples
    fi

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
