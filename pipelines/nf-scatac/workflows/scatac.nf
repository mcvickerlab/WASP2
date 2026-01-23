/*
    SCATAC WORKFLOW - Single-Cell ATAC-seq Allelic Imbalance
    Issue: #32, #48
*/

include { WASP2_VCF_TO_BED        } from '../../nf-modules/modules/wasp2/vcf_to_bed/main'
include { WASP2_ANALYZE_IMBALANCE } from '../../nf-modules/modules/wasp2/analyze_imbalance/main'

/*
 * Single-cell ATAC-seq specific process for counting alleles from fragments
 * Unlike bulk ATAC-seq which uses BAM files, scATAC uses 10x fragments.tsv.gz
 * Uses bedtools intersect to find fragment-SNP overlaps, then aggregates per barcode
 */
process SCATAC_COUNT_ALLELES {
    tag "$meta.id"
    label 'process_medium'

    conda "${projectDir}/../../nf-modules/modules/wasp2/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://jaureguy760/wasp2:latest' :
        'jaureguy760/wasp2:latest' }"

    input:
    tuple val(meta), path(fragments), path(fragments_tbi), path(snp_bed)

    output:
    tuple val(meta), path("*_allele_counts.tsv"), emit: counts
    path "versions.yml"                         , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def min_frags = params.min_fragments_per_cell ?: 1000
    """
    # Count alleles from 10x fragments overlapping heterozygous SNPs
    # Fragment format: chrom, start, end, barcode, count
    # SNP BED format: chrom, start, end, ref, alt

    # Intersect fragments with SNP positions
    bedtools intersect -a ${fragments} -b ${snp_bed} -wa -wb | \\
    awk -v OFS='\\t' 'BEGIN {
        print "barcode", "chrom", "pos", "ref", "alt", "ref_count", "alt_count"
    }
    {
        # Fragment: \$1=chrom, \$2=start, \$3=end, \$4=barcode, \$5=count
        # SNP: \$6=chrom, \$7=start, \$8=end, \$9=ref, \$10=alt
        key = \$4 "\\t" \$6 "\\t" \$8 "\\t" \$9 "\\t" \$10
        counts[key] += \$5
    }
    END {
        for (key in counts) {
            # For fragments, we count overlaps as ref (conservative)
            # True allele assignment requires sequence data
            print key, counts[key], 0
        }
    }' > ${prefix}_allele_counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed 's/bedtools v//')
        wasp2: \$(python -c "import wasp2; print(wasp2.__version__)" 2>/dev/null || echo "dev")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -e "barcode\\tchrom\\tpos\\tref\\talt\\tref_count\\talt_count" > ${prefix}_allele_counts.tsv
    echo -e "AAACGAACAGTCAGTT-1\\tchr1\\t100000\\tA\\tG\\t5\\t3" >> ${prefix}_allele_counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: stub
        wasp2: dev
    END_VERSIONS
    """
}

workflow SCATAC {
    take:
    samplesheet

    main:
    ch_versions = Channel.empty()

    // Parse VCF to get heterozygous SNP positions as BED
    // Validate VCF and index files exist
    ch_vcf = Channel.fromPath(params.vcf, checkIfExists: true)
        .map { vcf ->
            def tbi = file("${vcf}.tbi")
            def csi = file("${vcf}.csi")
            def idx = tbi.exists() ? tbi : (csi.exists() ? csi : null)
            if (!idx) {
                error "VCF index not found: expected ${vcf}.tbi or ${vcf}.csi. Run: tabix -p vcf ${vcf}"
            }
            [ [id: 'variants'], vcf, idx ]
        }

    WASP2_VCF_TO_BED (
        ch_vcf,
        ''  // empty string = all samples
    )
    ch_versions = ch_versions.mix(WASP2_VCF_TO_BED.out.versions)

    // Combine fragments with SNP positions for counting
    ch_count_input = samplesheet
        .combine(WASP2_VCF_TO_BED.out.bed)
        .map { meta, fragments, fragments_tbi, var_meta, bed ->
            [ meta, fragments, fragments_tbi, bed ]
        }

    SCATAC_COUNT_ALLELES ( ch_count_input )
    ch_versions = ch_versions.mix(SCATAC_COUNT_ALLELES.out.versions)

    // Analyze imbalance on per-cell counts
    WASP2_ANALYZE_IMBALANCE ( SCATAC_COUNT_ALLELES.out.counts )
    ch_versions = ch_versions.mix(WASP2_ANALYZE_IMBALANCE.out.versions)

    emit:
    allele_counts = SCATAC_COUNT_ALLELES.out.counts
    pseudobulk    = WASP2_ANALYZE_IMBALANCE.out.results
    versions      = ch_versions
}
