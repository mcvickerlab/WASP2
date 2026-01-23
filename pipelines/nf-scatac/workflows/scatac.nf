/*
    SCATAC WORKFLOW - Single-Cell ATAC-seq Allelic Imbalance
    Issue: #32
*/

include { WASP2_COUNT_ALLELES     } from '../../nf-modules/modules/wasp2_count_alleles/main'
include { WASP2_ANALYZE_IMBALANCE } from '../../nf-modules/modules/wasp2_analyze_imbalance/main'
include { VCF_TO_BED              } from '../../nf-modules/modules/vcf_to_bed/main'

workflow SCATAC {
    take:
    samplesheet

    main:
    ch_versions = Channel.empty()

    ch_vcf = Channel.fromPath(params.vcf)
        .map { vcf -> [ [id: 'variants'], vcf, file("${vcf}.tbi") ] }

    VCF_TO_BED ( ch_vcf )
    ch_versions = ch_versions.mix(VCF_TO_BED.out.versions)

    ch_count_input = samplesheet
        .combine(VCF_TO_BED.out.bed)
        .map { meta, fragments, fragments_tbi, var_meta, bed ->
            [ meta, fragments, fragments_tbi, bed ]
        }

    WASP2_COUNT_ALLELES ( ch_count_input )
    ch_versions = ch_versions.mix(WASP2_COUNT_ALLELES.out.versions)

    WASP2_ANALYZE_IMBALANCE ( WASP2_COUNT_ALLELES.out.counts )
    ch_versions = ch_versions.mix(WASP2_ANALYZE_IMBALANCE.out.versions)

    emit:
    allele_counts  = WASP2_COUNT_ALLELES.out.counts
    pseudobulk     = WASP2_ANALYZE_IMBALANCE.out.results
    cell_clusters  = Channel.empty()
    versions       = ch_versions
    multiqc_report = Channel.empty()
}
