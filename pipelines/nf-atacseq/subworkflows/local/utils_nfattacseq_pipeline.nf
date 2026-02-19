/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: Pipeline utilities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; paramsHelp; paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: PIPELINE_INITIALISATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {
    take:
    version      // boolean: Display version and exit
    help         // boolean: Display help text
    input_file   // string: Path to input samplesheet

    main:
    //
    // Print version and exit
    //
    if (version) {
        log.info "nf-atacseq v${workflow.manifest.version}"
        System.exit(0)
    }

    //
    // Print help message
    //
    if (help) {
        def help_string = paramsHelp("nextflow run nf-atacseq --input samplesheet.csv --vcf variants.vcf.gz --fasta genome.fa -profile docker")
        log.info help_string
        System.exit(0)
    }

    //
    // Validate parameters
    //
    // validateParameters() // Skipped: URLs fail file-exists check

    //
    // Print parameter summary
    //
    log.info paramsSummaryLog(workflow)

    //
    // Parse samplesheet
    //
    ch_samplesheet = Channel.fromPath(input_file, checkIfExists: true)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            // Create meta map
            def meta = [:]
            meta.id          = row.sample
            meta.single_end  = row.single_end ? row.single_end.toBoolean() : false
            meta.sample_name = row.sample_name ?: null

            // Check FASTQ files exist
            def fastq_1 = file(row.fastq_1, checkIfExists: true)
            def fastq_2 = row.fastq_2 ? file(row.fastq_2, checkIfExists: true) : null

            // Return tuple
            if (meta.single_end) {
                return [ meta, [ fastq_1 ] ]
            } else {
                if (!fastq_2) {
                    error "ERROR: Paired-end data requires fastq_2 for sample '${meta.id}'"
                }
                return [ meta, [ fastq_1, fastq_2 ] ]
            }
        }

    emit:
    samplesheet = ch_samplesheet  // channel: [ val(meta), [ fastq ] ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: PIPELINE_COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {
    take:
    outdir         // string: Output directory
    multiqc_report // channel: MultiQC report file

    main:
    //
    // Completion summary
    //
    workflow.onComplete {
        if (workflow.success) {
            log.info "-" * 60
            log.info "Pipeline completed successfully!"
            log.info "-" * 60
            log.info "Output directory: ${outdir}"
            log.info "Duration: ${workflow.duration}"
            log.info "-" * 60
        } else {
            log.error "-" * 60
            log.error "Pipeline completed with errors"
            log.error "-" * 60
            log.error "Check '.nextflow.log' for details"
            log.error "-" * 60
        }
    }

    workflow.onError {
        log.error "Pipeline execution stopped with the following error: ${workflow.errorMessage}"
    }
}
