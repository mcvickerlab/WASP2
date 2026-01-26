/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: Pipeline utilities for nf-rnaseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; paramsHelp; paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-schema'

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
        log.info "nf-rnaseq v${workflow.manifest.version}"
        System.exit(0)
    }

    //
    // Print help message
    //
    if (help) {
        def help_string = paramsHelp("nextflow run nf-rnaseq --input samplesheet.csv --vcf variants.vcf.gz --star_index /path/to/star -profile docker")
        log.info help_string
        System.exit(0)
    }

    //
    // Validate parameters
    //
    validateParameters()

    //
    // Print parameter summary
    //
    log.info paramsSummaryLog(workflow)

    //
    // Parse and validate input samplesheet
    // Uses toList() to collect all rows before processing to enable thread-safe duplicate detection
    //
    ch_samplesheet = Channel
        .fromPath(input_file, checkIfExists: true)
        .splitCsv(header: true)
        .toList()
        .flatMap { rows ->
            // Check for empty samplesheet
            if (rows.size() == 0) {
                exit 1, "ERROR: Samplesheet is empty (no data rows found)"
            }

            // Validate required columns exist (check first row)
            def first_row = rows[0]
            if (!first_row.containsKey('sample'))  { exit 1, "ERROR: Samplesheet missing 'sample' column" }
            if (!first_row.containsKey('fastq_1')) { exit 1, "ERROR: Samplesheet missing 'fastq_1' column" }

            // Thread-safe duplicate detection: collect all sample IDs and find duplicates
            def sample_ids = rows.collect { it.sample?.trim() ?: '' }
            def duplicates = sample_ids.findAll { id -> id && sample_ids.count(id) > 1 }.unique()
            if (duplicates) {
                exit 1, "ERROR: Duplicate sample ID(s) found: ${duplicates.join(', ')}"
            }

            rows  // Emit rows for individual processing
        }
        .map { row ->
            // Empty values check
            if (!row.sample || row.sample.trim() == '') {
                exit 1, "ERROR: Empty sample ID found in samplesheet"
            }
            if (!row.fastq_1 || row.fastq_1.trim() == '') {
                exit 1, "ERROR: Empty fastq_1 path for sample: ${row.sample}"
            }

            // Sample ID validation - no spaces allowed (provides specific error message)
            if (row.sample =~ /\s/) {
                exit 1, "ERROR: Sample ID '${row.sample}' contains spaces"
            }

            // Sample ID validation - alphanumeric, underscore, hyphen only
            if (!(row.sample =~ /^[a-zA-Z0-9_-]+$/)) {
                exit 1, "ERROR: Sample ID '${row.sample}' contains invalid characters (use alphanumeric, underscore, hyphen only)"
            }

            // FASTQ file extension validation - must be gzipped (case-insensitive)
            if (!(row.fastq_1 =~ /(?i)\.(fq|fastq)\.gz$/)) {
                exit 1, "ERROR: fastq_1 for sample '${row.sample}' must be gzipped (.fq.gz or .fastq.gz)"
            }
            if (row.fastq_2 && row.fastq_2.trim() != '' && !(row.fastq_2 =~ /(?i)\.(fq|fastq)\.gz$/)) {
                exit 1, "ERROR: fastq_2 for sample '${row.sample}' must be gzipped (.fq.gz or .fastq.gz)"
            }

            // FASTQ file existence validation
            def fq1 = file(row.fastq_1)
            if (!fq1.exists()) {
                exit 1, "ERROR: fastq_1 file not found for sample '${row.sample}': ${row.fastq_1}"
            }

            def fq2 = null
            if (row.fastq_2 && row.fastq_2.trim() != '') {
                fq2 = file(row.fastq_2)
                if (!fq2.exists()) {
                    exit 1, "ERROR: fastq_2 file not found for sample '${row.sample}': ${row.fastq_2}"
                }
            }

            def meta = [id: row.sample, single_end: fq2 == null, sample: row.sample]
            def fastqs = fq2 ? [fq1, fq2] : [fq1]
            tuple(meta, fastqs)
        }

    //
    // Prepare VCF channel with index
    //
    ch_vcf = Channel.empty()
    if (params.vcf) {
        ch_vcf = Channel.fromPath(params.vcf)
            .map { vcf ->
                def vcf_index = file("${vcf}.tbi")
                if (!vcf_index.exists()) {
                    vcf_index = file("${vcf}.csi")
                }
                if (!vcf_index.exists()) {
                    exit 1, "ERROR: VCF index not found. Expected ${vcf}.tbi or ${vcf}.csi"
                }
                tuple([id: 'reference'], vcf, vcf_index)
            }
    }

    emit:
    samplesheet = ch_samplesheet  // channel: [ val(meta), [ fastq ] ]
    vcf         = ch_vcf          // channel: [ val(meta), path(vcf), path(vcf_index) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: PIPELINE_COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {
    take:
    outdir    // string: Output directory
    versions  // channel: versions channel

    main:
    //
    // Collect and deduplicate versions
    //
    versions
        .unique()
        .collectFile(name: 'software_versions.yml', storeDir: "${outdir}/pipeline_info")

    //
    // Completion summary
    //
    workflow.onComplete {
        if (workflow.success) {
            log.info "-" * 60
            log.info " WASP2 RNA-seq ASE Pipeline Complete!"
            log.info "-" * 60
            log.info " Output directory: ${outdir}"
            log.info " Duration: ${workflow.duration}"
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
