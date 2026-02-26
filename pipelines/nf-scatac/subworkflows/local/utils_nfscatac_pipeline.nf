/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Utility subworkflows for nf-scatac pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {
    take:
    version
    help
    validate_params
    input

    main:
    if (version) {
        log.info "nf-scatac v${workflow.manifest.version}"
        System.exit(0)
    }

    if (help) {
        log.info helpMessage()
        System.exit(0)
    }

    // Parse samplesheet with scATAC-specific meta map fields
    // Supports optional barcodes, peaks, and BAM columns
    // BAM input enables true allele-specific counting with ref/alt/hap layers
    ch_samplesheet = Channel.fromPath(input, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            if (!row.sample) {
                error "Samplesheet error: 'sample' column is required. Found columns: ${row.keySet()}"
            }
            if (!row.fragments && !row.cellranger_dir && !row.bam) {
                error "Samplesheet error for '${row.sample}': provide 'fragments', 'cellranger_dir', or 'bam'"
            }

            def meta = [
                id: row.sample,
                single_end: false,
                cell_barcode_tag: row.barcode_tag ?: 'CB',
                chemistry: row.chemistry ?: '10x-atac-v2',
                has_bam: row.bam as boolean || row.cellranger_dir as boolean
            ]

            // Resolve fragments file path (optional when BAM is provided)
            def fragments = file('NO_FILE_FRAGS')
            def fragments_tbi = file('NO_FILE_FRAGS_TBI')
            if (row.fragments) {
                fragments = file(row.fragments, checkIfExists: true)
                fragments_tbi = file("${fragments}.tbi", checkIfExists: true)
            } else if (row.cellranger_dir) {
                def frag_path = "${row.cellranger_dir}/outs/fragments.tsv.gz"
                if (file(frag_path).exists()) {
                    fragments = file(frag_path, checkIfExists: true)
                    fragments_tbi = file("${frag_path}.tbi", checkIfExists: true)
                }
            }

            // Optional: BAM file for true allele-specific counting
            def bam = file('NO_FILE_BAM')
            def bai = file('NO_FILE_BAI')
            if (row.bam && row.bam.trim()) {
                bam = file(row.bam, checkIfExists: true)
                // Try common BAI naming conventions: .bam.bai and .bai
                def bai_path1 = file("${bam}.bai")
                def bai_path2 = file("${bam}".replaceAll(/\.bam$/, '.bai'))
                if (bai_path1.exists()) {
                    bai = bai_path1
                } else if (bai_path2.exists()) {
                    bai = bai_path2
                } else {
                    error "Samplesheet error for '${row.sample}': BAM index not found. Tried: ${bai_path1}, ${bai_path2}"
                }
            } else if (row.cellranger_dir) {
                def bam_path = "${row.cellranger_dir}/outs/possorted_bam.bam"
                if (file(bam_path).exists()) {
                    bam = file(bam_path, checkIfExists: true)
                    bai = file("${bam_path}.bai", checkIfExists: true)
                }
            }

            // Optional: cell barcode whitelist file
            def barcodes = row.barcodes && row.barcodes.trim()
                ? file(row.barcodes, checkIfExists: true)
                : file('NO_FILE_BARCODES')

            // Optional: peak BED file for restricting analysis to peak regions
            def peaks = row.peaks && row.peaks.trim()
                ? file(row.peaks, checkIfExists: true)
                : file('NO_FILE_PEAKS')

            [ meta, fragments, fragments_tbi, barcodes, peaks, bam, bai ]
        }

    emit:
    samplesheet = ch_samplesheet  // channel: [ val(meta), path(fragments), path(fragments_tbi), path(barcodes), path(peaks), path(bam), path(bai) ]
}

workflow PIPELINE_COMPLETION {
    take:
    outdir
    multiqc_report

    main:
    log.info """
    =============================================================
     nf-scatac COMPLETE
    =============================================================
    Output: ${outdir}

    Results:
      - allele_counts/: Per-cell allele counts at het SNPs
      - anndata/: AnnData H5AD files for scverse ecosystem
      - cell_qc/: Cell QC metrics
      - imbalance/: Allelic imbalance analysis results (pseudo-bulk)
    """.stripIndent()
}

def helpMessage() {
    return """
    nf-scatac - Single-Cell ATAC-seq Allelic Imbalance Pipeline

    Usage:
        nextflow run nf-scatac -profile docker --input samplesheet.csv --vcf variants.vcf.gz

    Required:
        --input       Samplesheet CSV (see format below)
        --vcf         Indexed VCF/BCF with heterozygous SNPs (phased recommended)

    Optional:
        --outdir                  Output directory [default: ./results]
        --min_fragments_per_cell  Minimum fragments per cell to include [default: 1000]
        --min_cells_per_snp       Minimum cells per SNP for pseudo-bulk [default: 3]
        --create_zarr             Also output Zarr format for GenVarLoader [default: false]
        --skip_pseudobulk         Skip pseudo-bulk aggregation [default: false]
        --skip_anndata            Skip AnnData H5AD creation [default: false]

    Samplesheet format:
        sample,fragments,cellranger_dir,bam,barcode_tag,chemistry,barcodes,peaks
        sample1,/path/to/fragments.tsv.gz,,,CB,10x-atac-v2,/path/to/barcodes.txt,/path/to/peaks.bed
        sample2,,/path/to/cellranger/output,,CB,10x-atac-v2,,
        sample3,,,/path/to/possorted.bam,CB,10x-atac-v2,/path/to/barcodes.txt,

    Column descriptions:
        sample         - Sample identifier (required)
        fragments      - Path to fragments.tsv.gz
        cellranger_dir - Path to CellRanger ATAC output (auto-detects fragments & BAM)
        bam            - Path to indexed BAM (enables true allele-specific counting)
        barcode_tag    - BAM tag for cell barcodes [default: CB]
        chemistry      - Library chemistry [default: 10x-atac-v2]
        barcodes       - Optional: file with valid cell barcodes (one per line)
        peaks          - Optional: BED file with peak regions to restrict analysis

    Input priority: At least one of fragments, cellranger_dir, or bam is required.
    When BAM is provided, true allele-specific counting is performed with ref/alt/hap layers.
    When only fragments are provided, overlap counting is performed (total counts only).

    Profiles:
        -profile docker        Run with Docker
        -profile singularity   Run with Singularity
        -profile conda         Run with Conda
        -profile test_stub     Run stub tests
        -profile test_real     Run integration tests with real data

    Output:
        results/
        ├── allele_counts/      # Per-cell allele counts
        ├── anndata/            # AnnData H5AD files (with ref/alt/hap layers if BAM provided)
        ├── zarr/               # Zarr files (if --create_zarr)
        ├── cell_qc/            # Cell QC metrics
        ├── pseudobulk/         # Pseudo-bulk aggregated counts
        ├── imbalance/          # Allelic imbalance results
        └── pipeline_info/      # Execution reports
    """.stripIndent()
}
