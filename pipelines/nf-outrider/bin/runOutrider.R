#!/usr/bin/env Rscript

#' Run OUTRIDER for Aberrant Expression Detection
#'
#' This script runs the OUTRIDER autoencoder for detecting
#' aberrant expression patterns in RNA-seq data.
#'
#' Based on nf-core/drop implementation with WASP2 enhancements.
#'
#' @param counts Path to count matrix TSV (genes x samples)
#' @param output_model Output path for trained OUTRIDER model (.rds)
#' @param output_results Output path for results table (.tsv)
#' @param padj Adjusted p-value cutoff for outlier calling
#' @param zscore Z-score cutoff for outlier calling
#' @param q Encoding dimension (NULL for auto-estimation)
#' @param iterations Maximum fitting iterations
#' @param convergence Convergence threshold
#' @param threads Number of parallel threads

suppressPackageStartupMessages({
    library(OUTRIDER)
    library(BiocParallel)
    library(data.table)
    library(argparse)
})

# Parse command line arguments
parser <- ArgumentParser(description = "Run OUTRIDER for aberrant expression detection")
parser$add_argument("--counts", required = TRUE, help = "Path to count matrix TSV (genes x samples)")
parser$add_argument("--output_model", default = "outrider_model.rds", help = "Output RDS file for trained model")
parser$add_argument("--output_results", default = "outrider_results.tsv", help = "Output TSV file for results")
parser$add_argument("--padj", type = "double", default = 0.05, help = "Adjusted p-value cutoff for outliers")
parser$add_argument("--zscore", type = "double", default = 2.0, help = "Z-score cutoff for outliers")
parser$add_argument("--q", type = "integer", default = NULL, help = "Encoding dimension (auto-estimate if not provided)")
parser$add_argument("--iterations", type = "integer", default = 15, help = "Maximum iterations")
parser$add_argument("--convergence", type = "double", default = 1e-5, help = "Convergence threshold")
parser$add_argument("--threads", type = "integer", default = 1, help = "Number of threads")

args <- parser$parse_args()

message("=== OUTRIDER Aberrant Expression Detection ===")
message(sprintf("Input: %s", args$counts))
message(sprintf("P-adj cutoff: %s", args$padj))
message(sprintf("Z-score cutoff: %s", args$zscore))

# Set parallel processing
if (args$threads > 1) {
    register(MulticoreParam(args$threads))
} else {
    register(SerialParam())
}

# Load count matrix with error handling
message("Loading count matrix...")
tryCatch({
    if (!file.exists(args$counts)) {
        stop(sprintf("Input file not found: %s", args$counts))
    }
    counts <- fread(args$counts, data.table = FALSE)
    if (nrow(counts) == 0) {
        stop(sprintf("Input file is empty: %s", args$counts))
    }
    if (ncol(counts) < 2) {
        stop(sprintf("Input file must have at least 2 columns (gene_id + samples): %s", args$counts))
    }
}, error = function(e) {
    message(sprintf("ERROR: Failed to load count matrix: %s", e$message))
    message("Expected format: TSV with gene_id in first column, samples in remaining columns")
    quit(status = 1)
})

rownames(counts) <- counts[[1]]
counts <- counts[, -1, drop = FALSE]

message(sprintf("Loaded: %d genes x %d samples", nrow(counts), ncol(counts)))

# Filter low-expressed genes
min_samples <- max(2, floor(ncol(counts) * 0.5))  # At least 50% of samples
row_sums <- rowSums(counts >= 10)
keep_genes <- row_sums >= min_samples
counts_filtered <- counts[keep_genes, , drop = FALSE]

message(sprintf("After filtering: %d genes", nrow(counts_filtered)))

# Validate sufficient genes for analysis
if (nrow(counts_filtered) < 10) {
    message(sprintf("ERROR: Too few genes (%d) for OUTRIDER analysis. Minimum 10 required.", nrow(counts_filtered)))
    quit(status = 1)
} else if (nrow(counts_filtered) < 100) {
    message(sprintf("WARNING: Only %d genes passed filtering. Results may be unreliable.", nrow(counts_filtered)))
    message("Consider relaxing the filter or using a larger sample size.")
}

# Create OutriderDataSet
message("Creating OutriderDataSet...")
ods <- OutriderDataSet(countData = as.matrix(counts_filtered))

# Filter samples with low coverage
ods <- filterExpression(ods, minCounts = TRUE, filterGenes = FALSE)

# Estimate size factors
message("Estimating size factors...")
ods <- estimateSizeFactors(ods)

# Determine encoding dimension (q)
if (is.null(args$q)) {
    message("Auto-estimating encoding dimension...")
    ods <- findEncodingDim(ods)
    q_val <- getBestQ(ods)
    message(sprintf("Optimal q: %d", q_val))
} else {
    q_val <- args$q
    message(sprintf("Using provided q: %d", q_val))
}

# Fit OUTRIDER model
message("Fitting OUTRIDER autoencoder...")
tryCatch({
    ods <- OUTRIDER(
        ods,
        q = q_val,
        controlData = TRUE,
        maxIterations = args$iterations,
        convergence = args$convergence
    )
}, error = function(e) {
    message(sprintf("ERROR: OUTRIDER model fitting failed: %s", e$message))
    if (grepl("singular|convergence", e$message, ignore.case = TRUE)) {
        message("This often indicates too few samples or highly correlated expression patterns.")
        message("Suggestions: (1) increase sample size, (2) reduce encoding dimension (--q)")
    }
    quit(status = 1, save = "no")
})

# Save trained model
message(sprintf("Saving model to: %s", args$output_model))
saveRDS(ods, args$output_model)

# Extract results
message("Extracting outlier results...")
tryCatch({
    res <- results(
        ods,
        padjCutoff = args$padj,
        zScoreCutoff = args$zscore,
        all = TRUE
    )
}, error = function(e) {
    message(sprintf("ERROR: Failed to extract OUTRIDER results: %s", e$message))
    quit(status = 1, save = "no")
})

# Convert to data.table and clean up
res_dt <- as.data.table(res)

# Add aberrant flag
res_dt[, aberrant := (padjust < args$padj) & (abs(zScore) > args$zscore)]

# Save results
message(sprintf("Saving results to: %s", args$output_results))
fwrite(res_dt, args$output_results, sep = "\t")

# Summary statistics
n_outliers <- sum(res_dt$aberrant, na.rm = TRUE)
n_genes <- uniqueN(res_dt$geneID)
n_samples <- uniqueN(res_dt$sampleID)

message("=== Summary ===")
message(sprintf("Total gene-sample pairs tested: %d", nrow(res_dt)))
message(sprintf("Aberrant expression events: %d", n_outliers))
message(sprintf("Genes with aberrant expression: %d / %d",
                uniqueN(res_dt[aberrant == TRUE]$geneID), n_genes))
message(sprintf("Samples with aberrant expression: %d / %d",
                uniqueN(res_dt[aberrant == TRUE]$sampleID), n_samples))

message("=== OUTRIDER Complete ===")
