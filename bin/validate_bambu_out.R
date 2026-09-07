#!/usr/bin/env Rscript

# Regenerate the Bambu PCA and heatmap plots after the transcriptome has been
# validated and curated by the metadata refinement steps. The raw plots produced
# by bambu.R are built from the unfiltered assembly, so they still include novel
# calls that were later discarded. This script subsets the SummarizedExperiment
# objects to the validated transcript and gene sets and re-plots them.

# --- Load necessary libraries ---

library(ggplot2)
library(bambu)
library(SummarizedExperiment)
library(optparse)

# --- Parse command-line arguments ---
option_list <- list(
  make_option(c("-s", "--se_rds"), type = "character", default = NULL,
              help = "Path to the transcript-level SummarizedExperiment RDS (se_multiSample.rds)", metavar = "character"),
  make_option(c("-g", "--se_gene_rds"), type = "character", default = NULL,
              help = "Path to the gene-level SummarizedExperiment RDS (seGene_multiSample.rds)", metavar = "character"),
  make_option(c("-t", "--counts_transcript"), type = "character", default = NULL,
              help = "Path to the validated transcript counts file", metavar = "character"),
  make_option(c("-c", "--counts_gene"), type = "character", default = NULL,
              help = "Path to the validated gene counts file", metavar = "character"),
  make_option(c("-o", "--outdir"), type = "character", default = ".",
              help = "Output directory", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$se_rds) || is.null(opt$se_gene_rds) ||
    is.null(opt$counts_transcript) || is.null(opt$counts_gene)) {
  print_help(opt_parser)
  stop("Missing required arguments.")
}

output_dir <- opt$outdir
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Load the SummarizedExperiment objects ---
se <- tryCatch({
  readRDS(opt$se_rds)
}, error = function(e) {
  stop(paste("Error reading transcript-level RDS:", e$message))
})

seGene <- tryCatch({
  readRDS(opt$se_gene_rds)
}, error = function(e) {
  stop(paste("Error reading gene-level RDS:", e$message))
})

# --- Read the validated ID sets ---
# Column 1 of the validated counts files holds the feature ID, with a header row.
read_validated_ids <- function(path, label) {
  ids <- tryCatch({
    tab <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                      check.names = FALSE, comment.char = "")
    as.character(tab[[1]])
  }, error = function(e) {
    stop(paste("Error reading", label, "counts file:", e$message))
  })

  ids <- unique(ids[!is.na(ids) & nzchar(ids)])
  if (length(ids) == 0) {
    stop(paste("No validated", label, "IDs found in", path))
  }
  ids
}

valid_tx_ids   <- read_validated_ids(opt$counts_transcript, "transcript")
valid_gene_ids <- read_validated_ids(opt$counts_gene, "gene")

# --- Subset the objects to the validated features ---
subset_se <- function(object, ids, label) {
  keep <- rownames(object) %in% ids
  n_keep <- sum(keep)

  if (n_keep == 0) {
    stop(sprintf(paste("None of the %d validated %s IDs matched the %s rownames in the RDS object.",
                       "Check that the counts file and the RDS come from the same Bambu run."),
                 length(ids), label, label))
  }

  cat(sprintf("Retained %d of %d quantified %ss after validation filtering.\n",
              n_keep, nrow(object), label))

  missing <- setdiff(ids, rownames(object))
  if (length(missing) > 0) {
    warning(sprintf("%d validated %s IDs were not present in the RDS object and were ignored.",
                    length(missing), label))
  }

  object[keep, , drop = FALSE]
}

se_filtered     <- subset_se(se, valid_tx_ids, "transcript")
seGene_filtered <- subset_se(seGene, valid_gene_ids, "gene")

# --- Identify the expressed features ----------------------------------------
# Bambu quantifies against the extended annotation, so the raw objects carry a
# row for every reference feature, including those with no reads in any sample.
# Those rows are dropped by the zero-count filter in filter_bambu_counts.sh, so
# measuring the yield against them would count the reference itself as a loss.
# The denominator below is the same set that filter keeps: features with at
# least one count in at least one sample.
expressed_rows <- function(object, label) {
  assay_names <- SummarizedExperiment::assayNames(object)
  if (length(assay_names) == 0) {
    stop(sprintf("The %s-level RDS object carries no assay to test for counts.", label))
  }
  counts_assay <- if ("counts" %in% assay_names) "counts" else assay_names[1]

  totals <- rowSums(SummarizedExperiment::assay(object, counts_assay), na.rm = TRUE)
  keep   <- totals > 0

  cat(sprintf("%d of %d %ss carry at least one %s count in at least one sample.\n",
              sum(keep), length(keep), label, counts_assay))
  keep
}

tx_expressed   <- expressed_rows(se, "transcript")
gene_expressed <- expressed_rows(seGene, "gene")

# --- Record the validation yield as a result --------------------------------
# A headline number for the run, so it is written as data rather than left to be
# scraped back out of the console log. Both the numerator and the denominator
# are restricted to the expressed features.
retained_expressed <- function(object, ids, expressed, label) {
  keep <- rownames(object) %in% ids & expressed
  silent <- sum(rownames(object) %in% ids & !expressed)
  if (silent > 0) {
    warning(sprintf(paste("%d validated %s IDs have no counts in any sample and were",
                          "excluded from the validation summary."), silent, label))
  }
  sum(keep)
}

# Columns state what each number is, not what the step did to it:
#   curated_set          kept by curation and carrying reads
#   total_quantified     at least one count in at least one sample
#   extended_annotations every row Bambu quantified against the extended annotation
#   recovery             curated_set / total_quantified
validation_summary <- data.frame(
  feature_type         = c("transcript", "gene"),
  curated_set          = c(retained_expressed(se, valid_tx_ids, tx_expressed, "transcript"),
                           retained_expressed(seGene, valid_gene_ids, gene_expressed, "gene")),
  total_quantified     = c(sum(tx_expressed), sum(gene_expressed)),
  extended_annotations = c(nrow(se), nrow(seGene)),
  stringsAsFactors     = FALSE
)
validation_summary$recovery <-
  validation_summary$curated_set / validation_summary$total_quantified

write.csv(validation_summary,
          file.path(output_dir, "validation_summary.csv"), row.names = FALSE)

# --- Save the filtered objects for downstream reuse ---
saveRDS(se_filtered, file = file.path(output_dir, "se_multiSample_validated.rds"))
saveRDS(seGene_filtered, file = file.path(output_dir, "seGene_multiSample_validated.rds"))

# --- Regenerate the plots ---
# Mirrors the plotBambu calls in bambu.R so the post-refinement figures are
# directly comparable with the raw ones.
plots <- list(
  heatmap_transcript_validated = plotBambu(se_filtered, type = "heatmap", group.variable = "groupVar"),
  heatmap_gene_validated       = plotBambu(seGene_filtered, type = "heatmap", group.variable = "groupVar"),
  pca_validated                = plotBambu(se_filtered, type = "pca"),
  pca_grouped_validated        = plotBambu(se_filtered, type = "pca", group.variable = "groupVar")
)

lapply(names(plots), function(x) {
  file_path <- file.path(output_dir, paste0(x, ".png"))
  png(file_path, width = 5, height = 5, units = "in", res = 600)
  print(plots[[x]])
  dev.off()
})

cat("Post-refinement plots written.\n")
