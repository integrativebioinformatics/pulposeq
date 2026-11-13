#!/usr/bin/env Rscript

# Sample-level PCA and expression heatmap, before and after curation, and the
# curation yield.
#
# Plots are transcript-level only, drawn with bambu::plotBambu so both versions
# use the same method (CPM, transcripts above the median standard deviation,
# log2(CPM + 1); PCA on those values, heatmap of the sample-by-sample Spearman
# correlation). The counts they come from follow the expression source:
#   - Bambu (default): Bambu's transcript object, every model it quantified, then
#     the curated transcripts
#   - oarfish (--oarfish_cpm given): oarfish's CPM table, every model oarfish
#     quantified, then the curated transcripts. CPM is computed once on all
#     models and subset, as in the curated tables.
#
# The yield (validation_summary.csv) is measured against the same source: a
# feature counts as quantified if it has any counts summed across all samples.

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
  make_option("--oarfish_counts", type = "character", default = NULL,
              help = "oarfish transcript counts, every model (oarfish_counts_transcript.txt)", metavar = "character"),
  make_option("--oarfish_counts_gene", type = "character", default = NULL,
              help = "oarfish gene counts, every gene (oarfish_counts_gene.txt)", metavar = "character"),
  make_option("--oarfish_cpm", type = "character", default = NULL,
              help = "oarfish transcript CPM, every model (oarfish_CPM_transcript.txt)", metavar = "character"),
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

has <- function(x) !is.null(x) && nzchar(x) && file.exists(x)
use_oarfish <- has(opt$oarfish_cpm)
if (use_oarfish && !(has(opt$oarfish_counts) && has(opt$oarfish_counts_gene))) {
  stop("--oarfish_cpm needs --oarfish_counts and --oarfish_counts_gene as well.")
}
expression_source <- if (use_oarfish) "oarfish" else "bambu"
cat(sprintf("Expression source: %s\n", expression_source))

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
                       "Check that the counts file and the RDS come from the same run."),
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

# --- Per-feature totals from the expression source ----------------------------
# Bambu: the counts assay of its objects. oarfish: its merged tables, which hold
# every model oarfish quantified.
read_counts_table <- function(path) {
  tab <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  m <- as.matrix(tab[, setdiff(names(tab), c("TXNAME", "GENEID")), drop = FALSE])
  rownames(m) <- tab[[1]]
  m
}

se_totals <- function(object, label) {
  assay_names <- SummarizedExperiment::assayNames(object)
  if (length(assay_names) == 0) {
    stop(sprintf("The %s-level RDS object carries no assay to test for counts.", label))
  }
  counts_assay <- if ("counts" %in% assay_names) "counts" else assay_names[1]
  setNames(rowSums(as.matrix(SummarizedExperiment::assay(object, counts_assay)), na.rm = TRUE),
           rownames(object))
}

if (use_oarfish) {
  tx_totals   <- rowSums(read_counts_table(opt$oarfish_counts), na.rm = TRUE)
  gene_totals <- rowSums(read_counts_table(opt$oarfish_counts_gene), na.rm = TRUE)
} else {
  tx_totals   <- se_totals(se, "transcript")
  gene_totals <- se_totals(seGene, "gene")
}

# --- Record the validation yield as a result --------------------------------
# A headline number for the run, written as data rather than left to be scraped
# back out of the console log.
#
# The denominator is every feature with any counts in the expression source.
# Bambu and oarfish quantify against the extended annotation, so their tables
# carry a row for every reference feature, including those with no reads;
# counting those would count the reference itself as a loss. The denominator is
# deliberately looser than PREFILTER's --min_total_counts rule, so the recovery
# figure also shows what that expression filter removed.
yield_row <- function(totals, ids, label) {
  expressed <- names(totals)[totals > 0]
  cat(sprintf("%d of %d %ss carry counts in %s.\n", length(expressed), length(totals),
              label, expression_source))
  silent <- setdiff(ids, expressed)
  if (length(silent) > 0) {
    warning(sprintf(paste("%d validated %s IDs have no %s counts in any sample and were",
                          "excluded from the validation summary."), length(silent), label,
                    expression_source))
  }
  c(curated_set = length(intersect(ids, expressed)), total_quantified = length(expressed))
}

# Columns state what each number is, not what the step did to it:
#   curated_set          curated and carrying counts in the expression source
#   total_quantified     any counts in the expression source, summed across samples
#   extended_annotations every row of Bambu's extended annotation
#   recovery             curated_set / total_quantified
#   expression_source    bambu or oarfish
tx_yield   <- yield_row(tx_totals, valid_tx_ids, "transcript")
gene_yield <- yield_row(gene_totals, valid_gene_ids, "gene")
validation_summary <- data.frame(
  feature_type         = c("transcript", "gene"),
  curated_set          = c(tx_yield[["curated_set"]], gene_yield[["curated_set"]]),
  total_quantified     = c(tx_yield[["total_quantified"]], gene_yield[["total_quantified"]]),
  extended_annotations = c(nrow(se), nrow(seGene)),
  stringsAsFactors     = FALSE
)
validation_summary$recovery <-
  validation_summary$curated_set / validation_summary$total_quantified
validation_summary$expression_source <- expression_source

write.csv(validation_summary,
          file.path(output_dir, "validation_summary.csv"), row.names = FALSE)

# --- Save the filtered objects for downstream reuse ---
saveRDS(se_filtered, file = file.path(output_dir, "se_multiSample_validated.rds"))
saveRDS(seGene_filtered, file = file.path(output_dir, "seGene_multiSample_validated.rds"))

# --- The objects the plots are drawn from -------------------------------------
if (use_oarfish) {
  cpm <- read_counts_table(opt$oarfish_cpm)

  # Sample groups come from Bambu's object. Bambu names its samples after the
  # BAM files ("<sample>.bam"); oarfish uses the bare sample name.
  bambu_info <- as.data.frame(SummarizedExperiment::colData(se))
  bambu_names <- if ("name" %in% names(bambu_info)) as.character(bambu_info$name) else colnames(se)
  bambu_names <- sub("\\.bam$", "", bambu_names)
  idx <- match(colnames(cpm), bambu_names)
  if (anyNA(idx)) {
    stop(sprintf("oarfish samples not found in Bambu's object: %s",
                 paste(colnames(cpm)[is.na(idx)], collapse = ", ")))
  }
  group <- if ("groupVar" %in% names(bambu_info)) bambu_info$groupVar[idx] else colnames(cpm)

  plot_all <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(CPM = cpm),
    colData = S4Vectors::DataFrame(name = colnames(cpm), groupVar = group,
                                   row.names = colnames(cpm))
  )
  plot_curated <- plot_all[rownames(plot_all) %in% valid_tx_ids, , drop = FALSE]
  cat(sprintf("Plotting oarfish CPM: %d models, %d curated.\n",
              nrow(plot_all), nrow(plot_curated)))
} else {
  plot_all     <- se
  plot_curated <- se_filtered
}

# --- Draw the plots -------------------------------------------------------------
# Every model first ("_all_models"), then the curated transcripts ("_validated").
plots <- list(
  pca_all_models                = plotBambu(plot_all, type = "pca"),
  pca_grouped_all_models        = plotBambu(plot_all, type = "pca", group.variable = "groupVar"),
  heatmap_transcript_all_models = plotBambu(plot_all, type = "heatmap", group.variable = "groupVar"),
  pca_validated                 = plotBambu(plot_curated, type = "pca"),
  pca_grouped_validated         = plotBambu(plot_curated, type = "pca", group.variable = "groupVar"),
  heatmap_transcript_validated  = plotBambu(plot_curated, type = "heatmap", group.variable = "groupVar")
)

lapply(names(plots), function(x) {
  file_path <- file.path(output_dir, paste0(x, ".png"))
  png(file_path, width = 5, height = 5, units = "in", res = 600)
  print(plots[[x]])
  dev.off()
})

cat(sprintf("Transcript-level plots written from %s counts.\n", expression_source))
