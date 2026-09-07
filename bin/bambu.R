#!/usr/bin/env Rscript

# --- Load necessary libraries ---

library(ggplot2)
library(Rsamtools)
library(bambu)
library(readr)
library(BiocParallel)
library(rtracklayer)
library(GenomicRanges)
library(optparse)

# --- Parse command-line arguments ---
option_list <- list(
  make_option(c("-g", "--genome"), type = "character", default = NULL,
              help = "Path to the genome FASTA file", metavar = "character"),
  make_option(c("-a", "--annotation"), type = "character", default = NULL,
              help = "Path to the GTF annotation file", metavar = "character"),
  make_option(c("-b", "--bamfiles"), type = "character", default = NULL,
              help = "Path to a file containing a list of BAM files, one per line", metavar = "character"),
  make_option(c("-s", "--sampleinfo"), type = "character", default = NULL,
              help = "Path to the sample information Excel file", metavar = "character"),
  make_option(c("-n", "--ncores"), type = "numeric", default = 1,
              help = "Number of cores to use for parallel processing", metavar = "numeric"),
  make_option(c("-d", "--ndr"), type = "character", default = NULL,
              help = "Novel Discovery Rate cut-off (optional)", metavar = "character"),
  make_option(c("-t", "--stranded"), type = "character", default = "false",
              help = "Whether the library is stranded (true/false)", metavar = "character"),
  make_option(c("-o", "--outdir"), type = "character", default = "output",
              help = "Output directory", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if required arguments are provided
if (is.null(opt$genome) ||
    is.null(opt$annotation) || is.null(opt$bamfiles) || is.null(opt$sampleinfo)) {
  print_help(opt_parser)
  stop("Missing required arguments.")
}

# --- Define file paths ---
genomeSequence <- opt$genome
gtf.file <- opt$annotation
bamFiles <- readLines(opt$bamfiles)  # Read BAM file paths from the input file
sample_info_file <- opt$sampleinfo
output_dir <- opt$outdir
ncores <- opt$ncores  # Use the specified number of cores

# Parse ndr parameter: convert "null" string to NULL, otherwise to numeric
if (is.null(opt$ndr) || opt$ndr == "null") {
  ndr_cutoff <- NULL
} else {
  ndr_cutoff <- as.numeric(opt$ndr)
}

# Parse stranded parameter: accept the usual truthy spellings
stranded_flag <- tolower(opt$stranded) %in% c("true", "t", "yes", "1")

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Prepare annotations ---
annotation <- tryCatch({
  prepareAnnotations(gtf.file)
}, error = function(e) {
  stop(paste("Error preparing annotations:", e$message))
})

# --- Load BAM files ---
totalData <- tryCatch({
  BamFileList(bamFiles)
}, error = function(e) {
  stop(paste("Error loading BAM files:", e$message))
})

# --- Run bambu ---
se.multiSample <- tryCatch({
  # Build bambu parameters
  bambu_params <- list(
    ncore = ncores,
    reads = totalData,
    annotations = annotation,
    genome = genomeSequence,
    stranded = stranded_flag,
    trackReads = TRUE,
    opt.discovery = list(min.txScore.singleExon = 0)
  )
  
  # Only add NDR if it's not NULL
  if (!is.null(ndr_cutoff)) {
    bambu_params$NDR <- ndr_cutoff
  }
  
  # Call bambu with the constructed parameters
  do.call(bambu, bambu_params)
}, error = function(e) {
  stop(paste("Error running bambu:", e$message))
})

# --- Add sample metadata ---
sample_info <- tryCatch({
  read.table(sample_info_file, sep="\t")
}, error = function(e) {
  stop(paste("Error reading sample information file:", e$message))
})

# basenames of BAM paths (in the order listed in bamFiles)
bam_basenames     <- basename(bamFiles)
bam_names_noext   <- sub('\\.bam(?:\\.gz)?$', '', bam_basenames, ignore.case = TRUE)

# determine which column in sample_info contains the BAM filename
# common case for sampinfo_samplesheet.tsv: column 2 is filename (col1 = group)
if (ncol(sample_info) >= 2) {
  sample_filenames <- as.character(sample_info[[2]])
} else {
  sample_filenames <- as.character(sample_info[[1]])
}
sample_basenames_noext <- sub('\\.bam(?:\\.gz)?$', '', basename(sample_filenames), ignore.case = TRUE)

# match sample rows to bamFiles by basename without extension
match_idx <- match(sample_basenames_noext, bam_names_noext)

# try fallback: match by full basename (with extension) for any NAs
na_rows <- which(is.na(match_idx))
if (length(na_rows) > 0) {
  fallback_idx <- match(basename(sample_filenames[na_rows]), bam_basenames)
  match_idx[na_rows] <- fallback_idx
}

if (all(is.na(match_idx))) {
  warning("Could not match any sample_info filenames to bamlist entries — keeping original sample_info order")
} else {
  # order sample_info rows by the matched index (unmatched rows will be placed at the end)
  ord <- order(ifelse(is.na(match_idx), Inf, match_idx), na.last = TRUE)
  sample_info <- sample_info[ord, , drop = FALSE]
}

# Optional sanity check
if (nrow(sample_info) != length(bamFiles)) {
  warning(sprintf("Number of sample_info rows (%d) differs from number of BAMs (%d). Check for missing or extra entries.",
                  nrow(sample_info), length(bamFiles)))
}

colData(se.multiSample)$group <- as.factor(sample_info[[1]])
colData(se.multiSample)$groupVar <- sample_info[[1]]

# --- Transcript-level analysis ---

# Convert to gene expression
seGene.multiSample <- transcriptToGeneExpression(se.multiSample)

# Add sample metadata to gene-level object
colData(seGene.multiSample)$group <- as.factor(sample_info[[1]])
colData(seGene.multiSample)$groupVar <- sample_info[[1]]


# --- Save SummarizedExperiment objects ---
saveRDS(se.multiSample, file = file.path(output_dir, "se_multiSample.rds"))
saveRDS(seGene.multiSample, file = file.path(output_dir, "seGene_multiSample.rds"))

# --- Save Bambu's official outputs ---
writeBambuOutput(se.multiSample, output_dir, prefix = "BambuOutput_")

# --- Save the new transcripts
transcript_annotations <- rtracklayer::import("BambuOutput_extended_annotations.gtf")

newtx_gtf <- transcript_annotations[grep("^BambuTx", transcript_annotations$transcript_id)]
rtracklayer::export(newtx_gtf, "bambu_novel_transcripts.gtf")

# --- Create and save plots ---

# Create a list of plots
plots <- list(
  heatmap_transcript = plotBambu(se.multiSample, type = "heatmap", group.variable="groupVar"),
  heatmap_gene = plotBambu(seGene.multiSample, type = "heatmap", group.variable="groupVar"),
  pca = plotBambu(se.multiSample, type = "pca"),
  pca_grouped = plotBambu(se.multiSample, type = "pca", group.variable = "groupVar")
)

# Use lapply to iterate and save each plot with its corresponding name
lapply(names(plots), function(x) {
  # Determine file extension based on plot type
  file_ext <- "png"

  # Save the plot
  file_path <- file.path(output_dir, paste0(x, ".", file_ext))

  # Use the appropriate device based on file extension
  png(file_path, width = 5, height = 5, units = "in", res = 600)  # Adjust dimensions and resolution as needed

  print(plots[[x]])  # Print the plot to the device
  dev.off()
})

# Function to reorder GTF data
reorder_gtf <- function(gtf_data) {
  return(gtf_data[order(gtf_data$transcript_id, gtf_data$type == "transcript", decreasing = TRUE)])
}

# Reorder the GTF data
sorted_novel <- reorder_gtf(newtx_gtf)
sorted_ext <- reorder_gtf(transcript_annotations)

# Write the reordered data to the output GTF file
rtracklayer::export(sorted_novel, "bambu_novel_transcripts.gtf")
rtracklayer::export(sorted_ext, "BambuOutput_extended_annotations.gtf")
