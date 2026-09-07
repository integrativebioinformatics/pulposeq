#!/usr/bin/env Rscript

Sys.setenv(HOME = tempdir())

# Load required libraries
suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(GenomicRanges)
    library(rtracklayer)
    library(optparse)
})

# Shared reference-GTF helpers, staged alongside this script by the module
source("gtf_annotation_utils.R")

# Define command line options
option_list <- list(
    make_option(c("--transcript_counts"), type="character", default=NULL,
                help="Path to Bambu transcript counts file", metavar="character"),
    make_option(c("--gene_counts"), type="character", default=NULL,
                help="Path to Bambu gene counts file", metavar="character"),
    make_option(c("--gtf_file"), type="character", default=NULL,
                help="Path to Bambu GTF annotations file", metavar="character"),
    make_option(c("--annotation"), type="character", default=NULL,
                help="Path to the reference annotation GTF (Ensembl or GENCODE)", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check required arguments
if (is.null(opt$transcript_counts) || is.null(opt$gene_counts) ||
    is.null(opt$gtf_file) || is.null(opt$annotation)) {
    print_help(opt_parser)
    stop("All input files must be specified.", call.=FALSE)
}

# Main analysis
cat("Starting transcript annotation analysis...\n")

# Read transcript counts
cat("Reading transcript counts...\n")
tx <- read_table(opt$transcript_counts, show_col_types = FALSE)

# Read all transcript metadata out of the reference annotation
reference <- read_reference_gtf(opt$annotation)
ref_tx <- reference$tx

# Identify which assembled transcripts are known, i.e. present in the reference.
# Bambu names novel transcripts BambuTx* and passes reference identifiers through
# unchanged, so membership in the supplied annotation is the definition of "known".
# Matching on whichever id form the counts use keeps Ensembl (unversioned) and
# GENCODE (versioned) references both working.
version_suffix <- any(tx$TXNAME %in% ref_tx$ensembl_transcript_id_version &
                      !tx$TXNAME %in% ref_tx$ensembl_transcript_id)

if (version_suffix) {
    id_col <- "ensembl_transcript_id_version"
} else {
    id_col <- "ensembl_transcript_id"
}

ens_ids <- tx$TXNAME[tx$TXNAME %in% ref_tx[[id_col]]]
cat(paste("Found", length(ens_ids), "known transcript IDs in the reference annotation\n"))

if (length(ens_ids) == 0) {
    stop(paste0("None of the assembled transcripts matched the reference annotation.\n",
                "  Check that --annotation is the same GTF that was supplied to Bambu."),
         call. = FALSE)
}

# Subset the reference metadata to the assembled known transcripts
ens_tx <- ref_tx[ref_tx[[id_col]] %in% ens_ids, ]

# Exon table for the same set
ref_exons <- reference$exons
ens_exons <- ref_exons[ref_exons[[id_col]] %in% ens_ids, ]

# Write transcriptome metadata
write.csv(ens_tx, "annotated_transcriptome_metadata.csv", row.names=FALSE)
cat("Written annotated_transcriptome_metadata.csv\n")

#' Write the metadata and exon-length tables for one biotype.
#'
#' Mirrors the previous per-biotype blocks, including the empty-file fallbacks so
#' downstream processes always have something to consume.
process_biotype <- function(biotype, metadata_file, exonlength_file, label) {
    cat(paste0("Processing ", label, "...\n"))

    subset_tx <- ens_tx[ens_tx$gene_biotype == biotype, ]
    subset_ids <- subset_tx[[id_col]]

    if (length(subset_ids) == 0) {
        write.csv(data.frame(), metadata_file, row.names=FALSE)
        write.csv(data.frame(), exonlength_file, row.names=FALSE)
        cat(paste0("No ", label, " found, created empty files\n"))
        return(subset_tx)
    }

    subset_exons <- ens_exons[ens_exons[[id_col]] %in% subset_ids, ]

    # Number of exons per transcript
    exon_counts <- exon_counts_per_transcript(subset_exons)

    subset_tx <- merge(subset_tx, exon_counts,
                       by.x="ensembl_transcript_id_version",
                       by.y="ensembl_transcript_id_version", all.x=TRUE)

    write.csv(subset_tx, metadata_file, row.names=FALSE)
    cat(paste("Written", metadata_file, "\n"))

    if (nrow(subset_exons) > 0) {
        exon_lengths <- data.frame(
            ensembl_transcript_id         = subset_exons$ensembl_transcript_id,
            ensembl_transcript_id_version = subset_exons$ensembl_transcript_id_version,
            ensembl_exon_id               = subset_exons$ensembl_exon_id,
            width                         = subset_exons$exon_chrom_end - subset_exons$exon_chrom_start + 1
        )
        write.csv(exon_lengths, exonlength_file, row.names=FALSE)
        cat(paste("Written", exonlength_file, "\n"))
    } else {
        write.csv(data.frame(), exonlength_file, row.names=FALSE)
    }

    subset_tx
}

ens_lnc <- process_biotype("lncRNA",
                           "annotated_lncRNAs_metadata.csv",
                           "annotated_lncRNAs_exonlength.csv",
                           "lncRNAs")

ens_pc <- process_biotype("protein_coding",
                          "annotated_protein-coding_metadata.csv",
                          "annotated_protein-coding_exonlength.csv",
                          "protein-coding transcripts")

# Export GTF and counts for annotated transcripts
cat("Processing GTF files and counts...\n")
gtf <- import(opt$gtf_file)

# get the IDs
tx_ids  <- ens_tx[[id_col]]
lnc_ids <- if (nrow(ens_lnc) > 0) ens_lnc[[id_col]] else character(0)
pc_ids  <- if (nrow(ens_pc) > 0) ens_pc[[id_col]] else character(0)

# Reference attributes to write into the exported GTFs, keyed by transcript id
known_attrs <- list(
    transcript_status   = rep("known", nrow(ens_tx)),
    gene_name           = ens_tx$external_gene_name,
    gene_biotype        = ens_tx$gene_biotype,
    transcript_name     = ens_tx$external_transcript_name,
    transcript_biotype  = ens_tx$transcript_biotype
)

# Export annotated transcriptome GTF
ann_tx_gtf <- subset(gtf, transcript_id %in% tx_ids)
ann_tx_gtf <- annotate_gtf(ann_tx_gtf, tx_ids, known_attrs)
export(ann_tx_gtf, "bambu_annotated_transcriptome.gtf")
cat("Written bambu_annotated_transcriptome.gtf\n")

# Export annotated transcript counts
ann_tx_counts <- subset(tx, TXNAME %in% tx_ids)
write.csv(ann_tx_counts, "bambu_annotated_transcriptome_tx_counts.csv", row.names = FALSE)
cat("Written bambu_annotated_transcriptome_tx_counts.csv\n")

# Export lncRNA GTF
if (length(lnc_ids) > 0) {
    ann_lnc_gtf <- subset(gtf, transcript_id %in% lnc_ids)
    ann_lnc_gtf <- annotate_gtf(ann_lnc_gtf, tx_ids, known_attrs)
    export(ann_lnc_gtf, "bambu_annotated_lncRNAs.gtf")
    cat("Written bambu_annotated_lncRNAs.gtf\n")
} else {
    # Create empty GTF file
    export(GRanges(), "bambu_annotated_lncRNAs.gtf")
    cat("No lncRNAs found, created empty bambu_annotated_lncRNAs.gtf\n")
}

# Export protein-coding GTF
if (length(pc_ids) > 0) {
    ann_pc_gtf <- subset(gtf, transcript_id %in% pc_ids)
    ann_pc_gtf <- annotate_gtf(ann_pc_gtf, tx_ids, known_attrs)
    export(ann_pc_gtf, "bambu_annotated_mRNAs.gtf")
    cat("Written bambu_annotated_mRNAs.gtf\n")
} else {
    # Create empty GTF file
    export(GRanges(), "bambu_annotated_mRNAs.gtf")
    cat("No protein-coding transcripts found, created empty bambu_annotated_mRNAs.gtf\n")
}

# Process gene counts
cat("Processing gene counts...\n")
gn <- read_table(opt$gene_counts, show_col_types = FALSE)
if (version_suffix) {
    gn_ids <- ens_tx$ensembl_gene_id_version
} else {
    gn_ids <- ens_tx$ensembl_gene_id
}

ann_gn_counts <- subset(gn, GENEID %in% gn_ids)
write.csv(ann_gn_counts, "bambu_annotated_transcriptome_gene_counts.csv", row.names = FALSE)
cat("Written bambu_annotated_transcriptome_gene_counts.csv\n")

cat("Transcript annotation analysis completed successfully!\n")
