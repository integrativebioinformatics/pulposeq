#!/usr/bin/env Rscript

Sys.setenv(HOME = tempdir())
Sys.setenv(BIOMART_CACHE = file.path(tempdir(), "biomart_cache"))

cat("HOME defined as:", Sys.getenv("HOME"), "\n")
cat("biomart cache:", Sys.getenv("BIOMART_CACHE"), "\n")

# Load required libraries
suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(biomaRt)
    library(GenomicRanges)
    library(rtracklayer)
    library(optparse)
    library(httr2)	
})

# Define command line options
option_list <- list(
    make_option(c("--transcript_counts"), type="character", default=NULL,
                help="Path to Bambu transcript counts file", metavar="character"),
    make_option(c("--gene_counts"), type="character", default=NULL,
                help="Path to Bambu gene counts file", metavar="character"),
    make_option(c("--gtf_file"), type="character", default=NULL,
                help="Path to Bambu GTF annotations file", metavar="character"),
    make_option(c("--ensembl_organism_dataset"), type="character", default=NULL,
                help="Ensembl organism genome dataset:", metavar="character"),
    make_option(c("--ensembl_version"), type="integer", default = NULL,
                help="Ensembl release version", metavar="integer")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check required arguments
if (is.null(opt$transcript_counts) || is.null(opt$gene_counts) || is.null(opt$gtf_file)) {
    print_help(opt_parser)
    stop("All input files must be specified.", call.=FALSE)
}

# Main analysis
cat("Starting transcript annotation analysis...\n")

# Read transcript counts
cat("Reading transcript counts...\n")
tx <- read_table(opt$transcript_counts, show_col_types = FALSE)

# Identify organism and release version
cat("Identifying organism and release version...\n")
organism <- opt$ensembl_organism_dataset
ensembl_version <- opt$ensembl_version

# Extract Ensembl transcript IDs
ens_ids <- tx$TXNAME[startsWith(tx$TXNAME, "ENS")]
cat(paste("Found", length(ens_ids), "Ensembl transcript IDs\n"))

# Apply biomaRt to gather metadata
cat("Connecting to Ensembl biomaRt...\n")

ensembl <- useEnsembl(biomart="genes", dataset=organism, version=ensembl_version)

attributes <- c(
    "chromosome_name",
    "gene_id", "gene_id_version",
    "transcript_id", "transcript_id_version",
    "external_transcript_name", "external_gene_name",
    "strand", "transcript_start", "transcript_end",
    "transcript_length", "gene_biotype", "transcript_biotype"
)

cat("Retrieving transcript metadata from biomaRt...\n")

# Verify Ensembl transcript IDs
if (any(grepl("\\.", ens_ids))) {
    version_suffix <- TRUE
    filters <- "transcript_id_version"
    ens_tx <- getBM(attributes = attributes, filters = filters, values = ens_ids, mart = ensembl)
} else {
    version_suffix <- FALSE
    filters <- "transcript_id"
    ens_tx <- getBM(attributes = attributes, filters = filters, values = ens_ids, mart = ensembl)
}

# Fix strand notation
ens_tx$strand <- ifelse(ens_tx$strand=="-1", "-", "+")

# Write transcriptome metadata
write.csv(ens_tx, "annotated_transcriptome_metadata.csv", row.names=FALSE)
cat("Written annotated_transcriptome_metadata.csv\n")

# Process lncRNAs
cat("Processing lncRNAs...\n")
ens_lnc <- ens_tx[ens_tx$gene_biotype=="lncRNA", ]
if (version_suffix == TRUE) {
    ens_lnc_ids <- ens_lnc$transcript_id_version
} else {
    ens_lnc_ids <- ens_lnc$transcript_id
}

if (length(ens_lnc_ids) > 0) {
    # Get exon information for lncRNAs
    exon_attributes <- c("chromosome_name", "transcript_id", "transcript_id_version", "exon_id_id", 
                        "exon_chrom_start", "exon_chrom_end")
    
    exon_ens_lnc <- getBM(attributes=exon_attributes, filters=filters, 
                         values=ens_lnc_ids, mart=ensembl)
    
    # Calculate the number of exons per transcript
    exon_counts <- exon_ens_lnc %>%
        group_by(transcript_id_version) %>%
        summarize(num_exons = n_distinct(exon_id_id), .groups = 'drop')
    
    ens_lnc <- merge(ens_lnc, exon_counts, by.x="transcript_id_version", 
                    by.y="transcript_id_version", all.x=TRUE)
    
    write.csv(ens_lnc, "annotated_lncRNAs_metadata.csv", row.names=FALSE)
    cat("Written annotated_lncRNAs_metadata.csv\n")
    
    # Calculate exon lengths for lncRNAs
    if (nrow(exon_ens_lnc) > 0) {
        exons_gr <- makeGRangesFromDataFrame(exon_ens_lnc,
                                           keep.extra.columns = TRUE,
                                           start.field = "exon_chrom_start",
                                           end.field = "exon_chrom_end",
                                           seqnames.field="chromosome_name")
        
        exon_lengths <- data.frame(
            transcript_id = exons_gr$transcript_id,
            transcript_id_version = exons_gr$transcript_id_version,
            exon_id_id = exons_gr$exon_id_id,
            width = width(exons_gr)
        ) 
        
        write.csv(exon_lengths, "annotated_lncRNAs_exonlength.csv", row.names=FALSE)
        cat("Written annotated_lncRNAs_exonlength.csv\n")
    } else {
        # Create empty file if no exons found
        write.csv(data.frame(), "annotated_lncRNAs_exonlength.csv", row.names=FALSE)
    }
} else {
    # Create empty files if no lncRNAs found
    write.csv(data.frame(), "annotated_lncRNAs_metadata.csv", row.names=FALSE)
    write.csv(data.frame(), "annotated_lncRNAs_exonlength.csv", row.names=FALSE)
    cat("No lncRNAs found, created empty files\n")
}

# Process protein-coding transcripts
cat("Processing protein-coding transcripts...\n")
ens_pc <- ens_tx[ens_tx$gene_biotype=="protein_coding", ]
if (version_suffix == TRUE) {
    ens_pc_ids <- ens_pc$transcript_id_version
} else {
    ens_pc_ids <- ens_pc$transcript_id
}

if (length(ens_pc_ids) > 0) {
    # Get exon information for protein-coding transcripts
    exon_attributes <- c("chromosome_name", "transcript_id", "transcript_id_version", "exon_id_id", 
                        "exon_chrom_start", "exon_chrom_end")
    
    exon_ens_pc <- getBM(attributes=exon_attributes, filters=filters, 
                        values=ens_pc_ids, mart=ensembl)
    
    # Calculate the number of exons per transcript
    exon_counts_pc <- exon_ens_pc %>%
        group_by(transcript_id_version) %>%
        summarize(num_exons = n_distinct(exon_id_id), .groups = 'drop')
    
    ens_pc <- merge(ens_pc, exon_counts_pc, by.x="transcript_id_version", 
                   by.y="transcript_id_version", all.x=TRUE)
    
    write.csv(ens_pc, "annotated_protein-coding_metadata.csv", row.names=FALSE)
    cat("Written annotated_protein-coding_metadata.csv\n")
    
    # Calculate exon lengths for protein-coding transcripts
    if (nrow(exon_ens_pc) > 0) {
        exons_gr_pc <- makeGRangesFromDataFrame(exon_ens_pc,
                                              keep.extra.columns = TRUE,
                                              start.field = "exon_chrom_start",
                                              end.field = "exon_chrom_end",
                                              seqnames.field="chromosome_name")
        
        exon_lengths_pc <- data.frame(
            transcript_id = exons_gr_pc$transcript_id,
            transcript_id_version = exons_gr_pc$transcript_id_version,
            exon_id_id = exons_gr_pc$exon_id_id,
            width = width(exons_gr_pc)
        ) 
        
        write.csv(exon_lengths_pc, "annotated_protein-coding_exonlength.csv", row.names=FALSE)
        cat("Written annotated_protein-coding_exonlength.csv\n")
    } else {
        # Create empty file if no exons found
        write.csv(data.frame(), "annotated_protein-coding_exonlength.csv", row.names=FALSE)
    }
} else {
    # Create empty files if no protein-coding transcripts found
    write.csv(data.frame(), "annotated_protein-coding_metadata.csv", row.names=FALSE)
    write.csv(data.frame(), "annotated_protein-coding_exonlength.csv", row.names=FALSE)
    cat("No protein-coding transcripts found, created empty files\n")
}

# Export GTF and counts for annotated transcripts
cat("Processing GTF files and counts...\n")
gtf <- import(opt$gtf_file)

# Test not removing version suffixes

# get the IDs
if (version_suffix == TRUE) {
    tx_ids <- ens_tx$transcript_id_version
    lnc_ids <- ens_lnc$transcript_id_version
    pc_ids <- ens_pc$transcript_id_version
} else {
    tx_ids <- ens_tx$transcript_id
    lnc_ids <- ens_lnc$transcript_id
    pc_ids <- ens_pc$transcript_id
}

# Export annotated transcriptome GTF
ann_tx_gtf <- subset(gtf, transcript_id %in% tx_ids)
export(ann_tx_gtf, "bambu_annotated_transcriptome.gtf")
cat("Written bambu_annotated_transcriptome.gtf\n")

# Export annotated transcript counts
ann_tx_counts <- subset(tx, TXNAME %in% tx_ids)
write.csv(ann_tx_counts, "bambu_annotated_transcriptome_tx_counts.csv", row.names = FALSE)
cat("Written bambu_annotated_transcriptome_tx_counts.csv\n")

# Export lncRNA GTF
if (length(lnc_ids) > 0) {
    ann_lnc_gtf <- subset(gtf, transcript_id %in% lnc_ids)
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
    gn_ids <- ens_tx$gene_id_version
} else {
    gn_ids <- ens_tx$gene_id
}

ann_gn_counts <- subset(gn, GENEID %in% gn_ids)
write.csv(ann_gn_counts, "bambu_annotated_transcriptome_gene_counts.csv", row.names = FALSE)
cat("Written bambu_annotated_transcriptome_gene_counts.csv\n")

cat("Transcript annotation analysis completed successfully!\n")