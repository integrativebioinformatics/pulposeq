#!/usr/bin/env Rscript

# PREFILTER: decide which transcripts every later step of the pipeline sees.
#
# Runs once, right after Bambu, and before gffread, gffcompare and CPC2. The kept
# transcript list it writes is the only one: the filtered count tables, the novel
# models sent to classification, and everything curated from them all follow it.
#
# Three rules, applied in this order so each removal is counted under one reason:
#
#   1. Unstranded novel models are removed. pulposeq never admits unoriented reads,
#      so a model without a strand cannot be placed in the sense/antisense
#      classification that novel lncRNA calls depend on.
#   2. Novel models shorter than 500 nt are removed. 500 nt is the consensus lower
#      bound for a long non-coding RNA; it is applied to every novel candidate so
#      the three novel categories stay comparable to one another.
#   3. Any transcript, reference or novel, whose counts summed across all samples
#      fall below --min_total_counts is removed. This is a lenient detection
#      filter: it asks whether there is any real evidence the transcript is
#      expressed, not whether it is suitable for differential expression, which is
#      left to the user's own downstream analysis.
#
# Rules 1 and 2 apply to novel models only. Reference transcripts are never removed
# for strand or length.
#
# The expression counts come from Bambu, or from oarfish when it runs. Which source
# decided is written to the filter summary.

suppressPackageStartupMessages(library(optparse))

MIN_NOVEL_LENGTH <- 500L
NOVEL_TX_PREFIX  <- "BambuTx"

option_list <- list(
    make_option("--counts_transcript", type = "character", default = NULL,
                help = "BambuOutput_counts_transcript.txt"),
    make_option("--counts_gene", type = "character", default = NULL,
                help = "BambuOutput_counts_gene.txt"),
    make_option("--cpm_transcript", type = "character", default = NULL,
                help = "BambuOutput_CPM_transcript.txt"),
    make_option("--full_length", type = "character", default = NULL,
                help = "BambuOutput_fullLengthCounts_transcript.txt"),
    make_option("--unique", type = "character", default = NULL,
                help = "BambuOutput_uniqueCounts_transcript.txt"),
    make_option("--novel_gtf", type = "character", default = NULL,
                help = "bambu_novel_transcripts.gtf: Bambu's novel models"),
    make_option("--expression_source", type = "character", default = "bambu",
                help = "Which counts decide expression: bambu or oarfish [default %default]"),
    make_option("--oarfish_counts", type = "character", default = NULL,
                help = "oarfish count table (TXNAME, GENEID, one column per sample); required with --expression_source oarfish"),
    make_option("--min_total_counts", type = "double", default = 5,
                help = "Minimum counts summed across all samples [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

required <- c("counts_transcript", "counts_gene", "cpm_transcript",
              "full_length", "unique", "novel_gtf")
missing_args <- required[vapply(required, function(x) is.null(opt[[x]]), logical(1))]
if (length(missing_args) > 0) {
    stop("Missing required arguments: ", paste(missing_args, collapse = ", "), call. = FALSE)
}
if (!opt$expression_source %in% c("bambu", "oarfish")) {
    stop("--expression_source must be 'bambu' or 'oarfish', got: ", opt$expression_source,
         call. = FALSE)
}
if (opt$expression_source == "oarfish" && is.null(opt$oarfish_counts)) {
    stop("--expression_source oarfish needs --oarfish_counts.", call. = FALSE)
}
if (is.na(opt$min_total_counts) || opt$min_total_counts < 0) {
    stop("--min_total_counts must be a number >= 0, got: ", opt$min_total_counts, call. = FALSE)
}

# --- Readers -------------------------------------------------------------------

#' Read a Bambu-style count table: identifier column(s) first, one column per
#' sample after. check.names = FALSE keeps sample names exactly as Bambu wrote them.
read_table_tsv <- function(path, id_cols) {
    tab <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
    absent <- setdiff(id_cols, names(tab))
    if (length(absent) > 0) {
        stop(sprintf("%s has no %s column. Found: %s", basename(path),
                     paste(absent, collapse = "/"), paste(names(tab), collapse = ", ")),
             call. = FALSE)
    }
    tab
}

sample_columns <- function(tab, id_cols) setdiff(names(tab), id_cols)

#' Write a table in the same layout Bambu uses: tab-separated, unquoted, no row names.
write_table_tsv <- function(tab, path) {
    utils::write.table(tab, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

#' Strand and spliced length of each novel model, read straight from the GTF text.
#'
#' Parsed as text rather than imported, so the filtered GTF written later keeps
#' Bambu's lines exactly as they were.
read_novel_models <- function(path) {
    lines   <- readLines(path)
    is_body <- nzchar(lines) & !startsWith(lines, "#")
    fields  <- strsplit(lines[is_body], "\t", fixed = TRUE)
    if (any(lengths(fields) < 9)) {
        stop(basename(path), " has lines with fewer than 9 tab-separated fields.", call. = FALSE)
    }
    type   <- vapply(fields, `[`, character(1), 3)
    start  <- as.integer(vapply(fields, `[`, character(1), 4))
    end    <- as.integer(vapply(fields, `[`, character(1), 5))
    strand <- vapply(fields, `[`, character(1), 7)
    attrs  <- vapply(fields, `[`, character(1), 9)

    has_tx <- grepl('transcript_id "', attrs, fixed = TRUE)
    tx_id  <- rep(NA_character_, length(attrs))
    tx_id[has_tx] <- sub('.*transcript_id "([^"]+)".*', "\\1", attrs[has_tx])

    is_exon <- type == "exon" & has_tx
    tx_length <- tapply(end[is_exon] - start[is_exon] + 1L, tx_id[is_exon], sum)
    tx_strand <- tapply(strand[has_tx], tx_id[has_tx], function(s) s[1])

    line_tx <- rep(NA_character_, length(lines))
    line_tx[is_body] <- tx_id

    list(lines = lines, is_body = is_body, line_tx = line_tx,
         length = tx_length, strand = tx_strand)
}

# --- Load ----------------------------------------------------------------------

TX_ID   <- c("TXNAME", "GENEID")
GENE_ID <- "GENEID"

cat("Reading Bambu tables...\n")
tx      <- read_table_tsv(opt$counts_transcript, TX_ID)
gene    <- read_table_tsv(opt$counts_gene, GENE_ID)
cpm     <- read_table_tsv(opt$cpm_transcript, TX_ID)
full    <- read_table_tsv(opt$full_length, TX_ID)
uniq    <- read_table_tsv(opt$unique, TX_ID)
models  <- read_novel_models(opt$novel_gtf)

is_novel <- startsWith(tx$TXNAME, NOVEL_TX_PREFIX)
novel_ids <- tx$TXNAME[is_novel]

not_in_gtf <- setdiff(novel_ids, names(models$strand))
if (length(not_in_gtf) > 0) {
    stop(sprintf("%d novel transcripts in the count table are missing from %s, e.g. %s. ",
                 length(not_in_gtf), basename(opt$novel_gtf), paste(head(not_in_gtf, 3), collapse = ", ")),
         "The two come from the same Bambu run and should match.", call. = FALSE)
}
no_exons <- setdiff(novel_ids, names(models$length))
if (length(no_exons) > 0) {
    stop(sprintf("%d novel transcripts have no exon records in %s, e.g. %s.",
                 length(no_exons), basename(opt$novel_gtf), paste(head(no_exons, 3), collapse = ", ")),
         call. = FALSE)
}

# --- Expression source ---------------------------------------------------------

if (opt$expression_source == "bambu") {
    source_tab <- tx
} else {
    cat("Reading oarfish counts...\n")
    source_tab <- read_table_tsv(opt$oarfish_counts, "TXNAME")
}
source_samples <- sample_columns(source_tab, TX_ID)
if (length(source_samples) == 0) {
    stop("The expression-source table has no sample columns.", call. = FALSE)
}
source_total <- rowSums(as.matrix(source_tab[, source_samples, drop = FALSE]))
names(source_total) <- source_tab$TXNAME

# A transcript the source never saw (with oarfish: one left out of its index) has
# no evidence of expression, so it counts as zero rather than being kept by default.
total <- unname(source_total[tx$TXNAME])
total[is.na(total)] <- 0

# --- Rules ---------------------------------------------------------------------

novel_strand <- unname(models$strand[tx$TXNAME])
novel_length <- unname(models$length[tx$TXNAME])

# Counts from Bambu and oarfish are model estimates, so a transcript backed by
# exactly 5 reads can come back as 4.9999999999999. The tolerance keeps
# floating-point noise, and the rounding of a table written to text, from deciding
# whether a transcript is kept.
COUNT_TOLERANCE <- 1e-6

unstranded <- is_novel & !(novel_strand %in% c("+", "-"))
too_short  <- is_novel & !unstranded & novel_length < MIN_NOVEL_LENGTH
too_short[is.na(too_short)] <- FALSE
low_expr   <- !unstranded & !too_short & total < opt$min_total_counts - COUNT_TOLERANCE
keep       <- !(unstranded | too_short | low_expr)

kept_tx    <- tx$TXNAME[keep]
kept_genes <- unique(tx$GENEID[keep])

# --- Write filtered tables -----------------------------------------------------
#
# File names match the ones the rest of the pipeline already reads. Counts, CPM and
# gene tables follow the kept list exactly. The full-length and unique tables keep
# their existing meaning: a row is a kept transcript with at least one such read.

has_reads <- function(tab) {
    m <- as.matrix(tab[, sample_columns(tab, TX_ID), drop = FALSE])
    rowSums(m != 0) > 0
}

write_table_tsv(tx[keep, , drop = FALSE], "BambuOutput_counts_transcript_filter.txt")
write_table_tsv(cpm[cpm$TXNAME %in% kept_tx, , drop = FALSE], "BambuOutput_CPM_transcript_filter.txt")
write_table_tsv(full[full$TXNAME %in% kept_tx & has_reads(full), , drop = FALSE],
                "BambuOutput_fullLengthCounts_transcript_filter.txt")
write_table_tsv(uniq[uniq$TXNAME %in% kept_tx & has_reads(uniq), , drop = FALSE],
                "BambuOutput_uniqueCounts_transcript_filter.txt")
write_table_tsv(gene[gene$GENEID %in% kept_genes, , drop = FALSE], "BambuOutput_counts_gene_filter.txt")

writeLines(kept_tx, "prefilter_kept_transcripts.txt")

# Novel models for classification: header lines plus every line of a kept model,
# copied as Bambu wrote them.
keep_line <- !models$is_body | (models$line_tx %in% kept_tx)
writeLines(models$lines[keep_line], "bambu_novel_transcripts_prefiltered.gtf")

# --- Summary -------------------------------------------------------------------

count_split <- function(flag) {
    c(reference = sum(flag & !is_novel), novel = sum(flag & is_novel), total = sum(flag))
}
summary_rows <- rbind(
    transcripts_in                 = count_split(rep(TRUE, nrow(tx))),
    removed_unstranded             = count_split(unstranded),
    removed_shorter_than_500nt     = count_split(too_short),
    removed_below_min_total_counts = count_split(low_expr),
    transcripts_kept               = count_split(keep)
)
summary_tab <- data.frame(step = rownames(summary_rows), summary_rows,
                          row.names = NULL, check.names = FALSE)
summary_tab <- rbind(summary_tab,
                     data.frame(step = c("genes_in", "genes_kept"),
                                reference = NA, novel = NA,
                                total = c(nrow(gene), sum(gene$GENEID %in% kept_genes))))

con <- file("filter_summary.tsv", open = "w")
writeLines(c(sprintf("# expression_source: %s", opt$expression_source),
             sprintf("# min_total_counts: %s", format(opt$min_total_counts)),
             sprintf("# min_novel_length_nt: %d", MIN_NOVEL_LENGTH),
             "# rules applied in order: unstranded, shorter_than_500nt, below_min_total_counts; each transcript is counted under the first rule it fails"),
           con)
utils::write.table(summary_tab, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "-")
close(con)

cat(sprintf("Expression source: %s, min_total_counts = %s\n",
            opt$expression_source, format(opt$min_total_counts)))
cat(sprintf("Transcripts: %d in, %d kept (reference %d, novel %d)\n",
            nrow(tx), sum(keep), sum(keep & !is_novel), sum(keep & is_novel)))
cat(sprintf("Removed: %d unstranded, %d shorter than %d nt, %d below min_total_counts\n",
            sum(unstranded), sum(too_short), MIN_NOVEL_LENGTH, sum(low_expr)))
cat(sprintf("Genes: %d in, %d kept\n", nrow(gene), sum(gene$GENEID %in% kept_genes)))
