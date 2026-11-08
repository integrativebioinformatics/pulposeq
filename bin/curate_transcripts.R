#!/usr/bin/env Rscript

# CURATE: write the final, curated transcript set.
#
# After PREFILTER and novel-transcript curation, one list of transcripts remains:
# every reference transcript PREFILTER kept, plus the novel models that passed
# class-code and biotype curation. Every table and annotation written here
# contains exactly that list, so Bambu's and oarfish's results line up row for row.
#
# Replaces validate_bambu_counts.sh and validate_bambu_gtf.sh/.awk. Their outputs
# keep the same names and the same rule:
#   - transcript tables: reference rows as PREFILTER kept them, plus novel rows
#     listed in the curated novel metadata
#   - gene table: genes with at least one curated transcript
#   - GTFs: the curated transcripts, and the subsets with at least one full-length
#     or unique read
# Tables and GTFs are subset as text, line by line, so every kept line is written
# exactly as it arrived.
#
# New here:
#   - transcript_read_support.tsv: per transcript, total Bambu reads, full-length
#     reads and unique reads, whether any read covers the full structure, and the
#     oarfish total when oarfish ran (plus its mean InfRV when it ran with bootstraps)
#   - oarfish's tables restricted to the same transcripts, and a
#     SummarizedExperiment with the bootstrap replicates (assays infRep1..N, the
#     layout fishpond expects)

suppressPackageStartupMessages(library(optparse))

option_list <- list(
    make_option("--novel_metadata", type = "character", default = NULL,
                help = "novel_transcripts_validated_metadata.csv (qry_id column)"),
    make_option("--counts_transcript", type = "character", default = NULL),
    make_option("--counts_gene", type = "character", default = NULL),
    make_option("--cpm_transcript", type = "character", default = NULL),
    make_option("--full_length", type = "character", default = NULL),
    make_option("--unique", type = "character", default = NULL),
    make_option("--extended_gtf", type = "character", default = NULL,
                help = "BambuOutput_extended_annotations.gtf"),
    make_option("--expression_source", type = "character", default = "bambu"),
    make_option("--oarfish_counts", type = "character", default = NULL),
    make_option("--oarfish_counts_gene", type = "character", default = NULL),
    make_option("--oarfish_cpm", type = "character", default = NULL),
    make_option("--oarfish_tpm", type = "character", default = NULL),
    make_option("--oarfish_bootstraps", type = "character", default = NULL)
)
opt <- parse_args(OptionParser(option_list = option_list))
required <- c("novel_metadata", "counts_transcript", "counts_gene", "cpm_transcript",
              "full_length", "unique", "extended_gtf")
missing_args <- required[vapply(required, function(x) is.null(opt[[x]]), logical(1))]
if (length(missing_args) > 0) {
    stop("Missing required arguments: ", paste(missing_args, collapse = ", "), call. = FALSE)
}
has <- function(x) !is.null(x) && nzchar(x) && file.exists(x)

NOVEL_TX_PREFIX <- "BambuTx"

# --- Line-level helpers ----------------------------------------------------------

first_field  <- function(lines) sub("\t.*$", "", lines)
second_field <- function(lines) sub("^[^\t]*\t([^\t]*).*$", "\\1", lines)

#' Keep the header plus the body lines whose first field passes `keep`.
subset_table_lines <- function(path, keep) {
    lines <- readLines(path)
    if (length(lines) == 0) return(lines)
    body <- lines[-1]
    c(lines[1], body[keep(first_field(body))])
}

# --- The curated transcript list -------------------------------------------------

novel_meta <- utils::read.csv(opt$novel_metadata, stringsAsFactors = FALSE, check.names = FALSE)
if (!"qry_id" %in% names(novel_meta)) {
    stop(basename(opt$novel_metadata), " has no qry_id column.", call. = FALSE)
}
curated_novel <- unique(as.character(novel_meta$qry_id))

keep_tx <- function(ids) !startsWith(ids, NOVEL_TX_PREFIX) | ids %in% curated_novel

out_tx <- c(
    counts_transcript = "BambuOutput_counts_transcript_validated.txt",
    cpm_transcript    = "BambuOutput_CPM_transcript_validated.txt",
    full_length       = "BambuOutput_fullLengthCounts_transcript_validated.txt",
    unique            = "BambuOutput_uniqueCounts_transcript_validated.txt"
)
validated <- list()
for (key in names(out_tx)) {
    validated[[key]] <- subset_table_lines(opt[[key]], keep_tx)
    writeLines(validated[[key]], out_tx[[key]])
}
ids_of <- function(lines) if (length(lines) > 1) first_field(lines[-1]) else character(0)

curated_tx    <- ids_of(validated$counts_transcript)
curated_genes <- unique(second_field(validated$counts_transcript[-1]))
gene_lines    <- subset_table_lines(opt$counts_gene, function(ids) ids %in% curated_genes)
writeLines(gene_lines, "BambuOutput_counts_gene_validated.txt")

# --- Curated GTFs ----------------------------------------------------------------
#
# The transcript identifier is read the way validate_bambu_gtf.awk read it: the
# text after the last "transcript_id", without an opening quote, up to the first
# quote, semicolon or whitespace. Comment lines and lines without a transcript_id
# are dropped.

gtf_lines <- readLines(opt$extended_gtf)
gtf_lines <- gtf_lines[!startsWith(gtf_lines, "#") & grepl("transcript_id", gtf_lines, fixed = TRUE)]
gtf_tx <- sub(".*transcript_id[ \t]*", "", gtf_lines)
gtf_tx <- sub("^[\"']", "", gtf_tx)
gtf_tx <- sub("[\"'; \t].*$", "", gtf_tx)

write_gtf_subset <- function(ids, path) {
    writeLines(gtf_lines[nzchar(gtf_tx) & gtf_tx %in% ids], path)
}
write_gtf_subset(curated_tx, "BambuOutput_annotations_validated.gtf")
write_gtf_subset(ids_of(validated$full_length), "BambuOutput_fullLength_validated.gtf")
write_gtf_subset(ids_of(validated$unique), "BambuOutput_uniquelyMapped_validated.gtf")

# --- Read support ----------------------------------------------------------------

read_numeric <- function(path) utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
row_totals <- function(tab) {
    m <- as.matrix(tab[, setdiff(names(tab), c("TXNAME", "GENEID")), drop = FALSE])
    setNames(rowSums(m), tab$TXNAME)
}
lookup <- function(totals, ids) {
    v <- unname(totals[ids])
    v[is.na(v)] <- 0
    v
}

counts_tab <- read_numeric(out_tx[["counts_transcript"]])
fl_total   <- row_totals(read_numeric(out_tx[["full_length"]]))
uq_total   <- row_totals(read_numeric(out_tx[["unique"]]))

support <- data.frame(
    TXNAME                  = counts_tab$TXNAME,
    GENEID                  = counts_tab$GENEID,
    transcript_status       = ifelse(startsWith(counts_tab$TXNAME, NOVEL_TX_PREFIX), "novel", "reference"),
    bambu_counts_total      = unname(row_totals(counts_tab)),
    full_length_reads_total = lookup(fl_total, counts_tab$TXNAME),
    unique_reads_total      = lookup(uq_total, counts_tab$TXNAME),
    stringsAsFactors        = FALSE
)
support$full_length_supported <- support$full_length_reads_total > 0

# --- oarfish, restricted to the same transcripts ---------------------------------

if (has(opt$oarfish_counts)) {
    oar_lines <- subset_table_lines(opt$oarfish_counts, function(ids) ids %in% curated_tx)
    missing_oar <- setdiff(curated_tx, ids_of(oar_lines))
    if (length(missing_oar) > 0) {
        stop(sprintf("%d curated transcripts have no oarfish counts, e.g. %s.", length(missing_oar),
                     paste(head(missing_oar, 3), collapse = ", ")), call. = FALSE)
    }
    writeLines(oar_lines, "oarfish_counts_transcript_validated.txt")
    if (has(opt$oarfish_cpm)) {
        writeLines(subset_table_lines(opt$oarfish_cpm, function(ids) ids %in% curated_tx),
                   "oarfish_CPM_transcript_validated.txt")
    }
    if (has(opt$oarfish_tpm)) {
        writeLines(subset_table_lines(opt$oarfish_tpm, function(ids) ids %in% curated_tx),
                   "oarfish_TPM_transcript_validated.txt")
    }
    if (has(opt$oarfish_counts_gene)) {
        writeLines(subset_table_lines(opt$oarfish_counts_gene, function(ids) ids %in% curated_genes),
                   "oarfish_counts_gene_validated.txt")
    }

    oar_tab <- read_numeric("oarfish_counts_transcript_validated.txt")
    oar_tab <- oar_tab[match(curated_tx, oar_tab$TXNAME), , drop = FALSE]
    support$oarfish_counts_total <- unname(row_totals(oar_tab))

    # SummarizedExperiment: counts plus one assay per bootstrap replicate
    samples <- setdiff(names(oar_tab), c("TXNAME", "GENEID"))
    counts  <- as.matrix(oar_tab[, samples, drop = FALSE])
    rownames(counts) <- oar_tab$TXNAME
    assays <- list(counts = counts)
    if (has(opt$oarfish_bootstraps)) {
        boot <- readRDS(opt$oarfish_bootstraps)
        n_rep <- unique(vapply(boot, ncol, integer(1)))
        if (length(n_rep) != 1) stop("Samples have different numbers of bootstrap replicates.", call. = FALSE)
        for (k in seq_len(n_rep)) {
            assays[[paste0("infRep", k)]] <- vapply(samples, function(s) boot[[s]][curated_tx, k],
                                                     numeric(length(curated_tx)))
            dimnames(assays[[paste0("infRep", k)]]) <- dimnames(counts)
        }
        # Inferential relative variance (Zhu et al. 2019), as fishpond::computeInfRV
        # computes it: per sample, max(var - mean, 0) / (mean + 5) + 0.01 over the
        # bootstrap replicates, then the mean across samples. It needs at least two
        # replicates for a variance.
        if (n_rep >= 2) {
            infrv <- vapply(samples, function(s) {
                b  <- boot[[s]][curated_tx, , drop = FALSE]
                mu <- rowMeans(b)
                v  <- rowSums((b - mu)^2) / (ncol(b) - 1)
                pmax(v - mu, 0) / (mu + 5) + 0.01
            }, numeric(length(curated_tx)))
            support$oarfish_mean_infrv <- unname(rowMeans(matrix(infrv, nrow = length(curated_tx))))
        }
    }
    se <- SummarizedExperiment::SummarizedExperiment(
        assays  = assays,
        rowData = S4Vectors::DataFrame(GENEID = oar_tab$GENEID,
                                       transcript_status = support$transcript_status,
                                       row.names = oar_tab$TXNAME),
        colData = S4Vectors::DataFrame(sample = samples, row.names = samples)
    )
    S4Vectors::metadata(se)$expression_source <- opt$expression_source
    saveRDS(se, "oarfish_validated_se.rds")
}

utils::write.table(support, "transcript_read_support.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("Curated transcripts: %d (reference %d, novel %d); genes: %d\n",
            nrow(support), sum(support$transcript_status == "reference"),
            sum(support$transcript_status == "novel"), length(gene_lines) - 1))
cat(sprintf("With at least one full-length read: %d of %d\n",
            sum(support$full_length_supported), nrow(support)))
if (has(opt$oarfish_counts)) cat("oarfish tables restricted to the same transcripts.\n")
