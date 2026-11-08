#!/usr/bin/env Rscript

# OARFISH_MERGE: combine oarfish's per-sample output into tables shaped like Bambu's.
#
# oarfish runs once per sample. PREFILTER, and everything after it, needs one table
# with a row per transcript and a column per sample, the layout Bambu already
# writes. This script builds that table and, alongside it:
#
#   - gene counts, summed from transcripts with Bambu's own transcript-to-gene map
#   - CPM for every library, and TPM for PacBio only. TPM divides by transcript
#     length, which over-corrects ONT data that already over-represents short
#     transcripts (Paniagua et al. 2026). Both are computed on every model in
#     oarfish's index, so library sizes are not distorted by later filtering.
#   - the bootstrap replicates of every sample, gathered into one R object
#   - a per-sample summary, including one Bambu-oarfish correlation as an alarm
#     for mismatched identifiers or samples
#
# Reads every *.quant, *.infreps.pq and *.meta_info.json file in the working
# directory. Output tables are internal; only the summaries are published.

suppressPackageStartupMessages(library(optparse))

option_list <- list(
    make_option("--bambu_counts", type = "character", default = NULL,
                help = "BambuOutput_counts_transcript.txt: transcript-to-gene map and the correlation check"),
    make_option("--index_summary", type = "character", default = NULL,
                help = "oarfish_index_summary.tsv from OARFISH_INDEX_FASTA"),
    make_option("--library", type = "character", default = NULL,
                help = "ONT_cDNA, ONT_DRS or PacBio; TPM is written for PacBio only")
)
opt <- parse_args(OptionParser(option_list = option_list))
for (arg in c("bambu_counts", "index_summary", "library")) {
    if (is.null(opt[[arg]])) stop("Missing required argument --", arg, call. = FALSE)
}

QUANT_COLS <- c("tname", "len", "num_reads")

# --- Per-sample quantifications --------------------------------------------------

quant_files <- sort(list.files(".", pattern = "\\.quant$"))
if (length(quant_files) == 0) stop("No *.quant files found.", call. = FALSE)
samples <- sub("\\.quant$", "", quant_files)

read_quant <- function(path) {
    q <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
    absent <- setdiff(QUANT_COLS, names(q))
    if (length(absent) > 0) {
        stop(sprintf("%s has no %s column. Found: %s. The oarfish output format may have changed.",
                     path, paste(absent, collapse = "/"), paste(names(q), collapse = ", ")),
             call. = FALSE)
    }
    q
}
quants <- setNames(lapply(quant_files, read_quant), samples)

# Every sample was quantified against the same index, so every file lists the same
# transcripts in the same order. Checked rather than assumed, because the bootstrap
# files carry no transcript names and rely on that order.
tx_names <- quants[[1]]$tname
tx_len   <- quants[[1]]$len
for (s in samples[-1]) {
    if (!identical(quants[[s]]$tname, tx_names)) {
        stop(sprintf("%s.quant lists different transcripts from %s.quant; all samples must share one index.",
                     s, samples[1]), call. = FALSE)
    }
}
counts <- vapply(quants, function(q) as.numeric(q$num_reads), numeric(length(tx_names)))
counts <- matrix(counts, ncol = length(samples), dimnames = list(tx_names, samples))

# --- Transcript-to-gene map, from Bambu -----------------------------------------

bambu <- utils::read.delim(opt$bambu_counts, check.names = FALSE, stringsAsFactors = FALSE)
tx2gene <- setNames(bambu$GENEID, bambu$TXNAME)
unmapped <- setdiff(tx_names, names(tx2gene))
if (length(unmapped) > 0) {
    stop(sprintf("%d oarfish transcripts are not in Bambu's table, e.g. %s. oarfish's index is built from Bambu's models, so every name should match.",
                 length(unmapped), paste(head(unmapped, 3), collapse = ", ")), call. = FALSE)
}
gene_id <- unname(tx2gene[tx_names])

# --- Tables in Bambu's layout ----------------------------------------------------

write_tsv <- function(tab, path) utils::write.table(tab, path, sep = "\t", quote = FALSE, row.names = FALSE)
as_bambu_layout <- function(m) data.frame(TXNAME = tx_names, GENEID = gene_id, m, check.names = FALSE)

write_tsv(as_bambu_layout(counts), "oarfish_counts_transcript.txt")

gene_counts <- rowsum(counts, gene_id, reorder = TRUE)
write_tsv(data.frame(GENEID = rownames(gene_counts), gene_counts, check.names = FALSE),
          "oarfish_counts_gene.txt")

lib_size <- colSums(counts)
cpm <- sweep(counts, 2, ifelse(lib_size > 0, lib_size, 1), "/") * 1e6
write_tsv(as_bambu_layout(cpm), "oarfish_CPM_transcript.txt")

if (opt$library == "PacBio") {
    rpk <- counts / (tx_len / 1000)
    rpk_total <- colSums(rpk)
    tpm <- sweep(rpk, 2, ifelse(rpk_total > 0, rpk_total, 1), "/") * 1e6
    write_tsv(as_bambu_layout(tpm), "oarfish_TPM_transcript.txt")
}

# --- Bootstrap replicates --------------------------------------------------------

infreps_files <- list.files(".", pattern = "\\.infreps\\.pq$")
n_bootstraps <- setNames(rep(0L, length(samples)), samples)
if (length(infreps_files) > 0) {
    if (!requireNamespace("arrow", quietly = TRUE)) {
        stop("Bootstrap files were found but the R package 'arrow' is not installed in this container.",
             call. = FALSE)
    }
    boot <- list()
    for (s in samples) {
        path <- paste0(s, ".infreps.pq")
        if (!file.exists(path)) {
            stop(sprintf("%s has no bootstrap file (%s) although other samples do.", s, path), call. = FALSE)
        }
        m <- as.matrix(as.data.frame(arrow::read_parquet(path)))
        if (nrow(m) != length(tx_names)) {
            stop(sprintf("%s has %d rows but %s.quant has %d transcripts.", path, nrow(m), s, length(tx_names)),
                 call. = FALSE)
        }
        storage.mode(m) <- "double"
        dimnames(m) <- list(tx_names, paste0("rep", seq_len(ncol(m))))
        boot[[s]] <- m
        n_bootstraps[s] <- ncol(m)
    }
    saveRDS(boot, "oarfish_bootstraps.rds")
}

# --- Run information from meta_info.json -----------------------------------------
#
# The JSON layout is oarfish's own and may change between versions, so every
# scalar field is flattened into a long table rather than picking named fields.

meta_files <- list.files(".", pattern = "\\.meta_info\\.json$")
meta_tab <- NULL
if (length(meta_files) > 0 && requireNamespace("jsonlite", quietly = TRUE)) {
    meta_rows <- lapply(meta_files, function(path) {
        fields <- tryCatch(unlist(jsonlite::fromJSON(path, simplifyVector = TRUE)),
                           error = function(e) NULL)
        if (is.null(fields) || length(fields) == 0) return(NULL)
        data.frame(sample = sub("\\.meta_info\\.json$", "", path),
                   field = names(fields), value = unname(as.character(fields)),
                   stringsAsFactors = FALSE)
    })
    meta_tab <- do.call(rbind, meta_rows)
    if (!is.null(meta_tab)) write_tsv(meta_tab, "oarfish_meta_info.tsv")
}

# --- Summary ---------------------------------------------------------------------

index <- utils::read.delim(opt$index_summary, stringsAsFactors = FALSE)
index_n <- function(category) {
    n <- index$transcripts[index$category == category]
    if (length(n) == 1) n else NA
}

# Bambu names its sample columns after the BAM files ("<sample>.bam"); oarfish's
# files are named after the sample alone. Matched without the extension.
bambu_columns <- setdiff(names(bambu), c("TXNAME", "GENEID"))
bambu_mat <- as.matrix(bambu[, bambu_columns, drop = FALSE])
rownames(bambu_mat) <- bambu$TXNAME
bambu_samples <- sub("\\.bam$", "", bambu_columns)
colnames(bambu_mat) <- bambu_samples
shared_tx <- intersect(tx_names, rownames(bambu_mat))

# Spearman correlation over transcripts either tool saw reads for. Leaving out the
# transcripts both called zero keeps the large unexpressed part of the annotation
# from inflating the figure.
spearman_vs_bambu <- vapply(samples, function(s) {
    if (!s %in% bambu_samples) return(NA_real_)
    o <- counts[shared_tx, s]
    b <- bambu_mat[shared_tx, s]
    seen <- o > 0 | b > 0
    if (sum(seen) < 3) return(NA_real_)
    suppressWarnings(stats::cor(o[seen], b[seen], method = "spearman"))
}, numeric(1))

# Read-level mapping outcome, from oarfish's discard table (0.10.x field names).
# NA when a field is absent, e.g. after a change in oarfish's output.
meta_value <- function(sample, field) {
    if (is.null(meta_tab)) return(NA_real_)
    v <- meta_tab$value[meta_tab$sample == sample & meta_tab$field == field]
    if (length(v) == 1) suppressWarnings(as.numeric(v)) else NA_real_
}
reads_mapped   <- vapply(samples, meta_value, numeric(1), field = "discard_table.valid_best_aln")
reads_unmapped <- vapply(samples, meta_value, numeric(1), field = "discard_table.no_mapping")
reads_filtered <- vapply(samples, meta_value, numeric(1), field = "discard_table.no_valid_aln")
mapping_rate   <- reads_mapped / (reads_mapped + reads_unmapped + reads_filtered)

summary_tab <- data.frame(
    sample                    = samples,
    reads_mapped              = reads_mapped,
    reads_unmapped            = reads_unmapped,
    reads_no_valid_alignment  = reads_filtered,
    mapping_rate              = round(mapping_rate, 4),
    assigned_reads            = round(lib_size, 2),
    transcripts_with_reads    = colSums(counts > 0),
    bootstraps                = n_bootstraps[samples],
    spearman_vs_bambu         = round(spearman_vs_bambu, 4),
    index_reference           = index_n("reference"),
    index_novel               = index_n("novel"),
    index_novel_unstranded_left_out = index_n("novel_unstranded_left_out"),
    check.names = FALSE
)
write_tsv(summary_tab, "oarfish_summary.tsv")

cat(sprintf("Merged %d samples x %d transcripts (%d reference, %d novel in the index)\n",
            length(samples), length(tx_names), index_n("reference"), index_n("novel")))
cat(sprintf("Bootstraps per sample: %s\n", paste(unique(n_bootstraps), collapse = ", ")))
cat(sprintf("Spearman vs Bambu: %s\n", paste(sprintf("%s=%.3f", samples, spearman_vs_bambu), collapse = ", ")))
low <- samples[!is.na(spearman_vs_bambu) & spearman_vs_bambu < 0.3]
if (length(low) > 0) {
    warning(sprintf("Low agreement with Bambu (Spearman < 0.3) for: %s. Check that sample names and transcript IDs match.",
                    paste(low, collapse = ", ")), call. = FALSE)
}
