#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(GenomicRanges)
    library(rtracklayer)
    library(optparse)
})

# Shared reference-GTF helpers, staged alongside this script by the module
source("gtf_annotation_utils.R")

option_list <- list(
    make_option(c("--annotations_gtf"), type="character", default=NULL,
                help="BambuOutput_annotations_validated.gtf", metavar="character"),
    make_option(c("--known_metadata"), type="character", default=NULL,
                help="annotated_transcriptome_metadata.csv", metavar="character"),
    make_option(c("--novel_metadata"), type="character", default=NULL,
                help="novel_transcripts_validated_metadata.csv", metavar="character"),
    make_option(c("--annotation"), type="character", default=NULL,
                help="Path to the reference annotation GTF (Ensembl or GENCODE)", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

required <- c("annotations_gtf",
              "known_metadata", "novel_metadata", "annotation")
missing <- required[vapply(required, function(x) is.null(opt[[x]]), logical(1))]
if (length(missing) > 0) {
    print_help(opt_parser)
    stop(paste("Missing required arguments:", paste(missing, collapse=", ")), call.=FALSE)
}

# --- Build the attribute lookup ------------------------------------------------

#' read.csv that tolerates the empty placeholder files the upstream steps write
#' when a biotype yielded no transcripts.
read_metadata <- function(path, label) {
    tab <- tryCatch(
        read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
        error = function(e) data.frame()
    )
    cat(sprintf("  %s: %d rows\n", label, nrow(tab)))
    tab
}

cat("Reading metadata tables...\n")
known <- read_metadata(opt$known_metadata, "known transcripts")
novel <- read_metadata(opt$novel_metadata, "novel transcripts")

reference <- read_reference_gtf(opt$annotation)
ref_gene_biotype <- reference$gene_biotype
rm(reference)

# The validated GTFs carry Bambu's transcript_id, which for known transcripts is
# whichever id form the reference annotation used. Key the lookup on both forms so
# either matches.
lookup <- list()

if (nrow(known) > 0) {
    known_attrs <- data.frame(
        key                    = known$ensembl_transcript_id_version,
        transcript_status      = "known",
        gene_name              = NA_character_,
        pulposeq_gene_biotype  = NA_character_,
        transcript_name        = NA_character_,
        pulposeq_transcript_biotype = NA_character_,
        class_code             = NA_character_,
        classification         = NA_character_,
        BambuTxClass           = NA_character_,
        BambuNDR               = NA_character_,
        compared_gene_id            = NA_character_,
        compared_gene_name          = NA_character_,
        compared_gene_biotype       = NA_character_,
        compared_transcript_id      = NA_character_,
        compared_transcript_name    = NA_character_,
        compared_transcript_biotype = NA_character_,
        stringsAsFactors       = FALSE
    )
    )
    # duplicate the rows under the unversioned key as well
    known_bare <- known_attrs
    known_bare$key <- known$ensembl_transcript_id
    lookup$known <- rbind(known_attrs, known_bare)
}

if (nrow(novel) > 0) {
    compared_gene <- as.character(novel$compared_gene_id)
    compared_gene[is.na(compared_gene) | compared_gene == "-" |
                      !nzchar(compared_gene)] <- NA_character_

    gene_biotype <- unname(ref_gene_biotype[compared_gene])
    gene_biotype[is.na(gene_biotype)] <- "novel"

    # Taken from the metadata, never re-derived from the prediction. The routing
    # upstream uses the prediction AND the reference gene's biotype together: a
    # non-coding call on a model that shares splice structure with a protein-coding
    # gene becomes novel_non_coding, not novel_lncRNA. An ifelse on the prediction
    # alone cannot represent that distinction and silently collapses the two, which
    # is exactly the bug this replaced. The fallback covers a metadata file written
    # before the column existed, and says so rather than pretending otherwise.
    tx_biotype <- if ("transcript_biotype" %in% names(novel)) {
        as.character(novel$transcript_biotype)
    } else {
        warning("novel metadata has no transcript_biotype column; falling back to ",
                "the coding prediction, which cannot distinguish novel_lncRNA from ",
                "novel_non_coding.")
        ifelse(novel$prediction == "coding", "novel_protein_coding", "novel_lncRNA")
    }

    col <- function(name) if (name %in% names(novel)) as.character(novel[[name]]) else NA_character_

    lookup$novel <- data.frame(
        key                    = novel$qry_id,
        transcript_status      = "novel",
        gene_name              = as.character(novel$gene_name),
        pulposeq_gene_biotype  = gene_biotype,
        transcript_name        = NA_character_,
        pulposeq_transcript_biotype = tx_biotype,
        class_code             = as.character(novel$class_code),
        classification         = as.character(novel$classification),
        BambuTxClass           = col("BambuTxClass"),
        BambuNDR               = col("BambuNDR"),
        compared_gene_id            = compared_gene,
        compared_gene_name          = col("compared_gene_name"),
        compared_gene_biotype       = col("compared_gene_biotype"),
        compared_transcript_id      = col("compared_transcript_id"),
        compared_transcript_name    = col("compared_transcript_name"),
        compared_transcript_biotype = col("compared_transcript_biotype"),
        stringsAsFactors       = FALSE
    )
}

if (length(lookup) == 0) {
    stop("Both metadata tables were empty; nothing to attach.", call.=FALSE)
}

attrs <- do.call(rbind, unname(lookup))
attrs <- attrs[!duplicated(attrs$key), ]
cat(sprintf("Attribute lookup built for %d transcripts\n", nrow(attrs)))

NOVEL_ATTR_COLS <- setdiff(names(attrs), "key")
KNOWN_ATTR_COLS <- "transcript_status"

#' Sort a GTF into a readable order: loci in genomic order, transcripts within a
#' locus, and each transcript's features beneath it.
#'
#' No gene rows are synthesised. Bambu emits transcript and exon only, and a gene
#' row invented here would be a pipeline artefact sitting in a file whose other
#' rows all come from a real annotation.
order_gtf <- function(gr) {
    # The full type list matters for annotations_final.gtf, which carries the
    # reference's CDS, UTR and codon rows. A type missing from this vector matches
    # NA, and NA sorts last, so those rows would be exiled to the end of each
    # transcript instead of sitting with the exons they belong to.
    GTF_TYPE_ORDER <- c("transcript", "exon", "CDS",
                        "five_prime_utr", "three_prime_utr", "UTR",
                        "start_codon", "stop_codon", "Selenocysteine")

    key_gene <- as.character(mcols(gr)$gene_id)
    key_gene[is.na(key_gene)] <- ""

    # Loci in genomic order, from the first base of each gene's earliest record.
    # There is no gene row to read an extent off, so it is derived.
    gene_start <- stats::ave(start(gr), key_gene, FUN = min)

    type_rank <- match(as.character(mcols(gr)$type), GTF_TYPE_ORDER)
    type_rank[is.na(type_rank)] <- length(GTF_TYPE_ORDER) + 1L

    key_tx <- as.character(mcols(gr)$transcript_id)
    key_tx[is.na(key_tx)] <- ""

    gr[order(as.character(seqnames(gr)), gene_start, key_gene,
             key_tx, type_rank, start(gr))]
}


# ---------------------------------------------------------------------------
# annotations_final.gtf -- the plotting annotation
# ---------------------------------------------------------------------------
#
# The same transcripts as the enriched validated GTF, but the KNOWN ones are taken
# from the reference annotation with every feature type intact, rather than from
# Bambu.
#
# Bambu emits transcript and exon rows only. A GTF built from those says where a
# transcript is but not which part of it codes, and plotgardener's plotTranscripts
# derives its thick-versus-thin rendering from the TxDb's CDS: with no CDS records
# every model, coding or not, draws as one uniform box. Reading the known
# transcripts from the reference restores CDS, UTR and codon features, so a
# protein-coding isoform renders with its coding region distinguishable from its
# untranslated ends.
#
# Novel transcripts keep their Bambu structure, because there is nothing else to
# use: Bambu does not call ORFs, so a novel model has no CDS to draw. They render as
# uniform boxes, which is the honest depiction -- the coding-potential prediction is
# a statement about the sequence, not a claim about where a start codon sits.
build_final_gtf <- function() {
    if (is.null(opt$annotation) || !file.exists(opt$annotation)) {
        cat("No reference annotation available; skipping annotations_final.gtf\n")
        return(invisible(NULL))
    }
    if (is.null(opt$annotations_gtf) || !file.exists(opt$annotations_gtf)) {
        cat("No validated GTF available; skipping annotations_final.gtf\n")
    }

        # The input, not an enriched copy: nothing is enriched any more. Status therefore
    # comes from the lookup rather than from an attribute on the file.
    validated <- import(opt$annotations_gtf)
    status    <- attrs$transcript_status[
        match(as.character(mcols(validated)$transcript_id), attrs$key)]

    novel_gr <- validated[!is.na(status) & status == "novel"]
    novel_gr <- annotate_gtf(novel_gr, attrs$key, as.list(attrs[NOVEL_ATTR_COLS]))

    known_id <- unique(as.character(
        mcols(validated)$transcript_id[!is.na(status) & status == "known"]))
    known_id <- known_id[!is.na(known_id)]

    cat(sprintf("Building annotations_final.gtf: %d known transcripts from the reference, %d novel from Bambu\n",
                length(known_id), length(unique(mcols(novel_gr)$transcript_id))))

    # Every feature type, which is the point of reading the reference at all.
    ref_all <- import(opt$annotation)
    ref_tx  <- as.character(mcols(ref_all)$transcript_id)

    # Matched on the bare identifier: Ensembl carries the version in a separate
    # attribute and GENCODE inline, and the validated GTF can name either form.
    keep     <- !is.na(ref_tx) & strip_version(ref_tx) %in% strip_version(known_id)
    ref_keep <- ref_all[keep]

    if (!length(ref_keep)) {
        warning("No reference records matched the validated known transcripts; ",
                "annotations_final.gtf will carry novel transcripts only.")
    } else {
        # transcript_status only. Everything else on a known record is the
        # reference's own, and stays that way.
        ref_keep <- annotate_gtf(ref_keep, attrs$key, as.list(attrs[KNOWN_ATTR_COLS]))
    }

    # The two sources carry different attribute columns, so both are widened to the
    # union before binding: c() on GRanges requires identical mcols. A reference row
    # ends up with NA in the pulposeq_* and compared_* columns and a novel row with
    # NA in level, tag and the rest; export omits NA attributes either way, so each
    # row is written with exactly the attributes that apply to it.
    harmonise <- function(gr, cols) {
        for (nm in cols) {
            if (!nm %in% colnames(mcols(gr))) mcols(gr)[[nm]] <- NA_character_
        }
        mcols(gr) <- mcols(gr)[, cols, drop = FALSE]
        gr
    }

    cols  <- union(colnames(mcols(ref_keep)), colnames(mcols(novel_gr)))
    final <- c(harmonise(ref_keep, cols), harmonise(novel_gr, cols))

    final <- order_gtf(final)
    
    export(final, "annotations_final.gtf")

    cat(sprintf("  wrote annotations_final.gtf (%d records; feature types: %s)\n",
                length(final),
                paste(sort(unique(as.character(mcols(final)$type))), collapse = ", ")))
    invisible(NULL)
}

build_final_gtf()

cat("Final annotation written.\n")
