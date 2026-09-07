#!/usr/bin/env Rscript

# Attach biotype and classification attributes to the validated GTFs produced by
# filter_bambu_gtf.sh.
#
# That script filters the Bambu extended annotation down to the transcripts that
# survived count validation, using awk. It is deliberately left untouched: Bambu's
# GTF carries only gene_id, transcript_id and exon_number, so the metadata has to
# come from elsewhere and adding it in awk would be fragile. This script
# post-processes the three outputs instead.
#
# Each validated GTF holds a mixture of known and novel transcripts, so both
# lookups are applied:
#   known -> the reference annotation, via annotated_transcriptome_metadata.csv
#   novel -> the gffcompare/RNAmining results, via novel_transcripts_validated_metadata.csv

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
    make_option(c("--fulllength_gtf"), type="character", default=NULL,
                help="BambuOutput_fullLength_validated.gtf", metavar="character"),
    make_option(c("--unique_gtf"), type="character", default=NULL,
                help="BambuOutput_uniquelyMapped_validated.gtf", metavar="character"),
    make_option(c("--known_metadata"), type="character", default=NULL,
                help="annotated_transcriptome_metadata.csv", metavar="character"),
    make_option(c("--novel_metadata"), type="character", default=NULL,
                help="novel_transcripts_validated_metadata.csv", metavar="character"),
    make_option(c("--annotation"), type="character", default=NULL,
                help="Path to the reference annotation GTF (Ensembl or GENCODE)", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

required <- c("annotations_gtf", "fulllength_gtf", "unique_gtf",
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
        key                = known$ensembl_transcript_id_version,
        transcript_status  = "known",
        gene_name          = known$external_gene_name,
        gene_biotype       = known$gene_biotype,
        transcript_name    = known$external_transcript_name,
        transcript_biotype = known$transcript_biotype,
        class_code             = NA_character_,
        classification         = NA_character_,
        # Bambu's txClassDescription is "annotation" for every reference transcript,
        # which says nothing. NA, like class_code.
        BambuTxClass           = NA_character_,
        BambuNDR               = NA_character_,
        # A known transcript IS the reference, so it has nothing to have been
        # compared against. NA rather than self-reference, matching class_code.
        ref_gene_id            = NA_character_,
        ref_gene_name          = NA_character_,
        ref_gene_biotype       = NA_character_,
        ref_transcript_id      = NA_character_,
        ref_transcript_name    = NA_character_,
        ref_transcript_biotype = NA_character_,
        stringsAsFactors       = FALSE
    )
    # duplicate the rows under the unversioned key as well
    known_bare <- known_attrs
    known_bare$key <- known$ensembl_transcript_id
    lookup$known <- rbind(known_attrs, known_bare)
}

if (nrow(novel) > 0) {
    ref_gene <- as.character(novel$ref_gene_id)
    ref_gene[is.na(ref_gene) | ref_gene == "-" | !nzchar(ref_gene)] <- NA_character_

    gene_biotype <- unname(ref_gene_biotype[ref_gene])
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
        gene_biotype           = gene_biotype,
        transcript_name        = NA_character_,
        transcript_biotype     = tx_biotype,
        class_code             = as.character(novel$class_code),
        classification         = as.character(novel$classification),
        BambuTxClass           = col("BambuTxClass"),
        BambuNDR               = col("BambuNDR"),
        ref_gene_id            = ref_gene,
        ref_gene_name          = col("ref_gene_name"),
        ref_gene_biotype       = col("ref_gene_biotype"),
        ref_transcript_id      = col("ref_id"),
        ref_transcript_name    = col("ref_transcript_name"),
        ref_transcript_biotype = col("ref_transcript_biotype"),
        stringsAsFactors       = FALSE
    )
}

if (length(lookup) == 0) {
    stop("Both metadata tables were empty; nothing to attach.", call.=FALSE)
}

attrs <- do.call(rbind, unname(lookup))
attrs <- attrs[!duplicated(attrs$key), ]
cat(sprintf("Attribute lookup built for %d transcripts\n", nrow(attrs)))

attr_cols <- setdiff(names(attrs), "key")

# --- Apply to each validated GTF ----------------------------------------------

#' Synthesise one `gene` feature per gene_id, spanning all of its records.
#'
#' Bambu emits only transcript and exon rows, and validate_bambu_gtf.sh drops any
#' line without a transcript_id, so the validated GTFs reach here with no gene
#' features at all. Everything in this pipeline derives gene extent from the
#' transcript rows and does not need them, but tools fed the published GTF do:
#' IGV cannot collapse a locus without a gene row, and several browsers and
#' converters assume the canonical gene/transcript/exon hierarchy.
add_gene_features <- function(gr) {
    gene_ids <- as.character(mcols(gr)$gene_id)
    ok <- !is.na(gene_ids) & nzchar(gene_ids)
    if (!any(ok)) {
        warning("No gene_id values found; no gene features were added.")
        return(gr)
    }

    sub <- gr[ok]
    gid <- gene_ids[ok]

    # range() over a split keeps seqname and strand with each group, so a gene whose
    # records somehow land on two contigs yields one row per fragment rather than a
    # single span across both. That should not happen, so say so if it does.
    spans <- unlist(range(split(sub, gid)), use.names = TRUE)
    if (anyDuplicated(names(spans))) {
        dups <- unique(names(spans)[duplicated(names(spans))])
        warning(sprintf("%d gene(s) span more than one sequence or strand; one gene row emitted per fragment.",
                        length(dups)))
    }

    # Gene-level values are constant within a gene: annotate_gtf keys on
    # transcript_id, and every transcript of a gene resolves to the same gene
    # metadata. So the first matching record is representative.
    md <- mcols(sub)[match(names(spans), gid), , drop = FALSE]

    # Blank the isoform-level fields. Carrying transcript_id or class_code onto a
    # gene row would assert that the locus has one isoform's identity.
    # `x[] <- NA` rather than `x <- NA`: the former keeps the column's type, the
    # latter would turn a character column logical and then c() below would refuse
    # to bind it against the original.
    for (nm in intersect(c("transcript_id", "transcript_name", "transcript_biotype",
                           "transcript_status", "exon_number", "exon_id",
                           "class_code", "classification", "score", "phase"),
                         colnames(md))) {
        md[[nm]][] <- NA
    }
    md$gene_id <- names(spans)

    # import() returns source and type as factors, neither carrying a "gene" or
    # "pulposeq" level -- assigning to them directly would silently yield NA and
    # export a typeless row. Character on both sides sidesteps the level juggling,
    # and export() coerces either way.
    mcols(gr)$type   <- as.character(mcols(gr)$type)
    mcols(gr)$source <- as.character(mcols(gr)$source)
    md$type   <- "gene"
    md$source <- "pulposeq"

    mcols(spans) <- md
    names(spans) <- NULL

    combined <- c(spans, gr)

    # Canonical GTF ordering:
    #
    #   gene
    #     transcript A
    #       exons and CDS of A
    #     transcript B
    #       exons and CDS of B
    #
    # Sorting on type before anything transcript-specific does NOT give this. It
    # gives the gene row, then every transcript of the gene, then every exon of
    # every transcript -- valid GTF, but a reader cannot walk one isoform without
    # scanning the whole locus, and the exons of A sit nowhere near A.
    #
    # bambu.R's reorder_gtf() solves the same problem with
    # order(transcript_id, type == "transcript", decreasing = TRUE), which groups
    # features under their transcript. That works on Bambu's output because it has
    # no gene rows; here it would strand them, since a gene row has no
    # transcript_id to sort on. Loci are keyed on gene_id instead, which every row
    # carries, and transcripts run alphabetically within their gene.
    #
    # The full type list matters for annotations_final.gtf, which carries the
    # reference's CDS, UTR and codon rows. A type missing from this vector matches
    # NA, and NA sorts last, so those rows would be exiled to the end of each
    # transcript instead of sitting with the exons they belong to.
    GTF_TYPE_ORDER <- c("gene", "transcript", "exon", "CDS",
                        "five_prime_utr", "three_prime_utr", "UTR",
                        "start_codon", "stop_codon", "Selenocysteine")

    # gene_id is the anchor, because it is the one identifier every row carries.
    # transcript_id cannot do this job: gene rows have none, so sorting on it strands
    # them at the end of the file, away from the gene they describe.
    key_gene   <- as.character(mcols(combined)$gene_id)
    gene_start <- start(spans)[match(key_gene, mcols(spans)$gene_id)]

    type_rank  <- match(as.character(mcols(combined)$type), GTF_TYPE_ORDER)
    # Anything unlisted keeps a stable place after the known types rather than
    # being scattered by NA.
    type_rank[is.na(type_rank)] <- length(GTF_TYPE_ORDER) + 1L

    # The empty string is what leads each locus with its gene row: it sorts before
    # any real identifier, and it keeps NA out of order(), which given several keys
    # can move any element carrying a missing value to the end whatever the earlier
    # keys say.
    key_tx <- as.character(mcols(combined)$transcript_id)
    key_tx[is.na(key_tx)] <- ""

    combined <- combined[order(as.character(seqnames(combined)),
                               gene_start, key_gene,
                               key_tx, type_rank, start(combined))]

    cat(sprintf("  added %d gene features\n", length(spans)))
    combined
}

enrich <- function(path) {
    if (is.null(path) || !file.exists(path)) {
        cat("Skipping missing file:", path, "\n")
        return(invisible(NULL))
    }

    gr <- import(path)
    if (length(gr) == 0) {
        cat("Empty GTF, re-exporting unchanged:", basename(path), "\n")
        export(gr, basename(path))
        return(invisible(NULL))
    }

    gr <- annotate_gtf(gr, attrs$key, as.list(attrs[attr_cols]))

    matched <- sum(!is.na(gr$transcript_status))
    cat(sprintf("%s: %d of %d records annotated\n",
                basename(path), matched, length(gr)))

    if (matched == 0) {
        warning(sprintf(paste0("No records in %s matched the metadata tables. ",
                               "Check that the transcript identifiers agree."),
                        basename(path)))
    }

    gr <- add_gene_features(gr)

    export(gr, basename(path))
    invisible(NULL)
}

cat("Enriching validated GTFs...\n")
enrich(opt$annotations_gtf)
enrich(opt$fulllength_gtf)
enrich(opt$unique_gtf)

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
    if (is.null(opt$annotations_gtf) || !file.exists(basename(opt$annotations_gtf))) {
        cat("No enriched validated GTF available; skipping annotations_final.gtf\n")
        return(invisible(NULL))
    }

    # The enriched copy written above, not the input, so the novel rows already
    # carry the pipeline's attributes.
    validated <- import(basename(opt$annotations_gtf))
    status    <- as.character(mcols(validated)$transcript_status)

    novel_gr <- validated[!is.na(status) & status == "novel"]
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
        ref_keep <- annotate_gtf(ref_keep, attrs$key, as.list(attrs[attr_cols]))
    }

    # The two sources carry different attribute columns, so both are widened to the
    # union before binding: c() on GRanges requires identical mcols.
    harmonise <- function(gr, cols) {
        for (nm in cols) {
            if (!nm %in% colnames(mcols(gr))) mcols(gr)[[nm]] <- NA_character_
        }
        mcols(gr) <- mcols(gr)[, cols, drop = FALSE]
        gr
    }

    cols  <- union(colnames(mcols(ref_keep)), colnames(mcols(novel_gr)))
    final <- c(harmonise(ref_keep, cols), harmonise(novel_gr, cols))

    final <- add_gene_features(final)
    export(final, "annotations_final.gtf")

    cat(sprintf("  wrote annotations_final.gtf (%d records; feature types: %s)\n",
                length(final),
                paste(sort(unique(as.character(mcols(final)$type))), collapse = ", ")))
    invisible(NULL)
}

build_final_gtf()

cat("Validated GTF enrichment completed successfully!\n")
