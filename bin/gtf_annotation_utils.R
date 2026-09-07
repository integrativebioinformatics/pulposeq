# Shared helpers for reading transcript metadata out of a reference annotation GTF.
#
# Everything needed is already present in the annotation the user supplies with
# --annotation, so it is read locally rather than queried from a server: no network
# dependency, and the metadata is guaranteed to match the annotation the assembly
# actually used.
#
# Both Ensembl and GENCODE GTFs are accepted. They differ in three ways that matter:
#
#   attribute names   Ensembl: gene_biotype / transcript_biotype
#                     GENCODE: gene_type / transcript_type
#   identifiers       Ensembl: unversioned, with separate gene_version /
#                              transcript_version / exon_version attributes
#                     GENCODE: versioned inline (ENSG00000290825.2), no *_version
#   sequence names    Ensembl: 1        GENCODE: chr1
#
# Sourced by ref_transcripts.R, novel_transcripts.R and annotation.R,
# each of which stages this file alongside itself.

suppressPackageStartupMessages({
    library(GenomicRanges)
    library(rtracklayer)
})

#' Remove a trailing version suffix from a stable identifier.
#'
#' ENST00000832824.1 -> ENST00000832824. Ensembl identifiers carry no inline
#' version, so they pass through untouched and this is safe to call without
#' branching on the annotation source.
strip_version <- function(ids) {
    sub("\\.[0-9]+$", "", ids)
}

#' Abort if any identifier repeats.
#'
#' Duplicate transcript or gene IDs silently corrupt the downstream merges by
#' multiplying rows rather than raising an error. The likeliest cause is a
#' GENCODE "ALL" build, whose alternate loci and assembly patches duplicate genes
#' across haplotypes; the CHR or PRI builds do not. Failing loudly here is far
#' better than a report with quietly inflated counts.
assert_unique_ids <- function(ids, what) {
    dups <- unique(ids[duplicated(ids)])
    if (length(dups) > 0) {
        stop(sprintf(
            paste0("Found %d duplicated %s in the reference annotation, e.g. %s.\n",
                   "  This usually means the annotation includes alternate loci or assembly patches\n",
                   "  (the GENCODE 'ALL' build). Use the CHR or PRI build instead."),
            length(dups), what, paste(utils::head(dups, 3), collapse = ", ")
        ), call. = FALSE)
    }
    invisible(TRUE)
}

#' Resolve the biotype column names for either annotation source.
detect_biotype_cols <- function(cols) {
    gene_col <- if ("gene_biotype" %in% cols) {
        "gene_biotype"
    } else if ("gene_type" %in% cols) {
        "gene_type"
    } else {
        NA_character_
    }

    tx_col <- if ("transcript_biotype" %in% cols) {
        "transcript_biotype"
    } else if ("transcript_type" %in% cols) {
        "transcript_type"
    } else {
        NA_character_
    }

    if (is.na(gene_col) || is.na(tx_col)) {
        stop(paste0("Could not find biotype attributes in the reference annotation. ",
                    "Expected gene_biotype/transcript_biotype (Ensembl) or ",
                    "gene_type/transcript_type (GENCODE)."), call. = FALSE)
    }

    list(gene = gene_col, transcript = tx_col)
}

#' Build both the versioned and unversioned form of an identifier.
#'
#' Ensembl supplies the unversioned id plus a separate version attribute; GENCODE
#' bakes the version into the id itself. Returns a list(bare=, versioned=).
build_id_forms <- function(ids, versions) {
    ids <- as.character(ids)

    if (!is.null(versions)) {
        # Ensembl: join the separate version attribute back on
        versions <- as.character(versions)
        versioned <- ifelse(is.na(versions) | !nzchar(versions),
                            ids, paste0(ids, ".", versions))
        list(bare = ids, versioned = versioned)
    } else {
        # GENCODE: the version is already part of the id
        list(bare = strip_version(ids), versioned = ids)
    }
}

#' Read a reference annotation GTF into the shape the pipeline needs.
#'
#' Only transcript and exon features are imported. A full GENCODE annotation is
#' several hundred megabytes and the CDS / UTR / codon records are a large share
#' of that, none of which is used here.
#'
#' Returns a list with:
#'   $tx            one row per transcript
#'   $exons         one row per exon
#'   $gene_biotype  named lookup, gene id (both forms) -> biotype, used to give
#'                  novel transcripts the biotype of the reference gene they
#'                  overlap
read_reference_gtf <- function(path) {
    cat("Reading reference annotation:", path, "\n")

    gtf <- rtracklayer::import(path, feature.type = c("transcript", "exon"))
    cols <- names(mcols(gtf))
    biotype_cols <- detect_biotype_cols(cols)

    cat("  biotype attributes:", biotype_cols$gene, "/", biotype_cols$transcript, "\n")

    has_gene_version <- "gene_version" %in% cols
    has_tx_version   <- "transcript_version" %in% cols
    has_exon_version <- "exon_version" %in% cols
    has_exon_id      <- "exon_id" %in% cols
    cat("  identifier style:", if (has_tx_version) "Ensembl (separate *_version)" else "GENCODE (inline version)", "\n")

    # Sequence names are written exactly as the annotation has them: "1" for
    # Ensembl, "chr1" for GENCODE. No prefix is stripped, for two reasons.
    #
    # Either source may also name a contig with a GRC accession such as
    # KI270728.1 or GL000009.2, where there is no prefix and no safe
    # transformation to apply, so any rule would have to special-case them.
    #
    # And the novel-transcript side reports seqnames straight from the
    # gffcompare GTF, which is never rewritten. Stripping here made the two
    # halves of the report disagree -- the annotated chromosome figure showed
    # "1" while the novel one showed "chr1" for the same GENCODE run.
    norm_seqnames <- function(gr) as.character(seqnames(gr))

    tx_gr   <- gtf[gtf$type == "transcript"]
    exon_gr <- gtf[gtf$type == "exon"]

    if (length(tx_gr) == 0) {
        stop("No transcript features found in the reference annotation.", call. = FALSE)
    }

    tx_ids   <- build_id_forms(tx_gr$transcript_id,
                               if (has_tx_version) tx_gr$transcript_version else NULL)
    gene_ids <- build_id_forms(tx_gr$gene_id,
                               if (has_gene_version) tx_gr$gene_version else NULL)

    assert_unique_ids(tx_ids$versioned, "transcript identifiers")

    # Mature transcript length is the sum of exon widths, which is what biomaRt's
    # transcript_length attribute reports -- not end - start.
    exon_tx_ids <- build_id_forms(exon_gr$transcript_id,
                                  if (has_tx_version) exon_gr$transcript_version else NULL)
    tx_length <- tapply(width(exon_gr), exon_tx_ids$versioned, sum)

    tx <- data.frame(
        chromosome_name               = norm_seqnames(tx_gr),
        ensembl_gene_id               = gene_ids$bare,
        ensembl_gene_id_version       = gene_ids$versioned,
        ensembl_transcript_id         = tx_ids$bare,
        ensembl_transcript_id_version = tx_ids$versioned,
        external_transcript_name      = if ("transcript_name" %in% cols) as.character(tx_gr$transcript_name) else NA_character_,
        external_gene_name            = if ("gene_name" %in% cols) as.character(tx_gr$gene_name) else NA_character_,
        strand                        = as.character(strand(tx_gr)),
        transcript_start              = start(tx_gr),
        transcript_end                = end(tx_gr),
        transcript_length             = as.integer(tx_length[tx_ids$versioned]),
        gene_biotype                  = as.character(mcols(tx_gr)[[biotype_cols$gene]]),
        transcript_biotype            = as.character(mcols(tx_gr)[[biotype_cols$transcript]]),
        stringsAsFactors              = FALSE
    )

    if (has_exon_id) {
        exon_ids <- build_id_forms(exon_gr$exon_id,
                                   if (has_exon_version) exon_gr$exon_version else NULL)$bare
    } else {
        # Fall back to a synthetic id so num_exons still counts correctly.
        exon_ids <- paste0(exon_tx_ids$versioned, ":", start(exon_gr), "-", end(exon_gr))
    }

    exons <- data.frame(
        chromosome_name               = norm_seqnames(exon_gr),
        ensembl_transcript_id         = exon_tx_ids$bare,
        ensembl_transcript_id_version = exon_tx_ids$versioned,
        ensembl_exon_id               = exon_ids,
        exon_chrom_start              = start(exon_gr),
        exon_chrom_end                = end(exon_gr),
        stringsAsFactors              = FALSE
    )

    # Lookups keyed on BOTH id forms, so they succeed whichever form gffcompare
    # reports in ref_id / ref_gene_id. Ensembl writes the version in a separate
    # attribute and GENCODE writes it inline, and a novel-transcript record can end
    # up carrying either.
    dual_key <- function(values, bare, versioned) {
        out <- c(setNames(values, bare), setNames(values, versioned))
        out[!duplicated(names(out))]
    }

    gene_biotype <- dual_key(tx$gene_biotype,
                             tx$ensembl_gene_id, tx$ensembl_gene_id_version)
    gene_name    <- dual_key(tx$external_gene_name,
                             tx$ensembl_gene_id, tx$ensembl_gene_id_version)

    # Transcript level, which the gene level cannot stand in for: a
    # nonsense_mediated_decay or retained_intron isoform of a protein_coding gene
    # has gene_biotype "protein_coding" and transcript_biotype something else
    # entirely, and it is the transcript a novel model was actually compared against.
    tx_biotype <- dual_key(tx$transcript_biotype,
                           tx$ensembl_transcript_id, tx$ensembl_transcript_id_version)
    tx_name    <- dual_key(tx$external_transcript_name,
                           tx$ensembl_transcript_id, tx$ensembl_transcript_id_version)

    # Gene strand, which the sense/antisense split of intronic novel models needs.
    # gffcompare's i code covers both orientations and they are not the same
    # finding: a model inside a same-strand intron cannot be told apart from a
    # fragment of that gene's pre-mRNA without independent evidence, while an
    # antisense one cannot be explained that way at all. Every transcript of a gene
    # shares the gene's strand, so the first occurrence per gene is the gene's.
    gene_strand <- dual_key(tx$strand,
                            tx$ensembl_gene_id, tx$ensembl_gene_id_version)

    cat(sprintf("  %d transcripts, %d exons, %d genes\n",
                nrow(tx), nrow(exons), length(unique(tx$ensembl_gene_id))))

    list(tx = tx, exons = exons,
         gene_biotype = gene_biotype, gene_name = gene_name, gene_strand = gene_strand,
         tx_biotype = tx_biotype, tx_name = tx_name)
}

#' Count distinct exons per transcript, keyed on the versioned transcript id.
exon_counts_per_transcript <- function(exons) {
    counts <- tapply(exons$ensembl_exon_id, exons$ensembl_transcript_id_version,
                     function(x) length(unique(x)))
    data.frame(
        ensembl_transcript_id_version = names(counts),
        num_exons                     = as.integer(counts),
        stringsAsFactors              = FALSE
    )
}

#' Attach metadata columns to a GRanges before export, so the written GTF carries
#' them as attributes. Values are matched by transcript id.
#'
#' `values` is a named list of vectors, each indexed by the transcript ids in `keys`.
annotate_gtf <- function(gr, keys, values) {
    idx <- match(as.character(gr$transcript_id), keys)
    for (nm in names(values)) {
        mcols(gr)[[nm]] <- values[[nm]][idx]
    }
    gr
}
