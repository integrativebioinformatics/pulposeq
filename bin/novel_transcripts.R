#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(rtracklayer)
    library(GenomicRanges)
    library(optparse)
})

# Shared reference-GTF helpers, staged alongside this script by the module
source("gtf_annotation_utils.R")

# Define command line options
option_list <- list(
    make_option(c("--bambu_gtf"), type="character", default=NULL,
                help="Path to bambu novel transcripts GTF file", metavar="character"),
    make_option(c("--annotation"), type="character", default=NULL,
                help="Path to the reference annotation GTF (Ensembl or GENCODE)", metavar="character"),
    make_option(c("--compared_gtf"), type="character", default=NULL,
                help="Path to compared transcriptome annotated GTF file", metavar="character"),
    make_option(c("--tmap_file"), type="character", default=NULL,
                help="Path to tmap results file", metavar="character"),
    make_option(c("--coding_predictions"), type="character", default=NULL,
                help="Path to the coding-potential table (CPC2 or RNAmining)", metavar="character"),
    make_option(c("--tx_counts"), type="character", default=NULL,
                help="Path to transcript counts file", metavar="character"),
    make_option(c("--gene_counts"), type="character", default=NULL,
                help="Path to gene counts file", metavar="character"),
    make_option(c("--tx_classes"), type="character", default=NULL,
                help="bambu_novel_tx_classes.csv: txClassDescription and NDR per novel transcript",
                metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check if all required arguments are provided
if (is.null(opt$bambu_gtf) || is.null(opt$compared_gtf) || is.null(opt$tmap_file) ||
    is.null(opt$coding_predictions) || is.null(opt$tx_counts) || is.null(opt$gene_counts) ||
    is.null(opt$annotation)) {
    print_help(opt_parser)
    stop("All input files must be specified.", call.=FALSE)
}

# Reference lookups, used to say what each novel transcript was compared against:
# the biotype, name and strand of the reference gene, and the biotype and name of
# the reference transcript itself. Kept as named vectors rather than the whole
# table, which is several hundred megabytes for a full GENCODE annotation.
reference <- read_reference_gtf(opt$annotation)
ref_gene_biotype <- reference$gene_biotype
ref_gene_name    <- reference$gene_name
ref_gene_strand  <- reference$gene_strand
ref_tx_biotype   <- reference$tx_biotype
ref_tx_name      <- reference$tx_name
rm(reference)

# Import GTF files
cat("Loading GTF files...\n")
gtf <- import(opt$bambu_gtf)
tx_gtf <- import(opt$compared_gtf)

tx_table <- as.data.frame(tx_gtf[tx_gtf$type == "transcript"], )

# Load tmap results
cat("Loading tmap results...\n")
tmap <- read_table(opt$tmap_file)

cat("Loading coding-potential predictions...\n")

#' Read a coding-potential table from either predictor.
#'
#' The two formats differ in more than layout. RNAmining writes a preamble and then
#' <id> <coding|non-coding> <score>, where the score is the probability of whichever
#' class it chose. CPC2 writes a "#ID" header with eight tab-separated columns,
#' spells the label "noncoding" without the hyphen, and its coding_probability is
#' the probability of CODING specifically -- so a non-coding row carries a low value
#' where RNAmining carries a high one. Passing either straight through would put two
#' different quantities in the same column.
#'
#' Both are normalised here: prediction in {coding, non-coding}, and coding_prob as
#' P(coding) whichever tool produced it. CPC2 columns are located by header name
#' rather than position, because it inserts ORF_Start when run with --ORF.
read_coding_predictions <- function(path) {
    lines <- readLines(path)
    hdr   <- grep("^#ID\t", lines)

    if (length(hdr)) {
        cols <- strsplit(sub("^#", "", lines[hdr[1]]), "\t", fixed = TRUE)[[1]]
        body <- lines[seq.int(hdr[1] + 1L, length(lines))]
        body <- body[nzchar(body)]
        if (!length(body)) {
            stop("CPC2 table ", path, " has a header but no rows.", call. = FALSE)
        }

        d <- read.table(text = body, sep = "\t", header = FALSE, quote = "",
                        comment.char = "", stringsAsFactors = FALSE)
        if (ncol(d) != length(cols)) {
            stop("CPC2 table ", path, " has ", ncol(d), " columns against a ",
                 length(cols), "-column header.", call. = FALSE)
        }
        names(d) <- cols

        missing <- setdiff(c("ID", "label", "coding_probability"), names(d))
        if (length(missing)) {
            stop("CPC2 table ", path, " is missing column(s): ",
                 paste(missing, collapse = ", "), ". Found: ",
                 paste(names(d), collapse = ", "), call. = FALSE)
        }

        cat(sprintf("Read %d CPC2 predictions\n", nrow(d)))
        return(data.frame(
            transcript_id    = as.character(d$ID),
            prediction       = ifelse(d$label == "coding", "coding", "non-coding"),
            coding_prob      = as.numeric(d$coding_probability),
            coding_predictor = "cpc2",
            stringsAsFactors = FALSE))
    }

    # RNAmining. The preamble line count has changed between versions, and skipping a
    # fixed number silently drops one transcript's call whenever that count is wrong,
    # so the body is identified by matching it instead.
    is_data <- grepl("^[^#[:space:]][^\t]*\t(coding|non-coding)\t", lines)
    if (!any(is_data)) {
        stop("No prediction rows found in ", path, ". Expected either a CPC2 table ",
             "with a #ID header, or RNAmining <id> <coding|non-coding> <score> lines.",
             call. = FALSE)
    }

    d <- read.table(text = lines[is_data], header = FALSE, sep = "\t",
                    stringsAsFactors = FALSE)
    colnames(d) <- c("transcript_id", "prediction", "score")

    cat(sprintf("Read %d RNAmining predictions (%d preamble lines skipped)\n",
                nrow(d), sum(!is_data)))
    data.frame(
        transcript_id    = as.character(d$transcript_id),
        prediction       = as.character(d$prediction),
        coding_prob      = ifelse(d$prediction == "coding", d$score, 1 - d$score),
        coding_predictor = "rnamining",
        stringsAsFactors = FALSE)
}

rnam <- read_coding_predictions(opt$coding_predictions)

# Select relevant information from gtf
tx_table <- dplyr::select(tx_table, seqnames, transcript_id, start, end, strand)

# Complete info dataframe
cat("Merging data...\n")
tx_info <- merge(tmap, tx_table, by.x="qry_id", by.y="transcript_id", all.x=TRUE)

# Attach the coding-potential calls
tx_info <- merge(tx_info, rnam, by.x="qry_id", by.y="transcript_id", all.x=TRUE)
tx_info$qry_gene_name    <- unname(ref_gene_name[tx_info$qry_gene_id])
tx_info$qry_gene_biotype <- unname(ref_gene_biotype[tx_info$qry_gene_id])
tx_info <- dplyr::select(tx_info, seqnames, qry_id, ref_id,
                  qry_gene_id, qry_gene_name, qry_gene_biotype, ref_gene_id,
                  class_code, strand, start, end, len, num_exons, prediction, coding_prob,
                  coding_predictor)

# Remove unstranded transcripts
tx_info <- tx_info[tx_info$strand != "*", ]

tx_info$exon <- ifelse(tx_info$num_exons == 1, 'mono-exonic', 'multi-exonic')

#' Reference gene id with gffcompare's placeholders normalised to NA. "-" means
#' no reference match, which is a different statement from a missing value.
normalise_ref_gene <- function(ids) {
    ids <- as.character(ids)
    ids[is.na(ids) | ids == "-" | !nzchar(ids)] <- NA_character_
    ids
}

# What the novel transcript was compared against, recorded on every novel record
# rather than only in the GTF.

tx_info$ref_gene_id <- normalise_ref_gene(tx_info$ref_gene_id)
tx_info$ref_id      <- normalise_ref_gene(tx_info$ref_id)

tx_info$ref_gene_biotype <- unname(ref_gene_biotype[tx_info$ref_gene_id])
tx_info$ref_gene_name    <- unname(ref_gene_name[tx_info$ref_gene_id])

# Transcript level as well as gene level, because the two disagree in exactly the
# cases worth looking at. A y-class novel transcript containing TARDBP-221 has
# ref_gene_biotype "protein_coding" and ref_transcript_biotype
# "nonsense_mediated_decay", and only the second says what it actually contains.
tx_info$ref_transcript_biotype <- unname(ref_tx_biotype[tx_info$ref_id])
tx_info$ref_transcript_name    <- unname(ref_tx_name[tx_info$ref_id])

# Orientation relative to the reference gene. NA where gffcompare reported no
# reference, which is the honest value for u: there is nothing to be sense or
# antisense to. x is antisense by construction and should always come out FALSE,
# which makes it a free check that the lookup is keyed correctly.
tx_info$ref_gene_strand    <- unname(ref_gene_strand[tx_info$ref_gene_id])
tx_info$same_strand_as_ref <- ifelse(is.na(tx_info$ref_gene_strand), NA,
                                     tx_info$strand == tx_info$ref_gene_strand)
# --- Naming: what the pipeline/Bambu decided vs what it was compared against ---------
#
names(tx_info)[names(tx_info) == "ref_gene_id"]            <- "compared_gene_id"
names(tx_info)[names(tx_info) == "ref_id"]                 <- "compared_transcript_id"
names(tx_info)[names(tx_info) == "ref_gene_name"]          <- "compared_gene_name"
names(tx_info)[names(tx_info) == "ref_transcript_name"]    <- "compared_transcript_name"
names(tx_info)[names(tx_info) == "ref_gene_biotype"]       <- "compared_gene_biotype"
names(tx_info)[names(tx_info) == "ref_transcript_biotype"] <- "compared_transcript_biotype"
names(tx_info)[names(tx_info) == "ref_gene_strand"]        <- "compared_gene_strand"
names(tx_info)[names(tx_info) == "same_strand_as_ref"]     <- "same_strand_as_compared"

# --- Bambu's own view of why each model is novel -------------------------------
#
# gffcompare and Bambu answer different questions about the same transcript, and
# neither subsumes the other. A class code says how the model sits relative to the
# reference it was matched against -- inside an intron, antisense, sharing a
# junction. txClassDescription says which PART of the model Bambu had to invent to
# build it: a new first exon, a new junction, a whole new gene. A `j` model whose
# only novelty is `newLastExon` is a different claim from a `j` model that is
# `allNew`, and the class code alone cannot separate them.
#
# One or more classes are packed into a colon-separated string, which is kept whole
# here. The report splits it into sets; splitting it into columns at this point
# would fix the vocabulary to whatever this Bambu version emits.
#
# BambuNDR is the novel discovery rate Bambu assigned to this transcript, not the
# run-level threshold: a per-model score, where the run parameter is the cutoff it
# was tested against. Carried per transcript so a candidate can be read against the
# threshold that admitted it.
#
# Both come from a small CSV the Bambu module writes while the SummarizedExperiment
# is still in memory. Reading the RDS back here would load every assay matrix to
# recover two columns.
tx_info$BambuTxClass <- NA_character_
tx_info$BambuNDR     <- NA_real_

if (!is.null(opt$tx_classes) && file.exists(opt$tx_classes)) {
    cat("Reading Bambu transcript classes...\n")
    bambu_meta <- read.csv(opt$tx_classes, stringsAsFactors = FALSE)

    idx <- match(tx_info$qry_id, as.character(bambu_meta$TXNAME))
    tx_info$BambuTxClass <- as.character(bambu_meta$txClassDescription)[idx]
    tx_info$BambuNDR     <- as.numeric(bambu_meta$NDR)[idx]

    cat(sprintf("  %d of %d novel transcripts carry a Bambu class\n",
                sum(!is.na(tx_info$BambuTxClass)), nrow(tx_info)))
} else {
    cat("No Bambu transcript classes supplied; BambuTxClass will be NA throughout.\n")
}

# Load counts
cat("Filtering transcripts with zero counts...\n")
tx_counts <- read_table(opt$tx_counts)
gene_counts <- read_table(opt$gene_counts)

# Remove novel transcripts that had 0 transcript counts in all samples
tx_info <- tx_info[tx_info$qry_id %in% tx_counts$TXNAME, ]

#' Human-readable label for a gffcompare class code. Kept in one place because
#' three outputs need the same mapping and copies of it drifted apart.
#'
#' The wording follows gffcompare's own definitions rather than paraphrasing them.
#' A shorter gloss reads better on a figure but invites the reader to reason from
#' the paraphrase instead of from what the tool actually asserts, and the two are
#' not always the same thing -- see the note on class code assignment in the
#' report's introduction.
CLASS_LABELS <- c(
    u = 'unknown or intergenic (u)',
    i = 'intronic (i)',
    x = 'antisense exonic overlap (x)',
    j = 'multi-exon SJ match (j)',
    k = 'contains reference (k)',
    o = 'sense exonic overlap (o)',
    y = 'contains reference within intron (y)',
    m = 'retained intron (m)',
    n = 'retained intron (n)'
)

#' Human-readable label for a class code, with i qualified by orientation.
#'
#' Containment within an intron covers two situations that warrant different
#' scrutiny. A model lying inside a same-strand intron cannot be told apart from a
#' fragment of that gene's unspliced precursor without independent evidence -- the
#' lncDACH1 problem -- while an antisense one is not explainable that way at all.
#' The orientation is a pulposeq addition, not part of gffcompare's definition, so
#' it is appended in parentheses rather than replacing the wording. Where the
#' orientation is unknown the label stays unqualified rather than guessing.
classify_class_code <- function(codes, same_strand = NULL) {
    codes <- as.character(codes)
    out   <- unname(CLASS_LABELS[codes])

    if (!is.null(same_strand)) {
        is_i      <- !is.na(codes) & codes == 'i' & !is.na(same_strand)
        out[is_i] <- ifelse(same_strand[is_i], 'sense intronic (i)',
                                               'antisense intronic (i)')
    }
    out
}

tx_info$classification <- classify_class_code(tx_info$class_code,
                                              tx_info$same_strand_as_compared)

# Save the metadata of all novel transcripts. Written after the classification so
# the full table carries it as well, including for the class codes that never
# become candidates.
cat("Saving novel transcripts metadata...\n")
write.csv(tx_info, file="novel_transcripts_metadata.csv", row.names = FALSE)

#' Build the attribute set written into the novel GTFs.
novel_attrs <- function(meta) {
    gene_biotype <- meta$compared_gene_biotype
    gene_biotype[is.na(gene_biotype)] <- "novel"

    list(
        transcript_status           = rep("novel", nrow(meta)),
        # pulposeq_ prefix: these are this pipeline's calls, not the annotation's.
        # A novel model sitting in an annotated gene would otherwise carry the
        # reference's gene_type and a pipeline gene_biotype on the same row with
        # nothing to say which is which.
        pulposeq_transcript_biotype = as.character(meta$pulposeq_transcript_biotype),
        pulposeq_gene_biotype       = gene_biotype,
        class_code                  = as.character(meta$class_code),
        classification              = as.character(meta$classification),
        BambuTxClass                = as.character(meta$BambuTxClass),
        BambuNDR                    = as.character(meta$BambuNDR),
        compared_gene_id            = meta$compared_gene_id,
        compared_gene_name          = as.character(meta$compared_gene_name),
        compared_gene_biotype       = as.character(meta$compared_gene_biotype),
        compared_transcript_id      = as.character(meta$compared_transcript_id),
        compared_transcript_name    = as.character(meta$compared_transcript_name),
        compared_transcript_biotype = as.character(meta$compared_transcript_biotype),
        gene_name                   = as.character(meta$qry_gene_name),
        qry_gene_biotype            = as.character(meta$qry_gene_biotype)
    )
}
#' Subset the Bambu GTF to a set of transcripts, attach the novel attributes and
#' write it out. Returns the exon records, which the exon-length summaries reuse.
write_novel_gtf <- function(meta, path) {
    ids   <- meta$qry_id
    tx    <- subset(gtf, type == "transcript" & transcript_id %in% ids)
    exons <- subset(gtf, type == "exon" & transcript_id %in% ids)

    out <- c(tx, exons)
    out <- out[order(out$transcript_id, out$type == "transcript", decreasing = TRUE)]
    out <- annotate_gtf(out, ids, novel_attrs(meta))
    export(out, path)

    invisible(exons)
}

# ---------------------------------------------------------------------------
# Routing novel models into biotypes
# ---------------------------------------------------------------------------
#
# Two groups, split by whether the model shares splice structure with the
# reference it was matched against.
#
# u, i and x share none of it. u has no reference at all, i lies wholly inside an
# intron, and x is exonic overlap on the opposite strand. For these the coding
# prediction stands on its own, so a non-coding call becomes novel_lncRNA whatever
# the reference gene happens to be: a sense-intronic or antisense lncRNA sitting
# inside a protein-coding locus is still an lncRNA, and both Ensembl and GENCODE
# annotate them as such.
#
# j, k, o, y, m and n all do share structure with a reference transcript, and
# there a non-coding call is much weaker evidence. An unproductive isoform of a
# protein-coding gene -- a retained intron, an NMD target, a truncated model --
# predicts non-coding too, and calling one a novel lncRNA asserts a new non-coding
# gene at a locus that already has a coding one. So lncRNA is claimed only where
# the reference gene is itself an lncRNA. Everything else in this group is filed by
# the prediction alone.
#
# novel_non_coding is deliberately not novel_lncRNA. It says the model has no
# coding potential without asserting what it is, which is as much as the evidence
# supports for a non-coding isoform of a coding gene.
#
# Note = and c never appear at all: only novel transcripts reach gffcompare, so a
# model matching or contained by a reference was already resolved as annotated.
INDEPENDENT_CODES <- c('u', 'i', 'x')
RELATED_CODES     <- c('j', 'k', 'o', 'y', 'm', 'n')
CANDIDATE_CODES   <- c(INDEPENDENT_CODES, RELATED_CODES)

# Reference gene biotypes that count as lncRNA when deciding whether a
# structurally-related model can be called one. Ensembl and GENCODE have
# consolidated on "lncRNA"; an annotation predating that merge would need its own
# names adding here.
LNCRNA_GENE_BIOTYPES <- c(
    "lncRNA"
)

is_lnc_ref  <- !is.na(tx_info$compared_gene_biotype) &
                   tx_info$compared_gene_biotype %in% LNCRNA_GENE_BIOTYPES
is_coding   <- !is.na(tx_info$prediction) & tx_info$prediction == 'coding'
is_noncod   <- !is.na(tx_info$prediction) & tx_info$prediction == 'non-coding'
independent <- tx_info$class_code %in% INDEPENDENT_CODES
related     <- tx_info$class_code %in% RELATED_CODES

biotype <- rep(NA_character_, nrow(tx_info))
biotype[(independent | related) & is_coding] <- "novel_protein_coding"
biotype[independent & is_noncod]             <- "novel_lncRNA"
biotype[related & is_noncod & is_lnc_ref]    <- "novel_lncRNA"
biotype[related & is_noncod & !is_lnc_ref]   <- "novel_non_coding"
tx_info$pulposeq_transcript_biotype <- biotype

# The 500 nt floor is the consensus lower bound for a long non-coding RNA. Applied
# to every candidate rather than only the lncRNA branch, so the three categories
# stay comparable to one another.
eligible <- tx_info[tx_info$class_code %in% CANDIDATE_CODES &
                        tx_info$len >= 500 &
                        !is.na(tx_info$pulposeq_transcript_biotype), ]

routed <- table(eligible$pulposeq_transcript_biotype)
cat(sprintf("Routed %d candidates: %s\n", nrow(eligible),
            paste(sprintf("%s=%d", names(routed), as.integer(routed)),
                  collapse = ", ")))

# Select novel lncRNA candidates
cat("Processing lncRNA candidates...\n")
new_lncRNAs <- eligible[eligible$pulposeq_transcript_biotype == "novel_lncRNA", ]

write.csv(new_lncRNAs, file="novel_lncRNAs_metadata.csv", row.names = FALSE)
new_lncRNAs_exons_gtf <- write_novel_gtf(new_lncRNAs, "novel_lncRNAs.gtf")

# Select novel protein-coding candidates
cat("Processing protein-coding candidates...\n")
new_mRNAs <- eligible[eligible$pulposeq_transcript_biotype == "novel_protein_coding", ]

write.csv(new_mRNAs, file="novel_protein-coding_metadata.csv", row.names = FALSE)
new_mRNAs_exons_gtf <- write_novel_gtf(new_mRNAs, "novel_protein-coding.gtf")

# Select the non-coding models that are not claimed as lncRNA
cat("Processing non-coding candidates...\n")
new_ncRNAs <- eligible[eligible$pulposeq_transcript_biotype == "novel_non_coding", ]

write.csv(new_ncRNAs, file="novel_non_coding_metadata.csv", row.names = FALSE)
write_novel_gtf(new_ncRNAs, "novel_non_coding.gtf")

# Save combined metadata. All three categories: the novel_non_coding models must
# not drop out of the validated GTF or the counts derived from it.
new_mRNA_lncRNA <- rbind(new_mRNAs, new_lncRNAs, new_ncRNAs)
write.csv(new_mRNA_lncRNA, "novel_transcripts_validated_metadata.csv", row.names = FALSE)

# Update novel transcripts GTF with the new classifications
write_novel_gtf(new_mRNA_lncRNA, "novel_transcripts_validated.gtf")

# Get exon lengths
cat("Calculating exon lengths...\n")
new_lncRNA_exon_len <- data.frame(
  transcript_id = new_lncRNAs_exons_gtf$transcript_id,
  exon_number = new_lncRNAs_exons_gtf$exon_number,
  width = width(new_lncRNAs_exons_gtf)
) 

write.csv(new_lncRNA_exon_len, "novel_lncRNA_exon_lengths.csv", row.names = FALSE)

new_mRNA_exon_len <- data.frame(
  transcript_id = new_mRNAs_exons_gtf$transcript_id,
  exon_number = new_mRNAs_exons_gtf$exon_number,
  width = width(new_mRNAs_exons_gtf)
) 

write.csv(new_mRNA_exon_len, "novel_protein-coding_exon_lengths.csv", row.names = FALSE)

# Save counts for novel mRNAs and lncRNAs
tx_ids <- new_mRNA_lncRNA$qry_id
gn_ids <- new_mRNA_lncRNA$qry_gene_id

novel_tx_counts <- subset(tx_counts, TXNAME %in% tx_ids)
write.csv(novel_tx_counts, "bambu_novel_tx_counts.csv", row.names = FALSE)

novel_gn_counts <- subset(gene_counts, GENEID %in% gn_ids)
write.csv(novel_gn_counts, "bambu_novel_gene_counts.csv", row.names = FALSE)

cat("Analysis completed successfully!\n")
