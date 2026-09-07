#!/usr/bin/env Rscript

# Genomic context figures, drawn with plotgardener. Two families of PNG -- one per
# selected gene, one per sense-intronic candidate -- each stacking:
#
#   title       the subject and the region actually plotted, chromosome included
#   chromosome  a bar marking where on the chromosome this window sits
#   structure   the gene's (or host's) annotated exon footprint, collapsed to one row
#   isoforms    every transcript in the window -- detected, novel, and annotated but
#               not detected -- coloured by biotype and labelled with its status
#   coverage    one track per sample, on a y range shared across samples
#   axis        genomic coordinates
#   legend      the biotypes this panel contains
#
# Two scope decisions, both deliberate.
#
# Only known genes carrying at least one novel isoform are drawn. A wholly novel locus
# has no annotated model to place it against, and that comparison is most of what makes
# the view worth looking at. No biotype filter applies: a novel lncRNA inside a known
# lncRNA gene competes for a panel on the same terms as one inside a protein-coding
# gene.
#
# Both tracks are drawn against the REFERENCE annotation, not only the validated set,
# and that is what makes a panel interpretable.
#
# The structure row is the gene's full annotated exon footprint, taken from
# --annotation. It has to be, because gffcompare assigns class codes against the whole
# annotation. Built from the validated set instead, the row omits the exons of
# undetected isoforms and a panel can contradict its own caption.
#
# The isoform rows carry three kinds of model: detected known transcripts, novel
# models, and annotated transcripts with no support in these samples.
# Undetected models keep their real biotype colour and carry "not detected" in
# the label, so they are never read as observed -- the coverage tracks below say what
# was actually seen.

suppressPackageStartupMessages({
  library(optparse)
  library(rtracklayer)
  library(plotgardener)
})

option_list <- list(
  make_option(c("-t", "--gtf"), type = "character", default = NULL,
              help = "Validated, enriched transcriptome GTF", metavar = "character"),
  make_option(c("-w", "--bigwigs"), type = "character", default = NULL,
              help = "Comma-separated bigWig files, one per sample", metavar = "character"),
  make_option(c("-s", "--names"), type = "character", default = NULL,
              help = "Comma-separated sample labels, same order as --bigwigs", metavar = "character"),
  make_option(c("-o", "--outdir"), type = "character", default = ".",
              help = "Output directory [default %default]", metavar = "character"),
  make_option(c("-c", "--counts"), type = "character", default = NULL,
              help = "Validated transcript counts matrix, one column per sample", metavar = "character"),
  make_option("--fl_counts", type = "character", default = NULL,
              help = "Validated full-length transcript counts matrix", metavar = "character"),
  make_option(c("-a", "--annotation"), type = "character", default = NULL,
              help = "Reference annotation GTF, for annotated gene and intron structure", metavar = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$gtf) || is.null(opt$bigwigs)) {
  stop("--gtf and --bigwigs are both required.", call. = FALSE)
}
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# The cap on how many panels the stratified selection draws, and the isoform count
# above which a locus stops being legible in one panel. Fixed rather than exposed:
# these are properties of what fits on a page, not choices a user needs to make.
#
# The cap rarely bites. Most strata are empty on any given run -- not every class
# code turns up among all three novel biotypes -- so the usual output is well under
# it.
MAX_PANELS <- 24L
MAX_TX     <- 12L

# How many sense-intronic candidates to draw. A model inside a host intron on the
# host's own strand is sequence-identical to part of that host's unspliced precursor,
# and no property of the assembly separates the two. The panel shows what each one
# actually looks like -- coverage flat across the whole intron and continuous with the
# flanking exons, or a discrete block with quiet intron either side -- and leaves the
# reading to the eye.
N_INTRONIC <- 5L

# A floor on total count, so a candidate with three reads does not take a panel from
# one with two thousand. Skipped if it would empty the set, as it would on a small run.
MIN_COUNTS <- 20L

# The drawing window for a candidate is the host intron it sits in, padded so
# the flanking exons are on screen. Without them the panel cannot show whether
# coverage stops at the intron boundary or runs straight through it, which is the only
# question it exists to answer.
INTRON_PAD_FRAC <- 0.15
INTRON_PAD_MIN  <- 1000L

# Isoforms are coloured by transcript biotype, and the legend under each panel lists
# the biotypes that panel actually contains. Novel biotypes take the saturated end of
# the palette so a novel model still reads as novel at a glance -- which is what the
# old novel-only highlight did -- while known biotypes stay muted. Anything unlisted
# falls back to grey rather than erroring: an annotation can always carry a biotype
# this list has not seen.
BIOTYPE_COLORS <- c(
  novel_lncRNA                   = "#D95F02",
  novel_protein_coding           = "#7570B3",
  novel_non_coding               = "#1B9E77",
  protein_coding                 = "#7FBC41",
  lncRNA                         = "#4393C3",
  retained_intron                = "#8DA0CB",
  nonsense_mediated_decay        = "#E7969C",
  protein_coding_CDS_not_defined = "#BDBDBD",
  processed_transcript           = "#D9D9D9",
  processed_pseudogene           = "#C7C7C7",
  unprocessed_pseudogene         = "#C7C7C7",
  "biotype not set"              = "#969696"
)
BIOTYPE_FALLBACK <- "#969696"

# Isoform track geometry, passed straight to plotTranscripts.
#
# boxHeight is raised from plotgardener's 2 mm default for contrast: a taller exon box
# against the same thin intron line is the thick-against-thin reading we want, and
# plotTranscripts draws UTR thinner than CDS by itself, so that comes free.
#
# plotTranscripts centres a label on the transcript's true midpoint and does not draw it
# when that falls outside the window. Gene panels never hit this,
# because their window spans min(start) to max(end) over the gene's own transcripts, so
# every midpoint is inside it by construction. Intronic panels are windowed on one host
# intron, so host isoforms running the length of the gene lose their names -- see the
# note on the isoform track in draw_panel().
TX_LABEL_SIZE <- 6
# A 6 pt label is about 0.083 in tall, so a 0.09 in band left it touching the model
# above it. The label band clears the text, and the row gap is trailing space so one
# row's model cannot run into the next row's label.
TX_LABEL_H    <- 0.13    # band above each model for its label
TX_ROW_GAP    <- 0.08    # clear space below a model, before the next row's label
TX_BOX_MM     <- 4       # exon box height, against plotgardener's 2 mm default
TX_SPACE_H    <- 0.2     # padding inside a model's own band, as a fraction of the box
TX_STROKE     <- 0.1     # transcript body outline, plotgardener's default

TX_MODEL_H <- (TX_BOX_MM / 25.4) * (1 + TX_SPACE_H)
TX_ROW_H   <- TX_LABEL_H + TX_MODEL_H + TX_ROW_GAP

# The isoform track calls plotTranscripts once per transcript rather than once per
# panel, and this is deliberate -- do not collapse it back into a single call.
#
# One call per panel is the obvious way to write it and was what ran first. plotgardener
# then owns the row packing and the labels, and it places a label at the transcript's
# midpoint and omits it when that midpoint falls outside the window. In a panel windowed
# on one host intron, every host isoform running the length of the gene has its midpoint
# off-window, so it goes unnamed: the ACOT7 panel drew ten models and named two.
#
# Filtering each call to a single transcript with transcriptFilter, at a row chosen
# here, moves both decisions to us. The label is then drawn separately and clamped
# inside the track, so an off-window midpoint still gets a name. plotgardener still
# renders the model, so UTR-vs-CDS thickness and exon structure are unchanged.
#
# The cost is one call per model instead of one per panel -- roughly a hundred across a
# typical run rather than ten. Both were rendered side by side on the same data before
# choosing this one.

biotype_color <- function(bt) {
  out <- unname(BIOTYPE_COLORS[bt])
  out[is.na(out)] <- BIOTYPE_FALLBACK
  out
}

# No OrgDb is needed, and that is a property of the assembly below rather than an
# oversight. plotgardener consults the OrgDb only to translate one identifier type
# into another: getExons() skips the lookup entirely when gene.id.column equals
# display.column, and only calls AnnotationDbi::select() when they differ. The
# defaults differ -- ENTREZID to SYMBOL -- which is why the documented examples all
# pass org.Hs.eg.db. Setting both to GENEID makes the lookup a no-op, and keeping it
# NULL is what lets this run for any organism. Bambu's novel genes have no entry in
# any OrgDb whatever the species.
ORGDB <- NULL

bw_files <- trimws(strsplit(opt$bigwigs, ",")[[1]])
bw_files <- bw_files[nzchar(bw_files)]

# Chromosome lengths for the ideogram, taken from the bigWig header rather than from
# an assembly package: no extra input, no network, and it is guaranteed to agree with
# the coverage being drawn beneath it. If it cannot be read the ideogram is skipped
# and the rest of the panel is unaffected.
chrom_lens <- tryCatch(
  GenomeInfoDb::seqlengths(rtracklayer::BigWigFile(bw_files[1])),
  error = function(e) NULL)
sample_names <- if (!is.null(opt$names)) {
  trimws(strsplit(opt$names, ",")[[1]])
} else {
  sub("\\.bw$", "", basename(bw_files))
}
if (length(sample_names) != length(bw_files)) {
  stop("--names must have the same number of entries as --bigwigs.", call. = FALSE)
}

# --- Read the annotation once -------------------------------------------------

cat("Reading validated annotation:", opt$gtf, "\n")
# CDS is read as well as exon, and that is what makes the isoform track legible.
#
# plotTranscripts renders a coding region thicker than its untranslated ends, and it
# takes that from the TxDb's CDS. Given exon rows alone it has nothing to work from
# and every model -- protein-coding or not -- draws as one uniform box. The known
# transcripts in --gtf carry the reference's CDS for exactly this reason.
#
# Novel models have no CDS and cannot: Bambu does not call ORFs. They draw as
# uniform boxes, which is the honest depiction, since the coding-potential
# prediction is a statement about the sequence and not a claim about where a start
# codon sits.
gtf <- rtracklayer::import(opt$gtf, feature.type = c("transcript", "exon", "CDS"))
md  <- S4Vectors::mcols(gtf)

col_or_na <- function(name) {
  if (name %in% colnames(md)) as.character(md[[name]]) else NA_character_
}

is_tx  <- as.character(md$type) == "transcript"
tx_gr  <- gtf[is_tx]
tx_md  <- S4Vectors::mcols(tx_gr)

tx_col <- function(name) {
  if (name %in% colnames(tx_md)) as.character(tx_md[[name]]) else NA_character_
}

tx <- data.frame(
  gene_id            = tx_col("gene_id"),
  gene_name          = tx_col("gene_name"),
  transcript_id      = tx_col("transcript_id"),
  transcript_name    = tx_col("transcript_name"),
  transcript_biotype = tx_col("transcript_biotype"),
  transcript_status  = tx_col("transcript_status"),
  class_code         = tx_col("class_code"),
  classification     = tx_col("classification"),
  # What gffcompare matched the model against. Needed because the novel model and
  # the gene it relates to are often not the same gene_id: Bambu gives an intronic
  # or antisense model its own gene, so a panel built from the model's own gene
  # alone would show it against nothing.
  ref_gene_id        = tx_col("ref_gene_id"),
  ref_gene_name      = tx_col("ref_gene_name"),
  ref_gene_biotype   = tx_col("ref_gene_biotype"),
  # The reference TRANSCRIPT a novel model was classified against, which decides
  # which undetected isoforms are worth drawing, and names the match on the label.
  ref_transcript_id   = tx_col("ref_transcript_id"),
  ref_transcript_name = tx_col("ref_transcript_name"),
  # Bambu's own account of which part of the model is new.
  BambuTxClass        = tx_col("BambuTxClass"),
  chrom              = as.character(GenomeInfoDb::seqnames(tx_gr)),
  strand             = as.character(BiocGenerics::strand(tx_gr)),
  start              = BiocGenerics::start(tx_gr),
  end                = BiocGenerics::end(tx_gr),
  stringsAsFactors   = FALSE
)
tx <- tx[!is.na(tx$gene_id) & nzchar(tx$gene_id), , drop = FALSE]

# Label each model by transcript_name where the annotation supplies one, falling
# back to the identifier. Bambu writes no transcript_name at all, so novel isoforms
# always fall back -- which is the desired result, since BambuTx is what you would
# search the GTF for. The biotype rides along because it is the thing a reader most
# wants to know about an isoform they have never seen before.
display_name <- ifelse(!is.na(tx$transcript_name) & nzchar(tx$transcript_name),
                       tx$transcript_name, tx$transcript_id)
biotype <- ifelse(is.na(tx$transcript_biotype) | !nzchar(tx$transcript_biotype),
                  "biotype not set", tx$transcript_biotype)
# plotTranscripts labels each model with its TXNAME, so the label has to BE the
# transcript_id as far as the TxDb is concerned -- one string, used on the figure and
# in the rewritten sub-GTF alike.
#
# Biotype is left out of it: every model is coloured by biotype and the legend names
# it, so repeating it on each label spends a third of the track width restating the
# palette. The GTF keeps transcript_biotype as its own attribute, so nothing is lost.
# The bare class code goes too -- the mapping to its wording is one-to-one, "intronic"
# says everything "i" does, and the letter is still on every row of the metadata CSV.
cls <- ifelse(is.na(tx$classification) | !nzchar(tx$classification), "",
              paste0(" | ", tx$classification))

# What the model was matched against, named on the label itself. A class code says
# a novel model shares a junction with "a reference transcript"; without naming
# which, the panel leaves the reader to work it out from the tracks. The reference
# transcript's name where the annotation gives one, its identifier otherwise.
ref_tx <- ifelse(!is.na(tx$ref_transcript_name) & nzchar(tx$ref_transcript_name),
                 tx$ref_transcript_name, tx$ref_transcript_id)
ref_mark <- ifelse(is.na(ref_tx) | !nzchar(ref_tx), "", paste0(" | ", ref_tx))

str_mark <- ifelse(tx$strand %in% c("+", "-"), paste0(" | ", tx$strand), "")

tx$label <- paste0(display_name, cls, ref_mark, str_mark)
# Kept as its own column because the figures colour by it and the legend lists it.
tx$biotype_key <- biotype

# --- Candidate selection ------------------------------------------------------

# "Known gene" is decided by whether any of its transcripts came from the reference,
# not by an identifier prefix. Prefixes differ between Ensembl, GENCODE and every
# non-model organism, and Bambu assigns novel isoforms at a known locus the reference
# gene_id anyway, so the status column is the only reliable signal.
by_gene <- split(seq_len(nrow(tx)), tx$gene_id)

genes <- do.call(rbind, lapply(names(by_gene), function(gid) {
  idx <- by_gene[[gid]]
  nm  <- tx$gene_name[idx]
  nm  <- nm[!is.na(nm) & nzchar(nm)]
  data.frame(
    gene_id     = gid,
    gene_name   = if (length(nm)) nm[1] else NA_character_,
    chrom       = tx$chrom[idx][1],
    start       = min(tx$start[idx]),
    end         = max(tx$end[idx]),
    n_tx        = length(idx),
    n_known     = sum(tx$transcript_status[idx] == "known", na.rm = TRUE),
    n_novel     = sum(tx$transcript_status[idx] == "novel", na.rm = TRUE),
    n_novel_lnc = sum(tx$transcript_biotype[idx] == "novel_lncRNA", na.rm = TRUE),
    biotypes    = paste(sort(unique(tx$transcript_biotype[idx])), collapse = ";"),
    stringsAsFactors = FALSE
  )
}))

# Coarsen the reference gene's biotype to the distinction that decides how a novel
# call should be read. The exact biotype still reaches the panel subtitle; this is
# only the stratifying key, and leaving it fine-grained would scatter single
# transcripts across a dozen near-empty strata.
LNCRNA_REF_BIOTYPES <- "lncRNA"

ref_class_of <- function(bt) {
  ifelse(is.na(bt) | !nzchar(bt),     "no reference",
  ifelse(bt %in% LNCRNA_REF_BIOTYPES, "lncRNA gene",
  ifelse(bt == "protein_coding",      "protein-coding gene", "other gene")))
}

slug <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "-", x)
  substr(gsub("(^-|-$)", "", x), 1, 60)
}

# One panel per stratum, rather than the five richest loci overall.
#
# The flat top-N ranked on novel lncRNA count, which is a popularity contest: it
# returned five variations on the same finding -- the most isoform-rich loci -- and
# never showed what an antisense lncRNA, or a retained intron inside a coding gene,
# actually looks like, because those loci are smaller. A reader cannot judge a class
# of call they have never been shown an example of.
#
# A stratum is (novel biotype, classification, reference gene class), so the panel
# set spans the ways a novel model can relate to the annotation. A sense-intronic
# lncRNA inside a protein-coding gene is a different claim from the same model
# inside an lncRNA gene, and that pair is the lncDACH1 argument in one figure each.
#
# Within a stratum the richest locus wins, so each panel is the most informative
# example of its kind.
#
# u is excluded: it is intergenic, so there is no reference gene to draw it against
# and the panel would be the novel model alone on an empty window.
novel_tx <- tx[tx$transcript_status %in% "novel" &
                 !is.na(tx$class_code) & !tx$class_code %in% "u", , drop = FALSE]

n_tx_of    <- setNames(genes$n_tx,    genes$gene_id)
n_known_of <- setNames(genes$n_known, genes$gene_id)
n_novel_of <- setNames(genes$n_novel, genes$gene_id)
lookup0 <- function(map, key) { v <- unname(map[key]); v[is.na(v)] <- 0L; as.integer(v) }

if (nrow(novel_tx)) {
  # A panel is the novel model's own gene plus the reference gene it was matched
  # against, which are frequently different and occasionally the same.
  ref_same <- !is.na(novel_tx$ref_gene_id) & novel_tx$ref_gene_id == novel_tx$gene_id
  ref_use  <- !is.na(novel_tx$ref_gene_id) & !ref_same

  novel_tx$panel_n_tx    <- lookup0(n_tx_of,    novel_tx$gene_id) +
    ifelse(ref_use, lookup0(n_tx_of,    novel_tx$ref_gene_id), 0L)
  novel_tx$panel_n_known <- lookup0(n_known_of, novel_tx$gene_id) +
    ifelse(ref_use, lookup0(n_known_of, novel_tx$ref_gene_id), 0L)
  novel_tx$panel_n_novel <- lookup0(n_novel_of, novel_tx$gene_id) +
    ifelse(ref_use, lookup0(n_novel_of, novel_tx$ref_gene_id), 0L)

  novel_tx$ref_class <- ref_class_of(novel_tx$ref_gene_biotype)
  novel_tx$stratum   <- paste(novel_tx$transcript_biotype, novel_tx$classification,
                              novel_tx$ref_class, sep = " | ")

  # A panel needs a known isoform to supply the context, and has to stay legible.
  novel_tx <- novel_tx[novel_tx$panel_n_known > 0 & novel_tx$panel_n_tx <= MAX_TX, ,
                       drop = FALSE]
}

cat(sprintf("%d novel models across %d strata are drawable within the %d-isoform limit\n",
            nrow(novel_tx), length(unique(novel_tx$stratum)), MAX_TX))

cand <- genes[0, , drop = FALSE]

if (nrow(novel_tx)) {
  novel_tx <- novel_tx[order(-novel_tx$panel_n_tx, -novel_tx$panel_n_novel,
                             novel_tx$transcript_id), , drop = FALSE]

  # Strata taken in order of how many models they hold, so if the cap does bite it
  # drops the thinnest evidence rather than whatever happened to sort last.
  strat_order <- names(sort(table(novel_tx$stratum), decreasing = TRUE))

  locus_of <- paste(novel_tx$gene_id, novel_tx$ref_gene_id, sep = "~")
  picked   <- integer(0)
  used     <- character(0)

  for (s in strat_order) {
    rows <- which(novel_tx$stratum == s)
    free <- rows[!locus_of[rows] %in% used]
    # Fall back to the best example even when its locus is already drawn: the strata
    # are what the figure set is meant to span, and a repeated locus under a
    # different heading still answers a different question about it.
    take   <- if (length(free)) free[1] else rows[1]
    used   <- c(used, locus_of[take])
    picked <- c(picked, take)
    if (length(picked) >= MAX_PANELS) break
  }

  sel <- novel_tx[picked, , drop = FALSE]

  cand <- do.call(rbind, lapply(seq_len(nrow(sel)), function(i) {
    r  <- sel[i, ]
    pg <- unique(c(r$gene_id, r$ref_gene_id))
    pg <- pg[!is.na(pg) & pg %in% genes$gene_id]
    g  <- genes[match(pg, genes$gene_id), , drop = FALSE]
    # A reference on another sequence is a broken match, not a panel. The model's own
    # gene is always on its own chromosome, so this never empties g.
    g  <- g[g$chrom == r$chrom, , drop = FALSE]

    data.frame(
      gene_id     = if (!is.na(r$ref_gene_id) && r$ref_gene_id %in% g$gene_id)
                      r$ref_gene_id else r$gene_id,
      gene_name   = if (!is.na(r$ref_gene_name) && nzchar(r$ref_gene_name))
                      r$ref_gene_name
                    else if (!is.na(r$ref_gene_id)) r$ref_gene_id else r$gene_id,
      chrom       = r$chrom,
      start       = min(g$start),
      end         = max(g$end),
      n_tx        = sum(g$n_tx),
      n_known     = sum(g$n_known),
      n_novel     = sum(g$n_novel),
      n_novel_lnc = sum(g$n_novel_lnc),
      biotypes    = paste(sort(unique(unlist(strsplit(g$biotypes, ";")))), collapse = ";"),
      panel_genes = paste(g$gene_id, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }))

  cand$stratum        <- sel$stratum
  cand$novel_biotype  <- sel$transcript_biotype
  cand$classification <- sel$classification
  cand$class_code     <- sel$class_code
  cand$ref_class      <- sel$ref_class
  cand$example_tx     <- sel$transcript_id
  cand$ref_gene_biotype <- sel$ref_gene_biotype
  cand$BambuTxClass   <- sel$BambuTxClass

  cand$label <- ifelse(is.na(cand$gene_name) | !nzchar(cand$gene_name),
                       cand$gene_id, cand$gene_name)
  cand$pad   <- ceiling((cand$end - cand$start) * 0.05)
  cand$win_s <- pmax(1, cand$start - cand$pad)
  cand$win_e <- cand$end + cand$pad

  cat(sprintf("Drawing %d stratified panels:\n%s\n", nrow(cand),
              paste0("  ", cand$stratum, "  -> ", cand$label, collapse = "\n")))
} else {
  message("No novel model outside class code u has a drawable reference context.")
}

# --- Sense-intronic candidates ------------------------------------------------

# The second selection, and a different unit: one transcript rather than one gene.
# Sense intronic is the one relationship position cannot resolve -- a model inside a
# host intron, on the host's own strand, is sequence-identical to part of that host's
# unspliced precursor. Nothing here scores them; they are drawn so the report can show
# one beside a clean call and leave the reading to the eye.
#
# `classification` is written by novel_transcripts.R from the class code and the strand
# comparison, and rides on the enriched GTF this script already reads, so the selection
# needs nothing from the alignments.
flagged <- data.frame()

sense_i <- which(grepl("^sense intronic", tx$classification))

if (length(sense_i)) {
  flagged <- data.frame(
    qry_id           = tx$transcript_id[sense_i],
    ref_gene_id      = tx$ref_gene_id[sense_i],
    ref_gene_biotype = tx$ref_gene_biotype[sense_i],
    class_code       = tx$class_code[sense_i],
    strand           = tx$strand[sense_i],
    BambuTxClass     = tx$BambuTxClass[sense_i],
    stringsAsFactors = FALSE
  )

  # Exon count is not carried on `tx`, and it is the one number that says whether a
  # candidate is fully unspliced.
  all_exons <- gtf[as.character(S4Vectors::mcols(gtf)$type) == "exon"]
  n_exons   <- table(as.character(S4Vectors::mcols(all_exons)$transcript_id))
  flagged$num_exons <- as.integer(n_exons[flagged$qry_id])

  # Per-sample counts. Not evidence about the candidate: this decides whether a panel
  # has anything to look at, and which coverage tracks the figure marks in bold.
  read_counts <- function(path) {
    if (is.null(path) || !file.exists(path)) return(NULL)
    df <- tryCatch(read.delim(path, header = TRUE, sep = "\t", check.names = FALSE),
                   error = function(e) NULL)
    if (is.null(df) || !nrow(df)) return(NULL)
    num <- df[, vapply(df, is.numeric, logical(1)), drop = FALSE]
    if (!ncol(num)) return(NULL)
    rownames(num) <- as.character(df[[1]])
    num
  }

  cnt <- read_counts(opt$counts)
  if (!is.null(cnt)) {
    m <- cnt[match(flagged$qry_id, rownames(cnt)), , drop = FALSE]
    flagged$counts_total       <- rowSums(m, na.rm = TRUE)
    flagged$samples_quantified <- rowSums(m > 0, na.rm = TRUE)
    flagged$samples_total      <- ncol(cnt)
    # Semicolon-separated, since a sample name may itself contain a comma.
    flagged$quantified_samples <- vapply(seq_len(nrow(m)), function(i) {
      hit <- colnames(m)[which(!is.na(m[i, ]) & m[i, ] > 0)]
      if (length(hit)) paste(hit, collapse = ";") else NA_character_
    }, character(1))
    flagged$counts_total[is.na(flagged$counts_total)]             <- 0
    flagged$samples_quantified[is.na(flagged$samples_quantified)] <- 0L
  } else {
    flagged$counts_total       <- NA_real_
    flagged$samples_quantified <- NA_integer_
    flagged$samples_total      <- NA_integer_
    flagged$quantified_samples <- NA_character_
  }

  # Membership of Bambu's full-length matrix: a zero/non-zero fact, not a threshold.
  fl <- read_counts(opt$fl_counts)
  flagged$full_length_support <- if (is.null(fl)) NA else flagged$qry_id %in% rownames(fl)

  if (nrow(flagged)) {
    enough <- !is.na(flagged$counts_total) & flagged$counts_total >= MIN_COUNTS
    if (any(enough)) flagged <- flagged[enough, , drop = FALSE]

    # Ranked by structure, not by a score. Mono-exonic first -- an unspliced model
    # inside an intron is the profile a precursor fragment has -- then models Bambu
    # never spanned end to end, then by count so a panel has coverage worth looking at.
    flagged <- flagged[order(-(flagged$num_exons %in% 1L),
                             -(flagged$full_length_support %in% FALSE),
                             -flagged$counts_total,
                             flagged$qry_id), , drop = FALSE]
    flagged <- head(flagged, N_INTRONIC)
  }
}

# --- Reference exon structure -------------------------------------------------

# The gene track is drawn from the REFERENCE annotation, not from the validated
# GTF, and that distinction decides whether a panel can be checked at all.
#
# gffcompare assigns every class code against the full annotation, including
# isoforms never detected in these samples. The validated GTF holds only what was
# detected. Building the gene track from the latter produces panels that
# contradict their own captions: a model coded `x` -- antisense exonic overlap --
# against an undetected isoform appears to sit in an empty intron, with no exon
# drawn anywhere near it, because the isoform carrying that exon was filtered out
# before the figure was made.
#
# Reduced per gene, so the track is the gene's full annotated exonic footprint and
# the class code refers to something visible.
bare_id <- function(x) sub("\\.[0-9]+$", "", as.character(x))

ref_exons_by_gene   <- list()
ref_introns_by_gene <- list()

wanted_genes <- unique(c(
  if (nrow(cand))    unlist(strsplit(cand$panel_genes, ";")) else character(0),
  if (nrow(flagged)) flagged$ref_gene_id else character(0)))
wanted_genes <- wanted_genes[!is.na(wanted_genes) & nzchar(wanted_genes)]

if (length(wanted_genes) && !is.null(opt$annotation) && file.exists(opt$annotation)) {
  ref_ex  <- rtracklayer::import(opt$annotation, feature.type = "exon")
  ref_gid <- as.character(S4Vectors::mcols(ref_ex)$gene_id)

  # Ensembl carries the gene version in a separate attribute and GENCODE inline, so
  # a novel record can name either form. Both sides are matched on the bare id.
  ref_ex <- ref_ex[bare_id(ref_gid) %in% bare_id(wanted_genes)]

  if (length(ref_ex)) {
    by_g <- split(ref_ex, bare_id(S4Vectors::mcols(ref_ex)$gene_id))
    ex_r <- GenomicRanges::reduce(by_g)
    ref_exons_by_gene <- as.list(ex_r)

    g_r <- unlist(range(by_g), use.names = TRUE)
    g_r <- g_r[!duplicated(names(g_r))]
    ref_introns_by_gene <- as.list(GenomicRanges::psetdiff(g_r, ex_r[names(g_r)]))
  }
}

if (nrow(flagged)) {
  # Window on the host intron the candidate sits in rather than the candidate
  # itself. A window drawn tight around a truncated fragment looks like a discrete
  # transcript no matter what it is; the flanking intron is where the difference
  # shows.
  tx_pos <- tx[match(flagged$qry_id, tx$transcript_id), c("chrom", "start", "end")]
  flagged$chrom <- tx_pos$chrom
  flagged$start <- tx_pos$start
  flagged$end   <- tx_pos$end

  # Drop anything the GTF cannot place.
  flagged <- flagged[!is.na(flagged$chrom), , drop = FALSE]
}

if (nrow(flagged)) {
  win <- t(vapply(seq_len(nrow(flagged)), function(i) {
    gid <- bare_id(flagged$ref_gene_id[i])
    ivs <- ref_introns_by_gene[[gid]]
    if (!is.null(ivs) && length(ivs)) {
      # The intron containing the candidate; if it straddles more than one, the
      # union of those it touches.
      hit <- which(BiocGenerics::start(ivs) <= flagged$end[i] &
                   BiocGenerics::end(ivs)   >= flagged$start[i])
      if (length(hit)) {
        i_s <- min(BiocGenerics::start(ivs)[hit])
        i_e <- max(BiocGenerics::end(ivs)[hit])
        # Padded past the intron edges: a window clipped exactly to the intron leaves
        # the flanking exons off screen, so read-through at the boundary is invisible.
        pad <- max(INTRON_PAD_MIN, ceiling((i_e - i_s) * INTRON_PAD_FRAC))
        return(c(i_s - pad, i_e + pad))
      }
    }
    # No usable intron: fall back to the candidate plus half its length either side,
    # which still shows whether coverage stops at the boundaries.
    pad <- ceiling((flagged$end[i] - flagged$start[i]) * 0.5)
    c(flagged$start[i] - pad, flagged$end[i] + pad)
  }, numeric(2)))

  flagged$win_s <- pmax(1, win[, 1])
  flagged$win_e <- win[, 2]
  flagged$ref_gene_name <- genes$gene_name[match(flagged$ref_gene_id, genes$gene_id)]
}

cat(sprintf("%d sense-intronic candidates selected for drawing\n", nrow(flagged)))

if (nrow(cand) == 0 && nrow(flagged) == 0) {
  message("Nothing to draw; writing empty candidate tables.")
  write.csv(genes[0, ], file.path(opt$outdir, "genomic_context_candidates.csv"),
            row.names = FALSE)
  write.csv(data.frame(), file.path(opt$outdir, "intronic_context_candidates.csv"),
            row.names = FALSE)
  quit(save = "no", status = 0)
}

# --- Annotation for plotgardener ----------------------------------------------

# plotTranscripts draws from a TxDb and labels each model with its TXNAME, so the label
# has to BE the transcript id as far as the TxDb is concerned. The GTF is rewritten
# accordingly -- but only over the handful of windows being drawn, which also keeps
# makeTxDbFromGFF off the whole-genome annotation. On this data that is a few thousand
# records instead of several million.
#
# Note plotTranscripts also has a transcriptFilter argument that takes transcript names
# to display. Restricting the GTF here rather than filtering per call is still the right
# way round: it is what keeps the TxDb build small, and one build is shared by every
# panel.
# Both selections share one TxDb: makeTxDbFromGFF over the flagged windows as well
# is a few hundred extra records, where building it twice would repeat the whole
# parse.
windows <- c(
  if (nrow(cand)) GenomicRanges::GRanges(
    cand$chrom, IRanges::IRanges(cand$win_s, cand$win_e)) else GenomicRanges::GRanges(),
  if (nrow(flagged)) GenomicRanges::GRanges(
    flagged$chrom, IRanges::IRanges(flagged$win_s, flagged$win_e)) else GenomicRanges::GRanges()
)
# --- Annotated transcripts at the drawn loci ----------------------------------
#
# The isoform track shows what was detected, and on its own that can leave a panel
# impossible to read. AK5 has an annotated antisense gene, AK5-AS1, sitting inside
# one of its introns. With none of its transcripts detected it never appeared, so a
# novel antisense model in that intron looked like it had arrived from nowhere --
# when the question worth asking is whether it is a new locus at all, or an
# undetected form of something already annotated. That question cannot be put
# without the annotation on the page.
#
# Annotated transcripts overlapping a drawn window are read in here, excluding any
# already in the validated set, which are drawn from their own record. Which of them
# a given panel actually draws is decided later, per panel, by whether a class code
# points at them. They keep their real biotype colour, so the palette still says
# what they are, and carry "not detected" in the label so they are never read as
# observed.
ref_draw_gr <- NULL
if (length(windows) && !is.null(opt$annotation) && file.exists(opt$annotation)) {
  # CDS here too, so an undetected annotated transcript draws with the same
  # thick-versus-thin structure as a detected one rather than as a flat bar.
  ref_all <- rtracklayer::import(opt$annotation,
                                 feature.type = c("transcript", "exon", "CDS"))
  ref_all <- IRanges::subsetByOverlaps(ref_all, windows)

  # Compared on the bare id: the validated GTF and the annotation can disagree on
  # whether the version suffix is carried inline, and a mismatch would draw the same
  # transcript twice -- once as detected and once as not.
  ref_txid <- as.character(S4Vectors::mcols(ref_all)$transcript_id)
  ref_all  <- ref_all[!is.na(ref_txid) &
                        !(bare_id(ref_txid) %in% bare_id(tx$transcript_id))]

  if (length(ref_all)) {
    ref_draw_gr <- ref_all

    rt     <- ref_all[as.character(S4Vectors::mcols(ref_all)$type) == "transcript"]
    rt_md  <- S4Vectors::mcols(rt)
    rcol   <- function(n) if (n %in% colnames(rt_md)) as.character(rt_md[[n]]) else NA_character_
    bt_col <- if ("transcript_biotype" %in% colnames(rt_md)) "transcript_biotype" else "transcript_type"

    ref_tx <- data.frame(
      gene_id            = rcol("gene_id"),
      gene_name          = rcol("gene_name"),
      transcript_id      = rcol("transcript_id"),
      transcript_name    = rcol("transcript_name"),
      transcript_biotype = rcol(bt_col),
      transcript_status  = "annotated",
      class_code         = NA_character_,
      classification     = NA_character_,
      ref_gene_id        = NA_character_,
      ref_gene_name      = NA_character_,
      ref_gene_biotype   = NA_character_,
      # These undetected reference transcripts are not novel models, so they have
      # no reference match and no Bambu class. The columns exist because `tx` has
      # them and rbind needs both frames to agree.
      ref_transcript_id   = NA_character_,
      ref_transcript_name = NA_character_,
      BambuTxClass        = NA_character_,
      chrom              = as.character(GenomeInfoDb::seqnames(rt)),
      strand             = as.character(BiocGenerics::strand(rt)),
      start              = BiocGenerics::start(rt),
      end                = BiocGenerics::end(rt),
      stringsAsFactors   = FALSE
    )

    ref_tx$label <- paste0(
      ifelse(!is.na(ref_tx$transcript_name) & nzchar(ref_tx$transcript_name),
             ref_tx$transcript_name, ref_tx$transcript_id),
      " | not detected",
      ifelse(ref_tx$strand %in% c("+", "-"), paste0(" | ", ref_tx$strand), ""))
    ref_tx$biotype_key <- ifelse(
      is.na(ref_tx$transcript_biotype) | !nzchar(ref_tx$transcript_biotype),
      "biotype not set", ref_tx$transcript_biotype)

    tx <- rbind(tx, ref_tx)
    cat(sprintf("Added %d annotated transcripts not present in the validated set\n",
                nrow(ref_tx)))
  }
}

# Which transcripts each panel is meant to contain -- this set, not the window alone,
# is what decides a panel's contents. Restricting it stops a neighbouring gene's
# isoforms from crowding out the ones the figure is about: a plain overlap query on
# one 84 kb window pulled in five genes and 31 transcripts to draw a caption that
# promised eight.
draw_ids <- character(0)

# What belongs in a panel is decided by the GENE and by the class codes, not by the
# window. The window is the gene padded by 5%, and that padding reaches into
# whatever sits next door.
#
# Three things qualify:
#
#   detected   every DETECTED transcript of the panel's genes
#   mine       novel models classified against one of those genes. Bambu gives each
#              intronic and antisense model its own gene id, so restricting to gene
#              id alone hides the pile -- at AK5, four single-exon antisense models
#              sit over one transcript, each in its own Bambu gene
#   referenced annotated transcripts, of any gene, that one of this panel's novel
#              models was classified against. Those are what a class code points
#              at, so a reader can check the code against the picture
#
# Nothing else undetected is drawn. Admitting a neighbour's transcripts because
# they fall inside the reference gene's span does not scale: at a 2 kb TEC locus
# overlapped by a 40-isoform lincRNA, all 40 qualified and the panel drew 2
# detected models beneath 40 undetected ones.
#
# Kept per panel rather than pooled, because a pooled set leaks between panels that
# happen to overlap.
panel_ids <- vector("list", nrow(cand))

if (nrow(cand)) {
  for (i in seq_len(nrow(cand))) {
    in_win <- tx$chrom == cand$chrom[i] &
              tx$end   >= cand$win_s[i] &
              tx$start <= cand$win_e[i]
    row_genes <- strsplit(cand$panel_genes[i], ";")[[1]]

    own  <- tx$gene_id %in% row_genes
    mine <- tx$transcript_status %in% "novel" &
              (tx$ref_gene_id %in% row_genes | own)

    # An undetected annotated transcript earns a row only when gffcompare
    # classified one of this panel's novel models against it. Those are the
    # transcripts a class code actually refers to, so a reader can check the code
    # against the picture; everything else undetected adds a row while saying
    # nothing the gene track above -- already the gene's full annotated footprint
    # -- does not show.
    #
    # The rule applies to every gene, not only the panel's own. Admitting a
    # neighbour's transcripts because they fall inside the reference gene's span
    # does not scale: at a 2 kb TEC locus overlapped by a 40-isoform lincRNA, all
    # 40 qualified and the panel drew 2 detected models under 40 undetected ones.
    ref_txids <- unique(tx$ref_transcript_id[in_win & mine])
    ref_txids <- bare_id(ref_txids[!is.na(ref_txids) & nzchar(ref_txids)])

    is_ann <- tx$transcript_status %in% "annotated"
    keep   <- in_win & ((own & !is_ann) | mine |
                        (is_ann & bare_id(tx$transcript_id) %in% ref_txids))

    panel_ids[[i]] <- unique(tx$transcript_id[keep])
    draw_ids <- c(draw_ids, panel_ids[[i]])
  }
}

flagged_ids <- vector("list", nrow(flagged))

if (nrow(flagged)) {
  for (i in seq_len(nrow(flagged))) {
    in_win <- tx$chrom == flagged$chrom[i] &
              tx$end   >= flagged$win_s[i] &
              tx$start <= flagged$win_e[i]
    # The host's isoforms for context, plus novel models classified against the host:
    # an intronic candidate carries its own Bambu gene id, so filtering on the host
    # gene alone would drop the transcript the figure exists to show.
    # %in% throughout rather than ==: ref_gene_id is NA for every known transcript,
    # and == would return NA, which indexes as NA and injects missing values.
    host <- tx$gene_id %in% flagged$ref_gene_id[i]
    mine <- tx$transcript_status %in% "novel" &
              (tx$ref_gene_id %in% flagged$ref_gene_id[i] | host)

    # Same rule as the stratified panels: an undetected annotated transcript is
    # drawn only where a novel model in this panel was classified against it.
    ref_txids <- unique(tx$ref_transcript_id[in_win & mine])
    ref_txids <- bare_id(ref_txids[!is.na(ref_txids) & nzchar(ref_txids)])

    is_ann <- tx$transcript_status %in% "annotated"
    keep   <- in_win & ((host & !is_ann) | mine |
                        (is_ann & bare_id(tx$transcript_id) %in% ref_txids))

    flagged_ids[[i]] <- unique(tx$transcript_id[keep])
    draw_ids <- c(draw_ids, flagged_ids[[i]])
  }
}
draw_ids <- unique(draw_ids[!is.na(draw_ids)])

# The validated and reference records carry different attribute sets, so both are
# cut down to the three fields the TxDb build actually reads before being combined.
# It also keeps the exported GTF small, which is the point of subsetting at all.
minimal_gtf <- function(gr) {
  out <- GenomicRanges::granges(gr)
  m   <- S4Vectors::mcols(gr)

  # phase rides along wherever the source has it. It is the reading frame of a CDS
  # row, and dropping it would export every CDS with a "." there. The drawing does
  # not use it -- plotTranscripts needs the ranges, not the frame -- but a GTF whose
  # CDS records have no phase is malformed, and makeTxDbFromGFF is entitled to say
  # so. Cheaper to keep than to find out which Bioconductor release starts caring.
  cols <- S4Vectors::DataFrame(
    type          = as.character(m$type),
    gene_id       = as.character(m$gene_id),
    transcript_id = as.character(m$transcript_id))
  if ("phase" %in% colnames(m)) cols$phase <- m$phase

  S4Vectors::mcols(out) <- cols
  out
}

sub_gtf <- minimal_gtf(IRanges::subsetByOverlaps(gtf, windows))
if (!is.null(ref_draw_gr)) {
  sub_gtf <- c(sub_gtf, minimal_gtf(ref_draw_gr))
}
row_tx_id <- as.character(S4Vectors::mcols(sub_gtf)$transcript_id)
# Gene-level rows carry no transcript_id and are kept regardless.
sub_gtf   <- sub_gtf[is.na(row_tx_id) | row_tx_id %in% draw_ids]
sub_tx_id <- as.character(S4Vectors::mcols(sub_gtf)$transcript_id)

# make.unique has to run over the transcript -> label mapping, one entry per
# transcript, NOT over GTF rows. sub_gtf carries one row per transcript AND one per
# exon, so uniquifying the row vector renames a transcript's own exons to "... #1",
# "... #2": they stop matching their parent transcript_id, every model loads into the
# TxDb with zero exons, and an 8-exon transcript draws as one featureless bar plus
# eight single-exon decoys. Uniquify the mapping, then apply it to every row.
uniq_tx        <- unique(sub_tx_id[!is.na(sub_tx_id)])
lab            <- tx$label[match(uniq_tx, tx$transcript_id)]
have_lab       <- !is.na(lab)
lab[have_lab]  <- make.unique(lab[have_lab], sep = " #")
lab[!have_lab] <- uniq_tx[!have_lab]
names(lab)     <- uniq_tx
relabel        <- unname(lab[sub_tx_id])
S4Vectors::mcols(sub_gtf)$transcript_id <- ifelse(is.na(relabel), sub_tx_id, relabel)

txdb_src <- file.path(opt$outdir, "genomic_context_regions.gtf")
rtracklayer::export(sub_gtf, txdb_src, format = "gtf")

cat(sprintf("Building TxDb from %d records across %d regions...\n",
            length(sub_gtf), nrow(cand)))
txdb <- if (requireNamespace("txdbmaker", quietly = TRUE)) {
  txdbmaker::makeTxDbFromGFF(txdb_src, format = "gtf")
} else {
  # makeTxDbFromGFF lived in GenomicFeatures before Bioconductor 3.19
  GenomicFeatures::makeTxDbFromGFF(txdb_src, format = "gtf")
}

pulposeq_assembly <- assembly(
  Genome         = "pulposeq_validated",
  TxDb           = txdb,
  OrgDb          = ORGDB,
  gene.id.column = "GENEID",
  display.column = "GENEID"
)

# --- Drawing ------------------------------------------------------------------

exons_gr  <- gtf[as.character(S4Vectors::mcols(gtf)$type) == "exon"]
exon_gid  <- as.character(S4Vectors::mcols(exons_gr)$gene_id)

#' Exon footprint for the gene track.
#'
#' The reference gene's annotated structure where the annotation has it, so the
#' track shows every exon a class code could have been assigned against -- not
#' only the exons of isoforms that happened to be detected. Falls back to the
#' validated GTF for a novel Bambu gene, which has no reference entry at all.
gene_structure <- function(gene_id) {
  key <- bare_id(gene_id)
  if (!is.na(key) && key %in% names(ref_exons_by_gene)) {
    return(ref_exons_by_gene[[key]])
  }
  GenomicRanges::reduce(exons_gr[exon_gid == gene_id])
}
exon_txid <- as.character(S4Vectors::mcols(exons_gr)$transcript_id)

#' Assign each model to a row, packing models that do not overlap onto the same one.
#'
#' A label sits above every model, so a row is claimed by whichever is wider, the model
#' or its label. The assignment returned here both draws the track and sizes it, which
#' is what makes the panel height exact. Sizing on the transcript count
#' alone leaves half the figure blank when models pack together; sizing on the bar
#' height alone clips the overflow, which is how a panel captioned "11 isoforms" came
#' to draw seven and lose three of its four novel candidates.
pack_rows <- function(starts, ends, labels, win_s, win_e, track_w, fontsize) {
  n <- length(starts)
  if (!n) return(integer(0))
  span    <- max(1, win_e - win_s)
  # Roughly the advance width of one character at this size, in inches, as bp.
  char_in <- fontsize * 0.0075
  lab_bp  <- nchar(labels) * char_in / track_w * span
  # Clamped to the window: a model running off the edge only competes for the space it
  # actually occupies on screen, and its label is centred on the visible part.
  vs  <- pmax(starts, win_s)
  ve  <- pmin(ends,   win_e)
  mid <- (vs + ve) / 2
  l   <- pmin(vs, mid - lab_bp / 2)
  r   <- pmax(ve, mid + lab_bp / 2)
  pad <- 0.02 * span
  row_end <- numeric(0)
  rows    <- integer(n)
  for (i in order(l)) {
    slot <- which(row_end + pad < l[i])
    k <- if (length(slot)) slot[1] else length(row_end) + 1L
    row_end[k] <- r[i]
    rows[i]    <- k
  }
  rows
}

#' A chromosome bar with the drawn window marked on it, and the zoom lines down to
#' the tracks.
#'
#' Deliberately not plotgardener's plotIdeogram(): that resolves cytobands through
#' AnnotationHub, which means a network call from a compute node that usually cannot
#' make one, and a failed run at the last step. This carries what the figure needs --
#' where on the chromosome the window sits -- and needs nothing but the chromosome
#' length, which the bigWig header already supplies.
draw_ideogram <- function(chrom, win_s, win_e, chrom_len,
                          x, y, w, h, panel_x0, panel_x1, zoom_y1) {
  plotRect(x = x, y = y, width = w, height = h, just = c("left", "top"),
           default.units = "inches", fill = "#F0F0F0", linecolor = "#BDBDBD")

  f0 <- max(0, min(1, win_s / chrom_len))
  f1 <- max(0, min(1, win_e / chrom_len))
  rx0 <- x + f0 * w
  # Widened to a floor, or a 20 kb window on a 250 Mb chromosome is a hairline.
  rx1 <- max(x + f1 * w, rx0 + 0.03)

  plotRect(x = rx0, y = y, width = rx1 - rx0, height = h, just = c("left", "top"),
           default.units = "inches", fill = "#C6E48B", linecolor = NA)
  plotSegments(x0 = (rx0 + rx1) / 2, y0 = y - 0.05,
               x1 = (rx0 + rx1) / 2, y1 = y + h + 0.05,
               default.units = "inches", linecolor = "#D62728", lwd = 1.2)
  plotText(label = sprintf("Chromosome %s", sub("^chr", "", chrom)),
           x = x + w, y = y + h + 0.07, just = c("right", "top"),
           fontsize = 8, fontcolor = "grey45", default.units = "inches")

  plotSegments(x0 = rx0, y0 = y + h + 0.03, x1 = panel_x0, y1 = zoom_y1,
               default.units = "inches", linecolor = "#CCCCCC", lty = 2, lwd = 0.8)
  plotSegments(x0 = rx1, y0 = y + h + 0.03, x1 = panel_x1, y1 = zoom_y1,
               default.units = "inches", linecolor = "#CCCCCC", lty = 2, lwd = 0.8)
}

#' Draw one region: chromosome position, gene structure, the isoform models over it,
#' and one coverage track per sample.
#'
#' Parameterised rather than branching internally, because the two callers differ
#' only in what they consider the structure row, which transcripts belong in the
#' panel, and what gets highlighted.
#'
#' @param structure_gr GRanges drawn as the single "gene" row above the isoforms
#' @param panel_tx     rows of `tx` this panel will draw; sizes the track and supplies
#'                     the biotype each model is coloured by
#' @param quantified   sample names that quantified the candidate this panel is about,
#'                     or NULL. Every sample's coverage is drawn either way -- seeing
#'                     the locus in the samples that did not call the model is the
#'                     point of the figure -- but the ones behind the read-level
#'                     numbers in the subtitle are named in bold.
draw_panel <- function(chrom, win_s, win_e, structure_gr, panel_tx,
                       title, subtitle, out_png, structure_label = "gene",
                       quantified = NULL) {

  # A shared y range across samples, so track heights compare directly instead of
  # each being rescaled to its own maximum.
  peaks <- vapply(bw_files, function(f) {
    sig <- rtracklayer::import.bw(
      f, which = GenomicRanges::GRanges(chrom, IRanges::IRanges(win_s, win_e)))
    if (length(sig) == 0) 0 else max(S4Vectors::mcols(sig)$score, na.rm = TRUE)
  }, numeric(1))
  ymax <- max(c(peaks, 1))

  margin   <- 0.5
  width    <- 8.0
  track_w  <- width - 2 * margin
  title_h  <- 0.34
  ideo_h   <- 0.17
  ideo_lab <- 0.22
  zoom_h   <- 0.40
  gene_h   <- 0.26
  # The row pitch is ours, not plotgardener's, because every model is placed here.
  row_h    <- TX_ROW_H
  sig_h    <- 0.55
  sig_lab_h <- 0.17        # band above each trace for its scale and sample name
  gap      <- 0.12
  label_h  <- 0.45
  legend_h <- 0.34

  # The same packing that draws the track sizes it, so the height is exact rather than
  # estimated -- no spare row, and nothing to clip.
  tx_rows <- pack_rows(panel_tx$start, panel_tx$end, panel_tx$label,
                       win_s, win_e, track_w, TX_LABEL_SIZE)
  n_rows  <- if (length(tx_rows)) max(tx_rows) else 1L
  # Exact, with no slack: every model sits at a row chosen here, so nothing can overflow
  # the track and nothing can be silently dropped from it.
  tx_h    <- max(0.4, n_rows * row_h)

  chrom_len <- if (!is.null(chrom_lens) && chrom %in% names(chrom_lens)) {
    as.numeric(chrom_lens[[chrom]])
  } else NA_real_
  show_ideo <- is.finite(chrom_len) && chrom_len > 0
  head_h    <- if (show_ideo) ideo_h + ideo_lab + zoom_h else 0

  # plotLegend lays its entries in a single row and divides the width equally between
  # them, so five entries including protein_coding_CDS_not_defined ran off the page.
  # Chunk into as many rows as the longest label needs, and draw one legend per row.
  bts       <- sort(unique(panel_tx$biotype_key))
  show_lgnd <- length(bts) > 0
  lgnd_rows <- if (show_lgnd) {
    # Calibrated against rendered panels: a 30-character biotype occupies ~1.4 in of
    # text at fontsize 7, so four entries fit across 7 in and five do not.
    per_row <- max(1L, floor(track_w / (max(nchar(bts)) * 0.048 + 0.30)))
    unname(split(bts, ceiling(seq_along(bts) / per_row)))
  } else list()
  lgnd_h    <- if (show_lgnd) length(lgnd_rows) * 0.20 + 0.14 else 0

  height  <- margin * 2 + title_h + head_h + gene_h + tx_h + gap * 3 +
             length(bw_files) * (sig_lab_h + sig_h + gap) + label_h + lgnd_h

  png(out_png, width = width, height = height, units = "in", res = 300)
  pageCreate(width = width, height = height, default.units = "inches",
             showGuides = FALSE)

  pars <- pgParams(
    chrom = chrom, chromstart = win_s, chromend = win_e,
    assembly = pulposeq_assembly,
    x = margin, width = width - 2 * margin, default.units = "inches"
  )

  # The region is named in full here, chromosome included, so a figure lifted out of
  # the report still says where it came from.
  plotText(
    label = sprintf("%s  |  %s:%s-%s", title, chrom,
                    format(win_s, big.mark = ","), format(win_e, big.mark = ",")),
    x = margin, y = margin, just = c("left", "top"),
    fontsize = 11, fontface = "bold", default.units = "inches"
  )
  plotText(
    label = subtitle,
    x = margin, y = margin + 0.17, just = c("left", "top"),
    fontsize = 8, fontcolor = "grey35", default.units = "inches"
  )

  y <- margin + title_h
  if (show_ideo) {
    draw_ideogram(chrom, win_s, win_e, chrom_len,
                  x = margin, y = y, w = track_w, h = ideo_h,
                  panel_x0 = margin, panel_x1 = margin + track_w,
                  zoom_y1  = y + ideo_h + ideo_lab + zoom_h - 0.04)
    y <- y + ideo_h + ideo_lab + zoom_h
  }

  if (length(structure_gr)) {
    plotRanges(
      params = pars, data = structure_gr,
      collapse = TRUE, fill = "#4D4D4D", linecolor = NA,
      y = y, height = gene_h
    )
    plotText(label = structure_label, x = margin - 0.06, y = y + gene_h / 2,
             just = c("right", "center"), fontsize = 7,
             fontcolor = "grey35", default.units = "inches")
  }

  y <- y + gene_h + gap

  # --- Isoform track ----------------------------------------------------------
  #
  # Colours come from `lab`, the same transcript -> label mapping the sub-GTF was
  # rewritten with, so a name that genuinely needed a uniquifying suffix still matches.
  # colorbyStrand is off: it would otherwise decide the colour of anything the
  # highlight table misses, and strand is already on the label.
  panel_lab <- unname(lab[panel_tx$transcript_id])
  panel_col <- biotype_color(panel_tx$biotype_key)

  # Genomic coordinate to page inches, clamped to the window.
  gx <- function(p) {
    margin + (min(max(p, win_s), win_e) - win_s) / max(1, win_e - win_s) * track_w
  }

  for (i in seq_len(nrow(panel_tx))) {
    lb <- panel_lab[i]
    if (is.na(lb)) next
    ry <- y + (tx_rows[i] - 1L) * row_h

    plotTranscripts(
      params = pars, y = ry + TX_LABEL_H, height = TX_MODEL_H,
      # Off, or plotgardener's own label lands on top of the one drawn below it and
      # every name appears twice.
      labels = NULL,
      # grid::unit rather than a bare unit(): grid is a base package and always
      # resolves, where relying on plotgardener to re-export it is an assumption with
      # no upside.
      boxHeight = grid::unit(TX_BOX_MM, "mm"), spaceHeight = TX_SPACE_H,
      stroke = TX_STROKE, limitLabel = FALSE,
      fill = c(BIOTYPE_FALLBACK, BIOTYPE_FALLBACK), colorbyStrand = FALSE,
      transcriptFilter     = lb,
      transcriptHighlights = data.frame(transcript = lb, color = panel_col[i],
                                        stringsAsFactors = FALSE)
    )

    # Centred on the visible extent, then clamped inside the track, so a model whose
    # midpoint falls outside the window still gets a name -- which is the whole reason
    # the labels are drawn here rather than left to plotTranscripts.
    lw <- nchar(lb) * TX_LABEL_SIZE * 0.0075
    lx <- min(max((gx(panel_tx$start[i]) + gx(panel_tx$end[i])) / 2,
                  margin + lw / 2),
              margin + track_w - lw / 2)
    plotText(label = lb, x = lx, y = ry, just = c("center", "top"),
             fontsize = TX_LABEL_SIZE, fontcolor = panel_col[i],
             default.units = "inches")
  }

  y <- y + tx_h + gap

  # Matched on the normalised name, because one side of the pair is a bigWig
  # basename and the other a BAM basename. Both descend from the same sample, but
  # only after the pipeline has renamed the file twice.
  norm <- function(x) gsub("[^a-z0-9]+", "", tolower(x))
  is_quant <- if (is.null(quantified) || !length(quantified)) {
    rep(NA, length(bw_files))
  } else {
    norm(sample_names) %in% norm(quantified)
  }

  for (i in seq_along(bw_files)) {
    # Scale and sample name in a band above the trace rather than inside it.
    # plotSignal's own scale = TRUE / label = draw both over the plotting area, where a
    # peak reaching ymax runs straight through them: "[0 - 121]" came out with the 1
    # struck through by the trace.
    # Named, not just bracketed: a bare "[0 - 130]" beside a transcript track reads as
    # though it might be a length in bp. It is the vertical scale -- reads stacked over
    # a single base -- and saying so costs five characters.
    # Rounded because bigWig scores are numeric: an unrounded ymax prints as 130.457.
    plotText(label = sprintf("depth [0 - %s]", format(round(ymax), scientific = FALSE,
                                                      trim = TRUE)),
             x = margin, y = y, just = c("left", "top"),
             fontsize = 7, fontcolor = "grey35", default.units = "inches")
    # Bold and black where the sample quantified the candidate, grey where it did
    # not. The trace itself is drawn identically in both, so the comparison the
    # figure exists to support is unaffected.
    plotText(label = sample_names[i], x = margin + track_w, y = y,
             just = c("right", "top"), fontsize = 8,
             fontface  = if (isTRUE(is_quant[i])) "bold" else "plain",
             fontcolor = if (isFALSE(is_quant[i])) "grey55" else "#1A1A1A",
             default.units = "inches")
    plotSignal(
      params = pars, data = bw_files[i],
      y = y + sig_lab_h, height = sig_h,
      range = c(0, ymax),
      linecolor = "#3B6FB6", fill = "#3B6FB6"
    )
    y <- y + sig_lab_h + sig_h + gap
  }

  plotGenomeLabel(
    chrom = chrom, chromstart = win_s, chromend = win_e,
    assembly = pulposeq_assembly,
    x = margin, y = y, length = track_w,
    default.units = "inches", scale = "bp", fontsize = 8
  )

  # Only the biotypes this panel contains, so the key never explains a colour that is
  # not on screen.
  if (show_lgnd) {
    per_row <- length(lgnd_rows[[1]])
    ly <- y + label_h
    for (r in lgnd_rows) {
      # Width proportional to the entries on this row, so a short final row keeps the
      # same swatch spacing as a full one instead of stretching across the page.
      plotLegend(
        legend = r, fill = biotype_color(r), border = FALSE,
        orientation = "h", fontsize = 7,
        x = margin, y = ly, width = track_w * length(r) / per_row, height = 0.18,
        just = c("left", "top"), default.units = "inches"
      )
      ly <- ly + 0.20
    }
  }

  pageGuideHide()
  dev.off()
  basename(out_png)
}

# --- Genes carrying novel isoforms --------------------------------------------

if (nrow(cand)) {
  cand$figure      <- NA_character_
  cand$n_annotated <- 0L
  for (i in seq_len(nrow(cand))) {
    row <- cand[i, ]

    # The set chosen for this panel specifically. draw_ids is the pooled union used
    # to size the shared TxDb; filtering on that instead would let an overlapping
    # panel's transcripts in.
    panel_tx <- tx[tx$transcript_id %in% panel_ids[[i]], , drop = FALSE]

    # Counted from what is actually drawn, not from the genes the panel was built
    # around. The caption has to describe the picture.
    n_drawn <- nrow(panel_tx)
    n_nov   <- sum(panel_tx$transcript_status == "novel",     na.rm = TRUE)
    n_ann   <- sum(panel_tx$transcript_status == "annotated", na.rm = TRUE)
    n_kno   <- n_drawn - n_nov - n_ann
    cand$n_tx[i]        <- n_drawn
    cand$n_novel[i]     <- n_nov
    cand$n_known[i]     <- n_kno
    cand$n_annotated[i] <- n_ann

    cat(sprintf("Drawing [%s] %s at %s:%d-%d (%d drawn: %d detected, %d novel, %d annotated only)\n",
                row$stratum, row$label, row$chrom, row$win_s, row$win_e,
                n_drawn, n_kno, n_nov, n_ann))

    ref_bt <- if (is.na(row$ref_gene_biotype) || !nzchar(row$ref_gene_biotype)) {
      row$ref_class
    } else row$ref_gene_biotype

    cand$figure[i] <- draw_panel(
      chrom        = row$chrom, win_s = row$win_s, win_e = row$win_e,
      # The gene's full ANNOTATED exon footprint, so every exon a class code could
      # have been assigned against is on screen -- including those of isoforms not
      # detected in these samples, which the isoform rows below cannot show.
      structure_gr = gene_structure(row$gene_id),
      panel_tx     = panel_tx,
      title        = row$label,
      # The stratum is the point of the panel, so it leads: this figure is here as
      # the example of its kind, not because the locus is remarkable.
      subtitle     = sprintf("%s, %s | reference: %s%s",
                             row$novel_biotype, row$classification, ref_bt,
                             if (is.na(row$BambuTxClass) || !nzchar(row$BambuTxClass)) ""
                               else sprintf(" | Bambu: %s", row$BambuTxClass)),
      out_png      = file.path(opt$outdir,
                               sprintf("genomic_context_%s_%s.png",
                                       slug(row$label), slug(row$stratum)))
    )
  }

  write.csv(cand[, c("stratum", "novel_biotype", "classification", "class_code",
                     "ref_class", "ref_gene_biotype", "BambuTxClass", "example_tx",
                     "gene_id", "gene_name", "label", "panel_genes", "chrom",
                     "start", "end", "win_s", "win_e", "n_tx", "n_known",
                     "n_novel", "n_annotated", "n_novel_lnc", "biotypes", "figure")],
            file.path(opt$outdir, "genomic_context_candidates.csv"), row.names = FALSE)
} else {
  # Same columns as the populated case. The report reads this file by name, and a
  # zero-row table with a different header is harder to handle than an empty one
  # with the right header.
  empty_cols <- c("stratum", "novel_biotype", "classification", "class_code",
                  "ref_class", "ref_gene_biotype", "BambuTxClass", "example_tx", "gene_id",
                  "gene_name", "label", "panel_genes", "chrom", "start", "end",
                  "win_s", "win_e", "n_tx", "n_known", "n_novel", "n_annotated", "n_novel_lnc",
                  "biotypes", "figure")
  write.csv(as.data.frame(setNames(replicate(length(empty_cols), character(0),
                                             simplify = FALSE), empty_cols)),
            file.path(opt$outdir, "genomic_context_candidates.csv"), row.names = FALSE)
}
cat("Wrote genomic_context_candidates.csv\n")

# --- Sense-intronic candidates ------------------------------------------------

if (nrow(flagged)) {
  flagged$figure <- NA_character_
  for (i in seq_len(nrow(flagged))) {
    row <- flagged[i, ]
    cat(sprintf("Drawing sense-intronic %s in %s intron at %s:%d-%d\n",
                row$qry_id, ifelse(is.na(row$ref_gene_name), row$ref_gene_id,
                                   row$ref_gene_name),
                row$chrom, row$win_s, row$win_e))

    # Everything in the window, not just the host's isoforms: an intronic candidate
    # carries its own Bambu gene id, so filtering on the host gene would leave other
    # novel models in the window out of the panel entirely. The candidate is in here
    # too, and appears only here, where the track label carries its identifier, class
    # code wording and strand.
    panel_tx <- tx[tx$transcript_id %in% flagged_ids[[i]], , drop = FALSE]

    host_lab <- if (is.na(row$ref_gene_name) || !nzchar(row$ref_gene_name)) {
      row$ref_gene_id
    } else {
      row$ref_gene_name
    }

    flagged$figure[i] <- draw_panel(
      chrom        = row$chrom, win_s = row$win_s, win_e = row$win_e,
      structure_gr = gene_structure(row$ref_gene_id),
      panel_tx     = panel_tx,
      title        = sprintf("%s in %s", row$qry_id, host_lab),
      subtitle     = sprintf(
        "sense intronic (%s) | %s | %s%s | %s counts%s",
        ifelse(is.na(row$ref_gene_biotype), "unknown biotype", row$ref_gene_biotype),
        if (isTRUE(row$num_exons == 1L)) "mono-exonic"
          else sprintf("%d exons", row$num_exons),
        if (isTRUE(row$full_length_support)) "full-length support"
          else if (isFALSE(row$full_length_support)) "no full-length support"
          else "full-length support not recorded",
        if (is.na(row$BambuTxClass) || !nzchar(row$BambuTxClass)) ""
          else sprintf(" | Bambu: %s", row$BambuTxClass),
        format(round(row$counts_total), big.mark = ","),
        if (is.na(row$samples_quantified)) "" else
          sprintf(" in %d/%d samples",
                  row$samples_quantified, row$samples_total)),
      out_png      = file.path(opt$outdir,
                               sprintf("intronic_context_%s.png", row$qry_id)),
      structure_label = "host",
      # Bold marks which samples the count came from, so it says nothing when every
      # sample quantified the candidate.
      quantified   = if (is.na(row$quantified_samples) ||
                         !nzchar(row$quantified_samples) ||
                         isTRUE(row$samples_quantified == row$samples_total)) NULL else
                       trimws(strsplit(row$quantified_samples, ";", fixed = TRUE)[[1]])
    )
  }

  write.csv(flagged[, c("qry_id", "ref_gene_id", "ref_gene_name", "ref_gene_biotype",
                        "chrom", "start", "end", "win_s", "win_e", "strand",
                        "class_code", "BambuTxClass", "num_exons",
                        "samples_quantified", "samples_total",
                        "quantified_samples", "counts_total",
                        "full_length_support",
                        "figure")],
            file.path(opt$outdir, "intronic_context_candidates.csv"), row.names = FALSE)
} else {
  write.csv(data.frame(), file.path(opt$outdir, "intronic_context_candidates.csv"),
            row.names = FALSE)
}
cat("Wrote intronic_context_candidates.csv\n")
