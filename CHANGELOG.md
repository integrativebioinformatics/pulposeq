# integrativebioinformatics/pulposeq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0 - [date]

Initial release of integrativebioinformatics/pulposeq, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- `PREFILTER` step (`bin/prefilter_transcripts.R`): the one place transcripts are filtered, right after Bambu and before classification. It removes unstranded novel models, novel models shorter than 500 nt, and any transcript whose counts summed across all samples fall below the new `--min_total_counts` parameter (default 5). It writes a per-rule `filter_summary.tsv`.
- `--min_total_counts` parameter.
- Optional oarfish quantification (`--run_oarfish`, `--oarfish_num_bootstraps`), using the nf-core `oarfish/readmode` module. oarfish quantifies every Bambu model except unstranded novel ones, on the same reads minimap2 aligns, with bootstrap replicates. When it runs, its counts replace Bambu's as the expression source for `PREFILTER`. New local modules `OARFISH_INDEX_FASTA` and `OARFISH_MERGE` (`bin/oarfish_merge.R`), and an `oarfish/` output folder with per-sample results and `oarfish_summary.tsv`.
- Resource settings for the `process_high` label, used by the nf-core oarfish module, which pulposeq did not define before.
- `CURATE` step (`bin/curate_transcripts.R`): writes the final transcript set. Every curated table and GTF, Bambu's and oarfish's, holds exactly the same transcripts. Adds `transcript_read_support.tsv` (Bambu total, full-length and unique reads per transcript, `full_length_supported`, and the oarfish total when oarfish ran) and, with oarfish, curated count, CPM, TPM (PacBio) and gene tables plus a `SummarizedExperiment` with bootstrap replicates (`infRep1..N`, fishpond layout).
- `oarfish_summary.tsv` compares oarfish with Bambu as two labelled diagnostics: detection overlap (detected by each, by both, Jaccard) and the Spearman correlation over transcripts both tools detected. It warns only for clearly broken runs.

### `Changed`

- The expression filter now keeps a transcript only if its counts summed across all samples reach `--min_total_counts` (default 5), replacing "non-zero in at least one sample". This changes results: fewer lowly expressed transcripts are kept.
- Classification (gffread, gffcompare, CPC2) now runs only on the novel models `PREFILTER` keeps, instead of every novel model.
- `novel_transcripts.R` no longer filters on strand, length or counts; it applies class codes and assigns biotypes.
- New output layout: everything from Bambu onwards is published under `transcriptome/`, in `bambu/` (raw), `oarfish/` (raw), `classification/` and `curated/` (the results, in `annotation/`, `counts/`, `metadata/` and `plots/`). Replaces `bambu/`, `bambu_filter/`, `bambu_validated/`, `gffcompare/`, `gffread/`, `coding_potential/`, `novel_transcripts/` and `ref_transcripts/`. The expression filter's intermediate `*_filter.txt` tables are no longer published.

### `Fixed`

### `Dependencies`

### `Deprecated`

### `Removed`

- `FILTER_BAMBU_COUNTS` module and `bin/filter_bambu_counts.sh`, replaced by `PREFILTER`.
- `VALIDATE_BAMBU_COUNTS` and `VALIDATE_BAMBU_GTF` modules, and `bin/validate_bambu_counts.sh`, `bin/validate_bambu_gtf.sh` and `bin/validate_bambu_gtf.awk`, replaced by `CURATE`, which writes byte-identical tables and GTFs.

- RNAmining coding-potential predictor, together with the `coding_potential_pred` and `organism` parameters. CPC2 is now the only coding-potential predictor.
