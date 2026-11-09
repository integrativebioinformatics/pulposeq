# integrativebioinformatics/pulposeq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0 - [date]

Initial release of integrativebioinformatics/pulposeq, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- `PREFILTER` step (`bin/prefilter_transcripts.R`): the one place transcripts are filtered, right after Bambu and before classification. It removes unstranded novel models, novel models shorter than 500 nt, and any transcript whose counts summed across all samples fall below the new `--min_total_counts` parameter (default 5). It writes a per-rule `filter_summary.tsv`.
- `--min_total_counts` parameter.
- Optional oarfish quantification (`--run_oarfish`, `--oarfish_num_bootstraps`), using the nf-core `oarfish/readmode` module. oarfish quantifies every Bambu model except unstranded novel ones, on the same reads minimap2 aligns, with bootstrap replicates. When it runs, its counts replace Bambu's as the expression source for `PREFILTER`. New local modules `OARFISH_INDEX_FASTA` and `OARFISH_MERGE` (`bin/oarfish_merge.R`), and an `oarfish/` output folder with per-sample results and `oarfish_summary.tsv`.
- Resource settings for the `process_high` label, used by the nf-core oarfish module, which pulposeq did not define before, checked against a test-run trace. `CURATE` gets 8 GB when oarfish runs (it loads the bootstrap replicates) and keeps 2 GB otherwise. The `medium` and `large` profiles get their own oarfish settings (`OARFISH_READMODE` 16 cpus/32 GB/6 h and 32 cpus/96 GB/12 h; `OARFISH_MERGE` and, with oarfish, `CURATE` 32 GB and 128 GB), sized from the oarfish paper's runtime and memory figures; without oarfish, `CURATE` keeps each profile's `process_low` memory.
- `CURATE` step (`bin/curate_transcripts.R`): writes the final transcript set. Every curated table and GTF, Bambu's and oarfish's, holds exactly the same transcripts. Adds `transcript_read_support.tsv` (Bambu total, full-length and unique reads per transcript, `full_length_supported`, and the oarfish total when oarfish ran) and, with oarfish, curated count, CPM, TPM (PacBio) and gene tables plus a `SummarizedExperiment` with bootstrap replicates (`infRep1..N`, fishpond layout).
- `oarfish_summary.tsv` compares oarfish with Bambu as two labelled diagnostics: detection overlap (detected by each, by both, Jaccard) and the Spearman correlation over transcripts both tools detected. It warns only for clearly broken runs.
- Report section *Bambu vs oarfish* (only when oarfish ran): per-sample agreement table, curated transcripts by which tool gave them reads (including those kept on oarfish's counts alone), a Bambu-vs-oarfish scatter, and quantification uncertainty from the bootstraps (mean InfRV per transcript, added to `transcript_read_support.tsv` as `oarfish_mean_infrv`).
- The report's run-parameters box states the expression source and `--min_total_counts`, and shows `filter_summary.tsv` per rule.

### `Changed`

- The expression filter now keeps a transcript only if its counts summed across all samples reach `--min_total_counts` (default 5), replacing "non-zero in at least one sample". This changes results: fewer lowly expressed transcripts are kept.
- Classification (gffread, gffcompare, CPC2) now runs only on the novel models `PREFILTER` keeps, instead of every novel model.
- `novel_transcripts.R` no longer filters on strand, length or counts; it applies class codes and assigns biotypes.
- Sample PCA and expression heatmaps are transcript-level only, drawn for every quantified model and for the curated set from the expression source: oarfish's counts when it ran, Bambu's otherwise. The gene-level heatmap is no longer drawn after curation. `VALIDATE_BAMBU_OUTPUTS` now writes `pca_all_models.png`, `pca_grouped_all_models.png` and `heatmap_transcript_all_models.png` alongside the `_validated` versions.
- `validation_summary.csv` measures the yield against the expression source and gains an `expression_source` column. With oarfish, curated transcripts that Bambu gave no reads are now counted instead of dropped with a warning.
- New output layout: everything from Bambu onwards is published under `transcriptome/`, in `bambu/` (raw), `oarfish/` (raw), `classification/` and `curated/` (the results, in `annotation/`, `counts/`, `metadata/` and `plots/`). Replaces `bambu/`, `bambu_filter/`, `bambu_validated/`, `gffcompare/`, `gffread/`, `coding_potential/`, `novel_transcripts/` and `ref_transcripts/`. The expression filter's intermediate `*_filter.txt` tables are no longer published.

### `Fixed`

- `VALIDATE_BAMBU_OUTPUTS` memory raised from 2 GB to 8 GB: it peaked at 1.7 GB on a chromosome 1 test run, leaving no headroom for a whole-genome annotation.
- Report: the novel read-support fraction figure joined the full-length and unique tables with an inner join. Those tables list only transcripts with at least one such read, so novel models with no full-length or no unique read were left out of the figure, hiding its low tail. They now count as 0.

### `Dependencies`

### `Deprecated`

### `Removed`

- `FILTER_BAMBU_COUNTS` module and `bin/filter_bambu_counts.sh`, replaced by `PREFILTER`.
- `VALIDATE_BAMBU_COUNTS` and `VALIDATE_BAMBU_GTF` modules, and `bin/validate_bambu_counts.sh`, `bin/validate_bambu_gtf.sh` and `bin/validate_bambu_gtf.awk`, replaced by `CURATE`, which writes byte-identical tables and GTFs.

- RNAmining coding-potential predictor, together with the `coding_potential_pred` and `organism` parameters. CPC2 is now the only coding-potential predictor.
