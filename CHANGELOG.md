# integrativebioinformatics/pulposeq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0 - [date]

Initial release of integrativebioinformatics/pulposeq, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- `PREFILTER` step (`bin/prefilter_transcripts.R`): the one place transcripts are filtered, right after Bambu and before classification. It removes unstranded novel models, novel models shorter than 500 nt, and any transcript whose counts summed across all samples fall below the new `--min_total_counts` parameter (default 5). It writes a per-rule `filter_summary.tsv`.
- `--min_total_counts` parameter.

### `Changed`

- The expression filter now keeps a transcript only if its counts summed across all samples reach `--min_total_counts` (default 5), replacing "non-zero in at least one sample". This changes results: fewer lowly expressed transcripts are kept.
- Classification (gffread, gffcompare, CPC2) now runs only on the novel models `PREFILTER` keeps, instead of every novel model.
- `novel_transcripts.R` no longer filters on strand, length or counts; it applies class codes and assigns biotypes.

### `Fixed`

### `Dependencies`

### `Deprecated`

### `Removed`

- `FILTER_BAMBU_COUNTS` module and `bin/filter_bambu_counts.sh`, replaced by `PREFILTER`.

- RNAmining coding-potential predictor, together with the `coding_potential_pred` and `organism` parameters. CPC2 is now the only coding-potential predictor.
