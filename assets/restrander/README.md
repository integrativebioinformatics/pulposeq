# Custom Restrander configurations

`--restrand_kit` covers the protocols Restrander ships presets for: `PCB109`, `PCB111`,
`PCB114`, `DCS109`, `DCS-LSK114`, `NEBNext` and `trimmed`. For anything else, write a
configuration and pass it with `--restrand_config`, which overrides `--restrand_kit`:

```bash
--restrand_config /path/to/my_protocol.json
```

`custom_config_example.json` in this directory is a working starting point. Copy it,
replace the primer sequences, and point `--restrand_config` at your copy. JSON has no
comment syntax, so the field reference is here rather than inline.

Whichever configuration a run resolves to is copied into `results/restrander/` as
`<sample>.restrander_config.json`, so the primer sequences behind every orientation call
stay with the results.

## Field reference

| Field | Meaning |
|---|---|
| `name` | Label for the configuration. Appears in Restrander's own output. |
| `description` | Free text. Worth writing properly — it is the only place the protocol is recorded. |
| `pipeline` | Ordered list of classification stages. Each read passes through them in turn and stops at the first that resolves its orientation. |
| `silent` | Suppresses Restrander's progress messages. |
| `exclude-unknowns` | Whether reads that could not be oriented are written to a separate `-unknowns.fastq.gz` instead of the main output. **Leave this `true`** — see below. |
| `error-rate` | Proportion of mismatches tolerated when matching a primer. `0.25` is the value every shipped preset uses. |

### The `poly` stage

Classifies on polyA/polyT tails.

| Field | Meaning |
|---|---|
| `tail-length` | Minimum run of A or T to count as a tail. |
| `search-size` | How many bases from each read end to search. |

### The `primer` stage

Classifies on the protocol's primer sequences. This is the stage to edit.

| Field | Meaning |
|---|---|
| `tso` | Strand-switching primer, also called SSP. Marks the 5' end of the transcript. |
| `rtp` | Reverse-transcription primer, also called VNP. Marks the 3' end. |
| `report-artefacts` | Detects reads with the wrong combination of ends — `TSO-TSO`, `RTP-RTP` — which indicate chimeras or mispriming. |

## Two things worth knowing

**Keep `exclude-unknowns` set to `true`.** With it enabled, reads whose orientation could
not be determined go to `<sample>.restranded-unknowns.fastq.gz` rather than the main
output. That matters because everything downstream then treats the reads as genuinely
oriented: minimap2 runs with `-uf` and Bambu with `stranded = TRUE`. A read with an
undetermined orientation surviving into that path would be assigned a coin-flip strand
that the rest of the pipeline is instructed to trust. The unknowns file is published and
its read count reported, so nothing disappears silently.

**A `primer` stage is what makes artefact detection possible.** `report-artefacts` lives
inside it, so a poly-only configuration — `trimmed`, or a custom one with the primer stage
removed — reports no artefact classes. It will also orient a smaller share of reads,
because a read has to carry a detectable polyA or polyT tail to be classified at all. If
your reads have primers, keep the primer stage.
