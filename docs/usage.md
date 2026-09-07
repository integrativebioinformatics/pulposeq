# Usage

> *Pipeline parameters are described below and should match the current pipeline schema/configuration.*

## Introduction

<!-- TODO integrative/bioinformatics: Add documentation about anything specific to running your pipeline. For general topics, please point to (and add to) the main nf-core website. -->

## Scope and assumptions

Two boundaries determine what pulposeq can find. Neither is a bug, and both are
easier to plan around than to discover in the results.

### Eukaryotic organisms

Splice-aware alignment, the intron-based gffcompare class codes the classification
depends on, and the coding-potential models all assume a eukaryotic
transcriptome. The pipeline has no meaningful behaviour on prokaryotic data.

Reference annotations must come from **Ensembl or GENCODE**. Other sources use
different attribute names for gene and transcript biotypes, and the pipeline reads
those attributes directly.

### Polyadenylated transcriptomes

This is a property of the **library**, not of the pipeline, but it decides what
can be discovered — so it belongs here rather than in the results.

Every standard long-read RNA kit selects for poly(A):

| Method | How it engages the transcript |
|---|---|
| ONT PCR-cDNA / direct cDNA | RT adapter with a poly-T overhang, or oligo-dT |
| ONT direct RNA (RNA004) | adapter ligated onto the poly(A) tail |
| PacBio Iso-Seq / Kinnex | oligo-dT primed reverse transcription |

Reverse transcriptase extends from a primer annealed near the template's 3' end,
so it needs a 3' handle, and in every standard kit the poly(A) tail is that handle.

Two consequences follow:

- **Non-polyadenylated transcripts are under-represented.** This includes a
  structurally distinct and well-characterised class of lncRNA — `MALAT1` and
  `NEAT1_2` are both triple-helix terminated rather than polyadenylated, as are
  replication-dependent histone mRNAs and sno-lncRNAs. Ribo-depleted short-read
  data is the usual way to find out what a poly(A)-selected library missed.
- **Oligo-dT can bind genomic A-runs, not only poly(A) tails.** When it does,
  reverse transcription starts from an internal position and the resulting cDNA
  is a truncated fragment that appears to end inside an intron. This is *internal
  priming*, and it produces exactly the profile of a candidate intronic lncRNA:
  unspliced, intronic, no ORF, a few hundred bases to a few kb. The pipeline does
  not filter the affected population, because filtering would remove real
  transcripts alongside artifacts. It surfaces it instead: such models are counted
  as mono-exonic sense-intronic candidates in the report, and each of the ones
  drawn is shown over its whole host intron in the genomic context panels, where
  internal priming reads as coverage continuous with the host rather than as a
  discrete block.

Neither constraint is fundamental to long-read sequencing. Protocols exist that
capture 3' ends independently of polyadenylation, and being free of oligo-dT they
are also free of internal priming — but they are published protocols rather than
off-the-shelf kits, and the pipeline assumes a standard library.

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use the parameter `--input` in the bash command to specify its location. It has to be a comma-separated file with 3 columns, and a header row as shown in the examples below.

``` bash
--input '[path to samplesheet file]'
```

### Full samplesheet

A final `samplesheet.csv` file consisting of single-end data may look something like the one below. This is for 6 samples, where we have 2 experimental groups and 3 replicates per group.

``` csv
sample,group,fastq
R1_H1975,H1975,home/user/R1_H1975.fastq.gz
R2_H1975,H1975,home/user/R2_H1975.fastq.gz
R3_H1975,H1975,home/user/R3_H1975.fastq.gz
R1_HCC827,HCC827,home/user/R1_HCC827.fastq.gz
R2_HCC827,HCC827,home/user/R2_HCC827.fastq.gz
R3_HCC827,HCC827,home/user/R3_HCC827.fastq.gz
```

| Column | Description |
|-------------|-----------------------------------------------------------|
| `sample` | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `group` | Experimental group name. For example: `treatment` vs `control` or `cell_line1` vs `cell_line2` |
| `fastq` | Full path to FastQ file for ONT or PacBio long-reads. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz". |

Another [example samplesheet] has been provided with the pipeline.

  [example samplesheet]: ../assets/samplesheet.csv

## Pipeline parameters

After seeting up the samplesheet, follow up to set the pipeline parameters.

> [!WARNING]
> Relative paths might cause issues! Always input double-quoted **FULL PATHS** (e.g. `"/full/path/to/file/file.extension"`).

> [!TIP]
> Stay up to date! Remember to always use the **lastest Ensembl releases** for reference genomes and annotations.

| Parameter | Description |
|----------------|--------------------------------------------------------|
| `input` | Full path to the `samplesheet.csv` file |
| `outdir` | Full path to the output directory where you want the results to be saved |
| `skip_nanocomp` | `NanoComp` reports to skip: `raw`, `restrander`, `chopper`, `minimap2` or `all`, comma-separated (see [below](#choosing-which-nanocomp-reports-to-run)) |
| `skip_chopper` | Skip filtering and trimming with `chopper`; the reads reach alignment exactly as supplied — or as Restrander left them, for unoriented ONT cDNA (`library: ONT_cDNA` with `stranded_library: false`) |
| `minqual` | Minimum read average Phred-score quality cut-off |
| `minlen` | Minimum read length (bp) |
| `maxlen` | Maximum read length (bp), optional — leave unset for no upper limit |
| `maxgc` | Maximum GC content (%) |
| `mingc` | Minimum GC content (%) |
| `headcrop` | Cut `x` number of bases from the beginning of the reads |
| `tailcrop` | Cut `y` number of bases from the end of the reads |
| `skip_minimap2` | Skip mapping/alignment to genome reference, and with it every later stage (runs MultiQC, then finishes pipeline execution) |
| `skip_bambu` | Skip transcriptome assembly with `bambu`, and with it classification, metadata refinement and the report |
| `reference` | Full path to a reference genome `FASTA` file from Ensembl |
| `annotation` | Full path to a reference annotation `GTF` file from Ensembl |
| `library` | Sequencing library type: `ONT_cDNA`, `ONT_DRS` or `PacBio` |
| `stranded_library` | Whether the reads **on disk** are already oriented (`true`/`false`) |
| `map_hq` | Use minimap2's `splice:hq` preset instead of `splice`. Unset follows the library default; not valid for `ONT_DRS` |
| `restrand_kit` | Sequencing kit whose primers Restrander should look for — **required** for unoriented ONT cDNA |
| `restrand_config` | Optional custom Restrander configuration JSON; overrides `restrand_kit` |
| `restrand_min_frac` | Minimum oriented fraction, per sample; below it the run stops (default `0.80`, minimum `0.5`) |
| `skip_class` | Skip transcriptome characterization; `bambu` still runs (runs MultiQC, then finishes pipeline execution) |
| `ndr` | Bambu's novel discovery rate. Leave unset for Bambu's own automatic choice (see [below](#novel-discovery-rate)) |
| `coding_potential_pred` | Coding-potential predictor: `cpc2` (default) or `rnamining` |
| `organism` | Organism scientific name (e.g. `"Homo_sapiens"`). Required by `rnamining`, ignored by `cpc2` |

### Choosing which NanoComp reports to run

NanoComp runs at up to four points, and `skip_nanocomp` names the ones to drop.
One parameter rather than a flag per report: they are the same tool answering the
same question at four moments of the same run, and what a run actually decides is
where read-level QC is worth its runtime.

| Value | Report skipped | Runs on |
|---|---|---|
| unset (default) | none — all four run | |
| `raw` | `NANOCOMP_RAW` | the reads as supplied |
| `restrander` | `NANOCOMP_REST` | the oriented reads (only exists when restranding runs) |
| `chopper` | `NANOCOMP_CHOPPER` | the filtered reads |
| `minimap2` | `NANOCOMP_MINIMAP2` | the alignments |
| `all` | all four | |

Values combine, comma-separated — `'raw,chopper'`, `'rest,minimap2'`,
`'raw,rest,chopper'`. `rest` is accepted as a synonym for `restrander`.

```bash
nextflow run integrativebioinformatics/pulposeq -profile docker --input samplesheet.csv --outdir results --skip_nanocomp 'raw,chopper'
```

Two things it deliberately does not do:

- **It never turns off the tool being reported on.** `--skip_nanocomp chopper` still
  filters; `--skip_nanocomp restrander` still orients. To skip the filtering itself,
  use `skip_chopper`, and note that this takes its report with it — there is nothing
  left to report on. Restranding has no skip at all (see
  [Restranding ONT cDNA libraries](#restranding-ont-cdna-libraries)).
- **It is not a resource dial in disguise.** `minimap2` is the report most worth
  dropping on large runs, since it is the one that reads whole BAMs — see
  [Measured resource requirements](#measured-resource-requirements) for the numbers.

The `restrander` report exists only for ONT cDNA the pipeline has to orient itself.
Naming it for a library that is already oriented is not an error — the run proceeds
and logs a warning saying the token had nothing to act on.

### Novel discovery rate {#novel-discovery-rate}

`ndr` is Bambu's threshold for calling a read class a novel transcript, and it has
the single largest effect on how many novel candidates exist. Leave it unset and
Bambu chooses one itself from the data; set it to a number between 0 and 1 to fix it.

Bambu reports an automatically chosen NDR only on stdout, so the pipeline parses it
back out and writes it to `bambu/bambu_run_metrics.csv` alongside what was requested.
The report prints both at the top of **Overall Assembly Statistics**, so a figure is
never read without knowing the threshold behind it.

A lower NDR admits fewer, better-supported novel models; a higher one admits more and
shifts the burden onto the evidence reported per candidate. Neither is a default worth
recommending sight-unseen — if you fix it, record why.

### Coding potential

Novel transcripts are evaluated by a
coding-potential call, so the predictor decides what the pipeline's main output
contains.

**`cpc2` is the default and needs no organism.**

`rnamining` remains selectable and requires `organism`, which picks a per-organism
model; the run stops rather than passing an empty one. 

Whichever runs, its output is published to `coding_potential/`, and every novel
record carries a `coding_predictor` column naming the tool and a `coding_prob`
column holding P(coding). Those are normalised across the two, because the tools
report different quantities natively — RNAmining scores the class it chose, CPC2
scores coding specifically.

### Library type and strandedness

The `library` and `stranded_library` parameters are **required** whenever alignment is not skipped. Together they set the `minimap2` alignment preset and Bambu's `stranded` argument — the two are deliberately driven from a single declaration of library chemistry, so they cannot be configured independently.

> [!IMPORTANT]
> `stranded_library` describes **the reads on disk**, not the library chemistry. ONT PCR-cDNA chemistry *is* strand-specific, but its basecalled reads are *not* oriented, because the protocol sequences either cDNA strand at random. Conflating the two is the most common way to get this wrong.

> [!IMPORTANT]
> **The pipeline processes stranded data only.** There is no unstranded path. Unoriented reads reaching alignment and quantification erase mono-exonic novel transcripts, push novel isoforms of known genes into the antisense class, and bleed reads between overlapping sense/antisense pairs — the three things an isoform-level lncRNA pipeline is least able to absorb. ONT cDNA is therefore either oriented before it arrives, or oriented by the pipeline.

| `library` | `stranded_library` | Restrander | `minimap2` preset (`map_hq` unset) | Bambu `stranded` |
|---|---|---|---|---|
| `ONT_DRS` | automatically `true` | no | `-ax splice -uf -k14` | `TRUE` |
| `PacBio` | automatically `true` | no | `-ax splice:hq -uf` | `TRUE` |
| `ONT_cDNA` | `true` | no | `-ax splice -uf` | `TRUE` |
| `ONT_cDNA` | `false` | **yes** | `-ax splice -uf` | `TRUE` |

If restranding does not reach `restrand_min_frac` on **every** sample, the run stops. It does not continue unstranded.

> [!IMPORTANT]
> A single execution must use one library type. Running PacBio and ONT samples together in the same run is not supported, as the difference in error profiles and library chemistry makes joint assembly difficult to interpret.

**PacBio.** The pipeline requires PacBio reads that have already been processed and stranded by PacBio's own independent standard workflows. pulposeq does not perform that preprocessing, and assumes it has already been done.

**ONT_DRS.** Direct RNA sequencing reads the native RNA molecule, so these libraries are considered automatically stranded. Setting `stranded_library: false` for `ONT_DRS` is ignored, and the pipeline emits a warning.

**ONT_cDNA.** This covers both PCR-cDNA and direct-cDNA protocols. Standard ONT cDNA library preparation sequences either the first or second cDNA strand, so read direction is roughly 50/50 with respect to the RNA strand. Set `stranded_library: false` for such reads — the pipeline then orients them itself with [Restrander](https://github.com/mritchielab/restrander), as described below. Set `stranded_library: true` only if you have already oriented the reads yourself, for example with [Pychopper](https://github.com/epi2me-labs/pychopper).

### Choosing the alignment preset with `map_hq`

`map_hq` selects between minimap2's two spliced presets. `splice:hq` is `splice` with three penalties changed:

| Option | `splice` | `splice:hq` | Penalises |
|---|---|---|---|
| `-B` | 2 | **4** | a mismatch |
| `-O` | 2,32 | **6,24** | opening a gap |
| `-C` | 9 | **5** | a junction that is not `GT-AG` |

Gap cost is `min(O1 + k×E1, O2 + k×E2)`, and both presets use `-E1,0`. So a gap of length *k* costs `min(2+k, 32)` under `splice` and `min(6+k, 24)` under `splice:hq`. Below 22 bases `splice:hq` is dearer, above it cheaper — short indels become harder to invent, introns easier to place.

Taken together, `splice:hq` **trusts the read**: it reads a discrepancy as real biology rather than basecalling error. That suits accurate data and misfires on noisy data, where it can assemble junctions out of errors.

**Leave `map_hq` unset** and the preset follows the library — `splice:hq` for PacBio, which is minimap2's documented Iso-Seq preset, and `splice` for ONT cDNA. Set it `true` or `false` to override for either platform.

Two reasons to set it explicitly:

- **`map_hq: true` for ONT cDNA.** R10.4.1 with `sup` basecalling is far more accurate than the chemistry `splice` was tuned for, and published work on this chemistry uses `splice:hq`. Worth evaluating rather than adopting untested: it changes junction placement, so it moves the class-code distribution and the mono-exonic count.
- **`map_hq: false` for PacBio.** Aligns both platforms identically, so the preset is not a confound when comparing their novel transcript calls. The pipeline warns, because `splice:hq` is the recommended preset for Iso-Seq, but the comparison is a legitimate reason to override it.

`map_hq` cannot be used with `ONT_DRS`: that preset sets `-k14` to cope with a high error rate, and `splice:hq` asserts the opposite, so the two cannot both be right about the same reads. Setting it is an error rather than a silently ignored value.

> [!NOTE]
> Changing `map_hq` changes the alignment command, so Nextflow re-runs alignment and everything downstream. Run it as a deliberate comparison, and not in the same run as a change to restranding — both move the mono-exonic count, in opposite directions, and you will not be able to attribute which did what.

### Restranding ONT cDNA libraries

When `library` is `ONT_cDNA` and `stranded_library` is `false`, the pipeline runs [Restrander](https://github.com/mritchielab/restrander) before filtering. It orients reads using poly(A)/poly(T) tails and the protocol's primer sequences, then reverse-complements the reverse-strand reads so that every read matches the orientation of the transcript it came from.

This matters more than it sounds. Without it, `minimap2` can only infer orientation from the splice signal, which it reports in the `ts:A` tag. Its own documentation is explicit about the limit:

> This tag is inferred from the GT-AG signal and is thus only available to spliced reads.
> — [minimap2 cookbook](https://github.com/lh3/minimap2/blob/master/cookbook.md)

So every unspliced read reaches Bambu with no transcript strand, and the curation step drops it. Spliced reads whose junctions are all non-canonical should be affected the same way, since there is no GT-AG signal to infer from — though minimap2 does not document that case explicitly, so treat it as expected rather than established.

The result is that mono-exonic novel transcripts disappear almost entirely, and correctly-stranded isoforms get misclassified as antisense. Poly(A) tails and primers are present on unspliced reads too, which is exactly why restranding recovers what splice-signal inference cannot.

> [!NOTE]
> The same documentation warns against the shortcut of simply declaring raw reads stranded: *"some intermediate reads are not stranded. For these reads, option `-uf` will lead to more errors."* Under `-ax splice` minimap2 defaults to `-ub`, searching both strands, which is correct for unoriented reads. `-uf` is only safe once the reads genuinely are oriented — which is what Restrander establishes.

Restranding is not optional for unoriented ONT cDNA. If the reads were already oriented outside the pipeline — with Pychopper, or with Restrander run before demultiplexing — declare that with `stranded_library: true` and the pipeline will trust it and skip this step.

No skip reaches it. Restrander is its own subworkflow, ahead of both QC and filtering, so `--skip_chopper` and `--skip_nanocomp` leave it untouched: the first changes which reads are filtered, the second only which reports are drawn. `--skip_nanocomp restrander` drops the NanoComp report on the oriented reads and nothing else.

#### Choosing a kit

`restrand_kit` is **required** and has no default, because the presets differ in their primer sequences and the wrong one produces a low orientation rate rather than an error.

| `restrand_kit` | Library preparation |
|---|---|
| `PCB109`, `PCB111`, `PCB114` | PCR-cDNA barcoding kits |
| `DCS109`, `DCS-LSK114` | Direct cDNA sequencing (no PCR) |
| `NEBNext` | NEBNext low-input / single-cell cDNA |
| `trimmed` | Any of the above, with the primers already removed |

Restrander also ships 10X Genomics presets. They are deliberately not accepted here — pulposeq is a bulk transcriptome pipeline and does not support single-cell data.

#### Barcoded, multiplexed libraries

**Demultiplex before running the pipeline.** The samplesheet takes one FASTQ per sample, so by the time pulposeq sees the data, one barcode is one row. Restrander then runs once per sample, independently, and never looks for barcodes at all.

Two things follow, and the second one catches people out.

**The kit is a property of the library preparation, not of demultiplexing.** A barcoded PCR-cDNA library is still a PCB library after demultiplexing. Selecting a `DCS` preset for it makes Restrander search for direct-cDNA primers that were never in those reads, orientation collapses, and the run silently falls back to unstranded.

**What demultiplexing does affect is whether the primers survived.** Demuxers trim to varying depths: some remove only the barcode, others take the surrounding adapter and primer with it. Restrander searches the first and last 200 bases of each read, so a barcode left in place is harmless — but a trimmed-away primer is invisible to every preset except `trimmed`.

Check rather than assume. Read the primer out of the preset that matches your chemistry:

```bash
docker run --rm emiyoshi/restrander:v1.1.3 cat /home/restrander/config/PCB114.json
```

Then count it, and its reverse complement, across the first 10,000 reads:

```bash
zcat sample.fastq.gz | head -40000 | awk 'NR%4==2' | grep -c "<tso sequence>"
```

Thousands of hits means the primers survived — use the chemistry preset. Near zero means use `trimmed`.

`trimmed` has two costs worth planning for. It classifies on poly(A)/poly(T) alone, so the orientation rate is lower than a primer-based preset can reach and may fall under `restrand_min_frac`. And artefact detection lives inside the primer stage, so `artefactStats` comes back empty — that is the configuration, not a failure.

#### Trimming parameters

Set `headcrop` and `tailcrop` to `0` when restranding. Restrander locates primers and poly(A)/poly(T) tails at exactly the read positions those crop. It runs before `chopper`, so cropping does not break it — but primer removal is Restrander's job, and the pipeline warns if either value is non-zero while restranding is active.

#### The success gate

Restranding is not assumed to have worked. Restrander reports how many reads it oriented, and the pipeline checks that fraction against `restrand_min_frac` (default `0.80`) **for every sample**. Any sample below the threshold stops the run, with an error listing every sample's fraction.

The check runs immediately after restranding rather than before Bambu, so a failure costs minutes rather than the hours Chopper and minimap2 would have spent first.

For reference, Restrander's own PCR-cDNA example data orients at **98.95%** with a matched preset. A run landing anywhere near `0.80` is telling you the kit or the trimming state is wrong, not that the threshold needs relaxing — which is why `restrand_min_frac` cannot be set below `0.5`, and warns below `0.75`.

When it fires, the fractions in the error message usually make the cause obvious:

- **all samples clustered just under the threshold** — the configuration doesn't match the reads. Check `restrand_kit` against your chemistry, and whether demultiplexing trimmed the primers.
- **one sample far below the rest** — that library has a problem. Excluding it is almost always better than lowering the threshold to accommodate it.

#### Why there is no unstranded fallback

Earlier versions continued with `-ax splice` and Bambu `stranded = FALSE` when strand could not be established. That behaviour has been removed, because the damage it causes is invisible in the outputs:

- **Mono-exonic novel transcripts disappear.** With no splice junctions there is nothing for minimap2 to infer orientation from, so these models get no strand and the curation step drops them. In matched runs over the same cell lines, an unstranded ONT run yielded **2** mono-exonic novel transcripts against PacBio's **462**.
- **Novel isoforms of known genes are misclassified as antisense.** GffCompare's `x` class means exonic overlap *on the opposite strand*, so every strand error converts a `j` into an `x`. The same matched runs show `x` inflated (13.7% against 9.8%) and `j` depleted by a corresponding amount (61.1% against 64.6%), while the strand-insensitive classes agree — `u` at 10.0% against 9.2%.
- **Quantification bleeds between overlapping transcripts.** Without strand, a read at a locus with sense and antisense transcription is compatible with both, so counts leak between them. This affects **annotated** transcripts too, not only novel ones, and it is worst exactly where a lncRNA pipeline cares most: every `*-AS1` lncRNA in the annotation is defined by overlapping a protein-coding gene on the opposite strand.

None of this raises a warning in the results. A transcript that was never called leaves no trace, and a misassigned read looks like any other read. Stopping the run is the only failure mode that is actually visible.

#### Reads that could not be oriented

Reads whose orientation cannot be determined are written to a separate `*-unknowns.fastq.gz` in `results/restrander/` instead of being passed on with an arbitrary strand. This is intentional: everything downstream treats the reads as genuinely oriented, so a read of unknown orientation surviving into that path would be assigned a coin-flip strand that the rest of the pipeline is instructed to trust.

Artefactual reads — those with the wrong combination of ends, reported as `TSO-TSO` or `RTP-RTP` — are counted in the report but **kept** in the main output. They are flagged, not discarded, so read counts are not silently perturbed and Bambu's automatically selected NDR is not shifted by the restranding step.

#### Custom configurations

For a protocol none of the presets cover, pass your own configuration with `restrand_config`, which overrides `restrand_kit`:

```bash
--restrand_config /path/to/my_protocol.json
```

`assets/restrander/` contains a working template and a field reference. Whichever configuration a run resolves to is copied into `results/restrander/` as `<sample>.restrander_config.json`, so the primer sequences behind every orientation call stay with the results.

### Ensembl vs GENCODE annotations

Transcript metadata is read **directly from the annotation `GTF` you supply** with `annotation`. Nothing is queried over the network, so runs are reproducible, work on offline compute nodes, and are unaffected by Ensembl archiving older releases or retiring service endpoints. It also guarantees the metadata describes exactly the annotation used for assembly, rather than whichever release a remote server happened to serve.

Both Ensembl and GENCODE annotations are accepted, and the differences between them are handled automatically:

| | Ensembl | GENCODE |
|---------------------|-------------------------------------------|--------------------------------------|
| Biotype attributes | `gene_biotype` / `transcript_biotype` | `gene_type` / `transcript_type` |
| Identifiers | unversioned, with separate `gene_version` / `transcript_version` | versioned inline (e.g. `ENSG00000290825.2`) |
| Sequence names | `1` | `chr1` |

Both identifier forms are always produced, so the `ensembl_transcript_id` and `ensembl_transcript_id_version` columns are populated whichever source you use.

Sequence names are **not** rewritten. `chromosome_name` reads `1` with an Ensembl annotation and `chr1` with GENCODE, exactly as the file you supplied names them. Stripping the prefix would have to special-case the GRC accessions (`KI270728.1`, `GL000009.2`) that both sources use for unplaced scaffolds, and it would put the metadata tables out of step with the novel-transcript outputs, which report the sequence names that reach them from the assembly untouched.

> [!IMPORTANT]
> With GENCODE, use the **CHR** or **PRI** annotation build. The **ALL** build additionally contains alternate loci, assembly patches and haplotypes, which duplicate gene and transcript identifiers. pulposeq checks for duplicated identifiers and stops with an error rather than silently double-counting them.

Regarding the pseudoautosomal regions (PAR) of chromosome Y: the gene annotation in these regions is identical between chromosomes X and Y, and until GENCODE release 43 the chrY copies were distinguished by an `_PAR_Y` suffix appended to their identifiers. From GENCODE release 44 (equivalent to Ensembl release 110) onward, the chrY PAR annotation carries its own distinct identifiers, so within the supported release range below this is not a concern. Additionally, the GENCODE GTF contains a number of attributes not present in the Ensembl GTF, but that does not interfere with the pipeline execution nor is required.

Check detailed information at the [GENCODE FAQ][1].

  [1]: https://www.gencodegenes.org/pages/faq.html

> [!WARNING]
> We highly recommend always using the latest stable release available, whether its GENCODE or Ensembl. The pipeline was only tested from Ensembl release 112 onwards. Due to reported differences between older versions, we recommend using, at least, from GENCODE v44 (Ensembl release 110) onwards. Versions older than that might cause issues with the transcriptome characterization steps.

## Running the pipeline {#running-the-pipeline}

The typical command for running the pipeline is as follows:

``` bash
nextflow run main.nf --input ./samplesheet.csv --outdir ./results --minqual [value] --reference [fasta] --annotation [gtf] --coding_potential_pred [cpc2/rnamining] --library [ONT_cDNA/ONT_DRS/PacBio] --stranded_library [true/false] -profile [test/medium/large],[container runtime: docker/singularity/apptainer]
```

Note that the pipeline will create the following files in your working directory:

``` bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications], other infrastructural tweaks (such as output directories), or module arguments (args).

  [tuning process resource specifications]: https://nf-co.re/docs/usage/configuration#tuning-workflow-resources

The above pipeline run specified with a params file in yaml format:

``` bash
nextflow run main.nf -profile docker -params-file params.yaml
```

with `params.yaml` containing:

``` yaml
input: './samplesheet.csv'
outdir: './results/'
<...>
```

> [!TIP]
> Follow the examples from the [test] and the [example run].

  [test]: ../test_data/testing.yml
  [example run]: ../examplerun.yml

### Updating the pipeline

``` bash
git clone https://github.com/integrativebioinformatics/pulposeq.git
```

When you run the above command, Git automatically clones the pipeline code from GitHub and stores it. When running the pipeline after this, it will always use this version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the commits in the pipeline:

``` bash
git fetch origin main
```

``` bash
git pull origin main
```

### Reproducibility

IIt is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [integrativebioinformatics/pulposeq releases page] and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

  [integrativebioinformatics/pulposeq releases page]: https://github.com/integrativebioinformatics/pulposeq/releases

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter files] to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

  [parameter files]: #running-the-pipeline

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a *single* hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Apptainer, Conda) - see below.

- `test`
  - Calibrated to the bundled chromosome 1 test data. Peaks at roughly 24 GB, so it runs on a workstation.
- `medium`
  - Full reference genome and annotation with modest sample sizes.
- `large`
  - Full reference with large ONT or PacBio samples, tens of GB each. Bambu alone can need several hundred GB; see [Resource requests](#resource-requests).
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)

Please only use [Conda](https://conda.io/docs/) as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, or Apptainer. Note that pulposeq does not support conda for local modules.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.


### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Resource allocations live in `conf/`. One file is always loaded and the profile you choose overrides it:

```
conf/base.config          always loaded, defines every label
    ↓ overridden by
conf/{test,medium,large}.config   only the labels each one defines
```

A label that a profile does not define keeps its `base.config` value. This is worth knowing, because
a process whose label no config matches falls back to the generic defaults **silently** — the usual
symptom is an out-of-memory kill (exit `137`) on a step you thought you had sized. To make unmatched
selectors report themselves, add the `debug` profile:

``` bash
nextflow run main.nf -profile large,apptainer,debug -params-file params.yml
```

#### Retries

`errorStrategy` resubmits on out-of-memory and related exit codes, with `maxRetries = 1` — so at most
one retry. Whether that retry asks for **more** than the first attempt depends on how the profile is
written:

``` groovy
memory = { 16.GB * task.attempt }   // escalates: 16 GB, then 32 GB
memory = 16.GB                      // does not: 16 GB, then 16 GB again
```

Fixed values are appropriate where the input is known and unchanging, such as the `test` profile.
For real data, `task.attempt` scaling lets you allocate close to the measured requirement and still
survive an unusual sample.

#### `resourceLimits`

`resourceLimits` caps every request **after** the `task.attempt` multiplier, and it clamps silently.
Two consequences:

- It must be at or above the largest request in the file, or that label is quietly cut back.
- If you rely on retry escalation, it must be at or above **twice** the largest first attempt, or the
  retry is clamped to the value that already failed.

#### Tuning for your own system

Prefer a site config passed with `-c` over editing the repository, so your settings survive
`git pull`:

``` groovy
// ~/.nextflow/config or site.config, used with -c
process {
    executor       = 'slurm'
    queue          = 'main'
    resourceLimits = [ cpus: 256, memory: 700.GB, time: 24.h ]

    withLabel:process_bambu { memory = 700.GB }   // override a single label
}

apptainer.cacheDir = '/path/to/apptainer_images'
```

### Measured resource requirements

Peak resident memory observed in real runs. Use these to judge whether a dataset fits your hardware.

| Process | test (chr1 subset) | medium (full reference, small samples) | large (8 ONT samples, 25–30 GB each) |
|---|---|---|---|
| `CHOPPER` | 7.5 GB | ~8 GB | **54.7 GB** |
| `MINIMAP2_ALIGN` | 16.3 GB (8 cpu) | 31.1 GB (20 cpu) | 56.1 GB (32 cpu) |
| `NANOCOMP` | 1.8 GB | 1.2 GB | 33.4 GB |
| `NANOCOMP_MINIMAP2` | 1.5 GB | 3.5 GB | 48.4 GB |
| `BAMBU` | 13.6 GB | 63.6 GB | **639.8 GB** |
| metadata refinement | 0.9 GB | ~4 GB | 4.1 GB |
| `RENDER_REPORT` | 1.0 GB | 1.3 GB | — |

Four rules let you predict your own case rather than guessing:

**Chopper is linear in input size**, at close to **0.38 × bytes read**, verified from 2.5 GB to
143 GB of input. A 30 GB gzipped FASTQ expands to roughly 130 GB, so expect a ~50 GB peak.

**Bambu grows by roughly 54 GB per additional sample** — 532 GB at 6 samples, 640 GB at 8 — because
it holds every sample in a single object.

> [!IMPORTANT]
> This puts a ceiling on samples per run. On a 768 GB node, Bambu becomes unschedulable somewhere
> around **10 samples**. It is a property of the assembly step, not something configuration can
> tune away: beyond that you need a larger-memory node, or to split the run.

**minimap2 and NanoComp scale with thread count**, since their memory is a fixed index plus
per-thread buffers. Raising `cpus` raises memory too — NanoComp went from 23.8 GB on 1 CPU to
33.4 GB on 3.

**The metadata refinement steps scale with the annotation, not the data.** They parse the reference
`GTF`, so a full human annotation costs ~4 GB whether the input is a single small sample or eight
large ones. A chromosome-1 test run understates them by roughly fourfold.

### SLURM and `--mem-per-cpu`

Some clusters schedule by memory per core rather than by total. Nextflow supports this:

``` groovy
executor.perCpuMemAllocation = true   // Nextflow 23.10+, SLURM only
```

With it enabled, jobs are submitted as `--mem-per-cpu <memory / cpus>` instead of `--mem <memory>`.
**The `memory` directive still means total memory** — Nextflow does the division at submission — so
there is no second set of numbers to maintain, and the setting is ignored by the local executor.

What changes is that `cpus` becomes load-bearing. Every memory value you raise also raises the
per-core figure, which must stay under your partition's `MaxMemPerCPU`:

``` bash
scontrol show config | grep -iE "MaxMemPerCPU|DefMemPerCPU"
sinfo -o "%P %c %m" | sort -u
```

Where a step needs more memory than that ratio allows, **raise `cpus` rather than lowering memory**.
Requesting cores a single-threaded tool will not use is the price of getting the memory at all, and
the cores are unusable by other jobs anyway once the node's memory is committed.

This belongs in your site config, not in the pipeline — it describes a SLURM installation rather
than pulposeq.

### Reading the resource reports

Every run writes `pipeline_info/execution_trace_*.txt`, which is more useful than the HTML report for
this purpose: it carries `rchar` and `wchar` alongside memory, so you can relate usage to input size.

A few traps worth knowing:

- **Size from `peak_rss`, never `peak_vmem`.** `RENDER_REPORT` reports over 1 TB of virtual memory
  against about 1 GB resident — address space that R and Quarto reserve without ever touching.
- **`peak_rss` is unreliable for piped commands.** `CHOPPER` runs `zcat | chopper | gzip` and has
  reported `0.00 GB` while genuinely using ~50 GB, because short-lived children escape Nextflow's
  sampler. For those, ask SLURM instead:

  ``` bash
  sacct -j <jobid> --format=JobID,JobName,MaxRSS,ReqMem,State,ExitCode
  ```

- **Nextflow reads high relative to SLURM**, by around 12% on the same task. Its figures are already
  conservative, so there is no need to add a large safety factor on top.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#update-tool-versions) section of the nf-core website.


### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#modifying-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time. Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory. We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

``` bash
NXF_OPTS='-Xms1g -Xmx16g'
```
