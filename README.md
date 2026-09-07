## integrativebioinformatics/pulposeq <img src="figures/logo.svg" align=right height="200px"/>
[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/integrativebioinformatics/pulposeq) [![Nextflow](https://img.shields.io/badge/version-%E2%89%A526.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/) [![nf-core template version](https://img.shields.io/badge/nf--core_template-4.0.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.0.2) [![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/) [![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

`pulposeq` is a Nextflow pipeline designed for isoform-level lncRNA discovery and characterization from long-read RNA-seq data. The workflow encompasses QC, mapping, transcriptome assembly and quantification, followed by a detailed final characterization of the entire transcriptome with particular emphasis on lncRNA structure and isoforms across known annotations and novel candidates.

For more details and further functionality, please refer to the [usage](docs/usage.md) and [output](docs/output.md) documentations.

**The workflow**

![pulposeq workflow](figures/pulposeq.drawio.svg)

We can describe each step of the workflow as follows:

1.  Quality control of reads ([NanoComp](https://github.com/wdecoster/nanocomp "wdecoster/nanocomp"))
2.  Recover read orientation, for ONT cDNA that is not already oriented ([Restrander](https://github.com/mritchielab/restrander "mritchielab/restrander"))
3.  Filtering and trimming ([chopper](https://github.com/wdecoster/chopper "wdecoster/chopper"))
4.  Mapping to a genome reference ([minimap2](https://github.com/lh3/minimap2 "lh3/minimap2") and [samtools](https://github.com/samtools/samtools "samtools"))
5.  Quality control of mapped reads ([NanoComp](https://github.com/wdecoster/nanocomp "wdecoster/nanocomp"))
6.  Per-base coverage tracks as `bigWig`, one per sample ([rtracklayer](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html) and [GenomicAlignments](https://bioconductor.org/packages/release/bioc/html/GenomicAlignments.html))
7.  Transcriptome Assembly ([Bambu](https://www.bioconductor.org/packages/release/bioc/html/bambu.html))
8.  Compare novel transcripts to the annotation reference ([GffCompare](https://github.com/gpertea/gffcompare "gpertea/gffcompare"))
9.  Convert novel transcripts `GTF` file to `FASTA` ([GffRead](https://github.com/gpertea/gffread "gpertea/gffread"))
10. Predict transcripts as protein-coding or non-coding ([CPC2](https://github.com/gao-lab/CPC2_standalone "gao-lab/CPC2_standalone") by default, or [RNAmining](https://gitlab.com/integrativebioinformatics/RNAmining "integrativebioinformatics/RNAmining"))
11. Gather all data from previous steps and generate informative and re-usable metadata `.csv` and `GTF` files for both novel and annotated transcripts, with biotypes read directly from the supplied reference annotation ([tidyverse](https://tidyverse.org/), [rtracklayer](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html), and [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html))
12. Restrict the count matrices and `GTF` files to the curated set, and attach biotype and classification attributes to the validated annotations
13. Regenerate Bambu's PCA and expression heatmaps from the curated transcriptome, so the sample-level view can be read beside the one built from the raw assembly
14. Draw genomic context figures over the coverage tracks — known genes carrying novel isoforms, and the sense-intronic candidates, which lie inside a host intron on the host's own strand and cannot be resolved from the assembly alone ([plotgardener](https://bioconductor.org/packages/release/bioc/html/plotgardener.html))
15. Provide a report and data visualization for the full transcriptome, with emphasis on lncRNAs ([Quarto](https://quarto.org/), [tidyverse](https://tidyverse.org/), [cowplot](https://cran.r-project.org/web/packages/cowplot/index.html), [scales](https://cran.r-project.org/web/packages/scales/index.html), etc)
16. Gather all possible QC information from the previous steps ([MultiQC](https://github.com/MultiQC/MultiQC "MultiQC"))

**What you get**

| Path | |
|---|---|
| `transcriptome_report/report.html` | The self-contained report, every figure embedded |
| `bambu_validated/` | Curated count matrices, annotations and the curation yield |
| `novel_transcripts/`, `ref_transcripts/` | Metadata `.csv` and `GTF` per category |
| `genomic_context/` | Coverage and transcript models drawn at selected loci |
| `coverage/*.bw` | Per-sample coverage as bigWig — two to three orders of magnitude smaller than the alignments, so this is the form of the data meant to leave the cluster. Load them in IGV to explore any region beyond the windows the pipeline chose to draw |

## Usage

pulposeq is compatible with **Ensembl or GENCODE** reference genomes and annotations. Transcript and gene biotypes are read directly from the annotation you supply, so no Ensembl release or BioMart dataset needs to be declared.

Coding potential is predicted with [CPC2](https://github.com/gao-lab/CPC2_standalone "gao-lab/CPC2_standalone") by default, which is organism-agnostic — so **the pipeline is no longer restricted to a fixed species list**. Any eukaryote with an Ensembl or GENCODE annotation is in scope; splice-aware alignment and the intron-based gffcompare class codes are what make the eukaryotic assumption, not the predictor.

[RNAmining](https://gitlab.com/integrativebioinformatics/RNAmining "integrativebioinformatics/RNAmining") remains selectable with `coding_potential_pred: "rnamining"`, and then requires an `organism` from its supported list:

> *Homo sapiens, Mus musculus, Danio rerio, Anolis carolinensis, Chrysemys picta bellii, Crocodylus porosus, Eptatretus burgeri, Gallus gallus, Latimeria chalumnae, Monodelphis domestica, Notechis scutatus, Ornithorhynchus anatinus, Petromyzon marinus, Rattus norvegicus, Sphenodon punctatus,* and *Xenopus tropicalis.*

> [!WARNING]
> RNAmining is currently under review and is **not recommended**. On this pipeline's test data it classified 86 of 100 GENCODE protein-coding transcripts as non-coding, where CPC2 misclassified 1. It is kept selectable so earlier runs can be reproduced and the cause investigated. See [Coding potential](docs/usage.md#coding-potential).

> [!WARNING]
> pulposeq requires Nextflow `>=26.04.0`, where the strict syntax parser is enabled by default. Make sure to setup appropriate configuration. See the current documentation at [Seqera Docs](https://docs.seqera.io/nextflow/strict-syntax).

### Nextflow setup and testing the pipeline

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data. The pipeline is compatible with Docker and Singularity/Apptainer.

You can run an example test by following the instructions:

Enter the `test_data` folder

``` bash
cd test_data
```

Download and unzip the reference `FASTA` and `GTF` files, and also download the fastq.gz files:

Make the file executable

``` bash
chmod +x download-ref-fastq.sh
```

Execute the script

``` bash
./download-ref-fastq.sh
```

Add YOUR full path for the samples in the `samplesheet.csv` ([file](test_data/samplesheet.csv)). For example, your full path for a sample could be:

`home/user/pulposeq/test_data/thesample.fastq.gz`

Go back to the main directory and execute the test with appropriate container runtime!

``` bash
cd ..
```

``` bash
nextflow run main.nf -profile test,[container runtime] -params-file test_data/testing.yml
```

Set the container runtime profile to `docker`, `singularity` or `apptainer`, according to your resources. pulposeq does **NOT** support nor recommend execution with `conda` environments for the local modules.

The `test` profile is calibrated to the bundled chromosome 1 data and peaks at roughly 24 GB, so it runs on a workstation. For real datasets use `medium` or `large`, and see [Resource requests](docs/usage.md#resource-requests) for measured requirements at each scale -- Bambu in particular can need several hundred GB.

In some cases, depending on your system's permissions and configuration regarding containers, you might need to set specific configurations before running. In that case, follow the instructions available at [Seqera Docs](https://docs.seqera.io/nextflow/reference/config).

> *As long as you use the supported container runtimes and do not modify the pipeline architecture, pulposeq ensures ***reproducibility*** independently of your system's specific setup.*

> [!WARNING]
> Please provide pipeline parameters via the CLI or input a `yaml` or `json`parameters file  via the Nextflow `-params-file` option (most recommended). Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration, ***except for parameters***; see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

You must declare your sequencing library chemistry with the `library` (`ONT_cDNA`, `ONT_DRS` or `PacBio`) and `stranded_library` (`true`/`false`) parameters. These set the `minimap2` alignment preset and Bambu's strandedness together, and a single execution must use one library type. PacBio libraries are expected to have been processed and stranded beforehand by PacBio's own standard workflows, and ONT direct RNA is oriented by construction.

Unoriented ONT cDNA (`library: "ONT_cDNA"` with `stranded_library: false`) is re-oriented **inside** the pipeline by Restrander, which then requires a `restrand_kit` matching your chemistry — there is no default, since the wrong preset yields a low orientation rate rather than an error. See [Library type and strandedness](docs/usage.md#library-type-and-strandedness) and [Restranding ONT cDNA libraries](docs/usage.md#restranding-ont-cdna-libraries) for the full details.


pulposeq publishes results in `symlink` mode by default. Every file in your output directory is a symbolic link pointing at the real file inside `work/`, rather than a copy of it. This costs no additional disk space and completes instantly, which matters when the outputs are large `BAM` and `GTF` files.

The `work/` directory holds your actual results alongside every intermediate execution file (temporary files, `.command.sh`, `.command.log`). It is also what makes `-resume` able to restart from the last successful step if a run fails.

> [!WARNING]
> With `symlink` publishing, **deleting `work/` destroys your results** — the links in the output directory are left pointing at files that no longer exist. Before removing `work/`, or before archiving or sharing an output directory, either dereference the links or re-run with `--publish_dir_mode copy`.
>
> To copy rather than link, either pass `--publish_dir_mode copy` on the command line, set `publish_dir_mode: copy` in your params file, or change the default in [`nextflow.config`](nextflow.config). Bear in mind that copying large datasets takes noticeably longer and doubles the disk they occupy. The available modes are documented in the [Nextflow documentation](https://docs.seqera.io/nextflow/reference/process#publishdir).
>
> To turn an existing symlinked output directory into real files, dereference it into a new location:
>
> ``` bash
> cp -rL results results_standalone
> ```


## Citations

If you use pulposeq in your research, please consider citing it.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) initative, and reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> The nf-core framework for community-curated bioinformatics pipelines.
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> Nat Biotechnol. 2020 Feb 13. doi: 10.1038/s41587-020-0439-x.

An extensive list of references for the tools incorporated by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

## Acknowledgments

#### Development & Contributions

The `pulposeq` pipeline was originally developed by **Bárbara Borges** ([\@borgessbarbara](https://github.com/borgessbarbara)). We extend our sincere thanks to **Lucas Freitas** ([\@lfreitasl](https://github.com/lfreitasl)), **Gleison Azevedo** ([\@gleisonm](https://github.com/gleisonm)) and **João Cavalcante** ([\@jvfe](https://github.com/jvfe)) for their significant contributions and assistance during development.

#### Supervision & Collaborations

This project was developed under the supervision of **Vinicius Maracaja-Coutinho** and **Thaís Gaudencio**, in collaboration with **Rodrigo Dalmolin** and the [Dalmolin Systems Biology Group](https://github.com/dalmolingroup)

#### Computational resources

This project was supported by the High-Performance Computing Center at UFRN (NPAD/UFRN) and the National Laboratory for High Performance Computing (NLHPC) (CCSS210001) at UChile.

#### Funding

CAPES (001), CNPq (MCTI/FNDCT 445067/2024-1), FONDECYT-ANID Postdoctorado (3250452), FONDECYT-ANID (1211731), FONDAP-ANID (15130011 and 1523A0008), and Anillo-ANID (ATE220016).

#### Laboratories & Institutions Involved

Bárbara Borges <sup>1,2,3</sup>, Lucas Freitas <sup>4,5</sup>, Gleison Azevedo <sup>1</sup>, João Cavalcante <sup>1</sup>, Rodrigo Dalmolin <sup>1</sup>, Thaís Gaudencio <sup>1,3</sup>, Vinicius Maracaja-Coutinho <sup>1,2,6,7</sup>

<details>

<summary>Affiliations</summary>

<sup>1</sup> Bioinformatics Multidisciplinary Environment, Instituto Metrópole Digital, Universidade Federal do Rio Grande do Norte, Natal, Brazil

<sup>2</sup> Unidad de Genómica Avanzada, Advanced Center for Chronic Diseases, Facultad de Ciencias Químicas y Farmacéuticas, Universidad de Chile, Santiago, Chile

<sup>3</sup> Artificial Intelligence Applications, Universidade Federal da Paraíba, João Pessoa, Brazil

<sup>4</sup> Faculty of Natural Sciences, Universität Hohenheim, Stuttgart, Germany

<sup>5</sup> Staatliches Museum für Naturkunde Stuttgart, Department of Biodiversity Monitoring, Stuttgart, Germany

<sup>6</sup> Laboratório de Medicina e Saúde Pública de Precisão, Instituto Gonçalo Moniz, Fundação Oswaldo Cruz, Fiocruz, Salvador, Brazil

<sup>7</sup> Instituto Nacional de Ciência e Tecnologia em Saúde Digital (INCT-DigiSaúde)

</details>

<p align="center">

<picture> <source media="(prefers-color-scheme: dark)" srcset="figures/institutional-logos-dark-theme.png" height=150> <img src="figures/institutional-logos-light.png" height="150"/> </picture>

</p>
