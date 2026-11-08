/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MULTIQC                           } from '../modules/nf-core/multiqc/main'
include { NOVEL_TRANSCRIPTS                 } from '../modules/local/metadata_refinement/novel_transcripts/main'
include { PREFILTER                         } from '../modules/local/metadata_refinement/prefilter/main'
include { CURATE                            } from '../modules/local/metadata_refinement/curate/main'
include { REF_TRANSCRIPTS                   } from '../modules/local/metadata_refinement/ref_transcripts/main'
include { ANNOTATION                        } from '../modules/local/metadata_refinement/annotation/main'
include { BAM_COVERAGE                      } from '../modules/local/bam_coverage/main'
include { GENOMIC_CONTEXT                   } from '../modules/local/genomic_context/main'
include { VALIDATE_BAMBU_OUTPUTS            } from '../modules/local/validate_bambu_out/main'
include { RENDER_REPORT                     } from '../modules/local/report/main'
include { RESTRANDING                       } from '../subworkflows/local/restranding'
include { QC_FILT                           } from '../subworkflows/local/qc'
include { ALIGNMENT                         } from '../subworkflows/local/alignment'
include { ASSEMBLY                          } from '../subworkflows/local/assembly'
include { CLASSIFICATION                    } from '../subworkflows/local/classification'
include { QUANTIFICATION                    } from '../subworkflows/local/quantification'
include { oarfishEnabled                    } from '../subworkflows/local/utils_nfcore_pulposeq_pipeline'
include { paramsSummaryMap                  } from 'plugin/nf-schema'
include { paramsSummaryMultiqc              } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML            } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText            } from '../subworkflows/local/utils_nfcore_pulposeq_pipeline'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PULPOSEQ {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions             = channel.empty()
    ch_multiqc_files        = channel.empty()
    ch_gtf_new_transcripts  = channel.empty()
    // The oarfish channels stay empty lists unless oarfish runs
    ch_oarfish_counts       = channel.value([])
    ch_oarfish_gene         = channel.value([])
    ch_oarfish_cpm          = channel.value([])
    ch_oarfish_tpm          = channel.value([])
    ch_oarfish_bootstraps   = channel.value([])

    // The three structural skips are nested, because each stage consumes the last
    // one's output: nothing downstream of minimap2 can run without alignments, and
    // nothing downstream of Bambu can run without its transcriptome. Skipping a
    // stage therefore skips everything after it -- there is no way to ask for the
    // classification of a transcriptome that was never assembled.
    def run_alignment       = !params.skip_minimap2
    def run_bambu           = run_alignment && !params.skip_bambu
    def run_classification  = run_bambu && !params.skip_class

    // oarfish is opt-in. It quantifies against Bambu's models, so it can only run
    // with Bambu; validateInputParameters stops the run if --run_oarfish is combined
    // with a skip that removes Bambu. When it runs, its counts -- not Bambu's --
    // decide which transcripts PREFILTER keeps.
    def run_oarfish         = run_bambu && oarfishEnabled()
    def expression_source   = run_oarfish ? 'oarfish' : 'bambu'

    // What the user asserted about the reads on disk. ONT_DRS is native RNA and
    // PacBio is oriented upstream by lima, so both are stranded by construction;
    // ONT_cDNA is whatever the user declares, and if they declare it unoriented the
    // pipeline restrands it rather than proceeding without strand.
    //
    // Note this describes the FASTQ, not the library chemistry. ONT PCR-cDNA
    // chemistry is strand-specific while its basecalled reads are not oriented,
    // which is why the two must not be conflated.
    def declared_stranded = params.library in ['ONT_DRS', 'PacBio'] ||
                            params.stranded_library?.toString()?.toLowerCase() == 'true'

    //
    // Orientation correction, before any QC or filtering
    //
    // Its own subworkflow, and outside every skip. Restranding is not quality
    // control -- it is what makes read strand equal RNA strand, which the alignment
    // preset and Bambu both assume. There is deliberately no flag that turns it off
    // for reads that need it.
    //
    RESTRANDING (
        ch_samplesheet,
        declared_stranded
    )
    ch_reads_for_alignment = RESTRANDING.out.reads
    ch_multiqc_files = ch_multiqc_files.mix(RESTRANDING.out.multiqc)
    ch_versions = ch_versions.mix(RESTRANDING.out.versions)

    //
    // Run QC workflow
    //
    // Takes the reads as supplied for the raw-read report, and the restranded ones
    // for filtering, so NanoComp still describes the input rather than the oriented
    // subset.
    //
    // Called unconditionally: reporting and filtering are gated separately inside,
    // by --skip_nanocomp and --skip_chopper, so that turning the reports off does
    // not also change which reads reach alignment.
    //
    QC_FILT (
        ch_samplesheet,
        RESTRANDING.out.reads
    )
    ch_reads_for_alignment = QC_FILT.out.filt_reads
    ch_multiqc_files = ch_multiqc_files.mix(QC_FILT.out.multiqc)
    ch_versions = ch_versions.mix(QC_FILT.out.versions)

    //
    // Run alignment workflow
    //
    if (run_alignment) {
        ALIGNMENT(ch_reads_for_alignment)

        // Empty when --skip_nanocomp drops the minimap2 report, so no gate here
        ch_multiqc_files = ch_multiqc_files.mix(ALIGNMENT.out.multiqc)
        ch_versions = ch_versions.mix(ALIGNMENT.out.versions)

        //
        // Per-sample coverage tracks. bigWig rather than BAM: plotgardener's
        // plotSignal cannot read BAM, and the tracks are two to three orders of
        // magnitude smaller than the alignments they come from, so they are also
        // the only form of this data that travels off the cluster.
        //
        BAM_COVERAGE (
            ALIGNMENT.out.bam.join(ALIGNMENT.out.index),
            file("${projectDir}/bin/bam_coverage.R", checkIfExists: true)
        )
        ch_versions = ch_versions.mix(BAM_COVERAGE.out.versions)
    }

    //
    // Transcriptome assembly with Bambu
    //
    if (run_bambu) {
        ASSEMBLY (
            ALIGNMENT.out.bam,
            params.reference,
            params.annotation
        )

        ch_versions = ch_versions.mix(ASSEMBLY.out.versions)

        //
        // Optional oarfish quantification, against every Bambu model and on the same
        // reads minimap2 aligned
        //
        if (run_oarfish) {
            QUANTIFICATION (
                ch_reads_for_alignment,
                ASSEMBLY.out.gtf_all_transcripts,
                params.reference,
                ASSEMBLY.out.transcript_counts
            )
            ch_oarfish_counts     = QUANTIFICATION.out.counts_transcript
            ch_oarfish_gene       = QUANTIFICATION.out.counts_gene
            ch_oarfish_cpm        = QUANTIFICATION.out.cpm_transcript
            // TPM exists for PacBio only, bootstraps only when requested: an empty
            // list stands in so curation is not left waiting for a missing file
            ch_oarfish_tpm        = QUANTIFICATION.out.tpm_transcript.ifEmpty([])
            ch_oarfish_bootstraps = QUANTIFICATION.out.bootstraps.ifEmpty([])
            ch_versions = ch_versions.mix(QUANTIFICATION.out.versions)
        }

        //
        // The one place transcripts are filtered. Runs whenever Bambu runs, before
        // classification, so gffread, gffcompare, CPC2 and everything curated after
        // them see the same transcript list. Unstranded and < 500 nt novel models are
        // removed here, and any transcript below --min_total_counts summed across all
        // samples -- counted by oarfish when it runs, by Bambu otherwise. Bambu's raw
        // output stays untouched.
        //
        PREFILTER (
            ASSEMBLY.out.transcript_counts,
            ASSEMBLY.out.gene_counts,
            ASSEMBLY.out.CPM,
            ASSEMBLY.out.full_length,
            ASSEMBLY.out.unique_counts,
            ASSEMBLY.out.gtf_new_transcripts,
            expression_source,
            ch_oarfish_counts,
            params.min_total_counts,
            file("${projectDir}/bin/prefilter_transcripts.R", checkIfExists: true)
        )
        ch_versions = ch_versions.mix(PREFILTER.out.versions)

        // Only the kept novel models go on to classification
        PREFILTER.out.novel_gtf
            .set { ch_gtf_new_transcripts }
    }

    if (run_classification) {
        CLASSIFICATION (
            ch_gtf_new_transcripts,
            params.annotation,
            params.reference
        )

        NOVEL_TRANSCRIPTS (
            ch_gtf_new_transcripts,
            CLASSIFICATION.out.annotated_gtf,
            CLASSIFICATION.out.tmap,
            CLASSIFICATION.out.predictions,
            PREFILTER.out.counts_transcript_filter,
            PREFILTER.out.counts_gene_filter,
            ASSEMBLY.out.bambu_tx_classes,
            params.annotation,
            file("${projectDir}/bin/novel_transcripts.R", checkIfExists: true),
            file("${projectDir}/bin/gtf_annotation_utils.R", checkIfExists: true)
        )
        ch_versions = ch_versions.mix(NOVEL_TRANSCRIPTS.out.versions)

        //
        // The final transcript set: reference transcripts PREFILTER kept, plus the
        // novel models that passed curation. Every table and GTF from here on --
        // Bambu's and, when it ran, oarfish's -- holds exactly these transcripts.
        //
        CURATE (
            NOVEL_TRANSCRIPTS.out.novel_combined_metadata,
            PREFILTER.out.counts_transcript_filter,
            PREFILTER.out.counts_gene_filter,
            PREFILTER.out.cpm_transcript_filter,
            PREFILTER.out.full_length_counts_transcript_filter,
            PREFILTER.out.unique_counts_transcript_filter,
            ASSEMBLY.out.gtf_all_transcripts,
            expression_source,
            ch_oarfish_counts,
            ch_oarfish_gene,
            ch_oarfish_cpm,
            ch_oarfish_tpm,
            ch_oarfish_bootstraps,
            file("${projectDir}/bin/curate_transcripts.R", checkIfExists: true)
        )
        ch_versions = ch_versions.mix(CURATE.out.versions)

        REF_TRANSCRIPTS (
            CURATE.out.counts_transcript_validated,
            CURATE.out.counts_gene_validated,
            ASSEMBLY.out.gtf_all_transcripts,
            params.annotation,
            file("${projectDir}/bin/ref_transcripts.R", checkIfExists: true),
            file("${projectDir}/bin/gtf_annotation_utils.R", checkIfExists: true)
        )
        ch_versions = ch_versions.mix(REF_TRANSCRIPTS.out.versions)

        //
        // The final annotation, for plotting. The validated GTFs are read here but
        // never rewritten -- they are published as CURATE wrote them.
        //
        ANNOTATION (
            CURATE.out.annotations_validated_gtf,
            REF_TRANSCRIPTS.out.transcriptome_metadata,
            NOVEL_TRANSCRIPTS.out.novel_combined_metadata,
            params.annotation,
            file("${projectDir}/bin/annotation.R", checkIfExists: true),
            file("${projectDir}/bin/gtf_annotation_utils.R", checkIfExists: true)
        )
        ch_versions = ch_versions.mix(ANNOTATION.out.versions)

        //
        // Genomic context figures: known genes carrying novel isoforms, plus the
        // sense-intronic candidates. A model inside a host intron on the host's own
        // strand is sequence-identical to part of that host's unspliced precursor,
        // so it is drawn on the whole intron and read by eye rather than scored.
        //
        // The final annotation, not the Bambu-derived one. Its known transcripts
        // come from the reference with CDS and UTR features intact, which is what
        // lets plotTranscripts draw a coding region thicker than its UTRs.
        GENOMIC_CONTEXT (
            ANNOTATION.out.final_gtf,
            BAM_COVERAGE.out.bigwig.map { _meta, bw -> bw }.collect(),
            CURATE.out.counts_transcript_validated,
            CURATE.out.full_length_counts_transcript_validated,
            params.annotation,
            file("${projectDir}/bin/genomic_context.R", checkIfExists: true)
        )
        ch_versions = ch_versions.mix(GENOMIC_CONTEXT.out.versions)

        //
        // Regenerate the Bambu PCA and heatmaps from the validated transcriptome
        //
        VALIDATE_BAMBU_OUTPUTS (
            ASSEMBLY.out.rds_transcript,
            ASSEMBLY.out.rds_gene,
            CURATE.out.counts_transcript_validated,
            CURATE.out.counts_gene_validated,
            file("${projectDir}/bin/validate_bambu_out.R", checkIfExists: true)
        )
        ch_versions = ch_versions.mix(VALIDATE_BAMBU_OUTPUTS.out.versions)

        RENDER_REPORT (
            CURATE.out.counts_gene_validated,
            CURATE.out.counts_transcript_validated,
            CURATE.out.full_length_counts_transcript_validated,
            CURATE.out.unique_counts_transcript_validated,
            REF_TRANSCRIPTS.out.transcriptome_metadata,
            REF_TRANSCRIPTS.out.protein_coding_metadata,
            REF_TRANSCRIPTS.out.lncrna_metadata,
            REF_TRANSCRIPTS.out.protein_coding_exonlength,
            REF_TRANSCRIPTS.out.lncrna_exonlength,
            NOVEL_TRANSCRIPTS.out.novel_combined_metadata,
            NOVEL_TRANSCRIPTS.out.novel_lncrna_exon_lengths,
            NOVEL_TRANSCRIPTS.out.novel_mrna_exon_lengths,
            NOVEL_TRANSCRIPTS.out.novel_non_coding_metadata,
            GENOMIC_CONTEXT.out.candidates,
            GENOMIC_CONTEXT.out.figures.ifEmpty([]),
            GENOMIC_CONTEXT.out.intronic_candidates,
            GENOMIC_CONTEXT.out.intronic_figures.ifEmpty([]),
            ASSEMBLY.out.pca,
            ASSEMBLY.out.pca_grouped,
            ASSEMBLY.out.h_gene,
            ASSEMBLY.out.h_transcript,
            VALIDATE_BAMBU_OUTPUTS.out.pca,
            VALIDATE_BAMBU_OUTPUTS.out.pca_grouped,
            VALIDATE_BAMBU_OUTPUTS.out.h_gene,
            VALIDATE_BAMBU_OUTPUTS.out.h_transcript,
            ASSEMBLY.out.bambu_metrics,
            VALIDATE_BAMBU_OUTPUTS.out.validation_summary,
            file("${projectDir}/bin/report.qmd", checkIfExists: true)
        )
        ch_versions = ch_versions.mix(RENDER_REPORT.out.versions)
    }
//
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'pulposeq_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    def multiqc_config_file                   = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    def multiqc_custom_config_file            = params.multiqc_config ? file(params.multiqc_config, checkIfExists: true) : []
    def multiqc_logo_file                     = params.multiqc_logo ? file(params.multiqc_logo, checkIfExists: true) : []
    def summary_params                        = paramsSummaryMap(workflow)
    def ch_workflow_summary                   = channel.value(paramsSummaryMultiqc(summary_params))
    def ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description                = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    def multiqc_configs = multiqc_custom_config_file ? [multiqc_config_file, multiqc_custom_config_file] : [multiqc_config_file]

    MULTIQC (
        ch_multiqc_files.collect().map { files ->
            [ [id: 'multiqc'], files, multiqc_configs, multiqc_logo_file, [], [] ]
        }
    )

    emit:
    multiqc_report = MULTIQC.out.report.map { _meta, report -> report }.toList()
    versions       = ch_versions

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
