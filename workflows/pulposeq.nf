/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MULTIQC                           } from '../modules/nf-core/multiqc/main'
include { NOVEL_TRANSCRIPTS                 } from '../modules/local/metadata_refinement/novel_transcripts/main'
include { FILTER_BAMBU_COUNTS               } from '../modules/local/metadata_refinement/filter_bambu_counts/main'
include { VALIDATE_BAMBU_GTF                } from '../modules/local/metadata_refinement/validate_bambu_gtf/main'
include { VALIDATE_BAMBU_COUNTS             } from '../modules/local/metadata_refinement/validate_bambu_counts/main'
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

    // The three structural skips are nested, because each stage consumes the last
    // one's output: nothing downstream of minimap2 can run without alignments, and
    // nothing downstream of Bambu can run without its transcriptome. Skipping a
    // stage therefore skips everything after it -- there is no way to ask for the
    // classification of a transcriptome that was never assembled.
    def run_alignment       = !params.skip_minimap2
    def run_bambu           = run_alignment && !params.skip_bambu
    def run_classification  = run_bambu && !params.skip_class

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

        ASSEMBLY.out.gtf_new_transcripts
            .set { ch_gtf_new_transcripts }

        ch_versions = ch_versions.mix(ASSEMBLY.out.versions)
    }

    if (run_classification) {
        CLASSIFICATION (
            ch_gtf_new_transcripts,
            params.annotation,
            params.reference
        )
        ch_versions = ch_versions.mix(CLASSIFICATION.out.versions)

        FILTER_BAMBU_COUNTS (
            ASSEMBLY.out.gene_counts,
            ASSEMBLY.out.transcript_counts,
            ASSEMBLY.out.CPM,
            ASSEMBLY.out.full_length,
            ASSEMBLY.out.unique_counts
        )

        NOVEL_TRANSCRIPTS (
            ch_gtf_new_transcripts,
            CLASSIFICATION.out.annotated_gtf,
            CLASSIFICATION.out.tmap,
            CLASSIFICATION.out.predictions,
            FILTER_BAMBU_COUNTS.out.counts_transcript_filter,
            FILTER_BAMBU_COUNTS.out.counts_gene_filter,
            // Bambu's txClassDescription lives only in the SummarizedExperiment,
            // so the RDS comes along to be read for it.
            ASSEMBLY.out.rds_transcript,
            params.annotation,
            file("${projectDir}/bin/novel_transcripts.R", checkIfExists: true),
            file("${projectDir}/bin/gtf_annotation_utils.R", checkIfExists: true)
        )
        ch_versions = ch_versions.mix(NOVEL_TRANSCRIPTS.out.versions)

        VALIDATE_BAMBU_COUNTS (
            NOVEL_TRANSCRIPTS.out.novel_combined_metadata,
            FILTER_BAMBU_COUNTS.out.counts_gene_filter,
            FILTER_BAMBU_COUNTS.out.counts_transcript_filter,
            FILTER_BAMBU_COUNTS.out.cpm_transcript_filter,
            FILTER_BAMBU_COUNTS.out.full_length_counts_transcript_filter,
            FILTER_BAMBU_COUNTS.out.unique_counts_transcript_filter
        )

        REF_TRANSCRIPTS (
            VALIDATE_BAMBU_COUNTS.out.counts_transcript_validated,
            VALIDATE_BAMBU_COUNTS.out.counts_gene_validated,
            ASSEMBLY.out.gtf_all_transcripts,
            params.annotation,
            file("${projectDir}/bin/ref_transcripts.R", checkIfExists: true),
            file("${projectDir}/bin/gtf_annotation_utils.R", checkIfExists: true)
        )
        ch_versions = ch_versions.mix(REF_TRANSCRIPTS.out.versions)

        VALIDATE_BAMBU_GTF (
            ASSEMBLY.out.gtf_all_transcripts,
            VALIDATE_BAMBU_COUNTS.out.counts_transcript_validated,
            VALIDATE_BAMBU_COUNTS.out.full_length_counts_transcript_validated,
            VALIDATE_BAMBU_COUNTS.out.unique_counts_transcript_validated
        )
        ch_versions = ch_versions.mix(VALIDATE_BAMBU_GTF.out.versions)

        //
        // Attach biotype and classification attributes to the validated GTFs
        //
        ANNOTATION (
            VALIDATE_BAMBU_GTF.out.annotations_validated_gtf,
            VALIDATE_BAMBU_GTF.out.fullLength_validated_gtf,
            VALIDATE_BAMBU_GTF.out.uniquelyMapped_validated_gtf,
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
            VALIDATE_BAMBU_COUNTS.out.counts_transcript_validated,
            VALIDATE_BAMBU_COUNTS.out.full_length_counts_transcript_validated,
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
            VALIDATE_BAMBU_COUNTS.out.counts_transcript_validated,
            VALIDATE_BAMBU_COUNTS.out.counts_gene_validated,
            file("${projectDir}/bin/validate_bambu_out.R", checkIfExists: true)
        )
        ch_versions = ch_versions.mix(VALIDATE_BAMBU_OUTPUTS.out.versions)

        RENDER_REPORT (
            VALIDATE_BAMBU_COUNTS.out.counts_gene_validated,
            VALIDATE_BAMBU_COUNTS.out.counts_transcript_validated,
            VALIDATE_BAMBU_COUNTS.out.full_length_counts_transcript_validated,
            VALIDATE_BAMBU_COUNTS.out.unique_counts_transcript_validated,
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
