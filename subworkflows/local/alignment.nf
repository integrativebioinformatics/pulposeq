//
// MODULES: MINIMAP2_ALIGN and NANOCOMP
//
include { MINIMAP2_ALIGN                  } from '../../modules/nf-core/minimap2/align/main'
include { NANOCOMP as NANOCOMP_MINIMAP2   } from '../../modules/nf-core/nanocomp/main'
include { nanocompSkips                   } from './utils_nfcore_pulposeq_pipeline'

/*
========================================================================================
    RUN ALIGNMENT WORKFLOW
========================================================================================
*/

workflow ALIGNMENT {
    take:
    reads

    main:
    ch_versions         = channel.empty()
    ch_bam              = channel.empty()
    ch_index            = channel.empty()
    ch_combined_mapping = channel.empty()
    ch_alignment_qc     = channel.empty()
    ch_reference        = channel.empty()

    // Building metamap for the reference
    channel
        .fromPath(params.reference)
        .map { file_path ->
            def basename = file_path.baseName
            [basename, file_path.toString()]
        }
        .collect()
        .set { ch_reference }

    // Alignment with the minimap2 module
    MINIMAP2_ALIGN (
        reads,
        ch_reference,
        params.bam_format,
        params.bam_index_extension,
        params.cigar_paf_format,
        params.cigar_bam
    )

    MINIMAP2_ALIGN.out.bam
        .set { ch_bam }

    // The .bai files were always produced and published, but ch_index was left as
    // the empty channel declared above, so nothing downstream could consume them.
    // An empty channel yields zero tasks rather than an error, so any process fed
    // from here would have been skipped silently.
    MINIMAP2_ALIGN.out.index
        .set { ch_index }

    // The alignment report is selected through --skip_nanocomp minimap2, alongside
    // the read-level ones, rather than through a flag of its own: it is the same
    // tool answering the same question one step later, and it is the most expensive
    // of the four, so it is the one most often worth dropping on large runs.
    if (!('minimap2' in nanocompSkips())) {

        ch_bam
            .collect { item -> item[1] }
            .map { filelist -> [[id: "All"], filelist] }
            .set { ch_combined_mapping }

        NANOCOMP_MINIMAP2 (
            ch_combined_mapping
        )

        ch_alignment_qc = ch_alignment_qc.mix(NANOCOMP_MINIMAP2.out.stats_txt.collect { item -> item[1] }.ifEmpty([]))

        ch_versions = ch_versions.mix(NANOCOMP_MINIMAP2.out.versions.ifEmpty(null))
    }

    emit:
    index    = ch_index
    bam      = ch_bam
    reference = ch_reference
    multiqc  = ch_alignment_qc
    versions = ch_versions
}