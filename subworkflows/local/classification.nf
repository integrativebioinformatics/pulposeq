//
// MODULE: nf-core modules
//
include { GFFCOMPARE } from '../../modules/nf-core/gffcompare/main'
include { GFFREAD    } from '../../modules/nf-core/gffread/main'
include { CPC2       } from '../../modules/nf-core/cpc2/main'

/*
========================================================================================
    RUN CLASSIFICATION WORKFLOW
========================================================================================
*/

workflow CLASSIFICATION {
    take:
    gtf         // path: GTF from transcript reconstruction
    annotation  // path/string: reference annotation GTF (params.annotation)
    reference   // path/string: reference genome FASTA (params.reference)

    main:
    //
    // Prepare inputs for nf-core/gffcompare
    // nf-core module expects: tuple(meta, gtfs), tuple(meta2, fasta, fai), tuple(meta3, reference_gtf)
    //
    def ch_gffcompare_gtf = gtf.map { gtf_file ->
        [ [id: 'gffcompare'], gtf_file ]
    }
    def ch_gffcompare_fasta = [ [id: 'genome'], [], [] ]  // not using fasta input
    def ch_gffcompare_ref   = [ [id: 'reference'], file(annotation) ]

    GFFCOMPARE(
        ch_gffcompare_gtf,
        ch_gffcompare_fasta,
        ch_gffcompare_ref
    )

    // Extract bare paths from tuple outputs for downstream compatibility
    ch_annotated_gtf = GFFCOMPARE.out.annotated_gtf.map { _meta, f -> f }
    ch_tmap          = GFFCOMPARE.out.tmap.map { _meta, f -> f }

    //
    // Prepare inputs for nf-core/gffread
    // nf-core module expects: tuple(meta, gff), path(fasta)
    //
    def ch_gffread_input = gtf.map { gtf_file ->
        [ [id: gtf_file.baseName], gtf_file ]
    }

    GFFREAD(
        ch_gffread_input,
        reference
    )

    // Extract bare fasta path for downstream compatibility
    ch_gffread_fasta = GFFREAD.out.gffread_fasta.map { _meta, f -> f }

    //
    // Coding potential prediction with CPC2: a table of one row per transcript with a
    // coding/non-coding call. novel_transcripts.R keys on its header rather than on
    // column position.
    //
    CPC2(
        GFFREAD.out.gffread_fasta
    )
    ch_predictions = CPC2.out.txt.map { _meta, f -> f }

    // GFFCOMPARE, GFFREAD and CPC2 all report versions through the `versions` topic
    // channel rather than a versions.yml emit, so this subworkflow emits none.

    emit:
    annotated_gtf = ch_annotated_gtf
    tmap          = ch_tmap
    gffread_fasta = ch_gffread_fasta
    predictions   = ch_predictions
}
