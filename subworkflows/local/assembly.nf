include { BAMBU } from '../../modules/local/bambu/main'

workflow ASSEMBLY {
    take:
    bams
    reference
    annotation

    main:
    ch_versions              = channel.empty()
    ch_bamlist               = channel.empty()
    ch_bambu_metrics         = channel.empty()
    ch_samp_info             = channel.empty()
    ch_reference             = channel.empty()
    ch_annotation            = channel.empty()
    ch_assembled_new_gtf     = channel.empty()
    ch_assembled_all_gtf     = channel.empty()
    ch_unique_counts         = channel.empty()
    ch_tx_counts             = channel.empty()
    ch_gene_counts           = channel.empty()
    ch_CPM                   = channel.empty()
    ch_full_length           = channel.empty()
    ch_pca                   = channel.empty()
    ch_pca_grouped           = channel.empty()
    ch_h_transcript          = channel.empty()
    ch_h_gene                = channel.empty()
    ch_rds_transcript        = channel.empty()
    ch_rds_gene              = channel.empty()

    // Setting channel for the reference
    ch_reference  = channel.fromPath(reference, checkIfExists: true)
    ch_annotation = channel.fromPath(annotation, checkIfExists: true)

    // Setting TSV file with sample information
    bams
        .map { meta, path -> meta.group + '\t' + path.getName() }
        .collectFile(name: 'sampinfo_samplesheet.tsv', newLine: true, sort: true)
        .set { ch_samp_info }

    // Setting up the BAM list
    bams
        .map { _meta, path -> path.toString() }
        .collectFile(name: 'bamlist.txt', newLine: true, sort: true)
        .set { ch_bamlist }

    BAMBU (
        ch_bamlist,
        ch_reference,
        ch_annotation,
        ch_samp_info,
    )

    BAMBU.out.gtf_new_transcripts
        .set { ch_assembled_new_gtf }

    BAMBU.out.gtf_all_transcripts
        .set { ch_assembled_all_gtf }

    BAMBU.out.gene_counts
        .set { ch_gene_counts }

    BAMBU.out.tx_counts
        .set { ch_tx_counts }

    BAMBU.out.CPM
        .set { ch_CPM }

    BAMBU.out.full_length
        .set { ch_full_length }

    BAMBU.out.unique_counts
        .set { ch_unique_counts }

    BAMBU.out.h_gene
        .set { ch_h_gene }

    BAMBU.out.h_transcript
        .set { ch_h_transcript }

    BAMBU.out.pca
        .set { ch_pca }

    BAMBU.out.pca_grouped
        .set { ch_pca_grouped }

    BAMBU.out.rds_transcript
        .set { ch_rds_transcript }

    BAMBU.out.rds_gene
        .set { ch_rds_gene }

    BAMBU.out.metrics
        .set { ch_bambu_metrics }

    ch_versions = ch_versions.mix(BAMBU.out.versions.ifEmpty(null))

    emit:
    versions            = ch_versions
    pca                 = ch_pca
    pca_grouped         = ch_pca_grouped
    rds_transcript      = ch_rds_transcript
    rds_gene            = ch_rds_gene
    h_gene              = ch_h_gene
    h_transcript        = ch_h_transcript
    CPM                 = ch_CPM
    full_length         = ch_full_length
    transcript_counts   = ch_tx_counts
    gene_counts         = ch_gene_counts
    gtf_new_transcripts = ch_assembled_new_gtf
    gtf_all_transcripts = ch_assembled_all_gtf
    unique_counts       = ch_unique_counts
    bamlist             = ch_bamlist
    samp_info           = ch_samp_info
    reference           = ch_reference
    annotation          = ch_annotation
    bambu_metrics       = ch_bambu_metrics
}