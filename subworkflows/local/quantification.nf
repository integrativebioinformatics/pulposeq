//
// oarfish quantification: optional, with --run_oarfish true
//
include { OARFISH_INDEX_FASTA } from '../../modules/local/oarfish_index_fasta/main'
include { OARFISH_READMODE    } from '../../modules/nf-core/oarfish/readmode/main'
include { OARFISH_MERGE       } from '../../modules/local/oarfish_merge/main'
include { oarfishSeqTech      } from './utils_nfcore_pulposeq_pipeline'

/*
========================================================================================
    RUN QUANTIFICATION WORKFLOW
========================================================================================
*/

//
// oarfish quantifies against every model Bambu assembled, before any filtering, so
// reads from models PREFILTER later removes stay on those models instead of being
// pushed onto the transcripts that are kept. The one exception is unstranded novel
// models, which oarfish's forward-strand filter cannot use reliably.
//
workflow QUANTIFICATION {
    take:
    reads          // channel: [ meta, fastq ] -- the same reads minimap2 aligns
    extended_gtf   // channel: BambuOutput_extended_annotations.gtf
    reference      // path or string: genome FASTA (params.reference)
    bambu_counts   // channel: BambuOutput_counts_transcript.txt

    main:
    ch_versions = channel.empty()

    OARFISH_INDEX_FASTA (
        extended_gtf,
        file(reference, checkIfExists: true)
    )
    ch_versions = ch_versions.mix(OARFISH_INDEX_FASTA.out.versions)

    // One index for every sample: value channels, so each sample's task reuses them
    def ch_annotated = OARFISH_INDEX_FASTA.out.annotated
        .map { fasta -> [ [id: 'annotated'], fasta ] }
        .first()

    // An empty novel FASTA (no novel models) becomes an empty input, and the module
    // then leaves out --novel
    def ch_novel = OARFISH_INDEX_FASTA.out.novel
        .map { fasta -> [ [id: 'novel'], fasta.size() > 0 ? fasta : [] ] }
        .first()

    def seq_tech = oarfishSeqTech()

    OARFISH_READMODE (
        reads.map { meta, fastq -> [ meta, fastq, seq_tech ] },
        ch_annotated,
        ch_novel
    )
    // OARFISH_READMODE reports its version through the `versions` topic channel

    OARFISH_MERGE (
        OARFISH_READMODE.out.quant.map     { _meta, file -> file }.toList(),
        OARFISH_READMODE.out.infreps.map   { _meta, file -> file }.toList(),
        OARFISH_READMODE.out.meta_info.map { _meta, file -> file }.toList(),
        bambu_counts,
        OARFISH_INDEX_FASTA.out.summary,
        params.library,
        file("${projectDir}/bin/oarfish_merge.R", checkIfExists: true)
    )
    ch_versions = ch_versions.mix(OARFISH_MERGE.out.versions)

    emit:
    counts_transcript = OARFISH_MERGE.out.counts_transcript   // all models: PREFILTER's expression source
    counts_gene       = OARFISH_MERGE.out.counts_gene
    cpm_transcript    = OARFISH_MERGE.out.cpm_transcript
    tpm_transcript    = OARFISH_MERGE.out.tpm_transcript
    bootstraps        = OARFISH_MERGE.out.bootstraps
    summary           = OARFISH_MERGE.out.summary
    versions          = ch_versions
}
