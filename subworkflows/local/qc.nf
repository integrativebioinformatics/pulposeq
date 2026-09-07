//
// MODULE: Installed directly from nf-core/modules
//

include { NANOCOMP as NANOCOMP_RAW          } from '../../modules/nf-core/nanocomp/main'
include { NANOCOMP as NANOCOMP_CHOPPER      } from '../../modules/nf-core/nanocomp/main'
include { CHOPPER                           } from '../../modules/nf-core/chopper/main'
include { nanocompSkips                     } from './utils_nfcore_pulposeq_pipeline'

/*
========================================================================================
    RUN QC_FILT WORKFLOW
========================================================================================
*/

workflow QC_FILT {
    take:
    // Two read channels rather than one. NanoComp's raw-read figures describe the
    // reads as they were supplied, while filtering operates on whatever RESTRANDING
    // produced -- which for unoriented ONT cDNA is a reverse-complemented set minus
    // the reads that could not be oriented. When restranding does not apply the two
    // are the same channel.
    raw_reads
    reads

    main:
    ch_versions          = channel.empty()
    ch_multiqc_raw       = channel.empty()
    ch_multiqc_filt      = channel.empty()
    ch_multiqc_all       = channel.empty()
    ch_combined_raw      = channel.empty()
    ch_combined_filtered = channel.empty()
    ch_reads             = reads

    // Reporting and filtering are gated separately: --skip_nanocomp drops a report
    // without touching the reads that reach alignment, while --skip_chopper changes
    // those reads. Filtering therefore still runs with every report turned off, and
    // every report still runs with filtering turned off.
    def nanocomp_skips = nanocompSkips()

    // Running nanocomp on raw reads
    if (!('raw' in nanocomp_skips)) {

        raw_reads
            .collect { item -> item[1] }
            .map { filelist -> [[id: "All"], filelist] }
            .set { ch_combined_raw }

        NANOCOMP_RAW(ch_combined_raw)

        ch_versions = ch_versions.mix(NANOCOMP_RAW.out.versions)

        // Generating a multiqc file for raw reads report
        ch_multiqc_raw = ch_multiqc_raw.mix(NANOCOMP_RAW.out.stats_txt.collect { item -> item[1] }.ifEmpty([]))

        ch_multiqc_all = ch_multiqc_all.mix(ch_multiqc_raw.ifEmpty([]))
    }

    // Putting conditional to whether run filtering on samples
    if (!params.skip_chopper) {

        CHOPPER(ch_reads, [])

        // Removed the broken CHOPPER.out.versions line entirely

        // Running quality check in filtered reads
        if (!('chopper' in nanocomp_skips)) {

            CHOPPER.out.fastq
                .collect { item -> item[1] }
                .map { filelist -> [[id: "All"], filelist] }
                .set { ch_combined_filtered }

            NANOCOMP_CHOPPER(ch_combined_filtered)

            ch_multiqc_filt = ch_multiqc_filt.mix(NANOCOMP_CHOPPER.out.stats_txt.collect { item -> item[1] }.ifEmpty([]))

            ch_multiqc_all = ch_multiqc_all.mix(ch_multiqc_filt.ifEmpty([]))
        }

        // Putting the output of chopper as new ch_reads
        ch_reads = CHOPPER.out.fastq
    }

    emit:
    filt_reads = ch_reads
    multiqc    = ch_multiqc_all
    versions   = ch_versions
}
