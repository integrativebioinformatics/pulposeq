//
// MODULE: Restrander, from the local modules
//

include { RESTRANDER                } from '../../modules/local/restrander/main'
include { NANOCOMP as NANOCOMP_REST } from '../../modules/nf-core/nanocomp/main'
include { nanocompSkips             } from './utils_nfcore_pulposeq_pipeline'

/*
========================================================================================
    RUN RESTRANDING WORKFLOW
========================================================================================
*/

//
// Fraction of reads Restrander could orient, read back from its JSON report.
// Reads it could not orient are counted under strandStats["?"], so the oriented
// fraction is the + and - counts over the total rather than one minus the other.
//
def orientedFraction(stats_file) {
    // Restrander prints its report on stdout, which the module redirects wholesale.
    // Slicing from the first brace to the last means any progress line the tool
    // writes alongside it does not turn a low orientation rate into a parse error.
    def raw   = stats_file.text
    def first = raw.indexOf('{')
    def last  = raw.lastIndexOf('}')
    if (first < 0 || last < first) {
        error("Could not find a JSON report in ${stats_file}. Restrander may have failed before writing its statistics.")
    }
    def report = new groovy.json.JsonSlurper().parseText(raw.substring(first, last + 1))
    def strand = report?.stats?.strandStats ?: [:]
    def total  = (report?.stats?.totalReads ?: 0) as long
    def placed = ((strand['+'] ?: 0) as long) + ((strand['-'] ?: 0) as long)
    total > 0 ? (placed / total) as double : 0d
}

//
// Orientation correction, deliberately its own subworkflow rather than a step inside
// QC.
//
// It is not quality control and it is not filtering: it is the step that makes read
// strand equal RNA strand, which everything downstream is built on. Keeping it
// separate is what allows --skip_chopper and --skip_nanocomp to remain freely usable
// while leaving no route to an unoriented alignment. Neither turns Restrander off:
// --skip_nanocomp restrander drops the report on its output, not the orientation.
//
// Ordering matters and is expressed by where this sits in the caller: Restrander
// locates primers and polyA/polyT tails at the read ends, exactly the positions
// Chopper's headcrop and tailcrop remove, so it must run first.
//
workflow RESTRANDING {
    take:
    reads
    declared_stranded    // what the user asserted about the reads on disk

    main:
    ch_versions      = channel.empty()
    ch_stats         = channel.empty()
    ch_multiqc       = channel.empty()
    ch_combined_rest = channel.empty()
    ch_reads         = reads

    def nanocomp_skips = nanocompSkips()

    if (params.library == 'ONT_cDNA' && !declared_stranded) {

        if (params.headcrop > 0 || params.tailcrop > 0) {
            log.warn "Restranding is active but headcrop=${params.headcrop} / tailcrop=${params.tailcrop} are non-zero. " +
                     "Primer and adapter removal is Restrander's job; cropping read ends afterwards is fine, but " +
                     "these values are normally 0 for restranded runs."
        }

        RESTRANDER (
            reads,
            params.restrand_config ? file(params.restrand_config, checkIfExists: true) : []
        )
        ch_versions = ch_versions.mix(RESTRANDER.out.versions)
        ch_stats    = RESTRANDER.out.stats

        // Every sample must clear the threshold. There is no unstranded fallback:
        // losing strand erases mono-exonic novel transcripts, pushes novel isoforms
        // of known genes into the antisense class, and bleeds reads between
        // overlapping sense/antisense pairs -- all of it landing on exactly what
        // this pipeline exists to measure. A warning would let that through into
        // results nobody re-reads the log for, so this stops instead.
        //
        // Checked here rather than before Bambu: Restrander takes minutes and
        // Chopper hours per sample, so failing now costs a coffee break rather than
        // most of a day.
        ch_reads = RESTRANDER.out.reads
            .join(RESTRANDER.out.stats)
            .map { meta, restranded, stats ->
                [ meta, restranded, orientedFraction(stats) ]
            }
            .toList()
            .flatMap { rows ->
                def failed = rows.findAll { row -> row[2] < params.restrand_min_frac }
                if (failed) {
                    def detail = rows
                        .collect { row -> "  ${row[0].id}: ${String.format('%.1f%%', row[2] * 100)}" }
                        .join('\n')
                    error(
                        "Restranding did not reach --restrand_min_frac " +
                        "(${String.format('%.1f%%', params.restrand_min_frac * 100)}) for " +
                        "${failed.size()} of ${rows.size()} samples, so the reads cannot be treated as oriented.\n" +
                        "${detail}\n" +
                        "The usual cause is a configuration that does not match the reads. If demultiplexing trimmed " +
                        "the primers, no kit preset can find them: use --restrand_kit trimmed, which orients on " +
                        "polyA/polyT alone. Otherwise check that --restrand_kit matches the library chemistry.\n" +
                        "If a single sample is far below the rest, that library likely has a problem and is better " +
                        "excluded than accommodated by lowering the threshold."
                    )
                }
                rows.collect { row -> [ row[0] + [restrand_frac: row[2]], row[1] ] }
            }

        // Read report on the oriented reads. Fed from the checked channel rather
        // than RESTRANDER.out.reads directly, so a run that fails the orientation
        // threshold stops without spending a NanoComp on reads it is about to
        // reject. It sits between the raw and the Chopper reports: what changed
        // here is which strand a read is on, not how many reads or how long they
        // are, so a difference against the raw report is reads Restrander dropped.
        if (!('restrander' in nanocomp_skips)) {

            ch_reads
                .collect { item -> item[1] }
                .map { filelist -> [[id: "All"], filelist] }
                .set { ch_combined_rest }

            NANOCOMP_REST(ch_combined_rest)

            ch_multiqc  = ch_multiqc.mix(NANOCOMP_REST.out.stats_txt.collect { item -> item[1] }.ifEmpty([]))
            ch_versions = ch_versions.mix(NANOCOMP_REST.out.versions)
        }

    } else {
        ch_reads = reads.map { meta, fastq -> [ meta + [restrand_frac: null], fastq ] }
    }

    emit:
    reads    = ch_reads
    stats    = ch_stats
    multiqc  = ch_multiqc
    versions = ch_versions
}
