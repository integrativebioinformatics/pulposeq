//
// Subworkflow with functionality specific to the integrativebioinformatics/pulposeq pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    def before_text = ""
    def after_text = ""
    before_text = """
-\033[2m----------------------------------------------------\033[0m-
\033[0;35m  integrativebioinformatics/pulposeq ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/','')}"}.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/integrativebioinformatics/pulposeq/blob/main/CITATIONS.md
"""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN(
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command,
        null
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Create channel from input file provided through params.input
    //
    ch_samplesheet = channel.fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    _outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    _multiqc_report  //  string: Path to MultiQC report

    main:

    //
    // Completion summary
    //
    workflow.onComplete {
        completionSummary(monochrome_logs)
    }

    workflow.onError {
        log.error("Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting")
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Stages whose NanoComp report --skip_nanocomp turns off.
//
// One parameter rather than one boolean per stage: the four reports are the same
// tool run at four points of the same pipeline, and what a run actually decides is
// where read-level QC is worth its runtime. Returns the set of stage names to skip,
// so call sites read as `if (!('raw' in nanocompSkips()))`.
//
// 'all' expands to every stage. Tokens may be separated by commas or whitespace,
// and 'rest' is accepted for 'restrander' because it is how the stage is usually
// spoken about and reads evenly next to 'raw' and 'chopper'.
//
def nanocompSkips() {
    def stages = ['raw', 'restrander', 'chopper', 'minimap2']

    def requested = nanocompSkipTokens()
        .collect { token -> token == 'rest' ? 'restrander' : token }
        .unique()

    def unknown = requested - (stages + 'all')
    if (unknown) {
        error(
            "--skip_nanocomp does not recognise ${unknown.collect { token -> "'${token}'" }.join(', ')}. " +
            "Valid values are ${stages.join(', ')} and all, alone or comma-separated " +
            "(for example --skip_nanocomp 'raw,chopper'). Leave it unset to run every report."
        )
    }

    return (requested.contains('all') ? stages : requested) as Set
}

//
// The raw tokens as the user wrote them, before 'all' is expanded. Only useful for
// telling "they asked for this stage by name" from "they asked for everything".
//
def nanocompSkipTokens() {
    return params.skip_nanocomp
        ?.toString()
        ?.toLowerCase()
        ?.tokenize(', \t')
        ?.findAll { token -> token } ?: []
}

//
// Whether Restrander runs: ONT cDNA that the user has not declared oriented.
//
def restrandingActive() {
    return params.library == 'ONT_cDNA' &&
           params.stranded_library?.toString()?.toLowerCase() != 'true'
}

//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    // NanoComp report selection. Parsed here as well as at the call sites so a typo
    // fails at startup rather than hours in, where the only symptom would be a
    // report quietly missing from the MultiQC page.
    nanocompSkips()

    if (nanocompSkipTokens().any { token -> token in ['rest', 'restrander'] } && !restrandingActive()) {
        def why = params.library == 'ONT_cDNA'
            ? "--stranded_library is true, so these reads are already oriented"
            : "--library ${params.library} is always oriented"

        log.warn(
            "--skip_nanocomp names the restrander report, but restranding is not running: ${why}. " +
            "That report only exists for ONT_cDNA reads the pipeline has to orient itself, so the " +
            "token has no effect."
        )
    }

    // Alignment requires reference genome and annotation
    if (!params.skip_minimap2) {
        if (!params.reference) {
            error("--reference must be provided when alignment is not skipped (skip_minimap2 = false)")
        }
        if (!params.annotation) {
            error("--annotation must be provided when alignment is not skipped (skip_minimap2 = false)")
        }

        // Library type drives the minimap2 preset and Bambu strandedness
        def valid_libraries = ['ONT_cDNA', 'ONT_DRS', 'PacBio']
        if (!params.library) {
            error("--library must be provided when alignment is not skipped (skip_minimap2 = false). Valid values: ${valid_libraries.join(', ')}")
        }
        if (!(params.library in valid_libraries)) {
            error("--library '${params.library}' is not valid. Valid values: ${valid_libraries.join(', ')}")
        }
        if (params.library in ['ONT_DRS', 'PacBio'] && params.stranded_library?.toString()?.toLowerCase() != 'true') {
            log.warn("--library ${params.library} is always stranded; --stranded_library false will be ignored.")
        }

        // Coding potential predictor. RNAmining carries a per-organism model and needs
        // to be told which, where CPC2 does not.
        def valid_predictors = ['cpc2', 'rnamining']
        if (!(params.coding_potential_pred in valid_predictors)) {
            error("--coding_potential_pred '${params.coding_potential_pred}' is not valid. Valid values: ${valid_predictors.join(', ')}")
        }
        if (params.coding_potential_pred == 'rnamining' && !params.organism) {
            error(
                "--organism must be provided when --coding_potential_pred rnamining. " +
                "RNAmining selects a per-organism model and has no usable default.\n" +
                "Supply --organism, or use the default --coding_potential_pred cpc2, which needs none."
            )
        }
        if (params.coding_potential_pred == 'cpc2' && params.organism) {
            log.warn("--organism is only used by RNAmining and will be ignored with --coding_potential_pred cpc2.")
        }

        // splice:hq is `splice` with -C5 -O6,24 -B4: it trusts the read, treating a
        // discrepancy as real biology rather than basecalling error. -k14 exists for
        // the opposite reason, so the two cannot both be right about the same data.
        if (params.map_hq != null && params.library == 'ONT_DRS') {
            error(
                "--map_hq cannot be used with --library ONT_DRS. splice:hq assumes a low error " +
                "rate, while the ONT_DRS preset sets -k14 to cope with a high one, so the two " +
                "cannot both be right about the same reads.\n" +
                "Drop --map_hq, or use --library ONT_cDNA if these are cDNA reads."
            )
        }
        // Left unset, the preset follows the library. Say which was chosen, since it
        // changes junction placement and therefore the novel transcript calls, and a
        // reader of the log should not have to infer it from the library type.
        if (params.map_hq == null) {
            log.info("--map_hq unset: using minimap2 " +
                     (params.library == 'PacBio' ? 'splice:hq' : 'splice') +
                     " for --library ${params.library}.")
        }
        else if (params.map_hq && params.library == 'PacBio') {
            log.info("--map_hq true for --library PacBio matches the default preset; no change.")
        }
        else if (!params.map_hq && params.library == 'PacBio') {
            log.warn("--map_hq false overrides PacBio's default splice:hq preset with plain splice. " +
                     "splice:hq is minimap2's documented preset for Iso-Seq reads.")
        }

        // Restranding applies only to ONT cDNA that is not already oriented. The
        // pipeline admits no unstranded path: losing strand erases mono-exonic novel
        // transcripts, depletes novel isoforms of known genes into the antisense
        // class, and bleeds reads between overlapping sense/antisense pairs. Such
        // reads are either oriented before they arrive or oriented here.
        def restrand_active = restrandingActive()

        if (restrand_active) {
            // The kit is deliberately not defaulted: the presets differ in their
            // TSO/RTP primer sequences, and the wrong one yields a low orientation
            // rate rather than an error.
            if (!params.restrand_kit && !params.restrand_config) {
                def kits = ['PCB109', 'PCB111', 'PCB114', 'DCS109', 'DCS-LSK114', 'NEBNext', 'trimmed']
                error(
                    "--restrand_kit must be provided for unoriented ONT_cDNA libraries. Valid values: ${kits.join(', ')}.\n" +
                    "Alternatively pass --restrand_config with a custom Restrander JSON.\n" +
                    "If the reads were already oriented outside the pipeline, set --stranded_library true instead."
                )
            }

            // A floor, not just a range check. Without one, --restrand_min_frac 0
            // would wave through reads Restrander could not orient while minimap2
            // and Bambu are both told the data is stranded -- worse than the
            // unstranded path this replaced, because it fails silently.
            def frac = params.restrand_min_frac as double
            if (frac < 0.5 || frac > 1) {
                error(
                    "--restrand_min_frac must be between 0.5 and 1, got ${params.restrand_min_frac}. " +
                    "Reads that Restrander could not orient must not reach an alignment and a quantification " +
                    "step that assume strand is known, so this threshold cannot be disabled."
                )
            }
            if (frac < 0.75) {
                log.warn "--restrand_min_frac is set to ${params.restrand_min_frac}. Restrander reaches ~0.99 on " +
                         "matched PCR-cDNA data, so a threshold this low usually means the kit or the trimming " +
                         "state is wrong rather than that the threshold needs relaxing."
            }
        }
        else if (params.restrand_kit || params.restrand_config) {
            // Reached only when one of these was actually passed: both default to
            // null, so their presence here means they were set for a library that
            // skips restranding. Warned about rather than quietly dropped -- the run
            // record should reflect what was asked for, and hiding the value would
            // leave the misunderstanding invisible.
            //
            // restrand_min_frac is deliberately not covered: it defaults to 0.8, so
            // there is no way to tell a user who passed 0.8 from one who passed
            // nothing at all.
            def why = params.library == 'ONT_cDNA'
                ? "--stranded_library is true, so these reads are already oriented"
                : "--library ${params.library} is always oriented"

            log.warn(
                "restrand_kit / restrand_config are set, but restranding is neither running nor applied: " +
                "${why}. Restranding runs only for ONT_cDNA reads that are not already oriented, so these " +
                "values are ignored."
            )
        }
    }

    // Filtering length checks
    if (params.maxlen && params.minlen && (params.maxlen as int) < (params.minlen as int)) {
        error("--maxlen (${params.maxlen}) cannot be less than --minlen (${params.minlen})")
    }

    // GC content range check
    if ((params.mingc as double) < 0 || (params.mingc as double) > 1) {
        error("--mingc must be between 0 and 1, got: ${params.mingc}")
    }
    if ((params.maxgc as double) < 0 || (params.maxgc as double) > 1) {
        error("--maxgc must be between 0 and 1, got: ${params.maxgc}")
    }
    if ((params.mingc as double) > (params.maxgc as double)) {
        error("--mingc (${params.mingc}) cannot be greater than --maxgc (${params.maxgc})")
    }

    // Crop values must be non-negative
    if ((params.headcrop as int) < 0) {
        error("--headcrop must be >= 0, got: ${params.headcrop}")
    }
    if ((params.tailcrop as int) < 0) {
        error("--tailcrop must be >= 0, got: ${params.tailcrop}")
    }

    // NDR must be between 0 and 1 if provided
    if (params.ndr != null) {
        if ((params.ndr as double) < 0 || (params.ndr as double) > 1) {
            error("--ndr must be between 0 and 1, got: ${params.ndr}")
        }
    }
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    def citation_text = [
            "Tools used in the workflow included:",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
