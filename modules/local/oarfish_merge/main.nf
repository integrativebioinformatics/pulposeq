process OARFISH_MERGE {
    tag "all samples"

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://itsiaguara/pulposeq:test':
        'docker.io/itsiaguara/pulposeq:test' }"

    input:
    path quants, stageAs: 'quant/*'           // every sample's <id>.quant
    path infreps, stageAs: 'infreps/*'        // every sample's <id>.infreps.pq; empty list with 0 bootstraps
    path meta_info, stageAs: 'meta_info/*'    // every sample's <id>.meta_info.json
    path bambu_counts                         // BambuOutput_counts_transcript.txt
    path index_summary                        // oarfish_index_summary.tsv
    val  library
    path r_script

    output:
    path "oarfish_counts_transcript.txt" , emit: counts_transcript
    path "oarfish_counts_gene.txt"       , emit: counts_gene
    path "oarfish_CPM_transcript.txt"    , emit: cpm_transcript
    path "oarfish_TPM_transcript.txt"    , emit: tpm_transcript, optional: true
    path "oarfish_bootstraps.rds"        , emit: bootstraps    , optional: true
    path "oarfish_meta_info.tsv"         , emit: meta_info     , optional: true
    path "oarfish_summary.tsv"           , emit: summary
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # The script reads every oarfish file from the working directory
    for f in quant/* infreps/* meta_info/*; do
        [ -e "\$f" ] && ln -s "\$f" "\$(basename "\$f")"
    done

    Rscript $r_script \\
        --bambu_counts ${bambu_counts} \\
        --index_summary ${index_summary} \\
        --library ${library} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | sed 's/R version //; s/ (.*//' | head -1)
        r-optparse: \$(Rscript -e "cat(as.character(packageVersion('optparse')))")
        r-arrow: \$(Rscript -e "cat(if (requireNamespace('arrow', quietly = TRUE)) as.character(packageVersion('arrow')) else 'not installed')")
    END_VERSIONS
    """

    stub:
    """
    touch oarfish_counts_transcript.txt oarfish_counts_gene.txt oarfish_CPM_transcript.txt oarfish_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | sed 's/R version //; s/ (.*//' | head -1)
    END_VERSIONS
    """
}
