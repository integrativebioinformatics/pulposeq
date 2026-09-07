process RNAMINING {
    tag 'Coding potential prediction'
    label 'process_low'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://samuelismael/rnamining:1.1.0-nextflow':
        'docker.io/samuelismael/rnamining:1.1.0-nextflow' }"

    input:
    path fasta

    output:
    path 'codings.txt'        , emit: coding
    path 'noncodings.txt'     , emit: noncoding
    path 'predictions.txt'    , emit: preds
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ? task.ext.args : "-organism_name ${params.organism} -prediction_type coding_prediction"
    """
    rnamining \\
            -f $fasta \\
            $args \\
            -output_folder ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rnamining: \$(rnamining --version 2>&1 | sed 's/^rnamining[[:space:]]*//')
    END_VERSIONS
    """

    stub:
    """
    touch codings.txt
    touch noncodings.txt
    cat <<-END_PREDICTIONS > predictions.txt
    # RNAmining predictions
    # ID	prediction	score
    # RNAmining 1.1.0
    # FASTA transcript predictions
    transcript_stub	coding	0.5
    END_PREDICTIONS
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rnamining: 1.1.0
    END_VERSIONS
    """
}
