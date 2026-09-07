process VALIDATE_BAMBU_GTF {
    label 'process_low'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://itsiaguara/pulposeq:test':
        'docker.io/itsiaguara/pulposeq:test' }"

    input:
    path gtf_file
    path counts_transcript
    path full_length_counts_transcript
    path unique_counts_transcript

    output:
    path "BambuOutput_annotations_validated.gtf"    , emit: annotations_validated_gtf
    path "BambuOutput_fullLength_validated.gtf"     , emit: fullLength_validated_gtf
    path "BambuOutput_uniquelyMapped_validated.gtf" , emit: uniquelyMapped_validated_gtf
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    validate_bambu_gtf.sh \\
        --gtf ${gtf_file} \\
        --awk_script "\$(which validate_bambu_gtf.awk)" \\
        --counts ${counts_transcript} \\
        --full_length ${full_length_counts_transcript} \\
        --unique ${unique_counts_transcript} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mawk: \$(awk --version | head -n1 | sed 's/mawk //; s/,.*//')
        bash: \$(bash --version | head -n1 | sed 's/GNU bash, version //; s/ .*//')
    END_VERSIONS
    """

    stub:
    """
    touch BambuOutput_annotations_validated.gtf
    touch BambuOutput_fullLength_validated.gtf
    touch BambuOutput_uniquelyMapped_validated.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mawk: \$(awk --version | head -n1 | sed 's/mawk //; s/,.*//')
        bash: \$(bash --version | head -n1 | sed 's/GNU bash, version //; s/ .*//')
    END_VERSIONS
    """
}
