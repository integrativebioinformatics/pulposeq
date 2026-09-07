process VALIDATE_BAMBU_COUNTS {
    label 'process_low'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://itsiaguara/pulposeq:test':
        'docker.io/itsiaguara/pulposeq:test' }"

    input:
    path metadata_csv
    path counts_gene
    path counts_transcript
    path cpm_transcript
    path full_length_counts_transcript
    path unique_counts_transcript

    output:
    path "BambuOutput_counts_gene_validated.txt"                 , emit: counts_gene_validated
    path "BambuOutput_counts_transcript_validated.txt"           , emit: counts_transcript_validated
    path "BambuOutput_CPM_transcript_validated.txt"              , emit: cpm_transcript_validated
    path "BambuOutput_fullLengthCounts_transcript_validated.txt" , emit: full_length_counts_transcript_validated
    path "BambuOutput_uniqueCounts_transcript_validated.txt"     , emit: unique_counts_transcript_validated
    path "versions.yml"                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    validate_bambu_counts.sh \\
        --metadata ${metadata_csv} \\
        --counts_gene ${counts_gene} \\
        --counts_transcript ${counts_transcript} \\
        --cpm_transcript ${cpm_transcript} \\
        --full_length ${full_length_counts_transcript} \\
        --unique ${unique_counts_transcript} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mawk: \$(awk --version | head -n1 | sed 's/GNU Awk //; s/,.*//')
        bash: \$(bash --version | head -n1 | sed 's/GNU bash, version //; s/ .*//')
    END_VERSIONS
    """

    stub:
    """
    touch BambuOutput_counts_gene_validated.txt
    touch BambuOutput_counts_transcript_validated.txt
    touch BambuOutput_CPM_transcript_validated.txt
    touch BambuOutput_fullLengthCounts_transcript_validated.txt
    touch BambuOutput_uniqueCounts_transcript_validated.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mawk: \$(awk --version | head -n1 | sed 's/GNU Awk //; s/,.*//')
        bash: \$(bash --version | head -n1 | sed 's/GNU bash, version //; s/ .*//')
    END_VERSIONS
    """
}
