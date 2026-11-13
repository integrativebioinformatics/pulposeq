process CURATE {
    label 'process_low'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://itsiaguara/pulposeq:test':
        'docker.io/itsiaguara/pulposeq:test' }"

    input:
    path novel_metadata                 // novel_transcripts_validated_metadata.csv
    path counts_transcript              // PREFILTER's tables
    path counts_gene
    path cpm_transcript
    path full_length_counts_transcript
    path unique_counts_transcript
    path extended_gtf                   // BambuOutput_extended_annotations.gtf
    val  expression_source              // 'bambu' or 'oarfish'
    path oarfish_counts                 // the oarfish inputs are empty lists unless oarfish ran
    path oarfish_counts_gene
    path oarfish_cpm
    path oarfish_tpm                    // empty for ONT libraries
    path oarfish_bootstraps             // empty with --oarfish_num_bootstraps 0
    path r_script

    output:
    path "BambuOutput_counts_transcript_validated.txt"           , emit: counts_transcript_validated
    path "BambuOutput_counts_gene_validated.txt"                 , emit: counts_gene_validated
    path "BambuOutput_CPM_transcript_validated.txt"              , emit: cpm_transcript_validated
    path "BambuOutput_fullLengthCounts_transcript_validated.txt" , emit: full_length_counts_transcript_validated
    path "BambuOutput_uniqueCounts_transcript_validated.txt"     , emit: unique_counts_transcript_validated
    path "BambuOutput_annotations_validated.gtf"                 , emit: annotations_validated_gtf
    path "BambuOutput_fullLength_validated.gtf"                  , emit: fullLength_validated_gtf
    path "BambuOutput_uniquelyMapped_validated.gtf"              , emit: uniquelyMapped_validated_gtf
    path "transcript_read_support.tsv"                           , emit: read_support
    path "oarfish_counts_transcript_validated.txt"               , emit: oarfish_counts_transcript, optional: true
    path "oarfish_counts_gene_validated.txt"                     , emit: oarfish_counts_gene      , optional: true
    path "oarfish_CPM_transcript_validated.txt"                  , emit: oarfish_cpm_transcript   , optional: true
    path "oarfish_TPM_transcript_validated.txt"                  , emit: oarfish_tpm_transcript   , optional: true
    path "oarfish_validated_se.rds"                              , emit: oarfish_se               , optional: true
    path "versions.yml"                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def oarfish_args = [
        oarfish_counts      ? "--oarfish_counts ${oarfish_counts}"           : '',
        oarfish_counts_gene ? "--oarfish_counts_gene ${oarfish_counts_gene}" : '',
        oarfish_cpm         ? "--oarfish_cpm ${oarfish_cpm}"                 : '',
        oarfish_tpm         ? "--oarfish_tpm ${oarfish_tpm}"                 : '',
        oarfish_bootstraps  ? "--oarfish_bootstraps ${oarfish_bootstraps}"   : ''
    ].findAll { arg -> arg }.join(' ')
    """
    Rscript $r_script \\
        --novel_metadata ${novel_metadata} \\
        --counts_transcript ${counts_transcript} \\
        --counts_gene ${counts_gene} \\
        --cpm_transcript ${cpm_transcript} \\
        --full_length ${full_length_counts_transcript} \\
        --unique ${unique_counts_transcript} \\
        --extended_gtf ${extended_gtf} \\
        --expression_source ${expression_source} \\
        ${oarfish_args} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | sed 's/R version //; s/ (.*//' | head -1)
        r-optparse: \$(Rscript -e "cat(as.character(packageVersion('optparse')))")
        bioconductor-summarizedexperiment: \$(Rscript -e "cat(as.character(packageVersion('SummarizedExperiment')))")
    END_VERSIONS
    """

    stub:
    """
    touch BambuOutput_counts_transcript_validated.txt BambuOutput_counts_gene_validated.txt
    touch BambuOutput_CPM_transcript_validated.txt BambuOutput_fullLengthCounts_transcript_validated.txt
    touch BambuOutput_uniqueCounts_transcript_validated.txt
    touch BambuOutput_annotations_validated.gtf BambuOutput_fullLength_validated.gtf BambuOutput_uniquelyMapped_validated.gtf
    touch transcript_read_support.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | sed 's/R version //; s/ (.*//' | head -1)
    END_VERSIONS
    """
}
