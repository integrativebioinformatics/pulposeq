process PREFILTER {
    label 'process_low'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://itsiaguara/pulposeq:test':
        'docker.io/itsiaguara/pulposeq:test' }"

    input:
    path counts_transcript
    path counts_gene
    path cpm_transcript
    path full_length_counts_transcript
    path unique_counts_transcript
    path novel_gtf
    val  expression_source       // 'bambu' or 'oarfish'
    path oarfish_counts          // [] unless expression_source is 'oarfish'
    val  min_total_counts
    path r_script

    output:
    path "BambuOutput_counts_transcript_filter.txt"           , emit: counts_transcript_filter
    path "BambuOutput_counts_gene_filter.txt"                 , emit: counts_gene_filter
    path "BambuOutput_CPM_transcript_filter.txt"              , emit: cpm_transcript_filter
    path "BambuOutput_fullLengthCounts_transcript_filter.txt" , emit: full_length_counts_transcript_filter
    path "BambuOutput_uniqueCounts_transcript_filter.txt"     , emit: unique_counts_transcript_filter
    path "bambu_novel_transcripts_prefiltered.gtf"            , emit: novel_gtf
    path "prefilter_kept_transcripts.txt"                     , emit: kept_transcripts
    path "filter_summary.tsv"                                 , emit: summary
    path "versions.yml"                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def oarfish_arg = expression_source == 'oarfish' ? "--oarfish_counts ${oarfish_counts}" : ''
    """
    Rscript $r_script \\
        --counts_transcript ${counts_transcript} \\
        --counts_gene ${counts_gene} \\
        --cpm_transcript ${cpm_transcript} \\
        --full_length ${full_length_counts_transcript} \\
        --unique ${unique_counts_transcript} \\
        --novel_gtf ${novel_gtf} \\
        --expression_source ${expression_source} \\
        ${oarfish_arg} \\
        --min_total_counts ${min_total_counts} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | sed 's/R version //; s/ (.*//' | head -1)
        r-optparse: \$(Rscript -e "cat(as.character(packageVersion('optparse')))")
    END_VERSIONS
    """

    stub:
    """
    touch BambuOutput_counts_transcript_filter.txt
    touch BambuOutput_counts_gene_filter.txt
    touch BambuOutput_CPM_transcript_filter.txt
    touch BambuOutput_fullLengthCounts_transcript_filter.txt
    touch BambuOutput_uniqueCounts_transcript_filter.txt
    touch bambu_novel_transcripts_prefiltered.gtf
    touch prefilter_kept_transcripts.txt
    touch filter_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | sed 's/R version //; s/ (.*//' | head -1)
        r-optparse: \$(Rscript -e "cat(as.character(packageVersion('optparse')))")
    END_VERSIONS
    """
}
