process VALIDATE_BAMBU_OUTPUTS {
    label 'process_post_refinement'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://itsiaguara/pulposeq:test':
        'docker.io/itsiaguara/pulposeq:test' }"

    input:
    path se_rds
    path se_gene_rds
    path counts_transcript_validated
    path counts_gene_validated
    path r_script

    output:
    path "pca_validated.png"                , emit: pca
    path "pca_grouped_validated.png"        , emit: pca_grouped
    path "heatmap_gene_validated.png"       , emit: h_gene
    path "heatmap_transcript_validated.png" , emit: h_transcript
    path "se_multiSample_validated.rds"     , emit: rds_transcript
    path "seGene_multiSample_validated.rds" , emit: rds_gene
    path "versions.yml"                     , emit: versions
    path "validation_summary.csv"           , emit: validation_summary

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    Rscript $r_script \\
        --se_rds ${se_rds} \\
        --se_gene_rds ${se_gene_rds} \\
        --counts_transcript ${counts_transcript_validated} \\
        --counts_gene ${counts_gene_validated} \\
        --outdir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | sed 's/R version //; s/ (.*//' | head -1)
        r-optparse: \$(Rscript -e "cat(as.character(packageVersion('optparse')))")
        r-ggplot2: \$(Rscript -e "cat(as.character(packageVersion('ggplot2')))")
        bioconductor-bambu: \$(Rscript -e "cat(as.character(packageVersion('bambu')))")
        bioconductor-summarizedexperiment: \$(Rscript -e "cat(as.character(packageVersion('SummarizedExperiment')))")
    END_VERSIONS
    """

    stub:
    """
    touch pca_validated.png
    touch pca_grouped_validated.png
    touch heatmap_gene_validated.png
    touch heatmap_transcript_validated.png
    touch se_multiSample_validated.rds
    touch seGene_multiSample_validated.rds
    printf 'feature_type,curated_set,total_quantified,extended_annotations,recovery\\ntranscript,1,2,4,0.5\\ngene,1,2,4,0.5\\n' > validation_summary.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | sed 's/R version //; s/ (.*//' | head -1)
        r-optparse: \$(Rscript -e "cat(as.character(packageVersion('optparse')))")
        r-ggplot2: \$(Rscript -e "cat(as.character(packageVersion('ggplot2')))")
        bioconductor-bambu: \$(Rscript -e "cat(as.character(packageVersion('bambu')))")
        bioconductor-summarizedexperiment: \$(Rscript -e "cat(as.character(packageVersion('SummarizedExperiment')))")
    END_VERSIONS
    """
}
