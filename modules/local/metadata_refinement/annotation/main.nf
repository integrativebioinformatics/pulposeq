process ANNOTATION {
    label 'process_low'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://itsiaguara/pulposeq:test':
        'docker.io/itsiaguara/pulposeq:test' }"

    input:
    path annotations_gtf
    path known_metadata
    path novel_metadata
    path annotation
    path r_script
    path gtf_utils

    output:
    // Known transcripts re-read from the reference with CDS and UTR intact, plus
    // the novel models. This is what the genomic context figures are drawn from.
    path "annotations_final.gtf"                    , emit: final_gtf
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    Rscript $r_script \\
        --annotations_gtf ${annotations_gtf} \\
        --known_metadata ${known_metadata} \\
        --novel_metadata ${novel_metadata} \\
        --annotation ${annotation} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | sed 's/R version //; s/ (.*//' | head -1)
        r-optparse: \$(Rscript -e "cat(as.character(packageVersion('optparse')))")
        bioconductor-genomicranges: \$(Rscript -e "cat(as.character(packageVersion('GenomicRanges')))")
        bioconductor-rtracklayer: \$(Rscript -e "cat(as.character(packageVersion('rtracklayer')))")
    END_VERSIONS
    """

    stub:
    """
    touch annotations_final.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | sed 's/R version //; s/ (.*//' | head -1)
        r-optparse: \$(Rscript -e "cat(as.character(packageVersion('optparse')))")
        bioconductor-genomicranges: \$(Rscript -e "cat(as.character(packageVersion('GenomicRanges')))")
        bioconductor-rtracklayer: \$(Rscript -e "cat(as.character(packageVersion('rtracklayer')))")
    END_VERSIONS
    """
}
