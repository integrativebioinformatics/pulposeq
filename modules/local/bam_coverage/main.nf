process BAM_COVERAGE {
    tag "$meta.id"
    // Not process_low. That label is for the small R metadata scripts, which hold
    // the annotation and peak around 4 GB. This reads alignments, and even chunked
    // per contig its peak tracks the largest chromosome's reads.
    label 'process_bam_coverage'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://itsiaguara/pulposeq:test':
        'docker.io/itsiaguara/pulposeq:test' }"

    input:
    // bam and bai arrive joined on meta, so a sample can never be paired with
    // another sample's index. Both are `path` so Nextflow stages them and binds
    // their real locations into the container.
    tuple val(meta), path(bam), path(bai)
    path r_script

    output:
    tuple val(meta), path("*.bw"), emit: bigwig
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    Rscript $r_script \\
        --bam ${bam} \\
        --prefix ${prefix} \\
        --outdir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | sed 's/R version //; s/ (.*//' | head -1)
        r-optparse: \$(Rscript -e "cat(as.character(packageVersion('optparse')))")
        bioconductor-rsamtools: \$(Rscript -e "cat(as.character(packageVersion('Rsamtools')))")
        bioconductor-genomicalignments: \$(Rscript -e "cat(as.character(packageVersion('GenomicAlignments')))")
        bioconductor-rtracklayer: \$(Rscript -e "cat(as.character(packageVersion('rtracklayer')))")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bw

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | sed 's/R version //; s/ (.*//' | head -1)
        r-optparse: \$(Rscript -e "cat(as.character(packageVersion('optparse')))")
        bioconductor-rsamtools: \$(Rscript -e "cat(as.character(packageVersion('Rsamtools')))")
        bioconductor-genomicalignments: \$(Rscript -e "cat(as.character(packageVersion('GenomicAlignments')))")
        bioconductor-rtracklayer: \$(Rscript -e "cat(as.character(packageVersion('rtracklayer')))")
    END_VERSIONS
    """
}
