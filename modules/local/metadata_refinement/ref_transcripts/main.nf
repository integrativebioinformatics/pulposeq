process REF_TRANSCRIPTS {
    label 'process_low'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://itsiaguara/pulposeq:test':
        'docker.io/itsiaguara/pulposeq:test' }"

    input:
    path transcript_counts
    path gene_counts
    path gtf_file
    path annotation
    path r_script
    path gtf_utils

    output:
    path "annotated_transcriptome_metadata.csv"          , emit: transcriptome_metadata
    path "annotated_lncRNAs_metadata.csv"                , emit: lncrna_metadata
    path "annotated_lncRNAs_exonlength.csv"              , emit: lncrna_exonlength
    path "annotated_protein-coding_metadata.csv"         , emit: protein_coding_metadata
    path "annotated_protein-coding_exonlength.csv"       , emit: protein_coding_exonlength
    path "bambu_annotated_transcriptome.gtf"             , emit: annotated_transcriptome_gtf
    path "bambu_annotated_transcriptome_tx_counts.csv"   , emit: annotated_tx_counts
    path "bambu_annotated_lncRNAs.gtf"                   , emit: annotated_lncrna_gtf
    path "bambu_annotated_mRNAs.gtf"                     , emit: annotated_mrna_gtf
    path "bambu_annotated_transcriptome_gene_counts.csv" , emit: annotated_gene_counts
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    """
    # Run the R script directly using the input path variable
    Rscript $r_script \\
        --transcript_counts ${transcript_counts} \\
        --gene_counts ${gene_counts} \\
        --gtf_file ${gtf_file} \\
        --annotation ${annotation} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | sed 's/R version //; s/ (.*//' | head -1)
        r-readr: \$(Rscript -e "cat(as.character(packageVersion('readr')))")
        r-optparse: \$(Rscript -e "cat(as.character(packageVersion('optparse')))")
        r-dplyr: \$(Rscript -e "cat(as.character(packageVersion('dplyr')))")
        bioconductor-genomicranges: \$(Rscript -e "cat(as.character(packageVersion('GenomicRanges')))")
        bioconductor-rtracklayer: \$(Rscript -e "cat(as.character(packageVersion('rtracklayer')))")
    END_VERSIONS
    """

    stub:
    """
    touch annotated_transcriptome_metadata.csv
    touch annotated_lncRNAs_metadata.csv
    touch annotated_lncRNAs_exonlength.csv
    touch annotated_protein-coding_metadata.csv
    touch annotated_protein-coding_exonlength.csv
    touch bambu_annotated_transcriptome.gtf
    touch bambu_annotated_transcriptome_tx_counts.csv
    touch bambu_annotated_lncRNAs.gtf
    touch bambu_annotated_mRNAs.gtf
    touch bambu_annotated_transcriptome_gene_counts.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | sed 's/R version //; s/ (.*//' | head -1)
        r-readr: \$(Rscript -e "cat(as.character(packageVersion('readr')))")
        r-optparse: \$(Rscript -e "cat(as.character(packageVersion('optparse')))")
        r-dplyr: \$(Rscript -e "cat(as.character(packageVersion('dplyr')))")
        bioconductor-genomicranges: \$(Rscript -e "cat(as.character(packageVersion('GenomicRanges')))")
        bioconductor-rtracklayer: \$(Rscript -e "cat(as.character(packageVersion('rtracklayer')))")
    END_VERSIONS
    """
}
