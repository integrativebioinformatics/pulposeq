process NOVEL_TRANSCRIPTS {
    label 'process_novel_transcripts'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://itsiaguara/pulposeq:test':
        'docker.io/itsiaguara/pulposeq:test' }"

    input:
    path bambu_gtf
    path compared_gtf
    path tmap_file
    path coding_predictions
    path tx_counts
    path gene_counts
    path se_rds
    path annotation
    path r_script
    path gtf_utils

    output:
    path "novel_transcripts_metadata.csv"          , emit: novel_transcripts_metadata
    path "novel_lncRNAs_metadata.csv"              , emit: novel_lncrnas_metadata
    path "novel_protein-coding_metadata.csv"       , emit: novel_mrnas_metadata
    path "novel_non_coding_metadata.csv"           , emit: novel_non_coding_metadata
    path "novel_transcripts_validated_metadata.csv"          , emit: novel_combined_metadata
    path "novel_lncRNAs.gtf"                       , emit: novel_lncrnas_gtf
    path "novel_protein-coding.gtf"                , emit: novel_mrnas_gtf
    path "novel_non_coding.gtf"                    , emit: novel_non_coding_gtf
    path "novel_transcripts_validated.gtf"         , emit: novel_gtf
    path "novel_lncRNA_exon_lengths.csv"           , emit: novel_lncrna_exon_lengths
    path "novel_protein-coding_exon_lengths.csv"   , emit: novel_mrna_exon_lengths
    path "bambu_novel_tx_counts.csv"          , emit: novel_tx_counts
    path "bambu_novel_gene_counts.csv"        , emit: novel_gene_counts
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    Rscript $r_script \\
        --bambu_gtf ${bambu_gtf} \\
        --compared_gtf ${compared_gtf} \\
        --tmap_file ${tmap_file} \\
        --coding_predictions ${coding_predictions} \\
        --tx_counts ${tx_counts} \\
        --gene_counts ${gene_counts} \\
        --annotation ${annotation} \\
        --se_rds ${se_rds} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | sed 's/R version //; s/ (.*//' | head -1)
        r-readr: \$(Rscript -e "cat(as.character(packageVersion('readr')))")
        r-optparse: \$(Rscript -e "cat(as.character(packageVersion('optparse')))")
        r-dplyr: \$(Rscript -e "cat(as.character(packageVersion('dplyr')))")
        bioconductor-rtracklayer: \$(Rscript -e "cat(as.character(packageVersion('rtracklayer')))")
        bioconductor-genomicranges: \$(Rscript -e "cat(as.character(packageVersion('GenomicRanges')))")
        bioconductor-summarizedexperiment: \$(Rscript -e "cat(as.character(packageVersion('SummarizedExperiment')))")
    END_VERSIONS
    """

    stub:
    """
    touch novel_transcripts_metadata.csv
    touch novel_lncRNAs_metadata.csv
    touch novel_protein-coding_metadata.csv
    touch novel_non_coding_metadata.csv
    touch novel_transcripts_validated_metadata.csv
    touch novel_lncRNAs.gtf
    touch novel_protein-coding.gtf
    touch novel_non_coding.gtf
    touch novel_transcripts_validated.gtf
    touch novel_lncRNA_exon_lengths.csv
    touch novel_protein-coding_exon_lengths.csv
    touch bambu_novel_tx_counts.csv
    touch bambu_novel_gene_counts.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | sed 's/R version //; s/ (.*//' | head -1)
        r-readr: \$(Rscript -e "cat(as.character(packageVersion('readr')))")
        r-optparse: \$(Rscript -e "cat(as.character(packageVersion('optparse')))")
        r-dplyr: \$(Rscript -e "cat(as.character(packageVersion('dplyr')))")
        bioconductor-rtracklayer: \$(Rscript -e "cat(as.character(packageVersion('rtracklayer')))")
        bioconductor-genomicranges: \$(Rscript -e "cat(as.character(packageVersion('GenomicRanges')))")
        bioconductor-summarizedexperiment: \$(Rscript -e "cat(as.character(packageVersion('SummarizedExperiment')))")
    END_VERSIONS
    """
}
