process RENDER_REPORT {
    label 'process_reports'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://itsiaguara/pulposeq:test':
        'docker.io/itsiaguara/pulposeq:test' }"

    input:
    // bambu outputs
    path counts_genes
    path counts_transcript
    path fl_counts_transcript
    path u_counts_transcript

    // metadata csv tables
    path transcriptome_meta
    path protein_coding_meta
    path lncrna_meta
    path protein_coding_exonlength
    path lncrna_exonlength
    path novel_transcriptome_meta
    path novel_lncrna_exonlength
    path novel_protein_coding_exonlength
    path novel_non_coding_meta

    // genomic context figures. The PNGs are staged rather than passed as
    // parameters: the candidate tables carry the filenames, and the report reads
    // them from the working directory.
    path genomic_context_candidates
    path genomic_context_figures
    path intronic_context_candidates
    path intronic_context_figures

    // raw bambu plots (pre-refinement assembly)
    path raw_pca
    path raw_pca_grouped
    path raw_heatmap_gene
    path raw_heatmap_transcript

    // regenerated plots from the validated transcriptome
    path post_pca
    path post_pca_grouped
    path post_heatmap_gene
    path post_heatmap_transcript

    // run metrics
    path bambu_metrics
    path validation_summary

    // quarto template
    path qmd_report

    output:
    path "*.html"                      , emit: report
    path "Bambu_assembly_summary.csv"  , emit: bambu_assembly_summary
    path "Bambu_lncRNA_PC_summary.csv" , emit: bambu_lncRNA_PC_summary
    path "versions.yml"                , emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export XDG_CACHE_HOME=/tmp/quarto_cache_home
    export XDG_DATA_HOME=/tmp/quarto_data_home

    quarto render $qmd_report \\
        -P counts_genes:${counts_genes} \\
        -P counts_transcript:${counts_transcript} \\
        -P fl_counts_transcript:${fl_counts_transcript} \\
        -P u_counts_transcript:${u_counts_transcript} \\
        -P transcriptome_meta:${transcriptome_meta} \\
        -P protein_coding_meta:${protein_coding_meta} \\
        -P lncrna_meta:${lncrna_meta} \\
        -P protein_coding_exonlength:${protein_coding_exonlength} \\
        -P lncrna_exonlength:${lncrna_exonlength} \\
        -P novel_transcriptome_meta:${novel_transcriptome_meta} \\
        -P novel_lncrna_exonlength:${novel_lncrna_exonlength} \\
        -P novel_protein_coding_exonlength:${novel_protein_coding_exonlength} \\
        -P novel_non_coding_meta:${novel_non_coding_meta} \\
        -P genomic_context_candidates:${genomic_context_candidates} \\
        -P intronic_context_candidates:${intronic_context_candidates} \\
        -P raw_pca:${raw_pca} \\
        -P raw_pca_grouped:${raw_pca_grouped} \\
        -P raw_heatmap_gene:${raw_heatmap_gene} \\
        -P raw_heatmap_transcript:${raw_heatmap_transcript} \\
        -P post_pca:${post_pca} \\
        -P post_pca_grouped:${post_pca_grouped} \\
        -P post_heatmap_gene:${post_heatmap_gene} \\
        -P post_heatmap_transcript:${post_heatmap_transcript} \\
        -P bambu_metrics:${bambu_metrics} \\
        -P validation_summary:${validation_summary} \\
        --to html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: \$(quarto --version | head -n1 | sed 's/Quarto //')
        r-base: \$(R --version | head -n1 | sed 's/R version //; s/ .*//')
        r-tidyverse: \$(Rscript -e "cat(as.character(packageVersion('tidyverse')))")
        r-cowplot: \$(Rscript -e "cat(as.character(packageVersion('cowplot')))")
        r-scales: \$(Rscript -e "cat(as.character(packageVersion('scales')))")
        r-RColorBrewer: \$(Rscript -e "cat(as.character(packageVersion('RColorBrewer')))")
        r-viridis: \$(Rscript -e "cat(as.character(packageVersion('viridis')))")
    END_VERSIONS
    """
    
    stub:
    """
    touch report.html
    touch Bambu_assembly_summary.csv
    touch Bambu_lncRNA_PC_summary.csv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: \$(quarto --version | head -n1 | sed 's/Quarto //')
        r-base: \$(R --version | head -n1 | sed 's/R version //; s/ .*//')
        r-tidyverse: \$(Rscript -e "cat(as.character(packageVersion('tidyverse')))")
        r-cowplot: \$(Rscript -e "cat(as.character(packageVersion('cowplot')))")
        r-scales: \$(Rscript -e "cat(as.character(packageVersion('scales')))")
        r-RColorBrewer: \$(Rscript -e "cat(as.character(packageVersion('RColorBrewer')))")
        r-viridis: \$(Rscript -e "cat(as.character(packageVersion('viridis')))")
    END_VERSIONS
    """
}
