process GENOMIC_CONTEXT {
    label 'process_genomic_context'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://itsiaguara/pulposeq:test':
        'docker.io/itsiaguara/pulposeq:test' }"

    input:
    path gtf_file      // validated, enriched transcriptome GTF
    path bigwigs       // every sample's coverage track, collected
    path counts        // validated transcript counts, one column per sample
    path fl_counts     // validated full-length transcript counts
    path annotation    // reference annotation, for annotated gene and intron structure
    path r_script

    output:
    path "genomic_context_*.png"          , emit: figures, optional: true
    path "genomic_context_candidates.csv" , emit: candidates
    // Sense-intronic candidates, windowed on the whole host intron so the coverage
    // either side of the candidate is visible.
    path "intronic_context_*.png"          , emit: intronic_figures, optional: true
    path "intronic_context_candidates.csv" , emit: intronic_candidates
    // The plotted regions as a small GTF. Built to feed makeTxDbFromGFF, but useful
    // on its own: load it in IGV beside the published bigWigs to explore the same
    // loci interactively.
    path "genomic_context_regions.gtf"    , emit: regions_gtf, optional: true
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def bw_list   = bigwigs.collect { bw -> bw.name }.join(',')
    def name_list = bigwigs.collect { bw -> bw.name.replaceAll(/\.bw$/, '') }.join(',')
    """
    Rscript $r_script \\
        --gtf ${gtf_file} \\
        --bigwigs ${bw_list} \\
        --names ${name_list} \\
        --counts ${counts} \\
        --fl_counts ${fl_counts} \\
        --annotation ${annotation} \\
        --outdir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | sed 's/R version //; s/ (.*//' | head -1)
        r-optparse: \$(Rscript -e "cat(as.character(packageVersion('optparse')))")
        bioconductor-plotgardener: \$(Rscript -e "cat(as.character(packageVersion('plotgardener')))")
        bioconductor-rtracklayer: \$(Rscript -e "cat(as.character(packageVersion('rtracklayer')))")
        bioconductor-txdbmaker: \$(Rscript -e "cat(tryCatch(as.character(packageVersion('txdbmaker')), error = function(e) 'not installed'))")
    END_VERSIONS
    """

    stub:
    """
    # Columns and figure names track what genomic_context.R actually writes. The
    # report reads these tables by column name, so a stub header that drifts from
    # the real one makes -stub pass while a real run fails.
    touch genomic_context_STUBGENE_stub-stratum.png
    touch genomic_context_regions.gtf
    printf 'stratum,novel_biotype,classification,class_code,ref_class,ref_gene_biotype,BambuTxClass,example_tx,gene_id,gene_name,label,panel_genes,chrom,start,end,win_s,win_e,n_tx,n_known,n_novel,n_annotated,n_novel_lnc,biotypes,figure\\n' > genomic_context_candidates.csv
    printf 'stub-stratum,novel_lncRNA,sense intronic (i),i,lncRNA,protein_coding,newWithin,BambuTxSTUB,ENSG00000000000,STUBGENE,STUBGENE,ENSG00000000000,chr1,1000,9000,950,9050,4,3,1,0,1,protein_coding;novel_lncRNA,genomic_context_STUBGENE_stub-stratum.png\\n' >> genomic_context_candidates.csv

    touch intronic_context_STUBTX.png
    printf 'qry_id,ref_gene_id,ref_gene_name,ref_gene_biotype,chrom,start,end,win_s,win_e,strand,class_code,BambuTxClass,num_exons,samples_quantified,samples_total,quantified_samples,counts_total,full_length_support,figure\\n' > intronic_context_candidates.csv
    printf 'BambuTxSTUB,ENSG00000000000,STUBGENE,protein_coding,chr1,3000,3800,2000,5000,+,i,newWithin,1,2,3,s1;s2,42,FALSE,intronic_context_STUBTX.png\\n' >> intronic_context_candidates.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | sed 's/R version //; s/ (.*//' | head -1)
        r-optparse: \$(Rscript -e "cat(as.character(packageVersion('optparse')))")
        bioconductor-plotgardener: \$(Rscript -e "cat(as.character(packageVersion('plotgardener')))")
        bioconductor-rtracklayer: \$(Rscript -e "cat(as.character(packageVersion('rtracklayer')))")
        bioconductor-txdbmaker: \$(Rscript -e "cat(tryCatch(as.character(packageVersion('txdbmaker')), error = function(e) 'not installed'))")
    END_VERSIONS
    """
}
