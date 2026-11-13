process OARFISH_INDEX_FASTA {
    tag "transcriptome"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.7--hdcf5f25_4' :
        'quay.io/biocontainers/gffread:0.12.7--hdcf5f25_4' }"

    input:
    path gtf        // BambuOutput_extended_annotations.gtf: every reference transcript + every novel model
    path genome     // reference genome FASTA

    output:
    path "oarfish_annotated.fasta"   , emit: annotated
    path "oarfish_novel.fasta"       , emit: novel
    path "oarfish_index_summary.tsv" , emit: summary
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Split Bambu's extended annotation into reference transcripts and novel
    // (BambuTx) models. Unstranded novel models are left out: oarfish runs with a
    // forward-strand filter, so a model without a strand cannot be quantified
    // reliably, and PREFILTER removes them anyway.
    """
    awk -F '\\t' '
        /^#/ { next }
        {
            id = ""
            if (match(\$9, /transcript_id "[^"]+"/)) {
                id = substr(\$9, RSTART + 15, RLENGTH - 16)
            }
            if (index(id, "BambuTx") == 1) {
                if (\$7 != "+" && \$7 != "-") { unstranded[id] = 1; next }
                novel[id] = 1
                print > "novel.gtf"
            } else {
                if (id != "") { reference[id] = 1 }
                print > "annotated.gtf"
            }
        }
        END {
            n_ref = 0; n_novel = 0; n_unstr = 0
            for (k in reference)  n_ref++
            for (k in novel)      n_novel++
            for (k in unstranded) n_unstr++
            print "category\\ttranscripts"                       > "oarfish_index_summary.tsv"
            print "reference\\t" n_ref                          > "oarfish_index_summary.tsv"
            print "novel\\t" n_novel                            > "oarfish_index_summary.tsv"
            print "novel_unstranded_left_out\\t" n_unstr        > "oarfish_index_summary.tsv"
        }
    ' ${gtf}

    touch annotated.gtf novel.gtf

    gffread -w oarfish_annotated.fasta -g ${genome} annotated.gtf

    # No novel models is possible (e.g. a very strict --ndr). oarfish then runs on
    # the reference alone.
    if [ -s novel.gtf ]; then
        gffread -w oarfish_novel.fasta -g ${genome} novel.gtf
    else
        touch oarfish_novel.fasta
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1)
    END_VERSIONS
    """

    stub:
    """
    touch oarfish_annotated.fasta oarfish_novel.fasta
    printf 'category\\ttranscripts\\nreference\\t0\\nnovel\\t0\\nnovel_unstranded_left_out\\t0\\n' > oarfish_index_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1)
    END_VERSIONS
    """
}
