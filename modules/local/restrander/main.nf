process RESTRANDER {
    tag "$meta.id"
    label 'process_restrander'

    // Restrander is in neither bioconda nor biocontainers, so there is no upstream
    // image to fall back on. Pinned to a tag rather than latest: the preset configs
    // ship inside the image and carry the primer sequences behind every orientation
    // call, so a silent change to them would change results.
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://emiyoshi/restrander:v1.1.3':
        'docker.io/emiyoshi/restrander:v1.1.3' }"

    input:
    tuple val(meta), path(reads)
    // Optional custom config. Pass [] to use the preset named by params.restrand_kit
    // from inside the container instead.
    path custom_config

    output:
    tuple val(meta), path("*.restranded.fastq.gz")          , emit: reads
    tuple val(meta), path("*.restrander_stats.json")        , emit: stats
    tuple val(meta), path("*.restranded-unknowns.fastq.gz") , emit: unknown, optional: true
    path "*.restrander_config.json"                         , emit: config_used
    path "versions.yml"                                     , emit: versions
    // The unknowns file is renamed in the script rather than matched by a glob.
    // Restrander inserts "-unknowns" before the FIRST dot of the output name, not
    // before the extension: out.fastq.gz becomes out-unknowns.fastq.gz, and
    // a.b.fastq.gz becomes a-unknowns.b.fastq.gz. Measured on v1.1.3.

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def kit      = params.restrand_kit
    def supplied = custom_config ? custom_config.name : ''
    // Restrander has no --version flag and prints no version anywhere -- running it
    // bare gives only a usage message. Taken from the pinned container tag instead;
    // keep the two in step.
    def version  = '1.1.3'
    """
    cfg="${supplied}"

    if [ -z "\$cfg" ]; then
        # The image publishes no documentation of its layout, so the preset is
        # located rather than assumed. Getting this wrong would not raise an error --
        # it would orient almost nothing and fall back to unstranded processing --
        # so it is worth failing loudly here instead.
        for d in /home/restrander/config /opt/restrander/config /usr/local/share/restrander/config /config; do
            if [ -f "\$d/${kit}.json" ]; then cfg="\$d/${kit}.json"; break; fi
        done
    fi

    if [ -z "\$cfg" ]; then
        cfg=\$(find / -type f -name "${kit}.json" 2>/dev/null | head -n 1 || true)
    fi

    if [ -z "\$cfg" ] || [ ! -f "\$cfg" ]; then
        echo "RESTRANDER: no config found for kit '${kit}' inside the container." >&2
        echo "Looked in /home/restrander/config (where emiyoshi/restrander:v1.1.3 keeps them)," >&2
        echo "/opt/restrander/config, /usr/local/share/restrander/config and /config," >&2
        echo "then searched the whole filesystem." >&2
        echo "Supply one from outside the image with --restrand_config /path/to/config.json." >&2
        exit 1
    fi

    # Copied into the task directory and published: the config carries the primer
    # sequences that determined every orientation call, so a run is not
    # reconstructable without it.
    cp "\$cfg" ${prefix}.restrander_config.json

    # Restrander takes three positional arguments and prints its statistics as JSON
    # on stdout, so the report is a redirect rather than a named output. Progress
    # messages go to stderr, so they do not contaminate it.
    #
    # Only reads whose orientation could not be determined are diverted to the
    # unknowns file. Artefactual reads -- TSO-TSO, RTP-RTP -- are counted in the
    # report but stay in the main output, so they are flagged rather than discarded
    # and read counts are not silently perturbed. Measured on v1.1.3: 20,000 reads
    # in, 210 unorientable diverted, 19,790 out, with 176 artefacts among them.
    restrander \\
        ${reads} \\
        ${prefix}.restranded.fastq.gz \\
        ${prefix}.restrander_config.json \\
        $args \\
        > ${prefix}.restrander_stats.json

    # Renamed to a predictable name. Restrander builds this filename by inserting
    # "-unknowns" before the first dot of the output, so the result depends on what
    # the sample id contains and cannot be matched by a fixed output glob. The loop
    # is a no-op when a sample has no unorientable reads.
    for u in *-unknowns*; do
        [ -e "\$u" ] || continue
        mv "\$u" ${prefix}.restranded-unknowns.fastq.gz
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        restrander: ${version}
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def version = '1.1.3'
    """
    echo '{}' > ${prefix}.restrander_config.json
    echo "" | gzip > ${prefix}.restranded.fastq.gz
    echo "" | gzip > ${prefix}.restranded-unknowns.fastq.gz

    # Shaped like the real report, because the workflow parses this file rather than
    # just passing it along: orientedFraction() in restranding.nf reads it to decide whether
    # the run may proceed at all, so an empty stub would break a stub run outright.
    #
    # These are the actual counts from restrander v1.1.3 on the 20,000-read PCB109
    # sample in the tool's own vignette, so the numbers a stub run reasons about are
    # the ones real data produces. (9526 + 10264) / 20000 = 0.9895, comfortably over
    # the 0.80 default, so a stub exercises the success path. Run with
    # --restrand_min_frac 0.99 to drive the fallback instead.
    cat <<-END_STATS > ${prefix}.restrander_stats.json
    {
        "stats": {
            "artefactStats": { "RTP-RTP": 47, "TSO-TSO": 129, "no artefact": 19824 },
            "strandStats": { "+": 9526, "-": 10264, "?": 210 },
            "totalReads": 20000
        }
    }
    END_STATS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        restrander: ${version}
    END_VERSIONS
    """
}
