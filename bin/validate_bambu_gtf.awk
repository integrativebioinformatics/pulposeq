BEGIN {
    while ((getline id < ids_file) > 0) { ids[id]=1 }
    close(ids_file)
}
{
    if ($0 ~ /^#/) next
    if ($0 !~ /transcript_id/) next

    line = $0
    sub(/.*transcript_id[ \t]*/, "", line)

    # Remove optional leading quote (double or single)
    sub(/^["']/, "", line)

    # transcript id ends at first quote, semicolon, or whitespace
    split(line, a, /["'; \t]/)
    tid = a[1]

    if (tid != "" && (tid in ids)) print $0
}