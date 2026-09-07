#!/usr/bin/env bash
# Filter out zero-count rows from Bambu count/CPM files.
# Removes any gene or transcript where all sample columns are 0.
#
# Usage: filter_bambu_counts.sh --counts_gene <file> --counts_transcript <file> \
#          --cpm_transcript <file> --full_length <file> --unique <file>

set -euo pipefail

# --- Argument Parsing ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --counts_gene) COUNTS_GENE="$2"; shift 2 ;;
        --counts_transcript) COUNTS_TRANSCRIPT="$2"; shift 2 ;;
        --cpm_transcript) CPM_TRANSCRIPT="$2"; shift 2 ;;
        --full_length) FULL_LENGTH="$2"; shift 2 ;;
        --unique) UNIQUE="$2"; shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

# --- Configuration ---
declare -a INPUT_FILES=(
    "$COUNTS_GENE"
    "$COUNTS_TRANSCRIPT"
    "$CPM_TRANSCRIPT"
    "$FULL_LENGTH"
    "$UNIQUE"
)

# --- Processing Functions ---

# Filter gene count file (checks cols 2-NF)
process_gene_file() {
    local input_file="$1"
    local output_file="$2"

    echo "Processing gene file: $(basename "$input_file")"

    awk '
    NR == 1 {
        print $0
        next
    }
    {
        has_counts = 0
        for (i = 2; i <= NF; i++) {
            if ($i != 0) {
                has_counts = 1
                break
            }
        }
        if (has_counts) {
            print $0
        }
    }' "$input_file" > "$output_file"
}

# Filter transcript files (checks cols 3-NF)
process_transcript_file() {
    local input_file="$1"
    local output_file="$2"

    echo "Processing transcript file: $(basename "$input_file")"

    awk '
    NR == 1 {
        print $0
        next
    }
    {
        has_counts = 0
        for (i = 3; i <= NF; i++) {
            if ($i != 0) {
                has_counts = 1
                break
            }
        }
        if (has_counts) {
            print $0
        }
    }' "$input_file" > "$output_file"
}

# --- Main Logic ---

echo "Starting zero-count filtering process..."

processed_files=0
for input_file in "${INPUT_FILES[@]}"; do
    if [ ! -f "$input_file" ] || [ ! -r "$input_file" ]; then
        echo "Warning: Input file '$input_file' not found or is not readable. Skipping."
        continue
    fi

    # Generate output filename
    base_name=$(basename "$input_file" .txt)
    output_file="${base_name}_filter.txt"

    # Determine file type and process accordingly
    if [[ "$input_file" == *"counts_gene"* ]]; then
        process_gene_file "$input_file" "$output_file"
    else
        process_transcript_file "$input_file" "$output_file"
    fi

    if [ -s "$output_file" ]; then
        echo "✓ filter file '$(basename "$output_file")' created successfully."
    else
        if [ -f "$output_file" ]; then
            echo "⚠ filter file '$(basename "$output_file")' created, but it might be empty or only contain the header."
        else
            echo "✗ Error: Failed to create filter file '$(basename "$output_file")'."
            exit 1
        fi
    fi

    processed_files=$((processed_files + 1))
    echo ""
done

if [ "$processed_files" -eq 0 ]; then
    echo "Error: No input files were processed. Exiting."
    exit 1
fi

echo "Zero-count filtering process completed."
echo ""
echo "Generated files:"
for input_file in "${INPUT_FILES[@]}"; do
    base_name=$(basename "$input_file" .txt)
    output_file="${base_name}_filter.txt"
    if [ -f "$output_file" ]; then
        echo "  - $(basename "$output_file")"
    fi
done