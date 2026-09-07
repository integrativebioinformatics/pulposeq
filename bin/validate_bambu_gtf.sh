#!/usr/bin/env bash
# Validate GTF file based on validated transcript IDs from count files
#
# Usage: validate_bambu_gtf.sh --gtf <file> --awk_script <file> \
#          --counts <file> --full_length <file> --unique <file>

set -euo pipefail

# --- Argument Parsing ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --gtf) GTF_FILE="$2"; shift 2 ;;
        --awk_script) AWK_SCRIPT="$2"; shift 2 ;;
        --counts) COUNTS_FILE="$2"; shift 2 ;;
        --full_length) FULL_LENGTH_FILE="$2"; shift 2 ;;
        --unique) UNIQUE_FILE="$2"; shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

# --- Configuration ---
declare -a VALIDATED_FILES=(
    "$COUNTS_FILE"
    "$FULL_LENGTH_FILE"
    "$UNIQUE_FILE"
)

declare -a OUTPUT_GTF_FILES=(
    "BambuOutput_annotations_validated.gtf"
    "BambuOutput_fullLength_validated.gtf"
    "BambuOutput_uniquelyMapped_validated.gtf"
)

# --- Validation Functions ---

check_files() {
    if [ ! -f "$GTF_FILE" ] || [ ! -r "$GTF_FILE" ]; then
        echo "Error: GTF file '$GTF_FILE' not found or is not readable."
        exit 1
    fi

    local found_file=false
    for file in "${VALIDATED_FILES[@]}"; do
        if [ -f "$file" ] && [ -r "$file" ]; then
            found_file=true
            break
        fi
    done

    if [ "$found_file" = false ]; then
        echo "Error: No validated transcript files found."
        exit 1
    fi
}

extract_transcript_ids() {
    local input_file="$1"
    local temp_file="$2"

    echo "Extracting transcript IDs from $(basename "$input_file")..."

    txname_col=$(head -1 "$input_file" | awk -F'\t' '{for(i=1;i<=NF;i++) if($i=="TXNAME") print i}')
    if [ -z "$txname_col" ]; then
        echo "Error: TXNAME column not found in header of '$input_file'."
        rm -f "$temp_file"
        return 1
    fi

    awk -v col="$txname_col" 'NR > 1 { print $col }' "$input_file" | sort -u > "$temp_file"

    local count
    count=$(wc -l < "$temp_file")
    echo "  Found $count unique transcript IDs"
}

subset_gtf() {
    local transcript_ids_file="$1"
    local output_gtf="$2"

    echo "Creating GTF subset: $(basename "$output_gtf")..."

    awk -v ids_file="$transcript_ids_file" -f "$AWK_SCRIPT" "$GTF_FILE" > "$output_gtf"

    local line_count
    line_count=$(grep -v "^#" "$output_gtf" | wc -l)
    echo "  GTF subset created with $line_count feature lines"
}

# --- Main Logic ---

echo "Starting GTF subsetting process..."

check_files

processed_files=0
for i in "${!VALIDATED_FILES[@]}"; do
    input_file="${VALIDATED_FILES[$i]}"
    output_gtf="${OUTPUT_GTF_FILES[$i]}"

    if [ ! -f "$input_file" ] || [ ! -r "$input_file" ]; then
        echo "Warning: Validated file '$input_file' not found or is not readable. Skipping."
        continue
    fi

    temp_ids_file="$(mktemp -p . temp_transcript_ids.XXXXXX)"

    if ! extract_transcript_ids "$input_file" "$temp_ids_file"; then
        echo "Error: Failed to extract transcript IDs from '$input_file'." >&2
        rm -f "$temp_ids_file"
        exit 1
    fi

    if [ ! -s "$temp_ids_file" ]; then
        echo "Warning: No transcript IDs found in '$input_file'. Skipping GTF subset creation."
        rm -f "$temp_ids_file"
        continue
    fi

    subset_gtf "$temp_ids_file" "$output_gtf"

    if [ -f "$output_gtf" ] && [ -s "$output_gtf" ]; then
        echo "✓ GTF subset '$(basename "$output_gtf")' created successfully."
    else
        echo "✗ Error: Failed to create GTF subset '$(basename "$output_gtf")'."
        rm -f "$temp_ids_file"
        exit 1
    fi

    rm -f "$temp_ids_file"
    processed_files=$((processed_files + 1))
    echo ""
done

if [ "$processed_files" -eq 0 ]; then
    echo "Error: No validated files were processed. Exiting."
    exit 1
fi

echo "GTF subsetting process completed."
echo ""
echo "Generated GTF files:"
for output_gtf in "${OUTPUT_GTF_FILES[@]}"; do
    if [ -f "$output_gtf" ]; then
        echo "  - $(basename "$output_gtf")"
    fi
done