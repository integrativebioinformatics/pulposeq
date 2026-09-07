#!/usr/bin/env bash
# Validate Bambu count files against novel transcripts metadata.
# Keeps only genes transcribing isoforms assembled by compatible RCs
#
# Usage: validate_bambu_counts.sh --metadata <file> --counts_gene <file> \
#          --counts_transcript <file> --cpm_transcript <file> \
#          --full_length <file> --unique <file>

set -euo pipefail

# --- Argument Parsing ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --metadata) METADATA_FILE="$2"; shift 2 ;;
        --counts_gene) COUNTS_GENE="$2"; shift 2 ;;
        --counts_transcript) COUNTS_TRANSCRIPT="$2"; shift 2 ;;
        --cpm_transcript) CPM_TRANSCRIPT="$2"; shift 2 ;;
        --full_length) FULL_LENGTH="$2"; shift 2 ;;
        --unique) UNIQUE="$2"; shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

# --- Configuration ---
SCRIPT_DIR="$(pwd)"
TEMP_TXID_FILE="${SCRIPT_DIR}/temp_valid_tx_ids.txt"
TEMP_GENEID_FILE="${SCRIPT_DIR}/temp_valid_gene_ids.txt"

declare -a TX_INPUT_FILES=(
    "$COUNTS_TRANSCRIPT"
    "$CPM_TRANSCRIPT"
    "$FULL_LENGTH"
    "$UNIQUE"
)

declare -a GENE_INPUT_FILES=(
    "$COUNTS_GENE"
)

# --- Validation Functions ---

check_metadata_file() {
    if [ ! -f "$METADATA_FILE" ] || [ ! -r "$METADATA_FILE" ]; then
        echo "Error: Metadata file '$METADATA_FILE' not found or is not readable."
        exit 1
    fi
}

extract_valid_tx_ids() {
    # Extract qry_id (column 2) from CSV, skip header
    awk -F',' 'NR > 1 {
        gsub(/"/, "", $2);
        print $2;
    }' "$METADATA_FILE" | sort -u > "$TEMP_TXID_FILE"
}

# --- Processing Functions ---

process_transcript_file() {
    local input_file="$1"
    local output_file="$2"

    echo "Processing transcript file: $(basename "$input_file")"

    awk -v tx_temp_file="$TEMP_TXID_FILE" '
    BEGIN {
        while ((getline line < tx_temp_file) > 0) {
            valid_ids[line] = 1
        }
        close(tx_temp_file)
    }
    NR == 1 {
        print $0
        next
    }
    {
        if ($1 ~ /^BambuTx/) {
            if ($1 in valid_ids) {
                print $0
            }
        } else {
            print $0
        }
    }' "$input_file" > "$output_file"
}

extract_valid_gene_ids() {
    local tx_validated_file="$1"

    echo "Extracting unique gene IDs from validated transcript file..."

    # Extract column 2 (gene IDs), skip header, sort and get unique values
    awk 'NR > 1 { print $2 }' "$tx_validated_file" | sort -u > "$TEMP_GENEID_FILE"

    if [ ! -s "$TEMP_GENEID_FILE" ]; then
        echo "Warning: No gene IDs were extracted from '$(basename "$tx_validated_file")'. Gene file will not be processed."
        return 1
    fi

    echo "✓ Extracted $(wc -l < "$TEMP_GENEID_FILE") unique gene IDs"
}

process_gene_file() {
    local input_file="$1"
    local output_file="$2"

    echo "Processing gene file: $(basename "$input_file")"

    awk -v gene_temp_file="$TEMP_GENEID_FILE" '
    BEGIN {
        while ((getline line < gene_temp_file) > 0) {
            valid_gene_ids[line] = 1
        }
        close(gene_temp_file)
    }
    NR == 1 {
        print $0
        next
    }
    {
        if ($1 in valid_gene_ids) {
            print $0
        }
    }' "$input_file" > "$output_file"
}

# --- Main Logic ---

echo "Starting validation process..."

check_metadata_file

echo "Extracting valid transcript IDs from metadata..."
extract_valid_tx_ids

if [ ! -s "$TEMP_TXID_FILE" ]; then
    echo "Error: No valid transcript IDs were extracted from metadata. '$TEMP_TXID_FILE' is empty."
    exit 1
fi

# Process each transcript input file
for input_file in "${TX_INPUT_FILES[@]}"; do
    if [ ! -f "$input_file" ] || [ ! -r "$input_file" ]; then
        echo "Warning: Input file '$input_file' not found or is not readable. Skipping."
        continue
    fi

    # Generate output filename
    base_name=$(basename "$input_file" _filter.txt)
    output_file="${SCRIPT_DIR}/${base_name}_validated.txt"

    # Prevent overwriting the input file
    if [ "$(realpath "$input_file")" == "$(realpath "$output_file")" ]; then
        echo "Error: Input and output filenames would be the same for '$input_file'. Skipping."
        continue
    fi

    process_transcript_file "$input_file" "$output_file"

    if [ -s "$output_file" ]; then
        echo "✓ Validated file '$(basename "$output_file")' created successfully."
    else
        if [ -f "$output_file" ]; then
            echo "⚠ Validated file '$(basename "$output_file")' created, but it might be empty or only contain the header."
        else
            echo "✗ Error: Failed to create validated file '$(basename "$output_file")'."
        fi
    fi

    echo ""
done

# Extract valid gene IDs from the validated transcript counts file
# and then process the gene counts file
tx_validated_file="${SCRIPT_DIR}/BambuOutput_counts_transcript_validated.txt"

if [ -f "$tx_validated_file" ] && [ -s "$tx_validated_file" ]; then
    if extract_valid_gene_ids "$tx_validated_file"; then
        for input_file in "${GENE_INPUT_FILES[@]}"; do
            if [ ! -f "$input_file" ] || [ ! -r "$input_file" ]; then
                echo "Warning: Gene input file '$input_file' not found or is not readable. Skipping."
                continue
            fi

            base_name=$(basename "$input_file" _filter.txt)
            output_file="${SCRIPT_DIR}/${base_name}_validated.txt"

            if [ "$(realpath "$input_file")" == "$(realpath "$output_file")" ]; then
                echo "Error: Input and output filenames would be the same for '$input_file'. Skipping."
                continue
            fi

            process_gene_file "$input_file" "$output_file"

            if [ -s "$output_file" ]; then
                echo "✓ Validated file '$(basename "$output_file")' created successfully."
            else
                if [ -f "$output_file" ]; then
                    echo "⚠ Validated file '$(basename "$output_file")' created, but it might be empty or only contain the header."
                else
                    echo "✗ Error: Failed to create validated file '$(basename "$output_file")'."
                fi
            fi

            echo ""
        done
    else
        echo "⚠ Skipping gene file processing due to empty gene ID extraction."
    fi
else
    echo "Error: Validated transcript file not found. Cannot extract gene IDs."
    exit 1
fi

# Clean up temporary files
rm -f "$TEMP_TXID_FILE" "$TEMP_GENEID_FILE"

echo "Validation and filtering process completed."
echo ""
echo "Generated files:"
for input_file in "${TX_INPUT_FILES[@]}"; do
    base_name=$(basename "$input_file" _filter.txt)
    output_file="${SCRIPT_DIR}/${base_name}_validated.txt"
    if [ -f "$output_file" ]; then
        echo "  - $(basename "$output_file")"
    fi
done
for input_file in "${GENE_INPUT_FILES[@]}"; do
    base_name=$(basename "$input_file" _filter.txt)
    output_file="${SCRIPT_DIR}/${base_name}_validated.txt"
    if [ -f "$output_file" ]; then
        echo "  - $(basename "$output_file")"
    fi
done