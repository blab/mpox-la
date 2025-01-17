#!/bin/bash

# Check if the input file is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input.tsv>"
    exit 1
fi

INPUT_FILE="$1"
OUTPUT_FILE="updated_metadata.tsv"

# Extract the header and determine the position of the "focus_areas" column
HEADER=$(head -n 1 "$INPUT_FILE")
FOCUS_COLUMN_INDEX=$(echo "$HEADER" | awk -F'\t' '{for(i=1;i<=NF;i++) if($i=="focus_areas") print i}')

# Check if "focus_areas" was found
if [ -z "$FOCUS_COLUMN_INDEX" ]; then
    echo "\"focus_areas\" column not found in the input file."
    exit 1
fi

# Use awk to rearrange the columns
awk -v focus_index="$FOCUS_COLUMN_INDEX" -F'\t' '
BEGIN { OFS = FS }
{
    # Print columns 1 to 3
    for (i = 1; i <= 3; i++) {
        printf "%s", $i
        if (i < 3) printf OFS
    }

    # Print the "focus_areas" column as the 6th column
    printf OFS $focus_index

    # Print the remaining columns, skipping the original "focus_areas" position
    for (i = 4; i <= NF; i++) {
        if (i != focus_index) {
            printf OFS $i
        }
    }

    # End the line
    printf "\n"
}' "$INPUT_FILE" > "$OUTPUT_FILE"

echo "The 'focus_areas' column has been moved to the 6th position in $OUTPUT_FILE."
