#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_file> <output_file>"
    exit 1
fi

input_file="$1"
output_file="$2"

# Initialize variables
node_counter=0
current_tree=""

# Function to process each tree
process_tree() {
    local tree_content="$1"
    local node_counter=0

    echo "$tree_content" | awk -v node_counter="$node_counter" '
    {
        while (match($0, /node0/)) {
            sub(/node0/, "node" node_counter)
            node_counter++
        }
        print
    }'
}

# Clear the output file
> "$output_file"

# Read the input file and process each tree
while IFS= read -r line || [ -n "$line" ]; do
    if [[ $line =~ ^tree\ STATE_ ]]; then
        # If we encounter a new tree, process the previous tree
        if [ -n "$current_tree" ]; then
            processed_tree=$(process_tree "$current_tree")
            echo "$processed_tree" >> "$output_file"
            current_tree=""
        fi
    fi
    current_tree+="$line"$'\n'
done < "$input_file"

# Process the last tree
if [ -n "$current_tree" ]; then
    processed_tree=$(process_tree "$current_tree")
    echo "$processed_tree" >> "$output_file"
fi

echo "Processing complete. Output written to $output_file."
