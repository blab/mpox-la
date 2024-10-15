#!/bin/bash

# Prompt the user to enter the path to the XML file
read -p "Enter the path to the XML file: " xml_file

# Check if the file exists
if [ ! -f "$xml_file" ]; then
    echo "File not found!"
    exit 1
fi

# Read each line of the XML file
while IFS= read -r line; do
    # Extract taxon ID from the line
    taxon_id=$(echo "$line" | sed -n 's/.*id="\([^"]*\)".*/\1/p')

    # Check if the line contains a taxon ID
    if [ ! -z "$taxon" ]; then
        # Extract the part before the second underscore
        part_before_second_underscore=$(echo "$taxon" | awk -F'_' '{print $1"_"$2}')

        # Print the result
        echo "Taxon ID: $taxon_id"
    fi
done < "$xml_file"
