#!/bin/bash

# Ensure the script is being run with bash
if [ -z "$BASH_VERSION" ]; then
  echo "This script requires bash to run."
  exit 1
fi

# Check if a file path is provided as an argument
if [ -z "$1" ]; then
  echo "Usage: $0 path_to_xml_file"
  exit 1
fi

# Extract the value of each taxon attribute
grep -o 'taxon="[^"]*"' "$1" | awk -F '"' '{print $2}' > taxon_values.txt

# Print all extracted taxon values
echo "Extracted taxon values:"
cat taxon_values.txt

# Find and print duplicates
echo
echo "Duplicate taxon values:"
sort taxon_values.txt | uniq -d -c | while read count taxon; do
  if [ "$count" -gt 1 ]; then
    echo "Taxon: $taxon (Count: $count)"
  fi
done

# Clean up
rm taxon_values.txt
