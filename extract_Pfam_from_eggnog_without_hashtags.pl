#!/bin/bash

# Directory containing the files
directory="path/curated"

# Create the processed directory if it doesn't exist
mkdir -p pfam

# Loop over all files in the directory
for file in "$directory"/*
do
  # Extract the base name of the file (without the path)
  basename=$(basename "$file")
  
  # Extract the 21st column and add the header
  awk -F'\t' 'BEGIN {print "Pfam"} {print $21}' "$file" > "$directory/pfam_$basename"
  
  # Move the processed file to the "processed" directory
  mv "$directory/pfam_$basename" "pfam/"
done

