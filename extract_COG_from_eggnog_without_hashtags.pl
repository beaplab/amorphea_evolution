#!/bin/bash

# Directory containing the files
directory="path/curated"

# Create the processed directory if it doesn't exist
mkdir -p COG

# Loop over all files in the directory
for file in "$directory"/*
do
  # Extract the base name of the file (without the path)
  basename=$(basename "$file")
  
  # Extract the 7th column and add the header
  awk -F'\t' 'BEGIN {print "COG_category"} {print $7}' "$file" > "$directory/COG_$basename"
  
  # Move the processed file to the "processed" directory
  mv "$directory/COG_$basename" "COG/"
done


