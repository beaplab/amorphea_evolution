library(tidyverse)
library(readr)
library(dplyr)
library(tidyr)

# Define the input and output directories (remove hashtags from eggnog annotations and extract_COG_from_eggnog_without_hashtags.pl beforehand)
input_directory <- "/home/agalvez/data/metabolic_amoeba/version4/data/amoebozoa_check/curated/COG"
output_directory <- "/home/agalvez/data/metabolic_amoeba/version4/data/amoebozoa_check/curated/processed"

# Create output directory if it doesn't exist
if (!dir.exists(output_directory)) {
  dir.create(output_directory)
}

# List all files in the input directory
files <- list.files(input_directory, pattern = "\\.fasta\\.annotations$", full.names = TRUE)

# Function to process a single file
process_file <- function(file) {
  # Extract file name and create output file path
  filename <- basename(file)
  output_file <- file.path(output_directory, paste0("processed_", filename))
  
  # Read sequences
  annotated_species <- read_tsv(file)
  
  # Filter data
  metabolic_potential <- annotated_species %>%
    filter(!COG_category %in% c("-", "S", "R")) |> 
    mutate(count = 1 / nchar(COG_category)) |> 
    separate_rows(COG_category, sep = "") |> 
    filter(COG_category != "") |> 
    group_by(COG_category) |> 
    summarise( total = sum(count)) |> 
    mutate(relab = (total / sum(total)*100))

  
  # Create dataframe
  dataframe_percentage_metabolic_potential <- metabolic_potential |> 
    select(COG_category, relab) |> 
    rename(Category = COG_category, Percentage.Freq = relab)
  
  # Save to file
  write_tsv(dataframe_percentage_metabolic_potential, output_file)
}

# Process all files
lapply(files, process_file)

# Generate table
# List all files in the input directory
freqs <- list.files(output_directory, pattern = "\\.annotations$", full.names = TRUE)

# Function to read a file and add a column for the file name
read_and_label_file <- function(file) {
  df <- read_tsv(file, show_col_types = FALSE)
  # Use the file name (without path and extension) as the column name
  col_name <- tools::file_path_sans_ext(basename(file))
  colnames(df)[colnames(df) == 'Percentage.Freq'] <- col_name
  return(df)
}

# Read and label all files
file_dfs <- lapply(freqs, read_and_label_file)

# Perform a left join on all data frames by 'Category'
combined_df <- reduce(file_dfs, full_join, by = "Category")

# Remove the prefix and suffix from column names
prefix <- "processed_COG_curated_"
suffix <- ".fasta"

# Update column names
colnames(combined_df) <- gsub(paste0("^", prefix), "", colnames(combined_df))
colnames(combined_df) <- gsub(paste0(suffix, "$"), "", colnames(combined_df))

transposed_combined_df <- as.data.frame(t(combined_df))
new_colnames <- as.character(transposed_combined_df[1, ])
transposed_combined_df <- transposed_combined_df[-1, ]
colnames(transposed_combined_df) <- new_colnames

# Save
# Convert row names to a column
transposed_combined_df <- transposed_combined_df %>%
  rownames_to_column(var = "RowName")

# Save the data frame to a TSV file
write_tsv(
  transposed_combined_df, 
  "/home/agalvez/data/metabolic_amoeba/version4/results/amoebozoa_check_tips_COG_percentage_freq.tsv"
)
