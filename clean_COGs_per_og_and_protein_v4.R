library(tidyverse)
library(readxl)
library(dplyr)
library(readr)
library(tidyr)
library(KEGGREST)
library(parallel)


#Read OG sequences
og_sequences_table <- read_tsv(file = "/home/agalvez/data/metabolic_amoeba/version4/data/general/formatted_Orthogroups.tsv")
#Read annotations per sequence
sequence_annotations <- read.table(file = "/home/agalvez/data/metabolic_amoeba/version4/data/general/COG_combined_annotations_v4.txt", sep = "\t", quote = "", header = FALSE, fill = TRUE, col.names = c("Gene","COG"))
#Read path to human
cog_human <- read.table("/home/agalvez/data/metabolic_amoeba/version4/data/general/translation_COG.csv", sep = ';', header = TRUE)


#Merge sequences per OG
og_sequence <- unite(og_sequences_table, All_sequences, Acancast:Ylip, sep = ",")


## Separate the All_sequences
og_sequence_unnested <- og_sequence %>%
  separate_rows(All_sequences, sep = ",") %>%
  mutate(All_sequences = trimws(All_sequences))  

## Join to get COG values
cog_per_og <- og_sequence_unnested %>%
  left_join(sequence_annotations, by = c("All_sequences" = "Gene")) %>%
  select(Orthogroup, COG) |> 
  mutate(COG = ifelse(is.na(COG), "S", COG)) |> 
  mutate(COG = ifelse(COG == "R", "S", COG)) |> 
  filter(COG != "S") |> 
  filter(COG != "-") |> 
  mutate(count = 1 / nchar(COG)) |> 
  separate_rows(COG, sep = "") |> 
  filter(COG != "") |>
  group_by(Orthogroup, COG) |> 
  summarise( total = sum(count)) |> 
  mutate(relab = (total / sum(total)*100)) |> 
  rename(count = total, percentage = relab)



# save translation table
write.csv(cog_per_og,"/home/agalvez/data/metabolic_amoeba/version4/results/translate_OG_to_function_COGs_v4.csv", row.names = FALSE)
