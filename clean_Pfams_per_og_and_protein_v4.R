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
sequence_annotations <- read.table(file = "/home/agalvez/data/metabolic_amoeba/version4/data/general/pfam_combined_annotations_v4.txt", sep = "\t", quote = "", header = FALSE, fill = TRUE, col.names = c("Gene","Pfam"))

#Merge sequences per OG
og_sequence <- unite(og_sequences_table, All_sequences, Acancast:Ylip, sep = ",")


## Separate the All_sequences
og_sequence_unnested <- og_sequence %>%
  separate_rows(All_sequences, sep = ",") %>%
  mutate(All_sequences = trimws(All_sequences)) 

## Join to get Pfam values
pfam_per_og <- og_sequence_unnested %>%
  left_join(sequence_annotations, by = c("All_sequences" = "Gene")) %>%
  select(Orthogroup, Pfam) |> 
  filter(!is.na(Pfam)) |> 
  filter(Pfam != "-") |> 
  mutate(count = 1 / (str_count(Pfam, ",") + 1)) |> 
  separate_rows(Pfam, sep = ",") |> 
  group_by(Orthogroup, Pfam) |> 
  summarise( total = sum(count)) |> 
  mutate(relab = (total / sum(total)*100)) |> 
  dplyr::rename(count = total, percentage = relab)

## Alternative Pfam clan values
clans <- read_tsv(file = "/home/agalvez/data/metabolic_amoeba/version4/data/general/Pfam-A.clans.tsv", col_names = c("Pfam accession", "clan accession", "Clan_ID", "Pfam ID", "Pfam description")) |>
  select(c("Clan_ID", "Pfam ID")) |> 
  mutate(`Clan_ID` = if_else(is.na(`Clan_ID`), `Pfam ID`, `Clan_ID`))

clan_per_og <- left_join(pfam_per_og, clans, by = c("Pfam"="Pfam ID")) |> 
  select(!c(Pfam)) |> 
  group_by(Orthogroup, Clan_ID) |> 
  summarise(percentage = sum(percentage))


check_clan_per_og <- clan_per_og |> 
  group_by(Orthogroup) |> 
  summarise(
    total_percentage = sum(percentage, na.rm = TRUE)
  )


# Save translation table
write.csv(pfam_per_og,"/home/agalvez/data/metabolic_amoeba/version4/results/translate_OG_to_function_Pfams_v4.csv", row.names = FALSE)
write.csv(clan_per_og,"/home/agalvez/data/metabolic_amoeba/version4/results/translate_OG_to_function_clans_Pfams_v4.csv", row.names = FALSE)

