library(tidyverse)
library(readxl)
library(dplyr)
library(readr)
library(tidyr)
library(KEGGREST)
library(parallel)
library(FactoMineR)
library(factoextra)
library(ggrepel)

# Read my COGs per species
percentages_per_species <- read_tsv("/home/agalvez/data/metabolic_amoeba/version4/results/tips_COG_percentage_freq_v4.tsv") |> 
  column_to_rownames(var = "RowName") |> 
  mutate(across(everything(), ~replace_na(., 0))) %>% 
  rownames_to_column(var = "RowName")

# Reading the data
tips_taxons <- read_tsv("/home/agalvez/data/metabolic_amoeba/version4/data/tips_classes_and_orders_ancestors.txt")
rownames <- read_tsv("/home/agalvez/data/metabolic_amoeba/version4/results/tips_COG_percentage_freq_v4.tsv") %>% 
  select("RowName")

# Left join by RowName in rownames and Notung in tips_taxons
tips_taxons <- left_join(rownames, tips_taxons, by = c("RowName" = "Notung"))

#  Amoebozoa cluster, excluding Dictdisc and Hetepall
amoebozoa_cluster <- tips_taxons %>%
  filter(Supergroup == "Amoebozoa", !RowName %in% c("Dictdisc", "Hetepall"))

#  The two species of interest
dictis_cluster <- tips_taxons %>%
  filter(RowName %in% c("Dictdisc", "Hetepall"))

#  Fungi cluster
fungi_cluster <- tips_taxons %>%
  filter(Supergroup == "Fungi")

# Filter data
amoebozoa_percentages <- percentages_per_species %>%
  filter(RowName %in% amoebozoa_cluster$RowName) %>%
  column_to_rownames(var = "RowName")

dictis_percentages <- percentages_per_species %>%
  filter(RowName %in% dictis_cluster$RowName) %>%
  column_to_rownames(var = "RowName")

fungi_percentages <- percentages_per_species %>%
  filter(RowName %in% fungi_cluster$RowName) %>%
  column_to_rownames(var = "RowName")

# Mean values
amoebozoa_means <- colMeans(amoebozoa_percentages)
amoebozoa_means_df <- data.frame(
  COG = names(amoebozoa_means),
  MeanValue = as.numeric(amoebozoa_means)
)

dictis_means <- colMeans(dictis_percentages)
dictis_means_df <- data.frame(
  COG = names(dictis_means),
  MeanValue = as.numeric(dictis_means)
)

fungi_means <- colMeans(fungi_percentages)
fungi_means_df <- data.frame(
  COG = names(fungi_means),
  MeanValue = as.numeric(fungi_means)
)

# Distances
# Make sure both dataframes are sorted by COG
fungi_means_sorted <- fungi_means_df %>% arrange(COG)
dictis_means_sorted <- dictis_means_df %>% arrange(COG)
amoebozoa_means_sorted <- amoebozoa_means_df %>% arrange(COG)

# Absolute difference to Fungi
distance_to_fungi <- fungi_means_sorted %>%
  mutate(DictValue = dictis_means_sorted$MeanValue,
         Distance = abs(MeanValue - DictValue)) %>%
  select(COG, Distance)

# Absolute difference to Amoebozoa
distance_to_amoebozoa <- amoebozoa_means_sorted %>%
  mutate(DictValue = dictis_means_sorted$MeanValue,
         Distance = abs(MeanValue - DictValue)) %>%
  select(COG, Distance)

# Add a label column
distance_to_fungi$Group <- "Fungi"
distance_to_amoebozoa$Group <- "Amoebozoa"

# Combine both
distance_combined <- bind_rows(distance_to_fungi, distance_to_amoebozoa)

# Pivot into wide format
distance_wide <- distance_combined %>%
  pivot_wider(names_from = Group, values_from = Distance)

# Determine which reference is closer and compute differences
distance_classified <- distance_wide %>%
  mutate(
    CloserTo = ifelse(Fungi < Amoebozoa, "Fungi", "Amoebozoa"),
    AbsDifference = abs(Fungi - Amoebozoa),                # absolute difference
    RelDifference = AbsDifference / pmax(Fungi, Amoebozoa) # relative difference
  )

distance_classified_filtered <- distance_classified %>%
  filter(RelDifference > 0.5, AbsDifference >= 0.1)

# Sum distances per group
total_distances <- distance_combined %>%
  group_by(Group) %>%
  summarise(TotalDistance = sum(Distance, na.rm = TRUE))




##################

#  Amoebozoa cluster, excluding Dictdisc and Hetepall
amoebozoa_cluster <- tips_taxons %>%
  filter(RowName %in% c("Acancast"))

#  The two species of interest
dictis_cluster <- tips_taxons %>%
  filter(RowName %in% c("Dictdisc", "Hetepall"))

#  Fungi cluster
fungi_cluster <- tips_taxons %>%
  filter(RowName %in% c("Amac"))

# Filter
amoebozoa_percentages <- percentages_per_species %>%
  filter(RowName %in% amoebozoa_cluster$RowName) %>%
  column_to_rownames(var = "RowName")

dictis_percentages <- percentages_per_species %>%
  filter(RowName %in% dictis_cluster$RowName) %>%
  column_to_rownames(var = "RowName")

fungi_percentages <- percentages_per_species %>%
  filter(RowName %in% fungi_cluster$RowName) %>%
  column_to_rownames(var = "RowName")

# Mean values
amoebozoa_means <- colMeans(amoebozoa_percentages)
amoebozoa_means_df <- data.frame(
  COG = names(amoebozoa_means),
  MeanValue = as.numeric(amoebozoa_means)
)

dictis_means <- colMeans(dictis_percentages)
dictis_means_df <- data.frame(
  COG = names(dictis_means),
  MeanValue = as.numeric(dictis_means)
)

fungi_means <- colMeans(fungi_percentages)
fungi_means_df <- data.frame(
  COG = names(fungi_means),
  MeanValue = as.numeric(fungi_means)
)

# Distances
# Make sure both dataframes are sorted by COG
fungi_means_sorted <- fungi_means_df %>% arrange(COG)
dictis_means_sorted <- dictis_means_df %>% arrange(COG)
amoebozoa_means_sorted <- amoebozoa_means_df %>% arrange(COG)

# Absolute difference to Fungi
distance_to_fungi <- fungi_means_sorted %>%
  mutate(DictValue = dictis_means_sorted$MeanValue,
         Distance = abs(MeanValue - DictValue)) %>%
  select(COG, Distance)

# Absolute difference to Amoebozoa
distance_to_amoebozoa <- amoebozoa_means_sorted %>%
  mutate(DictValue = dictis_means_sorted$MeanValue,
         Distance = abs(MeanValue - DictValue)) %>%
  select(COG, Distance)

# Add a label column
distance_to_fungi$Group <- "Fungi"
distance_to_amoebozoa$Group <- "Amoebozoa"

# Combine both
distance_combined <- bind_rows(distance_to_fungi, distance_to_amoebozoa)

# Pivot into wide format
distance_wide <- distance_combined %>%
  pivot_wider(names_from = Group, values_from = Distance)

# Determine which reference is closer and compute differences
distance_classified <- distance_wide %>%
  mutate(
    CloserTo = ifelse(Fungi < Amoebozoa, "Fungi", "Amoebozoa"),
    AbsDifference = abs(Fungi - Amoebozoa),                # absolute difference
    RelDifference = AbsDifference / pmax(Fungi, Amoebozoa) # relative difference
  )

distance_classified_filtered2 <- distance_classified %>%
  filter(RelDifference > 0.5, AbsDifference >= 0.1)

# Sum distances per group
total_distances2 <- distance_combined %>%
  group_by(Group) %>%
  summarise(TotalDistance = sum(Distance, na.rm = TRUE))



#Consistency
# Ensure both tables only have COG and CloserTo columns
df1 <- distance_classified_filtered %>% select(COG, CloserTo)
df2 <- distance_classified_filtered2 %>% select(COG, CloserTo)

# Find COGs that belong to the same CloserTo group in both tables
consistent_cogs <- inner_join(df1, df2, by = "COG") %>%
  filter(CloserTo.x == CloserTo.y) %>%
  select(COG, CloserTo = CloserTo.x)

inconsistent_cogs <- inner_join(df1, df2, by = "COG") %>%
  filter(CloserTo.x != CloserTo.y) %>%
  select(COG, CloserTo1 = CloserTo.x, CloserTo2 = CloserTo.y)


# Q,A Fungi-like
# Z,E,D,C Amoebozoa-like



