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


# CA WITH METAZOA AND FUNGI, NOT SHOWN IN VISUALIZATION!

# Read my COGs per species
percentages_per_species <- read_tsv("/home/agalvez/data/metabolic_amoeba/version4/results/tips_COG_percentage_freq_v4.tsv") |> 
  column_to_rownames(var = "RowName") |> 
  mutate(across(everything(), ~replace_na(., 0)))

gt_percentages_per_species <- read_tsv("/home/agalvez/data/metabolic_amoeba/version4/results/amoebozoa_check_tips_COG_percentage_freq.tsv") |> 
  column_to_rownames(var = "RowName") |> 
  mutate(across(everything(), ~replace_na(., 0))) # Add extra Amoebozoa genomes and transcriptomes

# Reading the data
tips_taxons <- read_tsv("/home/agalvez/data/metabolic_amoeba/version4/data/tips_classes_and_orders_ancestors.txt")
rownames <- read_tsv("/home/agalvez/data/metabolic_amoeba/version4/results/tips_COG_percentage_freq_v4.tsv") %>% 
  select("RowName")

gt_tips_taxons <- read_tsv("/home/agalvez/data/metabolic_amoeba/version4/data/amoebozoa_check/short_amoebozoa_check_tips.txt") # Add extra Amoebozoa genome and transcriptome tags
gt_rownames <- read_tsv("/home/agalvez/data/metabolic_amoeba/version4/results/amoebozoa_check_tips_COG_percentage_freq.tsv") %>% 
  select("RowName")

# Left join by RowName in rownames and Notung in tips_taxons
tips_taxons <- left_join(rownames, tips_taxons, by = c("RowName" = "Notung"))

rows_to_remove1 <- tips_taxons |>
  filter(Supergroup %in% c("Outgroup","Metazoa","Fungi")) |>
  pull(RowName)

rows_to_remove2 <- tips_taxons |>
  filter(RowName %in% c('MicrA134', 'Squajapo', 'Ovaldese', 'Mayocant', 'ThecSK13', 'Echibisp', 'Malacali', 'Brevanat', 'Mycagemm', 'Tychacut', 'Pellcata', 'Nutohowe', 'Colltric', 'Multmedi')) |>
  pull(RowName)

rows_to_remove <- c(rows_to_remove1, rows_to_remove2)

gt_tips_taxons <- left_join(gt_rownames, gt_tips_taxons, by = c("RowName" = "Notung"))

gt_rows_to_remove <- gt_tips_taxons |> 
  filter(is.na(Supergroup)) |> 
  pull(RowName)

# Merge both tables
og_amorphea_COG_percentages <- percentages_per_species |> 
  rownames_to_column(var = "rowname") |> 
  filter(!rowname %in% rows_to_remove) |> 
  column_to_rownames(var = "rowname")

gt_amorphea_COG_percentages <- gt_percentages_per_species |> 
  rownames_to_column(var = "rowname") |> 
  filter(!rowname %in% gt_rows_to_remove) |> 
  column_to_rownames(var = "rowname")

amorphea_COG_percentages <- rbind(og_amorphea_COG_percentages, gt_amorphea_COG_percentages)

# Perform Correspondence Analysis
res.ca <- CA(amorphea_COG_percentages, graph = FALSE)

# Extract eigenvector coordinates (column points from CA)
eigenvectors <- data.frame(res.ca$col$coord)
eigenvectors$variable <- rownames(eigenvectors)

# Visualize the proportion of variance explained by each dimension
fviz_screeplot(res.ca, addlabels = TRUE, ylim = c(0, 70))

# Biplot of rows and columns
fviz_ca_biplot(res.ca, repel = TRUE)

# Plot row points
fviz_ca_row(res.ca, repel = TRUE)

# Plot column points
fviz_ca_col(res.ca, repel = TRUE)

# Add colors Fungi Metazoa Outgroups Amoebozoa 
tip_groups <- read.table("/home/agalvez/data/metabolic_amoeba/version4/data/tips_classes_and_orders_ancestors.txt", header = TRUE, sep = "\t")

# Merge CA results with tip_groups data
ca_df <- data.frame(res.ca$row$coord)
ca_df$Species <- rownames(ca_df)
merged_df <- left_join(ca_df, tip_groups, by = c("Species" = "Notung"))
merged_no_out <- merged_df |> 
  filter(!Supergroup %in% c("Fungi", "Metazoa")) |> 
  mutate(Species = recode(Species,
                          "Acast0RNA" = "AcancastRNA-redundans", 
                          "AcastRNA" = "AcancastRNA-cdhit", 
                          "VvermdnaT1" = "VermvermDNA-T1",
                          "VvermdnaT2" = "VermvermDNA-T2"))


# Function to create convex hull for each group
compute_hull <- function(df) {
  df[chull(df$Dim.1, df$Dim.2), ]
}

# Apply convex hull calculation to each Supergroup
hull_data <- merged_no_out %>%
  group_by(Supergroup) %>%
  do(compute_hull(.))

species_to_label <- merged_no_out |> 
  filter(Species %in% c("AcancastRNA-redundans", "AcancastRNA-cdhit", "VermvermDNA-T1", "VermvermDNA-T2", "Acancast", "Vermverm")
)
  
points_without_label <- merged_no_out |>
  filter(!Species %in% c("AcancastRNA-redundans", "AcancastRNA-cdhit", "VermvermDNA-T1", "VermvermDNA-T2", "Acancast", "Vermverm"))

final_ca <- ggplot() + geom_polygon(data = hull_data, aes(x = Dim.1, y = Dim.2, fill = Supergroup, color = Supergroup, group = Supergroup), 
                                    alpha = 0.1, linetype = "dashed") +   # Add convex hulls for each Supergroup with custom fill colors
  geom_point(data = species_to_label, aes(x = Dim.1, y = Dim.2, color = Supergroup, shape = Ancestor), size = 3, alpha = 1) + 
  # Points without labels (lower alpha)
  geom_point(data = points_without_label, aes(x = Dim.1, y = Dim.2, color = Supergroup, shape = Ancestor), size = 3, alpha = 0.3) +   geom_text_repel(data = species_to_label, aes(x = Dim.1, y = Dim.2, label = Species, color = Supergroup), 
                  size = 4.5, max.overlaps = Inf, segment.alpha = 0.3) +   # Use ggrepel to add species names, matching label color with point color, and make labels smaller
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +   # Add dashed lines for axis
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = c("Fungi" = "blue", "Metazoa" = "red", "Olympicamoeba" = "orange2", "AmoebozoaRR" = "purple", "AmoebozoaRC" = "purple2", "AmoebozoaGP" = "red4", "AmoebozoaGE" = "pink4", "Amoebozoa" = "purple4")) +  # Customize fill for the convex hulls and color for points and labels
  scale_color_manual(values = c("Amorphea" = "black", "Fungi" = "blue", "Metazoa" = "red",
                                "Olympicamoeba" = "orange2", "AmoebozoaRR" = "purple", "AmoebozoaRC" = "purple2", "AmoebozoaGP" = "red4", "AmoebozoaGE" = "pink4", "Amoebozoa" = "purple4")) +   # Customize point and label colors
  scale_shape_manual(values = c("No" = 16, "Uramoebozoan" = 17, "Urfungi" = 17, "Urmetazoan" = 17, "Uramorphean" = 17)) +
  theme_minimal() +
  labs(x = "Dimension 1 (34.8%)",
       y = "Dimension 2 (15.7%)") 

# THIS CONTAINS METAZOA AND FUNGI! THEY ARE JUST NOT PLOTTED!
ggsave("/home/agalvez/data/metabolic_amoeba/version4/results/final_figures/pre_inkscape_check_amoebozoa_ca.pdf", plot = final_ca, width = 10, height = 6)


################### Test placement Amoebozoa species per Busco value.
busco <- read.table("/home/agalvez/data/metabolic_amoeba/version4/data/general/busco.txt", header = FALSE)
colnames(busco) <- c("Species", "Busco")

# Sort the data frame by Busco in descending order
busco <- busco[order(-busco$Busco), ]

# Create a new column for the grades
busco$Grade <- NA  # Initialize the Grade column

# Assign grades based on the position in the sorted table
busco$Grade <- cut(
  busco$Busco,
  breaks = c(-Inf, 20, 40, 60, 80, Inf),  # Adjusted breaks for Busco values
  labels = c("E", "D", "C", "B", "A"),  # Grades from lowest to highest
  right = FALSE
)
species_to_filter <- busco |> 
  filter(Grade %in% c("D", "E"))

species_vector <- species_to_filter$Species   # This set of species was previously removed from our analysis
####################

#Species above 80% BUSCO
high_busco <- read.table("/home/agalvez/data/metabolic_amoeba/version4/data/general/busco.txt", header = FALSE)
colnames(high_busco) <- c("Species", "Busco")

high_busco <- high_busco |> 
  filter(Busco> 80)

high_busco_vector <- as.vector(high_busco$Species)
