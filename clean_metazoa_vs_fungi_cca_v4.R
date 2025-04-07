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

#DO NOT PLOT OUTGROUPS? ADD SHAPES TO ANCESTORS!

# Read my COGs per species 
percentages_per_species <- read_tsv("/home/agalvez/data/metabolic_amoeba/version3/results/tips_COG_percentage_freq_v3.tsv") |> 
  column_to_rownames(var = "RowName") |> 
  mutate(across(everything(), ~replace_na(., 0))) # We use version 3 as it still has Rirr

mf_percentages_per_species <- read_tsv("/home/agalvez/data/metabolic_amoeba/version4/results/metazoa_vs_fungi_tips_COG_percentage_freq.tsv") |> 
  column_to_rownames(var = "RowName") |> 
  mutate(across(everything(), ~replace_na(., 0)))


# Reading the data
tips_taxons <- read_tsv("/home/agalvez/data/metabolic_amoeba/version4/data/tips_classes_and_orders_ancestors.txt")
rownames <- read_tsv("/home/agalvez/data/metabolic_amoeba/version3/results/tips_COG_percentage_freq_v3.tsv") %>% 
  select("RowName")

mf_tips_taxons <- read_tsv("/home/agalvez/data/metabolic_amoeba/version3/data/metazoa_vs_fungi/metazoa_vs_fungi_tips.txt")
mf_rownames <- read_tsv("/home/agalvez/data/metabolic_amoeba/version4/results/metazoa_vs_fungi_tips_COG_percentage_freq.tsv") %>% 
  select("RowName")

# Left join by RowName in rownames and Notung in tips_taxons
tips_taxons <- left_join(rownames, tips_taxons, by = c("RowName" = "Notung"))

rows_to_remove <- tips_taxons |> 
  filter(Supergroup %in% c("Outgroup","Amoebozoa","Olympicamoeba")) |> 
  pull(RowName)

mf_tips_taxons <- left_join(mf_rownames, mf_tips_taxons, by = c("RowName" = "Notung"))

mf_rows_to_remove <- mf_tips_taxons |> 
  filter(Supergroup %in% c("Outgroup","Amoebozoa","Olympicamoeba")) |> 
  pull(RowName)

# Merge both tables
og_amorphea_COG_percentages <- percentages_per_species |> 
  rownames_to_column(var = "rowname") |> 
  filter(!rowname %in% rows_to_remove) |> 
  column_to_rownames(var = "rowname")

mf_amorphea_COG_percentages <- mf_percentages_per_species |> 
  rownames_to_column(var = "rowname") |> 
  filter(!rowname %in% mf_rows_to_remove) |> 
  column_to_rownames(var = "rowname")

amorphea_COG_percentages <- rbind(og_amorphea_COG_percentages, mf_amorphea_COG_percentages)

# DO Correspondence analyses were done in R67 with FactoMiner package 
# - PCA NOT APPROPIATE for Compositional Data
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
merged_no_out <- merged_df %>% filter(Supergroup != "Outgroup")


# Function to create convex hull for each group
compute_hull <- function(df) {
  df[chull(df$Dim.1, df$Dim.2), ]
}

# Apply convex hull calculation to each Supergroup
hull_data <- merged_no_out %>%
  group_by(Supergroup) %>%
  do(compute_hull(.))

final_ca <- ggplot() + geom_polygon(data = hull_data, aes(x = Dim.1, y = Dim.2, fill = Supergroup, color = Supergroup, group = Supergroup), 
                                    alpha = 0.1, linetype = "dashed") +   # Add convex hulls for each Supergroup with custom fill colors
  geom_segment(data = eigenvectors, aes(x = 0, y = 0, xend = Dim.1, yend = Dim.2),
               arrow = arrow(length = unit(0.3, "cm")), color = "black", alpha = 0.4) +  # Add eigenvectors as arrows
  geom_text(data = eigenvectors, aes(x = Dim.1, y = Dim.2, label = variable), vjust = 1.5, color = "black") +  # Add species points with smaller size and color based on Supergroup
  geom_point(data = merged_no_out, aes(x = Dim.1, y = Dim.2, color = Supergroup, shape = Ancestor), size = 2) + 
  geom_text_repel(data = merged_no_out, aes(x = Dim.1, y = Dim.2, label = Species, color = Supergroup), 
                  size = 2.5, max.overlaps = Inf, segment.alpha = 0.3) +   # Use ggrepel to add species names, matching label color with point color, and make labels smaller
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +   # Add dashed lines for axis
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = c("Fungi" = "blue", "Metazoa" = "red", "FungiRNA" = "cyan3", "MetazoaRNA" = "red4")) +  # Customize fill for the convex hulls and color for points and labels
  scale_color_manual(values = c("Fungi" = "blue", "Metazoa" = "red", "FungiRNA" = "cyan3", "MetazoaRNA" = "red4")) +   # Customize point and label colors
  scale_shape_manual(values = c("No" = 16)) +
  theme_minimal() +
  labs(title = "CA Biplot with Eigenvectors, Custom Convex Hull Colors, and Supergroup Colors",
       x = "Dimension 1 (58%)",
       y = "Dimension 2 (12.3)")

ggsave("/home/agalvez/data/metabolic_amoeba/version4/results/final_figures/pre_inkscape_metazoa_vs_fungi_ca.pdf", plot = final_ca, width = 10, height = 7)

