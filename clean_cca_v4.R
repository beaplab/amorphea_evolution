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
  mutate(across(everything(), ~replace_na(., 0)))

# Reading the data
    tips_taxons <- read_tsv("/home/agalvez/data/metabolic_amoeba/version4/data/tips_classes_and_orders_ancestors.txt")
    rownames <- read_tsv("/home/agalvez/data/metabolic_amoeba/version4/results/tips_COG_percentage_freq_v4.tsv") %>% 
      select("RowName")

# Left join by RowName in rownames and Notung in tips_taxons
    tips_taxons <- left_join(rownames, tips_taxons, by = c("RowName" = "Notung"))

    rows_to_remove1 <- tips_taxons |>
      filter(Supergroup %in% c("Outgroup")) |>
      pull(RowName)

    rows_to_remove2 <- tips_taxons |>
      filter(RowName %in% c('MicrA134', 'Squajapo', 'Ovaldese', 'Mayocant', 'ThecSK13', 'Echibisp', 'Malacali', 'Brevanat', 'Mycagemm', 'Tychacut', 'Pellcata', 'Nutohowe', 'Colltric', 'Multmedi')) |>
      pull(RowName)

    rows_to_remove <- c(rows_to_remove1, rows_to_remove2) # Remove Outgroups and low Busco Amoebozoa

    # Merge both tables
amorphea_COG_percentages <- percentages_per_species |> 
  rownames_to_column(var = "rowname") |> 
  filter(!rowname %in% rows_to_remove) |> 
  column_to_rownames(var = "rowname")


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
merged_no_out <- merged_df


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
  geom_point(data = merged_no_out, aes(x = Dim.1, y = Dim.2, color = Supergroup), size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +   # Add dashed lines for axis
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = c("Fungi" = "blue", "Metazoa" = "red", "Amoebozoa" = "purple")) +  # Customize fill for the convex hulls and color for points and labels
  scale_color_manual(values = c("Amorphea" = "black", "Fungi" = "blue", "Metazoa" = "red",
                                "Olympicamoeba" = "orange2", "Amoebozoa" = "purple")) + # Olympicamoeba is a former candidate name for Apostamoebida
  theme_minimal() +
  labs(x = "Dimension 1 (45.5%)",
       y = "Dimension 2 (15.1%)")

final_names <- ggplot() + geom_polygon(data = hull_data, aes(x = Dim.1, y = Dim.2, fill = Supergroup, color = Supergroup, group = Supergroup),
                                    alpha = 0.1, linetype = "dashed") +   # Add convex hulls for each Supergroup with custom fill colors
  geom_point(data = merged_no_out, aes(x = Dim.1, y = Dim.2, color = Supergroup), size = 2) +
  geom_text_repel(data = merged_no_out, aes(x = Dim.1, y = Dim.2, label = Species, color = Supergroup),
                  size = 2.5, max.overlaps = Inf, segment.alpha = 0.3) +   # Use ggrepel to add species names, matching label color with point color, and make labels smaller
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +   # Add dashed lines for axis
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = c("Fungi" = "blue", "Metazoa" = "red", "Amoebozoa" = "purple")) +  # Customize fill for the convex hulls and color for points and labels
  scale_color_manual(values = c("Amorphea" = "black", "Fungi" = "blue", "Metazoa" = "red",
                                "Olympicamoeba" = "orange2", "Amoebozoa" = "purple")) +
  theme_minimal() +
  labs(x = "Dimension 1 (45.5%)",
       y = "Dimension 2 (15.1%)")


final_ca_points <- ggplot() +
  geom_segment(data = eigenvectors, aes(x = 0, y = 0, xend = Dim.1, yend = Dim.2),
               arrow = arrow(length = unit(0.3, "cm")), color = "black", alpha = 0.4) +  # Add eigenvectors as arrows
  geom_text(data = eigenvectors, aes(x = Dim.1, y = Dim.2, label = variable), vjust = 1.5, color = "black") +  # Add species points with smaller size and color based on Supergroup
  geom_point(data = merged_no_out, aes(x = Dim.1, y = Dim.2, color = Supergroup), size = 4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +   # Add dashed lines for axis
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = c("Fungi" = "blue", "Metazoa" = "red", "Amoebozoa" = "purple")) +  # Customize fill for the convex hulls and color for points and labels
  scale_color_manual(values = c("Amorphea" = "black", "Fungi" = "blue", "Metazoa" = "red", 
                                "Olympicamoeba" = "orange2", "Amoebozoa" = "purple")) +
theme_minimal() +
  labs(x = "Dimension 1 (45.5%)",
       y = "Dimension 2 (15.1%)")

ggsave("/home/agalvez/data/metabolic_amoeba/version4/results/final_figures/pre_inkscape_hull_all_ca.pdf", plot = final_ca, width = 10, height = 8)
ggsave("/home/agalvez/data/metabolic_amoeba/version4/results/final_figures/pre_inkscape_points_all_ca.pdf", plot = final_ca_points, width = 10, height = 8)
ggsave("/home/agalvez/data/metabolic_amoeba/version4/results/final_figures/pre_inkscape_names_all_ca.pdf", plot = final_names, width = 10, height = 8)


# Function to create convex hull for each group
compute_hull <- function(df) {
  df[chull(df$Dim.3, df$Dim.4), ]
}

# Apply convex hull calculation to each Supergroup
hull_data_3_4 <- merged_no_out %>%
  group_by(Supergroup) %>%
  do(compute_hull(.))

# Plot for dimensions 3 and 4
final_ca_3_4 <- ggplot() +
  #geom_polygon(data = hull_data_3_4, aes(x = Dim.3, y = Dim.4, fill = Supergroup, color = Supergroup, group = Supergroup), 
   #            alpha = 0.1, linetype = "dashed") +
  geom_point(data = merged_no_out, aes(x = Dim.3, y = Dim.4, color = Supergroup, shape = Ancestor), size = 4) + 
  geom_segment(data = eigenvectors, aes(x = 0, y = 0, xend = Dim.3, yend = Dim.4),
               arrow = arrow(length = unit(0.3, "cm")), color = "black", alpha = 0.8) +
  geom_text(data = eigenvectors, aes(x = Dim.3, y = Dim.4, label = variable), vjust = 1.5, color = "black") +
  # geom_text_repel(data = merged_no_out, aes(x = Dim.3, y = Dim.4, label = Species, color = Supergroup), 
  #                size = 3, max.overlaps = Inf) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = c("Fungi" = "blue", "Metazoa" = "red", "Amoebozoa" = "purple")) +
  scale_color_manual(values = c("Amorphea" = "black", "Fungi" = "blue", "Metazoa" = "red", 
                                "Olympicamoeba" = "orange", "Amoebozoa" = "purple")) +
  scale_shape_manual(values = c("No" = 16, "Uramoebozoan" = 17, "Urfungi" = 17, "Urmetazoan" = 17, "Uramorphean" = 17)) +
  theme_minimal() +
  labs(title = "CA Biplot: Dimension 3 vs Dimension 4",
       x = "Dimension 3 (9.7%)",
       y = "Dimension 4 (8.3%)")


# Add Amoebozoa class info
# Function to create convex hull for each group
compute_hull <- function(df) {
  df[chull(df$Dim.1, df$Dim.2), ]
}

# Apply convex hull calculation to each Supergroup
hull_data_class <- merged_no_out %>%
  group_by(Class) %>%
  do(compute_hull(.))

final_ca_class <- ggplot() +
  geom_point(data = merged_no_out, aes(x = Dim.1, y = Dim.2, color = Class, shape = Class), size = 3) + 
  geom_segment(data = eigenvectors, aes(x = 0, y = 0, xend = Dim.1, yend = Dim.2),
               arrow = arrow(length = unit(0.3, "cm")), color = "black", alpha = 0.4) +  # Add eigenvectors as arrows
  geom_text(data = eigenvectors, aes(x = Dim.1, y = Dim.2, label = variable), vjust = 1.5, color = "black") +  # Add species points with smaller size and color based on Supergroup
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +   # Add dashed lines for axis
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Fungi" = "blue", "Metazoa" = "red", "Tubulinea" = "purple", "Discosea" = "purple4", "Evosea" = "pink3")) +   # Customize point and label colors
  scale_shape_manual(values = c("Fungi" = 19, "Metazoa" = 20, "Tubulinea" = 18, "Discosea" = 17, "Evosea" = 15)) +
  theme_minimal() +
  labs(x = "Dimension 1 (45.5%)",
       y = "Dimension 2 (15.1%)")

ggsave("/home/agalvez/data/metabolic_amoeba/version4/results/final_figures/pre_inkscape_class_and_others.pdf", plot = final_ca_class, width = 10, height = 8)

