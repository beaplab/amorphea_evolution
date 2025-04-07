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
library(ggnewscale)


# Read ancestors
percentages_ancestors <- read_csv("/home/agalvez/data/metabolic_amoeba/version4/results/relevant_ancestors_v4.csv") |> 
  column_to_rownames(var = "File") |> 
  mutate(across(everything(), ~replace_na(., 0)))

# Read my COGs per species
percentages_per_species <- read_tsv("/home/agalvez/data/metabolic_amoeba/version4/results/tips_COG_percentage_freq_v4.tsv") |> 
  column_to_rownames(var = "RowName") |> 
  mutate(across(everything(), ~replace_na(., 0)))

# Reading the data
tips_taxons <- read_tsv("/home/agalvez/data/metabolic_amoeba/version4/data/tips_classes_and_orders_ancestors.txt")

rownames1 <- read_tsv("/home/agalvez/data/metabolic_amoeba/version4/results/tips_COG_percentage_freq_v4.tsv") %>% 
  select("RowName")
rownames2 <- read_csv("/home/agalvez/data/metabolic_amoeba/version4/results/relevant_ancestors_v4.csv") %>% 
  select("File") %>% 
  rename("RowName"="File")

rownames <- rbind(rownames1, rownames2)
  
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
amorphea_COG_percentages <- rbind(percentages_per_species,percentages_ancestors) |> 
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
merged_no_out <- merged_df %>% filter(Supergroup != c("Outgroup"))


# Function to create convex hull for each group
compute_hull <- function(df) {
  df[chull(df$Dim.1, df$Dim.2), ]
}

# Apply convex hull calculation to each Supergroup
hull_data <- merged_no_out %>%
  group_by(Supergroup) %>%
  do(compute_hull(.))

point_alpha <- c(
  "Fungi" = 0.15,
  "Olympicamoeba" = 0.15,
  "Metazoa" = 0.15,
  "Amoebozoa" = 0.15,
  "Ancestor-Bpp" = 1.0,
  "Ancestor-Notung" = 1.0
)


final_ca_supl <- ggplot() +
  # Add convex hulls for each Supergroup with custom fill colors
  geom_polygon(data = hull_data, aes(x = Dim.1, y = Dim.2, fill = Supergroup, color = Supergroup, group = Supergroup),
               alpha = 0.1, linetype = "dashed") +
  # Add species points with custom color, shape, and alpha
  geom_point(data = merged_no_out, aes(x = Dim.1, y = Dim.2, color = Supergroup, shape = Ancestor, alpha = Supergroup),
             size = 2) +
  # Add species names using ggrepel for specific Supergroups
  geom_text_repel(
    data = merged_no_out %>%
      filter(Supergroup %in% c("Ancestor-Bpp", "Ancestor-Notung")),
    aes(x = Dim.1, y = Dim.2, label = Species, color = Supergroup),
    size = 2, max.overlaps = Inf, segment.alpha = 0.3
  ) +
  # Add dashed lines for axis
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  # Customize fill for convex hulls
  scale_fill_manual(values = c("Fungi" = "blue", "Metazoa" = "red", "Amoebozoa" = "purple",
                               "Ancestor-Bpp" = "darkgreen", "Ancestor-Notung" = "springgreen3")) +
  # Customize point and label colors
  scale_color_manual(values = c("Amorphea" = "black", "Fungi" = "blue", "Metazoa" = "red",
                                "Olympicamoeba" = "orange2", "Amoebozoa" = "purple",
                                "Ancestor-Bpp" = "darkgreen", "Ancestor-Notung" = "springgreen3")) +
  # Customize alpha values
  scale_alpha_manual(values = point_alpha) +
  # Customize shapes
  scale_shape_manual(values = c("No" = 16, "Yes" = 17)) +
  # Minimal theme
  theme_minimal() +
  # Add title and axis labels
  labs(x = "Dimension 1 (45.4%)",
       y = "Dimension 2 (13.1%)")


ggsave("/home/agalvez/data/metabolic_amoeba/version4/results/final_figures/pre_inkscape_ancestors_bothmethods2_ca.pdf", plot = final_ca_supl, width = 10, height = 8)


##### Only Notung

merged_no_out <- merged_no_out %>%
  mutate(
    Modified_Species = case_when(
      Species == "Opisthokonta(Bpp)" ~ "Opisthokonta",
      Species == "Amoebozoa(Notung)" ~ "Amoebozoa",
      Species == "Metazoa(Bpp)" ~ "Metazoa",
      Species == "Choanozoa(Bpp)" ~ "Choanozoa",
      Species == "Holozoa(Bpp)" ~ "Holozoa",
      Species == "Obazoa(Bpp)" ~ "Obazoa",
      Species == "Discosea(Bpp)" ~ "Discosea",
      Species == "Evosea(Bpp)" ~ "Evosea",
      Species == "Tubulinea(Bpp)" ~ "Tubulinea",
      Species == "Amoebozoa(Bpp)" ~ "Amoebozoa",
      Species == "Amorphea(Bpp)" ~ "Amorphea",
      Species == "Fungi(Bpp)" ~ "Fungi",
      Species == "Nucletmycea(Bpp)" ~ "Nucetmycea",
      Species == "Obazoa(Notung)" ~ "Obazoa", 
      Species == "Opisthokonta(Notung)" ~ "Opisthokonta",
      Species == "Tubulinea(Notung)" ~ "Tubulinea",
      Species == "Nucletmycea(Notung)" ~ "Nucletmycea",
      Species == "Holozoa(Notung)" ~ "Holozoa",
      Species == "Discosea(Notung)" ~ "Discosea",
      Species == "Evosea(Notung)" ~ "Evosea",
      Species == "Choanozoa(Notung)" ~ "Choanozoa",
      Species == "Fungi(Notung)" ~ "Fungi",
      Species == "Metazoa(Notung)" ~ "Metazoa",
      Species == "Amorphea(Notung)" ~ "Amorphea",   
      TRUE ~ Species  # Keep other labels unchanged
    )
  )

# Define alpha values for hulls and points
hull_alpha <- c(
  "Fungi" = 0.05, 
  "Metazoa" = 0.05, 
  "Amoebozoa" = 0.05, 
  "Ancestor-Bpp" = 0.15, 
  "Ancestor-Notung" = 0.15
)

point_alpha <- c(
  "Fungi" = 0.15, 
  "Olympicamoeba" = 0.15, 
  "Metazoa" = 0.15, 
  "Amoebozoa" = 0.15, 
  "Ancestor-Bpp" = 1.0, 
  "Ancestor-Notung" = 1.0
)

# Create the plot
final_ca <- ggplot() + 
  # Convex hulls with specific alpha values
  geom_polygon(
    data = hull_data, 
    aes(x = Dim.1, y = Dim.2, fill = Supergroup, color = Supergroup, group = Supergroup, alpha = Supergroup), 
    linetype = "dashed"
  ) +  
  scale_alpha_manual(values = hull_alpha, guide = "none") +  # Alpha for hulls
  
  # Reset the alpha scale for points
  new_scale("alpha") +  
  
  # Species points with different alpha values
  geom_point(
    data = merged_no_out, 
    aes(x = Dim.1, y = Dim.2, color = Supergroup, shape = Ancestor, alpha = Supergroup), 
    size = 2
  ) +  
  scale_alpha_manual(values = point_alpha) +  # Alpha for points
  
  # Labels only for Ancestor-Bpp and Ancestor-Notung species
  geom_text_repel(
    data = merged_no_out %>% 
      filter(Supergroup %in% c("Ancestor-Bpp", "Ancestor-Notung")),
    aes(x = Dim.1, y = Dim.2, label = Modified_Species, color = Supergroup), 
    size = 2.5, max.overlaps = Inf, segment.alpha = 0.3, 
    force = 2, box.padding = 0.5, point.padding = 0.3
  ) +  
  
  # Dashed axis lines
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  
  # Customize fill colors for convex hulls
  scale_fill_manual(values = c(
    "Fungi" = "blue", "Metazoa" = "red", "Amoebozoa" = "purple", 
    "Ancestor-Bpp" = "darkgreen", "Ancestor-Notung" = "springgreen4"
  )) +  
  
  # Customize point and label colors
  scale_color_manual(values = c(
    "Amorphea" = "black", "Fungi" = "blue", "Metazoa" = "red", 
    "Olympicamoeba" = "orange2", "Amoebozoa" = "purple", 
    "Ancestor-Bpp" = "darkgreen", "Ancestor-Notung" = "springgreen4"
  )) +  
  
  # Customize shapes for Ancestor points
  scale_shape_manual(values = c("No" = 16, "Yes" = 17)) +
  
  theme_minimal() +
  labs(x = "Dimension 1 (45.7%)", y = "Dimension 2 (14.6%)"
  )

ggsave("/home/agalvez/data/metabolic_amoeba/version4/results/final_figures/pre_inkscape_ancestors_notung_ca.pdf", plot = final_ca, width = 10, height = 8)
