library(tidyverse)
library(readxl)
library(dplyr)
library(readr)
library(tidyr)
library(KEGGREST)
library(parallel)
library(patchwork)
library(viridis)
library(cowplot)

# Read ancestors
percentages_ancestors <- read_csv("/home/agalvez/data/metabolic_amoeba/version4/results/relevant_ancestors_v4.csv") |> 
  column_to_rownames(var = "File") |> 
  mutate(across(everything(), ~replace_na(., 0)))

# Keep only NOTUNG inferences
tips_taxons <- read_tsv("/home/agalvez/data/metabolic_amoeba/version4/data/tips_classes_and_orders_ancestors.txt")

rownames <- read_csv("/home/agalvez/data/metabolic_amoeba/version4/results/relevant_ancestors_v4.csv") %>% 
  select("File") %>% 
  rename("RowName"="File")

tips_taxons <- left_join(rownames, tips_taxons, by = c("RowName" = "Notung"))

rows_to_remove <- tips_taxons |>
  filter(Supergroup %in% c("Ancestor-Bpp")) |> 
  pull(RowName)

percentages_ancestors <- percentages_ancestors |> 
  rownames_to_column(var = "rowname") |> 
  filter(!rowname %in% rows_to_remove) |> 
  column_to_rownames(var = "rowname")

# Calculate percentage changes per node in reference to ancestral node content
# Amoebozoa
amoebozoa <- ((percentages_ancestors["Amoebozoa(Notung)", ] - percentages_ancestors["Amorphea(Notung)", ]) /
              percentages_ancestors["Amorphea(Notung)", ]) * 100
rownames(amoebozoa) <- "Amoebozoa(changes)"

# Discosea
discosea <- ((percentages_ancestors["Discosea(Notung)", ] - percentages_ancestors["Amoebozoa(Notung)", ]) /
               percentages_ancestors["Amoebozoa(Notung)", ]) * 100
rownames(discosea) <- "Discosea(changes)"

# Tubulinea
tubulinea <- ((percentages_ancestors["Tubulinea(Notung)", ] - percentages_ancestors["Amoebozoa(Notung)", ]) /
                percentages_ancestors["Amoebozoa(Notung)", ]) * 100
rownames(tubulinea) <- "Tubulinea(changes)"

# Evosea
evosea <- ((percentages_ancestors["Evosea(Notung)", ] - percentages_ancestors["Amoebozoa(Notung)", ]) /
             percentages_ancestors["Amoebozoa(Notung)", ]) * 100
rownames(evosea) <- "Evosea(changes)"

# Obazoa
obazoa <- ((percentages_ancestors["Obazoa(Notung)", ] - percentages_ancestors["Amorphea(Notung)", ]) /
             percentages_ancestors["Amorphea(Notung)", ]) * 100
rownames(obazoa) <- "Obazoa(changes)"

# Opisthokonta
opisthokonta <- ((percentages_ancestors["Opisthokonta(Notung)", ] - percentages_ancestors["Obazoa(Notung)", ]) /
                   percentages_ancestors["Obazoa(Notung)", ]) * 100
rownames(opisthokonta) <- "Opisthokonta(changes)"

# Nucletmycea
nucletmycea <- ((percentages_ancestors["Nucletmycea(Notung)", ] - percentages_ancestors["Opisthokonta(Notung)", ]) /
                  percentages_ancestors["Opisthokonta(Notung)", ]) * 100
rownames(nucletmycea) <- "Nucletmycea(changes)"

# Fungi
fungi <- ((percentages_ancestors["Fungi(Notung)", ] - percentages_ancestors["Nucletmycea(Notung)", ]) /
            percentages_ancestors["Nucletmycea(Notung)", ]) * 100
rownames(fungi) <- "Fungi(changes)"

# Holozoa
holozoa <- ((percentages_ancestors["Holozoa(Notung)", ] - percentages_ancestors["Opisthokonta(Notung)", ]) /
              percentages_ancestors["Opisthokonta(Notung)", ]) * 100
rownames(holozoa) <- "Holozoa(changes)"

# Choanozoa
choanozoa <- ((percentages_ancestors["Choanozoa(Notung)", ] - percentages_ancestors["Holozoa(Notung)", ]) /
                percentages_ancestors["Holozoa(Notung)", ]) * 100
rownames(choanozoa) <- "Choanozoa(changes)"

# Metazoa
metazoa <- ((percentages_ancestors["Metazoa(Notung)", ] - percentages_ancestors["Choanozoa(Notung)", ]) /
              percentages_ancestors["Choanozoa(Notung)", ]) * 100
rownames(metazoa) <- "Metazoa(changes)"

# Metazoa_from_Amorphea
metazoa_from_amorphea <- ((percentages_ancestors["Metazoa(Notung)", ] - percentages_ancestors["Amorphea(Notung)", ]) /
                            percentages_ancestors["Amorphea(Notung)", ]) * 100
rownames(metazoa_from_amorphea) <- "Metazoa_from_Amorphea(changes)"

# Fungi_from_Amorphea
fungi_from_amorphea <- ((percentages_ancestors["Fungi(Notung)", ] - percentages_ancestors["Amorphea(Notung)", ]) /
                          percentages_ancestors["Amorphea(Notung)", ]) * 100
rownames(fungi_from_amorphea) <- "Fungi_from_Amorphea(changes)"


# Bind new rows to the original dataframe
percentages_ancestors <- rbind(percentages_ancestors, amoebozoa, discosea, tubulinea, evosea, obazoa, opisthokonta, nucletmycea, fungi, holozoa, choanozoa, metazoa, fungi_from_amorphea, metazoa_from_amorphea)

# Keep only changes
changes <- percentages_ancestors[grep("\\(changes\\)", rownames(percentages_ancestors)), ]

changes$Node <- rownames(changes)

changes_long <- pivot_longer(changes, 
                             cols = -Node, 
                             names_to = "Variable", 
                             values_to = "Value")


changes_long$Node <- str_replace(changes_long$Node, "\\(changes\\)", "")


# Select relevant categories checking changes superior to 20 %.
# Used to decide threshold of relevant COGs (we selected 20%)
ggplot(changes_long, aes(x = Value)) +
     geom_density(fill = "lightblue", alpha = 0.5) +
     theme_minimal() +
     labs(title = "Density Plot of 'value'",
         x = "Value",
         y = "Density")

chosen_nodes <- c("Amoebozoa", "Fungi_from_Amorphea", "Metazoa_from_Amorphea")

detect_changes_in_paths <-changes_long %>%
  filter(Node %in% chosen_nodes, Value > 20 | Value < -20)

relevant_COGs <- unique(detect_changes_in_paths$Variable)

# AMOEBOZOA: No relevant changes in composition from Amorphea. 
# METAZOA: Expansion of T, Z and W. Reduction of Y, H, J, M and N.
# FUNGI: Expansion of K, A, B and Y. Reduction of N, M, V and W.

changes_long <- changes_long |>
  filter(Variable %in% relevant_COGs)


# Generate one plot per changes in each node
label_mapping <- c(
  "J" = "J = Translation, ribosomal structure and biogenesis",  
  "A" = "A = RNA processing and modification", 
  "K" = "K = Transcription", 
  "L" = "L = Replication, recombination and repair", 
  "B" = "B = Chromatin structure and dynamics", 
  "D" = "D = Cell cycle control, cell division, chromosome partitioning", 
  "Y" = "Y = Nuclear structure", 
  "V" = "V = Defense mechanisms", 
  "T" = "T = Signal transduction mechanisms", 
  "M" = "M = Cell wall/membrane/envelope biogenesis", 
  "N" = "N = Cell motility", 
  "Z" = "Z = Cytoskeleton", 
  "W" = "W = Extracellular structures", 
  "U" = "U = Intracellular trafficking, secretion, and vesicular transport", 
  "O" = "O = Posttranslational modification, protein turnover, chaperones", 
  "C" = "C = Energy production and conversion", 
  "G" = "G = Carbohydrate transport and metabolism", 
  "E" = "E = Amino acid transport and metabolism", 
  "F" = "F = Nucleotide transport and metabolism", 
  "H" = "H = Coenzyme transport and metabolism", 
  "I" = "I = Lipid transport and metabolism", 
  "P" = "P = Inorganic ion transport and metabolism", 
  "Q" = "Q = Secondary metabolites biosynthesis, transport and catabolism", 
  "R" = "R = General function prediction only", 
  "S" = "S = Function unknown"
)

plots_list <- list()

for (node_name in unique(changes_long$Node)) {
  node_data <- changes_long[changes_long$Node == node_name, ]
  
  p <- ggplot(node_data, aes(x = Variable, y = Value, fill = Variable)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +  # Set alpha for bars
    labs(
      # title = paste(node_name), 
      x = "COG category", 
      y = "Change in composition (%)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
    scale_fill_manual(values = viridis::viridis(length(label_mapping)),  # Use viridis colors for custom labels
                      labels = label_mapping) +  # Adjust legend labels
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5),  # Center the title
          legend.title = element_blank()) + # Hide the legend title
    geom_hline(yintercept = 0, color = "black", linetype = "solid") + # Add a horizontal line at y = 0
    ylim(-20, 30)

  plots_list[[node_name]] <- p
}

plots_list[[4]]

# Loop through the list and save each plot
# for (node_name in names(plots_list)) {
#   plot <- plots_list[[node_name]]
# 
#   # Save as PNG
#   ggsave(paste0("/home/agalvez/data/metabolic_amoeba/version4/results/gains_losses/", node_name, "_percentages_increment.png"), plot = plot, width = 22, height = 10)
# }

# Combined paths
combined_paths_plot <- plot_grid(plots_list[[1]] + theme(legend.position = "none"), plots_list[[12]] + theme(legend.position = "none"), plots_list[[13]] + theme(legend.position = "none"),
                                 labels = c("Amoebozoa", "Fungi", "Metazoa"), ncol = 1)  # Arrange in a single column
ggsave("/home/agalvez/data/metabolic_amoeba/version4/results/final_figures/pre_inkscape_differential_cogs_combined_paths_%increment_plot.pdf", combined_paths_plot, width = 18, height = 10)  # Adjust dimensions as needed

combined_tax2_paths_plot <- plot_grid(plots_list[[2]] + theme(legend.position = "none"), plots_list[[3]] + theme(legend.position = "none"), plots_list[[4]] + theme(legend.position = "none"),
                                      labels = c("Discosea", "Tubulinea", "Evosea"), ncol = 1)  # Arrange in a single column
ggsave("/home/agalvez/data/metabolic_amoeba/version4/results/final_figures/pre_inkscape_differential_cogs_tax2_paths_%increment_plot.pdf", combined_tax2_paths_plot, width = 18, height = 10)  # Adjust dimensions as needed