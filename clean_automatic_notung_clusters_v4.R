library(tidyverse)
library(readxl)
library(dplyr)
library(readr)
library(tidyr)
library(KEGGREST)
library(parallel)
library(patchwork)
library(ape)     
library(pheatmap) 
library(cluster)
library(factoextra)
library(dynamicTreeCut)
library(RColorBrewer)
library(colorspace)
library(ggpubr) 
library(svglite)

# Read Orthogroup(OG) counts
notung_raw <- read_tsv(file = "/home/agalvez/data/metabolic_amoeba/version4/data/general/ancestors/batch_file_v4.txt.SpeciesTree_rooted_node_labels.txt.geneCount.txt")
notung_raw$`GT name` <- sub("_tree\\.txt$", "", notung_raw$`GT name`)

n_columns <- colnames(notung_raw)[grep("^N\\d", colnames(notung_raw))]

notung_raw <- notung_raw |>
  select(!c(n_columns))

notung_raw <- as.data.frame(notung_raw) |> 
  column_to_rownames(var = "GT name")

# Read annotations
Pfam_OG <- read.table("/home/agalvez/data/metabolic_amoeba/version4/results/translate_OG_to_function_clans_Pfams_v4.csv", sep = ',', header = TRUE)|> 
  rename(Pfam = Clan_ID)

# Keep only Metazoa, Fungi and Amoebozoa with high buscos
columns_to_remove <- c('MicrA134', 'Squajapo', 'Ovaldese', 'Mayocant', 'ThecSK13', 'Echibisp', 'Malacali', 'Brevanat', 'Mycagemm', 'Tychacut', 'Pellcata', 'Nutohowe', 'Colltric', 'Multmedi')

high_busco_vector <- c("Echiexud", "Flabcita", "Idiovort", "Tadh",    
                       "Arcevulg", "Spom", "Hetepall", "Tmel",    
                       "Physpoly", "Amac", "Schivulg", "Pbla", "Dpul",    
                       "Lbic", "Hvul", "Dictdisc",
                       "Nvec", "Bisobp66", "Scom", "Skow", "Ylip",    
                       "Balamand", "Cochminu", "Spur", "Cang", "Bflo",    
                       "Protfung", "Lgig", "Rhizsaxo", "Vermverm",
                       "Vsepbp79", "Crev", "Ccor", "Luaphula", "Soliirre",
                       "Scil", "Hsap", "Coprprot", "Vermanta", "Spro", "Ccin", "Nirr", "Spun",    
                       "Masteilh", "Umay", "Ctel", "Aque",
                       "Mver", "Gpro", "Cavoapop", "Ncra", "Mlei")

pre_high_notung <- notung_raw |> 
  select(-all_of(columns_to_remove))

high_notung <- pre_high_notung |> 
  select(all_of(high_busco_vector))

# Filter absent Orthogroups (OG) in 90% of species
filtered_high_notung <- high_notung[rowMeans(high_notung == 0) < 0.9, ]

# Normalize counts
normalized_high_notung <- log10(filtered_high_notung + 0.001)  # Apply log transformation

# Load taxonomic information 
tip_groups <- read.table("/home/agalvez/data/metabolic_amoeba/version4/data/tips_classes_and_orders_ancestors.txt", header = TRUE, sep = "\t")

class <- tip_groups |> 
  select(c("Notung", "Class")) |> 
  column_to_rownames(var = "Notung")

class <- class[rownames(class) %in% colnames(normalized_high_notung), , drop = FALSE]

# Cluster OG counts 
pheatmap_high_counts <- pheatmap(normalized_high_notung,
                                 show_rownames = FALSE,
                                 show_colnames = TRUE,
                                 cluster_rows = TRUE,
                                 cluster_cols = TRUE,
                                 clustering_distance_cols = "euclidean",  
                                 clustering_method_cols = "complete",
                                 clustering_distance_rows = "manhattan",  # good for expression data
                                 clustering_method_rows = "complete",
                                 annotation_col = class)

# Extract the hierarchical clustering results from rows
row_dendrogram <- pheatmap_high_counts$tree_row  

# Cut the dendrogram to generate clusters
height_threshold <- 85  # Adjust based on the dendrogram
row_clusters <- cutree(row_dendrogram, h = height_threshold)
num_clusters <- length(unique(row_clusters))  # Count detected clusters

# Order clusters based on appearance in the dendrogram
ordered_clusters <- factor(row_clusters, levels = unique(row_clusters[row_dendrogram$order]))

# Apply color gradient based on order
cluster_colors <- colorRampPalette(brewer.pal(12, "Set3"))(num_clusters)
ordered_colors <- setNames(cluster_colors, levels(ordered_clusters))

# Create annotation data
cluster_annotation <- data.frame(Cluster = ordered_clusters)
rownames(cluster_annotation) <- rownames(normalized_high_notung)
annotation_row_colors <- list(Cluster = ordered_colors)

# Pheatmap to visualize clusters
clustered_pheatmap <- pheatmap(normalized_high_notung,
                               show_rownames = FALSE,
                               show_colnames = TRUE,
                               cluster_rows = TRUE,
                               cluster_cols = TRUE,
                               clustering_distance_cols = "euclidean",
                               clustering_method_cols = "complete",
                               clustering_distance_rows = "manhattan",  
                               clustering_method_rows = "complete",
                               annotation_col = class,
                               annotation_row = cluster_annotation,
                               annotation_colors = annotation_row_colors)

ggsave(
  filename = "/home/agalvez/data/metabolic_amoeba/version4/results/final_figures/pre_inkscape_dendrogram_cca.png",  # File name
  plot = clustered_pheatmap,           # Plot object
  width = 14,                     # Width in inches
  height = 13,                     # Height in inches
  dpi = 700  )


# Plot legend separately (as we have so many clusters that it does not fit the previous graph)
check_legend <- barplot(rep(1, length(ordered_colors)), 
                        col = rev(ordered_colors),  # Reverse color order
                        border = "black", 
                        names.arg = rev(names(ordered_colors)),  # Reverse labels
                        horiz = TRUE,  # Flip bars to horizontal
                        las = 1,       # Rotate text to be readable
                        cex.names = 0.4, 
                        main = "Cluster Colors")


# Select and merge clusters to fit Supergroups
panamorphea <- cluster_annotation |> 
  filter(Cluster %in% c(4,28,15,5,2,1,3,41,23,62,29,72,58))

metazoa <- cluster_annotation |> 
  filter(Cluster %in% c(12,44,95))

fungi <- cluster_annotation |> 
  filter(Cluster%in% c(34,71,56,96,78,94))

amoebozoa <- cluster_annotation |> 
  filter(Cluster %in% c(37,7,70,31,89,74,47,35,6,20,40))  

# Combine all groups into one dataframe fro visualization purposes
check_clusters_df <- bind_rows(panamorphea, metazoa, fungi, amoebozoa) |> select(Cluster)
check_clusters_df$Cluster <- as.character(check_clusters_df$Cluster)

annotation_row_colors <- list(
  Cluster = c(
    setNames(rep("cyan", length(c(4,28,15,5,2,1,3,41,23,62,29,72,58))), c(4,28,15,5,2,1,3,41,23,62,29,72,58)),
    setNames(rep("red", length(c(12,44,95))), c(12,44,95)),
    setNames(rep("green", length(c(34,71,56,96,78,94))), c(34,71,56,96,78,94)),
    setNames(rep("purple", length(c(37,7,70,31,89,74,47,35,6,20,40))), c(37,7,70,31,89,74,47,35,6,20,40))
  )
)

# Check selected clusters
check_clusters_pheatmap <- pheatmap(normalized_high_notung,
                                   show_rownames = FALSE,
                                   show_colnames = TRUE,
                                   cluster_rows = TRUE,
                                   cluster_cols = TRUE,
                                   clustering_distance_cols = "euclidean",
                                   clustering_method_cols = "complete",
                                   clustering_distance_rows = "manhattan",  # better for expression data
                                   clustering_method_rows = "complete",
                                   annotation_col = class,
                                   annotation_row = check_clusters_df,
                                   annotation_colors = annotation_row_colors)

ggsave(
  filename = "/home/agalvez/data/metabolic_amoeba/version4/results/final_figures/pre_inkscape_check_dendrogram_cca.png",  # File name
  plot = check_clusters_pheatmap,           # Plot object
  width = 14,                     # Width in inches
  height = 13,                     # Height in inches
  dpi = 700  )

# Detect Orthogroups in clusters
panamorphea <- panamorphea |> rownames_to_column(var = "Orthogroup")
metazoa <- metazoa |> rownames_to_column(var = "Orthogroup")
fungi <- fungi |> rownames_to_column(var = "Orthogroup")
amoebozoa <- amoebozoa |> rownames_to_column(var = "Orthogroup")


# Detect clans
preparing_data <- function(df){
  
  dataset <- df |> 
    left_join(Pfam_OG, by = "Orthogroup") |> 
    group_by(Pfam) |> 
    summarize(Score = sum(percentage, na.rm = TRUE)) |> 
    filter(!is.na(Pfam)) |>
    mutate(percent = (Score / sum(Score))*100) 
 
  return(dataset)
  
}

sel.pfams.panamo <- preparing_data(panamorphea)
sel.pfams.metazoa <- preparing_data(metazoa)
sel.pfams.fungi <- preparing_data(fungi) 
sel.pfams.amoebozoa <- preparing_data(amoebozoa)


# Selection of relevant clans
top_hits_f <- function(df1,df2,df3,df4){

  dataset1 <- df1 |>
    arrange(desc(Score)) |>
    slice_head(n = 20) |>
    pull(Pfam)

  dataset2 <- df2 |>
    arrange(desc(Score)) |>
    slice_head(n = 20) |>
    pull(Pfam)

  dataset3 <- df3 |>
    arrange(desc(Score)) |>
    slice_head(n = 20) |>
    pull(Pfam)

  dataset4 <- df4 |>
    arrange(desc(Score)) |>
    slice_head(n = 20) |>
    pull(Pfam)

  genes <- unique(c(dataset1,dataset2,dataset3,dataset4))

  return(genes)

}

gene_list <- top_hits_f(sel.pfams.amoebozoa,sel.pfams.fungi,sel.pfams.metazoa,sel.pfams.panamo)


# Keep only relevant clans
preparing_relevant <- function(df){
  
  dataset <- df |> 
    left_join(Pfam_OG, by = "Orthogroup") |> 
    group_by(Pfam) |> 
    summarize(Score = sum(percentage, na.rm = TRUE)) |> 
    filter(!is.na(Pfam)) |>
    mutate(percent = (Score / sum(Score))*100) |>
    filter(Pfam %in% gene_list)
  return(dataset)
  
}

rel.pfams.panamo <- preparing_relevant(panamorphea)
rel.pfams.metazoa <- preparing_relevant(metazoa)
rel.pfams.fungi <- preparing_relevant(fungi) 
rel.pfams.amoebozoa <- preparing_relevant(amoebozoa)

# Complete rel.pfams with 0 to missing clans
add_missing_pfams <- function(df, gene_list) {
  missing_pfams <- setdiff(gene_list, df$Pfam)
  
  if (length(missing_pfams) > 0) {
    missing_df <- data.frame(Pfam = missing_pfams, Score = 0, percent = 0)
    
    df <- bind_rows(df, missing_df)
  }
  
  return(df)
}

rel.pfams.panamo <- add_missing_pfams(rel.pfams.panamo, gene_list)
rel.pfams.metazoa <- add_missing_pfams(rel.pfams.metazoa, gene_list)
rel.pfams.fungi <- add_missing_pfams(rel.pfams.fungi, gene_list)
rel.pfams.amoebozoa <- add_missing_pfams(rel.pfams.amoebozoa, gene_list)

# Combine rel.pfams and keep track of which Supergroup they represent
bind_data <- bind_rows(
  rel.pfams.panamo %>% mutate(Group = "Panamorphea"),
  rel.pfams.amoebozoa %>% mutate(Group = "Amoebozoa"),
  rel.pfams.fungi %>% mutate(Group = "Fungi"),
  rel.pfams.metazoa %>% mutate(Group = "Metazoa"))


# Plot a heatmap representing relative presences of each clan in each Supergroup gene cluster
heatmap_plot <- ggplot(bind_data, aes(x = Group, y = Pfam, fill = percent)) +
  geom_tile(color = "black") +  # Add black borders around tiles
  scale_fill_gradient(low = "#f7fbff", high = "red4", name = "Percent") +  # Color scale for percent
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    axis.text.y = element_text(size = 10),             # Adjust y-axis label size
    panel.grid = element_blank()                       # Remove grid lines
  ) +
  labs(
    x = "Group",
    y = "Pfam"
  ) +
  coord_fixed()

ggsave(
  filename = "/home/agalvez/data/metabolic_amoeba/version4/results/final_figures/pre_inkscape_dendrogram_automatic_function_overrepresented_pfams_clan_annotation.pdf",  # File name
  plot = heatmap_plot,           # Plot object
  width = 16,                     # Width in inches
  height = 20,                     # Height in inches
  dpi = 300                       # Resolution
)

# Select Pfam clans that are overrepresented for each Supergroup
supergroup_tables <- split(bind_data, bind_data$Group)
fungi_table <- supergroup_tables[["Fungi"]]
panamo_table <- supergroup_tables[["Panamorphea"]]
metazoa_table <- supergroup_tables[["Metazoa"]]
amoebozoa_table <- supergroup_tables[["Amoebozoa"]]

differential_clans <- bind_data %>%
  group_by(Pfam) %>%
  arrange(desc(percent), .by_group = TRUE) %>% # Sort within each Pfam by percent
  summarise(
    Top_Group = first(Group),
    Top_Percent = first(percent),
    Second_Percent = nth(percent, 2, default = 0), # Handle cases with <2 entries
    Difference = Top_Percent - Second_Percent, 
    Ratio = ifelse(Second_Percent == 0, Inf, Top_Percent / Second_Percent) ) |> 
  filter(Top_Percent >= 0.5 & Ratio >= 1.4)  # 40% increase arbitrary 

write.csv(differential_clans, "/home/agalvez/data/metabolic_amoeba/version4/results/final_figures/post_inkscape/to_paper/definitive/differentially_expressed_clan_per_group.csv", row.names = FALSE)

