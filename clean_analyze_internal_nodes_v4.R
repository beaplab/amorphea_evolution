library(tidyverse)
library(readxl)
library(dplyr)
library(readr)
library(tidyr)
library(KEGGREST)
library(parallel)
library(patchwork)


# Read notung gene count
notung_raw <- read_tsv(file = "/home/agalvez/data/metabolic_amoeba/version4/data/general/ancestors/batch_file_v4.txt.SpeciesTree_rooted_node_labels.txt.geneCount.txt")
notung_raw$`GT name` <- sub("_tree\\.txt$", "", notung_raw$`GT name`)
notung_raw <- notung_raw |>
  select("GT name", starts_with("N")) |> 
  select(!c("Nutohowe", "Ncra", "Nirr", "Nvec", "NolaAFSM"))

# Read Bppancestor gene count
bpp_raw <- read_tsv(file = "/home/agalvez/data/metabolic_amoeba/version4/data/general/ancestors/input_to_bpp.txt.sites_1")
bpp_raw$Sites <- gsub("\\[|\\]", "", bpp_raw$Sites)
bpp_raw$Sites <- as.numeric(bpp_raw$Sites) - 1
bpp_raw$Sites <- sprintf("%07d", bpp_raw$Sites)
bpp_raw$Sites <- paste0("OG", bpp_raw$Sites)
bpp_raw <- bpp_raw[c("Sites", grep("^prob", names(bpp_raw), value = TRUE))]
bpp_raw <- bpp_raw[c(names(bpp_raw)[!grepl("0$", names(bpp_raw))])]

clean_column_name <- function(col_name) {
  if (str_starts(col_name, "prob.")) {
    col_name <- str_replace(col_name, "prob.", "")  # Remove "prob."
    col_name <- str_split(col_name, "\\.", simplify = TRUE)[1]  # Remove ".1" and keep the first part
  }
  return(col_name)
}

names(bpp_raw) <- sapply(names(bpp_raw), clean_column_name)


# Read annotations 
COG_OG <- read.table("/home/agalvez/data/metabolic_amoeba/version4/results/translate_OG_to_function_COGs_v4.csv", sep = ',', header = TRUE) |> 
  select(-count)

cog_human <- read.table("/home/agalvez/data/metabolic_amoeba/version4/data/general/translation_COG.csv", sep = ';', header = TRUE)



## COG
### Notung
calculate_cog_ponderation_notung <- function(node, notung_raw, COG_OG, cog_human) {
  column_data <- notung_raw %>%
    select(`GT name`, all_of(node)) %>%
    left_join(COG_OG, by = c(`GT name` = "Orthogroup")) %>%
    filter(!!sym(node) != 0)
  
  final_data_notung_cog <- column_data %>%
    mutate(ponderation = !!sym(node) * percentage) %>%
    select(ponderation, COG) %>%
    filter(COG != "S") |> 
    group_by(COG) %>%
    summarize(total_count = sum(ponderation)) %>%
    left_join(cog_human, by = "COG") 

  
  all_counts <- sum(final_data_notung_cog$total_count)
  
  final_data_notung_cog <- final_data_notung_cog %>%
    mutate(percentage = (total_count / all_counts) * 100)
  
  return(final_data_notung_cog)
  
}

plot_cog_ponderation_notung <- function(df){
  
  barplot <- df %>%
    ggplot(aes(x = percentage, y = reorder(Description, percentage), fill = General)) +
    geom_bar(stat = "identity") +
    labs(x = "Percentage", y = "Functional annotation (COGs)",
         fill = "General function", title = paste(node)) +
    theme_classic() +
    theme(legend.position = "bottom")
  
  return(barplot)
  
}

nodes_to_process <- colnames(notung_raw)[-1]  # Excluding the first column `GT name`

output_dir <- "/home/agalvez/data/metabolic_amoeba/version4/results/notung/COG"

barplots <- list()

for (node in nodes_to_process) {
  df <- calculate_cog_ponderation_notung(node, notung_raw, COG_OG, cog_human)
  plot <- plot_cog_ponderation_notung(df)
  barplots[[node]] <- plot  # Store the plot in the list
  filename <- paste0(output_dir, "/", node)
  write_csv(x = df, file = str_c(filename , ".csv"))
  ggsave(str_c(filename , "_barplot.png"), plot = plot, width = 14, height = 6)
  cat("Saved plot for column:", node, "\n")
}

print(barplots[["N11"]])

### Bppancestor
calculate_cog_ponderation_bpp<- function(node, bpp_raw, COG_OG, cog_human) {
  column_data <- bpp_raw %>%
    select(Sites, all_of(node)) %>%
    left_join(COG_OG, by = c(Sites = "Orthogroup")) %>%
    filter(!!sym(node) != 0)
  
  
  final_data_bpp_cog <- column_data %>%
    mutate(ponderation = !!sym(node) * percentage) %>%
    select(ponderation, COG) %>%
    filter(COG != "S") |> 
    group_by(COG) %>%
    summarize(total_count = sum(ponderation)) %>%
    left_join(cog_human, by = "COG")

  
  all_counts <- sum(final_data_bpp_cog$total_count)
  
  final_data_bpp_cog <- final_data_bpp_cog %>%
    mutate(percentage = (total_count / all_counts) * 100)
  
  return(final_data_bpp_cog)
  
}

plot_cog_ponderation_bpp <- function(df){
  
  barplot <- df %>%
    ggplot(aes(x = percentage, y = reorder(Description, percentage), fill = General)) +
    geom_bar(stat = "identity") +
    labs(x = "Percentage", y = "Functional annotation (COGs)",
         fill = "General function", title = paste(node)) +
    theme_classic() +
    theme(legend.position = "bottom")
  
  return(barplot)
  
}

nodes_to_process <- colnames(bpp_raw)[-1]  # Excluding the first column `Sites`

output_dir <- "/home/agalvez/data/metabolic_amoeba/version4/results/bpp/COG"

barplots <- list()

for (node in nodes_to_process) {
  df <- calculate_cog_ponderation_bpp(node, bpp_raw, COG_OG, cog_human)
  plot <- plot_cog_ponderation_bpp(df)
  barplots[[node]] <- plot  # Store the plot in the list
  filename <- paste0(output_dir, "/", node)
  write_csv(x = df, file = str_c(filename , ".csv"))
  ggsave(str_c(filename , "_barplot.png"), plot = plot, width = 14, height = 6)
  cat("Saved plot for column:", node, "\n")
}

print(barplots[["224"]])


# Path to manually selected ancestors directory (from the two previously generated directories per Notung and Bpp)
ancestors_dir <- "/home/agalvez/data/metabolic_amoeba/version4/results/relevant_ancestors"

# Obtener una lista de todos los archivos CSV en el directorio
files <- list.files(ancestors_dir, pattern = "*.csv", full.names = TRUE)

# Función para procesar cada archivo
process_ancestor <- function(file) {
  # Leer el archivo
  data <- read_csv(file)
  
  # Extraer el nombre del archivo sin la extensión
  file_name <- tools::file_path_sans_ext(basename(file))
  
  # Crear una fila con el nombre del archivo y los porcentajes de cada COG
  data %>%
    select(COG, percentage) %>%
    pivot_wider(names_from = COG, values_from = percentage) %>%
    mutate(File = file_name) %>%
    select(File, everything())
}

# Aplicar la función a todos los archivos y combinarlos en una tabla
combined_data <- lapply(files, process_ancestor) %>%
  bind_rows()

# Change names to human readable (This depends on your phylogeny)
combined_data <- combined_data |> 
  mutate(File = recode(File,
                       `224` = 'Amoebozoa(Bpp)',
                       `N11` = 'Amoebozoa(Notung)',
                       `172` = 'Discosea(Bpp)',
                       `N20` = 'Discosea(Notung)',
                       `205` = 'Evosea(Bpp)',
                       `N21` = 'Evosea(Notung)',
                       `223` = 'Tubulinea(Bpp)',
                       `N15` = 'Tubulinea(Notung)',
                       `225` = 'Amorphea(Bpp)',
                       `N8` = 'Amorphea(Notung)',
                       `101` = 'Metazoa(Bpp)',
                       `N39` = 'Metazoa(Notung)',
                       `49` = 'Fungi(Bpp)',
                       `N35` = 'Fungi(Notung)',
                       `109` = 'Opisthokonta(Bpp)',
                       `N12` = 'Opisthokonta(Notung)',
                       `225` = 'Amorphea(Bpp)',
                       `N8` = 'Amorphea(Notung)',
                       `117` = 'Obazoa(Bpp)',
                       `N10` = 'Obazoa(Notung)',
                       `108` = 'Holozoa(Bpp)',
                       `N17` = 'Holozoa(Notung)',
                       `107` = 'Choanozoa(Bpp)',
                       `N27` = 'Choanozoa(Notung)',
                       `59` = 'Nucletmycea(Bpp)',
                       `N16` = 'Nucletmycea(Notung)',
                       ))

# Bpp     Notung  Name
# 224     N11     Amoebozoa
# 101     N39     Metazoa
# 49      N35     Fungi
# 109     N12     Opisthokonta
# 225     N8      Amorphea
# -       -       -
# 117     N10     Obazoa
# 108     N17     Holozoa
# 107     N27     Choanozoa
# 59      N16     Nucletmycea
# -       -       -
# 223     N15     Tubulinea
# 205     N21     Evosea
# 172     N20     Discosea



# Save
write_csv(combined_data, "/home/agalvez/data/metabolic_amoeba/version4/results/relevant_ancestors_v4.csv")
