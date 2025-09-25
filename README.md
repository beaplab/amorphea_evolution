# amorphea_evolution

A curated collection of input data and analysis scripts used in the study:  
**_“The evolution of gene functional repertoire in Amorphea: Divergent strategies across Amoebozoa, Fungi, and Metazoa”_**

---

## 📜 Script Usage Overview

This repository documents all required **input files** and the **execution order** of scripts within our gene function evolution pipeline. Scripts are stored in Github, while input files are stored in Figshare (10.6084/m9.figshare.28741229).

---

## ⚙️ Pipeline Phases

### 🛠️ Phase 0: Input Files to Pipeline

| File | Description |
|------|-------------|
| `busco.txt` | Species BUSCO values |
| `short_amoebozoa_check_tips.txt` | Metadata for added Amoebozoa genomes/transcriptomes |
| `amoebozoa_check_curated_COG.tar.xz` | COGs per protein for Amoebozoa additions |
| `batch_file_v4.txt.SpeciesTree_rooted_node_labels.txt.geneCount.txt` | Notung-based gene gain/loss inference |
| `input_to_bpp.txt.sites_1` | Input for Bppancestor ancestral inference |
| `translation_COG.csv` | Abbreviations and descriptions of COG categories |
| `tips_COG_percentage_freq_v3.tsv` | Intermediate COG composition (Rirr excluded) |
| `general_curated_COG.tar.xz` | COGs per protein in core dataset species |
| `formatted_Orthogroups.tsv` | Orthogroup composition data |
| `COG_combined_annotations_v4.txt` | COG annotations for all proteins |
| `random_data.csv` | Randomly simulated species COG compositions |
| `metazoa_vs_fungi_tips.txt` | Transcriptome metadata for Metazoa and Fungi |
| `Pfam-A.clans.tsv` | Translation table from Pfam domains to clans |
| `pfam_combined_annotations_v4.txt` | Pfam annotations for all proteins |
| `script_inputs_and_outputs.csv` | Table of scripts with input/output file mappings |

---

### 🔍 Phase 1: Data Extraction & Preprocessing

1. `extract_COG_from_eggnog_without_hashtags.pl`  
   _Extracts raw COG annotations from EggNOG data._

2. `extract_Pfam_from_eggnog_without_hashtags.pl`  
   _Extracts Pfam clan information for use in downstream analyses._

---

### 📊 Phase 2: COG Processing

3. `clean_COG_per_species_and_protein_v4.R`  
   _Generates the main COG frequency matrix for downstream analyses._

4. `clean_COGs_per_og_and_protein_v4.R`  
   _Maps Orthogroups to COGs for ancestral inference._

5. `clean_analyze_internal_nodes_v4.R`  
   _Processes inferred ancestral nodes for comparative analysis._

---

### 🔄 Phase 3: Genome–Transcriptome Comparisons

6. `clean_amoebozoa_genome_transcriptome_check_COG_per_species_and_protein_v4.R`  
   _Genome–transcriptome comparison for Amoebozoa._

7. `clean_metazoa_vs_fungi_COG_per_species_and_protein_v4.R`  
   _Genome–transcriptome comparison for Metazoa vs. Fungi._

8. `clean_metazoa_vs_fungi_cca_v4.R`  
   _Generates Supplementary Figure 5 for Metazoa–Fungi._

9. `clean_amoebozoa_genome_transcriptome_check_cca_v4.R`  
   _Generates Supplementary Figure 5 for Amoebozoa._

---

### 🧠 Phase 4: Ancestor & Correspondence Analyses (CA)

10. `clean_cca_v4.R`  
    _Main CA plots (Figure 4, Supplementary Figure 3)._

11. `clean_ancestors_cca_v4.R`  
    _Performs CA on ancestral node data._

12. `clean_amoebozoa_cca_v4.R`  
    _Amoebozoa-specific CA (Supplementary Figures 6–7)._

13. `clean_percentual_gains_losses.R`  
    _Quantifies gene gain/loss at internal nodes (Supplementary Figure 8)._

14. `clean_metazoa_vs_amoebozoa_cca_v4.R`  
    _Generates Supplementary Figure 4 (Metazoa–Amoebozoa CA)._

15. `clean_random_check_cca_v4.R`  
    _Analyzes random data for comparison (Supplementary Figure 10)._

16. `clean_amoebozoa_vs_fungi_cca_v4.R`  
    _Generates Supplementary Figure 4 (Amoebozoa–Fungi CA)._

17. `clean_check_dictis_v4.R`  
    _Finds COG categories that are more similar to either Fungi or Amoebozoa in dictyostelids._

---

### 🔥 Phase 5: Pfam Clan Analyses

18. `clean_Pfams_per_og_and_protein_v4.R`  
    _Maps Pfam clans to orthogroups for clustering analysis._

19. `clean_automatic_notung_clusters_v4.R`  
    _Creates Pfam-based heatmaps (e.g., Figure 7)._
    
---   


## 🧾 Citation

If you use this repository, please cite the corresponding paper:  
> _The evolution of gene functional repertoire in Amorphea: Divergent strategies across Amoebozoa, Fungi, and Metazoa_

---

## 📫 Contact

For questions or contributions, feel free to reach out alexgalvez38@gmail.com or alex.galvez@ibe.upf-csic.es.

---

