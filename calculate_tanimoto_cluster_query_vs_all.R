#install.packages("rcdk")
#install.packages('tidyverse')
#install.packages('fingerprint')
library(rcdk)
library(fingerprint)
library(tidyverse)

# Load data
tsne_data <- readRDS("data/tsne_data.rds")
lipinski_df <- readRDS("data/lipinski.rds")

# Cluster compounds
tsne_data <- tsne_data %>%
  mutate(Compounds_clean = trimws(tolower(Compounds)))

set.seed(42)
tsne_data$cluster <- as.factor(kmeans(tsne_data[, c("V1", "V2")], centers = 10)$cluster)

# Merge SMILES
df <- tsne_data %>%
  left_join(lipinski_df %>% select(Compounds, SMILES), by = "Compounds") %>%
  filter(!is.na(SMILES))

# Function to compute similarity between one query and many targets
compute_tanimoto_against_all <- function(query_name, query_smiles, target_names, target_smiles, cluster_id) {
  # Parse query molecule
  query_mol <- parse.smiles(query_smiles)[[1]]
  query_fp <- get.fingerprint(query_mol, type = "extended")
  
  # Parse and fingerprint target molecules
  target_mols <- parse.smiles(target_smiles)
  target_fps <- lapply(target_mols, get.fingerprint, type = "extended")
  
  # Compute similarity of query_fp vs all target_fps
  sim_scores <- fp.sim.matrix(list(query_fp), target_fps, method = "tanimoto")[1, ]
  
  # Return long-form tibble
  tibble(
    Query = query_name,
    Target = target_names,
    Cluster = cluster_id,
    Tanimoto = round(sim_scores, 3)
  )
}


# Compute scores
tanimoto_scores_all <- list()

for (clus in unique(df$cluster)) {
  cat("Processing cluster:", clus, "\n")
  sub_df <- df %>% filter(cluster == clus)
  n <- nrow(sub_df)
  
  for (i in seq_len(n)) {
    query_name <- sub_df$Compounds[i]
    query_smiles <- sub_df$SMILES[i]
    target_names <- sub_df$Compounds
    target_smiles <- sub_df$SMILES
    
    scores <- compute_tanimoto_against_all(query_name, query_smiles, target_names, target_smiles, clus)
    tanimoto_scores_all[[length(tanimoto_scores_all) + 1]] <- scores
  }
}

# Combine and save
tanimoto_long <- bind_rows(tanimoto_scores_all)
saveRDS(tanimoto_long, "data/tanimoto_cluster_query_vs_all.rds")
