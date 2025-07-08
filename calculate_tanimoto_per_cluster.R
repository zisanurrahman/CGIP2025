#install.packages("rcdk")
#install.packages('tidyverse')
#install.packages('fingerprint')
library(rcdk)
library(tidyverse)
library(fingerprint)

# Load data
tsne_data <- readRDS("data/tsne_data.rds")
lipinski_df <- readRDS("data/lipinski.rds")

# Add cluster
df <- tsne_data %>%
  mutate(Compounds_clean = trimws(tolower(Compounds)))

set.seed(42)
df$cluster <- as.factor(kmeans(df[, c("V1", "V2")], centers = 10)$cluster)

# Join SMILES
compound_smiles <- lipinski_df %>% select(Compounds, SMILES)
df <- df %>% left_join(compound_smiles, by = c("Compounds" = "Compounds"))

# Keep only compounds with SMILES
df <- df %>% filter(!is.na(SMILES))

# Precompute Tanimoto scores by cluster
tanimoto_scores <- df %>%
  group_by(cluster) %>%
  group_modify(~ {
    # Parse SMILES into molecules
    mols <- parse.smiles(.x$SMILES)
    fps <- lapply(mols, get.fingerprint, type = "extended")
    
    # Use the first compound in the cluster as reference
    ref_fp <- fps[[1]]
    sim_scores <- sapply(fps, function(fp) fp.sim.matrix(list(ref_fp, fp), method = "tanimoto")[1, 2])
    
    .x %>%
      mutate(TanimotoScore = round(sim_scores, 3))
  }) %>%
  ungroup() %>%
  select(Compounds, cluster, TanimotoScore)

# Save as .rds
saveRDS(tanimoto_scores, "data/tanimoto_scores_by_cluster.rds")