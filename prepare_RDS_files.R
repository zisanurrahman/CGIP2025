# create data directory if not exists
dir.create("data", showWarnings = FALSE)

# Download and convert to .rds
volcano_csv <- "https://zenodo.org/records/13117525/files/unique_drugs_filtered_sorted_volcano_shiny_filtered2.csv"
lipinski_csv <- "https://zenodo.org/records/14994629/files/lipinski_parameters_250308.csv"
umap_csv <- "https://zenodo.org/records/13117525/files/unique_drugs_umap_clusters.csv"

volcano_df <- readr::read_csv(volcano_csv)
lipinski_df <- readr::read_csv(lipinski_csv)
umap_df <- readr::read_csv(umap_csv)

saveRDS(volcano_df, file = "data/volcano.rds")
saveRDS(lipinski_df, file = "data/lipinski.rds")
saveRDS(umap_df, file = "data/umap.rds")
pwd
getwd()
