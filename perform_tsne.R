# Install and load required packages
library(Rtsne)
library(tidyverse)

# Read the CSV file
input_filename <- "unique_drugs_filtered_sorted_with_mutant_count_BCAL_count_updt.csv"
data <- read.csv(input_filename)
data
# Specify the column names for t-SNE
tsne_columns <- c("Locus_tag", "Compounds")

# Drop rows with missing values in the specified columns
selected_data <- na.omit(data[tsne_columns])

# Separate multiple values in the "Locus_tag" column into individual rows
selected_data <- selected_data %>%
  separate_rows(Locus_tag, sep = ",\\s*", convert = TRUE)

# Create a binary matrix indicating the presence or absence of each "Compounds" value for each "Locus_tag"
binary_matrix <- table(selected_data$Compounds, selected_data$Locus_tag) > 0
binary_matrix <- as.matrix(binary_matrix)

# Remove duplicate rows
numeric_data_locus <- unique(binary_matrix)

# Ensure numeric data for t-SNE
numeric_data_locus <- as.matrix(numeric_data_locus)

# Perform t-SNE
set.seed(42)
tsne_result_locus <- Rtsne(numeric_data_locus, dims = 2, perplexity = 10, verbose = TRUE)

# Convert t-SNE result to a data frame
tsne_data <- as.data.frame(tsne_result_locus$Y)
tsne_data
tsne_data$Compounds <- rownames(numeric_data_locus)
tsne_data
# Left join tsne_data with selected_data based on the "Compounds" column
merged_data <- left_join(tsne_data, selected_data, by = "Compounds")
merged_data
write.csv(merged_data,file='tsne_data_merged_file.csv', row.names=FALSE)