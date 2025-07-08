#install.packages('jsonlite')
library(jsonlite)

json_path <- "data/compound_images.json"
compound_images <- fromJSON(json_path)
saveRDS(compound_images, "data/compound_images.rds")
