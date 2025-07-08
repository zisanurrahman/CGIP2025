# Optimized Shiny App using .rds files for fast loading
library(shiny)
library(ggplot2)
library(plotly)
library(DT)
library(dplyr)
library(readr)
#install.packages(c("ggforce", "scales"))
library(ggforce)
library(scales)

# === Load pre-converted .rds files ===
df <- reactiveVal(readRDS("data/volcano.rds"))
lipinski_df <- reactiveVal(readRDS("data/lipinski.rds"))
umap_df <- reactiveVal(readRDS("data/umap.rds"))
compound_images <- readRDS("data/compound_images.rds")
tsne_data <- readRDS("data/tsne_data.rds")

# === UI ===
ui <- fluidPage(
  titlePanel(HTML("<h2 style='text-align: center;'>Explore a Chemical-Genetic Interactions Profile of <i>B. cenocepacia</i> K56-2</h2>")),
  
  wellPanel(
    tags$div(
      tags$p("This app accesses data from: Rahman ASMZ et al. 202X. XXX. The complete dataset contains over 3 million chemical-genetic interactions generated from screening of a CRISPRi-based essential gene knockdown mutant library against 5000 compounds."),
      tags$p("For clarity and visualization purposes, only the compounds exhibiting strong interactions (Z-Score ≤ -5) with at least one essential gene knockdown mutant are displayed here."),
      tags$p("Special Note: t-SNE plot here projects the high-dimensional (609 knockdown mutants × mutant response) CGIPs for each compound onto two dimensions. The analysis concatenates the features (knockdown mutants and their response to each compound) into a single list for each compound, ensuring all profiles are padded to the same length. t-SNE is applied to reduce the high-dimensional data to two dimensions, creating a 2D representation of each compound profile. The 2D t-SNE results are then used as input for a K-means clustering algorithm. K-means partitioned the data into a predefined number of clusters (in this case, 10), grouping compounds based on the similarity of their 2D representations. The selected compound is highlighted in red and its corresponding cluster—determined via k-means clustering in the t-SNE space—is shaded with the same color to reflect its localized embedding.")
    )
  ),
  
  fluidRow(
    column(3, selectInput("compound", "Select Compound:", choices = NULL)),
    column(4, sliderInput("zscore", "Filter by Z-Score:",
                          min = -25, max = 10, value = c(-25, 10), step = 0.5))
  ),
 
   fluidRow(
    column(4, tags$div(style = "text-align: center;",
                       tags$img(src = "https://zenodo.org/records/14994629/files/heatmap_with_dendrogram_all5K.png",
                                style = "max-width: 100%; height: auto;"))),
    column(5, tags$div(style = "text-align: center;",
                       tags$img(src = "https://zenodo.org/records/14994629/files/Number_of_Compounds_vs_Percentage_of_Mutants_mean_240415_Compressed.png",
                                style = "max-width: 100%; height: auto;"))),
    column(3, h3("Most Depleted Mutants"), DTOutput("topMutants2"))
  )
,
  
  fluidRow(
    column(3, h3("Selected Compound"), uiOutput("lipinskiImage")),
    column(4, plotlyOutput("volcanoPlot", height = "320px")),
    column(5, plotOutput("umapPlot", height = "320px"))
  ),
  
  fluidRow(
    column(12, h3("Computed Chemical Properties"), DTOutput("lipinskiTable"))
  )
)

# === Server ===
server <- function(input, output, session) {
  
  observe({
    updateSelectInput(session, "compound", choices = unique(df()$Compounds))
  })
  
  filtered_data <- reactive({
    req(input$compound)
    df() %>%
      filter(Compounds == input$compound,
             Z.Score_ESet_mean >= input$zscore[1],
             Z.Score_ESet_mean <= input$zscore[2])
  })
  
  top_mutants_data <- reactive({
    filtered_data() %>%
      arrange(Z.Score_ESet_mean) %>%
      head(10) %>%
      select(Mutants, Z.Score_ESet_mean, COG2, Gene) %>%
      rename("Z-Score" = Z.Score_ESet_mean, "COG" = COG2)
  })
  
  output$topMutants2 <- renderDT({
    top_mutants_data()
  }, options = list(pageLength = 5))
  
  output$volcanoPlot <- renderPlotly({
    data <- filtered_data() %>%
      mutate(color = case_when(
        Z.Score_ESet_mean <= -1 ~ "red",
        Z.Score_ESet_mean >= 1 ~ "blue",
        TRUE ~ "grey"
      ))
    
    p <- ggplot(data, aes(x = Z.Score_ESet_mean, y = -log10(modified_padj),
                          color = color,
                          text = paste("Gene:", Gene, "<br>Mutants:", Mutants))) +
      geom_point(alpha = 0.7) +
      scale_color_identity() +
      labs(title = paste("Volcano Plot for", input$compound),
           x = "Z-Score", y = "-log10(Padj)") +
      theme_minimal() + theme(legend.position = "none")
    
    ggplotly(p, tooltip = "text") %>% layout(height = 320, showlegend = FALSE)
  })
  
  output$lipinskiImage <- renderUI({
    img_data <- compound_images[[input$compound]]
    if (!is.null(img_data)) {
      tags$img(src = img_data, height = "150px")
    } else {
      tags$p("Image not available.")
    }
  })
  
  output$lipinskiTable <- renderDT({
    data <- lipinski_df() %>% filter(Compounds == input$compound)
    if (nrow(data) > 0) {
      data %>% select(-SMILES)
    }
  }, options = list(dom = 't'))
  
  output$umapPlot <- renderPlot({
    req(input$compound)
    
    library(ggplot2)
    library(dplyr)
    library(scales)
    
    df <- tsne_data
    df$Compounds_clean <- trimws(tolower(df$Compounds))
    selected_compound <- trimws(tolower(input$compound))
    
    # Cluster
    set.seed(42)
    kmeans_result <- kmeans(df[, c("V1", "V2")], centers = 10)
    df$cluster <- as.factor(kmeans_result$cluster)
    
    # Define groupings and color mapping
    color_mapping <- c("rpoB" = "blue", "gyrB" = "red", "30s rRNA" = "orange", "dcw" = "purple")
    
    # Identify Fill group of selected compound
    selected_fill <- df %>% filter(Compounds_clean == selected_compound) %>% pull(Fill) %>% .[1]
    selected_color <- if (selected_fill %in% names(color_mapping)) color_mapping[selected_fill] else "darkred"
    
    # Set color and alpha per point
    df$point_color <- "grey"
    df$point_alpha <- 0.5
    df$point_size <- 2.5
    
    df$point_color[df$Compounds_clean == selected_compound] <- selected_color
    df$point_alpha[df$Compounds_clean == selected_compound] <- 1
    df$point_size[df$Compounds_clean == selected_compound] <- 2.5
    
    # Plot base
    p <- ggplot(df, aes(x = V1, y = V2)) +
      geom_point(aes(color = point_color, alpha = point_alpha), size = 2.5) +
      scale_color_identity() +
      scale_alpha_identity() +
      labs(x = "t-SNE Dimension I", y = "t-SNE Dimension II") +
      theme_minimal() +
      theme(
        legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_text(size = 16)
      )
    
    # Determine selected compound's cluster
    selected_cluster <- df %>%
      filter(Compounds_clean == selected_compound) %>%
      pull(cluster) %>% .[1]
    
    # Draw hulls
    for (cl in unique(df$cluster)) {
      cluster_data <- df %>% filter(cluster == cl)
      if (nrow(cluster_data) > 2) {
        hull_indices <- chull(cluster_data$V1, cluster_data$V2)
        hull_data <- cluster_data[hull_indices, ]
        
        if (cl == selected_cluster) {
          hull_fill <- alpha(selected_color, 0.15)
          hull_border <- selected_color
        } else {
          hull_fill <- alpha("grey", 0.12)
          hull_border <- "grey"
        }
        
        p <- p + geom_polygon(data = hull_data, aes(x = V1, y = V2),
                              fill = hull_fill, color = hull_border, linewidth = 0.4)
      }
    }
    
    p
  })
  
  
}

# === Launch App ===
shinyApp(ui = ui, server = server)
