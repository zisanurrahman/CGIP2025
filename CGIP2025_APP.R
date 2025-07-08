# Optimized Shiny App using .rds files for fast loading
library(shiny)
library(ggplot2)
library(plotly)
library(DT)
library(dplyr)
library(readr)

# === Load pre-converted .rds files ===
df <- reactiveVal(readRDS("data/volcano.rds"))
lipinski_df <- reactiveVal(readRDS("data/lipinski.rds"))
umap_df <- reactiveVal(readRDS("data/umap.rds"))

# === UI ===
ui <- fluidPage(
  titlePanel(HTML("<h2 style='text-align: center;'>Explore a Chemical-Genetic Interactions Profile of <i>B. cenocepacia</i> K56-2</h2>")),
  
  wellPanel(
    tags$div(
      tags$p("This app accesses data from: Rahman ASMZ et al. 202X. XXX. The complete dataset contains over 3 million chemical-genetic interactions generated from screening of a CRISPRi-based essential gene knockdown mutant library against 5000 compounds."),
      tags$p("For clarity and visualization purposes, only the compounds exhibiting strong interactions (Z-Score ≤ -5) with at least one essential gene knockdown mutant are displayed here."),
      tags$p("Special Note: UMAP (Uniform Manifold Approximation and Projection) plot here projects the 1827-dimensional (609 knockdown mutants × mutant response, GO annotation, and COG) CGIPs for each compound onto two dimensions. The analysis concatenates the features (knockdown mutants and their response, GO annotation, and COG) into a single list for each compound, ensuring all profiles are padded to the same length. UMAP is applied to reduce the high-dimensional data to 2 dimensions, creating a 2D representation of each compound profile. The 2D UMAP results are then used as input for a K-means clustering algorithm. K-means partitioned the data into a predefined number of clusters (in this case, 10), grouping compounds based on the similarity of their 2D representations.")
    )
  ),
  
  fluidRow(
    column(3, selectInput("compound", "Select Compound:", choices = NULL)),
    column(4, sliderInput("zscore", "Filter by Z-Score:",
                          min = -25, max = 10, value = c(-25, 10), step = 0.5))
  ),
  
  fluidRow(
    column(4, tags$img(src = "https://zenodo.org/records/14994629/files/heatmap_with_dendrogram_all5K.png",
                       width = "100%", height = "350px")),
    column(5, tags$img(src = "https://zenodo.org/records/14994629/files/Number_of_Compounds_vs_Percentage_of_Mutants_mean_240415_Compressed.png",
                       width = "102%", height = "340px")),
    column(3, h3("Most Depleted Mutants"), DTOutput("topMutants2"))
  ),
  
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
    data <- lipinski_df() %>% filter(Compounds == input$compound)
    if (nrow(data) > 0) {
      smiles <- data$SMILES[1]
      image_url <- paste0("https://cactus.nci.nih.gov/chemical/structure/", URLencode(smiles), "/image")
      tags$img(src = image_url, height = "150px")
    }
  })
  
  output$lipinskiTable <- renderDT({
    data <- lipinski_df() %>% filter(Compounds == input$compound)
    if (nrow(data) > 0) {
      data %>% select(-SMILES)
    }
  }, options = list(dom = 't'))
  
  output$umapPlot <- renderPlot({
    ggplot(umap_df(), aes(x = UMAP1, y = UMAP2, color = as.factor(Cluster))) +
      geom_point(alpha = 0.5) +
      geom_point(data = umap_df() %>% filter(Compounds == input$compound),
                 aes(x = UMAP1, y = UMAP2),
                 color = "black", size = 4, shape = 21, fill = "yellow") +
      labs(x = "UMAP 1", y = "UMAP 2", color = "Cluster") +
      theme_minimal() +
      theme(legend.position = c(0.5,0.6),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 14),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank())
  })
}

# === Launch App ===
shinyApp(ui = ui, server = server)
