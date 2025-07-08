# ================================
# Shiny App: B. cenocepacia CGIP Viewer (Two Tabs)
# Author: Zisanur Rahman
# ================================

# --- Load Required Libraries ---

#install.packages("ggtext") 
library(shiny)
library(ggplot2)
library(plotly)
library(DT)
library(dplyr)
library(readr)
library(ggforce)
library(scales)
library(ggtext)

# === Load pre-converted .rds files ===
# Optimized Shiny App using .rds files for fast loading
df <- reactiveVal(readRDS("data/volcano.rds"))
lipinski_df <- reactiveVal(readRDS("data/lipinski.rds"))
umap_df <- reactiveVal(readRDS("data/umap.rds"))
compound_images <- readRDS("data/compound_images.rds")
tsne_data <- readRDS("data/tsne_data.rds")
compound_counts <- readRDS("data/compound_counts_per_gene.rds")
# === UI ===
ui <- fluidPage(
  titlePanel(HTML("<h2 style='text-align: center;'>Explore a Chemical-Genetic Interactions Profile of <i>B. cenocepacia</i> K56-2</h2>")),
  
  tabsetPanel(
    # --- Tab 1: Compound to Mutant ---
    tabPanel("Compound → Mutant",
             wellPanel(
               tags$div(
                 tags$p(
                   "This app accesses data from: Rahman ASMZ et al. 202X. XXX. The complete dataset contains over 3 million chemical-genetic interactions generated from ",
                   tags$a(href = "https://www.cell.com/cell-reports/fulltext/S2211-1247(24)01318-4", 
                          target = "_blank", 
                          "CIMPLE-Seq"),
                   " screening of an essential gene knockdown mutant library against 5000 compounds."
                 ),
                 tags$p("For clarity and visualization purposes, only the compounds exhibiting strong interactions (Z-Score ≤ -5) with at least one essential gene knockdown mutant are displayed here."),
                 tags$p("Special Note: t-SNE plot here projects the high-dimensional (609 knockdown mutants × mutant response) CGIPs for each compound onto two dimensions. The analysis concatenates the features (knockdown mutants and their response to each compound) into a single list for each compound, ensuring all profiles are padded to the same length. t-SNE is applied to reduce the high-dimensional data to two dimensions, creating a 2D representation of each compound profile. The 2D t-SNE results are then used as input for a K-means clustering algorithm. K-means partitioned the data into a predefined number of clusters (in this case, 10), grouping compounds based on the similarity of their 2D representations. The selected compound is highlighted in red and its corresponding cluster—determined via k-means clustering in the t-SNE space—is shaded with the same color to reflect its localized embedding.")
               )
             ),
             
             fluidRow(
               column(3, selectInput("compound", "Select Compound:", choices = NULL)),
               column(4, sliderInput("zscore", "Filter by Z-Score:", min = -25, max = 10, value = c(-25, 10), step = 0.5))
             ),
             
             fluidRow(
               column(4, tags$div(style = "text-align: center;",
                                  tags$img(src = "https://zenodo.org/records/14994629/files/heatmap_with_dendrogram_all5K.png",
                                           style = "max-width: 100%; height: auto;"))),
               column(5, tags$div(style = "text-align: center;",
                                  tags$img(src = "https://zenodo.org/records/14994629/files/Number_of_Compounds_vs_Percentage_of_Mutants_mean_240415_Compressed.png",
                                           style = "max-width: 100%; height: auto;"))),
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
    ),
    
    # --- Tab 2: t-SNE Cluster Details ---
    tabPanel("t-SNE Cluster Details",
             fluidPage(
               br(),
               fluidRow(
                 column(12,
                        h4("Other Compounds in the Same t-SNE Cluster"),
                        DTOutput("clusterCompoundTable"))
               )
             )
    ),
    
    # --- Tab 3: Mutant to Compound ---
    tabPanel("Mutant → Compound",
             fluidPage(
               br(),
               fluidRow(column(4, selectInput("Gene", "Select Target Gene:", choices = NULL))),
               fluidRow(column(12, uiOutput("compoundSectionTitle"))),
               uiOutput("compoundCountText"),
               fluidRow(column(12, plotlyOutput("geneCompoundPlot", height = "450px"))),
               fluidRow(column(12, DTOutput("mutantCompoundTable")))
             )
    )
  )
)



  # === Server ===
  server <- function(input, output, session) {
  
  observe({
    updateSelectInput(session, "compound", choices = unique(df()$Compounds))
    updateSelectInput(session, "Gene", choices = unique(df()$Gene))  # Capital 'M' to match UI
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

    clickedCompound <- reactiveVal(NULL)
    
    observeEvent(event_data("plotly_click", source = "genePlot"), {
      clickedCompound(event_data("plotly_click", source = "genePlot")$key)
    })
    
    
    output$clickedCompoundImage <- renderUI({
      req(clickedCompound())
      img_data <- compound_images[[clickedCompound()]]
      if (!is.null(img_data)) {
        tags$img(src = img_data, height = "200px", style = "border:1px solid #ccc; padding:5px;")
      } else {
        tags$p("Structure not available.")
      }
    })
    
    
    # Global reactive
    gene_affected_compounds <- reactive({
      req(input$Gene)
      
      df() %>%
        filter(Gene == input$Gene, Z.Score_ESet_mean <= -5) %>%
        arrange(Compounds, Z.Score_ESet_mean) %>%
        group_by(Compounds) %>%
        slice_min(Z.Score_ESet_mean, n = 1) %>%
        ungroup() %>%
        distinct(Compounds, Mutants, .keep_all = TRUE) %>%  # ensure uniqueness
        mutate(`Z-Score` = round(Z.Score_ESet_mean, 2)) %>%
        select(Compound = Compounds,
               Mutant = Mutants,
               `Z-Score`,
               Gene,
               COG = COG2)
    })
    
    output$mutantCompoundTable <- renderDT({
      datatable(
        gene_affected_compounds(),
        selection = 'single',
        options = list(
          pageLength = 5,
          columnDefs = list(list(className = 'dt-center', targets = "_all"))  # Center headers
        ),
        rownames = FALSE,
        class = 'cell-border stripe'
      ) %>% formatStyle(
        columns = colnames(gene_affected_compounds()),
        `text-align` = 'center'  # Center cell values
      )
    })
    
    gene_compound_plot_data <- reactive({
      req(input$Gene)
      
      summary_df <- compound_gene_summary %>% filter(Gene == input$Gene)
      stats_df <- compound_gene_stats %>% filter(Gene == input$Gene)
      
      list(
        data = summary_df,
        affected_count = stats_df$AffectedCompounds,
        total_compounds = stats_df$TotalCompounds,
        percent = stats_df$PercentAffected
      )
    })
    
    output$geneCompoundPlot <- renderPlotly({
      info <- gene_compound_plot_data()
      data <- info$data
      label <- paste0(info$affected_count, " compounds affect ", input$Gene, " (Z ≤ -5)")
      
      p <- ggplot(data, aes(x = reorder(Compound, Z.Score), y = Z.Score, key = Compounds,  # <== for plotly click
                            text = paste("Compound:", Compound, "<br>Mutant:", Mutant, "<br>Z-Score:", Z.Score))) +
        geom_point(color = "steelblue", size = 2.5, alpha = 0.85) +
        #geom_hline(aes(yintercept = -5, linetype = ""), color = "red") +
        annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, label = label, size = 5, fontface = "bold") +
        labs(x = NULL, y = "Z-Score", linetype = "Z-Score = -5") +
        theme_minimal() +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              panel.grid = element_blank(),
              axis.title = element_text(size = 14),
              legend.position = "top")
      
      ggplotly(p, tooltip = "text")
    })
    
    observeEvent(input$mutantCompoundTable_rows_selected, {
      selected_row <- input$mutantCompoundTable_rows_selected
      if (!is.null(selected_row)) {
        selected_compound <- gene_affected_compounds()$Compound[selected_row]
        updateSelectInput(session, "compound", selected = selected_compound)
      }
    })
    
    output$compoundSectionTitle <- renderUI({
      req(input$Gene)
      HTML(paste0("<h3>Chemical-Genetic Interaction Profile of <i>", input$Gene, "</i> </h3>"))
    })
    
    
    
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

    affected_compound_count <- reactive({
      req(input$Gene)
      count <- compound_counts %>%
        filter(Gene == input$Gene) %>%
        pull(AffectedCompounds)
      
      if (length(count) == 0) count <- 0
      count
    })
    
    gene_compound_plot_data <- reactive({
      req(input$Gene)
      
      # Use all interactions for the selected gene
      data <- df() %>% 
        filter(Gene == input$Gene) %>%
        group_by(Compounds) %>%
        slice_min(Z.Score_ESet_mean, n = 1) %>%
        ungroup() %>%
        mutate(Z.Score = round(Z.Score_ESet_mean, 2))
      
      # Count how many compounds are significant (Z ≤ -5)
      total_compounds <- n_distinct(df()$Compounds)
      affected_count <- n_distinct(data$Compounds[data$Z.Score <= -5])
      percent <- round((affected_count / total_compounds) * 100, 2)
      
      list(
        data = data,
        affected_count = affected_count,
        total_compounds = total_compounds,
        percent = percent
      )
    })
    
    output$geneCompoundPlot <- renderPlotly({
      info <- gene_compound_plot_data()
      data <- info$data
      label <- paste0(info$affected_count, " compounds affect ", input$Gene, " (Z ≤ -5)")
      
      p <- ggplot(data, aes(x = reorder(Compounds, Z.Score), y = Z.Score,
                            text = paste("Compound:", Compounds, "<br>Mutant:", Mutants, "<br>Z-Score:", Z.Score))) +
        geom_point(color = "darkblue", size = 1, alpha = 0.1) +
        geom_hline(aes(yintercept = -5, linetype = "Z-Score = -5"), color = "red", linewidth = 0.2, linetype = "dashed") +
        annotate("text", x = Inf, y = Inf,
                 label = "<span style='color:red;'>---</span> Z-Score = -5",
                 hjust = 0.5, vjust = 2, size = 5, family = "sans", parse = FALSE)+
        annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
                 label = label, size = 5, fontface = "bold") +
        labs(x = NULL, y = "Z-Score", linetype = "") +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid = element_blank(),
          axis.title = element_text(size = 14),
          plot.caption = ggtext::element_markdown()  # <- enable styled text
        )
      
      
      ggplotly(p, tooltip = "text")
    })
    
    
    output$compoundCountText <- renderUI({
      req(input$Gene)
      count <- affected_compound_count()
      
      HTML(paste0("<h4><b>", count, "</b> compounds affect <i>", input$Gene, "</i> (Z-Score ≤ -5)</h4>"))
    })
    
    
    df <- tsne_data
    df$Compounds_clean <- trimws(tolower(df$Compounds))
    selected_compound <- trimws(tolower(input$compound))
    
    
    # Cluster
    set.seed(42)
    kmeans_result <- kmeans(df[, c("V1", "V2")], centers = 15)
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
  

  clusterCompoundData <- reactive({
    req(input$compound)
    
    df <- tsne_data
    df$Compounds_clean <- trimws(tolower(df$Compounds))
    selected_compound <- trimws(tolower(input$compound))
    
    set.seed(42)
    df$cluster <- as.factor(kmeans(df[, c("V1", "V2")], centers = 15)$cluster)
    selected_cluster <- df %>% filter(Compounds_clean == selected_compound) %>% pull(cluster) %>% .[1]
    
    df %>%
      filter(cluster == selected_cluster) %>%
      distinct(Compounds, .keep_all = TRUE) %>%  # <== This ensures uniqueness
      select(Compound = Compounds, Cluster = cluster)
  })
  
  
  output$clusterCompoundTable <- renderDT({
    clusterCompoundData() %>%
      datatable(options = list(pageLength = 10), rownames = FALSE)
  })
  
    
}

  
# === Launch the App ===
shinyApp(ui = ui, server = server)
