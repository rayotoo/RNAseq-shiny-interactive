library(shiny)
library(ggplot2)
library(DESeq2)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggrepel)
library(DT)
library(shinyjs)
library(shinydashboard)
library(GO.db) # added GO.db

# Load required data
countdata <- read.table("final_counts.txt", header = TRUE, skip = 1, row.names = 1)
colnames(countdata) <- substr(colnames(countdata), 1, 17)
colnames(countdata) <- gsub("^\\.\\.", "", colnames(countdata))
countdata <- countdata[, -c(1:5)]

metadata <- read.delim("metadata_2.txt", row.names = 1)
metadata$sampleid <- row.names(metadata)
metadata <- metadata[match(colnames(countdata), metadata$sampleid), ]

# Remove outlier
metadata <- metadata[rownames(metadata) != 'JvRad_Rhipp_204', ]
countdata <- countdata[, !names(countdata) %in% 'JvRad_Rhipp_204']

# UI and Server setup
ui <- dashboardPage(
  dashboardHeader(title = "RNA-seq Explorer"),
  dashboardSidebar(
    width = 250,
    sliderInput("pval", "Adjusted p-value threshold:", min = 0, max = 0.1, value = 0.05, step = 0.01),
    sliderInput("logfc", "Absolute log2 Fold Change threshold:", min = 0, max = 5, value = 1, step = 0.1),
    numericInput("go_n", "Number of GO terms to display:", min = 1, max = 30, value = 10),
    actionButton("runAnalysis", "Run Analysis"),
    br(),
    tags$div(id = "analysis_message", ""),
    downloadButton("downloadFilteredData", "Filtered Data"),
    downloadButton("downloadVolcanoPlot", "Volcano Plot"),
    downloadButton("downloadHeatmapPlot", "Heatmap Plot"),
    downloadButton("downloadPCAPlot", "PCA Plot"),
    downloadButton("downloadGOPlot_BP", "GO BP Plot"),
    downloadButton("downloadGOPlot_MF", "GO MF Plot"),
    downloadButton("downloadGOPlot_CC", "GO CC Plot"),
    
    # Developer Info in Sidebar
    tags$div(
      style = "padding: 10px; font-size: 0.9em; color: #444;",
      tags$h4("Developed by", style = "font-weight: bold; margin-bottom: 5px;"),
      tags$p("Raymond Anan Otoo"),
      tags$p("Email: ",
             tags$a(href = "mailto:rotoo@omicsanalyticsgroup.com", "rotoo@omicsanalyticsgroup.com", style = "color: #007bff;", target = "_blank")
      ),
      tags$p(
        tags$a(href = "https://github.com/rayotoo", "GitHub Profile", target = "_blank", style = "color: #007bff; text-decoration: none;"),
        br(),
        tags$a(href = "https://www.linkedin.com/in/raymondotoo", "LinkedIn Profile", target = "_blank", style = "color: #007bff; text-decoration: none;")
      )
    )
    
  ),
  dashboardBody(
    useShinyjs(),
    tabsetPanel(
      id = "mainTabset",
      tabPanel("Volcano Plot", plotOutput("volcanoPlot")),
      tabPanel("Heatmap", plotOutput("heatmapPlot")),
      tabPanel("PCA Plot", plotOutput("pcaPlot")),
      tabPanel("GO Analysis",
               tabsetPanel(
                 tabPanel("Biological Process",
                          selectInput("bp_display_type", "Display Type", choices = c("Table", "Plot")),
                          conditionalPanel(
                            condition = "input.bp_display_type == 'Table'",
                            dataTableOutput("goTable_BP")
                          ),
                          conditionalPanel(
                            condition = "input.bp_display_type == 'Plot'",
                            plotOutput("goPlot_BP")
                          )
                 ),
                 tabPanel("Molecular Function",
                          selectInput("mf_display_type", "Display Type", choices = c("Table", "Plot")),
                          conditionalPanel(
                            condition = "input.mf_display_type == 'Table'",
                            dataTableOutput("goTable_MF")
                          ),
                          conditionalPanel(
                            condition = "input.mf_display_type == 'Plot'",
                            plotOutput("goPlot_MF")
                          )
                 ),
                 tabPanel("Cellular Component",
                          selectInput("cc_display_type", "Display Type", choices = c("Table", "Plot")),
                          conditionalPanel(
                            condition = "input.cc_display_type == 'Table'",
                            dataTableOutput("goTable_CC")
                          ),
                          conditionalPanel(
                            condition = "input.cc_display_type == 'Plot'",
                            plotOutput("goPlot_CC")
                          )
                 )
               )
      )
    )
  )
)

server <- function(input, output, session) {
  # Initially disable tabs
  shinyjs::disable(selector = paste0("#mainTabset li:nth-child(", 2:5, ") a"))
  
  # Hide all download buttons initially
  hideAllDownloads <- function() {
    shinyjs::hide("downloadFilteredData")
    shinyjs::hide("downloadVolcanoPlot")
    shinyjs::hide("downloadHeatmapPlot")
    shinyjs::hide("downloadPCAPlot")
    shinyjs::hide("downloadGOPlot_BP")
    shinyjs::hide("downloadGOPlot_MF")
    shinyjs::hide("downloadGOPlot_CC")
  }
  hideAllDownloads()
  
  # Function to enable tabs
  enableTabs <- function() {
    shinyjs::enable(selector = paste0("#mainTabset li:nth-child(", 2:5, ") a"))
  }
  
  # Reactive expression for filtered data
  filtered <- reactive({
    req(analysis_results()$res_df)
    analysis_results()$res_df %>% filter(padj < input$pval, abs(log2FoldChange) > input$logfc)
  })
  
  # Store the results of the analysis
  analysis_results <- reactiveVal(NULL)
  
  observeEvent(input$runAnalysis, {
    # Input validation
    if (is.null(input$pval) || is.null(input$logfc)) {
      showModal(modalDialog(
        title = "Error",
        "Please provide values for p-value and log2 fold change thresholds.",
        easyClose = TRUE
      ))
      return()
    }
    
    # Create a progress bar
    progress <- Progress$new(session, min = 0, max = 100)
    on.exit(progress$close())
    
    # Update progress and message
    updateProgress <- function(value = NULL, message = NULL, detail = NULL) {
      if (!is.null(value)) progress$set(value = value)
      if (!is.null(message)) progress$set(message = message, detail = detail)
    }
    
    # Initial message
    updateProgress(message = "Running analysis...", detail = "Initializing")
    
    # DESeq2 Analysis
    updateProgress(value = 10, message = "Running DESeq2...", detail = "Calculating")
    tryCatch({
      ddsMat <- DESeqDataSetFromMatrix(countData = countdata, colData = metadata, design = ~Group)
      ddsMat <- DESeq(ddsMat, quiet = TRUE)
      res <- results(ddsMat)
      res_df <- as.data.frame(res)
      res_df$gene <- rownames(res_df)
      res_df <- na.omit(res_df)
    }, error = function(e) {
      showModal(modalDialog(
        title = "Error",
        paste("Error during DESeq2 analysis:", e$message),
        easyClose = TRUE
      ))
      return(NULL)
    })
    
    if (is.null(res_df)) return()
    
    # Convert to rlog for PCA and Heatmap
    updateProgress(value = 40, message = "Calculating rlog transformation...", detail = "Transforming data")
    tryCatch({
      ddsMat_rlog <- rlog(ddsMat, blind = FALSE)
    }, error = function(e) {
      showModal(modalDialog(
        title = "Error",
        paste("Error during rlog transformation:", e$message),
        easyClose = TRUE
      ))
      return(NULL)
    })
    
    if (is.null(ddsMat_rlog)) return()
    
    # Store results
    analysis_results(list(
      res_df = res_df,
      ddsMat_rlog = ddsMat_rlog
    ))
    
    # GO Analysis
    updateProgress(value = 60, message = "Performing GO analysis...", detail = "Calculating")
    
    significant_genes <- res_df %>%
      filter(padj < input$pval)
    
    if (!"entrez" %in% colnames(significant_genes)) {
      tryCatch({
        significant_genes$entrez <- mapIds(
          x = org.Mm.eg.db,
          keys = rownames(significant_genes),
          column = "ENTREZID",
          keytype = "SYMBOL",
          multiVals = "first"
        )
      }, error = function(e) {
        showModal(modalDialog(
          title = "Error",
          paste("Error during Entrez ID mapping:", e$message),
          easyClose = TRUE
        ))
        return(NULL)
      })
    }
    
    if (is.null(significant_genes$entrez)) return()
    
    significant_genes_entrez <- subset(significant_genes, is.na(entrez) == FALSE)
    
    gene_matrix <- significant_genes_entrez$log2FoldChange
    names(gene_matrix) <- significant_genes_entrez$entrez
    
    tryCatch({
      go_analysis_results <- list(
        BP = enrichGO(
          gene = names(gene_matrix),
          OrgDb = org.Mm.eg.db,
          readable = TRUE,
          ont = "BP",
          pvalueCutoff = 0.05,
          qvalueCutoff = 0.10
        ),
        MF = enrichGO(
          gene = names(gene_matrix),
          OrgDb = org.Mm.eg.db,
          readable = TRUE,
          ont = "MF",
          pvalueCutoff = 0.05,
          qvalueCutoff = 0.10
        ),
        CC = enrichGO(
          gene = names(gene_matrix),
          OrgDb = org.Mm.eg.db,
          readable = TRUE,
          ont = "CC",
          pvalueCutoff = 0.05,
          qvalueCutoff = 0.10
        )
      )
      analysis_results(c(analysis_results(), list(go_analysis_results = go_analysis_results)))
      
    }, error = function(e) {
      showModal(modalDialog(
        title = "Error",
        paste("Error during GO analysis:", e$message),
        easyClose = TRUE
      ))
      return(NULL)
    })
    
    if (is.null(analysis_results()$go_analysis_results)) return()
    
    # Update progress
    updateProgress(value = 100, message = "Analysis complete!", detail = "Rendering plots")
    
    # Enable the other tabs
    enableTabs()
    showAllDownloads()
    
    # Show a message
    output$analysis_message <- renderText({
      "Analysis Complete!  Please check the other tabs."
    })
    
  })
  
  # Function to show all download buttons using shinyjs
  showAllDownloads <- function() {
    shinyjs::show("downloadFilteredData")
    shinyjs::show("downloadVolcanoPlot")
    shinyjs::show("downloadHeatmapPlot")
    shinyjs::show("downloadPCAPlot")
    shinyjs::show("downloadGOPlot_BP")
    shinyjs::show("downloadGOPlot_MF")
    shinyjs::show("downloadGOPlot_CC")
  }
  
  # Render plots and tables
  output$volcanoPlot <- renderPlot({
    req(analysis_results()$res_df)
    res_df <- analysis_results()$res_df
    ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
      geom_point(aes(color = padj < input$pval & abs(log2FoldChange) > input$logfc)) +
      theme_minimal() +
      scale_color_manual(values = c("gray", "red")) +
      labs(title = "Volcano Plot", x = "log2 Fold Change", y = "-log10 Adjusted p-value")
  })
  
  output$heatmapPlot <- renderPlot({
    req(analysis_results()$res_df, analysis_results()$ddsMat_rlog)
    
    res_df <- analysis_results()$res_df
    ddsMat_rlog <- analysis_results()$ddsMat_rlog
    
    significant_genes <- res_df %>%
      filter(padj < input$pval, abs(log2FoldChange) > input$logfc) %>%
      rownames()
    
    filtered_genes <- significant_genes[!grepl("^Gm", significant_genes)]
    
    if (length(filtered_genes) > 0) {
      mat <- assay(ddsMat_rlog[filtered_genes, ])
      if (nrow(mat) > 40) {
        mat <- mat[1:40, ]
      }
      annotation_col <- data.frame(
        Group = factor(colData(ddsMat_rlog)$Group),
        Replicate = factor(colData(ddsMat_rlog)$Replicate),
        row.names = colData(ddsMat_rlog)$sampleid
      )
      
      ann_colors <- list(
        Group = c(sham = "lightblue", ir = "purple"),
        Replicate = c(Rep1 = "red", Rep2 = "green", Rep3 = "blue", Rep4 = "forestgreen", Rep5 = "black")
      )
      
      pheatmap(
        mat = mat,
        color = colorRampPalette(c("blue", "white", "red"))(255),
        scale = "row",
        annotation_col = annotation_col,
        annotation_colors = ann_colors,
        fontsize = 8,
        show_colnames = TRUE,
        main = "Heatmap of Top 40 Significant Genes"
      )
    } else {
      plot.new()
      text(
        x = 0.5,
        y = 0.5,
        "No significant genes to display with current thresholds.",
        cex = 1,
        col = "red",
        adj = c(0.5, 0.5)
      )
    }
  })
  
  
  output$pcaPlot <- renderPlot({
    req(analysis_results()$ddsMat_rlog)
    ddsMat_rlog <- analysis_results()$ddsMat_rlog
    
    pcaData <- plotPCA(ddsMat_rlog, intgroup = "Group", ntop = 500, returnData = TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    pcaData$sample <- row.names(pcaData)
    
    palette_colors <- c("red", "blue")
    
    ggplot(pcaData, aes(x = PC1, y = PC2, color = Group, label = sample)) +
      geom_point(size = 5, alpha = 0.9, shape = 16) +
      geom_text_repel(size = 5, box.padding = 0.6, max.overlaps = 15) +
      scale_color_manual(values = palette_colors) +
      labs(
        title = "Principal Component Analysis (PCA)",
        subtitle = "Top 500 Most Variable Genes",
        x = paste0("Principal Component 1 (", percentVar[1], "%)"),
        y = paste0("Principal Component 2 (", percentVar[2], "%)"),
        color = "Experimental Group"
      ) +
      theme_minimal(base_size = 18) +
      theme(
        legend.position = "top",
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(size = 14),
        plot.title = element_text(face = "bold", size = 20),
        plot.subtitle = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18)
      )
  })
  
  
  output$goTable_BP <- renderDataTable({
    req(analysis_results()$go_analysis_results)
    go_enrich_df <- as.data.frame(analysis_results()$go_analysis_results$BP)
    if (nrow(go_enrich_df) == 0) {
      data.frame(Message = "No significant GO terms with current p-value threshold.")
    } else {
      datatable(go_enrich_df, options = list(pageLength = 10))
    }
  })
  
  output$goTable_MF <- renderDataTable({
    req(analysis_results()$go_analysis_results)
    go_enrich_df <- as.data.frame(analysis_results()$go_analysis_results$MF)
    if (nrow(go_enrich_df) == 0) {
      data.frame(Message = "No significant GO terms with current p-value threshold.")
    } else {
      datatable(go_enrich_df, options = list(pageLength = 10))
    }
  })
  
  output$goTable_CC <- renderDataTable({
    req(analysis_results()$go_analysis_results)
    go_enrich_df <- as.data.frame(analysis_results()$go_analysis_results$CC)
    if (nrow(go_enrich_df) == 0) {
      data.frame(Message = "No significant GO terms with current p-value threshold.")
    } else {
      datatable(go_enrich_df, options = list(pageLength = 10))
    }
  })
  
  output$goPlot_BP <- renderPlot({
    req(analysis_results()$go_analysis_results)
    go_enrich <- analysis_results()$go_analysis_results$BP
    if (!is.null(go_enrich) && nrow(as.data.frame(go_enrich)) > 0) {
      barplot(go_enrich,
              drop = TRUE,
              showCategory = input$go_n,
              title = "GO Biological Pathways",
              font.size = 8)
    } else {
      plot.new()
      text(x = 0.5, y = 0.5, "No significant GO terms to display.",
           cex = 1, col = "red", adj = c(0.5, 0.5))
    }
  })
  
  output$goPlot_MF <- renderPlot({
    req(analysis_results()$go_analysis_results)
    go_enrich <- analysis_results()$go_analysis_results$MF
    if (!is.null(go_enrich) && nrow(as.data.frame(go_enrich)) > 0) {
      barplot(go_enrich,
              drop = TRUE,
              showCategory = input$go_n,
              title = "GO Molecular Functions",
              font.size = 8)
    } else {
      plot.new()
      text(x = 0.5, y = 0.5, "No significant GO terms to display.",
           cex = 1, col = "red", adj = c(0.5, 0.5))
    }
  })
  
  output$goPlot_CC <- renderPlot({
    req(analysis_results()$go_analysis_results)
    go_enrich <- analysis_results()$go_analysis_results$CC
    if (!is.null(go_enrich) && nrow(as.data.frame(go_enrich)) > 0) {
      barplot(go_enrich,
              drop = TRUE,
              showCategory = input$go_n,
              title = "GO Cellular Components",
              font.size = 8)
    } else {
      plot.new()
      text(x = 0.5, y = 0.5, "No significant GO terms to display.",
           cex = 1, col = "red", adj = c(0.5, 0.5))
    }
  })
  
  # Download Handlers
  output$downloadFilteredData <- downloadHandler(
    filename = "filtered_data.csv",
    content = function(file) {
      req(analysis_results()$res_df)
      write.csv(filtered(), file)
    }
  )
  
  output$downloadVolcanoPlot <- downloadHandler(
    filename = "volcano_plot.png",
    content = function(file) {
      req(analysis_results()$res_df)
      ggsave(file, plot = output$volcanoPlot, device = "png", width = 8, height = 6, units = "in", dpi = 300)
    }
  )
  
  output$downloadHeatmapPlot <- downloadHandler(
    filename = "heatmap_plot.png",
    content = function(file) {
      req(analysis_results()$res_df, analysis_results()$ddsMat_rlog)
      
      res_df <- analysis_results()$res_df
      ddsMat_rlog <- analysis_results()$ddsMat_rlog
      
      significant_genes <- res_df %>%
        filter(padj < input$pval, abs(log2FoldChange) > input$logfc) %>%
        rownames()
      
      filtered_genes <- significant_genes[!grepl("^Gm", significant_genes)]
      
      if (length(filtered_genes) > 0) {
        png(file, width = 8, height = 6, units = "in", res = 300) # Open png
        
        mat <- assay(ddsMat_rlog[filtered_genes, ])
        if (nrow(mat) > 40) {
          mat <- mat[1:40, ]
        }
        annotation_col = data.frame(
          Group = factor(colData(ddsMat_rlog)$Group),
          Replicate = factor(colData(ddsMat_rlog)$Replicate),
          row.names = colData(ddsMat_rlog)$sampleid
        )
        
        ann_colors = list(
          Group = c(sham = "lightblue", ir = "purple"),
          Replicate = c(Rep1 = "red", Rep2 = "green", Rep3 = "blue", Rep4 = "forestgreen", Rep5 = "black")
        )
        
        pheatmap(
          mat = mat,
          color = colorRampPalette(c("blue", "white", "red"))(255),
          scale = "row",
          annotation_col = annotation_col,
          annotation_colors = ann_colors,
          fontsize = 8,
          show_colnames = T,
          main = "Heatmap of Top 40 Significant Genes"
        )
        dev.off() # Close the png
      }
    }
  )
  
  output$downloadPCAPlot <- downloadHandler(
    filename = "pca_plot.png",
    content = function(file) {
      req(analysis_results()$ddsMat_rlog)
      ggsave(file, plot = output$pcaPlot, device = "png", width = 8, height = 6, units = "in", dpi = 300)
    }
  )
  
  output$downloadGOPlot_BP <- downloadHandler(
    filename = "go_bp_plot.png",
    content = function(file) {
      req(analysis_results()$go_analysis_results)
      ggsave(file, plot = output$goPlot_BP(), device = "png", width = 8, height = 6, units = "in", dpi = 300)
    }
  )
  
  output$downloadGOPlot_MF <- downloadHandler(
    filename = "go_mf_plot.png",
    content = function(file) {
      req(analysis_results()$go_analysis_results)
      ggsave(file, plot = output$goPlot_MF(), device = "png", width = 8, height = 6, units = "in", dpi = 300)
    }
  )
  
  output$downloadGOPlot_CC <- downloadHandler(
    filename = "go_cc_plot.png",
    content = function(file) {
      req(analysis_results()$go_analysis_results)
      ggsave(file, plot = output$goPlot_CC(), device = "png", width = 8, height = 6, units = "in", dpi = 300)
    }
  )
}

shinyApp(ui = ui, server = server)
