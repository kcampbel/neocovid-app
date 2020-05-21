library(shiny)
library(shinydashboard)
library(DT)
library(data.table)
library(tidyverse)

##### ~~~~~ UI ~~~~~ #####
ui <- dashboardPage(skin = "purple",
  dashboardHeader(title = "neoCOVID Explorer"),
  dashboardSidebar(
    width = 200,
    sidebarMenu(id = "sidebar", 
      menuItem("Virus", tabName = "virus", icon = icon("ellipsis-v")),
      div( id = 'sidebar_virus',
           conditionalPanel("input.sidebar == 'virus'",
                            selectInput("proteinInput", "Select a Protein", choices = c("Show All"), selected = "Show All")
                            )
           ),
      
      menuItem("Epitope", tabName = "epitope", icon = icon("flag")),
      div( id = 'sidebar_epitope',
           conditionalPanel("input.sidebar == 'epitope'"
                            # selectizeInput("epitopeInput", "Epitope of Interest", choices = NULL, selected = NULL, multiple = TRUE),
                            
                            )
           )
      ),
      
      menuItem("HLA", tabName = "hla",icon = icon("fingerprint")),
      
      menuItem("Population", tabName = "populations", icon = icon("globe"))
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "virus",
              box(title = "Viral View", width = 12,
                  plotlyOutput("mainPlotly")
                  ),
              box(title = "Epitopes", width = 12,
                  dataTableOutput("peptideDataTable")
                  )
              ),
      tabItem(tabName = "epitope",
              fluidRow(
                box(h1("Explore SARS-CoV-2 Epitopes"),
                    br(),
                    em("If your epitope is not found, the closest string match will be displayed."),
                    br(),
                    textInput("epitopeInput", "Find Epitope", value = ""),
                    width = 6, height = 300),
                box(width = 6, height = 300)
              ),
              box(title = "Filtered Epitope", width = 12,
                  dataTableOutput("epitopeDataTable")
              )
      )
    )
  )
)

##### ~~~~~ SERVER ~~~~~ #####
# Define server logic required to draw a histogram ----
server <- function(input, output, session) {
  load('~/neocovid-app/app/data/forApp.Rda')
  
  ##### Viral View #####
  # Get protein name options
  observe({
    updateSelectInput(
      session,
      "proteinInput",
      choices = c("Show All", unique(as.character(proteinSequences3$Protein))),
      selected = "Show All"
    )
  })
  
  # Subset to selected protein
  subsetProteins <- reactive({
    if(input$proteinInput == "Show All"){
      subsetProteins <- proteinSequences3
    } else {
      subsetProteins <- filter(proteinSequences3, Protein == input$proteinInput)
    }
    return(subsetProteins)
  })
  
  # Output spatial plot
  output$mainPlotly <- renderPlotly({
    plotData <- subsetProteins()
    plotEpisData <- plotEpis
    
    basePlot <- ggplot(data = plotData, aes(fill = Gene)) +
      geom_hline(yintercept = -1, colour = 'grey80', size = 2) +
      # Proteins
      geom_rect(aes(xmin = polar_start, xmax = polar_end, ymin = -1.5, ymax = -0.5)) +
      scale_fill_manual(values = sample(colors, 12)) +
      # Chains
      geom_rect(aes(xmin = chain_polar_start, xmax = chain_polar_end, ymin = -1.75, ymax = -0.25), fill = 'grey80', alpha = 0.3) +
      # Peptides
      geom_rect(data = plotEpisData, aes(xmin = peptide.polar_start, xmax = peptide.polar_end, ymin = y_new, ymax = y_new + 0.5)) +
      #
      scale_x_continuous(limits = c(min(plotData$polar_start), max(plotData$polar_end))) + 
      scale_y_continuous(limits = c(-8,max(plotEpis$y_new)+1)) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    #
    ggplotly(basePlot, source = "virusPlot") %>% layout(xaxis = list(rangeslider = list(type = "date")))
  })
  
  # Subset predictions to those in view (whether it's full, a single protein, or zoomed in)
  filterPredictions <- reactive({
    zoom <- event_data("plotly_relayout", "virusPlot")
    # if plot just rendered, event_data is NULL
    # if user double clicks for autozoom, then zoom$xaxis.autorange is TRUE
    # if user resizes page, then zoom$width is pixels of plot width
    if(input$proteinInput == "Show All" & (is.null(zoom) || names(zoom[1]) %in% c("xaxis.autorange", "width"))) {
      # xlim <- "default of plot"
      return(summPredictions)
    } else if(input$proteinInput != "Show All") {
      return(filter(summPredictions, Protein == input$proteinInput))
    } else {
      xmin <- zoom$`xaxis.range[0]`
      xmax <- zoom$`xaxis.range[1]`
      get_peptides <- filter(plotEpis, peptide.polar_start > xmin & peptide.polar_end < xmax)$Peptide
      return(filter(summPredictions, Peptide %in% get_peptides))
    }
  })
  # Output predictions table
  output$peptideDataTable <- renderDataTable(filterPredictions(), options = list(scrollX = TRUE, scrollY = TRUE))
  
  ##### Epitope View #####
  subsetEpitopes <- reactive({
    if(input$epitopeInput == "") {
      return(summPredictions)
    } else if(input$epitopeInput %in% unique(summPredictions$Peptide)) {
      return(filter(summPredictions, Peptide == input$epitopeInput))
    } else {
      sub_string <- grep(input$epitopeInput, unique(summPredictions$Peptide), value = TRUE)
      get_dist <- lapply(unique(summPredictions$Peptide), function(y){ adist(input$epitopeInput, y) }) %>% unlist
      closest_peptides <- unique(summPredictions$Peptide)[which(get_dist<3)]
      return(filter(summPredictions, Peptide %in% union(sub_string, closest_peptides)))
    }
  })
  
  output$epitopeDataTable <- renderDataTable(subsetEpitopes(), options = list(scrollX = TRUE, scrollY = TRUE))
  
}

##### ~~~~~ Run App ~~~~~ #####
shinyApp(ui = ui, server = server)




