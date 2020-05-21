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
      menuItem("Virus", tabName = "virus"),
      div( id = 'sidebar_virus',
           conditionalPanel("input.sidebar == 'virus'",
                            selectInput("proteinInput", "Select a Protein", choices = c("Show All"), selected = "Show All")
                            )
           ),
      
      menuItem("Epitope", tabName = "epitope"),
      
      menuItem("HLA", tabName = "hla"),
      
      menuItem("Population", tabName = "populations", icon = icon("globe"))
    )
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
              )
    )
  )
)

##### ~~~~~ SERVER ~~~~~ #####
# Define server logic required to draw a histogram ----
server <- function(input, output, session) {
  load('~/neocovid-app/app/data/forApp.Rda')
  output$test <- renderDataTable(predictions)
  
  observe({
    updateSelectInput(
      session,
      "proteinInput",
      choices = c("Show All", unique(as.character(proteinSequences3$Protein))),
      selected = "Show All"
    )
  })
  
  subsetProteins <- reactive({
    if(input$proteinInput == "Show All"){
      subsetProteins <- proteinSequences3
    } else {
      subsetProteins <- filter(proteinSequences3, Protein == input$proteinInput)
    }
    return(subsetProteins)
  })
  
  
  subsetEpis <- reactive({
    if(input$proteinInput == "Show All"){
      subsetEpis <- plotEpis
    } else {
      subsetEpis <- filter(plotEpis, Protein == input$proteinInput)
    }
    return(subsetEpis)
  })
  output$peptideDataTable <- renderDataTable(subsetEpis(), options = list(scrollX = TRUE, scrollY = TRUE))

  output$mainPlotly <- renderPlotly({
    plotData <- subsetProteins()
    plotEpisData <- subsetEpis()
    
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
      scale_y_continuous(limits = c(-8,max(plotEpis$y_new)+1))
    #
    ggplotly(basePlot) %>% layout(xaxis = list(rangeslider = list(type = "date")))
  })
  
}

##### ~~~~~ Run App ~~~~~ #####
shinyApp(ui = ui, server = server)




