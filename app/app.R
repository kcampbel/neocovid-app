library(shiny)
library(shinydashboard)
library(DT)
library(data.table)
library(plotly)
library(tidyverse)
library(leaflet)
library(wesanderson)
# 
# load('~/neocovid-app/app/data/forApp.Rda')
load(file = url("https://github.com/kcampbel/neocovid-app/raw/master/app/data/forApp.Rda"))


##### ~~~~~ UI ~~~~~ #####
ui <- dashboardPage(skin = "purple",
                    dashboardHeader(title = "neoCOVID Explorer"),
                    ##### ~~~~~ SIDEBAR ~~~~~ ######
                    dashboardSidebar(
                      tags$head(tags$style(HTML('.content-wrapper { height: 3000px !important;}'))),
                      sidebarMenu(id = "Sidebar",
                                  menuItem("Home", tabName = "home", icon = icon("searchengin")),
                                  menuItem("Putative Epitopes", tabName = "epitopes", icon = icon("flag")),
                                  menuItem("HLA", tabName = "hla", icon = icon("fingerprint")),
                                  menuItem("Populations", tabName = "populations", icon = icon("globe"))
                      )
                    ),
                    ##### ~~~~~ BODY ~~~~~ ######
                    dashboardBody(
                      tabItems(
                        tabItem(tabName = "home",
                                h1("Select Features"),
                                "Use this page to set filters for viral features (e.g. genes, proteins), epitopes, HLA types, or populations of interest. Other pages will depict filtered results specified here.",
                                br(),
                                box(title = "Virus", width = 12, collapsible = TRUE, collapsed = FALSE,
                                    div(style="display: inline-block;vertical-align:top; width: 300px;", selectInput("geneInput", "Select a Gene(s)", choices = levels(plotProteins$Gene), selected = NULL, multiple = TRUE)),
                                    div(style="display: inline-block;vertical-align:top; width: 100px;",HTML("<br>")),
                                    div(style="display: inline-block;vertical-align:top; width: 300px;", selectInput("proteinInput", "Select a Protein(s)", choices = levels(plotProteins$Protein), selected = NULL, multiple = TRUE))
                                    ),
                                box(title = "Putative Epitopes", width = 12, collapsible = TRUE, collapsed = FALSE,
                                    div(style="display: inline-block;vertical-align:top; width: 300px;", sliderInput("episizeInput", "Epitope Size", min = 8, max = 30, value = c(8,15))),
                                    div(style="display: inline-block;vertical-align:top; width: 100px;",HTML("<br>")),
                                    div(style="display: inline-block;vertical-align:top; width: 300px;", textInput("epitopeInput", "Search for Epitope", value = ""))
                                    ),
                                box(title = "HLA", width = 12, collapsible = TRUE, collapsed = FALSE,
                                    div(style="display: inline-block;vertical-align:top; width: 300px;", selectInput("hlageneInput", "Select HLA Gene(s)", choices = unique(gsub("(.+)\\*\\d+.*", "\\1", unlist(strsplit(unique(predictions$`HLA Allele`), split = "/|-")), perl = T)), selected = c('A','B','C'), multiple = TRUE)),
                                    div(style="display: inline-block;vertical-align:top; width: 300px;", selectizeInput("hlaInput", "Search HLA", choices = unlist(strsplit(unique(predictions$`HLA Allele`), split = "-|/")), selected = c('A*68:44','B*07:02','C*08:02'), multiple = TRUE)),
                                    div(style="display: inline-block;vertical-align:top; width: 100px;",HTML("<br>")),
                                    div(style="display: inline-block;vertical-align:top; width: 300px;", checkboxInput("hlafamInput", "Show peptides from HLA proteins within family"))
                                    ),
                                box(title = "Populations", width = 12, collapsible = TRUE, collapsed = FALSE,
                                    div(style="display: inline-block;vertical-align:top; width: 300px;", selectizeInput("countryInput", "Select a Country(s)", choices = unique(world@data$NAME), selected = NULL, multiple = TRUE)),
                                    div(style="display: inline-block;vertical-align:top; width: 100px;",HTML("<br>")),
                                    div(style="display: inline-block;vertical-align:top; width: 300px;", numericInput("hlafreqInput", "Minimum Population Allelic Frequency", 0, min = 0, max = 1))
                                    )
                        ),
                        tabItem(tabName = "epitopes",
                                # tags$style(type = "text/css", "#map {height: calc(100vh - 80px) !important;}"),
                                # dataTableOutput('test'),
                                box("N HLA per Predicted Epitope", width = 6, collapsible = TRUE, collapsed = FALSE, dataTableOutput("nHLA")),
                                # box("Previously Published", width = 6, collapsible = TRUE, collapsed = FALSE, dataTableOutput("publishedEpis")),
                                box("All Filtered Predicted Epitopes", width = 12, collapsible = TRUE, collapsed = FALSE, dataTableOutput("epitopeOutput"))
                        ),
                        tabItem(tabName = "hla", 
                                # tags$style(type = "text/css", "#map {height: calc(100vh - 80px) !important;}")
                        ),
                        tabItem(tabName = "populations", 
                                tags$style(type = "text/css", "#map {height: calc(100vh - 80px) !important;}"),
                                leafletOutput("worldmap"),
                                box("Population Allelic Frequencies", width = 12, collapsible = TRUE, collapsed = FALSE, dataTableOutput('populationAlleles'))
                        )
                      )
                    )
)

##### ~~~~~ SERVER ~~~~~ #####
server <- function(input, output, session) {
  
  # geneInput proteinInput episizeInput epitopeInput hlageneInput hlaInput hlafamInput countryInput hlafreqInput
  getProteins <- reactive({
    if(length(input$geneInput)>0) {
      if(length(input$proteinInput)>0) {
        filt <- plotProteins[Gene %in% input$geneInput | Protein %in% input$proteinInput]
      } else {
        filt = plotProteins[Gene %in% input$geneInput]
      }
    } else {
      if(length(input$proteinInput)>0) {
        filt <- plotProteins[Protein %in% input$proteinInput]
      } else {
        filt = plotProteins
      }
    }
    return(filt)
  })
  #
  getPopulations <- reactive({
    if(length(input$countryInput)>0){
      return(populationFrequencies[Country %in% input$countryInput & `Allele Frequency`>=input$hlafreqInput])
    } else {
      return(populationFrequencies[`Allele Frequency`>=input$hlafreqInput])
    }
  })
  #
  getEpitopes <- reactive({
    # Filter by Epitope size
    filter <- predictions[nchar(Peptide)>=input$episizeInput[1] & nchar(Peptide)<=input$episizeInput[2]]
    
    # Filter by corresponding Gene/Protein
    # locusFilter <- merge.data.table(getProteins(), sizeFilter, all.y = FALSE, by = c('Gene','Protein'))
    filter <- filter[Gene %in% getProteins()$Gene | Protein %in% getProteins()$Protein]
    
    # Filter by HLA Features
    if(length(input$hlageneInput)>0){
      filter[, tmp := ifelse(any(unlist(hlagene) %in% input$hlageneInput), 1, 0), by = c('Peptide','HLA Allele')]
      filter <- filter[tmp==1 | hlagene == input$hlageneInput]
    } 
    #
    if(!is.null(input$hlaInput)) {
      filter[, tmp := ifelse(any(unlist(hlaalleles) %in% input$hlaInput), 1, 0), by = c('Peptide','HLA Allele')]
      if(input$hlafamInput == TRUE) {
        hlafamilies = gsub("([A-Z]*\\d*\\*\\d+):\\d+", "\\1", input$hlaInput)
        filter[, tmp2 := ifelse(any(unlist(hlafam) %in% hlafamilies), 1, 0), by = c('Peptide','HLA Allele')]
        # hlafamFilter = locusFilter[hlafam %in% hlafamilies | (grepl("-|/", `HLA Allele`) & any(unlist(hlafam) %in% hlafamilies))]
        filter <- filter[tmp==1 | tmp2 == 1]
      } else {
        filter <- filter[tmp == 1]
      }
    }
    
    # Filter by HLA features (population)
    filter[, tmp := ifelse(any(unlist(hlaalleles) %in% getPopulations()$`HLA Allele`) | any(unlist(hlafam) %in% getPopulations()$`HLA Allele`), 1, 0), by = c('Peptide','HLA Allele')]
    filter[tmp == 1]
    
    if(length(input$epitopeInput)>0){
      filter[, dist := list(lapply(input$epitopeInput, function(epi){ adist(Peptide, epi) } )), by = c('Peptide','HLA Allele')]
      filter[, epiInPep := ifelse(any(grepl(paste0(input$epitopeInput, collapse = "|"), Peptide)), 1, 0), by = c('Peptide','HLA Allele')]
      filter[, pepInEpi := ifelse(any(grepl(Peptide, input$epitopeInput)), 1, 0), by = c('Peptide','HLA Allele')]
      filter <- filter[any(dist<3) | epiInPep == 1 | pepInEpi == 1]
    }
    return(filter)
    })
  
  ##### ~~~~~ Epitopes ~~~~~ #####
  output$epitopeOutput <- renderDataTable(getEpitopes(), options = list(scrollX = TRUE, scrollY = TRUE))
  
  # N HLA per Epitope
  nHLA <- reactive({
    nHLA <- getEpitopes()[, .(`N HLA Alleles` = .N, `HLA Alleles` = list(`HLA Allele`)), by = c('Peptide', 'Gene','peptide.polar_start')]
    nHLA <- nHLA[order(-`N HLA Alleles`)]
    
    return(nHLA)
  })
  output$nHLA <- renderDataTable(nHLA()[,c('Peptide','N HLA Alleles', 'HLA Alleles')], options = list(scrollX = TRUE, scrollY = TRUE))
  #
  # findPublished <- reactive({
  #   publishedEpis[, index := 1:.N]
  #   # publishedEpis[, calcdist := min(unlist(lapply(nHLA()$Peptide, function(epi){ adist(Peptide, epi) } ))), by = c('index')]
  #   publishedEpis[, epiInPep := ifelse(any(grepl(paste0(nHLA()$Peptide, collapse = "|"), Peptide)), 1, 0), by = c('index')]
  #   publishedEpis[, pepInEpi := ifelse(any(grepl(Peptide, nHLA()$Peptide)), 1, 0), by = c('index')]
  #   publishedEpis <- filter[epiInPep == 1 | pepInEpi == 1]
  #   return(publishedEpis[,c('Peptide','Protein','HLAA Restriction','Source','Label')])
  # })
  # output$publishedEpis <- renderDataTable(findPublished(), options = list(scrollX = TRUE, scrollY = TRUE))
  
  # output$epitopesPlot <- renderPlotly({
  #   
  #   base <- ggplot(data = nHLA(), aes(fill = Gen))
  # })
  
  ##### ~~~~~ HLA ~~~~~ #####
  
  ##### ~~~~~ POPULATIONS ~~~~~ #####
  autoZoomMap <- reactive({
    if(length(input$countryInput)>0){
      world@data %>% filter(NAME %in% input$countryInput)
    } else {
      world@data %>% filter(NAME == "United States")
    }
  })
  #
  output$worldmap <- renderLeaflet({
    leaflet(world) %>%
      addProviderTiles(providers$CartoDB.Positron,
                       options = providerTileOptions(noWrap = TRUE)) %>%
      addTiles() %>% 
      # setView(42, 16, 4)
      setView(lng = autoZoomMap()$LON[1], lat = autoZoomMap()$LAT[1], zoom = 3) #%>%
      # addCircleMarkers(~LON, ~LAT,
      #                  fillColor = ~mypalette(PICK),
      #                  fillOpacity = 0.7,
      #                  color = "white",
      #                  radius = 8,
      #                  stroke = FALSE,
      #                  label = mytext,
      #                  labelOptions = labelOptions( style = list("font-weight" = "normal",
      #                                                            padding = "3px 8px"),
      #                                               textsize = "13px",
      #                                               direction = "auto")
      # )
  })
  output$populationAlleles <- renderDataTable(getPopulations(), options = list(scrollX = TRUE, scrollY = TRUE))
}

##### ~~~~~ Run App ~~~~~ #####
shinyApp(ui = ui, server = server)




