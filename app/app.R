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
                                box(width = 12, plotlyOutput("epitopesPlot")),
                                box("N HLA per Predicted Epitope", width = 4, collapsible = TRUE, collapsed = FALSE, dataTableOutput("nHLA")),
                                box("Previously Published", width = 8, collapsible = TRUE, collapsed = FALSE, dataTableOutput("publishedEpis")),
                                box("All Filtered Predicted Epitopes", width = 12, collapsible = TRUE, collapsed = FALSE, dataTableOutput("epitopeOutput"))
                        ),
                        tabItem(tabName = "hla",
                                box(width = 12, plotlyOutput("hlaPlot")),
                                box("Top 10 HLA Types with the most predicted epitopes", width = 6, plotlyOutput('nEpiPlot')),
                                box("N Predicted Epitopes per HLA Allele", width = 6, dataTableOutput("nEpi"))
                                # tags$style(type = "text/css", "#map {height: calc(100vh - 80px) !important;}")
                        ),
                        tabItem(tabName = "populations", 
                                tags$style(type = "text/css", "#map {height: calc(100vh - 80px) !important;}"),
                                box(leafletOutput("worldmap"), width = 6),
                                box("Cover Set Solution", width = 6, plotlyOutput("setCoverPlot")), 
                                box("Population Allelic Frequencies", width = 8, collapsible = TRUE, collapsed = FALSE, dataTableOutput('populationAlleles')),
                                box("Cover Sets", width = 4, dataTableOutput("coverSetData"))
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
      # return(populationFrequencies[Country %in% input$countryInput & `Allele Frequency`>=input$hlafreqInput])
      return(populationFrequencies[.(input$countryInput)][`Allele Frequency`>=input$hlafreqInput])
    } else {
      return(populationFrequencies[`Allele Frequency`>=input$hlafreqInput])
    }
  })
  #
  getEpitopes <- reactive({
    # Filter by Epitope size
    filter <- predictions[nchar(Peptide)>=input$episizeInput[1] & nchar(Peptide)<=input$episizeInput[2]]
    
    # Filter by corresponding Gene/Protein
    filter <- filter[Gene %in% getProteins()$Gene | Protein %in% getProteins()$Protein]
    
    # Filter by HLA Features
    if(length(input$hlageneInput)>0){
      setkey(filter, hlagene)
      filter = filter[.(input$hlageneInput)]
    } 
    #
    if(length(input$hlaInput)>0) {
      setkey(filter, hlaalleles)
      filter = filter[.(input$hlaInput)]
      if(input$hlafamInput == TRUE) {
        setkey(filter, hlafam)
        hlafamilies = gsub("([A-Z]*\\d*\\*\\d+):\\d+", "\\1", input$hlaInput)
        filter = filter[.(hlafamilies)]
      }
    }
    
    # Filter by HLA features (population)
    filter = filter[hlaalleles %in% getPopulations()$`HLA Allele` | hlafam %in% getPopulations()$`HLA Allele`]

    if(length(input$epitopeInput)>0){
      filter[, dist := list(lapply(input$epitopeInput, function(epi){ adist(Peptide, epi) } )), by = c('Peptide','HLA Allele')]
      filter[, epiInPep := ifelse(any(grepl(paste0(input$epitopeInput, collapse = "|"), Peptide)), 1, 0), by = c('Peptide','HLA Allele')]
      filter[, pepInEpi := ifelse(any(grepl(Peptide, input$epitopeInput)), 1, 0), by = c('Peptide','HLA Allele')]
      filter <- filter[any(dist<3) | epiInPep == 1 | pepInEpi == 1]
    }
    return(filter)
    })
  
  ##### ~~~~~ Epitopes ~~~~~ #####
  output$epitopeOutput <- renderDataTable(getEpitopes()[,c('Peptide','Gene','Protein','Position','HLA Allele','Class',
                                                           'Median Affinity (nM)','Best Score (nM)','Predicted Stability (NetMHCstabpan)',
                                                           'Half Life','Stability Rank')],
                                          extensions = 'Buttons', 
                                          options = list(scrollX = TRUE, scrollY = TRUE,
                                                         dom = 'Blfrtip',
                                                         buttons = c('copy', 'csv')))
  
  # N HLA per Epitope
  nHLA <- reactive({
    nHLA <- getEpitopes()[, .(`N HLA Alleles` = .N, `HLA Alleles` = list(`HLA Allele`)), by = c('Peptide', 'Gene','peptide.polar_start')]
    nHLA <- nHLA[order(-`N HLA Alleles`)]
    return(nHLA)
  })
  output$nHLA <- renderDataTable(nHLA()[,c('Peptide','N HLA Alleles', 'HLA Alleles')],
                                 extensions = 'Buttons', 
                                 options = list(scrollX = TRUE, scrollY = TRUE,
                                                dom = 'Blfrtip',
                                                buttons = c('copy', 'csv')))
  #
  findPublished <- reactive({
    setkey(comparePublished, Peptide)
    comparePublished[.(nHLA()$Peptide)][!is.na(Comparison)]
    return(comparePublished)
  })
  output$publishedEpis <- renderDataTable(findPublished()[,c('Peptide','Published Epitope','Comparison','HLA Restriction','Type','Source')],
                                          extensions = 'Buttons', 
                                          options = list(scrollX = TRUE, scrollY = TRUE,
                                                         dom = 'Blfrtip',
                                                         buttons = c('copy', 'csv')))
  #
  output$epitopesPlot <- renderPlotly({
    plotProteins[, fillGene := ifelse(Gene %in% getProteins()$Gene, Gene, NA)]
    plotProteins$fillGene = factor(plotProteins$fillGene, levels = levels(plotProteins$Gene))
    forPlot <- nHLA()#merge.data.table(nHLA(), getEpitopes(), by = 'Peptide')
    forPlot$Gene = factor(forPlot$Gene, levels = levels(plotProteins$Gene))
    #
    mycolors <- sample(colors, 12)
    base <- ggplot(data = plotProteins, aes(fill = Gene)) +
      geom_hline(yintercept = -1, colour = 'grey80', size = 2) +
      # Proteins
      geom_rect(aes(xmin = polar_start, xmax = polar_end, ymin = -1.5, ymax = -0.5)) +
      scale_fill_manual(values = mycolors, na.value = 'grey80') +
      # Chains
      geom_rect(aes(xmin = chain_polar_start, xmax = chain_polar_end, ymin = -1.75, ymax = -0.25), fill = 'grey80', alpha = 0.3) +
      # Peptides
      geom_col(data = forPlot, aes(x = peptide.polar_start, y = `N HLA Alleles`, fill = Gene)) +
      # scale_colour_manual(values = mycolors, na.value = 'grey80') +
      #
      scale_x_continuous(limits = c(min(plotProteins$polar_start), max(plotProteins$polar_end))) + 
      scale_y_continuous(limits = c(-8,max(forPlot$`N HLA Alleles`)+5)) +
      labs(y = "N HLA Alleles", x = "Postion in Viral Proteome")
    #
    ggplotly(base, source = 'epiPlot') %>%
      layout(xaxis = list(rangeslider = list(type = 'date')))
  })
  
  ##### ~~~~~ HLA ~~~~~ #####
  # N HLA per Epitope
  nEpis <- reactive({
    nEpis <- getEpitopes()[, .(`N Peptides` = .N, `Peptides` = list(`Peptide`)), by = c('HLA Allele')]
    nEpis <- nEpis[order(-`N Peptides`)]
    return(nEpis)
  })
  output$nEpi <- renderDataTable(nEpis()[,c('HLA Allele','N Peptides', 'Peptides')],
                                 extensions = 'Buttons', 
                                 options = list(scrollX = TRUE, scrollY = TRUE,
                                                dom = 'Blfrtip',
                                                buttons = c('copy', 'csv')))
  output$hlaPlot <- renderPlotly({
    plotProteins[, fillGene := ifelse(Gene %in% getProteins()$Gene, Gene, NA)]
    plotProteins$fillGene = factor(plotProteins$fillGene, levels = levels(plotProteins$Gene))
    # forPlot <- nHLA()#merge.data.table(nHLA(), getEpitopes(), by = 'Peptide')
    forPlot <- getEpitopes()
    forPlot$Gene = factor(forPlot$Gene, levels = levels(plotProteins$Gene))
    #
    mycolors <- sample(colors, 12)
    base <- ggplot(data = plotProteins, aes(fill = Gene)) +
      geom_hline(yintercept = -25, colour = 'grey80', size = 2) +
      # Proteins
      geom_rect(aes(xmin = polar_start, xmax = polar_end, ymin = -40, ymax = -10)) +
      scale_fill_manual(values = mycolors, na.value = 'grey80', guide = FALSE) +
      # Chains
      geom_rect(aes(xmin = chain_polar_start, xmax = chain_polar_end, ymin = -50, ymax = 0), fill = 'grey80', alpha = 0.3) +
      # Peptides
      # geom_col(data = forPlot, aes(x = peptide.polar_start, y = `N HLA Alleles`, fill = Gene)) +
      geom_point(data = forPlot, aes(x = peptide.polar_start, y = 500-`Median Affinity (nM)`, colour = Gene, shape = hlagene), alpha = 0.8, size = 1.5) +
      scale_color_manual(values = mycolors, na.value = 'grey80', guide = FALSE) +
      #
      scale_x_continuous(limits = c(min(plotProteins$polar_start), max(plotProteins$polar_end))) +
      scale_y_continuous(limits = c(-50,550), breaks = c(0, 100, 200, 300, 400, 500),
                         labels = c(500,400,300,200,100,0)) +
      labs(y = "Predicted Binding Affinity (nM)", x = "Position in Viral Proteome",
           colour = "HLA Gene")
    #
    ggplotly(base, source = 'hlaPlot') %>%
      layout(xaxis = list(rangeslider = list(type = 'date'))) %>%
      hide_legend()
  })
  
  output$nEpiPlot <- renderPlotly({
    forPlot <- nEpis()
    forPlot$`HLA Allele` <- factor(forPlot$`HLA Allele`, levels = forPlot$`HLA Allele`)
    base <- ggplot(data = forPlot[1:min(nrow(forPlot), 10),], 
           aes(x = `HLA Allele`, y = `N Peptides`)) +
      geom_col(fill = '#686de0', colour = 'black') +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    ggplotly(base, source = 'nEpiPlot')
  })
  
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
      setView(lng = autoZoomMap()$LON[1], lat = autoZoomMap()$LAT[1], zoom = 2) #%>%
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
  #
  output$setCoverPlot <- renderPlotly({
    setkey(coverSets, World_Country)
    forPlot <- coverSets[.(unique(getPopulations()$Country))]
    base <- ggplot(data = forPlot, aes(x = order, y = perc, group = country)) +
      geom_line(color= "#D1495B", size = 1) + ##d1495b
      geom_point(size = 1, color="#2E4057", group = 1.5) +
      #ylim(0, 100) +
      theme_minimal() +
      labs(#title = c,
           #subtitle = "Set Cover Solutions",
           x = "Epitope Rank",
           y = "Cumulative Individuals Covered \n (Percent)")
      # theme(plot.title = element_text(face = "bold", hjust = 0.5),
            # plot.subtitle = element_text(hjust = 0.5))`
    ggplotly(base, source = 'coverSetPlot')
  })
  coverSetList <- reactive({
    setkey(coverSets, World_Country)
    filteredSets <- coverSets[.(unique(getPopulations()$Country))]
    filteredSets2 <- filteredSets[, .(`N Peptides` = length(unique(Epitopes_order)),
                    Peptides = list(unique(Epitopes_order))), by = c('World_Country')]
    colnames(filteredSets2) <- gsub('World_', '', colnames(filteredSets2))
    # final <- filteredSets2[-order(`N Peptides`)]
    return(filteredSets2)
  })
  output$coverSetData <- renderDataTable(coverSetList()[,c('Country','N Peptides','Peptides')],
                                         extensions = 'Buttons', 
                                         options = list(scrollX = TRUE, scrollY = TRUE,
                                                        dom = 'Blfrtip',
                                                        buttons = c('copy', 'csv')))
}

##### ~~~~~ Run App ~~~~~ #####
shinyApp(ui = ui, server = server)




