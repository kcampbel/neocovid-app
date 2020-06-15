library(shiny)
library(shinydashboard)
library(DT)
library(data.table)
library(plotly)
library(tidyverse)
library(leaflet)
library(wesanderson)

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
      #
      menuItem("Putative Epitopes", tabName = "epitope", icon = icon("flag")),
      div( id = 'sidebar_epitope',
           conditionalPanel("input.sidebar == 'epitope'"
                            # selectizeInput("epitopeInput", "Epitope of Interest", choices = NULL, selected = NULL, multiple = TRUE),

                            )
           ),
      # 
      menuItem("HLA", tabName = "hla",icon = icon("fingerprint")),
      div( id = 'sidebar_hla',
           conditionalPanel("input.sidebar == 'hla'",
                            # textInput("hlaInput", "Search HLA", value = "HLA-A*01:01"),
                            selectizeInput("hlaInput", "Search HLA", choices = c("A*68:44"), selected = "A*68:44"),
                            selectInput("hlaProteinInput", "Select a Protein", choices = c("Show All"), selected = "Show All")
                            )
           ),
      #
      menuItem("Population", tabName = "populations", icon = icon("globe")),
      div( id = 'sidebar_populations',
           conditionalPanel("input.sidebar == 'populations'",
                            selectizeInput("countryInput", "Country(s)", choices = c("United States"), selected = "United States", multiple = TRUE),
                            numericInput("allelicFreq", "Minimum Population Allelic Frequency", 0.05, min = 0, max = 1)
                            )
           )
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "virus",
              box(title = "Viral View", width = 12,
                  plotlyOutput("mainPlotly")
                  ),
              box(title = "Predicted Epitopes", width = 12,
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
                box(width = 6, height = 300
                    # textOutput("iedbTextOutput")
                    )
              ),
              box(title = "Predicted Epitopes", width = 12,
                  dataTableOutput("epitopeDataTable")
              )
      ),
      tabItem(tabName = "hla",
              box(title = "HLA View", width = 12, height = 600,
                  plotlyOutput("hlaPlotly")
                  ),
              box(title = "Predicted Epitopes", width = 12,
                  checkboxInput("otherHlaFamilyMembers", "Show peptides from other family members"),
                  dataTableOutput("hlaPeptideTable")
                  # dataTableOutput("test")
                  )
              ),
      tabItem(tabName = "populations",
              box(leafletOutput("worldmap"), width = 12),
              # box(dataTableOutput("alleleN"), width = 6),
              box(dataTableOutput("alleleBreakdown"), width = 12)
      )
    )
  )
)

##### ~~~~~ SERVER ~~~~~ #####
# Define server logic required to draw a histogram ----
server <- function(input, output, session) {
  # load('~/neocovid-app/app/data/forApp.Rda')
  load(file = url("https://github.com/kcampbel/neocovid-app/raw/master/app/data/forApp.Rda"))
  # download.file(url = "https://github.com/kcampbel/neocovid-app/blob/master/app/data/forApp.Rda", "myfile")
  # load("myfile")
  
  
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
      scale_y_continuous(limits = c(-8,max(plotEpis$y_new)+1)) + theme_bw() +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
      facet_grid(Class ~ ., scales = 'free_y')
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
  
  # subsetIedb <- reactive({
  #   if(input$epitopeInput == "") {
  #     return("No epitopes found in IEDB.")
  #   } else {
  #     sub_string <- grep(input$epitopeInput, uniqueIedbAnt, value = TRUE)
  #     get_dist <- lapply(uniqueIedbAnt, function(y){ adist(input$epitopeInput, y) }) %>% unlist
  #     closest_peptides <- uniqueIedbAnt[which(get_dist<3)]
  #     message <- paste0(
  #       ifelse(length(sub_string)>0, paste0("Epitope ", input$epitopeInput, " found in IEDB.\n"), ""),
  #       "There were ", length(union(sub_string, get_dist), " epitopes that contain \'", input$epitopeInput,"\'.\n"),
  #       stringr::str_wrap(paste0(union(sub_string, get_dist), collapse = ", "), 50)
  #     )
  #     return(message)
  #   }
  # })
  
  # output$iedbTextOutput <- renderText(subsetIedb())
  
  output$epitopeDataTable <- renderDataTable(subsetEpitopes(), options = list(scrollX = TRUE, scrollY = TRUE))
  
  ##### HLA View #####
  # Get HLA options
  observe({
    updateSelectizeInput(
      session,
      "hlaInput",
      choices = unique(summPredictions$`HLA Allele`),
      selected = "A*68:44"
    )
  })
  
  # Get protein name options
  observe({
    updateSelectInput(
      session,
      "hlaProteinInput",
      choices = c("Show All", unique(as.character(proteinSequences3$Protein))),
      selected = "Show All"
    )
  })

  # Get HLA protein family
  hlaFamily <- reactive({
    getHlaFamily <- gsub("([A-C]\\*\\d+):\\d+", "\\1", input$hlaInput)
    return(getHlaFamily)
  })

  # Filter to HLA family members
  filtHlaFamily <- reactive({
    filtHlaFamily <- summPredictions2 %>% ungroup %>%
      mutate(HlaFamily = gsub("([A-C]\\*\\d+):\\d+", "\\1", `HLA Allele`)) %>%
      filter(HlaFamily == hlaFamily() & `HLA Allele` != input$hlaInput)
    return(filtHlaFamily)
  })

  # Subset to selected protein
  subsetToProtein <- reactive({
    if(input$hlaProteinInput == "Show All"){
      subsetToProtein <- proteinSequences3
    } else {
      subsetToProtein <- filter(proteinSequences3, Protein == input$hlaProteinInput)
    }
    return(subsetToProtein)
  })

  output$hlaPlotly <- renderPlotly({
    plotData <- subsetToProtein()

    summHlaFamily <- filtHlaFamily() %>% group_by(Peptide, peptide.polar_start) %>% summarise(n = n())
    filterPreds <- filter(summPredictions2, `HLA Allele` == input$hlaInput)

    basePlot <- ggplot(data = NULL) +
      geom_hline(data = NULL, yintercept = -25, colour = 'grey80', size = 2) +
      # Peptides from the same HLA family members
      geom_segment(data = summHlaFamily, aes(x = peptide.polar_start, xend = peptide.polar_start,
                                             y = 510, yend = 590, colour = n)) +
      scale_colour_gradientn(colours = wes_palette('Zissou1', 50, type = 'continuous')) +
      # Chains
      geom_rect(data = plotData, aes(fill = Gene, xmin = chain_polar_start, xmax = chain_polar_end,
                                              ymin = -50, ymax = 0),
                fill = 'grey80', alpha = 0.3) +
      # Proteins
      geom_rect(data = plotData, aes(fill = Gene, xmin = polar_start, xmax = polar_end,
                                              ymin = -40, ymax = -10)) +
      scale_fill_manual(values = sample(colors, 12)) +
      # Peptides from specific HLA type
      geom_col(data = filterPreds, aes(fill = Gene, x = peptide.polar_start,
                                       y = 500-`Median Affinity (nM)`), width = 5) +
      #
      scale_x_continuous(limits = c(min(plotData$polar_start), max(plotData$polar_end))) +
      scale_y_continuous(limits = c(-50,600), breaks = c(0, 100, 200, 300, 400, 500, 550),
                         labels = c(500,400,300,200,100,0, paste0('Other ',hlaFamily()," Alleles"))) +
      labs(y = "Predicted Binding Affinity (nM)", x = "Position in SARS-CoV-2",
           colour = "N Other Alleles") +
      theme_bw() + 
      theme()

    ggplotly(basePlot, source = "hlaPeptidePlot", height = 550) %>% 
      layout(xaxis = list(rangeslider = list(type = "date")))
  })
  
  # Subset predictions to those in view (whether it's full, a single protein, or zoomed in)
  filterHlaPredictions <- reactive({
    zoom <- event_data("plotly_relayout", "hlaPeptidePlot")
    # if plot just rendered, event_data is NULL
    # if user double clicks for autozoom, then zoom$xaxis.autorange is TRUE
    # if user resizes page, then zoom$width is pixels of plot width
    # filterPreds <- filter(summPredictions2, `HLA Allele` == pickHLA)
    # filterPredsPlusHlaFamily <- filter(summPredictions2, `HLA Allele` == pickHLA | `HLA Allele` %in% filtHlaFamily$`HLA Allele`)
    # 
    if(input$hlaProteinInput == "Show All" & (is.null(zoom) || names(zoom[1]) %in% c("xaxis.autorange", "width"))) {
      if(input$otherHlaFamilyMembers == TRUE) {
        return(filter(summPredictions, `HLA Allele` == input$hlaInput | `HLA Allele` %in% filtHlaFamily()$`HLA Allele`))
      } else {
        return(filter(summPredictions, `HLA Allele` == input$hlaInput))
      }
    } else if(input$hlaProteinInput != "Show All") {
      if(input$otherHlaFamilyMembers == TRUE) {
        return(filter(summPredictions, Protein == input$hlaProteinInput & (`HLA Allele` == input$hlaInput | `HLA Allele` %in% filtHlaFamily()$`HLA Allele`)))
      } else {
        return(filter(summPredictions, Protein == input$hlaProteinInput & `HLA Allele` == input$hlaInput))
      }
    } else {
      xmin <- zoom$`xaxis.range[0]`
      xmax <- zoom$`xaxis.range[1]`
      get_peptides <- filter(plotEpis, peptide.polar_start > xmin & peptide.polar_end < xmax)$Peptide
      if(input$otherHlaFamilyMembers == TRUE) {
        return(filter(summPredictions, Peptide %in% get_peptides & (`HLA Allele` == input$hlaInput | `HLA Allele` %in% filtHlaFamily()$`HLA Allele`)))
      } else {
        return(filter(summPredictions, Peptide %in% get_peptides & `HLA Allele` == input$hlaInput))
      }
    }
  })
  # Output predictions table
  output$hlaPeptideTable <- renderDataTable(filterHlaPredictions(), options = list(scrollX = TRUE, scrollY = TRUE))

  ##### Population View #####
  mypalette <- colorFactor( palette=colors, domain=world@data$PICK, na.color = "transparent" )
  mytext <- world@data$NAME %>% lapply(htmltools::HTML)
  #
  observe({
    updateSelectizeInput(
      session,
      "countryInput",
      choices = unique(filter(world@data, !is.na(PICK))$NAME),
      selected = "United States"
    )
  })
  #
  getPopulationBreakdown <- reactive({
    countryPopFreq %>%
      mutate(Gene = paste0("HLA-", Locus)) %>%
      filter(Country %in% input$countryInput & `Allele Frequency`>=input$allelicFreq)
    # fulldf %>% 
    #   mutate(Country = World_Country, Gene = paste0("HLA-", Locus)) %>%
    #   dplyr::select(Region, Country, Population, `Sample Size`,
    #                 Gene, Allele, Superfamily, 
    #                 `Allele Frequency`,
    #                 `% of individuals that have the allele`) %>%
    #   filter(Country %in% input$countryInput & `Allele Frequency`>=input$allelicFreq)
  })
  #
  getAlleleBreakdown <- reactive({
    uniqueAlleles <- unique(getPopulationBreakdown()$`HLA Allele`)
    
  })
  
  #
  autoZoomMap <- reactive({
    world@data %>% filter(NAME %in% input$countryInput)
  })
  #
  output$worldmap <- renderLeaflet({
    leaflet(world) %>%
      addProviderTiles(providers$CartoDB.Positron,
      options = providerTileOptions(noWrap = TRUE)) %>%
      addTiles() %>% 
      setView(lng = autoZoomMap()$LON[1], lat = autoZoomMap()$LAT[1], zoom = 3) %>%
      addCircleMarkers(~LON, ~LAT,
                       fillColor = ~mypalette(PICK), 
                       fillOpacity = 0.7, 
                       color = "white", 
                       radius = 8, 
                       stroke = FALSE,
                       label = mytext,
                       labelOptions = labelOptions( style = list("font-weight" = "normal",
                                                                 padding = "3px 8px"),
                                                    textsize = "13px",
                                                    direction = "auto")
                       )
  })
  
  # output$alleleN <- renderDataTable( getAlleleBreakdown() )
  output$alleleBreakdown <- renderDataTable( getPopulationBreakdown() )
  
}

##### ~~~~~ Run App ~~~~~ #####
shinyApp(ui = ui, server = server)




