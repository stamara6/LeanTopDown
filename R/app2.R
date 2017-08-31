#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(reshape2)
library(ggcorrplot)
library(shinyjs)
library(logging)
library(shinycssloaders)
library(shinydashboard)

ui <- dashboardPage(

          dashboardHeader(title = "LeanTopDown"),

          dashboardSidebar(

            sidebarMenu(

              conditionalPanel(condition = "input.conditionedPanels == 'mir'"
                                       #checkboxInput("masslists", "Process Masslists", FALSE),
                                       ),

              conditionalPanel(condition = "input.conditionedPanels == 'cor'",
                                       #checkboxInput("masslists2", "Process Masslists", FALSE),
                                       actionLink("selectall2","Select All")),

              conditionalPanel(condition = "input.conditionedPanels == 'spa'",
                                       checkboxInput("masslists", "Process Masslists", FALSE),
                                       actionLink("selectall3","Select All")),

              conditionalPanel(condition = "input.conditionedPanels == 'frm'",
                                       checkboxInput("facet", "Facetting", FALSE),
                                       textInput("seq", "Sequence", ""),
                                       checkboxInput("septerm", "Separate Termini Processing", TRUE),
                                       numericInput("modnum", "Number of Modifications", 0, 0, 5, 1),
                                       checkboxInput("bound", "Binding of Masslist", TRUE),
                                       numericInput("ppm", "Accuracy", 3.5, 0.5, 20, 0.5)),

              conditionalPanel(condition = "input.conditionedPanels == 'frf'",
                                       fileInput("freqfile", "Frequent Flyers File", FALSE)),

              conditionalPanel(condition = "input.conditionedPanels == 'zpl'",
                                       fileInput("masslists2", "Upload Masslists", TRUE),
                                       checkboxInput("facet2", "Facetting", FALSE),
                                       textInput("seq2", "Sequence", ""),
                                       checkboxInput("septerm2", "Separate Termini Processing", TRUE),
                                       numericInput("modnum2", "Number of Modifications", 0, 0, 5, 1),
                                       checkboxInput("bound2", "Binding of Masslist", TRUE),
                                       numericInput("ppm2", "Accuracy", 3.5, 0.5, 20, 0.5)),

              conditionalPanel(condition = "input.conditionedPanels == 'frs'",
                                       fileInput("masslists3", "Upload Masslists", TRUE),
                                       checkboxInput("facet3", "Facetting", FALSE),
                                       textInput("seq3", "Sequence", ""),
                                       checkboxInput("septerm3", "Separate Termini Processing", TRUE),
                                       numericInput("modnum3", "Number of Modifications", 0, 0, 5, 1),
                                       checkboxInput("bound3", "Binding of Masslist", TRUE),
                                       numericInput("ppm3", "Accuracy", 3.5, 0.5, 20, 0.5)),

              conditionalPanel(condition = "input.conditionedPanels == 'erp'",
                                       fileInput("spectra", "Upload Spectra", TRUE),
                                       checkboxInput("facet4", "Facetting", FALSE),
                                       textInput("seq4", "Sequence", ""),
                                       checkboxInput("septerm4", "Separate Termini Processing", TRUE),
                                       numericInput("modnum4", "Number of Modifications", 0, 0, 5, 1),
                                       checkboxInput("bound4", "Binding of Masslist", TRUE),
                                       numericInput("ppm4", "Accuracy", 3.5, 0.5, 20, 0.5)),
              actionButton("plot", "Plot"),
              uiOutput("mzRange"),
              fileInput("fileIn", "Upload Masslists", TRUE),
              actionLink("selectall","Select All"),
              checkboxGroupInput("selectFiles", "Select Spectra", c("1" = "1"))


              )),
              dashboardBody(
                     tabsetPanel(
                       tabPanel("mirRaw Plot", box(withSpinner(plotOutput("mirPlot", hover = "plotHover",
                                                          brush = brushOpts("plotBrush", resetOnNew = TRUE),
                                                          dblclick = "plotDC"))), value = "mir"),
                       tabPanel("Correlation Plot", box(withSpinner(plotOutput("corPlot"))), value = "cor"),
                       tabPanel("Spectrum Annotation", box(withSpinner(plotOutput("spectrum", hover = "spechover",
                                                                  brush = brushOpts("specBrush", resetOnNew = TRUE),
                                                                  dblclick = "specDC"))), value = "spa"),
                       tabPanel("Fragment Maps", box(withSpinner(plotOutput("fragmap", hover = "fmhover",
                                                            brush = brushOpts("fmBrush", resetOnNew = TRUE),
                                                            dblclick = "fmDC"))), value = "frm"),
                       tabPanel("Frequent Flyers", box(withSpinner(plotOutput("freqfls", hover = "ffhover",
                                                            brush = brushOpts("ffBrush", resetOnNew = TRUE),
                                                            dblclick = "ffDC"))), value = "frf"),
                       tabPanel("Z Plot", box(withSpinner(plotOutput("zplot", hover = "zphover",
                                                            brush = brushOpts("zpBrush", resetOnNew = TRUE),
                                                            dblclick = "zpDC"))), value = "zpl"),
                       tabPanel("Fragmentation Stats", box(withSpinner(plotOutput("fragstat", hover = "fshover",
                                                            brush = brushOpts("fsBrush", resetOnNew = TRUE),
                                                            dblclick = "fsDC"))), value = "frs"),
                       tabPanel("Energy Resolved Tracing", box(withSpinner(plotOutput("erplot", hover = "erhover",
                                                            brush = brushOpts("erBrush", resetOnNew = TRUE),
                                                            dblclick = "erDC"))), value = "erp"),
                       id = "conditionedPanels"

                     ),
                     verbatimTextOutput("info")
              )


)

#basicConfig()

#options(shiny.error = browser)


server <- function(input, output, session) {

  df <- reactive({

    inFiles <- input$fileIn

    df <- list()
    if (is.null(inFiles))
      return(NULL)
    for (i in seq_along(inFiles$datapath)) {
      if(!grepl(".txt|.masslist|.csv",inFiles$datapath[i])) next

      df[[i]] <- read.table(inFiles$datapath[i], header = TRUE, sep = "\t")
      if("Intensity" %in% colnames(df[[i]]))df[[i]]$modInt <- 100/max(df[[i]]$Intensity)*df[[i]]$Intensity

    }
    df <- df[!sapply(df, is.null)]
    names(df) <- inFiles$name[grepl(".txt|.masslist|.csv", inFiles$name)]
    if(input$masslists == TRUE) {
      inFiles <<- input$fileIn
      df <- prepare_masslists(shiny = TRUE)
    }

    df

  })

  df2 <- reactive({

    inFiles <- input$fileIn2

    df <- list()
    if (is.null(inFiles))
      return(NULL)
    for (i in seq_along(inFiles$datapath)) {
      if(!grepl(".txt|.masslist|.csv",inFiles$datapath[i])) next

      df[[i]] <- read.table(inFiles$datapath[i], header = TRUE, sep = "\t")
      if("Intensity" %in% colnames(df[[i]]))df[[i]]$modInt <- 100/max(df[[i]]$Intensity)*df[[i]]$Intensity

    }
    df <- df[!sapply(df, is.null)]
    names(df) <- inFiles$name[grepl(".txt|.masslist|.csv", inFiles$name)]

    df

  })


  observe({

    #browser()

    spectra <- names(df())
    #spectra2 <- names(df2())

    updateCheckboxGroupInput(session, "selectFiles", "Select Spectra", choices = spectra,
                       selected = "")

    #updateCheckboxGroupInput(session, "selectFiles", "Select Spectra", choices = spectra2,
     #                        selected = "")

  })



  observe({

    if(input$selectall == 0) return(NULL)
    else if (input$selectall%%2 == 0)
    {
      updateCheckboxGroupInput(session, "spectra", "Select Spectra", choices = names(df()))
    } else
    {
      updateCheckboxGroupInput(session, "spectra", "Select Spectra", choices = names(df()), selected = names(df()))
    }
  })

  ranges <- reactiveValues(x = NULL, y = NULL)



  observeEvent(input$plotDC, {
    brush <- input$plotBrush
    if(!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)
    } else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
  })

  output$mzRange <- renderUI({
    data <- df()
    data <- do.call("rbind", data)
    #if(exists(data$Mz)) colnames(data)[colnames(data) == "Mz"] <- "m.z"
    sliderInput("range", "m/z Range (for correlation):",
                min = min(data$m.z), max = max(data$m.z), value = c(min(data$m.z), max(data$m.z)))
  })

  output$mzRange2 <- renderUI({
    data <- df2()
    data <- do.call("rbind", data)
    #if(exists(data$Mz)) colnames(data)[colnames(data) == "Mz"] <- "m.z"
    sliderInput("range", "m/z Range (for correlation):",
                min = min(data$m.z), max = max(data$m.z), value = c(1000,5000))
  })

  chosenspectra <- eventReactive(input$plot, {
    subspectra <- input$selectFiles

  })

  output$mirPlot <- renderPlot({
    #browser()
    data <- df()
    if(is.null(ranges$x)) {ranges$x <- c(input$range[1], input$range[2])
    ranges$y <- c(-100, 100)}
    subspectra <- chosenspectra()
    raw1 <- data[[subspectra[1]]]
    raw2 <- data[[subspectra[2]]]




    #colnames(raw1)[grepl("m.z|Mz|M.z|M.Z")] <- "Mz"
    ggplot(raw1, aes(x = m.z, y = Intensity/max(Intensity)*100)) +
      theme_bw() +
      geom_histogram(stat = "identity", color = "black") +
      geom_histogram(data = raw2, aes(x = m.z, y = -Intensity/max(Intensity)*100), stat = "identity", color = "red") +
      geom_hline(yintercept = 0) +

      scale_y_continuous(limits = ranges$y) +
      scale_x_continuous(limits = ranges$x) +
      xlab("m/z") +
      ylab("Normalized Intensities")
    }
  )

  output$info <- renderText({

    mzInt <- function(e) {
    if(is.null(e)) return("NULL\n")
      paste0("m/z = ", e$x, "\n")
    }

    mzInt(input$plotHover)
  })

  output$corPlot <- renderPlot({

    data <- df()
    subspectra <- chosenspectra()

    ranges <- seq(input$range[1], input$range[2], by = 0.1)

    subdataList <- lapply(subspectra, function(x)
    {
      df <- data[[x]]
      #df$exp <- names(data[x])
      df <- tapply(df$Intensity, cut(df$m.z, ranges), sum, na.rm = TRUE)
    })

    subdataDf <- do.call("cbind", subdataList)
    colnames(subdataDf) <- subspectra

    corMatrix <- cor(subdataDf, use = "pairwise.complete.obs", method = "pearson")

    ggcorrplot(corMatrix, method = "square", hc.order = TRUE, lab = TRUE)

  },
  width = 1000,
  height = 1000)

  output$spectrum <- renderPlot({

    data <- df()
    subspectra <- chosenspectra()
#browser()
    inFiles <- input$fileIn
    gList <<- data[subspectra]
    ranges$y <<- c(0, 100)
    LeanTopDown::annotate_spectrum(data = gList, ranges = ranges)
  })

  output$fragmap <- renderPlot({

    #data <- chosenspectra()
    inFiles <- input$fileIn[input$fileIn$name %in% chosenspectra(),]
    seq <- input$seq
    rm("data_merged", "Data_list", pos = ".GlobalEnv")
    #gList <- data
    #browser()
    Data_list <<- LeanTopDown::import_masslists(bound = input$bound, mod.number = input$modnum, shiny = TRUE, files = inFiles)
    data_merged <<- LeanTopDown::prepare_data_frame_term_sep_aa_all(sep_term = input$septerm,
                                                                        mod.number = input$modnum,
                                                                        ppm = input$ppm, tint = FALSE, seq = seq)
    LeanTopDown::fragment_map(seq = seq)

  })

  output$freqfls <- renderPlot({

    file <- input$freqfile
    LeanTopDown::frqfls(file = file$datapath)

  })


  output$zplot <- renderPlot(
    {
      data <- chosenspectra()

    }
  )

}



shinyApp(ui = ui, server = server)

