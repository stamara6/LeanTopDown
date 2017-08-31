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

ui <- tagList(fluidPage(useShinyjs(),
                        theme = "bootstrap.css",
            includeScript("./text.js"),
            titlePanel("LeanTopDown"),
            fluidRow(column(4, checkboxInput("facet", "Facetting", FALSE),
                               textInput("seq", "Sequence", ""),
                               checkboxInput("septerm", "Separate Termini Processing", TRUE),
                               numericInput("modnum", "Number of Modifications", 0, 0, 5, 1),
                               checkboxInput("bound", "Binding of Masslist", TRUE)),
                     column(4, numericInput("ppm", "Accuracy", 3.5, 0.5, 20, 0.5))),
            fluidRow(
              column(4,
                     wellPanel(
                       tags$div(class="form-group shiny-input-container",
                                tags$div(tags$label("File input")),
                                tags$div(tags$label("Choose folder", class="btn btn-primary",
                                                    tags$input(id = "fileIn", webkitdirectory = TRUE, type = "file", style="display: none;", onchange="pressed()"))),
                                tags$label("No folder chosen", id = "noFile"),
                                tags$div(id="fileIn_progress", class="progress progress-striped active shiny-file-input-progress",
                                         tags$div(class="progress-bar")
                                )
                       ),
                       uiOutput("mzRange"),

                       checkboxInput("masslists", "Process Masslists", FALSE),

                       actionLink("selectall","Select All"),

                       uiOutput("selectFiles")
                     )
              ),
              column(8,
                     tabsetPanel(
                       tabPanel("mirRaw Plot", plotOutput("mirPlot", hover = "plotHover",
                                                          brush = brushOpts("plotBrush", resetOnNew = TRUE), dblclick = "plotDC")),
                       tabPanel("Correlation Plot", plotOutput("corPlot")),
                       tabPanel("Spectrum Annotation", plotOutput("spectrum", hover = "spechover",
                                                                  brush = brushOpts("specBrush", resetOnNew = TRUE),
                                                                  dblclick = "specDC")),
                       tabPanel("Fragment Maps", plotOutput("fragmap", hover = "fmhover",
                                                            brush = brushOpts("fmBrush", resetOnNew = TRUE),
                                                            dblclick = "fmDC")),
                       tabPanel("Frequent Flyers", plotOutput("freqfls", hover = "ffhover",
                                                            brush = brushOpts("ffBrush", resetOnNew = TRUE),
                                                            dblclick = "ffDC")),
                       tabPanel("Z Plot", plotOutput("zplot", hover = "zphover",
                                                            brush = brushOpts("zpBrush", resetOnNew = TRUE),
                                                            dblclick = "zpDC")),
                       tabPanel("Fragmentation Stats", plotOutput("fragstat", hover = "fshover",
                                                            brush = brushOpts("fsBrush", resetOnNew = TRUE),
                                                            dblclick = "fsDC")),
                       tabPanel("Energy Resolved Tracing", plotOutput("erplot", hover = "erhover",
                                                            brush = brushOpts("erBrush", resetOnNew = TRUE),
                                                            dblclick = "erDC"))

                     ),
                     verbatimTextOutput("info")
              )
            )
),
HTML("<script type='text/javascript' src='getFolders.js'></script>")
)



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


  output$selectFiles <- renderUI({
    spectra <- names(df())
    checkboxGroupInput("spectra", "Select Spectra", spectra)
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
    sliderInput("range", "m/z Range (for correlation):",
                min = min(data$`m.z|Mz`), max = max(data$`m.z|Mz`), value = c(1000,5000))
  })

  output$mirPlot <- renderPlot({
    data <- df()
    subspectra <- input$spectra
    raw1 <- data[[subspectra[1]]]
    raw2 <- data[[subspectra[2]]]


    if(is.null(ranges$x)) {ranges$x <- c(input$range[1], input$range[2])
                           ranges$y <- c(-100, 100)}
    colnames(raw1)[grepl("m.z|Mz|M.z|M.Z")] <- "Mz"
    ggplot(raw1, aes(x = Mz, y = modInt)) +
      theme_bw() +
      geom_histogram(stat = "identity", color = "black") +
      geom_histogram(data = raw2, aes(x = Mz, y = -modInt), stat = "identity", color = "red") +
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
    subspectra <- input$spectra

    ranges <- seq(input$range[1], input$range[2], by = 0.1)

    subdataList <- lapply(subspectra, function(x)
    {
      df <- data[[x]]
      #df$exp <- names(data[x])
      df <- tapply(df$modInt, cut(df$Mz, ranges), sum, na.rm = TRUE)
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
    inFiles <- input$fileIn
    gList <<- data
    LeanTopDown::annotate_spectrum(data = gList)

  })

  output$fragmap <- renderPlot({

    data <- df()
    inFiles <- input$fileIn
    seq <- input$seq
    gList <- data
    Data_list <<- LeanTopDown::import_masslists(bound = input$bound, mod.number = input$modnum, shiny = TRUE)
    data_merged <<- LeanTopDown::prepare_data_frame_term_sep_aa_all(sep_term = input$septerm,
                                                                        mod.number = input$modnum,
                                                                        ppm = input$ppm)
    LeanTopDown::fragment_map()

  })

}

shinyApp(ui = ui, server = server)

