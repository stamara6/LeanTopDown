#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

source("helpers.R", local = TRUE)
source("spectrum annotator.R", local = TRUE)
source("freqflyers.R", local = TRUE)
source("frgmntstats.R", local = TRUE)


library(shiny)
library(ggplot2)
library(reshape2)
library(ggcorrplot)
library(shinyjs)
library(logging)
library(shinycssloaders)
library(shinydashboard)
library(LeanTopDown)
library(gridExtra)

ui <- dashboardPage(

          dashboardHeader(title = "LeanTopDown"),

          dashboardSidebar(

            sidebarMenu(

              conditionalPanel(condition = "input.conditionedPanels == 'mir'"
                                       #checkboxInput("masslists", "Process Masslists", FALSE),
                                       ),

              conditionalPanel(condition = "input.conditionedPanels == 'cor'"
                                       #checkboxInput("masslists2", "Process Masslists", FALSE),
                                       ),

              conditionalPanel(condition = "input.conditionedPanels == 'spa'"
                                       ),

              conditionalPanel(condition = "input.conditionedPanels == 'frm'"
                                       ),

              conditionalPanel(condition = "input.conditionedPanels == 'frf'",
                                       fileInput("freqfile", "Frequent Flyers File", FALSE)),

              conditionalPanel(condition = "input.conditionedPanels == 'frs'"
                                       ),

              conditionalPanel(condition = "input.conditionedPanels == 'erp'"
                                       ),
              actionButton("plot", "Plot"),
              downloadButton("saveplot", "Save Plot"),
              radioButtons("filetype", "Save as:",
                           choices = c("pdf", "png")),
              textInput("filename", "File name:", ""),
              uiOutput("mzRange"),
              checkboxInput("masslists", "Process Masslists", FALSE),
              fileInput("fileIn", "Upload Masslists", TRUE),
              actionLink("selectall","Select All"),
              checkboxGroupInput("selectFiles", "Select Spectra", c("1" = "1"))


              )),
              dashboardBody(
                     tabsetPanel(
                       tabPanel("mirRaw Plot", box(withSpinner(plotOutput("mirPlot", hover = "plotHover",
                                                          brush = brushOpts("plotBrush", resetOnNew = TRUE),
                                                          dblclick = "plotDC"))),
                                                value = "mir"),
                       tabPanel("Correlation Plot", box(withSpinner(plotOutput("corPlot"))),  value = "cor"),
                       tabPanel("Spectrum Annotation", box(withSpinner(plotOutput("spectrum", hover = "spechover",
                                                                  brush = brushOpts("specBrush", resetOnNew = TRUE),
                                                                  dblclick = "specDC"))),

                                box(checkboxInput("facet", "Facetting", FALSE)),
                                value = "spa"),
                       tabPanel("Fragment Maps", box(withSpinner(plotOutput("fragmap", hover = "fmhover",
                                                            brush = brushOpts("fmBrush", resetOnNew = TRUE),
                                                            dblclick = "fmDC"))),


                                box(checkboxInput("facet2", "Facetting", FALSE),
                                    textInput("seq", "Sequence", "MKSVITTVVSAADAAGRFPSNSDLESIQGNIQRSAARLEAAEKLAGNHEAVVKEAGDACFAKYAYLKNPGEAGENQEKINKCYRDVDHYMRLVNYCLVVGGTGPLDEWGIAGAREVYRTLNLPTSAYVASIAYTRDRLCVPRDMSAQAGVEFSAYLDYLINALS"),
                                    #checkboxInput("septerm", "Separate Termini Processing", TRUE),
                                    numericInput("modnum", "Number of Modifications", 0, 0, 5, 1),
                                    #checkboxInput("bound", "Binding of Masslist", TRUE),
                                    numericInput("ppm", "Accuracy", 3.5, 0.5, 20, 0.5)),value = "frm"),
                       tabPanel("Frequent Flyers", box(withSpinner(plotOutput("freqfls", hover = "ffhover",
                                                            brush = brushOpts("ffBrush", resetOnNew = TRUE),
                                                            dblclick = "ffDC"))),
                                                                                       value = "frf"),
                       tabPanel("Fragmentation Stats", box(withSpinner(plotOutput("fragtypes", hover = "fthover",
                                                            brush = brushOpts("ftBrush", resetOnNew = TRUE),
                                                            dblclick = "ftDC"))),
                                                       box(withSpinner(plotOutput("zplot", hover = "zphover",
                                                            brush = brushOpts("zpBrush", resetOnNew = TRUE),
                                                            dblclick = "zpDC"))),

                                box(checkboxInput("facet3", "Facetting", FALSE),
                                    textInput("seq3", "Sequence", "MKSVITTVVSAADAAGRFPSNSDLESIQGNIQRSAARLEAAEKLAGNHEAVVKEAGDACFAKYAYLKNPGEAGENQEKINKCYRDVDHYMRLVNYCLVVGGTGPLDEWGIAGAREVYRTLNLPTSAYVASIAYTRDRLCVPRDMSAQAGVEFSAYLDYLINALS"),
                                    numericInput("modnum2", "Number of Modifications", 0, 0, 5, 1),
                                    numericInput("ppm2", "Accuracy", 3.5, 0.5, 20, 0.5)),
                                value = "frs"),
                       tabPanel("Energy Resolved Tracing", box(withSpinner(plotOutput("erplot", hover = "erhover",
                                                            brush = brushOpts("erBrush", resetOnNew = TRUE),
                                                            dblclick = "erDC"))),

                                box(fileInput("spectra", "Upload Spectra", TRUE),
                                    checkboxInput("facet4", "Facetting", FALSE),
                                    textInput("seq4", "Sequence", ""),
                                    checkboxInput("septerm4", "Separate Termini Processing", TRUE),
                                    numericInput("modnum4", "Number of Modifications", 0, 0, 5, 1),
                                    checkboxInput("bound4", "Binding of Masslist", TRUE),
                                    numericInput("ppm4", "Accuracy", 3.5, 0.5, 20, 0.5)),
                                value = "erp"),
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

   # browser()

    spectra <- names(df())

    updateCheckboxGroupInput(session, "selectFiles", "Select Spectra", choices = spectra,
                       selected = "")

  })



  observe({

    if(input$selectall == 0) return(NULL)
    else if (input$selectall%%2 == 0)
    {
      updateCheckboxGroupInput(session, "selectFiles", "Select Spectra", choices = names(df()))
    } else
    {
      updateCheckboxGroupInput(session, "selectFiles", "Select Spectra", choices = names(df()), selected = names(df()))
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

  # output$mzRange2 <- renderUI({
  #   data <- df2()
  #   data <- do.call("rbind", data)
  #   #if(exists(data$Mz)) colnames(data)[colnames(data) == "Mz"] <- "m.z"
  #   sliderInput("range", "m/z Range (for correlation):",
  #               min = min(data$m.z), max = max(data$m.z), value = c(1000,5000))
  # })


  inputPlotMir <- function() {
    mirplot
  }
  inputPlotCor <- function() {
    corplot
  }
  inputPlotSpa <- function() {
    spaplot
  }
  inputPlotFrm <- function() {
    frmplot
  }
  inputPlotFrf <- function() {
    frfplot
  }
  inputPlotFrs <- function() {
    grid.arrange(zplot, frtplot)
  }
  inputPlotErp <- function() {
    erpplot
  }

#browser()
  output$saveplot <- downloadHandler(filename = function() {if(input$filename != "") {fn <- paste0("_", input$filename)} else {fn <- ""}
                                                            paste0(input$conditionedPanels,
                                                                   fn,".",
                                                                   input$filetype)},
                                     content = function(file) {
                                       if(input$conditionedPanels == "mir") inputPlot <<- inputPlotMir()
                                       else if(input$conditionedPanels == "cor") inputPlot <<- inputPlotCor()
                                       else if(input$conditionedPanels == "spa") inputPlot <<- inputPlotSpa()
                                       else if(input$conditionedPanels == "frm") inputPlot <<- inputPlotFrm()
                                       else if(input$conditionedPanels == "frf") inputPlot <<- inputPlotFrf()
                                       else if(input$conditionedPanels == "frs") inputPlot <<- inputPlotFrs()
                                       else if(input$conditionedPanels == "erp") inputPlot <<- inputPlotErp()
                                       ggsave(file, plot = inputPlot, device = input$filetype)})

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
    mirplot <<- ggplot(raw1, aes(x = m.z, y = Intensity/max(Intensity)*100)) +
      theme_bw() +
      geom_histogram(stat = "identity", color = "black") +
      geom_histogram(data = raw2, aes(x = m.z, y = -Intensity/max(Intensity)*100), stat = "identity", color = "red") +
      geom_hline(yintercept = 0) +

      scale_y_continuous(limits = ranges$y) +
      scale_x_continuous(limits = ranges$x) +
      xlab("m/z") +
      ylab("Normalized Intensities")

    print(mirplot)
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

    corplot <<- ggcorrplot(corMatrix, method = "square", hc.order = TRUE, lab = TRUE)
    print(corplot)
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
    spaplot <<- annotate_spectrum(data = gList, ranges = ranges)
    print(spaplot)
  })

  output$fragmap <- renderPlot({

    #data <- chosenspectra()
    inFiles <- input$fileIn[input$fileIn$name %in% chosenspectra(),]
    seq <- input$seq
    rm("data_merged", "Data_list", pos = ".GlobalEnv")
    #gList <- data
    #browser()
    Data_list <<- import_masslists(bound = TRUE, mod.number = input$modnum, shiny = TRUE, files = inFiles)
    data_merged <<- prepare_data_frame_term_sep_aa_all(sep_term = TRUE,
                                                                        mod.number = input$modnum,
                                                                        ppm = input$ppm, tint = FALSE, seq = seq)
    frmplot <<- fragment_map(seq = seq)
    print(frmplot)
  })

  output$freqfls <- renderPlot({

    file <- input$freqfile
    frfplot <<- frqfls(file = file$datapath)
    frfplot
  })



  output$zplot <- renderPlot(
    {
      inFiles <- input$fileIn[input$fileIn$name %in% chosenspectra(),]
      seq <- input$seq
      rm("data_merged", "Data_list", pos = ".GlobalEnv")

      #browser()
      Data_list <<- import_masslists(bound = TRUE, mod.number = input$modnum, shiny = TRUE, files = inFiles)
      data_merged <<- prepare_data_frame_term_sep_aa_all(sep_term = TRUE,
                                                         mod.number = input$modnum,
                                                         ppm = input$ppm, tint = FALSE, seq = seq)
      if(length(unique(data_merged$Charge)) == 1) return()
      zplot <<- fragment_charges(data = data_merged)

      zplot

    }
  )


  output$fragtypes <- renderPlot(
    {
      inFiles <- input$fileIn[input$fileIn$name %in% chosenspectra(),]
      seq <- input$seq
      rm("data_merged", "Data_list", pos = ".GlobalEnv")

      #browser()
      Data_list <<- import_masslists(bound = TRUE, mod.number = input$modnum, shiny = TRUE, files = inFiles)
      data_merged <<- prepare_data_frame_term_sep_aa_all(sep_term = TRUE,
                                                         mod.number = input$modnum,
                                                         ppm = input$ppm, tint = FALSE, seq = seq)

      frtplot <<- fragment_types(data = data_merged, ppm = input$ppm)

      frtplot

    }
  )

  output$tscplot <- renderPlot(
    {
      inFiles <- input$fileIn[input$fileIn$name %in% chosenspectra(),]
      seq <- input$seq
      rm("data_merged", "Data_list", pos = ".GlobalEnv")

      #browser()
      Data_list <<- import_masslists(bound = TRUE, mod.number = input$modnum, shiny = TRUE, files = inFiles)
      data_merged <<- prepare_data_frame_term_sep_aa_all(sep_term = TRUE,
                                                         mod.number = input$modnum,
                                                         ppm = input$ppm, tint = FALSE, seq = seq)

      tscplot <<- termini_sequence_coverage(data = data_merged, ppm = input$ppm)

      tscplot

    }
  )

}



shinyApp(ui = ui, server = server)

