shinyServer(
  function(input, output, session) {
    
    observeEvent(input$run, {
      updateNavbarPage(session, "master",
                       selected = "QC Results"
      )
    })
    
    ns <- eventReactive(input$run, {
      
      if (input$hk == "") {
        hk.genes <- NULL
      } else {
        hk.genes <- gsub(" ", "", strsplit(input$hk, split=",")[[1]])
      }
      
      nanostringData <- processNanostringData(input$expr$datapath,
                                              bgType = "t.test", bgPVal = input$bgP,
                                              housekeeping = hk.genes,
                                              includeQC = TRUE)
      
      nanostringDataBG <- processNanostringData(input$expr$datapath,
                                                bgType = "threshold", bgThreshold = 2, bgProportion = 0.5,
                                                housekeeping = hk.genes,
                                                includeQC = TRUE)
      
      colnames(nanostringData$exprs) <- 
        colnames(nanostringData$exprs.raw) <-
        names(nanostringData$pc.scalefactors) <-
        names(nanostringData$hk.scalefactors) <-
        colnames(nanostringData$qc) <- 
        rownames(nanostringDataBG$bg.stats) <-
        gsub(".*\\/|\\.RCC", "", colnames(nanostringData$exprs))
      
      if (is.null(input$phen$datapath)) {
        groups <- gsub("_.*", "", colnames(nanostringData$exprs))
      } else {
        groups <- as.character(read.table(file = input$phen$datapath, sep = "", as.is = TRUE))
      }
      
      if (input$basePhen == "") {
        base.group <- groups[1]
      } else {
        base.group <- input$basePhen
      }
      
      limmaResults <- runLimmaAnalysis(nanostringData, groups, base.group)
      
      ns <- list(dat = nanostringData,
                 bg.stats = nanostringDataBG$bg.stats,
                 deRes = limmaResults,
                 base.group = base.group)
      
      if (!is.null(input$gsDb$datapath)) {
        ns$gsRes <- mHG.genesets.multi(limmaResults, input$gsDb$datapath,
                                       t_thresh = input$setThresh, numReq_inSet = input$minSize, numReq_overall = 10)
      }
      
      return(ns)
    })
    
    posQC <- reactive({ positiveQC(ns()) })
    negQC <- reactive({ negativeQC(ns()) })
    hkQC <- reactive({ housekeepingQC(ns()) })
    pcaPlot <- reactive({ plotPCA(ns()) })
    deResults <- reactive({ deRes(ns()) })
    
    output$posTab <- renderTable({posQC()$tab})
    output$posPlot <- renderPlot({posQC()$plt})
    output$negTab <- renderTable({negQC()$tab})
    output$negPlot <- renderPlotly({negQC()$plt})
    output$hkTab <- renderTable({hkQC()$tab})
    output$hkPlot1 <- renderPlotly({hkQC()$plt1})
    output$hkPlot2 <- renderPlotly({hkQC()$plt2})
    output$pcaPlot <- renderPlotly({pcaPlot()})
    output$deCounts <- renderTable({deResults()$summary})
    output$deTab <- renderDT({deResults()$de}, rownames = TRUE)
    
    
    output$gsUI <- renderUI({
      #inputPanel(
          fluidRow(
      #      column(12,
                 selectInput("gsComp", label = "Comparison:",
                             choices = names(ns()$gsRes)),
      #    ),
      #    column(12,
                 radioButtons("gsDir", label = "Direction:",
                              choices = c("Higher in Group" = 1,
                                          "Higher in Base" = 2), selected = 1),
                 sliderInput("gsJac", label = "Clustering Threshold (Jaccard Index):",
                             min = 0, max = 1, value = 0.5, step = 0.05),
      #    ),
      #    column(12,
                 numericInput("gsQthresh", label = "q-value threshold:", value = 0.05, step = 0.05),
                 numericInput("gsClust", label = "Cluster to plot:", value = 1, step = 1)
      #    )
        )
      #)
    })

      leadingEdge <- reactive({
      mHG.to.LEdge.multi(ns()$gsRes, ns()$deRes,
                         t.thresh = 2, q.thresh = input$gsQthresh, input$gsDb$datapath)
    })

      groupedGenesets <- reactive({    
        groupByDistance(ns()$gsRes[[input$gsComp]][[as.numeric(input$gsDir)]],
                        leadingEdge()[[input$gsComp]][[as.numeric(input$gsDir)]],
                        join.threshold = (1-input$gsJac))
    })

    output$gsTab <- renderDT({
      datatable(groupedGenesets(),
                rownames = FALSE,
                options = list(autoWidth = TRUE,
                               columnDefs = list(list(width = '200px', targets = "_all")))) %>%
        formatStyle('best', target = 'row',
                    color = styleEqual(c("", "x"), c('grey', 'black')),
                    fontSize = '12px')
    }, rownames = FALSE
    )
    
    output$gsHM <- renderPlotly({
      plotlyHeatmap(ns(), groupedGenesets(), leadingEdge(), input$gsClust, input$gsComp, input$gsDir)
    })

    
  }
)