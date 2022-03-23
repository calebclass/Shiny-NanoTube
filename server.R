shinyServer(
  function(input, output, session) {
    options(shiny.maxRequestSize=100*1024^2)
    
    hide(id = "normTxt")
    hide(id = "hk")
    hide(id = "bgP")
    hide(id = "gseaTxt")
    hide(id = "minSize")
    
    observeEvent(input$adv, {
      toggle(id = "normTxt")
      toggle(id = "hk")
      toggle(id = "bgP")
      toggle(id = "gseaTxt")
      toggle(id = "minSize")
    })
    
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
                                              includeQC = FALSE)
      
      nanostringDataBG <- processNanostringData(input$expr$datapath,
                                                bgType = "threshold", bgThreshold = 2, bgProportion = 0.5,
                                                housekeeping = hk.genes,
                                                includeQC = TRUE,
                                                output.format = "list")
      
      
      file_input <- input$expr$datapath
      file.extension <- substr(file_input[1], 
                               (nchar(file_input[1])-3), nchar(file_input[1]))
      if(file.extension %in% c(".zip",".ZIP")){
        file_input <- unzip(file_input)
        file_input <-read_merge_rcc(file_input)
        nanoTableData <- file_input
        output$nanoTable <- renderTable(file_input)
      }
      
      colnames(nanostringData) <- 
        colnames(nanostringDataBG$exprs.raw) <-
        names(nanostringDataBG$pc.scalefactors) <-
        names(nanostringDataBG$hk.scalefactors) <-
        rownames(nanostringDataBG$qc) <- 
        rownames(nanostringDataBG$bg.stats) <-
        gsub(".*\\/|\\.RCC", "", colnames(nanostringData))
      
      if (is.null(input$phen$datapath)) {
        groups <- gsub("_.*", "", colnames(nanostringData))
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
                 dat.list = nanostringDataBG,
                 deRes = limmaResults,
                 base.group = base.group)
      
      if (!is.null(input$gsDb$datapath)) {
        ns$gsRes <- limmaToFGSEA(limmaResults, input$gsDb$datapath,
                                 min.set = input$minSize)
      }
      
      return(ns)

    })
    
    posQC <- reactive({ positiveQC(ns()$dat.list) })
    negQC <- reactive({ negativeQC(ns()$dat.list) })
    hkQC <- reactive({ housekeepingQC(ns()$dat.list) })
    pcaPlot <- reactive({ plotPCA(ns()) })
    deResults <- reactive({ deRes(ns()) })
    ####
    canoPlot <- reactive({deVolcano(ns()$deRes) +geom_hline(yintercept = 2, linetype = "dashed", colour = "darkred") +
        geom_vline(xintercept = 0.5, linetype = "dashed", colour = "darkred") +
        geom_vline(xintercept = -0.5, linetype = "dashed", colour = "darkred")})
    ###
    
    #####

    output$posTab <- renderTable({posQC()$tab})
    output$posPlot <- renderPlot({posQC()$plt})
    output$negTab <- renderTable({negQC()$tab}, rownames = TRUE)
    output$negPlot <- renderPlotly({negQC()$plt})
    output$hkTab <- renderTable({hkQC()$tab})
    output$hkPlot1 <- renderPlotly({hkQC()$plt1})
    output$hkPlot2 <- renderPlotly({hkQC()$plt2})
    output$pcaPlot <- renderPlotly({pcaPlot()})
    output$deCounts <- renderTable({deResults()$summary},
                                   rownames = TRUE, round = 0)
    output$deTab <- renderDT({deResults()$de}, rownames = TRUE)
    
    ####
    output$DEdownload <- downloadHandler(
      filename = function() {"de.csv"},
      content = function(file) {
        write.csv(deResults(), file)
      }
    )
 
    ###
    output$canoPlot <- renderPlotly({canoPlot()})
    ###
    output$NANOdownload <- downloadHandler(
      filename = function() {"nanoTable.csv"},
      content = function(file) {
        write.csv(fileInput(), file)
      }
    )
    
    
    
    
   ####
    
    
   
  
    
    output$gsUI <- renderUI({
      #inputPanel(
      fluidRow(
        #      column(12,
        selectInput("gsComp", label = "Comparison:",
                    choices = names(ns()$gsRes)),
        #    ),
        #    column(12,
        
        sliderInput("gsJac", label = "Clustering Threshold (Jaccard Index):",
                    min = 0, max = 1, value = 0.5, step = 0.05),
        #    ),
        #    column(12,
        numericInput("gsQthresh", label = "q-value threshold (must be between minumum q and 1):", value = 1, step = 0.05),
        numericInput("gsClust", label = "Cluster to plot:", value = 1, step = 1)
        #    )
      )
      #)
    })
    
    leadingEdge <- reactive({
      fgseaToLEdge(ns()$gsRes, 
                   cutoff.type = "padj",
                   cutoff = input$gsQthresh)
    })
    
    groupedGenesets <- reactive({    
      groupFGSEA(ns()$gsRes[[input$gsComp]],
                 leadingEdge()[[input$gsComp]],
                 join.threshold = (1-input$gsJac),
                 returns = "signif")
    })
    
   
    
    
    
    
    ##########
    output$gsTab <- renderDT({
      #### q value threshold must be less than 1
      #### q value must be greater than or equal to the lowest value in the column
      #### Be able to print a minimum q
      
      datatable(groupedGenesets()[,1:9],
                rownames = FALSE,
                options = list(autoWidth = TRUE,
                               columnDefs = list(list(width = '200px', targets = "_all")))) %>%
        formatStyle('Cluster.Max', target = 'row',
                    color = styleEqual(c("", "x"), c('grey', 'black')),
                    fontSize = '12px')
    }, rownames = FALSE
    )
    
    output$gsHM <- renderPlotly({
      plotlyHeatmap(ns(), groupedGenesets()[,1:9], leadingEdge(), input$gsClust, input$gsComp, input$gsDir)
    })
    ############
    
    
    
    
    
    
    hide(id="jobSetup1")
    hide(id="jobSetup2")
    hide(id="jobSetup3")
    hide(id="jobSetup4")
    
    hide(id = "nanoStringData")
    hide(id ="sampleGroupClassifiers")
    hide(id = "geneSetData")
    
    hide(id="qResults1")
    hide(id="aResults1")
    
    observeEvent(input$setup, {
      toggle(id = "jobSetup1")
      toggle(id =  "jobSetup2")
      toggle(id =  "jobSetup3")
      toggle(id =  "jobSetup4")
      toggle(id = "nanoStringData")
      toggle(id = "sampleGroupClassifiers")
      toggle(id = "geneSetData")
      
    })
    observeEvent(input$qcResults, {
      toggle(id = "qResults1")
    })
    observeEvent(input$anaResults, {
      toggle(id = "aResults1")
    })
    
  }
)