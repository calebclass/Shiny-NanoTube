shinyServer(
  function(input, output, session) {
    options(shiny.maxRequestSize=100*1024^2)
    
    sample_info <- reactive({
      req(input$phen)
      read.csv(input$phen$datapath, header = TRUE,
               sep = ",")
    })
    
    output$phenCol_input <- renderUI({
      req(sample_info())
      selectInput(inputId = "phenCol",
                  label = "Group Column",
                  choices = colnames(sample_info()))
    })
    
    output$basePhen_input <- renderUI({
      req(sample_info(), input$phenCol)
      selectInput(inputId = "basePhen",
                  label = "Base Group",
                  choices = unique(sample_info()[,input$phenCol]))
    })
    
    merged_info <- eventReactive(input$check, {
      req(input$expr, input$phen)
      nanostringData <- processNanostringData(input$expr$datapath,
                                              sampleTab = input$phen$datapath,
                                              groupCol = input$phenCol,
                                              normalization = "none",
                                              includeQC = FALSE,
                                              output.format = "list")
    })
    
    output$merged_info <- renderDataTable({
      req(merged_info())
      
      checkTable <- datatable(cbind(data.frame(Filename = colnames(merged_info()$exprs), 
                                     Group = merged_info()$groups),
                          merged_info()$samples),
                          rownames = FALSE,
                          options = list(
                            autoWidth = TRUE, scrollX = TRUE,
                            columnDefs = list(list(
                              width = "125px", targets = "_all"
                            )),
                            dom = 'tpB',
                            lengthMenu = list(c(5, 15,-1), c('5', '15', 'All')),
                            pageLength = 10
                          ))  # thanks to https://stackoverflow.com/questions/57946206/how-to-resize-a-datatable-in-order-to-fit-it-in-a-box-for-shinydashboard
    })
    
#    observeEvent(input$run, {
#      updateNavbarPage(session, "master",
#                       selected = "QC Results"
#      )
#    })
    
    ns <- eventReactive(input$run, {
      
      withProgress(message = "Please wait", value = 0, {
        
        if (input$hk == "") {
          hk.genes <- NULL
        } else {
          hk.genes <- gsub(" ", "", strsplit(input$hk, split=",")[[1]])
        }
        
        incProgress(1/4, detail = "Processing Data")
        
        nanostringData <- processNanostringData(input$expr$datapath,
                                                sampleTab = input$phen$datapath,
                                                groupCol = input$phenCol,
                                                bgType = "t.test", bgPVal = input$bgP,
                                                housekeeping = hk.genes,
                                                includeQC = FALSE)
        
        incProgress(1/4, detail = "Normalizing Data")
        
        nanostringDataBG <- processNanostringData(input$expr$datapath,
                                                  bgType = "threshold", bgThreshold = 2, bgProportion = 0.5,
                                                  housekeeping = hk.genes,
                                                  includeQC = TRUE,
                                                  output.format = "list")
        
        # Need to run this again, just to get Gene statistics vs. Background.
        # This will be updated in NanoTube R package
        nanostringDataBG2 <- processNanostringData(input$expr$datapath,
                                                sampleTab = input$phen$datapath,
                                                groupCol = input$phenCol,
                                                bgType = "t.test", bgPVal = input$bgP,
                                                housekeeping = hk.genes,
                                                includeQC = FALSE,
                                                output.format = "list")
        
        
        #      file_input <- input$expr$datapath
        #      file.extension <- substr(file_input[1], 
        #                               (nchar(file_input[1])-3), nchar(file_input[1]))
        #      if(file.extension %in% c(".zip",".ZIP")){
        #        file_input <- unzip(file_input)
        #        file_input <-read_merge_rcc(file_input)
        #        nanoTableData <- file_input
        #        output$nanoTable <- renderTable(file_input)
        #      }
        
        colnames(nanostringData) <- 
          colnames(nanostringDataBG$exprs.raw) <-
          names(nanostringDataBG$pc.scalefactors) <-
          names(nanostringDataBG$hk.scalefactors) <-
          rownames(nanostringDataBG$qc) <- 
          rownames(nanostringDataBG$bg.stats) <-
          gsub(".*\\/|\\.RCC", "", colnames(nanostringData))
        
        #      if (is.null(input$phen$datapath)) {
        #        groups <- gsub("_.*", "", colnames(nanostringData))
        #      } else {
        #        groups <- as.character(read.table(file = input$phen$datapath, sep = "", as.is = TRUE))
        #      }
        
        #      if (input$basePhen == "") {
        #        base.group <- groups[1]
        #      } else {
        base.group <- input$basePhen
        #      }
        
        incProgress(1/4, detail = "Analyzing Diff. Expr.")
        
        limmaResults <- runLimmaAnalysis(nanostringData, NULL, base.group)
        
        ns <- list(dat = nanostringData,
                   dat.list = nanostringDataBG,
                   deRes = limmaResults,
                   base.group = base.group,
                   gene.stats = nanostringDataBG2$gene.stats)
        
        if (!is.null(input$gsDb$datapath)) {
          incProgress(1/6, detail = "Analyzing Gene Sets")
          ns$gsRes <- limmaToFGSEA(limmaResults, input$gsDb$datapath,
                                   min.set = input$minSize)
        }
        
        return(ns)
        
      })
    })
    
    output$numSamps <- renderText({ 
      req(ns())
      paste0("All done! Analyzed ", ncol(ns()$dat), " samples.\nUse the other tabs to view results.") 
      })
    
    posQC <- reactive({ prepPosOutputs(positiveQC(ns()$dat.list)) })
    negQC <- reactive({ negativeQC(ns()$dat.list, interactive.plot = FALSE) })
    hkQC <- reactive({ housekeepingQC(ns()$dat.list) })
    pcaPlot <- reactive({ plotPCA(ns()) })
    deResults <- reactive({ deRes(ns(), input$summaryQ) })
    ####
    
    canoPlot <- reactive({deVolcanoInt(limmaResults = ns()$deRes,
                                       plotContrast = input$volComp) +
        geom_hline(yintercept =  -log10(input$volcanoHorLineInput), linetype =  "dashed", colour = 'darkred') +
        geom_vline(xintercept = input$volcanoVertLineInput, linetype = "dashed", colour = "darkred") +
        geom_vline(xintercept = -input$volcanoVertLineInput, linetype = "dashed", colour = "darkred")})
    
    ###
    
    #####

    output$posTab <- renderDataTable({ posQC()$DT })
    output$posPlot <- renderPlotly({ posQC()$plotly })
    output$negTab <- renderDataTable({datatable(negQC()$tab, rownames = TRUE,
                                                options = list(
                                                  columnDefs = list(list(className = 'dt-center', targets = 0:5))
                                                ))})
    output$negGenes <- renderDataTable({prepNegGenes(ns())})
    output$negSummary <- renderTable ({summarizeNegQC(ns())},
                                      rownames = FALSE, colnames = FALSE)
    output$negPlot <- renderPlotly({ggplotly(negQC()$plt,
                                             height = 120 + nrow(negQC()$tab) * 15 )})
    output$hkTab <- renderDataTable({datatable(hkQC()$tab, rownames = FALSE) })
    output$hkPlot1 <- renderPlotly({hkQC()$plt1})
    output$hkPlot2 <- renderPlotly({hkQC()$plt2})
    output$pcaPlot <- renderPlotly({
      req(ns())
      pcaPlot()})
    output$deCounts <- renderTable({
      req(ns())
      deResults()$summary},
      rownames = TRUE, round = 0)
    output$deTab <- renderDT({req(ns())
      deResults()$de}, 
      rownames = TRUE)
    
    ####
    output$DEdownload <- downloadHandler(
      filename = function() {"de.csv"},
      content = function(file) {
        write.csv(deResults(), file)
      }
    )
    
    ###
    output$canoPlot <- renderPlotly({req(ns())
      canoPlot()})
    ###
    output$NANOdownload <- downloadHandler(
      filename = function() {"nanoTable.csv"},
      content = function(file) {
        write.csv(fileInput(), file)
      }
    )
    
    
    
    
   ####
    
    
    output$volUI <- renderUI({
        req(ns())
        selectInput("volComp", label = "Comparison (vs. Base Group):",
                    choices = colnames(ns()$deRes)[!(colnames(ns()$deRes) %in% 
                                                      c("Intercept", "(Intercept)"))]
      )
    })
   
  
    
    output$gsUI <- renderUI({
      #inputPanel(
      req(ns()$gsRes)
      column(width = 10,
        #      column(12,
        selectInput("gsComp", label = "Comparison (vs. Base Group):",
                    choices = names(ns()$gsRes)),
        #    ),
        #    column(12,
        
        sliderInput("gsJac", label = "Clustering Threshold (Jaccard Index):",
                    min = 0, max = 1, value = 0.5, step = 0.05),
        #    ),
        #    column(12,
        numericInput("gsQthresh", label = "q-value threshold:", value = 1, step = 0.05),
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
      tab <- groupFGSEA(ns()$gsRes[[input$gsComp]],
                 leadingEdge()[[input$gsComp]],
                 join.threshold = (1-input$gsJac),
                 returns = "signif")
    })
    
   
    
    
    
    
    ##########
    output$gsTab <- renderDT({
      #### q value threshold must be less than 1
      #### q value must be greater than or equal to the lowest value in the column
      #### Be able to print a minimum q
      req(groupedGenesets())
      
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
      req(groupedGenesets())
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