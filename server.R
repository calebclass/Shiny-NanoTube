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
      if (!input$phenModel) {
        selectInput(inputId = "phenCol",
                    label = "Group Column",
                    choices = colnames(sample_info()))
      }
    })
    
    output$basePhen_input <- renderUI({
      req(sample_info(), input$phenCol)
      if (!input$phenModel) {
        selectInput(inputId = "basePhen",
                    label = "Base Group",
                    choices = unique(sample_info()[,input$phenCol]))
      }
    })
    
    phenCol <- reactive({
      if (input$phenModel) {
        NULL
      } else {
        input$phenCol
      }
    })
    
    merged_info <- eventReactive(input$check, {
      req(input$expr, input$phen)
      nanostringData <- processNanostringData(input$expr$datapath,
                                              sampleTab = input$phen$datapath,
                                              groupCol = phenCol(),
                                              normalization = "none",
                                              includeQC = FALSE,
                                              output.format = "list")
    })
    
    output$merged_info <- renderDataTable({
      req(merged_info())
      
      if (input$phenModel) {
        tab1 <- data.frame(Filename = colnames(merged_info()$exprs))
      } else {
        tab1 <- data.frame(Filename = colnames(merged_info()$exprs), 
                           Group = merged_info()$groups)
      }
      
      checkTable <- datatable(cbind(tab1,
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
                                                groupCol = phenCol(),
                                                normalization = input$normMethod,
                                                bgType = "t.test", bgPVal = input$bgP,
                                                housekeeping = hk.genes,
                                                n.unwanted = input$nUnwanted,
                                                RUVg.drop = input$RUVgDrop,
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
                                                groupCol = phenCol(),
                                                bgType = "t.test", bgPVal = input$bgP,
                                                housekeeping = hk.genes,
                                                includeQC = FALSE,
                                                output.format = "list")
        

        
        colnames(nanostringData) <- 
          colnames(nanostringDataBG$exprs.raw) <-
          names(nanostringDataBG$pc.scalefactors) <-
          names(nanostringDataBG$hk.scalefactors) <-
          rownames(nanostringDataBG$qc) <- 
          rownames(nanostringDataBG$bg.stats) <-
          gsub(".*\\/|\\.RCC", "", colnames(nanostringData))
        
        incProgress(1/4, detail = "Analyzing Diff. Expr.")
        
        if (input$phenModel) {
          # Design matrix has had 2 extra columns added -- these are removed here.
          base.group <- "Intercept"
          design.mat <- pData(nanostringData)[,2:(ncol(pData(nanostringData))-1)]
          limmaResults <- runLimmaAnalysis(nanostringData, design = design.mat)
        } else {
          base.group <- input$basePhen
          limmaResults <- runLimmaAnalysis(nanostringData, groups = NULL, base.group)
        }
        
        ns <- list(dat = nanostringData,
                   dat.list = nanostringDataBG,
                   deRes = limmaResults,
                   base.group = base.group,
                   gene.stats = nanostringDataBG2$gene.stats)
        
        if (!is.null(input$gsDb$datapath)) {
          incProgress(1/6, detail = "Analyzing Gene Sets")
          ns$gsRes <- limmaToFGSEA(limmaResults, input$gsDb$datapath,
                                   min.set = input$minSize)
        } else if (input$gsReactome) {
          incProgress(1/6, detail = "Analyzing Gene Sets")
          ns$gsRes <- limmaToFGSEA(limmaResults, "data/ReactomePathways.gmt",
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
    hkQC <- reactive({ housekeepingQC(ns()$dat.list, plotType = input$boxplotType) })
    boxHeight <- function() {
      hkQC()$pltHeight
    }
    pcaPlot <- reactive({ plotPCA(ns()) })
    deResults <- reactive({ deRes(ns(), input$summaryQ) })
    ####
    
    canoPlot <- reactive({deVolcanoInt(limmaResults = ns()$deRes,
                                       plotContrast = input$volComp) +
        geom_hline(yintercept =  -log10(input$volcanoHorLineInput), linetype =  "dashed", colour = 'darkred') +
        geom_vline(xintercept = input$volcanoVertLineInput, linetype = "dashed", colour = "darkred") +
        geom_vline(xintercept = -input$volcanoVertLineInput, linetype = "dashed", colour = "darkred")})
    
    ###
    
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
    output$normPlot1 <- renderPlot({hkQC()$plt1},
                                 height = boxHeight)
    output$normPlot2 <- renderPlot({hkQC()$plt2},
                                 height = boxHeight)
    output$hkPlot1 <- renderPlot({hkQC()$j1},
                                   height = boxHeight)
    output$hkPlot2 <- renderPlot({hkQC()$j2},
                                   height = boxHeight)
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
      filename = function() {"DE_table.csv"},
      content = function(file) {
        write.csv(deResults()$de$x$data, file, 
                  row.names = FALSE)
      }
    )
    
    output$GSdownload <- downloadHandler(
      filename = function() {"GS_table.csv"},
      content = function(file) {
        write.csv(groupedGenesets()[,1:9], file, 
                  row.names = FALSE)
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
        numericInput("gsClust", label = "Cluster for heatmap:", value = 1, step = 1)
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
      genesetsOut <- groupedGenesets()
      genesetsOut[,2:6] <- signif(genesetsOut[,2:6], digits = 3)
      
      datatable(genesetsOut[,1:9],
                rownames = FALSE,
                options = list(autoWidth = FALSE,
                               scrollX = TRUE,
                               columnDefs = list(list(width = '200px', targets = "_all"),
                                                 list(className = 'dt-center', targets = 1:8)))) %>%
        formatStyle('Cluster.Max', target = 'row',
                    color = styleEqual(c("", "x"), c('grey', 'black')),
                    fontSize = '12px')  # Round 2:6 (1:5 in python)
    }, rownames = FALSE
    )
    
    output$gsHM <- renderPlotly({
      req(groupedGenesets())
      plotlyHeatmap(ns(), groupedGenesets()[,1:9], leadingEdge(), input$gsClust, input$gsComp, input$gsDir)
    })
    ############
    

  
    
  }
)