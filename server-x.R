# server.R
library(iDINGO)
library(ggplot2)
library(plyr)
library(stringr)
source("helpers.R")

shinyServer(
  
  function(input, output, session) {
    options(shiny.maxRequestSize = 30*1024^2)
    
    observeEvent(input$run, {
      updateNavbarPage(session, "master",
                       selected = "Visualize Results"
      )
    })
    
    platforms <- eventReactive(input$run, {
      return(c(paste("1", input$exprName1, sep=ifelse(input$exprName1 == "", yes = "", no = "-")),
               paste("2", input$exprName2, sep=ifelse(input$exprName2 == "", yes = "", no = "-")),
               paste("3", input$exprName3, sep=ifelse(input$exprName3 == "", yes = "", no = "-"))
      ))
    })
    
    fit <- eventReactive(input$run, {
      
      # Read in datasets
      dataset <- read.table(file = input$expr$datapath, sep = "", row.names = 1)
      
      if (class(input$expr2$datapath) == "character") {
        dataset2 <- read.compare.datasets(file = input$expr2$datapath, prev.dataset = dataset)
      } else {
        dataset2 <- NULL
        Y2 <- NULL
      }
      
      if (class(input$expr3$datapath) == "character") {
        dataset3 <- read.compare.datasets(file = input$expr3$datapath, prev.dataset = dataset)
      } else {
        dataset3 <- NULL
        Y3 <- NULL
      }
      
      phenotypes <- read.table(file = input$phen$datapath, sep = "")
      
      # Remove rows with zero variance
      dataset <- dataset[apply(dataset, 1, var) > 0,]
      try(dataset2 <- dataset2[apply(dataset2, 1, var) > 0,], silent = TRUE)
      try(dataset3 <- dataset3[apply(dataset3, 1, var) > 0,], silent = TRUE)
      
      # Set up DINGO
      Y1 <- t(dataset)
      try(Y2 <- t(dataset2), silent = TRUE)
      try(Y3 <- t(dataset3), silent = TRUE)
      x <- as.numeric(unlist(phenotypes))
      
      fit <- try(idingo(dat=Y1, dat2=Y2, dat3=Y3, x=x, plats=platforms(), diff.score=TRUE, B=input$numBoot, verbose = TRUE, cores = input$numCore))
      
      return(fit)
      
    })
    
    cor.data <- reactive({ makeCorTable(fit(), input$threshType, input$thresh) })
    
    output$network <- renderVisNetwork({
      plotNetwork(fit(), threshold = input$thresh, thresh.type = input$threshType, layout = input$layout)
    })
    
    output$table <- renderTable({
      signif.pairs <- cor.data()[cor.data()$ds.col != "grey",]
      hub.count <- count(c(as.character(signif.pairs$gene1),
                           as.character(signif.pairs$gene2)))
      hub.count[,1] <- gsub(".*\\|\\|", "", hub.count[,1])
      makeHubTable(hub.count, input$hubs)
      
    }, digits = 3)
    
    output$correlations <- renderPlot({
      
      if (sum(cor.data()$ds.col != "grey") > 0) {
        ggplot(data = cor.data(), mapping = aes(x = R2, y = R1, col = ds.col)) +
          geom_abline(slope = 1, intercept = 0, colour = "grey") +
          geom_point(data = subset(cor.data(), ds.col == "grey"), aes(x = R2, y = R1, col = ds.col)) +
          geom_point(data = subset(cor.data(), ds.col != "grey"), aes(x = R2, y = R1, col = ds.col)) +
          scale_colour_manual(values = c("blue", "grey", "red")) +
          geom_hline(yintercept = 0, linetype = "longdash", colour = "black") +
          geom_vline(xintercept = 0, linetype = "longdash", colour = "black") +
          xlab("Partial Correlation (Group 1)") + ylab("Partial Correlation (Group 2)") +
          theme_bw() + theme(legend.position = "none")
      }
    })
    
    output$summary1 <- renderText({
      paste("Total pairs possible:", length(fit()$R1))
    })
    
    output$summary2 <- renderText({
      signif.num <- dim(cor.data()[cor.data()$ds.col != "grey",])[1]
      paste("Pairs meeting threshold: ", signif.num, " (",
            signif(100 * signif.num / length(fit()$R1), 3), "%)", sep = "")
    })
    
    output$download <- downloadHandler(
      filename = function() { paste(input$dlname, ".", input$fileformat, sep="") },
      content = function(file) {
        if (input$fileformat == "rds") saveRDS(fit(), file)
        else if (input$fileformat == "txt") write.table(cor.data()[,1:6], file, sep = "\t",
                                                        row.names = FALSE, quote = FALSE)
        else write.csv(cor.data()[,1:6], file, row.names = FALSE, quote = FALSE)
      }
    )
    
  }
)
