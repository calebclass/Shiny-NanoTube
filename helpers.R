# Add logo to dashboard header:
header <- dashboardHeader()
anchor <- tags$a(href='',
                 tags$img(src='NanoTube-LogoBlue2.jpg', height='50', width='144'))

header$children[[2]]$children <- tags$div(
  anchor,
  class = 'name')

######################

summarizeNegQC <- function(ns) {
  negSummary <- data.frame(nms = c("Pass", "Removed"),
                           counts = c(paste0(sum(ns$gene.stats$pass), " Targets"),
                                      paste0(sum(!ns$gene.stats$pass), " Targets")))
  
  return(negSummary)
}

prepNegGenes <- function(ns) {
  geneTab <- ns$gene.stats
  geneTab[,1:2] <- signif(geneTab[,1:2], digits = 2)
  
  return(datatable(geneTab,
                   options = list(
                     columnDefs = list(list(className = 'dt-center', targets = 0:2))
                   )))
}

housekeepingQC <- function(ns, plotType = "RLE") {
  hk.tab <- data.frame(Sample = names(ns$dat.list$hk.scalefactors),
                       `Scale Factor` = round(ns$dat.list$hk.scalefactors, 2))
  
  if (plotType == "RLE") {
    
    boxplot.dat <- log2(ns$dat.list$exprs.raw[grep("endogenous", ns$dat.list$dict.raw$CodeClass, ignore.case=TRUE),]+0.5)
    box1 <- ruv::ruv_rle(t(boxplot.dat), ylim = c(-2, 2)) +
      ylab("Relative log expression") +
      theme(axis.text.x = element_text(size = 12),
            axis.title.x = element_text(size = 12),
            axis.text.y = element_text(size = 14)) +
      ylab("Relative log expression") + 
      scale_x_discrete(limits = colnames(ns$dat.list$exprs.raw)) +
      coord_flip()
    
    boxplot.dat <- log2(exprs(ns$dat)[grep("endogenous", fData(ns$dat)$CodeClass, ignore.case=TRUE),]+0.5)
    box2 <- ruv::ruv_rle(t(boxplot.dat), ylim = c(-2, 2)) +
      theme(axis.text.x = element_text(size = 12),
            axis.title.x = element_text(size = 12),
            axis.text.y = element_text(size = 14)) +
      ylab("Relative log expression") + 
      scale_x_discrete(limits = colnames(exprs(ns$dat))) +
      coord_flip()
  
  } else {
    
    boxplot.dat <- as.data.frame(log2(ns$dat.list$exprs.raw+0.5))
    boxplot.dat$CodeClass <- ns$dat.list$dict.raw$CodeClass
    boxplot.df <- reshape::melt(boxplot.dat, "CodeClass")
    
    box1 <- ggplot() +
      geom_boxplot(data=boxplot.df[grep("endogenous", boxplot.df$CodeClass, ignore.case=TRUE),], 
                   aes(x=variable, y=value), fill="white") + 
      theme_bw() + xlab("") + ylab("log2(counts+0.5)") +
      theme(axis.text.x = element_text(size = 12),
            axis.title.x = element_text(size = 12),
            axis.text.y = element_text(size = 14)) +
      coord_flip()
    
    
    boxplot.dat <- log2(exprs(ns$dat)[grep("endogenous", fData(ns$dat)$CodeClass, ignore.case=TRUE),]+0.5)
    boxplot.df <- reshape::melt(boxplot.dat)
    
    box2 <- ggplot() +
      geom_boxplot(data=boxplot.df, aes(x=X2, y=value), fill="white") + 
      theme_bw() + xlab("") + ylab("log2(counts+0.5)") +
      theme(axis.text.x = element_text(size = 12),
            axis.title.x = element_text(size = 12),
            axis.text.y = element_text(size = 14)) +
      coord_flip()
    
  }
  
  # Housekeeping gene plot (raw data)
  hk.dat <- log2(ns$dat.list$exprs.raw[grep("housekeep", ns$dat.list$dict.raw$CodeClass, ignore.case=TRUE),]+0.5)
  hk.medians <- apply(hk.dat, 1, median)
  hk.dat <- hk.dat - hk.medians
  
  rownames(hk.dat) <- ns$dat.list$dict.raw$Name[grep("housekeep", ns$dat.list$dict.raw$CodeClass, ignore.case=TRUE)]
  hk.plot <- reshape::melt(hk.dat)
  jitter1 <- ggplot(data = hk.plot, aes(x = X1, y = value)) +
    geom_boxplot() +
    geom_jitter(height = 0, width = 0.25) +
    theme_bw() + xlab("") + ylab("RLE") +
    theme(axis.text.x = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          axis.text.y = element_text(size = 14)) +
    coord_flip()
    
  # Housekeeping gene plot (normalized data)
  hk.dat <- log2(ns$dat.list$exprs[grep("housekeep", ns$dat.list$dict$CodeClass, ignore.case=TRUE),]+0.5)
  hk.medians <- apply(hk.dat, 1, median)
  hk.dat <- hk.dat - hk.medians
  
  rownames(hk.dat) <- ns$dat.list$dict$Name[grep("housekeep", ns$dat.list$dict$CodeClass, ignore.case=TRUE)]
  hk.plot <- reshape::melt(hk.dat)
  jitter2 <- ggplot(data = hk.plot, aes(x = X1, y = value)) +
    geom_boxplot() +
    geom_jitter(height = 0, width = 0.25) +
    theme_bw() + xlab("") + ylab("RLE") +
    theme(axis.text.x = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          axis.text.y = element_text(size = 14)) +
    coord_flip()
  
  return(list(tab = hk.tab,
              plt1 = box1,
              plt2 = box2,
              j1 = jitter1,
              j2 = jitter2,
              pltHeight = 150 + nrow(hk.tab) * 15,
              jitterHeight = 150 + length(grep("housekeep", ns$dat.list$dict$CodeClass, ignore.case=TRUE)) * 15))
}


plotPCA <- function(ns) {
  pca.dat <- log2(exprs(ns$dat)[grepl("endogenous", fData(ns$dat)$CodeClass, ignore.case=TRUE) &
                                 rowSums(exprs(ns$dat) == 0) == 0,]+0.5)
  
  gp <- pData(ns$deRes$eset)$group
  
  pca <- prcomp(t(pca.dat),
                center = TRUE, scale = TRUE)
  pca.ly <- as.data.frame(pca$x[,1:2])
  pca.ly$sample <- row.names(pca$x)
  if (!is.null(gp)) pca.ly$group <- gp
  
  perc.var <- pca$sdev^2 / sum(pca$sdev^2) * 100
  
  layout.x <- layout.y <- list(
    showline = TRUE,
    showticklabels = TRUE,
    showgrid = FALSE,
    zeroline = FALSE,
    linecolor = toRGB("black"),
    linewidth = 1
  )
  
  layout.x$title <- paste0("PC1 (", signif(perc.var[1], 2), "% var)")
  layout.y$title <- paste0("PC2 (", signif(perc.var[2], 2), "% var)")
  
  if (is.null(gp)) {
    plt <- plot_ly(data = pca.ly, x = ~PC1, y = ~PC2,
                   text = ~paste0("Sample: ", sample),
                   type = "scatter", mode = "markers") %>%
      layout(xaxis = layout.x, yaxis = layout.y)
  } else {
    plt <- plot_ly(data = pca.ly, x = ~PC1, y = ~PC2,
                   text = ~paste0("Sample: ", sample), color = ~group,
                   type = "scatter", mode = "markers") %>%
      layout(xaxis = layout.x, yaxis = layout.y)
  }
  
  return(plt)
}


deRes <- function(ns, summaryQ) {
  diffExpr.tab <- rbind(colSums(ns$deRes$q.value < summaryQ & ns$deRes$coefficients > 0),
                        colSums(ns$deRes$q.value < summaryQ & ns$deRes$coefficients < 0))
  diffExpr.tab <- sapply(as.data.frame(diffExpr.tab[,!(colnames(diffExpr.tab) %in% c("Intercept", "(Intercept)"))]),
                         as, "integer")
  
  
  if (ncol(diffExpr.tab) == 1) {
    colnames(diffExpr.tab) <- ""
    gp <- colnames(ns$deRes$q.value)[!(colnames(ns$deRes$q.value) %in% c("Intercept", "(Intercept)"))]
    rownames(diffExpr.tab) <- c(paste0("logFC > 0 for ", gp),
                                 paste0("logFC < 0 for ", gp))
  } else {
    rownames(diffExpr.tab) <- c("logFC > 0", 
                                "logFC < 0")
  }
  
  
  #knitr::kable(diffExpr.tab[,-1], padding=2, format="html", 
  #             col.names = ifelse(ncol(diffExpr.tab)==2, yes = "", no = NA)) %>%
  #  kableExtra::kable_styling(full_width = F)
  

  diffExpr.full <- as.data.frame(makeDiffExprFile(ns$deRes, returns = "stats"))
  
  brks.fc <- seq(-2, 2, 4/100)
  ln <- length(brks.fc) + 1
  cols.fc <- paste0("rgb(",
                    c(round(seq(40, 255, length.out = ln/2), 0), rep(255, times = ln/2)), ",",
                    c(round(seq(40, 255, length.out = ln/2), 0), round(seq(255, 40, length.out = ln/2), 0)), ",",
                    c(rep(255, times = ln/2), round(seq(255, 40, length.out = ln/2), 0)), ")")
  brks.qv <- seq(0, 0.1, 0.1/100)
  cols.qv <- paste0("rgb(255,", round(seq(40, 255, length.out = length(brks.qv) + 1), 0), ",255)")
  dtable <- datatable(diffExpr.full, 
                      options = list(columnDefs = list(list(className = 'dt-center', targets = "_all")),
                                     autoWidth = TRUE, scrollX = TRUE,
                                     #columnDefs = list(list(width = "125px", targets = "_all")),
                                     dom = 'tpB',
                                     lengthMenu = list(c(5, 15,-1), c('5', '15', 'All')),
                                     pageLength = 10)) %>%
    formatStyle(names(diffExpr.full)[substr(names(diffExpr.full), 1, 6) == "Log2FC"], 
                backgroundColor = styleInterval(brks.fc, cols.fc)) %>%
    formatStyle(names(diffExpr.full)[substr(names(diffExpr.full), 1, 5) == "q-val"], 
                backgroundColor = styleInterval(brks.qv, cols.qv))

  return(list(summary = diffExpr.tab,
              de = dtable))
  
}


plotlyHeatmap <- function(ns, groupedGenesets, leadingEdge, gsClust, gsComp, gsDir) {
  
  if (is.null(ns$deRes$eset$group)) {
    dat.scaled <- exprs(ns$deRes$eset) -
      apply(exprs(ns$deRes$eset), 1, median)
  } else {
    dat.scaled <- exprs(ns$deRes$eset) -
      apply(exprs(ns$deRes$eset)[,ns$deRes$eset$group == 
                                   levels(ns$deRes$eset$group)[1]], 
            1, median)
  }
  
  # Force all data between -3 & 3 (log-scale)
  hm.max <- 3
  dat.scaled[dat.scaled > hm.max] <- hm.max
  dat.scaled[dat.scaled < -hm.max] <- -hm.max
  
  # Colors for leading edge heatmap
  acols <- c("grey95", "black")
  names(acols) <- c(0, 1)
  
  # Set up axes
  ay1 <- ay2 <- ay <- ax1 <- ax <- list(title = "", ticks = "")
  ax$tickangle <- ax1$tickangle <- 20
  ax$tickfont <- list(size = 8)
  
  # Include genes from selected gene set cluster (gsClust)
  subClust <- rownames(groupedGenesets)[groupedGenesets$Cluster == gsClust]
  annot_row <- leadingEdge[[gsComp]]
  annot_row <- annot_row[,subClust]
  
  if (length(subClust) > 1) {
    annot_row <- as.data.frame(annot_row)[rowSums(annot_row) >= 1,]
    annot_row.df <- sapply(annot_row, as.numeric)
    rownames(annot_row.df) <- rownames(annot_row)
    
    genesKeep <- rownames(dat.scaled) %in% rownames(annot_row.df)
  } else {
    genesKeep <- rownames(dat.scaled) %in% names(annot_row)[annot_row==1]
  }
  
  #Cluster genes
  genes <- labels(as.dendrogram(hclust(dist(dat.scaled[genesKeep,]), method="ward.D")))
  
  dat.hm <- dat.scaled[genes,]
  
  # Replace any spaces with periods in group names (automatically done by other functions)
  sample.groups <- gsub(" ", ".", ns$dat$groups)
  base.group <- gsub(" ", ".", ns$base.group)
  
  title.y <- 1.15
  margin.plots <- list(b = 100, t = 80)
  
  # Split to Base & Comparison group
  if (is.null(ns$deRes$eset$group)) {
    base.group <- ""
    dat.hmB <- dat.hm
    samp_namesB <- colnames(dat.scaled)
  } else {
    dat.hmB <- dat.hm[, sample.groups == base.group]
    dat.hmC <- dat.hm[, sample.groups == gsComp]
    samp_namesB <- colnames(dat.scaled)[sample.groups == base.group]
  }

  dat.hmB <- reshape::melt(dat.hmB)
  dat.hmB$X1 <- factor(dat.hmB$X1, levels = genes)
  dat.hmB$X2 <- factor(dat.hmB$X2, levels = samp_namesB)
  colnames(dat.hmB)[3] <- "expression"
  
  if (length(subClust) > 1) {
    annot_row.hm <- reshape::melt(annot_row.df)
    annot_row.hm$X1 <- factor(annot_row.hm$X1, levels = genes)
    annot_row.hm$X2 <- factor(annot_row.hm$X2, levels = colnames(annot_row.df))
    ay1$showticklabels <- FALSE
  }
  
  # Heatmap for Base Group (or all data if design matrix was used)
  pHM_Base <- plot_ly(data = dat.hmB, x = ~X2, y = ~X1, z = ~expression, xgap = 2, ygap = 2, 
                 type = "heatmap", colors = colorRamp(c("blue", "white", "red")), alpha = 1,
                 zauto = FALSE, zmax = hm.max, zmin = -hm.max) %>%
    layout(xaxis = ax1, yaxis = ay1, margin = margin.plots) %>%
#    add_annotations(text = gsub("\\.", " ", base.group), 
#                    x = dat.hmB$X2[1], y = as.numeric(dat.hmB$X1[nrow(dat.hmB)]) + title.offset,
#                    showarrow = FALSE)
    add_annotations(text = gsub("\\.", " ", base.group),
                    x = dat.hmB$X2[1], y = title.y,
                    xref = "x", yref = "paper", showarrow = FALSE, 
                    font = list(size = 20), xanchor = "left")
  
  # Hetmap for Comparison Group
  if (!is.null(ns$deRes$eset$group)) {
    dat.hmC <- reshape::melt(dat.hmC)
    dat.hmC$X1 <- factor(dat.hmC$X1, levels = genes)
    dat.hmC$X2 <- factor(dat.hmC$X2, levels = colnames(dat.scaled)[sample.groups == gsComp])
    colnames(dat.hmC)[3] <- "expression"
    ay2$showticklabels <- FALSE
    
    rel.width.B <- sum(sample.groups == base.group) / sum(sample.groups %in% c(base.group, gsComp))
    rel.width.C <- 1 - rel.width.B
    
    # Heatmap for Comparison Group
    pHM_Comp <- plot_ly(data = dat.hmC, x = ~X2, y = ~X1, z = ~expression, xgap = 2, ygap = 2, 
                        type = "heatmap", colors = colorRamp(c("blue", "white", "red")), alpha = 1,
                        zauto = FALSE, zmax = hm.max, zmin = -hm.max, showscale = FALSE) %>%
      layout(xaxis = ax1, yaxis = ay2, margin = margin.plots) %>%
      add_annotations(text = gsub("\\.", " ", gsComp),
                      x = dat.hmC$X2[1], y = title.y,
                      xref = "x", yref = "paper", showarrow = FALSE,
                      font = list(size = 20), xanchor = "left")
  }
  
  # Leading Edge Heatmap (If more than 1 gene set included)
  if (length(subClust) > 1) {
    
    pCB <- plot_ly(data = annot_row.hm, x = ~X2, y = ~X1, z = ~value, xgap = 2, ygap = 2,
                   type = "heatmap", colors = colorRamp(c("grey95", "black")), alpha = 1, 
                   showscale = FALSE) %>%
      layout(xaxis = ax, yaxis = ay, margin = margin.plots)
#      add_annotations(text = " ",
#                      x = 0, y = title.y,
#                      xref = "paper", yref = "paper", showarrow = FALSE)
    
    if (is.null(ns$deRes$eset$group)) {
      rel.width.ledge <- 0.8*ncol(annot_row) / (ncol(annot_row) + ncol(dat.hmB))
      return(subplot(pCB, pHM_Base, nrows = 1, widths = c(rel.width.ledge, 1) / (1+rel.width.ledge)))
    } else {
      rel.width.ledge <- 0.8*ncol(annot_row) / (ncol(annot_row) + 
                                                  sum(sample.groups %in% c(base.group, gsComp))) 
      return(subplot(pCB, pHM_Base, pHM_Comp, nrows = 1, widths = c(rel.width.ledge, rel.width.B, rel.width.C) / (1+rel.width.ledge)))
    }
  } else {
    if (is.null(ns$deRes$eset$group)) {
      return(pHM_Base)
    } else {
      return(subplot(pHM_Base, pHM_Comp, nrows = 1, widths = c(rel.width.B, rel.width.C)))
    }
    
  }
  
}


deVolcanoInt <- function(limmaResults, 
                       plotContrast = NULL, y.var = c("p.value", "q.value")) {
  
  # Bind local variables
  log2FC <- log10p <- NULL
  
  # Identify contrast if not provided
  if (is.null(plotContrast)) {
    plotContrast <-
      colnames(limmaResults$coefficients)[
        which(!(colnames(limmaResults) %in% 
                  c("Intercept", "(Intercept)")))[1]]
    
    cat("\n'plotContrast' not provided, setting it to", plotContrast, "\n")
  }
  
  # Set up data frame for plot
  df <- data.frame(log2FC = limmaResults$coefficients[,plotContrast],
                   log10p = -log10(limmaResults[[y.var[1]]][,plotContrast]),
                   name = limmaResults$genes$Name)
  
  
  plt <- ggplot(df, aes(x = log2FC, y = log10p, label = name)) +
    geom_point() +
    xlab("log2(Fold Change)") +
    ylab(paste0("-log10(", substr(y.var[1], 1, 1), ")")) +
    theme_bw()
  
  
  
  return(plt)
  
}

prepPosOutputs <- function(posQC) {
  posQC$plotly <- ggplotly(posQC$plt,
                           height = nrow(posQC$tab)*90) 
  
  posQC$DT <- posQC$tab
  posQC$DT$Scale.Factor <- signif(posQC$DT$Scale.Factor, digits = 3)
  posQC$DT$R.squared <- signif(posQC$DT$R.squared, digits = 3)
  posQC$DT <- datatable(posQC$DT, rownames = FALSE)
  
  return(posQC)
}