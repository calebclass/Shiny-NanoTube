housekeepingQC <- function(ns) {
  hk.tab <- data.frame(Sample = names(ns$hk.scalefactors),
                       `Scale Factor` = round(ns$hk.scalefactors, 2))
  
  boxplot.dat <- as.data.frame(log2(ns$exprs.raw+0.5))
  boxplot.dat$CodeClass <- ns$dict.raw$CodeClass
  boxplot.df <- reshape::melt(boxplot.dat, "CodeClass")
  
  b1 <- ggplot() +
    geom_boxplot(data=boxplot.df[boxplot.df$CodeClass == "Endogenous",], 
                 aes(x=variable, y=value), fill="grey70") + 
    theme_classic() + xlab("") + ylab("log2(counts+0.5)") +
    ggtitle("Pre-Normalization Data") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  ggplotly(b1) %>% 
    layout(margin = list(b=90))
  
  boxplot.dat <- log2(ns$exprs[ns$dict$CodeClass == "Endogenous",]+0.5)
  boxplot.df <- reshape::melt(boxplot.dat)
  
#  samp.df <- data.frame(sample = rep(colnames(boxplot.dat),times=2),
#                        scaleFactor = c(rep("Positive Control", times=ncol(boxplot.dat)), rep("Housekeeping", times=ncol(boxplot.dat))),
#                        dat = c(ns$pc.scalefactors,
#                                hk.scaleFactor = ns$hk.scalefactors))
  
  b2 <- ggplot() +
    geom_boxplot(data=boxplot.df, aes(x=X2, y=value), fill="grey70") + 
    theme_classic() + xlab("") + ylab("log2(counts+0.5)") +
    ggtitle("Normalized Data")
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  ggplotly(b2) %>% 
    layout(margin = list(b=90))
  
  return(list(tab = hk.tab,
              plt1 = b1,
              plt2 = b2))
}


plotPCA <- function(ns) {
  pca.dat <- log2(exprs(ns$dat)[fData(ns$dat)$CodeClass == "Endogenous" &
                                 rowSums(exprs(ns$dat) == 0) == 0,]+0.5)
  
  pca <- prcomp(t(pca.dat),
                center = TRUE, scale = TRUE)
  pca.ly <- as.data.frame(pca$x[,1:2])
  pca.ly$sample <- row.names(pca$x)
  pca.ly$group <- pData(ns$deRes$eset)$group
  
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
  
  plt <- plot_ly(data = pca.ly, x = ~PC1, y = ~PC2,
                 text = ~paste0("Sample: ", sample), color = ~group,
                 type = "scatter", mode = "markers") %>%
    layout(xaxis = layout.x, yaxis = layout.y)
  
  return(plt)
}


deRes <- function(ns) {
  diffExpr.tab <- rbind(colSums(ns$deRes$q.value < 0.05 & ns$deRes$coefficients > 0),
                        colSums(ns$deRes$q.value < 0.05 & ns$deRes$coefficients < 0))
  diffExpr.tab <- as.data.frame(diffExpr.tab[,!(colnames(diffExpr.tab) %in% c("Intercept", "(Intercept)"))])
  
  if (ncol(diffExpr.tab) == 1) {
    colnames(diffExpr.tab) <- ""
    rownames(diffExpr.tab) <- c(paste0("Higher in ", colnames(ns$deRes$q.value)[!(colnames(ns$deRes$q.value) %in% c("Intercept", "(Intercept)"))]),
                                 paste0("Higher in ", ns$base.group))
  } else {
    rownames(diffExpr.tab) <- c("Higher in Group", 
                                paste0("Higher in ", ns$base.group))
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
                      options = list(columnDefs = list(list(className = 'dt-center', targets = "_all")))) %>%
    formatStyle(names(diffExpr.full)[substr(names(diffExpr.full), 1, 6) == "Log2FC"], 
                backgroundColor = styleInterval(brks.fc, cols.fc)) %>%
    formatStyle(names(diffExpr.full)[substr(names(diffExpr.full), 1, 5) == "q-val"], 
                backgroundColor = styleInterval(brks.qv, cols.qv))

  return(list(summary = diffExpr.tab,
              de = dtable))
  
}


plotlyHeatmap <- function(ns, groupedGenesets, leadingEdge, gsClust, gsComp, gsDir) {
  diffExpr <- makeDiffExprFile(ns$deRes, filename = NULL)
  dat.scaled <- diffExpr[,-(1:(4*(ncol(ns$deRes$t)-1)))]
  
  hm.max <- 3
  dat.scaled[dat.scaled > hm.max] <- hm.max
  dat.scaled[dat.scaled < -hm.max] <- -hm.max
  
  acols <- c("grey95", "black")
  names(acols) <- c(0, 1)
  
  ay1 <- ay <- ax1 <- ax <- list(title = "", ticks = "")
  ax$tickangle <- ax1$tickangle <- 20
  ax$tickfont <- list(size = 8)
  

    subClust <- rownames(groupedGenesets)[groupedGenesets$Cluster == gsClust]
    annot_row <- leadingEdge[[gsComp]][[as.numeric(gsDir)]]
    annot_row <- annot_row[,subClust]
    
    if (length(subClust) > 1) {
      annot_row <- as.data.frame(annot_row)[rowSums(annot_row) >= 1,]
      annot_row.df <- sapply(annot_row, as.numeric)
      #annot_row.df <- as.data.frame(annot_row.df)
      rownames(annot_row.df) <- rownames(annot_row)
      
      genesKeep <- rownames(dat.scaled) %in% rownames(annot_row.df)
    } else {
      genesKeep <- rownames(dat.scaled) %in% names(annot_row)[annot_row==1]
    }
    
    #Cluster genes
    genes <- labels(as.dendrogram(hclust(dist(dat.scaled[genesKeep,]), method="ward.D")))
    
    dat.hm <- dat.scaled[genes,]
    dat.hm <- reshape::melt(dat.hm)
    dat.hm$X1 <- factor(dat.hm$X1, levels = genes)
    dat.hm$X2 <- factor(dat.hm$X2, levels = colnames(dat.scaled))
    colnames(dat.hm)[3] <- "expression"
    
    if (length(subClust) > 1) {
      annot_row.hm <- reshape::melt(annot_row.df)
      annot_row.hm$X1 <- factor(annot_row.hm$X1, levels = genes)
      annot_row.hm$X2 <- factor(annot_row.hm$X2, levels = colnames(annot_row.df))
      ay1$showticklabels <- FALSE
    }
    
    #acols <- list()
    #for (i in colnames(annot_row.df)) {
    #  acols[[i]] <- c("grey95", "black")
    #  names(acols[[i]]) <- c("0", "1")
    #}
    rel.width <- 0.8*ncol(annot_row) / (ncol(annot_row)+ncol(dat.scaled))
    
    pHM <- plot_ly(data = dat.hm, x = ~X2, y = ~X1, z = ~expression, xgap = 2, ygap = 2, 
                   type = "heatmap", colors = colorRamp(c("blue", "white", "red")), alpha = 1,
                   zauto = FALSE, zmax = hm.max, zmin = -hm.max) %>%
      layout(xaxis = ax1, yaxis = ay1, margin = list(b = 100))
    
    if (length(subClust) > 1) {
      pCB <- plot_ly(data = annot_row.hm, x = ~X2, y = ~X1, z = ~value, xgap = 2, ygap = 2,
                     type = "heatmap", colors = colorRamp(c("grey95", "black")), alpha = 1, 
                     showscale = FALSE) %>%
        layout(xaxis = ax, yaxis = ay, margin = list(b = 100))
      return(subplot(pCB, pHM, nrows = 1, widths = c(rel.width, 1-rel.width)))
    } else {
      return(pHM)
    }
    
}