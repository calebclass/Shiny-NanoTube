# Shiny-NanoTube

## Installation
1. Download and install [RStudio (and R)](https://www.rstudio.com/products/rstudio/download/)
2. Install the [NanoTube package](http://www.bioconductor.org/packages/release/bioc/html/NanoTube.html) from Bioconductor:
```{r}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("NanoTube")
```
3. Install other required libraries
```{r}
install.packages(c("shinyBS", "shinyjs", "plotly", "DT"))
```
4. Open 'ui.R' in Shiny-NanoTube, and click 'Run'!
