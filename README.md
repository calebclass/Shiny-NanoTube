# Shiny-NanoTube
<a href="https://zenodo.org/badge/latestdoi/334240182"><img src="https://zenodo.org/badge/334240182.svg" alt="DOI"></a>

## Installation
1. Download and install [RStudio (and R)](https://www.rstudio.com/products/rstudio/download/)
2. Install the [NanoTube package](http://www.bioconductor.org/packages/release/bioc/html/NanoTube.html) from Bioconductor:
```{r}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("NanoTube")
```
3. Install other required and recommended libraries.
```{r}
install.packages(c("shinyBS", "shinyjs", "plotly", "DT"))
BiocManager::install("qusage")
```
4. Open 'ui.R' in Shiny-NanoTube, and click 'Run'!

## A first analysis
1. Download example data sets from the 'data' folder.
2. Run the Shiny-NanoTube app (either from RStudio or the [web page](research.butler.edu/nanotube)).
3. In the **Setup** tab, load 'GSE132946.zip' in the **NanoString data** field, and load 'SampleData-GSE132946.csv' in the **Sample info table** field. A **Group Column** menu will appear, and open that menu to select 'outcome'. 
4. Click the **Check Samples** button to ensure that the expression data and sample metadata are merged correctly.  The 'Filename' and 'name' columns should match.
5. (Optional) Load a gmt file in the **Gene set database** field.
6. Click **Analyze Data**.
7. View results in the **QC Results**, **Differential Expression**, and **Gene Set Analysis** (optional) tabs.
8. If you're having trouble, check the **Help** page, raise an **Issue** on GitHub, or fill out [the form here](https://research.butler.edu/caleb-class-lab/nanotube/).

