library(shinyBS)
library(Biobase)
library(limma)
library(plotly)
library(ggplot2)
library(DT)
library(NanoTube)
source("helpers.R")

shinyUI(navbarPage("Shiny-NanoTube", id="master",
                   tabPanel("Job Setup",
                            fluidPage(
                              
                              h2("Data Entry"),
                              
                              fileInput("expr",
                                        label = "NanoString data",
                                        multiple = FALSE),
                              bsTooltip("expr",
                                        "Either a folder containing .RCC files, or an expression matrix in a .csv or .txt file.",
                                        placement = "bottom", trigger = "hover", options = NULL),
                              
                              fluidRow(
                                column(4,
                                       fileInput("phen",
                                                 label = "Sample group classifiers",
                                                 multiple = FALSE),
                                       bsTooltip("phen",
                                                 "This file should contain a group classifier for each sample, in the same order as in the expression dataset. If not provided, groups will be interpreted from sample file names (see Help).",
                                                 placement = "bottom", trigger = "hover", options = NULL)
                                ),
                                
                                column(4,
                                       textInput("basePhen",
                                                 label = "Base group",
                                                 value = ""),
                                       bsTooltip("basePhen",
                                                 "The group against which other groups will be compared (such as the Control group). If not provided, the first group will be used.",
                                                 placement = "bottom", trigger = "hover", options = NULL)
                                )
                                
                                
                              ),
                              
                              fileInput("gsDb",
                                        label = "Gene set database (Optional)",
                                        multiple = FALSE),
                              bsTooltip("gsDb",
                                        "A gene set database file, either in .gmt format or an .rds file containing an R-format list of gene sets",
                                        placement = "bottom", trigger = "hover", options = NULL),
                              
                              actionButton("run",
                                           label = "Analyze Data"),
                              #submitButton("Analyze Data"),
                              
                              h3("Advanced Options"),
                              ### Add action button for conditional input
                              
                              h4("Normalization Options"),
                              
                              textInput("hk",
                                        label = "Housekeeping Genes",
                                        value = "ALAS1, ABCF1, TBP, PPIA, TUBB"),
                              bsTooltip("hk",
                                        "Optional: A list of housekeeping genes present in the input data, separated by commas. If not provided, housekeeping genes will be identified as marked in the input file."),
                              
                              numericInput("bgP",
                                           label = "Negative Control Threshold (t test p-value)",
                                           value = 0.05,
                                           min = 0.00001,
                                           max = 2),
                              bsTooltip("bgP",
                                        "Expression threshold (vs. negative control genes) for inclusion, in the form of a p-value from a 2-sample t test (see Help). To include all genes in analysis, set to 2."),
                              
                              h4("Gene Set Analysis Options"),
                              
                              numericInput("minSize",
                                           label = "Min size (exclude smaller sets)",
                                           value = 5,
                                           min = 0,
                                           max = 500),
                              bsTooltip("minSize",
                                        "Minimum size of gene sets for inclusion in pathway analysis (only genes included in NanoString data are counted)."),
                              
                              numericInput("setThresh",
                                           label = "t threshold",
                                           value = 2,
                                           min = 0,
                                           max = 20),
                              bsTooltip("setThresh",
                                        "Minimum t-statistic for mHG significance cutoff (see Help).")
                              
                            )),
                   
                   tabPanel("QC Results",
                            navbarPage("QC", id = "qc",
                              tabPanel("Positive Controls",
                                       fluidPage(
                                         fluidRow(
                                           column(8,
                                                  plotOutput("posPlot", width = "100%")
                                           ),
                                           
                                           column(4,
                                                  tableOutput("posTab")
                                           )))),
                              
                              tabPanel("Negative Controls",
                                       fluidPage(
                                         tableOutput("negTab"),
                                         
                                         plotlyOutput("negPlot")
                                       )),
                              
                              tabPanel("Housekeeping",
                                       fluidPage(
                                         fluidRow(
                                           column(8,
                                                  plotlyOutput("hkPlot1", width = "100%"),
                                                  plotlyOutput("hkPlot2", width = "100%")
                                           ),
                                           
                                           column(4,
                                                  tableOutput("hkTab")
                                           ))))
                            )
                            
                   ),
                   

                   tabPanel("Analysis Results",
                            navbarPage("DE", id = "de",
                                       tabPanel("PCA",
                                                fluidPage(
                                                  plotlyOutput("pcaPlot")
                                                )
                                       ),
                                       
                                       tabPanel("Differential Expression",
                                                fluidPage(
                                                  tableOutput("deCounts"),
                                                  DTOutput("deTab")
                                                )
                                       ),
                                       
                                       #conditionalPanel(
                                       #  condition = "input.gsDb !== null",
                                       tabPanel("Gene Set Analysis",
                                                fluidPage(

                                                    #inputPanel(
                                                      uiOutput("gsUI"),
                                                    #),
                                                    DTOutput("gsTab"),
                                                    
                                                    plotlyOutput("gsHM")
                                                  )
                                                  
                                                )
                                                
                                       #)
                            
                                       
                            )
                   ),
                   
                   tabPanel("Help!",
                            fluidPage(

                            )
                            
                   )
)
)

