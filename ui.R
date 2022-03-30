library(NanoTube)
library(shinyBS)
library(shinyjs)
library(plotly)
library(DT)

source("helpers.R")
##Add gradient scroll 
##Add sticky header on scroll
#Add something interactive where person can see the files after they are uploaded
shinyUI(fluidPage(id="formatting",
                  
                  navbarPage(div(img(src="NanoTube-Logo.png", 
                                     style="float:right", height = 40, width = 110), ""),
                             tags$head(
                               tags$style(HTML('.navbar-nav > li > a, .navbar-brand {
                            padding-top:4px !important; 
                            padding-bottom:4px !important;
                            line-height: 40px !important;
                            height: 40px;
                            font-size: 15px;
                            }
                           .navbar {min-height:50px !important;}'))),
                             selected = "Welcome",
                             id = "master",
                             ####
                             #includeScript(),
                             useShinyjs(),
                             ##CSS here    
                             tags$head(
                               tags$link(rel = "stylesheet", type = "text/css", href = "styling.css")
                             ),
                             tags$style(type="text/css", "body {padding-top: 70px;}"),
                             
                             tabPanel("Welcome",
                                      fluidPage(
                                        h1("Welcome to the NanoTube!"),
                                        
                                        
                                        includeCSS("www/styling2.css"),
                                        includeHTML("www/landing.html"),
                                        includeScript("www/styling2.js")
                                        
                                       # img(src="www/cophs_horiz_4cp_1-4w.jpg", align = "left")
                                      )
                                      ),
                             
                             tabPanel("Job Setup",
                                      
                                      fluidPage(
                                        
                                        h2("Data Entry"),
                                        
                                        fluidRow(
                                          ######################
                                          column(12,
                                                 fileInput("expr",
                                                           label = "NanoString data",
                                                           multiple = FALSE),
                                                 bsTooltip("expr",
                                                           "Either a folder containing .RCC files, or an expression matrix in a .csv or .txt file.",
                                                           placement = "bottom", trigger = "hover", options = NULL)
                                          )
                                       
                                          
                                          
                                          ######################
                                        ),
                                        fluidRow(
                                          column(4,
                                                 fileInput("phen",
                                                           label = "Sample info table",
                                                           multiple = FALSE),
                                                 bsTooltip("phen",
                                                           "This should be a CSV file, containing sample information.",
                                                           placement = "bottom", trigger = "hover", options = NULL)
                                          ),
                                          column(12,
                                                 DTOutput("merged_info"))
                                        ),
                                        
                                        fluidRow(
                                          column(4,
                                                 uiOutput("phenCol_input"))
                                        ),
                                        
                                        fluidRow(
                                          column(4,
                                                 uiOutput("basePhen_input"))
                                          
                                          
                                        ),
                                        
                                        fileInput("gsDb",
                                                  label = "Gene set database (Optional)",
                                                  multiple = FALSE),
                                        
                                        bsTooltip("gsDb",
                                                  "A gene set database file, either in .gmt format or an .rds file containing an R-format list of gene sets",
                                                  placement = "bottom", trigger = "hover", options = NULL),
                                        
                                        actionButton("check",
                                                     label = "Check Samples"),
                                        
                                        actionButton("run",
                                                     label = "Analyze Data"),
                                        #submitButton("Analyze Data"),
                                        
                                        br(), br(),
                                        
                                        actionLink("adv", "Advanced Options"),
                                        
                                        h4("Normalization Options", id = "normTxt"),
                                        
                                        textInput("hk",
                                                  label = "Housekeeping Genes",
                                                  value = ""),
                                        bsTooltip("hk",
                                                  "Optional: A list of housekeeping genes present in the input data, separated by commas. If not provided, housekeeping genes will be identified as marked in the input file."),
                                        
                                        numericInput("bgP",
                                                     label = "Negative Control Threshold (t test p-value)",
                                                     value = 0.05,
                                                     min = 0.00001,
                                                     max = 2),
                                        bsTooltip("bgP",
                                                  "Expression threshold (vs. negative control genes) for inclusion, in the form of a p-value from a 2-sample t test (see Help). To include all genes in analysis, set to 2."),
                                        
                                        h4("Gene Set Analysis Options", id = "gseaTxt"),
                                        
                                        numericInput("minSize",
                                                     label = "Min size (exclude smaller sets)",
                                                     value = 5,
                                                     min = 0,
                                                     max = 500),
                                        bsTooltip("minSize",
                                                  "Minimum size of gene sets for inclusion in pathway analysis (only genes included in NanoString data are counted)."),
                                        
                                        
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
                                                                     br(),br(),
                                                                     plotlyOutput("hkPlot2", width = "100%")
                                                              ),
                                                              
                                                              column(4,
                                                                     tableOutput("hkTab")
                                                              ))))
                                      )
                                      
                             ),
                             
                             tabPanel("Analysis Results",
                                      navbarPage("Analysis", id = "de",
                                                 tabPanel("PCA",
                                                          fluidPage(
                                                            plotlyOutput("pcaPlot")
                                                          )
                                                 ),
                                                 
                                                 #####
                                                 tabPanel("Volcano",
                                                          fluidPage(
                                                            h4("Volcano Plot"),
                                                            
                                                            numericInput('volcanoVertLineInput','log2(FC) cutoff', value = 0, min = 0, max = 4),
                                                            numericInput('volcanoHorLineInput','p-val cutoff', value = 0.05, min = 0, max = 1),
                                                            
                                                            plotlyOutput("canoPlot")
                                                          )
                                                          
                                                 ),
                                                 #######
                                                 tabPanel("Differential Expression",
                                                          fluidPage(
                                                            h4("Differentially expressed genes (q < 0.05)"),
                                                            tableOutput("deCounts"),
                                                            
                                                            h4("Full Results"),
                                                            DTOutput("deTab"),
                                                            
                                                            downloadButton("DEdownload","Download Table")
                                                            #####
                                                            
                                                          )
                                                 ),
                                                 ####
                                                 tabPanel("Nanostring data table",
                                                          fluidPage(
                                                            h4("datatable"),
                                                            tableOutput('nanoTable'),
                                                            downloadButton("NANOdownload","Download Table"),
                                                            
                                                          )),
                                                 #####https://shiny.rstudio.com/reference/shiny/1.2.0/showTab.html
                                                 #conditionalPanel(
                                                 #  condition = "input.gsDb != null",
                                                 
                            
                                                 
                                                 
                                                 #)
                                                 
                                                 
                                      )
                             ),
                             
                             ####
                             tabPanel("Gene Set Analysis",
                                      navbarPage("GSA", id = "gsa",
                                                 tabPanel("Gene Set Analysis",
                                                          fluidPage(
                                                            
                                                            #inputPanel(
                                                            
                                                            uiOutput("gsUI"),
                                                            #),
                                                            DTOutput("gsTab"),
                                                            
                                                            plotlyOutput("gsHM"),
                                                            
                                                            includeHTML("www/gsa.html"),
                                                          )
                                                          
                                                 )
                                      )
                             ),
                             #####
                             tabPanel("Help",
                                      fluidPage(
                                        includeCSS("www/styling2.css"),
                                        includeHTML("www/help.html"),
                                        includeScript("www/styling2.js"),
                                      )
                                      
                                      ),
                             
                  )
)

)