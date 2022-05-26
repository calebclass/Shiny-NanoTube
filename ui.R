library(NanoTube)
library(shiny)
library(shinydashboard)
library(shinyBS)
library(shinyjs)
library(plotly)
library(DT)

source("helpers.R")
##Add gradient scroll 
##Add sticky header on scroll
#Add something interactive where person can see the files after they are uploaded
dashboardPage(skin = "blue",
              
              header,  # defined in 'helpers.R'
              
              dashboardSidebar(
                sidebarMenu(
                  menuItem("Welcome", tabName = "Welcome", icon = icon("flag")),
                  menuItem("Setup", tabName = "setup", icon = icon("th")),
                  menuItem("QC Results", tabName = "QCres", icon = icon("chart-line")),
                  menuItem("Differential Expression", tabName = "AnalysisRes", icon = icon("chart-bar")),
                  menuItem("Gene Set Analysis", tabName = "GSA", icon = icon("chart-bar")),
                  menuItem("Help", tabName = "Help", icon = icon("question-circle"))
                )
              ),
              
              dashboardBody(
                tabItems(
                  
                  tabItem(tabName = "Welcome",
                          fluidRow(
                            includeCSS("www/styling2.css"),
                            
                            box(
                              
                              width = 10,
                              
                              withTags({
                                div(class="norm", checked = NA,
                                    h1("Welcome to the NanoTube!"),
                                    img(src="flowchart.jpg", class = "resize", align = "left"),
                                    
                                    p("NanoTube performs data processing, quality control, normalization and analysis on NanoString gene expression data."),
                                    b("Click on the Setup tab to get started."),
                                    
                                    br(),
                                    br(),
                                    a(href='http://www.bioconductor.org/packages/release/bioc/html/NanoTube.html', b("Check out the NanoTube package on Bioconductor!")),
                                    p("This R package provides additional normalization and analysis options for NanoString nCounter data."),
                                    
                                    
                                    h2("Basic Features"),
                                    h3("Data Processing"),
                                    p("nCounter data are input as raw RCC files or CSV files (which possibly came from the RCC Collector tool). An additional sample information table is then loaded to allow comparisons between groups."),
                                    h3("Normalization"),
                                    p("This application performs manufacturer-recommended normalization steps, including positive and housekeeping normalization, as well as the removal of target genes found to have expression levels below 'background' (estimated from the negative control gene expression)."),
                                    h3("Analysis"),
                                    p("nCounter data are input as raw RCC files or CSV files (possibly from the RCC Collector tool), which can then be combined with a sample information table and saved as an ExpressionSet."),
                                    h3("Visualization"),
                                    p("nCounter data are input as raw RCC files or CSV files (possibly from the RCC Collector tool), which can then be combined with a sample information table and saved as an ExpressionSet."),
                                    
                                    h2("Citation"),
                                    p("If you use the NanoTube in your work, please cite our conference paper:"),
                                    b("Class CA, Bristow CA, Do K-A (2021). Easy NanoString Gene Expression Analysis with the Nanotube. The FASEB Journal 36(S1)."),
                                    
                                    h2("License"),
                                    p("The NanoTube is provided with the GNU General Public License (GPL-3), and without warranty."),
                                    a(href="https://github.com/calebclass/Shiny-NanoTube/blob/master/LICENSE", b("GPL-3 License for NanoTube")),
                                    br(),
                                    br(),
                                    br(),
                                    
                                    img(src="cophs_horiz_4cp_1-4w.jpg", class = "resize2", align = "left")
                                )
                              })
                            )
                          )
                          
                  ),
                  
                  tabItem(tabName = "setup",
                          
                          fluidRow(
                            column(width = 6,
                                   
                                   box(title = "Data Entry", width = NULL,
                                       
                                       fluidPage(
                                         
                                         fileInput("expr",
                                                   label = "NanoString data",
                                                   multiple = FALSE),
                                         bsTooltip("expr",
                                                   "Either a folder containing .RCC files, or an expression matrix in a .csv or .txt file.",
                                                   placement = "bottom", trigger = "hover", options = NULL),
                                         
                                         ######################
                                         
                                         fileInput("phen",
                                                   label = "Sample info table",
                                                   multiple = FALSE),
                                         bsTooltip("phen",
                                                   "This should be a CSV file, containing sample information.",
                                                   placement = "bottom", trigger = "hover", options = NULL),
                                         
                                         
                                         
                                         # Read columns in Sample info table, asks user which column corresponds to "Group"
                                         fluidRow(
                                           column(6,
                                                  uiOutput("phenCol_input"))
                                         ),
                                         
                                         # Of the groups in the "Group" column, which one is the group against which all others will be compared (the control group, for example)
                                         fluidRow(
                                           column(6,
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
                                         
                                         br(),
                                         verbatimTextOutput("numSamps"),
                                         br()
                                         
                                       )),
                                   
                                   box(title = "Advanced Options", width = NULL, 
                                       collapsible = TRUE, collapsed = TRUE,
                                       
                                       h4("Normalization Options", id = "gseaTxt"),
                                       
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
                                                 "Minimum size of gene sets for inclusion in pathway analysis (only genes included in NanoString data are counted).")
                                       
                                       
                                   )
                            ),
                            
                            column(width = 6,
                                   
                                   box(width = NULL,
                                       dataTableOutput("merged_info")
                                   )
                            ))
                  ),
                  
                  tabItem(tabName = "QCres",
                          navbarPage("QC", id = "qc",
                                     tabPanel("Positive Controls",
                                              box(
                                                column(width = 12, plotlyOutput("posPlot", width = "100%", height = "auto")),
                                                width = 8
                                              ),
                                              
                                              box(
                                                dataTableOutput("posTab"),
                                                width = 4
                                              )),
                                     
                                     tabPanel("Negative Controls",
                                              fluidRow(
                                                
                                                column(width = 4,
                                                       box(
                                                         tableOutput("negSummary"), width = NULL,
                                                         title = "Endogenous Targets vs. Negative Controls"
                                                       ),
                                                       
                                                       box(
                                                         dataTableOutput("negGenes"), width = NULL,
                                                         title = "Target Statistics"
                                                       )
                                                       
                                                ),
                                                
                                                column(width = 8,
                                                       box(
                                                         dataTableOutput("negTab"), width = NULL,
                                                         title = "Sample Statistics"
                                                       ),
                                                       
                                                       box(
                                                         plotlyOutput("negPlot", height = "auto"), width = NULL,
                                                         title = "Negative Target Counts"
                                                       )
                                                )
                                                
                                                
                                              )),
                                     
                                     tabPanel("Housekeeping",
                                              box(
                                                plotlyOutput("hkPlot1", width = "100%", height = "auto"),
                                                br(),br(),
                                                plotlyOutput("hkPlot2", width = "100%", height = "auto"),
                                                width = 8
                                              ),
                                              
                                              box(
                                                dataTableOutput("hkTab"),
                                                width = 4
                                              )))
                  ),
                  
                  tabItem(tabName = "AnalysisRes",
                          fluidRow(
                            column(width = 5,
                                   box(
                                     title = "PCA", width = NULL, 
                                     plotlyOutput("pcaPlot")
                                   ),
                                   
                                   box(
                                     title = "Volcano Plot", width = NULL,
                                     
                                     uiOutput("volUI"),
                                     numericInput('volcanoVertLineInput','log2(FC) cutoff', value = 0, min = 0, max = 10),
                                     numericInput('volcanoHorLineInput','p-val cutoff', value = 0.05, min = 0, max = 1),
                                     
                                     plotlyOutput("canoPlot")
                                   )
                            ),
                            column(width = 7,
                                   box(title = "Summary", width = NULL,
                                       column(width = 4, 
                                              numericInput('summaryQ', 'q-val cutoff', value = 0.05, min = 0, max = 1)),
                                       tableOutput("deCounts")),
                                   
                                   box(
                                     title = "Full Results", width = NULL,
                                     
                                     DTOutput("deTab"),
                                     downloadButton("DEdownload","Download Table")
                                   ))
                          )
                          
                          ####
                          #tabPanel("Nanostring data table",
                          #        fluidPage(
                          #          h4("datatable"),
                          #          tableOutput('nanoTable'),
                          #          downloadButton("NANOdownload","Download Table")
                          #          
                          #        ))
                          #####https://shiny.rstudio.com/reference/shiny/1.2.0/showTab.html
                          #conditionalPanel(
                          #  condition = "input.gsDb != null",
                          #)
                          
                          
                          
                  ),
                  
                  ####
                  tabItem(tabName = "GSA",
                          fluidRow(
                            column(width = 4,
                                   box(title = "GSEA Visualization Options", width = NULL,
                                       uiOutput("gsUI")
                                   )),
                            
                            column(width = 8,
                                   box(width = NULL,
                                       
                                       DTOutput("gsTab"),
                                   ),
                                   
                                   box(width = NULL,
                                       plotlyOutput("gsHM"),
                                   )
                            )
                          )
                          
                          
                          
                          
                  ),
                  #####
                  tabItem("Help",
                          fluidRow(
                            includeCSS("www/styling2.css"),
                            
                            box(title = "Setup", width = 10, 
                                collapsible = TRUE, collapsed = TRUE,
                            includeHTML("www/help_setup.html"),
                            ),
                            box(title = "Quality Control", width = 10, 
                                collapsible = TRUE, collapsed = TRUE,
                                includeHTML("www/help_qc.html"),
                            ),
                            box(title = "Analysis", width = 10, 
                                collapsible = TRUE, collapsed = TRUE,
                                includeHTML("www/help_de.html"),
                            ),
                            box(title = "Gene Set Analysis", width = 10, 
                                collapsible = TRUE, collapsed = TRUE,
                                includeHTML("www/help_gsa.html"),
                            )
                          )
                          
                  )
                  
                )
              )
)


