#------------------------------------------------------------------------------
packages = c("shiny","rstudioapi","tidyverse","bslib","shinyWidgets","DT","shinycustomloader","shinyBS","fst","data.table","shinydashboard","shinydashboardPlus",
             "survival","survminer","circlize","broom","EGAnet","ComplexHeatmap","qgraph","plyr","DBI","RSQLite","MoffittFunctions","shinyjs","gbm","factoextra",
             "ggcorrplot","igraph","GSVA") # rstudioapi, not needed
# 
# packages = c("shiny","tidyverse","DT","data.table","heatmaply","ComplexHeatmap","qgraph","rstudioapi",
#              "fst","survival","survminer","shinycustomloader","circlize","shinyWidgets","EGAnet","broom",
#              "shinydashboard","shinydashboardPlus","bslib","shinyBS","igraph","MoffittFunctions","WGCNA","AdjKM.CIF","stringi")
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    # some packages must be installed in seperately
    if (x %in% c("ComplexHeatmap","DESeq2","edgeR","clusterProfiler","enrichplot","pathview","EBSeq")) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install(x)
      
    } 
    #else if(x=="shinyTree"){
    # devtools::install_github("shinyTree/shinyTree")
    #}
    else {
      install.packages(x, dependencies = TRUE,repos = "http://cran.us.r-project.org")
    }
  }
  library(x, character.only = TRUE)
})

############################## Import useful functions #################################

folder_main = "/CGPA_IOtherapy"

setwd(folder_main)
source("pan_cancer_dashboard.R")
source("functions.R")
source("KM_surv.R")
source("Adj_KM_surv.R")
source("GHI.R")


genesets_pca = fread("www/genesets_pca.csv")




dark = bs_theme(version = 3,bg="black",fg = "white",warning = "#FFD300" )%>%
  bs_add_rules(sass::sass_file("www/style/style.scss"))




ui <- fluidPage(
  theme = dark,
  useShinydashboard(),
  div(class="navbar1",
      navbarPage(title = "CGPA: Immnotherapy",
                 
                 fluid = TRUE, 
                 collapsible = TRUE,
                 tabPanel("Single-gene Discovery",
                          br(),
                          tabsetPanel(type="pills",id="mainnav",
                                      tabPanel("PAN-CANCER SUMMARY",
                                               br(),
                                               fluidRow(
                                                 splitLayout(
                                                   cellWidths = c("80%","17%"), 
                                                   br(),
                                                   textInput("input_gene", NULL,value = "EGFR") 
                                                 )
                                                 
                                               ),
                                               fluidRow(
                                                 column(7,
                                                        div(id="box1",
                                                            box(
                                                              title="Prognostic Marker Summary",collapsible = TRUE, 
                                                              status="warning",
                                                              solidHeader = TRUE,
                                                              width = 12,height=600,
                                                              (div(style='max-width: 100%; height: 500px; width: auto;overflow-y: scroll;',
                                                                   div(style="text-align: center;",
                                                                       # DTOutput("test")
                                                                       plotOutput("forest")
                                                                   ),
                                                                   br(),
                                                                   div(style="text-align: left;",
                                                                       textOutput("summary_forest1"),
                                                                       textOutput("summary_forest2"))
                                                                   
                                                              )) 
                                                            )
                                                        )
                                                        
                                                 ),
                                                 column(5, pan_ui("box4","Gene Expression Profile",600,"cir_bar")),
                                                 column(7,
                                                        DTOutput("test"),
                                                        pan_ui2("box4","Kaplan-Meier Plots",530,"KM_general")),
                                                 column(5, pan_ui("box4","PPI network (STRING)",600,"string_net"))
                                               )
                                               
                                      ),
                                      tabPanel("MULTIVARIABLE ANALYSIS",
                                               br(),
                                               fluidRow(
                                                 column(3,
                                                        selectizeInput("study_type","Study type",
                                                                       choices=    c("GEO159067 (HNSC)","GSE135222 (NSCLC)","GSE179730 (COSCC)","GSE221733 (NSCLC)","GSE241876 (TNBC)",
                                                                                     "IMvigor (BLCA)","Javelin101 (KIRC)","Kallisto (BLCA)","Harmonized (Melanoma)"),multiple=F,selected="IMvigor (BLCA)"),
                                                        
                                                        sel_input_ui("input_gene_dash","Select a gene"),
                                                        sel_input_ui("survtype","OS or PFS/PFI",choices=NULL,multiple=FALSE),
                                                        sel_input_ui("km_type","KM or Adjusted KM",c("KM","Adjusted KM")),
                                                        sel_input_ui("cutoff_multi","KM cutoff",c("median","optimal","quartile")),
                                                        conditionalPanel(
                                                          condition="input.km_type=='Adjusted KM'",
                                                          style = "display: flex; align-items: flex-start;",  # Align items in the center
                                                          sel_input_ui("adj_cov_multi","Adjusted covariates",choices=NULL,multiple=TRUE)
                                                        )
                                                        
                                                 ),
                                                 column(9,
                                                        fluidRow(
                                                          column(9,br(),fluidRow(column(12,uiOutput("survival_dash")),
                                                                                 DTOutput("testt")
                                                          ))
                                                        )
                                                 )
                                               )
                                               
                                               
                                      ),
                                      tabPanel("GENE-HALLMARK INTERACTION",
                                               br(),
                                               fluidRow(
                                                 column(3,
                                                        selectizeInput("study_type_ghi","Study type",
                                                                       choices=    c("GEO159067 (HNSC)","GSE135222 (NSCLC)","GSE179730 (COSCC)","GSE221733 (NSCLC)","GSE241876 (TNBC)",
                                                                                     "IMvigor (BLCA)","Javelin101 (KIRC)","Kallisto (BLCA)","Harmonized (Melanoma)"),multiple=F,selected="IMvigor (BLCA)"),
                                                        sel_input_ui("gene_single_ghi","Choose a gene",choices=NULL,multiple=F),
                                                        
                                                        selectizeInput("survtype_ghi","OS or PFS/PFI",choices=NULL,multiple=F),
                                                        sel_input_ui("hallmark_match","Choose hallmarks",choices = NULL,multiple=T),
                                                        p("Click 'Run' button to get/update the figure"),
                                                        actionButton("search_single_ghi","Run")
                                                 ),
                                                 column(9,
                                                        div(id="box4",
                                                            box(
                                                              title="Gene-hallmark interaction model (z-score)",collapsible = TRUE, status="warning",
                                                              solidHeader = TRUE,
                                                              width = 8,height=800,
                                                              dropdown(
                                                                p("GHI full model: Hazard ~ a x Gene + b x Hallmark + c x Gene X Hallmark"),
                                                                p("GHI partial model: Hazard ~ a x Gene + c x Gene x Hallmark"),
                                                                br(),
                                                                p("z score of the interaction term between the gene and hallmark"),
                                                                DTOutput("GHI_zscore"),
                                                                circle = TRUE, status = "warning", icon = icon("cogs"), width = "600px",
                                                                tooltip = tooltipOptions(title = "Click to see more")
                                                              ),
                                                              tags$style(HTML('#sw-content-dropdown, .sw-dropdown-in {background-color: gray;}')),
                                                              #DTOutput("test_match")
                                                              
                                                              p("Interaction model"),
                                                              div(withLoader(plotOutput("edge_plot",width="80%",height="400px"),type="html",loader="loader1"),align="center")
                                                              
                                                            )
                                                            
                                                        )
                                                 )
                                               )
                                      )
                          )
                 ),
                 tabPanel("Gene Pair Discovery",
                          br(),
                          fluidRow(
                            column(2,
                                   br(),br(),
                                   selectizeInput("study_type_dual","Study type",
                                                  choices=   c("GEO159067 (HNSC)","GSE135222 (NSCLC)","GSE179730 (COSCC)","GSE221733 (NSCLC)","GSE241876 (TNBC)",
                                                               "IMvigor (BLCA)","Javelin101 (KIRC)","Kallisto (BLCA)","Harmonized (Melanoma)"),multiple=F,selected="IMvigor (BLCA)"),
                                   
                                   sel_input_ui("survtype_dual","OS or PFS/PFI",choices=NULL),
                                   selectizeInput("gene_dual","Select two genes",choices=NULL,multiple=T,options = list(maxItems = 2)),
                                   sel_input_ui("KM_cutoff","KM cutoff",choices=c("quartile","optimal","median"),multiple=F,selected = "median")
                                   # actionButton("run_dual", "Run")
                            ),
                            column(10,
                                   tabsetPanel(type="pills",id="mainnav", 
                                               br(),
                                               tabPanel("UNIVARIABLE MODEL",
                                                        column(5,
                                                               div(id="box1",
                                                                   box(
                                                                     title="KM plot",collapsible = TRUE, status="warning",
                                                                     solidHeader = TRUE,
                                                                     width = 12,height=800,
                                                                     withLoader(plotOutput("uni_km",height="280px"),type="html",loader="loader1"),
                                                                     br(),
                                                                     p("Log rank test for the selected genes"),
                                                                     div(DTOutput("uni_km_table"), style = "font-size:80%")  
                                                                   )
                                                               )
                                                        ),
                                                        column(5,
                                                               div(id="box1",
                                                                   box(
                                                                     title="Univariable Cox model",collapsible = TRUE, status="warning",
                                                                     solidHeader = TRUE,
                                                                     width = 12,height=800,
                                                                     plotOutput("blankUni1",height="50px"),
                                                                     plotOutput("uni_cox_plot1",width="85%",height="150px"),
                                                                     plotOutput("blankUni2",height="60px"),
                                                                     plotOutput("uni_cox_plot2",width="85%",height="150px"),
                                                                     plotOutput("blankUni3",height="85px")
                                                                   )
                                                               )
                                                        )
                                                        
                                                        
                                                        
                                               ),
                                               
                                               
                                               tabPanel("TWO GENE MODEL",
                                                        
                                                        column(5,
                                                               div(id="box1",
                                                                   box(
                                                                     title="KM plot",collapsible = TRUE, status="warning",
                                                                     solidHeader = TRUE,
                                                                     width = 12,height=800,
                                                                     br(),
                                                                     withLoader(plotOutput("multi_km",height="350px",width="60%"),type="html",loader="loader1"),
                                                                     br(),
                                                                     p("Log rank test for the selected genes"),
                                                                     div(DTOutput("multi_km_table"), style = "font-size:90%") 
                                                                     
                                                                   )
                                                               )
                                                        ),
                                                        column(5,
                                                               div(id="box1",
                                                                   box(
                                                                     title="Pairwise Log-rank test",collapsible = TRUE, status="warning",
                                                                     solidHeader = TRUE,
                                                                     width = 12,height=400,
                                                                     br(),
                                                                     plotOutput("blankpairlog1",height="15px"),
                                                                     div(DTOutput("pairwise_logrank"), style = "font-size:100%"),
                                                                     plotOutput("blankpairlog2",height="22px")
                                                                   )
                                                               )
                                                        )
                                                        
                                                        ,
                                                        column(5,
                                                               div(id="box1",
                                                                   box(
                                                                     title="Multivariable Cox model",collapsible = TRUE, status="warning",
                                                                     solidHeader = TRUE,
                                                                     width = 12,height=400,
                                                                     plotOutput("multi_cox",height="200px")
                                                                   )
                                                               )
                                                        )
                                                        
                                               ),
                                               tabPanel("INTERACTION MODEL",
                                                        
                                                        column(5,
                                                               div(id="box1",
                                                                   box(
                                                                     title="KM plot",collapsible = TRUE, status="warning",
                                                                     solidHeader = TRUE,
                                                                     width = 12,height=800,
                                                                     br(),
                                                                     plotOutput("blankinteractionkm1",height="50px"),
                                                                     plotOutput("km_interaction",height="340px"),
                                                                     plotOutput("blankinteractionkm2",height="50px")
                                                                     
                                                                   )
                                                               )
                                                        ),
                                                        
                                                        column(5,
                                                               div(id="box1",
                                                                   box(
                                                                     title="Multivariable Cox model with interaction",collapsible = TRUE, status="warning",
                                                                     solidHeader = TRUE,
                                                                     width = 12,height=800,
                                                                     br(),
                                                                     plotOutput("blankinteraction1",height="120px"),
                                                                     DTOutput("cox_interaction"),
                                                                     plotOutput("blankinteraction2",height="180px")
                                                                     
                                                                   )
                                                               )
                                                        )
                                               ),
                                               tabPanel("GENE RATIO MODEL",
                                                        column(5,
                                                               div(id="box1",
                                                                   box(
                                                                     title="KM plot",collapsible = TRUE, status="warning",
                                                                     solidHeader = TRUE,
                                                                     width = 12,height=800,
                                                                     plotOutput("km_ratio",height="350px",width="80%"),
                                                                     br(),
                                                                     p("Log-rank test for the gene ratio"),
                                                                     div(DTOutput("ratio_km_table"), style = "font-size:80%")  
                                                                     
                                                                   )
                                                               )
                                                        ),
                                                        
                                                        column(5,
                                                               div(id="box1",
                                                                   box(
                                                                     title="Cox model with gene ratio",collapsible = TRUE, status="warning",
                                                                     solidHeader = TRUE,
                                                                     width = 12,height=800,
                                                                     br(),
                                                                     plotOutput("blankRatio1",height="130px"),
                                                                     plotOutput("cox_ratio",height="200px"),
                                                                     plotOutput("blankRatio2",height="150px")
                                                                     
                                                                   )
                                                               )
                                                        )
                                               ),
                                               tabPanel("CORRELATION",
                                                        column(5,
                                                               div(id="box1",
                                                                   box(
                                                                     title="Correlation",collapsible = TRUE, status="warning",
                                                                     solidHeader = TRUE,
                                                                     width = 12,height=800,
                                                                     plotOutput("Scatterplot_dual",height="450px")
                                                                     
                                                                   )
                                                               )
                                                        )
                                               )
                                   )  
                            )
                            
                          )
                          
                 ),
                 tabPanel("Multi-gene Panel",
                         
                          sidebarLayout(
                                        sidebarPanel(width=2,
                                               br(),
                                               selectizeInput("study_type_multi","Study type",choices=c("GEO159067 (HNSC)","GSE135222 (NSCLC)","GSE179730 (COSCC)","GSE221733 (NSCLC)","GSE241876 (TNBC)",
                                                                                                        "IMvigor (BLCA)","Javelin101 (KIRC)","Kallisto (BLCA)","Harmonized (Melanoma)"),
                                                              multiple=F,selected="IMvigor (BLCA)"),
                                               selectizeInput("pathway","Select MSigDB pathway",
                                                              choices = c("None","MsigDB HALLMARK","MSigDB C2 (curated)","MSigDB C5 (ontology)","MSigDB C6 (oncogenic)","MSigDB C7 (immunologic)","Others"),selected = "MsigDB HALLMARK",multiple=FALSE),
                                               selectizeInput("geneset",label="Select a geneset",choices = NULL,multiple = FALSE),
                                               selectizeInput("gene_multi", label="Or customize gene list", 
                                                              choices =NULL,multiple=TRUE),
                                               selectizeInput("surv_multi","OS or PFI",choices = NULL),
                                               # wellPanel(
                                               #   useShinyjs(),
                                               #   div(
                                               #     fileInput('file', 
                                               #               span(
                                               #                 list(HTML("<p><abbr title='The uploaded data has to be one column csv file with geneNames on each row'>Or, upload your genelist (csv file)...</abbr></p>"))
                                               #               ),
                                               #               accept = c(
                                               #                 '.csv'
                                               #               )),
                                               #     style="font-size:80%;"
                                               #   ),
                                               #   actionButton('reset', 'Reset the file'),#Clear input dataset#
                                               #   tags$head(
                                               #     tags$style(HTML('#reset{padding:8px;font-size:80%}'))
                                               #   )
                                               #   
                                               # )
                                             
                                               ),
                                        mainPanel(width = 10,
                                               br() ,
                                            tabsetPanel(type="pills",id="mainnav",
                                               tabPanel("PROGNOSTIC RANKING",
                                                                 br(),
                                                                 
                                                                 fluidRow(column(9,
                                                                                 br(),
                                                                                 fluidRow(   
                                                                                   column(6,div(id="box2",
                                                                                                box(
                                                                                                  title="Waterfall plot for univariate Cox",collapsible = TRUE, status="warning",
                                                                                                  solidHeader = TRUE,
                                                                                                  dropdown(
                                                                                                    radioButtons("waterfall_menu","Exclude non-significant markers?",c("Yes","No"),"No"),
                                                                                                    # actionButton("search_net", "Run"),
                                                                                                    circle = TRUE, status = "warning", icon = icon("gear"), width = "300px",
                                                                                                    tooltip = tooltipOptions(title = "Click to see inputs !")
                                                                                                  ),
                                                                                                  width = 12,height=500,
                                                                                                  div(withLoader(plotOutput("waterfall_plot",width="80%"),type="html",loader="loader1"),align="center"),
                                                                                                  DTOutput("test_multi")
                                                                                                  
                                                                                                  
                                                                                                )
                                                                                   )),
                                                                                   column(6,div(id="box2",
                                                                                                box(
                                                                                                  title="Gradient boosting importance score",collapsible = TRUE, status="warning",
                                                                                                  solidHeader = TRUE,
                                                                                                  dropdown(
                                                                                                    numericInput("gb_menu","Input a threshold to exclude low influence markers",0),
                                                                                                    # actionButton("search_net", "Run"),
                                                                                                    circle = TRUE, status = "warning", icon = icon("gear"), width = "300px",
                                                                                                    tooltip = tooltipOptions(title = "Click to see inputs !")
                                                                                                  ),
                                                                                                  tags$style(HTML('#sw-content-dropdown, .sw-dropdown-in {background-color: gray;}')),
                                                                                                  width = 12,height=500,
                                                                                                  div(withLoader(plotOutput("gbm_plot",width="80%"),type="html",loader="loader1"),align="center")
                                                                                                
                                                                                                  
                                                                                                )
                                                                                   ))
                                                                                   
                                                                                 )
                                                                 )
                                                                 
                                                                 )
                                                        
                                                        ),
                                               tabPanel("JOINT SIGNATURE",
                                                        br(),
                                                        
                                                        fluidRow(
                                                          column(12,
                                                                 tabsetPanel(type="pills",id="mainnav_surv",
                                                                             tabPanel("Average score",
                                                                                      br(),
                                                                                      uiOutput("average_km_ui"),
                                                                                      
                                                                                      fluidRow(
                                                                                        column(9,br(),div(id="box2",
                                                                                                          box(
                                                                                                            title="Average score KM plot",collapsible = TRUE, status="warning",
                                                                                                            solidHeader = TRUE,
                                                                                                            width = 12,height=500,
                                                                                                            dropdown(
                                                                                                              downloadButton("download_average_KM", "Download"),
                                                                                                              # actionButton("search_net", "Run"),
                                                                                                              circle = TRUE, status = "warning", icon = icon("gear"), width = "300px",
                                                                                                              tooltip = tooltipOptions(title = "Click to download")
                                                                                                            ),
                                                                                                            conditionalPanel(
                                                                                                              condition="input.search_average_joint>0",
                                                                                                              div(withLoader(plotOutput("average_km_plot",width="60%"),type="html",loader="loader1"),align="center")
                                                                                                              #DTOutput("test_path")
                                                                                                            ),
                                                                                                            conditionalPanel(
                                                                                                              condition="input.search_average_joint==0",
                                                                                                              br(),
                                                                                                              p("Click 'Run' button to get the figure")
                                                                                                            )
                                                                                                            
                                                                                                            
                                                                                                          )
                                                                                        )
                                                                                        ),
                                                                                        column(9,br(),div(id="box2",
                                                                                                          box(
                                                                                                            title="Average score table",collapsible = TRUE, status="warning",
                                                                                                            solidHeader = TRUE,
                                                                                                            width = 12,height=500,
                                                                                                            dropdown(
                                                                                                              downloadButton("download_average_table", "Download"),
                                                                                                              # actionButton("search_net", "Run"),
                                                                                                              circle = TRUE, status = "warning", icon = icon("gear"), width = "300px",
                                                                                                              tooltip = tooltipOptions(title = "Click to download")
                                                                                                            ),
                                                                                                            conditionalPanel(
                                                                                                              condition="input.search_average_joint>0",
                                                                                                              div(withLoader(DTOutput("average_km_table"),type="html",loader="loader1"),align="center")
                                                                                                              #DTOutput("test_path")
                                                                                                            ),
                                                                                                            conditionalPanel(
                                                                                                              condition="input.search_average_joint==0",
                                                                                                              br(),
                                                                                                              p("Click 'Run' button to get the figure")
                                                                                                            )
                                                                                                            
                                                                                                            
                                                                                                          )
                                                                                        )
                                                                                        )
                                                                                        
                                                                                      )
                                                                                      
                                                                             ),
                                                                             tabPanel("ssGSEA",
                                                                                      br(),
                                                                                      uiOutput("gsea_km_ui"),
                                                                                      
                                                                                      fluidRow(
                                                                                        column(
                                                                                          9,br(),div(id="box2",
                                                                                                     box(
                                                                                                       title="ssGSEA KM plot",collapsible = TRUE, status="warning",
                                                                                                       solidHeader = TRUE,
                                                                                                       width = 12,height=500,
                                                                                                       dropdown(
                                                                                                         downloadButton("download_gsea_KM", "Download"),
                                                                                                         # actionButton("search_net", "Run"),
                                                                                                         circle = TRUE, status = "warning", icon = icon("gear"), width = "300px",
                                                                                                         tooltip = tooltipOptions(title = "Click to download")
                                                                                                       ),
                                                                                                       conditionalPanel(
                                                                                                         condition="input.search_gsea>0",
                                                                                                         div(withLoader(plotOutput("gsea_km_plot",width="60%"),type="html",loader="loader1"),align="center")
                                                                                                         #DTOutput("test_path")
                                                                                                         
                                                                                                       ),
                                                                                                       conditionalPanel(
                                                                                                         condition="input.search_gsea==0",
                                                                                                         br(),
                                                                                                         p("Click 'Run' button to get the figure")
                                                                                                       )
                                                                                                       
                                                                                                       
                                                                                                     )
                                                                                          )),
                                                                                        column(
                                                                                          9,br(),div(id="box2",
                                                                                                     box(
                                                                                                       title="ssGSEA Table",collapsible = TRUE, status="warning",
                                                                                                       solidHeader = TRUE,
                                                                                                       width = 12,height=500,
                                                                                                       dropdown(
                                                                                                         downloadButton("download_gsea_table", "Download"),
                                                                                                         # actionButton("search_net", "Run"),
                                                                                                         circle = TRUE, status = "warning", icon = icon("gear"), width = "300px",
                                                                                                         tooltip = tooltipOptions(title = "Click to download")
                                                                                                       ),
                                                                                                       conditionalPanel(
                                                                                                         condition="input.search_gsea>0",
                                                                                                         div(withLoader(DTOutput("gsea_table"),type="html",loader="loader1"),align="center")
                                                                                                         #DTOutput("test_path")
                                                                                                         
                                                                                                       ),
                                                                                                       conditionalPanel(
                                                                                                         condition="input.search_gsea==0",
                                                                                                         br(),
                                                                                                         p("Click 'Run' button to get the figure")
                                                                                                       )
                                                                                                       
                                                                                                       
                                                                                                     )
                                                                                          ))
                                                                                        
                                                                                        
                                                                                      )
                                                                                      
                                                                             )
                                                                 ),
                                                                 
                                                          )
                                                          
                                                        )
                                               ),
                                               tabPanel("SUBNETWORK",
                                                        br(),
                                                        fluidRow(
                                                          # column(2,
                                                          #        div(style = "display: flex; align-items:  flex-start;",  # Ensure items are vertically centered
                                                          #            tags$div("Include anchor genes?", style = "margin-right: 5px; font-size: 16px; white-space: nowrap; line-height: 1.5;margin-top:5px;"),  # Keep label styling simple
                                                          #            selectizeInput("anchor_not", NULL, c("No", "Yes"), "No", options = list(width = '100%'))
                                                          #        )
                                                          #        # radioButtons("anchor_not","Include anchor genes?",choices=c("No","Yes"),"No"),
                                                          #        
                                                          # )
                                                          column(5, uiOutput("subnet_ui"))
                                                        ),
                                                        fluidRow(
                                                          column(9,
                                                                 fluidRow(
                                                                   column(12,
                                                                          br(),
                                                                          div(id="box2",
                                                                              box(
                                                                                title="Subnetwork",collapsible = TRUE, status="warning",
                                                                                solidHeader = TRUE,
                                                                                width = 12,height="130%",
                                                                                conditionalPanel(
                                                                                  condition = "input.search_subnet>0",
                                                                                  div(style='overflow-x: scroll;overflow-y: scroll;',
                                                                                      withLoader(plotOutput("subnetwork",width="80%",height="600px"),type="html",loader="loader1"),align="center")
                                                                                  #DTOutput("test"))
                                                                                ),
                                                                                conditionalPanel(
                                                                                  condition="input.search_subnet==0",
                                                                                  br(),
                                                                                  p("Click 'Run' button to get the figure")
                                                                                )
                                                                              ))
                                                                   ),
                                                                   column(12,
                                                                          br(),
                                                                          div(id="box2",
                                                                              box(
                                                                                title="forest plot for each module",collapsible = TRUE, status="warning",
                                                                                solidHeader = TRUE,
                                                                                width = 12,height="130%",
                                                                                conditionalPanel(
                                                                                  condition = "input.search_subnet>0",
                                                                                  div(style='overflow-x: scroll;overflow-y: scroll;',
                                                                                      withLoader(plotOutput("subnetwork_forest",width="30%",height="400px"),type="html",loader="loader1"),align="center")
                                                                                  
                                                                                ),
                                                                                conditionalPanel(
                                                                                  condition="input.search_subnet==0",
                                                                                  br(),
                                                                                  p("Click 'Run' button to get the figure")
                                                                                )
                                                                              ))
                                                                   )
                                                                 )
                                                                 
                                                          )
                                                          
                                                        )),
                                               tabPanel("MULTI-GENE CORRELATION",
                                                        br(),
                                                        fluidRow(
                                                          #column(1),
                                                          column(8,div(id="box1",
                                                                       box(
                                                                         title="Correlation Network (Click the 'gear' button to select threshold)",collapsible = TRUE,status="warning",
                                                                         solidHeader = TRUE,
                                                                         width = 12,height=800,
                                                                         dropdown(
                                                                           numericInput("p_network","Choose a minimum correlation threshold (click 'update' once complete)",value=0.3,min=0,max=0.9),
                                                                           # actionButton("search_net", "Update"),
                                                                           circle = TRUE, status = "warning", icon = icon("gear"), width = "300px",
                                                                           tooltip = tooltipOptions(title = "Click to see inputs !")
                                                                         ),
                                                                         tags$style(HTML('#sw-content-dropdown, .sw-dropdown-in {background-color: gray;}')),
                                                                         div(withLoader(plotOutput("corr_network_plot",width="80%",height="500px"),type="html",loader="loader1"),align="center")
                                                                       )
                                                          )),
                                                          column(3),
                                                          column(1),
                                                          column(8,div(id="box1",
                                                                       box(
                                                                         title="PCA plot",collapsible = TRUE, status="warning",
                                                                         solidHeader = TRUE,
                                                                         width = 12,height=800,
                                                                         div(withLoader(plotOutput("pca_plot",width="80%",height="500px"),type="html",loader="loader1"),align="center")
                                                                         
                                                                       )
                                                                       
                                                                       
                                                          )
                                                          ),
                                                          column(3),
                                                          column(1),
                                                          column(8,div(id="box1",
                                                                       box(
                                                                         title="Correlogram (Spearman correlation)",collapsible = TRUE,status="warning",
                                                                         solidHeader = TRUE,
                                                                         width = 12,height=800,
                                                                         div(withLoader(plotOutput("corr_plot",width="80%",height="500px"),type="html",loader="loader1"),align="center")
                                                                       )
                                                          )),
                                                          column(3)
                                                          
                                                        )
                                               )
                                               )
                                      )
                                    
                          )
                 ),
                 tabPanel("Data Download",
                          DTOutput("data_ICI_info")
                          )
      ))
  
)

server <- function(input, output, session) {
  observe({
    query = parseQueryString(session$clientData$url_search)
    cat('input_gene ', query$geneid, '\n')
    
    if (!is.null(query[['geneid']])) {
      updateTextInput(session, "input_gene", "",value = toupper(query[['geneid']]))
    }
    
  })
  
  geneName = reactive({
    read_fst("www/data/geneNames.fst")
  })
  
  pheno = reactive({
    dt = read_fst(paste0(folder_main,"www/data/dt_ICI.fst"),columns=c("age","gender","clinical_stage","alcohol","smoking","hpv.status","immunotherapy.line",
                                                                      "anatomical.location","study","OS.time","OS","PFS.time","PFS","FMOne_mutation_burden_per_MB","Neoantigen_burden_per_MB",
                                                                      "Hypoxia_ssgsea","ISG.RS_ssgsea","IFN.Hallmark_ssgsea","CYT","CTL"
    ))
    dt
  }) ####################### CHANGE!!!!
  
  
  #=================================================================================
  # pan-cancer dashboard 
  #=================================================================================
  #--------------------------------------------------------------------------------
  # Forest plot #
  #------------------------------------------------------------------------------
  dt_forest = reactive({
    dt = read_fst(paste0(folder_main,"www/data/dt_ICI.fst"),columns=c("ID","OS.time","OS",
                                                                      "PFS.time","PFS","study",input$input_gene))
    names(dt) = make.names(names(dt))
    dt
  })
  
  os_forest_res = reactive({
    forest_table_server(input,output,session,input_dt = dt_forest(),surv_type="OS",sel_gene = make.names(input$input_gene))
  })
  pfs_forest_res = reactive({
    forest_table_server(input,output,session,input_dt = dt_forest(),surv_type="PFS",sel_gene = make.names(input$input_gene))
  })
  
  os_forest_plt = reactive({
    forest_fig(input,output,session,res=os_forest_res(),surv_type="OS")
    
  })
  
  pfs_forest_plt = reactive({
    forest_fig(input,output,session,res=pfs_forest_res(),surv_type="PFS")
    
  })
  
  # output$test = renderDT({
  #   pfs_forest_res()
  # })
  output$forest = renderPlot({
    tryCatch({
      p1 = os_forest_plt()
      p2 = pfs_forest_plt()
      ggarrange(p1,p2,ncol=2)
    }
    ,error=function(e){
      ""
    }
    )
    
  })
  
  output$summary_forest1 = renderText({
    tryCatch({
      pan_server_text(os_forest_res(),type="OS")
      
    },error = function(e){
      ""
    })
  }) 
  
  output$summary_forest2 = renderText({
    tryCatch({
      pan_server_text(pfs_forest_res(),type="PFS")
      
    },error = function(e){
      ""
    })
  
  })
  
  #--------------------------------------------------------------------------------
  # KM plot for significant genes from univariable cox, PUDATE: SHOW ALL KM not matter sig or not
  #------------------------------------------------------------------------------
  output$KM_general = renderUI({
    div(style='max-width: 100%; height: 500px; width: auto;overflow-y: scroll;',
        withLoader(plotOutput("km_optimal",width="70%",height="800px"),type="html",loader="loader1"),align="center")
  })
  

  
  km_plt_optimal = reactive({
  #  tryCatch({
    dt_sub  = dt_forest()
    covariates = os_forest_res()$study[!is.na(os_forest_res()$pval)]
    dt_sub_km = subset(dt_sub,dt_sub$study%in%covariates)
    
      groups <- dt_sub_km %>% group_split(study)
      dt_sub_km$OS.time = as.numeric(dt_sub_km$OS.time)
      dt_sub_km$OS = as.numeric(dt_sub_km$OS)
      
      # Generate the plots
      plots <- lapply(1:length(groups), function(x) {
        generate_km_plot(groups[[x]], make.names(input$input_gene), unique(groups[[x]]$study),cutoff=input$KM_cutoff_pan)
      })
      
      num_plots <- length(plots)
      num_cols <- ceiling(sqrt(num_plots))
      num_rows <- ceiling(num_plots / num_cols)
      
      # Arrange the plots
      ggarrange(plotlist = plots, ncol = num_cols, nrow = num_rows)
    # }
    # ,error=function(e){
    #   ""            
    # }
    # )
    
  })
  
  # km_plt_optimal = reactive({
  #   tryCatch({
  #     dt_sub  = dt_forest()
  #     covariates = os_forest_res()$study[os_forest_res()$pval<0.05]
  #     
  #     dt_sub_km = subset(dt_sub,dt_sub$study%in%covariates)
  #     
  #     groups <- dt_sub_km %>% group_split(study)
  #     dt_sub_km$OS.time = as.numeric(dt_sub_km$OS.time)
  #     dt_sub_km$OS = as.numeric(dt_sub_km$OS)
  #     
  #     # Generate the plots
  #     plots <- lapply(1:length(groups), function(x) {
  #       generate_km_plot(groups[[x]], make.names(input$input_gene), unique(groups[[x]]$study))
  #     })
  #     
  #     num_plots <- length(plots)
  #     num_cols <- ceiling(sqrt(num_plots))
  #     num_rows <- ceiling(num_plots / num_cols)
  #     
  #     # Arrange the plots
  #     ggarrange(plotlist = plots, ncol = num_cols, nrow = num_rows)
  #   }
  #   ,error=function(e){
  #     ggplot() +
  #       theme_void() +
  #       theme(
  #         panel.background = element_rect(fill = "white", color = NA),
  #         plot.background = element_rect(fill = "white", color = NA)
  #       )
  #     
  #   }
  #   )
  #   
  #   
  # })
  
  # output$test = renderDT({
  #   km_plt_optimal()
  # })
  output$km_optimal = renderPlot({
    km_plt_optimal()
  })
  
  #-------------------------------------------------------------------------------
  # Circular plot
  #-------------------------------------------------------------------------------
  output$cir_bar = renderUI({
    div(style='max-width: 100%; height: 500px; width: auto;overflow-y: scroll;',
        withLoader(plotOutput("circular_bar_plt",width="90%",height="500px"),type="html",loader="loader1"),align="center")
  })
  
  circular_bar = reactive({
    cir_bar_func(input$input_gene,input_dt = dt_forest())
  })
  
  output$circular_bar_plt = renderPlot({
    tryCatch({
      circular_bar()
    },error = function(e){
      ""
    })

  })
  
  ##################################### string network ######################################
  src= reactive({
    pan_server_fig(box_id="box4",input_gene = input[["input_gene"]])
  })
  
  output$string_net<-  renderUI({
    tags$img(src = src(),alt="string network",width=1000,height=400)
  }) 
  
  #===============================================================================
  # Multivariable analysis
  #===============================================================================
  
  observe({
    # List the CSV files in the 'data' directory
    data_dir <- paste0(folder_main,"www/data/geneNames_study/")
    files <- list.files(data_dir, pattern = "genenames_.*\\.csv$", full.names = F)
    
    input_study =  input$study_type # "GEO159067 (HNSC)"
    study_name_map = map_study_names(input_study)
    geneNames = fread(paste0(data_dir,study_name_map,".csv"))
    sel_input_server(input,output,session,"input_gene_dash","Select a gene",geneNames$geneNames,input$input_gene)  
  })
  
  # show clinical covariates
  observe({
    dt = pheno()
    filtered_data <- dt %>% filter(study == input$study_type)
    # filtered_data = data.frame(filtered_data[,1:(ncol(filtered_data)-1)])
    # Identify covariates that are not all missing
    covariates <- colnames(filtered_data)[apply(filtered_data, 2, function(col) any(!is.na(col)))]
    
    # Exclude non-covariates
    non_covariates <- c("ID", "PFS.time", "PFS", "OS.time", "OS", "study","CYT","Hypoxia_ssgsea","ISG.RS_ssgsea","IFN.Hallmark_ssgsea")
    covariates <- setdiff(covariates, non_covariates)
    
    # Update selectizeInput choices
    updateSelectizeInput(session, "adj_cov_multi", choices = covariates,select = covariates[1], server = TRUE)
  })
  
  # show survtype
  observe({
    dt = pheno()
    filtered_data <- dt %>% filter(study == input$study_type)
    filtered_data = data.frame(filtered_data[,c("OS","PFS")])
    covariates <- colnames(filtered_data)[apply(filtered_data, 2, function(col) any(!is.na(col)))]
    
    # Update selectizeInput choices
    updateSelectizeInput(session, "survtype", choices = covariates, server = TRUE)
  })
  
  ############################# unadjusted KM and cox ##############################
  ########## KM plot #############
  output$survival_dash = renderUI({
    
    div(
      survival_UI(headings=input$study_type,output_info="survival_km",
                  more_cox= "morecox"),
      style="display: inline-block; ")
    
  })
  
  input_data_unadjsurv <- reactive({
    dt = read_fst(paste0(folder_main,"www/data/dt_ICI.fst"),columns=c( "age","gender","clinical_stage","alcohol","smoking","hpv.status","immunotherapy.line",
                                                                       "anatomical.location","study","OS.time","OS","PFS.time","PFS","FMOne_mutation_burden_per_MB","Neoantigen_burden_per_MB",
                                                                       input$input_gene_dash))
    names(dt) = make.names(names(dt))
    dt
  })
  
  input_data_adj = reactive({
    sel_gene = make.names(input$input_gene_dash)
    dt = input_data_unadjsurv()
    pheno =pheno()
    
    dt_study = subset(dt,dt$study==input$study_type)
    pheno_study = subset(pheno,pheno$study==input$study_type)
    
    if(input$survtype=="OS"){
      dt_study = data.frame(dt_study[,c("OS.time","OS",sel_gene)])
    }else{
      dt_study = data.frame(dt_study[,c("PFS.time","PFS",sel_gene)])
      
    }
    names(dt_study) = c("time","status","gene")
    
    adj_dt = survival_adj_data(input$cutoff_multi,dt_study,pheno_study)
    return(adj_dt)
  })
  
  ###### KM plot ###########
  unadj_km = reactive({
    tryCatch({
      if(input$km_type=="KM"){
        tryCatch({
          dt = input_data_unadjsurv()
          dt_study = subset(dt,dt$study==input$study_type)
          if(input$survtype=="OS"){
            dt_study = data.frame(dt_study[,c("OS.time","OS",make.names(input$input_gene_dash))])
          }else if(input$survtype=="PFS/PFI"){
            dt_study = data.frame(dt_study[,c("PFS.time","PFS",make.names(input$input_gene_dash))])
          }
          names(dt_study) = c("time","status","gene")
          
          survival_server1(geneName=make.names(input$input_gene_dash),cutoff=input$cutoff_multi,data=dt_study)
          
        },error = function(e){
          "Create an empty image with 'please input a valid geneName'"
        })
      }else{
        tryCatch({
          
          survival_adj_server1(inputdata=input_data_adj(),adj_cov=input$adj_cov_multi)
          
        },error = function(e){
          "Create an empty image with 'please input a valid geneName'"
        })
      }
      
    },error=function(e){
      ""
    })
  })
  
  output$survival_km = renderPlot({
    unadj_km()
  })
  
  ########## coxph result ##############
  
  unadj_cox = reactive({
    # tryCatch({
    if(input$km_type=="KM"){
      #  tryCatch({
      dt = data.frame(input_data_unadjsurv())
      dt_study = subset(dt,dt$study==input$study_type)
      if(input$survtype=="OS"){
        dt_study = data.frame(dt_study[,c("OS.time","OS",make.names(input$input_gene_dash))])
      }else if(input$survtype=="PFS/PFI"){
        dt_study = data.frame(dt_study[,c("PFS.time","PFS",make.names(input$input_gene_dash))])
      }
      names(dt_study) = c("time","status","gene")
      
      res =  survival_server2(geneName=make.names(input$input_gene_dash),input_data=dt_study)
      names(res) = c("Gene","HR (95% CI)","p_value","n(events)")
      
      return(data.frame(res,check.names=F))
      # },error = function(e){
      #   "Create an empty image with 'please input a valid geneName'"
      # })
    }else{
      input_data_adj()
      tryCatch({
        res = survival_adj_server2(inputdata=input_data_adj(),adj_cov=input$adj_cov_multi)
        res$Variable[nrow(res)] = make.names(input$input_gene_dash)
        res
      },error = function(e){
        "Create an empty image with 'please input a valid geneName'"
      })
    }
    
    # },error=function(e){
    #   ""
    # })
  })
  
  output$morecox = renderDT({
    res = unadj_cox()
    DT:::datatable(
      data.frame(res),rownames = FALSE,colnames=FALSE,options =  list(scrollX=TRUE,ordering=F,dom="ft",pageLength=1000,paging= F,scrollCollapse=T,scrollY="50vh",
                                                                      initComplete = JS(
                                                                        "function(settings, json) {",
                                                                        "$(this.api().table().header()).css({'background-color': '#4B4B4B', 'color': '#fff'});",
                                                                        "}"),
                                                                      headerCallback = JS(
                                                                        "function(thead, data, start, end, display){",
                                                                        "  $(thead).css({'background-color': '#4B4B4B', 'color': '#fff'});",
                                                                        "}")
      ),
      container = tags$table(
        class="stripe row-border hover",
        tags$thead(tags$tr(lapply(colnames(res), tags$th)))
      )
    )
  })

  #===============================================================================
  # Gene-hallmark interaction
  #===============================================================================
  observe({
    # List the CSV files in the 'data' directory
    data_dir <- paste0(folder_main,"www/data/geneNames_study/")
    files <- list.files(data_dir, pattern = "genenames_.*\\.csv$", full.names = F)
    
    input_study =  input$study_type_ghi # "GEO159067 (HNSC)"
    study_name_map = map_study_names(input_study)
    geneNames = fread(paste0(data_dir,study_name_map,".csv"))
    sel_input_server(input,output,session,"gene_single_ghi","Select a gene",geneNames$geneNames,input$input_gene)  
  })
  
  # show survtype
  observe({
    dt = pheno()
    filtered_data <- dt %>% filter(study == input$study_type_ghi)
    filtered_data = data.frame(filtered_data[,c("OS","PFS")])
    covariates <- colnames(filtered_data)[apply(filtered_data, 2, function(col) any(!is.na(col)))]
    
    # Update selectizeInput choices
    updateSelectizeInput(session, "survtype_ghi", choices = covariates, server = TRUE)
  })
  
  # Hallmarks options
  observe({
    dt = pheno()
    filtered_data <- dt %>% filter(study == input$study_type_ghi)
    
    covariates <- colnames(filtered_data)[apply(filtered_data, 2, function(col) any(!is.na(col)))]
    
    # Exclude non-covariates
    non_covariates <- c("ID", "PFS.time", "PFS", "OS.time", "OS", "study","anatomical.location","hpv.status")
    covariates <- setdiff(covariates, non_covariates)
    covariates = c(covariates,"AR")
    # Update selectizeInput choices
    updateSelectizeInput(session, "hallmark_match", choices = covariates,select = covariates, server = TRUE)
  })
  
  ############ data for GHI model #############
  # GE data
  dt_match = reactive({
    if(length(input$study_type_ghi)>1){
      validate("Please select one cancer type at a time and click 'Run' button")
    }
    cov_list = unique(c(make.names(input$gene_single_ghi),make.names(input$hallmark_match)))
    dt = read_fst(paste0("www/data/","dt_ICI.fst"),columns=c("study","ID","OS.time","OS", "PFS.time","PFS",cov_list))
    dt = subset(dt,dt$study==input$study_type_ghi)
    dt
  })
  
  ############ edge plot ##############
  res_edge = eventReactive(input$search_single_ghi,{
    if(length(input$study_type_ghi)>1){
      validate("Please select one cancer type at a time")
    }
    
    ghi_model(cov=input$hallmark_match,input_gene=input$gene_single_ghi,hallmark_dt = dt_match(),response_match=input$survtype_ghi,tumor_purity_adj="No")
  },ignoreNULL = F)
  
  output$GHI_zscore = renderDT({
    res = t(res_edge())
    DT:::datatable(
      data.frame(res),options = list(searchable =F)
    )    %>%
      formatSignif(columns = c("Partial.model","Full.model"), digits = 3)
  })
  
  output$edge_plot = renderPlot({
    par(cex=0.5)
    chordDiagram( res_edge(), annotationTrack = "grid",  transparency = 0.6,
                  link.visible=abs(res_edge())>1.96,
                  preAllocateTracks = list(track.height = 0.2))
    
    circos.track(track.index = 1, panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, cex=1.5,
                  facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
    }, bg.border = NA) # here set bg.border to NA is important
  })
  #==================================================================================
  ##################################################################################
  ####################### TWO-gene ###############################
  ##################################################################################
  #==================================================================================
  #------------------------------------------------------------------------
  ################### User input ################################
  #------------------------------------------------------------------------
  observe({
    # List the CSV files in the 'data' directory
    data_dir <- paste0(folder_main,"www/data/geneNames_study/")
    files <- list.files(data_dir, pattern = "genenames_.*\\.csv$", full.names = F)
    
    input_study =  input$study_type_ghi # "GEO159067 (HNSC)"
    study_name_map = map_study_names(input_study)
    geneNames = fread(paste0(data_dir,study_name_map,".csv"))
    sel_input_server(input,output,session,"gene_dual","Select two genes",geneNames$geneNames,c("PIK3CA","TP53"))  
  })
  
  # show survtype
  observe({
    dt = pheno()
    filtered_data <- dt %>% filter(study == input$study_type_dual)
    filtered_data = data.frame(filtered_data[,c("OS","PFS")])
    covariates <- colnames(filtered_data)[apply(filtered_data, 2, function(col) any(!is.na(col)))]
    
    # Update selectizeInput choices
    updateSelectizeInput(session, "survtype_dual", choices = covariates, server = TRUE)
  })
  
  #--------------------------------------------------------------------------------
  ###################### Univariate cox model ############################
  #-------------------------------------------------------------------------------
  raw_dt = reactive({
    req(input$gene_dual, input$study_type_dual)
    dt = lapply(1:2,function(x){
      
      dt = read_fst(paste0(folder_main,"www/data/dt_ICI.fst"),
                    columns = c("OS","OS.time","PFS","PFS.time","study",input$gene_dual[x]))
      
      names(dt)[ncol(dt)] = make.names(input$gene_dual[x])
      dt = subset(dt,dt$study==input$study_type_dual)
      dt
    })
    
    ###################### calculate ratio  ###############################
    dt = data.frame(do.call(cbind,dt))
    
    dt$ratio = dt[,make.names(input$gene_dual[1])]/dt[,make.names(input$gene_dual[2])]
    #dt$ratio  = log2( dt$ratio +1)
    #dt[,make.names(input$gene_dual)] = log2(dt[,make.names(input$gene_dual)]+1)
    #dt$OS.time = dt$OS.time/30.417
    #dt$PFI.time= dt$PFI.time/30.417
    dt = expss::apply_labels(dt,ratio = paste0((input$gene_dual[1]),"/",(input$gene_dual[2])))
    
    
    ###################### calculate ratio group ###############################
    if(input$KM_cutoff=="median"){
      dt$ratio_group =  findInterval(dt$ratio,median(dt$ratio,na.rm=T))   
      dt$ratio_group = ifelse(dt$ratio_group==1,"high","low")
      
    } else if(input$KM_cutoff=="quartile"){
      dt$ratio_group =  findInterval(dt$ratio,quantile(dt$ratio,na.rm=T))   
      dt$ratio_group = ifelse(dt$ratio_group%in%c(2,3),NA,ifelse(dt$ratio_group%in%c(4,5),"high","low"))
      
    } else if(input$KM_cutoff=="optimal"){
      res.cut <- surv_cutpoint(dt,time = "OS.time", event = "OS", variables="ratio")$cutpoint[[1]]
      dt$ratio_group  = ifelse(dt$ratio >=res.cut ,"high","low") 
    }
    
    dt = expss::apply_labels(dt,ratio_group = paste0((input$gene_dual[1]),"/",(input$gene_dual[2])))
    
    dt
  })
  
  
  
  uni_cox_fig1 = reactive({
    
    req(input$gene_dual, input$survtype_dual)
    
    if(is.null(input$run_dual)){
      uni_cox = lapply(1:2,function(x){
        formula =  as.formula(paste('Surv(OS.time, OS)~',make.names(input$gene_dual[x])))
        fit = coxph(formula,data=raw_dt())
        
        
        
        forestmodel::forest_model(fit)
      })
      
      
    }else{
      if(input$survtype_dual=="OS"){
        uni_cox = lapply(1:2,function(x){
          formula =  as.formula(paste('Surv(OS.time, OS)~',make.names(input$gene_dual[x])))
          fit = coxph(formula,data=raw_dt())
          
          
          
          forestmodel::forest_model(fit)
        })
        
      }else if (input$survtype_dual=="PFS"){
        uni_cox = lapply(1:2,function(x){
          formula =  as.formula(paste('Surv(PFS.time, PFS)~',make.names(input$gene_dual[x])))
          fit = coxph(formula,data=raw_dt())
          forestmodel::forest_model(fit)
          
        })
        
      }
    }
    
    
    
  })
  
  
  output$uni_cox_plot1 = renderPlot({
    # tryCatch({
    uni_cox_fig1()[[1]]
    
    # },error=function(e){
    #   ""
    # })
  })
  
  output$uni_cox_plot2 = renderPlot({
    #  tryCatch({
    uni_cox_fig1()[[2]]
    # },error= function(e){
    #   ""
    # })
  })
  
  #===============================================================================
  # unadj KM and logrank test #
  #=============================================================================== 
  input_data_KM_dual <- reactive({
    dt = read_fst(paste0(folder_main,"www/data/dt_ICI.fst"),columns=c( "age","gender","clinical_stage","alcohol","smoking","hpv.status","immunotherapy.line",
                                                                       "anatomical.location","study","OS.time","OS","PFS.time","PFS","FMOne_mutation_burden_per_MB","Neoantigen_burden_per_MB",input$gene_dual))
    names(dt) = make.names(names(dt))
    dt
  })
  
  KM_input_dt = reactive({
    sel_gene = make.names(input$gene_dual)
    dt = input_data_KM_dual()
    pheno =pheno()
    
    dt_study = subset(dt,dt$study==input$study_type_dual)
    pheno_study = subset(pheno,pheno$study==input$study_type_dual)
    
    if(input$survtype_dual=="OS"){
      dt_study1 = data.frame(dt_study[,c("OS.time","OS",sel_gene[1])])
      dt_study2 = data.frame(dt_study[,c("OS.time","OS",sel_gene[2])])
    }else{
      dt_study1 = data.frame(dt_study[,c("PFS.time","PFS",sel_gene[1])])
      dt_study2 = data.frame(dt_study[,c("PFS.time","PFS",sel_gene[2])])
    }
    names(dt_study1) = c("time","status","gene")
    names(dt_study2) = c("time","status","gene")
    
    adj_dt1 = survival_adj_data(input$KM_cutoff,dt_study1,pheno_study)
    adj_dt2 = survival_adj_data(input$KM_cutoff,dt_study2,pheno_study)
    
    names(adj_dt1)[3] = sel_gene[1]
    names(adj_dt2)[3] = sel_gene[2]
    
    
    KM_dt = data.frame(cbind(adj_dt1[,c("time","status",sel_gene[1])],adj_dt2[,sel_gene[2]]))
    names(KM_dt)[4] = sel_gene[2]
    KM_dt$gene_combine = paste0(KM_dt[[sel_gene[1]]],"_",KM_dt[[sel_gene[2]]])
    return(KM_dt)
  })
  
  # KM plot
  uni_km_plot = reactive({
    req(input$gene_dual, input$survtype_dual)
    
    km = lapply(1:2,function(x){
      fit <- eval(parse(text = paste0("survfit(Surv(time,status) ~ ", make.names(input$gene_dual[x]), ", data = KM_input_dt())")))
      p = ggsurvplot(fit,data = KM_input_dt(),pval=T,legend.title=input$gene_dual[x],legend.labs=c("high","low"),palette = c("#DF8F44FF","#374E55FF"),size=1)
      p
    })
    
    
    km
    ggarrange(km[[1]]$plot,km[[2]]$plot,ncol=2)
  })
  
  output$uni_km = renderPlot({
    tryCatch({
      uni_km_plot()
      
    },error = function(e){
      ""
    })
  })
  
  # KM table
  uni_km_table_res = reactive({
    tryCatch({
      
      OS_KM_fit_table <- purrr::map_dfr( make.names(input$gene_dual), run_pretty_km_output, model_data = KM_input_dt(), time_in = "time", 
                                         event_in = "status", event_level = '1') %>% 
        select(Group, Level, everything())
      data.frame(OS_KM_fit_table)
    },error = function(e){
      data.table()
    })
    
  })
  
  output$uni_km_table = renderDT({
    tryCatch({
      uni_km_table_res()
      DT:::datatable(
        uni_km_table_res(),rownames = FALSE,colnames=F, options = list(dom='t',scrollX=TRUE,ordering=F,
                                                                       initComplete = JS(
                                                                         "function(settings, json) {",
                                                                         "$(this.api().table().header()).css({'background-color': '#4B4B4B', 'color': '#fff'});",
                                                                         "}"),
                                                                       headerCallback = JS(
                                                                         "function(thead, data, start, end, display){",
                                                                         "  $(thead).css({'background-color': 'white', 'color': '#4B4B4B'});",
                                                                         "}")
        ),
        container = tags$table(
          class="stripe row-border hover",
          tags$thead(tags$tr(lapply(colnames(uni_km_table_res()), tags$th)))
        )
      )
    },error = function(e){
      ""
    })
    
  })
  
  #------------------------------------------------------------------------------
  #################### Gene ratio cox ######################################
  #------------------------------------------------------------------------------
  ratio_cox = reactive({
    #eventReactive(input$run_dual,{
    if(input$survtype_dual=="OS"){
      formula =  as.formula(paste('Surv(OS.time, OS)~ratio'))
    }else{
      formula =  as.formula(paste('Surv(PFS.time, PFS)~ratio'))
    }
    fit = coxph(formula,data= raw_dt())
  })
  
  output$cox_ratio = renderPlot({
    tryCatch({
      forestmodel::forest_model(ratio_cox())
      
    },error = function(e){
      ""
    })
  })
  
  #------------------------------------------------------------------------------
  #################### Gene ratio km ######################################
  #------------------------------------------------------------------------------
  ratio_km = reactive({
    # eventReactive(input$run_dual,{
    if(input$survtype_dual=="OS"){
      formula =  as.formula(paste('Surv(OS.time, OS)~ratio_group'))
    }else{
      formula =  as.formula(paste('Surv(PFS.time, PFS)~ratio_group'))
    }
    fit = survfit(formula,data = raw_dt())
    fit$call$formula <- formula
    
    fit
  })
  
  output$km_ratio = renderPlot({
    tryCatch({
      p = ggsurvplot(ratio_km(),data = raw_dt(),
                     pval=T,legend.title=paste0(input$gene_dual[1],"/",input$gene_dual[2]),
                     legend.labs=c("high","low"),palette = c("#DF8F44FF","#374E55FF"),size=1)
      p$plot    
    },error=function(e){
      ""
    })
    
    
  })
  
  # KM table
  ratio_km_table_res = reactive({
    tryCatch({
      # eventReactive(input$run_dual,{
      if(input$survtype_dual=="OS"){
        time = "OS.time";status = "OS"
      }else{
        time = "PFS.time";status="PFS"
      }
      OS_KM_fit_table <- purrr::map_dfr("ratio_group", run_pretty_km_output, model_data = raw_dt(), time_in = time, 
                                        event_in = status, event_level = '1') %>% 
        select(Group, Level, everything())
      data.frame(OS_KM_fit_table)
    },error = function(e){
      data.table()
    })
    
  })
  
  output$ratio_km_table = renderDT({
    DT:::datatable(
      ratio_km_table_res(),rownames = FALSE,colnames=F,options = list(dom='t',scrollX=TRUE,ordering=F,
                                                                      initComplete = JS(
                                                                        "function(settings, json) {",
                                                                        "$(this.api().table().header()).css({'background-color': '#4B4B4B', 'color': '#fff'});",
                                                                        "}"),
                                                                      headerCallback = JS(
                                                                        "function(thead, data, start, end, display){",
                                                                        "  $(thead).css({'background-color': '#4B4B4B', 'color': '#fff'});",
                                                                        "}")
      ),
      container = tags$table(
        class="stripe row-border hover",
        tags$thead(tags$tr(lapply(colnames(ratio_km_table_res()), tags$th)))
      )
    )
  })
  
  #-------------------------------------------------------------------------------
  ################# Mutivariable cox ####################################
  #-------------------------------------------------------------------------------
  
  multi_cox =reactive({
    #eventReactive(input$run_dual,{
    if(input$survtype_dual=="OS"){
      formula = as.formula(paste0("Surv(OS.time,OS)~",paste(make.names(input$gene_dual),collapse = "+")))
    }else{
      formula = as.formula(paste0("Surv(PFS.time,PFS)~",paste(make.names(input$gene_dual),collapse = "+")))
    }
    fit = coxph(formula,data = raw_dt())
    fit
  })
  
  output$multi_cox = renderPlot({
    tryCatch({
      forestmodel::forest_model(multi_cox())
      
    },error=function(e){
      ""
    })
  })
  
  #-------------------------------------------------------------------------------
  ################ KM with 4 groups ########################
  #-------------------------------------------------------------------------------
  ########## KM plot ###########
  
  cutoff_dt_multi = reactive({
    na.omit(KM_input_dt())
  })
  
  multi_km_plot = reactive({
    #eventReactive(input$run_dual,{
    formula = as.formula("Surv(time,status)~gene_combine")
    fit = survfit(formula,data = cutoff_dt_multi())
    fit$call$formula <- formula
    if(length(unique(cutoff_dt_multi()$gene_combine))==2){
      p = ggsurvplot(fit,data = cutoff_dt_multi(),pval=T,legend.title=paste0(input$gene_dual[1],"+",input$gene_dual[2]),
                     legend.labs=c("high+high","low+low"),palette = c("#DF8F44FF","#374E55FF"))
    }else{
      p = ggsurvplot(fit,data = cutoff_dt_multi(),pval=T,legend.title=paste0(input$gene_dual[1],"+",input$gene_dual[2]),
                     legend.labs=c("high+high","high+low","low+high","low+low"),palette = c("#DF8F44FF","gray","#4682B433","#374E55FF"))+
        guides(colour = guide_legend(nrow = 2))
    }
    
    p$plot
  })
  
  output$multi_km = renderPlot({
    tryCatch({
      multi_km_plot()
      
    },error = function(e){
      ""
    })
  })
  
  ###############KM 4-group log-rank results ###############
  multi_km_table_res = reactive({
    tryCatch({
      
      OS_KM_fit_table <- purrr::map_dfr( "gene_combine", run_pretty_km_output, model_data = cutoff_dt_multi(), time_in = "time",
                                         event_in = "status", event_level = '1') %>%
        select(Group, Level, everything())
      data.frame(OS_KM_fit_table)
    },error = function(e){
      data.table()
    })
    
  })
  
  output$multi_km_table = renderDT({
    DT:::datatable(
      multi_km_table_res(),rownames = FALSE,colnames=F,options = list(dom='t',scrollX=TRUE,ordering=F,
                                                                      initComplete = JS(
                                                                        "function(settings, json) {",
                                                                        "$(this.api().table().header()).css({'background-color': '#4B4B4B', 'color': '#fff'});",
                                                                        "}"),
                                                                      headerCallback = JS(
                                                                        "function(thead, data, start, end, display){",
                                                                        "  $(thead).css({'background-color': '#4B4B4B', 'color': '#fff'});",
                                                                        "}")
      ),
      container = tags$table(
        class="stripe row-border hover",
        tags$thead(tags$tr(lapply(colnames(multi_km_table_res()), tags$th)))
      )
    )
  })
  
  ###############KM 4-group pairwise log-rank results ###############
  pairwise_km = reactive({
    tryCatch({
      #eventReactive(input$run_dual,{
      formula = as.formula("Surv(time,status)~gene_combine")
      
      res = pairwise_survdiff(formula,data=cutoff_dt_multi())
      res = data.frame(broom::tidy(res))
      res$p.value = ifelse(res$p.value>=0.001,round(res$p.value,3),
                           formatC(res$p.value,format = "e",digits=3))
      names(res)[1] = input$gene_dual[1]
      names(res)[2] = input$gene_dual[2]
      res
    },error = function(e){
      data.table()
    })
    
  })
  
  output$pairwise_logrank = renderDT({
    
    DT:::datatable(
      pairwise_km(),rownames = FALSE,colnames=F,options = list(dom='t',scrollX=TRUE,ordering=F,
                                                               initComplete = JS(
                                                                 "function(settings, json) {",
                                                                 "$(this.api().table().header()).css({'background-color': '#4B4B4B', 'color': '#fff'});",
                                                                 "}"),
                                                               headerCallback = JS(
                                                                 "function(thead, data, start, end, display){",
                                                                 "  $(thead).css({'background-color': '#4B4B4B', 'color': '#fff'});",
                                                                 "}")
      ),
      container = tags$table(
        class="stripe row-border hover",
        tags$thead(tags$tr(lapply(colnames(pairwise_km()), tags$th)))
      )
    )
  })
  #-------------------------------------------------------------------------------
  ################### Multivariable model with interaction ###############
  #-------------------------------------------------------------------------------
  ########### cox with interaction #############
  cox_interaction = reactive({
    tryCatch({
      #eventReactive(input$run_dual,{
      if(input$survtype_dual=="OS"){
        formula = as.formula(paste0("Surv(OS.time,OS)~",paste(make.names(input$gene_dual),collapse = "*")))
      }else{
        formula = as.formula(paste0("Surv(PFS.time,PFS)~",paste(make.names(input$gene_dual),collapse = "*")))
      }
      cox_fit = coxph(formula,data=raw_dt())
      res = tidy(cox_fit,exponentiate = T,conf.int = T)
      res = data.frame(res[,c("term","estimate","conf.low","conf.high","p.value")])
      
      res$p.value = ifelse(res$p.value>=0.001,round(res$p.value,3),
                           formatC(res$p.value,format = "e",digits=3))
      
      for(i in 2:4){
        res[,i] = round(res[,i],3)
      }
      names(res) = c("term","HR","conf.low", "conf.high", "p.value")
      res
    },error =function(e){
      data.table()
    })
    
  })
  
  output$cox_interaction = renderDT({
    DT:::datatable(
      cox_interaction(),rownames = FALSE,colnames=F,options = list(dom='t',scrollX=TRUE,ordering=F,
                                                                   initComplete = JS(
                                                                     "function(settings, json) {",
                                                                     "$(this.api().table().header()).css({'background-color': '#4B4B4B', 'color': '#fff'});",
                                                                     "}"),
                                                                   headerCallback = JS(
                                                                     "function(thead, data, start, end, display){",
                                                                     "  $(thead).css({'background-color': '#4B4B4B', 'color': '#fff'});",
                                                                     "}")
      ),
      container = tags$table(
        class="stripe row-border hover",
        tags$thead(tags$tr(lapply(colnames(cox_interaction()), tags$th)))
      )
    )
  })
  
  ########### KM with interaction #############
  # high and low for gene1 data
  
  km_interaction_hi = reactive({
    dt = KM_input_dt()
    dt1 = subset(dt,dt[[input$gene_dual[1]]]=="high")
    
    dt1
  })
  
  km_interaction_lo = reactive({
    dt = KM_input_dt()
    dt1 = subset(dt,dt[[input$gene_dual[1]]]=="low")
    dt1
    
  })
  
  # KM plot 
  km_interation_plot = reactive({
    #eventReactive(input$run_dual,{
    
    formula = as.formula(paste0("Surv(time,status)~",make.names(input$gene_dual[2])))
    fit1 = survfit(formula,data = km_interaction_hi())
    fit2 = survfit(formula,data = km_interaction_lo())
    fit1$call$formula <- formula
    fit2$call$formula <- formula
    
    p1 = ggsurvplot(fit1,data = km_interaction_hi(),pval=T,legend.title=input$gene_dual[2],legend.labs=c("high","low"),palette = c("#DF8F44FF","#374E55FF")) + ggtitle(paste0(input$gene_dual[1],":High"))
    p2 = ggsurvplot(fit2,data = km_interaction_lo(),pval=T,legend.title=input$gene_dual[2],legend.labs=c("high","low"),palette = c("#DF8F44FF","#374E55FF"))+ ggtitle(paste0(input$gene_dual[1],":Low"))
    ggarrange(p1$plot,p2$plot)
  })
  
  # output$test = renderDT({
  #   km_interaction_hi()
  # })
  output$km_interaction = renderPlot({
    tryCatch({
      km_interation_plot()  
    },error=function(e){
      ""
    })
    
  })
  
  ########################## scatterplot #######################################
  scatter_plot = reactive({
    #eventReactive(input$run_dual,{
    test = raw_dt()
    test_plt = data.frame(test[,c(make.names(input$gene_dual[1]),make.names(input$gene_dual[2]))])
    names(test_plt) = c("gene1","gene2")
    # for(i in 1:2){
    #   test_plt[,i] = log2(test_plt[,i]+1)
    # }
    
    pearson_corr <- round(cor(test_plt$gene1, test_plt$gene2, method = "pearson"), 2)
    spearman_corr <- round(cor(test_plt$gene1, test_plt$gene2, method = "spearman"), 2)
    
    scatter_plot <- ggscatter(test_plt, x = "gene1", y = "gene2",
                              add = "reg.line",                                 # Add regression line
                              conf.int = TRUE,                                  # Add confidence interval
                              add.params = list(color = "blue",
                                                fill = "lightgray")
    )+
      xlab(input$gene_dual[1])+
      ylab(input$gene_dual[2])
    
    # Add annotations for Pearson and Spearman correlations
    y_pos_pearson <- max(test_plt$gene2) - (max(test_plt$gene2) - min(test_plt$gene2)) * 0.1
    y_pos_spearman <- y_pos_pearson - (max(test_plt$gene2) - min(test_plt$gene2)) * 0.1
    
    scatter_plot <- scatter_plot +
      annotate("text", x = min(test_plt$gene1), y = y_pos_pearson, 
               label = paste("Pearson:", pearson_corr), 
               hjust = 0, vjust = 1, size = 5) +
      annotate("text", x = min(test_plt$gene1), y = y_pos_spearman, 
               label = paste("Spearman:", spearman_corr), 
               hjust = 0, vjust = 1, size = 5)
    # Print the plot
    scatter_plot
    
    
  })
  
  output$Scatterplot_dual = renderPlot({
    tryCatch({
      scatter_plot()
      
    },error=function(e){
      ""
    })
  })  
  #================================================================================
  ################################################################################
  ####################### multi-gene panel ###################################
  #================================================================================
  observe({
    dt = pheno()
    filtered_data <- dt %>% filter(study == input$study_type_multi)
    filtered_data = data.frame(filtered_data[,c("OS","PFS")])
    covariates <- colnames(filtered_data)[apply(filtered_data, 2, function(col) any(!is.na(col)))]
    
    # Update selectizeInput choices
    updateSelectizeInput(session, "surv_multi", choices = covariates, server = TRUE)
  }) 
  
  
  All_genes <- reactive({
    req(input$study_type_multi, input$pathway)
    data_dir <- paste0(folder_main,"www/data/geneNames_study/")
    
    study_name_map <- map_study_names(input$study_type_multi)
    gene_file_path <- paste0(data_dir, study_name_map, ".csv")
    
    print(paste("Checking gene file path:", gene_file_path))  # Debugging output
    
    tryCatch({
      if (file.exists(gene_file_path)) {
        geneNames <- fread(gene_file_path)
        all_genes_in_study <- geneNames$geneNames
        print(paste("Loaded genes for study:", input$study_type_multi))  # Debugging output
        all_genes_in_study
      } else {
        showNotification(paste("Gene file for study", input$study_type_multi, "does not exist."), type = "error")
        NULL
      }
    }, error = function(e) {
      showNotification(paste("Error reading gene file for study", input$study_type_multi, ":", e$message), type = "error")
      NULL
    })
  })
  
  rv = reactiveValues(data=NULL)
  
  observe({
    tryCatch({
      req(All_genes())
      
      sel_pathway <- input$pathway
      pathway_file_path <- paste0(folder_main, "www/pathway/", sel_pathway, ".fst")
      
      print(paste("Attempting to read pathway file:", pathway_file_path))  # Debugging output
      
      if (file.exists(pathway_file_path)) {
        print("File exists. Reading the file.")  # Debugging output
        genesets_pca <- read_fst(pathway_file_path)
        rv$genesets_pca = genesets_pca
      } else {
        showNotification(paste("Pathway file for", sel_pathway, "does not exist."), type = "error")
        NULL
      }
    },error = function(e){
      ""
    })
    
  })
  
  
  
  observe({
    req(input$study_type_multi, input$pathway)
    
    genesets_pca = rv$genesets_pca
    geneNames = All_genes()
    # 
    genesets_list <- lapply(genesets_pca, function(col) {
      col[!is.na(col)]
    })
    # 
    valid_genesets <- sapply(genesets_list, function(geneset) {
      all(geneset %in% geneNames)
    })
    # 
    filtered_genesets = names(genesets_pca)[valid_genesets]
    updateSelectizeInput(session,"geneset",choices = filtered_genesets,selected = filtered_genesets[1],server = T)
    
  })
  
  
  
  
  
  observe({
    req(rv$genesets_pca, All_genes(), input$geneset)
    
    genesets_pca <- rv$genesets_pca
    geneNames <- All_genes()
    
    print(paste("Updating gene_multi for geneset:", input$geneset))  # Debugging output
    
    if (input$geneset == "None") {
      updateSelectizeInput(session, "gene_multi", choices = geneNames, server = TRUE)
    } else if (is.null(input$pathway)) {
      updateSelectizeInput(session, "gene_multi", choices = geneNames, server = TRUE)
    } else {
      sel_genes <- genesets_pca[[input$geneset]][!is.na(genesets_pca[[input$geneset]])]
      updateSelectizeInput(session, "gene_multi", choices = geneNames, selected = sel_genes, server = TRUE)
    }
  })
  
  #--------------------------------------------------------------------------------
  # import data
  #--------------------------------------------------------------------------------
  data_multi = reactive({
    dt = read_fst(paste0(folder_main,"www/data/dt_ICI.fst"),columns=c("study","OS.time","OS","PFS.time","PFS",(input$gene_multi)))
    dt = subset(dt,dt$study==input$study_type_multi)
    names(dt) = make.names(names(dt))
    # if(is.null(input$file)){
    #   dt = read_fst(paste0(folder_main,"www/data/dt_ICI.fst"),columns=c("study","OS.time","OS","PFS.time","PFS",make.names(input$gene_multi)))
    #   dt = subset(dt,dt$study==input$study_type_multi)
    #   names(dt) = make.names(names(dt))
    # }else{
    #   req(input$file)
    #   infile <- input$file
    #   if (is.null(infile)){
    #     return(NULL)
    #   }
    #   genes = read.csv(input$file$datapath)
    #   
    #   sel_genes = unique(c(make.names(genes[,1]),make.names(input$gene_multi)))
    #   dt = read_fst(paste0(folder_main,"www/data/dt_ICI.fst"),columns=c("study","OS.time","OS","PFS.time","PFS",sel_genes))
    #   dt = subset(dt,dt$study==input$study_type_multi)
    #   names(dt) = make.names(names(dt))
    # }
    dt
  })
  
  
  # observeEvent(input$reset,{
  #   sel_genes = unique(make.names(input$gene_multi))
  #   
  #   rv$data = read_fst(paste0(folder_main,"www/data/dt_ICI.fst"),columns=c("study","OS.time","OS","PFS.time","PFS",make.names(input$gene_multi)))
  #   
  #   reset('file')
  # })
  
  #---------------------------------------------------------------------------------------
  ###################### Grdient boosting ############################################
  #---------------------------------------------------------------------------------------
  # output$test_multi = renderDT({
  #   gbm_figure()
  #   })
  
  gbm_figure = reactive({
    
    dt =   data_multi()
    
    sel_genes = make.names(names(dt)[c(6:ncol(dt))])
    
    
    if(input$surv_multi=="OS"){
      fit_formula = as.formula(paste("Surv(OS.time,OS)~", paste0(paste(sel_genes, collapse=" + "))))
      dt = na.omit(dt[,c("OS.time","OS",sel_genes)])
      
    }else if(input$surv_multi=="PFS"){
      fit_formula = as.formula(paste("Surv(PFS.time,PFS)~", paste0(paste(sel_genes, collapse=" + "))))
      dt = na.omit(dt[,c("PFS.time","PFS",sel_genes)])
      
    }
    gbm1 = gbm(fit_formula,data = dt,distribution="coxph")
    best.iter <- gbm.perf(gbm1,method="OOB") # returns test set estimate of best number of trees
    #summary(gbm1,n.trees=best.iter,cBars = length(input$gene_gbm))
    dt_p = summary(gbm1)
    dt_p = dt_p[order(dt_p$rel.inf,decreasing=T),]
    
    # geneName = geneName_mrna()
    # geneName$make.names = make.names(geneName$geneName)
    # dt_p = merge(geneName,dt_p,by.x="make.names",by.y="row.names")
    #dt_p = dt_p[!duplicated(dt_p$geneName),]
    #rownames(dt_p) = dt_p$geneName
    #dt_p = dt_p[,-c(1:2)]
    #dt_p$var=rownames(dt_p)
    dt_p = subset(dt_p,dt_p$rel.inf>=input$gb_menu)
    
    ggplot(dt_p,aes(x=reorder(var, rel.inf),y=rel.inf,fill=rel.inf))+
      geom_bar(stat="identity")+
      coord_flip()+
      theme_classic()+
      xlab("Relative influence")+
      ylab("")
  })
  
  
  
  
  output$gbm_plot = renderPlot({
    tryCatch({
      gbm_figure()
    },error=function(e){
      ""
    })
    
  })
  
  ############################ waterfall plot #############################
  #    cox_waterfall_plot = eventReactive(input$search_gbm,{
  cox_waterfall_plot = reactive({
    dt =   data_multi()
    
    sel_genes = make.names(names(dt)[c(6:ncol(dt))])
    
    # if(is.null(input$file)){
    #   #sel_genes =  make.names(input$gene)
    #   if(input$anchor_not=="No"){
    #     sel_genes = make.names(names(dt)[c(6:ncol(dt))])
    #   }else{
    #     length_anchor = length(input$anchor_gene)
    #     sel_genes = make.names(names(dt)[c((5+length_anchor+1):ncol(dt))])
    #     
    #   }
    # }else{
    #   if(input$anchor_not=="No"){
    #     sel_genes = make.names(names(dt)[c(6:ncol(dt))])
    #   }else{
    #     length_anchor = length(input$anchor_gene)
    #     sel_genes = make.names(names(dt)[c((5+length_anchor+1):ncol(dt))])
    #     
    #   }
    # }
    
    
    if(input$surv_multi=="OS"){
      dt = na.omit(dt[,c("OS.time","OS",sel_genes)])
      
      res = lapply(1:length(sel_genes),function(x){
        fit_formula = as.formula(paste("Surv(OS.time,OS)~", sel_genes[x]))
        fit = coxph(fit_formula,data=dt)
        res = data.frame(summary(fit)$coef)
        names(res) = c("coef","HR","se","z","p.value")
        res$gene = paste0(rownames(res),res$mark)
        res
      })
    }else if(input$surv_multi=="PFS"){
      dt = na.omit(dt[,c("PFS.time","PFS",sel_genes)])
      
      res = lapply(1:length(sel_genes),function(x){
        fit_formula = as.formula(paste("Surv(PFS.time,PFS)~", sel_genes[x]))
        fit = coxph(fit_formula,data=dt)
        res = data.frame(summary(fit)$coef)
        names(res) = c("coef","HR","se","z","p.value")
        res$gene = paste0(rownames(res),res$mark)
        res
      })
    }
    
    res = do.call(rbind,res)
    res$logHR = ifelse(res$coef>0,"> 0","<= 0")
    
    if(input$waterfall_menu=="Yes"){
      res = subset(res,res$p.value<=0.05)
      res$label_text = ""
    }else{
      res = res
      res$label_text = ifelse(res$p.value<0.1&res$p.value>0.05,".",ifelse(res$p.value<=0.05&res$p.value>0.01,"*",
                                                                          ifelse(res$p.value<=0.01,"**","")))
    }
    
    ggplot(res, aes(x=reorder(gene, desc(HR)), y=coef, fill=logHR, color=logHR)) +
      theme_classic() %+replace%
      theme(axis.line.x = element_blank(), axis.text.x = element_text(face="bold",angle=90), axis.ticks.x = element_blank(),
            axis.title.y = element_text(face="bold",angle=90)) +
      # coord_cartesian(ylim = c(-1,1))+
      geom_bar(stat="identity", width=0.7, position = position_dodge(width=0.4))+
      xlab("")+
      ylab("log HR")+
      scale_color_manual(values=c("gray","#DF8F44FF"))+
      scale_fill_manual(values=c("gray","#DF8F44FF"))+
      geom_text(aes(label=label_text), vjust=-1.5,color="black",nudge_x = 0,
                size = 5, fontface = "bold", family = "Fira Sans")
    
  })
  
  output$waterfall_plot = renderPlot({
    tryCatch({
      cox_waterfall_plot()
    },error=function(e){
      ""
    })
  })    
  
  #--------------------------------------------------------------------------------
  ################# Joint signature (SSGSEA) #####################################
  #--------------------------------------------------------------------------------
  output$average_km_ui = renderUI({
    div(
      fluidRow(
        
        column(2,
               div(style = "display: flex; align-items: flex-start;",  # Ensure items are vertically centered
                   tags$div("KM cutoff",style = "margin-right: 5px; font-size: 16px; white-space: nowrap; line-height: 1;margin-top:5px;"),  # Keep label styling simple
                   selectizeInput("gsea_km_cutoff_av", NULL, c("median", "optimal", "quartile"), "median", options = list(width = '100%'))
               )
        ),
        
        column(2,  # Keep the button aligned
               actionButton("search_average_joint", "Click to Run", style = "width: 50%; height: 38px; background-color: #FFD700; color: black; border: none;")
        )
      )
      
      
    )
    
  })
  
  output$gsea_km_ui = renderUI({
    div(
      fluidRow(
        
        column(2,
               div(style = "display: flex; align-items: flex-start;",  # Ensure items are vertically centered
                   tags$div("KM cutoff",style = "margin-right: 5px; font-size: 16px; white-space: nowrap; line-height: 1;margin-top:5px;"),  # Keep label styling simple
                   selectizeInput("gsea_km_cutoff", NULL, c("median", "optimal", "quartile"), "median", options = list(width = '100%'))
               )
        ),
        column(2,  # Keep the button aligned
               actionButton("search_gsea", "Click to Run", style = "width: 50%; height: 38px; background-color: #FFD700; color: black; border: none;")
        )
      )
      
    )
    
  })
  
  average_table = eventReactive(input$search_average_joint,{
    
    
    dt =  data_multi()
    dt_sub = subset(dt,dt$study==input$study_type_multi)
    dt_sub_gene = dt_sub[,-c(1:5)]
    
    dt_sub$value = rowMeans(dt_sub_gene,na.rm=T)
    dt_sub = data.frame(cbind(dt_sub[,1:5],"value"=dt_sub$value))
    dt_sub
  })
  
  average_table_group = reactive({
    dt_sub = average_table()
    if(input$gsea_km_cutoff_av=="median"){
      dt_sub$cut = findInterval(dt_sub$value,median(dt_sub$value,na.rm=T))
      dt_sub$cut = ifelse(is.na(dt_sub$cut),NA,ifelse(dt_sub$cut ==1,"high","low"))
    }else if(input$gsea_km_cutoff_av=="quartile"){
      dt_sub$cut = findInterval(dt_sub$value, quantile(dt_sub$value,na.rm=T))
      dt_sub$cut = ifelse(dt_sub$cut%in%c(2,3),NA,ifelse(dt_sub$cut%in%c(4,5),"high","low"))
    }else if(input$gsea_km_cutoff_av=="optimal"){
      if(input$gsea_km_type_av=="OS"){
        res.cut <- surv_cutpoint(dt_sub,time = "OS.time", event = "OS", variables="value")
      }else{
        res.cut <- surv_cutpoint(dt_sub,time = "PFI.time", event = "PFI", variables="value")
      }
      res.cat <- surv_categorize(res.cut)
      dt_sub = data.frame(cbind(dt_sub,"cut"=res.cat$value))
    }
    dt_sub
  })
  # output$test_path = renderDT({
  #   average_table_group()
  # })
  average_km_plot  = reactive({
    
    dt = average_table_group()
    
    #  if(input$gsea_cov_av=="None"){
    if(input$surv_multi=="OS"){
      fit = survfit(Surv(OS.time,OS)~cut,data = dt)
    }else{
      fit = survfit(Surv(PFI.time,PFI)~cut,data = dt)
    }
    p = ggsurvplot(fit,data = dt,pval=TRUE,legend.title="Median",legend.labs=c("high","low"),palette = c("#DF8F44FF","#374E55FF"))
    p = p$plot
    # }else{
    #   dtt = data.frame(na.omit(dt))
    #   if(input$gsea_km_type_av=="OS"){
    #     if(input$gsea_cov_av=="Age + Sex"){
    #       res = adjusted_KM(data=dtt,time='OS.time',status="OS",group="cut",covlist=c("age","gender"),stratified_cox = "Yes",reference_group="G&B")
    #     }else if(input$gsea_cov_av=="Age + Sex + tumor_Purity"){
    #       res = adjusted_KM(data=dtt,time='OS.time',status="OS",group="cut",covlist=c("age","gender","purity"),stratified_cox = "Yes",reference_group="G&B")
    #     }
    #   }else{
    #     if(input$gsea_cov_av=="Age + Sex"){
    #       res = adjusted_KM(data=dtt,time='PFI.time',status="PFI",group="cut",covlist=c("age","gender"),stratified_cox = "Yes",reference_group="G&B")
    #     }else if(input$gsea_cov_av=="Age + Sex + tumor_Purity"){
    #       res = adjusted_KM(data=dtt,time='PFI.time',status="PFI",group="cut",covlist=c("age","gender","purity"),stratified_cox = "Yes",reference_group="G&B")
    #     }
    #   }
    #
    #   p =  ggplot(res,aes(x=time,y = prob, group =class))+
    #     geom_step(aes(color = class),size=1.5)+
    #     theme_classic()+
    #     ylim(c(0,1))+
    #     theme(legend.position="top")+
    #     scale_color_manual(values=c("#DF8F44FF", "#374E55FF"))
    #
    # }
    p
  })
  output$average_km_plot = renderPlot({
    average_km_plot()
  })
  
  output$average_km_table = renderDT({
    res = average_table()
    #names(res)[2] = "average_value"
    #res = data.frame(res[,c("Row.names","value","purity","OS","OS.time","PFI","PFI.time","age","gender")])
    #names(res)[2]="average_value"
    DT:::datatable(
      res,
      rownames = FALSE,colnames=FALSE,
      options = list(scrollX=TRUE,ordering=F,dom="ft",pageLength=1000,paging= F,scrollCollapse=T,scrollY="50vh",
                     initComplete = JS(
                       "function(settings, json) {",
                       "$(this.api().table().header()).css({'background-color': '#4B4B4B', 'color': '#fff'});",
                       "}"),
                     headerCallback = JS(
                       "function(thead, data, start, end, display){",
                       "  $(thead).css({'background-color': '#4B4B4B', 'color': '#fff'});",
                       "}")
      ),
      container = tags$table(
        class="stripe row-border hover",
        tags$thead(tags$tr(lapply(colnames(res), tags$th)))
        
      )
    )
  })
  
  observe({
    output[["download_average_table"]]<-downloadHandler(
      filename = function() {
        paste0(input$cancer_type,"_joint_average",'.csv', sep = '')
        
      },
      content = function(file){
        fwrite(average_table(),file,row.names = T)
      }
    )
    
  })
  
  observe({
    output[["download_average_KM"]]<-downloadHandler(
      filename = function() {
        paste0(input$surv_multi,"_joint_average_KM",'.png', sep = '')
        
      },
      content = function(file){
        ggsave(file, average_km_plot(),type="cairo-png",width = 12, height = 12,units="cm")
      }
    )
  })
  
  #===========================================================================================
  ######################## SSGSEA ##################################################
  #===========================================================================================
  gsea_data = reactive({
    dt = read_fst(paste0(folder_main,"www/data/dt_ICI.fst"))
    dt = subset(dt,dt$study==input$study_type_multi)
    names(dt) = make.names(names(dt))
    dt
  })
  
  gsea_table = eventReactive(input$search_gsea,{
    validate(need(input$study_type_multi !="GEO159067 (HNSC)", "This study is not table to perform ssGSEA"))
    
    validate(need(input$search_gsea !="", "Please wait while calculations are running....."))
    
    progress = shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message="Estimating ssGSEA scores for the gene sets",
                 detail = "This will take ~1 min")
    
    progress$inc(1/5, message=paste("Import data"))
    
    dt = gsea_data()
    
    dt_sub = dt
    dt_sub = data.frame(dt_sub[,c(1,22:ncol(dt_sub))]) ############### might change !!!!!!!!!!!!!!!
    rownames(dt_sub) = dt_sub$ID
    dt_sub = data.frame(dt_sub[,-1])
    
    #dt_sub = na.omit(dt_sub)
    
    dt_sub = data.frame(t(dt_sub))
    #  list(c(input$gene_gsea))
    GE_matrix <- as.matrix(dt_sub)
    GE_matrix = na.omit(GE_matrix)
    
    
    
    gene_gsea = data_multi()
    gene_gsea = gene_gsea[,-c(1:5)]
    
    
    progress$inc(2/5, message=paste("Estimating ssGSEA scores for the gene sets"))
    gsva_H <- gsva(expr= GE_matrix, list(names(gene_gsea)), method="ssgsea")
    progress$inc(1/5, message=paste("Combining with clinical info"))
    gsva_H = data.frame(t(gsva_H))
    names(gsva_H)="value"
    gsva_H = merge(gsva_H,dt[,c("ID","OS.time","OS","PFS.time","PFS")],by.x="row.names",by.y="ID")
    gsva_H 
  })
  
  gsea_table_process = reactive({
    dt = gsea_table()
    if(input$gsea_km_cutoff=="median"){
      dt$cut = findInterval(dt$value,median(dt$value,na.rm=T))
      dt$cut = ifelse(is.na(dt$cut),NA,ifelse(dt$cut ==1,"high","low"))
    }else if(input$gsea_km_cutoff=="quartile"){
      dt$cut = findInterval(dt$value, quantile(dt$value,na.rm=T))
      dt$cut = ifelse(dt$cut%in%c(2,3),NA,ifelse(dt$cut%in%c(4,5),"high","low"))
    }else if(input$gsea_km_cutoff=="optimal"){
      if(input$gsea_km_type=="OS"){
        res.cut <- surv_cutpoint(dt,time = "OS.time", event = "OS", variables="value")
      }else{
        res.cut <- surv_cutpoint(dt,time = "PFI.time", event = "PFI", variables="value")
      }
      res.cat <- surv_categorize(res.cut)
      dt = data.frame(cbind(dt,"cut"=res.cat$value))
    }
    dt
  })
  
  # output$test_path = renderDT({
  #   gsea_table_process()})
  
  gsea_km_plot = reactive({
    
    dt = gsea_table_process()
    
    # if(input$gsea_cov=="None"){
    if(input$surv_multi=="OS"){
      fit = survfit(Surv(OS.time,OS)~cut,data = dt)
    }else{
      fit = survfit(Surv(PFS.time,PFS)~cut,data = dt)
    }
    p = ggsurvplot(fit,data = dt,pval=TRUE,legend.title="Median",legend.labs=c("high","low"),palette = c("#DF8F44FF","#374E55FF"))
    p = p$plot
    # }else{
    #   dtt = data.frame(na.omit(dt))
    #   if(input$gsea_km_type=="OS"){
    #     if(input$gsea_cov=="Age + Sex"){
    #       res = adjusted_KM(data=dtt,time='OS.time',status="OS",group="cut",covlist=c("age","gender"),stratified_cox = "Yes",reference_group="G&B")
    #     }else if(input$gsea_cov=="Age + Sex + tumor_Purity"){
    #       res = adjusted_KM(data=dtt,time='OS.time',status="OS",group="cut",covlist=c("age","gender","purity"),stratified_cox = "Yes",reference_group="G&B")
    #     }
    #   }else{
    #     if(input$gsea_cov=="Age + Sex"){
    #       res = adjusted_KM(data=dtt,time='PFI.time',status="PFI",group="cut",covlist=c("age","gender"),stratified_cox = "Yes",reference_group="G&B")
    #     }else if(input$gsea_cov=="Age + Sex + tumor_Purity"){
    #       res = adjusted_KM(data=dtt,time='PFI.time',status="PFI",group="cut",covlist=c("age","gender","purity"),stratified_cox = "Yes",reference_group="G&B")
    #     }
    #   }
    #   
    #   p =  ggplot(res,aes(x=time,y = prob, group =class))+
    #     geom_step(aes(color = class),size=1.5)+
    #     theme_classic()+
    #     ylim(c(0,1))+
    #     theme(legend.position="top")+
    #     scale_color_manual(values=c("#DF8F44FF", "#374E55FF"))
    #   
    # }
    p
  })
  
  output$gsea_km_plot = renderPlot({
    gsea_km_plot()
  })
  
  output$gsea_table = renderDT({
    res = gsea_table()
    #   names(res)[2] = "ssGSEA_value"
    #  res
    DT:::datatable(
      res,
      rownames = FALSE,colnames=FALSE,
      options = list(scrollX=TRUE,ordering=F,dom="ft",pageLength=1000,paging= F,scrollCollapse=T,scrollY="50vh",
                     initComplete = JS(
                       "function(settings, json) {",
                       "$(this.api().table().header()).css({'background-color': '#4B4B4B', 'color': '#fff'});",
                       "}"),
                     headerCallback = JS(
                       "function(thead, data, start, end, display){",
                       "  $(thead).css({'background-color': '#4B4B4B', 'color': '#fff'});",
                       "}")
      ),
      container = tags$table(
        class="stripe row-border hover",
        tags$thead(tags$tr(lapply(colnames(res), tags$th)))
        
      )
    )
  })
  
  observe({
    output[["download_gsea_table"]]<-downloadHandler(
      filename = function() {
        paste0(input$cancer_type,"_joint_ssGSEA",'.csv', sep = '')
        
      },
      content = function(file){
        fwrite(gsea_table(),file,row.names = T)
      }
    )
    
  })
  
  observe({
    output[["download_gsea_KM"]]<-downloadHandler(
      filename = function() {
        paste0(input$cancer_type,"_joint_ssGSEA_KM",'.png', sep = '')
        
      },
      content = function(file){
        ggsave(file, gsea_km_plot(),type="cairo-png",width = 12, height = 12,units="cm")
      }
    )
  })
  
  #----------------------------------------------------------------------------------------------------
  ######################################### EGA subnetwork #################################################
  #----------------------------------------------------------------------------------------------------
  output$subnet_ui = renderUI({
    div(
      
      column(5,
             actionButton("search_subnet", "Click to Run", style = "width: 50%; height: 38px; background-color: #FFD700; color: black; border: none;")
      )
      
    )
  })
  
  subnetwork_res = eventReactive(input$search_subnet,{
    dt = data_multi()
    dt_sub = subset(dt,dt$study==input$study_type_multi)
    dt_sub_ega = dt_sub[,-c(1:5)]
    
    ega = EGA(dt_sub_ega,model = "glasso",plot.EGA = T)
    
  },ignoreNULL = F)
  
  output$subnetwork = renderPlot({
    subnetwork_res()
  })
  
  ###################### forest plot #########################
  forest_res = eventReactive(input$search_subnet,{
    dt = data_multi()
    
    dt_sub = subset(dt,dt$study==input$study_type_multi)
    dt_sub_ega = dt_sub[,-c(1:5)]
    
    res = subnetwork_res()$dim.variables
    res = na.omit(res)
    dt_sub_ega_t = data.frame(t(dt_sub_ega))
    res_m = merge(res,dt_sub_ega_t,by="row.names")
    # for(i in 4:ncol(res_m)){
    #   res_m[,i] = log2(res_m[,i]+1)
    # }
    res_mean = lapply(4:ncol(res_m),function(x){
      res =  res_m %>%
        group_by(dimension) %>% 
        summarise_at(colnames(res_m)[x], funs(mean(., na.rm=TRUE)))
      res[,-1]
    })
    
    res_final = do.call(cbind,res_mean)
    res_final_t = t(res_final)
    
    
    
    res_final_t = data.frame("OS.time"=dt_sub$OS.time,"OS"=dt_sub$OS,res_final_t)
    
    covariates = names(res_final_t)[-c(1:2)]
    univ_formulas <- sapply(covariates,
                            function(x) as.formula(paste('Surv(OS.time, OS)~', x)))
    univ_models <- lapply( univ_formulas, function(x){coxph(x, data = res_final_t)})
    univ_results <- lapply(univ_models,
                           function(x){ 
                             x <- summary(x)
                             p.value<-signif(x$wald["pvalue"], digits=2)
                             wald.test<-signif(x$wald["test"], digits=2)
                             beta<-signif(x$coef[1], digits=2);#coeficient beta
                             HR <-signif(x$coef[2], digits=2);#exp(beta)
                             HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                             HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                             res<-c("beta"=beta, "HR"=HR, "HR_lower"=HR.confint.lower, "HR_upper"=HR.confint.upper,"pval"=p.value)
                             
                             return(res)
                             #return(exp(cbind(coef(x),confint(x))))
                           })
    univ_results = data.frame(do.call(rbind,univ_results))
    if(nrow(univ_results)==1){
      univ_results$names=c("group1")
    }else{
      univ_results$names=paste0("group",rownames(univ_results))
    }
    univ_results$symbol = ifelse(univ_results$beta>=0,"HR>=1","HR<1")
    
    p = ggplot(univ_results, aes(y=names, x=HR, xmin=HR_lower, xmax=HR_upper,color = symbol))+
      #Add data points and color them black
      #  geom_point(shape=15,size = 3)+
      geom_pointrange(shape=15,size = 1.5,fatten =1.5)+
      #add the CI error bars
      geom_errorbarh(height=.2)+
      #Specify the limits of the x-axis and relabel it to something more meaningful
      #  scale_x_continuous(limits=c(-2,2), name='Standardized Mean Difference (d)')+
      #Give y-axis a meaningful label
      ylab('')+
      xlab("OS (HR)")+
      #Add a vertical dashed line indicating an effect size of zero, for reference
      geom_vline(xintercept=1, color='black', linetype='dotted',size=1)+
      theme_classic()+
      #scale_color_brewer(palette="Set1")+
      scale_color_manual(values=c("#374E55FF", "#DF8F44FF"))+
      theme(legend.position='none',axis.line = element_line(colour = "gray"),
            axis.text.y = element_text(face="bold", color="gray",size=15),
            axis.text.x = element_text(face="bold", color="gray",size=15),
            axis.title.x = element_text(size = 15),
            axis.ticks.x = element_blank()
      )
    p
  })
  
  output$subnetwork_forest = renderPlot(forest_res()) 
  
  #-------------------------------------------------------------------
  ############# tab1: correlations #############################
  #-------------------------------------------------------------------
  dt_pca = reactive({
    tryCatch({
      dt = data_multi()
      #dt = read_fst(paste0(folder_main,"www/pca/lnc_mrna_clin.fst"),columns = c("type","EGFR","HOTAIR","CD8A","IFNG.AS1"))
      
      # dt = read_fst("C:/Users/4467777/Desktop/TCGA_shiny/TCGA_prognositic_app/www/pca/mRNA_clin.fst",columns = c("type",input$genePCA))
      
      dt = dt[,-c(1:5)]
      
      
      dt
    },error = function(e){
      ""
    })
    
  })
  ########################## PCA plot ############################################
  pca_plot = reactive({
    tryCatch({
      res.pca <- prcomp(as.matrix(na.omit(dt_pca())), scale = TRUE)
      p = fviz_pca_var(res.pca,arrowsize = 2, pointsize = 1.5,circlesize = 0.8,
                       col.var = "contrib", # Color by contributions to the PC
                       gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                       repel = TRUE     # Avoid text overlapping
      )+
        ggtitle("")
      p        
    },error = function(e){
      ""
    })
    
  })
  output$pca_plot = renderPlot({
    pca_plot()
  })
  
  ############################ corrplot ##########################################
  #corr_plot= eventReactive(input$search_pca,{
  corr_plot = reactive({
    tryCatch({
      corr = round(cor(dt_pca(),method = "spearman",use = 'pairwise.complete.obs'),3)
      
      if(nrow(corr)>15){
        label=FALSE
      }else{
        label=TRUE
      }
      #p.mat <- cor_pmat(dt)
      p = ggcorrplot(
        corr, hc.order = TRUE, type = "lower", outline.col = "white",
        ggtheme = ggplot2::theme_classic(),lab=label,
        colors = c("#6D9EC1", "white", "#E46726")
      )+
        theme(axis.line = element_line(colour = "gray"),
              axis.ticks.x = element_blank()
              
        )
      p 
    },error = function(e){
      ""
    })
    
  })
  
  
  output$corr_plot = renderPlot({
    corr_plot()
  })
  
  ############################## Network plot ####################################
  # net_data = eventReactive(input$search_pca,{
  net_data = reactive({
    tryCatch({
      adjm = cor(dt_pca(),method = "spearman",use = 'pairwise.complete.obs')
      #validate(input$search_net)
      adjm[abs(adjm)<input$p_network] <- 0
      
      adjm        
    },error=function(e){
      ""
    })
    
  })
  
  output$corr_network_plot = renderPlot({
    # tryCatch({
    network<- graph_from_adjacency_matrix(net_data(), weighted=T, mode="undirected", diag=F)
    
    par(mar=c(0,0,0,0))
    plot(network,
         vertex.size=28,# Size of the node (default is 15)
         vertex.label.family="Helvetica", #or "Times"
         vertex.label.cex=0.9,
         # vertex.color=my_color,
         edge.curved=0.9)
    # },error=function(e){
    #   ""
    # })
    
    
  },bg="grey13")  
  
############################ Output page #####################################
  data_ICI_info = reactive({
    dt = fread(paste0(folder_main,"www/data_ICI_info.csv"))
  })
  
output$data_ICI_info = renderDT({
  res = data_ICI_info()
  DT:::datatable(
    res,
    rownames = FALSE,colnames=FALSE,
    options = list(scrollX=TRUE,ordering=F,dom="ft",pageLength=1000,paging= F,scrollCollapse=T,scrollY="80vh",
                   initComplete = JS(
                     "function(settings, json) {",
                     "$(this.api().table().header()).css({'background-color': '#4B4B4B', 'color': '#fff'});",
                     "}"),
                   headerCallback = JS(
                     "function(thead, data, start, end, display){",
                     "  $(thead).css({'background-color': '#4B4B4B', 'color': '#fff'});",
                     "}")
    ),
    container = tags$table(
      class="stripe row-border hover",
      tags$thead(tags$tr(lapply(colnames(res), tags$th)))
      
    )
  )
  
})
}

shinyApp(ui, server)
