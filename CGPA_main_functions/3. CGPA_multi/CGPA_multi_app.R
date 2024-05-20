################################
# multi-gene panel app
###############################
packages = c("shiny","rstudioapi","data.table","tidyverse","bslib","shinyWidgets","shinydashboard","shinydashboardPlus","shinyjs","shinyBS","DT",
             "shinycustomloader","fst","ggcorrplot","igraph","factoextra","gbm","survival","survminer","AdjKMCIF","GSVA","EGAnet") # rstudioapi not neededd

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

folder_main = dirname(getSourceEditorContext()$path)
source(paste0(folder_main,"functions.R"))
genesets_pca = fread("www/genesets_pca.csv")
options(shiny.maxRequestSize = 500*1024^2)  # 500MB in bytes


# light = bs_theme(version = 4, bootswatch = "cosmo", bg="#E3E3E3",fg = "black",warning = "#FFD300") %>%
#     bs_add_rules(sass::sass_file("www/style/style_light.scss"))
dark = bs_theme(version = 3,bg="black",fg = "#fff")%>%
    bs_add_rules(sass::sass_file("www/style/style.scss"))

ui= fluidPage(
  useShinyjs(),
    theme = dark,
    tags$style("
              body {
    -moz-transform: scale(1, 1); /* Moz-browsers */
    zoom: 0.8; /* Other non-webkit browsers */
    zoom: 80%; /* Webkit browsers */
              }

              "),
    useShinydashboard(),
    
    div(class="navbar1",
        navbarPage(title = "CGPA PLUS",
                             
                             fluid = TRUE, 
                             collapsible = TRUE,
                         
                         
                            tabPanel("Multi-gene Panel Discovery",
                                     br(),
                             sidebarLayout(
                               sidebarPanel(width=2,
                                            br(),br(),
                                            selectizeInput("cancer_type","Cancer type",c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA",
                                                                                         "GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC",
                                                                                         "LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ",
                                                                                         "SARC","SKCM","STAD","TGCT","THCA","UCEC","UCS","UVM"),multiple=FALSE,
                                                           selected=c("HNSC")),
                                            selectizeInput("geneset",label="Select a pre-defined geneset",choices = c("customize your own gene"="None",genesets_pca$genesets),selected  = "None",multiple = FALSE),
                                            selectizeInput("gene", label="Or customize gene list", 
                                                           choices =NULL,multiple=TRUE),

                                            wellPanel(
                                              useShinyjs(),
                                              div(
                                              fileInput('file', 
                                                        span(
                                                          list(HTML("<p><abbr title='The uploaded data has to be one column csv file with geneNames on each row'>Or, upload your genelist (csv file)...</abbr></p>"))
                                                        ),
                                                        accept = c(
                                                          '.csv'
                                                        )),
                                              style="font-size:80%;"
                                              ),
                                              actionButton('reset', 'Reset the file'),#Clear input dataset#
                                              tags$head(
                                                tags$style(HTML('#reset{padding:8px;font-size:80%}'))
                                              )
                                              
                                            ),
                                            # actionButton("search_pca", "Update"),
                                            # tags$head(
                                            #     tags$style(HTML('#search_pca{border-color: #b58900;} #search_pca:hover{background-color:#b58900;} '))
                                            # )
                               ),
                               mainPanel(width=10,
                                         tabsetPanel(type="pills",id="mainnav",

                                           tabPanel("PROGNOSTIC RANKING",
                                                    br(),
                                                    
                                                   uiOutput("gbm_surv_ui"),
                                                      
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
                                                                                    div(withLoader(plotOutput("waterfall_plot",width="80%"),type="html",loader="loader1"),align="center")
                                                                                    # 
                                                                                    # conditionalPanel(
                                                                                    #   condition="input.search_gbm>0",
                                                                                    #   div(withLoader(plotOutput("waterfall_plot",width="80%"),type="html",loader="loader1"),align="center")
                                                                                    #   #DTOutput("test_path")
                                                                                    #   
                                                                                    # ),
                                                                                    # conditionalPanel(
                                                                                    #   condition="input.search_gbm==0",
                                                                                    #   br(),
                                                                                    #   p("Click 'Run' button to get the figure")
                                                                                    # )
                                                                                    
                                                                                    
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
                                                                                    # 
                                                                                    # conditionalPanel(
                                                                                    #   condition="input.search_gbm>0",
                                                                                    #   div(withLoader(plotOutput("gbm_plot",width="80%"),type="html",loader="loader1"),align="center")
                                                                                    #   #DTOutput("test_path")
                                                                                    #   
                                                                                    # ),
                                                                                    # conditionalPanel(
                                                                                    #   condition="input.search_gbm==0",
                                                                                    #   br(),
                                                                                    #   p("Click 'Run' button to get the figure")
                                                                                    # )
                                                                                    
                                                                                    
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
                                                      column(2,
                                                             div(style = "display: flex; align-items:  flex-start;",  # Ensure items are vertically centered
                                                                 tags$div("Include anchor genes?", style = "margin-right: 5px; font-size: 16px; white-space: nowrap; line-height: 1.5;margin-top:5px;"),  # Keep label styling simple
                                                                 selectizeInput("anchor_not", NULL, c("No", "Yes"), "No", options = list(width = '100%'))
                                                             )
                                                            # radioButtons("anchor_not","Include anchor genes?",choices=c("No","Yes"),"No"),
                                                            
                                                      ),
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
                                         ))
                                )
                            ),
                   tabPanel("Custom Portal",
                            br(),
                            
                            
                            sidebarLayout(
                              sidebarPanel(width=2,
                                           br(),br(),
                                           wellPanel(
                                             div(fileInput('file_omics',
                                                           span(
                                                             list(HTML("<p><abbr title='The uploaded data has to be a csv file with the 1st column as geneNames, and the remaining columns are samples'>Upload your omic file (csv file)...</abbr></p>"))
                                                           ),
                                                           accept = c(
                                                             '.csv'
                                                           )),
                                                 style="font-size:100%;"
                                             ),
                                             div(fileInput('file_pheno',
                                                           span(
                                                             list(HTML("<p><abbr title='Please upload your phenotype file in CSV format. The first column should contain sample identifiers, followed by columns with clinical information such as OS.time, OS, PFI.time, PFI, among others'>Upload your phenotype file (csv file)....</abbr></p>"))
                                                           ),
                                                           accept = c(
                                                             '.csv'
                                                           )),
                                                 style="font-size:100%;"
                                             ),
                                             actionButton("load_omics", "View omics Example"),
                                             bsModal("omics_example_window", "", "load_omics", size = "large",
                                                     withLoader(
                                                       
                                                       DTOutput("omics_example"),type="html",loader="loader1")),
                                             
                                             div(style = "padding-top: 20px;"),
                                             actionButton("load_pheno", "View phenotype Example"),
                                             bsModal("pheno_example_window", "", "load_pheno", size = "large",
                                                     withLoader(DTOutput("pheno_example"),type="html",loader="loader1"))
                                           ),
                                           
                                           wellPanel(
                                             useShinyjs(),
                                             style = "background:gray",
                                             selectizeInput("gene_custom", label="Input the interested genes", 
                                                            choices =NULL,multiple=TRUE),
                                             
                                             tags$head(
                                               tags$style(HTML('#reset{padding:8px;font-size:80%}'))
                                             )
                                             
                                           )
                                           
                              ),
                              mainPanel(width=10,
                                        tabsetPanel(type="pills",id="mainnav",
                                                    tabPanel("VIEW DATA",
                                                             br(),
                                                             div(id="box1",
                                                                 box(
                                                                   title = "View the top uploaded dataset",collapsible = TRUE,status="warning",
                                                                   solidHeader = TRUE,
                                                                   width = 12,height=800,
                                                                   div(style='max-width: 100%; height: 500px; width: auto;overflow-y: scroll;',
                                                                       withLoader(DTOutput("upload_omics_view"),type="html",loader="loader1")
                                                                       
                                                                   )
                                                                 )  
                                                             )
                                                    ),
                                                   
                                                    tabPanel("PROGNOSTIC RANKING",
                                                             fluidRow(column(9,
                                                                             br(),
                                                                             fluidRow(
                                                                             uiOutput("gbm_surv_ui_custom"))
                                                                             ,
                                                                             br(),
                                                                             fluidRow(
                                                                               column(6,div(id="box2",
                                                                                            box(
                                                                                              title="Waterfall plot for univariate Cox",collapsible = TRUE, status="warning",
                                                                                              solidHeader = TRUE,
                                                                                              dropdown(
                                                                                                radioButtons("waterfall_menu_custom","Exclude non-significant markers?",c("Yes","No"),"No"),
                                                                                                # actionButton("search_net", "Run"),
                                                                                                circle = TRUE, status = "warning", icon = icon("gear"), width = "300px",
                                                                                                tooltip = tooltipOptions(title = "Click to see inputs !")
                                                                                              ),
                                                                                              width = 12,height=500,
                                                                                              div(withLoader(plotOutput("waterfall_plot_custom",width="80%"),type="html",loader="loader1"),align="center")
                                                                                              
                                                                                            )
                                                                               )),
                                                                               column(6,div(id="box2",
                                                                                            box(
                                                                                              title="Gradient boosting importance score",collapsible = TRUE, status="warning",
                                                                                              solidHeader = TRUE,
                                                                                              dropdown(
                                                                                                numericInput("gb_menu_custom","Input a threshold to exclude low influence markers",0),
                                                                                                # actionButton("search_net", "Run"),
                                                                                                circle = TRUE, status = "warning", icon = icon("gear"), width = "300px",
                                                                                                tooltip = tooltipOptions(title = "Click to see inputs !")
                                                                                              ),
                                                                                              tags$style(HTML('#sw-content-dropdown, .sw-dropdown-in {background-color: gray;}')),
                                                                                              width = 12,height=500,
                                                                                              div(withLoader(plotOutput("gbm_plot_custom",width="80%"),type="html",loader="loader1"),align="center")
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
                                                                                       #hidden(div(id = "dataUI", uiOutput("average_km_ui_custom"))),
                                                                                         uiOutput("average_km_ui_custom"),
                                                                                           
                                                                                           fluidRow(
                                                                                             column(9,br(),div(id="box2",
                                                                                                               box(
                                                                                                                 title="Average score KM plot",collapsible = TRUE, status="warning",
                                                                                                                 solidHeader = TRUE,
                                                                                                                 width = 12,height=500,
                                                                                                                 dropdown(
                                                                                                                   downloadButton("download_average_km_custom", "Download"),
                                                                                                                   # actionButton("search_net", "Run"),
                                                                                                                   circle = TRUE, status = "warning", icon = icon("gear"), width = "300px",
                                                                                                                   tooltip = tooltipOptions(title = "Click to download")
                                                                                                                 ),
                                                                                                                 
                                                                                                                 conditionalPanel(
                                                                                                                   condition="input.search_average_joint_custom>0",
                                                                                                                   div(withLoader(plotOutput("average_km_plot_custom",width="60%"),type="html",loader="loader1"),align="center"),
                                                                                                                   DTOutput("test_path")
                                                                                                                 ),
                                                                                                                 conditionalPanel(
                                                                                                                   condition="input.search_average_joint_custom==0",
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
                                                                                                                   downloadButton("download_average_table_custom", "Download"),
                                                                                                                   # actionButton("search_net", "Run"),
                                                                                                                   circle = TRUE, status = "warning", icon = icon("gear"), width = "300px",
                                                                                                                   tooltip = tooltipOptions(title = "Click to download")
                                                                                                                 ),
                                                                                                                 
                                                                                                                 conditionalPanel(
                                                                                                                   condition="input.search_average_joint_custom>0",
                                                                                                                   div(withLoader(DTOutput("average_table_custom_output"),type="html",loader="loader1"),align="center"),
                                                                                                                 ),
                                                                                                                 conditionalPanel(
                                                                                                                   condition="input.search_average_joint_custom==0",
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
                                                                                           uiOutput("gsea_km_ui_custom"),
                                                                                           fluidRow(
                                                                                             column(
                                                                                               9,br(),div(id="box2",
                                                                                                          box(
                                                                                                            title="ssGSEA KM plot",collapsible = TRUE, status="warning",
                                                                                                            solidHeader = TRUE,
                                                                                                            width = 12,height=500,
                                                                                                            dropdown(
                                                                                                              downloadButton("download_gsea_km_custom", "Download"),
                                                                                                              # actionButton("search_net", "Run"),
                                                                                                              circle = TRUE, status = "warning", icon = icon("gear"), width = "300px",
                                                                                                              tooltip = tooltipOptions(title = "Click to download")
                                                                                                            ),
                                                                                                            
                                                                                                            conditionalPanel(
                                                                                                              condition="input.search_gsea_custom>0",
                                                                                                              div(withLoader(plotOutput("gsea_km_plot_custom",width="60%"),type="html",loader="loader1"),align="center")
                                                                                                              #DTOutput("test_path")
                                                                                                              
                                                                                                            ),
                                                                                                            conditionalPanel(
                                                                                                              condition="input.search_gsea_custom==0",
                                                                                                              br(),
                                                                                                              p("Click 'Run' button to get the figure")
                                                                                                            )
                                                                                                            
                                                                                                            
                                                                                                          )
                                                                                               )),
                                                                                             column(
                                                                                               9,br(),div(id="box2",
                                                                                                          box(
                                                                                                            title="ssGSEA table",collapsible = TRUE, status="warning",
                                                                                                            solidHeader = TRUE,
                                                                                                            width = 12,height=500,
                                                                                                            dropdown(
                                                                                                              downloadButton("download_gsea_table_custom", "Download"),
                                                                                                              # actionButton("search_net", "Run"),
                                                                                                              circle = TRUE, status = "warning", icon = icon("gear"), width = "300px",
                                                                                                              tooltip = tooltipOptions(title = "Click to download")
                                                                                                            ),
                                                                                                            
                                                                                                            conditionalPanel(
                                                                                                              condition="input.search_gsea_custom>0",
                                                                                                              div(withLoader(DTOutput("gsea_table_custom",width="60%"),type="html",loader="loader1"),align="center")
                                                                                                              #DTOutput("test_path")
                                                                                                              
                                                                                                            ),
                                                                                                            conditionalPanel(
                                                                                                              condition="input.search_gsea_custom==0",
                                                                                                              br(),
                                                                                                              p("Click 'Run' button to get the figure")
                                                                                                            )
                                                                                                            
                                                                                                            
                                                                                                          )
                                                                                               ))
                                                                                             # column(3,br(),
                                                                                             #        uiOutput("gsea_km_ui_custom")
                                                                                             # )
                                                                                             
                                                                                           )
                                                                                           
                                                                                  )
                                                                      ),
                                                                      
                                                               )
                                                               
                                                             )
                                                    ),
                                                    tabPanel("SUBNETWORK",
                                                             br(),
                                                             fluidRow(
                                                               uiOutput("subnet_ui_custom")
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
                                                                                       condition = "input.search_subnet_custom>0",
                                                                                       div(style='overflow-x: scroll;overflow-y: scroll;',
                                                                                           withLoader(plotOutput("subnetwork_custom",width="80%",height="600px"),type="html",loader="loader1"),align="center")
                                                                                       #DTOutput("test"))
                                                                                     ),
                                                                                     conditionalPanel(
                                                                                       condition="input.search_subnet_custom==0",
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
                                                                                       condition = "input.search_subnet_custom>0",
                                                                                       div(style='overflow-x: scroll;overflow-y: scroll;',
                                                                                           withLoader(plotOutput("subnetwork_forest_custom",width="30%",height="400px"),type="html",loader="loader1"),align="center")
                                                                                       
                                                                                     ),
                                                                                     conditionalPanel(
                                                                                       condition="input.search_subnet_custom==0",
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
                                                                                numericInput("p_network_custom","Choose a minimum correlation threshold (click 'update' once complete)",value=0.3,min=0,max=0.9),
                                                                                # actionButton("search_net", "Update"),
                                                                                circle = TRUE, status = "warning", icon = icon("gear"), width = "300px",
                                                                                tooltip = tooltipOptions(title = "Click to see inputs !")
                                                                              ),
                                                                              tags$style(HTML('#sw-content-dropdown, .sw-dropdown-in {background-color: gray;}')),
                                                                              div(withLoader(plotOutput("corr_network_plot_custom",width="80%",height="500px"),type="html",loader="loader1"),align="center")
                                                                            )
                                                               )),
                                                               column(3),
                                                               column(1),
                                                               column(8,div(id="box1",
                                                                            box(
                                                                              title="PCA plot",collapsible = TRUE, status="warning",
                                                                              solidHeader = TRUE,
                                                                              width = 12,height=800,
                                                                              div(withLoader(plotOutput("pca_plot_custom",width="80%",height="500px"),type="html",loader="loader1"),align="center")
                                                                              
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
                                                                              div(withLoader(plotOutput("corr_plot_custom",width="80%",height="500px"),type="html",loader="loader1"),align="center")
                                                                            )
                                                               )),
                                                               column(3)
                                                               
                                                             )
                                                    )
                                        ))
                            )
                   ),
                   tabPanel("Help",
                            mainPanel(
                              column(2),
                              column(10,
                                     br(),br(),
                                     uiOutput("tab_instructions")
                              )
                              
                            )
                   )
                            
                             # inputs =div(
                             #   class = "custom-control custom-switch", 
                             #   tags$input(
                             #     id = "light_mode", type = "checkbox", class = "custom-control-input",
                             #     onclick = HTML("Shiny.setInputValue('light_mode', document.getElementById('dark_mode').value);")
                             #   ),
                             #   tags$label(
                             #     "Light mode", `for` = "light_mode", class = "custom-control-label"
                             #   )
                             # )
                             )
        )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    
  query = reactive({
    query = parseQueryString(session$clientData$url_search)
    query
  })
  
  # observe({
  #   query=query()
  #   if (!is.null(query[['geneid']])) {
  #     updateSelectizeInput(session,"geneset",choices = c("None",genesets_pca$genesets),selected  = "None")
  #   }else{
  #     updateSelectizeInput(session,"geneset",choices = c("None",genesets_pca$genesets),selected  = genesets_pca$genesets[1])
  #     
  #   }
  # })
  
  observe({
    query=query()
    #updateSelectizeInput(session,"geneset",choices = c("None",genesets_pca$genesets),selected  = "None")
    
    if (!is.null(query[['geneid']])) {
      if(input$geneset=="None"){
        genes = strsplit(query[['geneid']],",")[[1]]
        updateSelectizeInput(session,"gene",choices = geneName_mrna()[,1],selected  = genes,server = T)
      }else{
        genes_sel = genesets_pca$gene[which(genesets_pca$genesets==input$geneset)]
        genes = strsplit(genes_sel,",")[[1]]
        updateSelectizeInput(session,"gene",choices = geneName_mrna()[,1],selected  = genes,server = T)
      }
    }else{
      if(input$geneset=="None"){
        updateSelectizeInput(session,"gene",choices = geneName_mrna()[,1],server = T)
      }else{
        # [which(genesets()$genesets==input$genesetPCA)]
        
        genes_sel = genesets_pca$gene[which(genesets_pca$genesets==input$geneset)]
        genes = strsplit(genes_sel,",")[[1]]
        updateSelectizeInput(session,"gene",choices = geneName_mrna()[,1],selected  = genes,server = T)
        
      }
    }
  })
  
  
  # observe({
  #     session$setCurrentTheme(
  #         if (isTRUE(input$light_mode)) light else dark
  #     )
  # })
  
  # geneName = reactive({
  #   read_fst(paste0(folder_main,"www/geneName.fst"))
  # })
  geneName_mrna = reactive({
    read_fst("www/geneName_mrna.fst")
  })
  # geneName_lnc = reactive({
  #   read_fst(paste0(folder_main,"www/geneName_lnc.fst"))
  # })
  
  ################################ common UI settings #######################################################  
  output$gbm_surv_ui = renderUI({
    # div(
    #   selectizeInput("gbm_surv_type","OS or PFI",c("OS","PFI"),"OS")
    # )
    # selectInput("gbm_surv_type", "OS or PFI", 
    #             choices = c("OS", "PFI"),
    #             style = "align-self: center; margin-right: 10px;")
    # 
    
    div(style = "display: flex; align-items: flex-start;",  # Adjust alignment to start at the top
        tags$div("OS or PFI ", style = "margin-right: 5px; margin-top: 5px; font-size: 16px;"),  # Adjust "Input 2" label position
        selectizeInput("gbm_surv_type",NULL,c("OS","PFI"),"OS")
    )
    
  })
  output$average_km_ui = renderUI({
    div(
      fluidRow(
        column(2,
               div(style = "display: flex; align-items:  flex-start;",  # Ensure items are vertically centered
                   tags$div("OS or PFI", style = "margin-right: 5px; font-size: 16px; white-space: nowrap; line-height: 1.5;margin-top:5px;"),  # Keep label styling simple
                   selectizeInput("gsea_km_type_av", NULL, c("OS", "PFI"), "OS", options = list(width = '100%'))
               )
        ),
        column(2,
               div(style = "display: flex; align-items: flex-start;",  # Ensure items are vertically centered
                   tags$div("KM cutoff",style = "margin-right: 5px; font-size: 16px; white-space: nowrap; line-height: 1;margin-top:5px;"),  # Keep label styling simple
                   selectizeInput("gsea_km_cutoff_av", NULL, c("median", "optimal", "quartile"), "median", options = list(width = '100%'))
               )
        ),
        column(3,
               div(style = "display: flex; align-items: flex-start;",  # Ensure items are vertically centered
                   tags$div("Adjust for covariates?", style = "margin-right: 5px; font-size: 16px; white-space: nowrap; line-height: 1.5;margin-top:5px;"),  # Keep label styling simple
                   selectizeInput("gsea_cov_av", NULL, c("None", "Age + Sex", "Age + Sex + tumor_Purity"), "None", options = list(width = '100%'))
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
               div(style = "display: flex; align-items:  flex-start;",  # Ensure items are vertically centered
                   tags$div("OS or PFI", style = "margin-right: 5px; font-size: 16px; white-space: nowrap; line-height: 1.5;margin-top:5px;"),  # Keep label styling simple
                   selectizeInput("gsea_km_type", NULL, c("OS", "PFI"), "OS", options = list(width = '100%'))
               )
        ),
        column(2,
               div(style = "display: flex; align-items: flex-start;",  # Ensure items are vertically centered
                   tags$div("KM cutoff",style = "margin-right: 5px; font-size: 16px; white-space: nowrap; line-height: 1;margin-top:5px;"),  # Keep label styling simple
                   selectizeInput("gsea_km_cutoff", NULL, c("median", "optimal", "quartile"), "median", options = list(width = '100%'))
               )
        ),
        column(3,
               div(style = "display: flex; align-items: flex-start;",  # Ensure items are vertically centered
                   tags$div("Adjust for covariates?", style = "margin-right: 5px; font-size: 16px; white-space: nowrap; line-height: 1.5;margin-top:5px;"),  # Keep label styling simple
                   selectizeInput("gsea_cov", NULL, c("None", "Age + Sex", "Age + Sex + tumor_Purity"), "None", options = list(width = '100%'))
               )
        ),
        column(2,  # Keep the button aligned
               actionButton("search_gsea", "Click to Run", style = "width: 50%; height: 38px; background-color: #FFD700; color: black; border: none;")
        )
      )
      
      
      
      
      
    )
    

  })
  
  output$subnet_ui = renderUI({
    div(
      fluidRow(
        column(5,
               conditionalPanel(
                 condition="input.anchor_not=='Yes'",
                 div(style = "display: flex; align-items: flex-start;",  # Ensure items are vertically centered
                     tags$div("Choose anchor genes", style = "margin-right: 5px; font-size: 16px; white-space: nowrap; line-height: 1.5;margin-top:5px;"),  # Keep label styling simple
                     selectizeInput("anchor_gene", NULL, c("PTPRC","CD3E","ACTA2","DCN",
                                                           "KRT17","KRT14"), c("PTPRC","KRT14","DCN"), options = list(width = '100%'))
                 )

                 
               )
               ),
        column(5,
               actionButton("search_subnet", "Click to Run", style = "width: 50%; height: 38px; background-color: #FFD700; color: black; border: none;")
               )
      )

     
    )
  })
  # observe({
  #   if(input$geneset=="None"){
  #     updateSelectizeInput(session,"gene",choices = geneName_mrna()[,1],server = T)
  #   }else{
  #     # [which(genesets()$genesets==input$genesetPCA)]
  # 
  #     genes_sel = genesets_pca$gene[which(genesets_pca$genesets==input$geneset)]
  #     genes = strsplit(genes_sel,",")[[1]]
  #     updateSelectizeInput(session,"gene",choices = geneName_mrna()[,1],selected  = genes,server = T)
  # 
  #   }
  # })
  # 
  rv = reactiveValues(data=NULL)
  observe({
    if(is.null(input$file)){
      if(input$anchor_not=="No"){
        rv$data = read_fst("www/pca/mRNA_clin.fst",columns = c("Row.names","purity","type","age","gender","type","OS","OS.time","PFI","PFI.time",make.names(input$gene)))
      }else if(input$anchor_not=="Yes"){
        rv$data = read_fst("www/pca/mRNA_clin.fst",columns = c("Row.names","purity","type","age","gender","type","OS","OS.time","PFI","PFI.time",make.names(input$anchor_gene),make.names(input$gene)))
        
      }
    }else{
      req(input$file)
      infile <- input$file
      if (is.null(infile)){
        return(NULL)
      }
      genes = read.csv(input$file$datapath)
      
      sel_genes = unique(c(make.names(genes[,1]),make.names(input$gene)))
      
      if(input$anchor_not=="No"){
        rv$data = read_fst("www/pca/mRNA_clin.fst",columns = c("Row.names","purity","type","age","gender","type","OS","OS.time","PFI","PFI.time",sel_genes))
      }else if(input$anchor_not=="Yes"){
        rv$data = read_fst("www/pca/mRNA_clin.fst",columns = c("Row.names","purity","type","age","gender","type","OS","OS.time","PFI","PFI.time",make.names(input$anchor_gene),sel_genes))
        
      }
    }
    
  })
  
  observeEvent(input$reset,{
    sel_genes = unique(make.names(input$gene))
    
    if(input$anchor_not=="No"){
      rv$data = read_fst("www/pca/mRNA_clin.fst",columns = c("Row.names","purity","type","age","gender","type","OS","OS.time","PFI","PFI.time",sel_genes))
    }else if(input$anchor_not=="Yes"){
      rv$data = read_fst("www/pca/mRNA_clin.fst",columns = c("Row.names","purity","type","age","gender","type","OS","OS.time","PFI","PFI.time",make.names(input$anchor_gene),sel_genes))
      
    }
    reset('file')
  })
  
  ############################################################################################################
  #-------------------------------------------------------------------
  ############# tab1: correlations #############################
  #-------------------------------------------------------------------
  dt_pca = reactive({
    tryCatch({
      dt = rv$data
      #dt = read_fst(paste0(folder_main,"www/pca/lnc_mrna_clin.fst"),columns = c("type","EGFR","HOTAIR","CD8A","IFNG.AS1"))
      
      # dt = read_fst("C:/Users/4467777/Desktop/TCGA_shiny/TCGA_prognositic_app/www/pca/mRNA_clin.fst",columns = c("type",input$genePCA))
      dt = subset(dt,dt$type==input$cancer_type)
      
      if(input$anchor_not=="No"){
        dt = dt[,-c(1:10)]
      }else if(input$anchor_not=="Yes"){
        length_anchor = length(input$anchor_gene)
        dt = dt[,-c(1:(10+length_anchor))]
      }
      
      dt = data.frame(t(dt))
      
      geneName = geneName_mrna()
      geneName$make.names = make.names(geneName$geneName)
      dt = merge(geneName,dt,by.x = "make.names",by.y = "row.names")
      dt = dt[!duplicated(dt$geneName),]
      
      rownames(dt) = dt$geneName
      dt = t(dt[,-c(1:2)])
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
    tryCatch({
      network<- graph_from_adjacency_matrix(net_data(), weighted=T, mode="undirected", diag=F)
      
      par(mar=c(0,0,0,0))
      plot(network,
           vertex.size=28,# Size of the node (default is 15)
           vertex.label.family="Helvetica", #or "Times"
           vertex.label.cex=0.9,
           # vertex.color=my_color,
           edge.curved=0.9)
    },error=function(e){
      ""
    })
    
    
  },bg="grey13")
  
  #---------------------------------------------------------------------------------------
  ###################### Grdient boosting ############################################
  #---------------------------------------------------------------------------------------
  
  
  gbm_figure = reactive({
    
    dt =   rv$data 
    dt = na.omit(dt)
    
    if(is.null(input$file)){
      #sel_genes =  make.names(input$gene)
      if(input$anchor_not=="No"){
        sel_genes = make.names(names(dt)[c(11:ncol(dt))])
      }else{
        length_anchor = length(input$anchor_gene)
        sel_genes = make.names(names(dt)[c((10+length_anchor+1):ncol(dt))])
        
      }
    }else{
      if(input$anchor_not=="No"){
        sel_genes = make.names(names(dt)[c(11:ncol(dt))])
      }else{
        length_anchor = length(input$anchor_gene)
        sel_genes = make.names(names(dt)[c((10+length_anchor+1):ncol(dt))])
        
      }
    }
    
    dt = subset(dt,dt$type==input$cancer_type)
    if(input$gbm_surv_type=="OS"){
      fit_formula = as.formula(paste("Surv(OS.time,OS)~", paste0(paste(sel_genes, collapse=" + "))))
    }else if(input$gbm_surv_type=="PFI"){
      fit_formula = as.formula(paste("Surv(PFI.time,PFI)~", paste0(paste(sel_genes, collapse=" + "))))
    }
    gbm1 = gbm(fit_formula,data = dt,distribution="coxph")
    best.iter <- gbm.perf(gbm1,method="OOB") # returns test set estimate of best number of trees
    #summary(gbm1,n.trees=best.iter,cBars = length(input$gene_gbm))
    dt_p = summary(gbm1)
    dt_p = dt_p[order(dt_p$rel.inf,decreasing=T),]
    
    geneName = geneName_mrna()
    geneName$make.names = make.names(geneName$geneName)
    dt_p = merge(geneName,dt_p,by.x="make.names",by.y="row.names")
    dt_p = dt_p[!duplicated(dt_p$geneName),]
    rownames(dt_p) = dt_p$geneName
    dt_p = dt_p[,-c(1:2)]
    dt_p$var=rownames(dt_p)
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
    dt =   rv$data 
    #sel_genes = make.names(input$gene_gbm)
    #dt = read_fst(paste0(folder_main,"www/pca/lnc_mrna_clin.fst"),columns = c("Row.names","type","OS","OS.time","PFI","PFI.time",sel_genes))
    
    if(is.null(input$file)){
      #sel_genes =  make.names(input$gene)
      if(input$anchor_not=="No"){
        sel_genes = make.names(names(dt)[c(11:ncol(dt))])
      }else{
        length_anchor = length(input$anchor_gene)
        sel_genes = make.names(names(dt)[c((10+length_anchor+1):ncol(dt))])
        
      }
    }else{
      if(input$anchor_not=="No"){
        sel_genes = make.names(names(dt)[c(11:ncol(dt))])
      }else{
        length_anchor = length(input$anchor_gene)
        sel_genes = make.names(names(dt)[c((10+length_anchor+1):ncol(dt))])
        
      }
    }
    
    dt = na.omit(dt)
    dt = subset(dt,dt$type==input$cancer_type)
    for(i in 7:ncol(dt)){
      dt[,i] = log2(dt[,i]+1)
    }
    if(input$gbm_surv_type=="OS"){
      #fit_formula = as.formula(paste("Surv(OS.time,OS)~", paste0(paste(sel_genes, collapse=" + "))))
      res = lapply(1:length(sel_genes),function(x){
        fit_formula = as.formula(paste("Surv(OS.time,OS)~", sel_genes[x]))
        fit = coxph(fit_formula,data=dt)
        res = data.frame(summary(fit)$coef)
        names(res) = c("coef","HR","se","z","p.value")
        res$gene = paste0(rownames(res),res$mark)
        res
      })
    }else if(input$gbm_surv_type=="PFI"){
      # fit_formula = as.formula(paste("Surv(PFI.time,PFI)~", paste0(paste(sel_genes, collapse=" + "))))
      res = lapply(1:length(sel_genes),function(x){
        fit_formula = as.formula(paste("Surv(PFI.time,PFI)~", sel_genes[x]))
        fit = coxph(fit_formula,data=dt)
        res = data.frame(summary(fit)$coef)
        names(res) = c("coef","HR","se","z","p.value")
        res$gene = paste0(rownames(res),res$mark)
        res
      })
    }
    # fit = coxph(fit_formula,data=dt)
    # 
    # res = data.frame(summary(fit)$coef)
    # names(res) = c("coef","HR","se","z","p.value")
    # res$mark = ifelse(res$p.value<0.05,"**","")
    # 
    # res$gene = paste0(rownames(res),res$mark)
    # res = res[order(res$HR,decreasing = T),]
    # range_res = range(res$HR)
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
  
  #===========================================================================================
  ######################## SSGSEA ##################################################
  #===========================================================================================
  gsea_table = eventReactive(input$search_gsea,{
    validate(need(input$search_gsea !="", "Please wait while calculations are running....."))
    
    progress = shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message="Estimating ssGSEA scores for the gene sets",
                 detail = "This will take ~1 min")
    
    progress$inc(1/5, message=paste("Import data"))
    
    dt = read_fst("www/pca/mRNA_clin.fst")
    
    dt_sub = subset(dt,dt$type==input$cancer_type)
    dt_sub = data.frame(dt_sub[,c(1,12:ncol(dt_sub))])
    rownames(dt_sub) = dt_sub$Row.names
    dt_sub = data.frame(dt_sub[,-1])
    dt_sub = data.frame(t(dt_sub))
    #  list(c(input$gene_gsea))
    GE_matrix <- as.matrix(dt_sub)
    
    gene_gsea = rv$data
    if(input$anchor_not=="No"){
      gene_gsea = gene_gsea[,-c(1:10)]
    }else if(input$anchor_not=="Yes"){
      length_anchor = length(input$anchor_gene)
      gene_gsea = gene_gsea[,-c(1:(10+length_anchor))]
    }
    
    progress$inc(2/5, message=paste("Estimating ssGSEA scores for the gene sets"))
    gsva_H <- gsva(expr= GE_matrix, list(names(gene_gsea)), method="ssgsea")
    progress$inc(1/5, message=paste("Combining with clinical info"))
    gsva_H = data.frame(t(gsva_H))
    names(gsva_H)="value"
    gsva_H = merge(gsva_H,dt[,c("Row.names","OS.time","OS","PFI.time","PFI","age","gender","purity")],by.x="row.names",by.y="Row.names")
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
    
    if(input$gsea_cov=="None"){
      if(input$gsea_km_type=="OS"){
        fit = survfit(Surv(OS.time,OS)~cut,data = dt)
      }else{
        fit = survfit(Surv(PFI.time,PFI)~cut,data = dt)
      }
      p = ggsurvplot(fit,data = dt,pval=TRUE,legend.title="Median",legend.labs=c("high","low"),palette = c("#DF8F44FF","#374E55FF"))
      p = p$plot
    }else{
      dtt = data.frame(na.omit(dt))
      if(input$gsea_km_type=="OS"){
        if(input$gsea_cov=="Age + Sex"){
          res = adjusted_KM(data=dtt,time='OS.time',status="OS",group="cut",covlist=c("age","gender"),stratified_cox = "Yes",reference_group="G&B")
        }else if(input$gsea_cov=="Age + Sex + tumor_Purity"){
          res = adjusted_KM(data=dtt,time='OS.time',status="OS",group="cut",covlist=c("age","gender","purity"),stratified_cox = "Yes",reference_group="G&B")
        }
      }else{
        if(input$gsea_cov=="Age + Sex"){
          res = adjusted_KM(data=dtt,time='PFI.time',status="PFI",group="cut",covlist=c("age","gender"),stratified_cox = "Yes",reference_group="G&B")
        }else if(input$gsea_cov=="Age + Sex + tumor_Purity"){
          res = adjusted_KM(data=dtt,time='PFI.time',status="PFI",group="cut",covlist=c("age","gender","purity"),stratified_cox = "Yes",reference_group="G&B")
        }
      }
      
      p =  ggplot(res,aes(x=time,y = prob, group =class))+
        geom_step(aes(color = class),size=1.5)+
        theme_classic()+
        ylim(c(0,1))+
        theme(legend.position="top")+
        scale_color_manual(values=c("#DF8F44FF", "#374E55FF"))
      
    }
    p
  })
  
  output$gsea_km_plot = renderPlot({
    gsea_km_plot()
  })
  
  output$gsea_table = renderDT({
    res = gsea_table()
    names(res)[2] = "ssGSEA_value"
    res
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
  
  #===========================================================================
  ################## average KM ###############################
  #=========================================================================
  average_table = eventReactive(input$search_average_joint,{
    
    
    # dt = read_fst(paste0(folder_main,"www/pca/mRNA_clin.fst"),columns=c("Row.names","purity","type","age","gender","OS.time","OS","PFI.time",
    #                                                                     "PFI",make.names(input$gene_gsea)))
    dt =  rv$data
    dt_sub = subset(dt,dt$type==input$cancer_type)
    if(input$anchor_not=="No"){
      dt_sub_gene = dt_sub[,-c(1:10)]
    }else if(input$anchor_not=="Yes"){
      length_anchor = length(input$anchor_gene)
      dt_sub_gene = dt_sub[,-c(1:(10+length_anchor))]
    }
    dt_sub$value = rowMeans(dt_sub_gene,na.rm=T)
    dt_sub = data.frame(cbind(dt_sub[,1:10],"value"=dt_sub$value))
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
    
    if(input$gsea_cov_av=="None"){
      if(input$gsea_km_type_av=="OS"){
        fit = survfit(Surv(OS.time,OS)~cut,data = dt)
      }else{
        fit = survfit(Surv(PFI.time,PFI)~cut,data = dt)
      }
      p = ggsurvplot(fit,data = dt,pval=TRUE,legend.title="Median",legend.labs=c("high","low"),palette = c("#DF8F44FF","#374E55FF"))
      p = p$plot
    }else{
      dtt = data.frame(na.omit(dt))
      if(input$gsea_km_type_av=="OS"){
        if(input$gsea_cov_av=="Age + Sex"){
          res = adjusted_KM(data=dtt,time='OS.time',status="OS",group="cut",covlist=c("age","gender"),stratified_cox = "Yes",reference_group="G&B")
        }else if(input$gsea_cov_av=="Age + Sex + tumor_Purity"){
          res = adjusted_KM(data=dtt,time='OS.time',status="OS",group="cut",covlist=c("age","gender","purity"),stratified_cox = "Yes",reference_group="G&B")
        }
      }else{
        if(input$gsea_cov_av=="Age + Sex"){
          res = adjusted_KM(data=dtt,time='PFI.time',status="PFI",group="cut",covlist=c("age","gender"),stratified_cox = "Yes",reference_group="G&B")
        }else if(input$gsea_cov_av=="Age + Sex + tumor_Purity"){
          res = adjusted_KM(data=dtt,time='PFI.time',status="PFI",group="cut",covlist=c("age","gender","purity"),stratified_cox = "Yes",reference_group="G&B")
        }
      }
      
      p =  ggplot(res,aes(x=time,y = prob, group =class))+
        geom_step(aes(color = class),size=1.5)+
        theme_classic()+
        ylim(c(0,1))+
        theme(legend.position="top")+
        scale_color_manual(values=c("#DF8F44FF", "#374E55FF"))
      
    }
    p
  })
  output$average_km_plot = renderPlot({
    average_km_plot()
  })
  
  output$average_km_table = renderDT({
    res = average_table()
    #names(res)[2] = "average_value"
    res = data.frame(res[,c("Row.names","value","purity","OS","OS.time","PFI","PFI.time","age","gender")])
    names(res)[2]="average_value"
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
        paste0(input$cancer_type,"_joint_average_KM",'.png', sep = '')
        
      },
      content = function(file){
        ggsave(file, average_km_plot(),type="cairo-png",width = 12, height = 12,units="cm")
      }
    )
  })
  #----------------------------------------------------------------------------------------------------
  ######################################### EGA subnetwork #################################################
  #----------------------------------------------------------------------------------------------------
  subnetwork_res = eventReactive(input$search_subnet,{
    dt = rv$data
    if(input$anchor_not=="No"){
      dt_sub = subset(dt,dt$type==input$cancer_type)
      dt_sub_ega = dt_sub[,-c(1:10)]
    }else{
      dt_sub = subset(dt,dt$type==input$cancer_type)
      length_anchor = length(input$anchor_gene)
      dt_sub_ega = dt_sub[,-c(1:10)]
      names(dt_sub_ega)[1:length_anchor]=paste0("***",names(dt_sub_ega)[1:length_anchor],"***")
    }
    ega = EGA(dt_sub_ega,model = "glasso",plot.EGA = T)
    
  },ignoreNULL = F)
  
  output$subnetwork = renderPlot({
    subnetwork_res()
  })
  
  ###################### forest plot #########################
  forest_res = eventReactive(input$search_subnet,{
    dt = rv$data
    
    dt_sub = subset(dt,dt$type==input$cancer_type)
    dt_sub_ega = dt_sub[,-c(1:10)]
    
    res = subnetwork_res()$dim.variables
    res = na.omit(res)
    dt_sub_ega_t = data.frame(t(dt_sub_ega))
    res_m = merge(res,dt_sub_ega_t,by="row.names")
    for(i in 4:ncol(res_m)){
      res_m[,i] = log2(res_m[,i]+1)
    }
    res_mean = lapply(4:ncol(res_m),function(x){
      res =  res_m %>%
        group_by(dimension) %>% 
        summarise_at(colnames(res_m)[x], funs(mean(., na.rm=TRUE)))
      res[,-1]
    })
    
    res_final = do.call(cbind,res_mean)
    res_final_t = t(res_final)
    
    dt = rv$data
    
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
    univ_results$names=paste0("group",rownames(univ_results))
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
  # output$test = renderDT({
  #  rv$data
  # })
  #-------------------------------------------------------------------------------
  ############### user mannual ###############
  #-------------------------------------------------------------------------------

  
  output$tab_instructions<-renderUI({
    list(
      div(img(align="center",src=paste0("images/","overall.jpg"),height="715px"), style="text-align: center;"),
      br(),

      h3("Prognostic Ranking"),
      div(img(align="center",src=paste0("images/","prognostic_ranking.jpg"),height="780px"), style="text-align: center;"),
      br(),
      h3("Joint Signature"),
      div(img(align="center",src=paste0("images/","joint_sig.jpg"),height="120px"), style="text-align: center;"),
      br(),
      h3("Subnetwork"),
      div(img(align="center",src=paste0("images/","subnetwork.jpg"),height="150px"), style="text-align: center;"),
      br(),
      h3("Multi-gene Correlation"),
      div(img(align="center",src=paste0("images/","multi-gene_corr.jpg"),height="88px"), style="text-align: center;")
      
      

      
    )
  })

  #---------------------------- data preparation ------------------------------------
  omics_example = reactive({
    fread(paste0(folder_main,"www/user_upload/omics_example.csv"))
    
  })
  pheno_example = reactive({
    fread(paste0(folder_main,"www/user_upload/phenotype_example.csv"))
  })
  omics <- reactive({
    infile <- input$file_omics
    if (is.null(infile)) {
      omics = fread(paste0(folder_main,"www/user_upload/omics.csv"))
    }else{
      omics= fread(infile$datapath)
    }
    omics = data.frame(omics)
    colnames(omics)[1] = "geneName"
    return(omics)
  })
  
  pheno = reactive({
    infile <- input$file_pheno
    if (is.null(infile)) {
      fread(paste0(folder_main,"www/user_upload/phenotype.csv"))
    } else{
      fread(infile$datapath)
    }
  })
  
  merged_dt = reactiveValues(data=NULL)
  observe({
    omics = data.frame(omics())
    colnames(omics)[1] = "geneName"
    rownames(omics) = omics$geneName
    omics = data.frame(omics[,-1])
    omics = data.frame(t(omics))
    
    pheno = data.frame(pheno())
    names(pheno)[1] = "Row.names"
    merged_dt$data = merge(pheno,omics,by.x="Row.names",by.y="row.names")
    
  })
  
  #------------------ view the merged uploaded dataset ----------------------------
  output$upload_omics_view = renderDT({
    data.table(merged_dt$data[1:500,1:100],options = list(scrollX = TRUE,
                                                          initComplete = JS(
                                                            "function(settings, json) {",
                                                            "$(this.api().table().header()).css({'background-color': '#4B4B4B', 'color': '#fff'});",
                                                            "}"),
                                                          headerCallback = JS(
                                                            "function(thead, data, start, end, display){",
                                                            "  $(thead).css({'background-color': '#4B4B4B', 'color': '#fff'});",
                                                            "}")
                                                          ))
    
    
  })
  
  
  #------------------ substract the geneNames -------------------------
  observe({
    updateSelectizeInput(session,"gene_custom",choices = omics()$geneName,selected = head(omics()$geneName),server = T)
  })
  
  ##################### click the button to show the example omics (two small dataset)#####################
  output$omics_example = renderDT({
    datatable(omics_example(),options = list(scrollX = TRUE))
  })
  
  output$pheno_example = renderDT({
    pheno_example()
  })
  
  ############################################################################################################
  #-------------------------------------------------------------------
  ############# tab1: correlations #############################
  #-------------------------------------------------------------------
  dt_pca_custom = reactive({
    tryCatch({
      omics = data.frame(omics())
      rownames(omics) = omics$geneName; omics=data.frame(omics[,-1])
      dt = data.frame(t(omics))
      dt = data.frame(dt[,make.names(input$gene_custom)])
      dt
    },error = function(e){
      ""
    })
    
  })
  # ########################## PCA plot ############################################
  pca_plot_custom = reactive({
    tryCatch({
      res.pca <- prcomp(as.matrix(na.omit(dt_pca_custom())), scale = TRUE)
      p = fviz_pca_var(res.pca,arrowsize = 2, pointsize = 1.5,circlesize = 0.8,
                       col.var = "contrib", # Color by contributions to the PC
                       gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                       repel = TRUE     # Avoid text overlapping
      )+
        ggtitle("")
      p
    },error = function(e){
      "Please select genes"
    })
    
  })
  output$pca_plot_custom = renderPlot({
    pca_plot_custom()
  })
  
  ############################## Network plot ####################################
  # net_data = eventReactive(input$search_pca,{
  net_data_custom = reactive({
    tryCatch({
      adjm = cor(dt_pca_custom(),method = "spearman",use = 'pairwise.complete.obs')
      #validate(input$search_net)
      adjm[abs(adjm)<input$p_network_custom] <- 0
      
      adjm        
    },error=function(e){
      ""
    })
    
  })
  
  output$corr_network_plot_custom = renderPlot({
    tryCatch({
      network<- graph_from_adjacency_matrix(net_data_custom(), weighted=T, mode="undirected", diag=F)
      
      par(mar=c(0,0,0,0))
      plot(network,
           vertex.size=28,# Size of the node (default is 15)
           vertex.label.family="Helvetica", #or "Times"
           vertex.label.cex=0.9,
           # vertex.color=my_color,
           edge.curved=0.9)
    },error=function(e){
      "Please select genes"
    })
    
    
  },bg="grey13")
  
  ############################ corrplot ##########################################
  #corr_plot= eventReactive(input$search_pca,{
  corr_plot_custom = reactive({
    tryCatch({
      corr = round(cor(dt_pca_custom(),method = "spearman",use = 'pairwise.complete.obs'),3)
      
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
  
  
  output$corr_plot_custom = renderPlot({
    corr_plot_custom()
  })
  
  #---------------------------------------------------------------------------------------
  ###################### Grdient boosting ############################################
  #---------------------------------------------------------------------------------------
  output$gbm_surv_ui_custom = renderUI({
    req(merged_dt$data) 
    
    if (!is.null(merged_dt$data)) {
      h4("Data is available")
    } else {
      h4("Data is not available")
    }
    
    dt =   merged_dt$data
    dt = na.omit(dt)
    

    
    # Initialize a variable to hold the result string
    result <- ""
    
    # Check for the presence of OS and OS.time columns
    os_columns_present <- all(c("OS", "OS.time") %in% names(dt))
    
    # Check for the presence of PFI and PFI.time columns
    pfi_columns_present <- all(c("PFI", "PFI.time") %in% names(dt))
    
    # Construct the result string based on the presence of column sets
    if (os_columns_present && pfi_columns_present) {
      result <- c("OS","PFI")
    } else if (os_columns_present) {
      result <- "OS"
    } else if (pfi_columns_present) {
      result <- "PFI"
    }else{
      NULL
    }

    default_selection <- if (length(result) > 0) result[1] else NULL
    
    if (length(result) > 0) {
      
      column(3,  # Adjusted from 2 for better spacing, adapt as needed
             style = "display: flex; align-items: flex-start;",  # Align items in the center
             tags$div("OS or PFI", style = "margin-right: 5px; font-size: 16px; white-space: nowrap; line-height: 1.5;margin-top:5px;"),  # Label
             selectizeInput("gbm_surv_type_custom", NULL, choices = result, selected = default_selection)
      )

    } else {
      # Provide a fallback UI element or message if there are no valid options
     h4("No valid OS or PFI columns present")
    }
    
  })

  
  gbm_figure_custom = reactive({
    
    dt =   merged_dt$data
    dt = na.omit(dt)
    
    sel_genes = make.names(input$gene_custom)
    selected_columns <- c()
    check_columns <- c("OS.time", "OS", "PFI.time", "PFI")
    for(col in check_columns) {
      if(col %in% names(dt)) {  # Check for column presence
        # Add column if it exists, you can also add additional checks here (e.g., for completeness)
        selected_columns <- c(selected_columns, col)
      }
    }
    
    dt <- dt[, c(selected_columns,sel_genes)]
    
    if(input$gbm_surv_type_custom=="OS"){
      fit_formula = as.formula(paste("Surv(OS.time,OS)~", paste0(paste(sel_genes, collapse=" + "))))
    }else if(input$gbm_surv_type_custom=="PFI"){
      fit_formula = as.formula(paste("Surv(PFI.time,PFI)~", paste0(paste(sel_genes, collapse=" + "))))
    }
    gbm1 = gbm(fit_formula,data = dt,distribution="coxph")
    best.iter <- gbm.perf(gbm1,method="OOB") # returns test set estimate of best number of trees
    #summary(gbm1,n.trees=best.iter,cBars = length(input$gene_gbm))
    dt_p = summary(gbm1)
    dt_p = dt_p[order(dt_p$rel.inf,decreasing=T),]
    
    dt_p$var=rownames(dt_p)
    dt_p = subset(dt_p,dt_p$rel.inf>=input$gb_menu_custom)
    ggplot(dt_p,aes(x=reorder(var, rel.inf),y=rel.inf,fill=rel.inf))+
      geom_bar(stat="identity")+
      coord_flip()+
      theme_classic()+
      xlab("Relative influence")+
      ylab("")
  })
  
  
  
  
  output$gbm_plot_custom = renderPlot({
    tryCatch({
      gbm_figure_custom()
      
    },error = function(e){
      ""
    })
  })
  
  ############################ waterfall plot #############################
  #    cox_waterfall_plot = eventReactive(input$search_gbm,{
  cox_waterfall_plot_custom = reactive({
    dt =   merged_dt$data
    
    sel_genes = make.names(input$gene_custom)
    
    dt = na.omit(dt)
    
    selected_columns <- c()
    check_columns <- c("OS.time", "OS", "PFI.time", "PFI")
    for(col in check_columns) {
      if(col %in% names(dt)) {  # Check for column presence
        # Add column if it exists, you can also add additional checks here (e.g., for completeness)
        selected_columns <- c(selected_columns, col)
      }
    }
    
    dt <- dt[, c(selected_columns,sel_genes)]
    
    
    dt[sel_genes] <- lapply(dt[sel_genes], function(x) log2(x + 1)) 
    
    
    if(input$gbm_surv_type_custom=="OS"){
      #fit_formula = as.formula(paste("Surv(OS.time,OS)~", paste0(paste(sel_genes, collapse=" + "))))
      res = lapply(1:length(sel_genes),function(x){
        fit_formula = as.formula(paste("Surv(OS.time,OS)~", sel_genes[x]))
        fit = coxph(fit_formula,data=dt)
        res = data.frame(summary(fit)$coef)
        names(res) = c("coef","HR","se","z","p.value")
        res$gene = paste0(rownames(res),res$mark)
        res
      })
    }else if(input$gbm_surv_type_custom=="PFI"){
      # fit_formula = as.formula(paste("Surv(PFI.time,PFI)~", paste0(paste(sel_genes, collapse=" + "))))
      res = lapply(1:length(sel_genes),function(x){
        fit_formula = as.formula(paste("Surv(PFI.time,PFI)~", sel_genes[x]))
        fit = coxph(fit_formula,data=dt)
        res = data.frame(summary(fit)$coef)
        names(res) = c("coef","HR","se","z","p.value")
        res$gene = paste0(rownames(res),res$mark)
        res
      })
    }
    
    res = do.call(rbind,res)
    res$logHR = ifelse(res$coef>0,"> 0","<= 0")
    
    if(input$waterfall_menu_custom=="Yes"){
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
  
  output$waterfall_plot_custom = renderPlot({
    tryCatch({
      cox_waterfall_plot_custom()
      
    },error = function(e){
      ""
    })
  })   
  
  #-------------------------------------------------------------------------------
  ######################### joint signature #################################
  #_------------------------------------------------------------------------------
  output$average_km_ui_custom <- renderUI({
    req(merged_dt$data)
    dt <- na.omit(merged_dt$data)  # Make sure data is available and NA values are omitted
    
    # Initialize a variable to hold the result string
    result <- character(0)  # It's better to initialize as an empty character vector
    
    # Check for the presence of OS and OS.time columns
    os_columns_present <- all(c("OS", "OS.time") %in% names(dt))
    
    # Check for the presence of PFI and PFI.time columns
    pfi_columns_present <- all(c("PFI", "PFI.time") %in% names(dt))
    
    # Construct the result string based on the presence of column sets
    if (os_columns_present && pfi_columns_present) {
      result <- c("OS","PFI")
    } else if (os_columns_present) {
      result <- "OS"
    } else if (pfi_columns_present) {
      result <- "PFI"
    }
    
    # Ensure there is at least one valid option for selectizeInput to prevent errors
    if (length(result) == 0) result <- c("None")
    
    div(
      fluidRow(
        column(2,  # Adjusted from 2 for better spacing, adapt as needed
               style = "display: flex; align-items: flex-start;",  # Align items in the center
               tags$div("OS or PFI", style = "margin-right: 5px; font-size: 16px; white-space: nowrap; line-height: 1.5;margin-top:5px;"),  # Label
               selectizeInput("gsea_km_type_av_custom", NULL, choices = result, selected = result[1])
        ),
        column(2,  # Adjusted from 2 for better spacing, adapt as needed
               style = "display: flex; align-items: flex-start;",  # Align items in the center
               tags$div("KM cutoff",style = "margin-right: 5px; font-size: 16px; white-space: nowrap; line-height: 1;margin-top:5px;"),  # Keep label styling simple # Label
               selectizeInput("gsea_km_cutoff_av_custom", NULL, choices = c("median", "optimal", "quartile"), selected = "median")
        ),
        column(2,  # Adjusted from 2 for consistency, adapt as needed
               actionButton("search_average_joint_custom", "Click to Run", style = "width: 50%; height: 38px; background-color: #FFD700; color: black; border: none;")
        )
      )
    )
  
    # div(
    #   selectizeInput("gsea_km_type_av_custom","OS or PFI",result,result[1]),
    #   selectizeInput("gsea_km_cutoff_av_custom","Choose a cutoff",c("median","optimal","quartile"),"median"),
    #   actionButton("search_average_joint_custom", "Run")
    # )
    
  })
  
  output$gsea_km_ui_custom = renderUI({
    req(merged_dt$data)
    dt =   merged_dt$data
    dt = na.omit(dt)
    
    # Initialize a variable to hold the result string
    result <- ""
    
    # Check for the presence of OS and OS.time columns
    os_columns_present <- all(c("OS", "OS.time") %in% names(dt))
    
    # Check for the presence of PFI and PFI.time columns
    pfi_columns_present <- all(c("PFI", "PFI.time") %in% names(dt))
    
    # Construct the result string based on the presence of column sets
    if (os_columns_present && pfi_columns_present) {
      result <- c("OS","PFI")
    } else if (os_columns_present) {
      result <- "OS"
    } else if (pfi_columns_present) {
      result <- "PFI"
    }
    if (length(result) == 0) result <- c("None")
    
    div(
      fluidRow(
        column(2,  # Adjusted from 2 for better spacing, adapt as needed
               style = "display: flex; align-items: flex-start;",  # Align items in the center
               tags$div("OS or PFI:", style = "margin-right: 5px; font-size: 16px; white-space: nowrap; line-height: 1.5;margin-top:5px;"),  # Label
               selectizeInput("gsea_km_type_custom", NULL, choices = result, selected = result[1])
        ),
        column(2,  # Adjusted from 2 for better spacing, adapt as needed
               style = "display: flex; align-items: flex-start;",  # Align items in the center
               tags$div("KM cutoff",style = "margin-right: 5px; font-size: 16px; white-space: nowrap; line-height: 1;margin-top:5px;"),  # Keep label styling simple # Label
               selectizeInput("gsea_km_cutoff_custom", NULL, choices = c("median", "optimal", "quartile"), selected = "median")
        ),
        column(2,  # Adjusted from 2 for consistency, adapt as needed
               actionButton("search_gsea_custom", "Click to Run", style = "width: 50%; height: 38px; background-color: #FFD700; color: black; border: none;")
        )
      )
    )
    
    # div(
    #   selectizeInput("gsea_km_type_custom","OS or PFI",result,result[1]),
    #   selectizeInput("gsea_km_cutoff_custom","Choose a cutoff",c("median","optimal","quartile"),"median"),
    #   actionButton("search_gsea_custom", "Run")
    # )
    
    
  })
  
  #===========================================================================================
  ######################## SSGSEA ##################################################
  #===========================================================================================
  gsea_table_custom = eventReactive(input$search_gsea_custom,{
    
    
    dt= data.frame(omics())
    
    rownames(dt) = dt$geneName
    dt = data.frame(dt[,-1])
    
    #  list(c(input$gene_gsea))
    GE_matrix <- as.matrix(dt)
    
    
    gene_gsea = input$gene_custom
    
    gsva_H <- gsva(expr= GE_matrix, list(gene_gsea), method="ssgsea")
    gsva_H = data.frame(t(gsva_H))
    names(gsva_H)="value"
    
    pheno = data.frame(pheno())
    gsva_H = merge(gsva_H,pheno,by.x="row.names",by.y="Row.names")
    gsva_H
  })
  
  gsea_table_process_custom = reactive({
    dt = gsea_table_custom()
    if(input$gsea_km_cutoff_custom=="median"){
      dt$cut = findInterval(dt$value,median(dt$value,na.rm=T))
      dt$cut = ifelse(is.na(dt$cut),NA,ifelse(dt$cut ==1,"high","low"))
    }else if(input$gsea_km_cutoff_custom=="quartile"){
      dt$cut = findInterval(dt$value, quantile(dt$value,na.rm=T))
      dt$cut = ifelse(dt$cut%in%c(2,3),NA,ifelse(dt$cut%in%c(4,5),"high","low"))
    }else if(input$gsea_km_cutoff_custom=="optimal"){
      if(input$gsea_km_type_custom=="OS"){
        res.cut <- surv_cutpoint(dt,time = "OS.time", event = "OS", variables="value")
      }else{
        res.cut <- surv_cutpoint(dt,time = "PFI.time", event = "PFI", variables="value")
      }
      res.cat <- surv_categorize(res.cut)
      dt = data.frame(cbind(dt,"cut"=res.cat$value))
    }
    dt
  })
  
  gsea_km_plot_custom = reactive({
    
    dt = gsea_table_process_custom()
    
    if(input$gsea_km_type_custom=="OS"){
      fit = survfit(Surv(OS.time,OS)~cut,data = dt)
    }else{
      fit = survfit(Surv(PFI.time,PFI)~cut,data = dt)
    }
    p = ggsurvplot(fit,data = dt,pval=TRUE,legend.title="Median",legend.labs=c("high","low"),palette = c("#DF8F44FF","#374E55FF"))
    p = p$plot
    
    p
  })
  
  output$gsea_km_plot_custom = renderPlot({
    gsea_km_plot_custom()
  })
  
  output$gsea_table_custom = renderDT({
    res = gsea_table_custom()
    res
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
    output[["download_gsea_table_custom"]]<-downloadHandler(
      filename = function() {
        paste0(input$cancer_type,"_ssGSEA",'.csv', sep = '')
        
      },
      content = function(file){
        fwrite(gsea_table(),file,row.names = T)
      }
    )
    
  })
  
  observe({
    output[["download_gsea_km_custom"]]<-downloadHandler(
      filename = function() {
        paste0(input$cancer_type,"_ssGSEA_KM",'.png', sep = '')
        
      },
      content = function(file){
        ggsave(file, gsea_km_plot_custom(),type="cairo-png",width = 12, height = 12,units="cm")
      }
    )
  })
  
  
  #===========================================================================
  ################## average KM ###############################
  #=========================================================================
  average_table_custom = eventReactive(input$search_average_joint_custom,{
    
    omics = data.frame(omics())
    omics = subset(omics,omics$geneName%in%input$gene_custom)
    rownames(omics) = omics[,1]
    omics = data.frame(omics[,-1])
    omics = data.frame(t(omics))
    dt_sub_gene = omics 
    dt_sub_rowmean= data.frame(rowMeans(dt_sub_gene,na.rm=T))
    
    pheno = data.frame(pheno())
    names(dt_sub_rowmean) = "value"
    dt_sub = merge(pheno,dt_sub_rowmean,by.x="Row.names",by.y="row.names")
    dt_sub
    
  })
  
  average_table_group_custom = reactive({
    dt_sub = average_table_custom()
    if(input$gsea_km_cutoff_av_custom=="median"){
      dt_sub$cut = findInterval(dt_sub$value,median(dt_sub$value,na.rm=T))
      dt_sub$cut = ifelse(is.na(dt_sub$cut),NA,ifelse(dt_sub$cut ==1,"high","low"))
    }else if(input$gsea_km_cutoff_av_custom=="quartile"){
      dt_sub$cut = findInterval(dt_sub$value, quantile(dt_sub$value,na.rm=T))
      dt_sub$cut = ifelse(dt_sub$cut%in%c(2,3),NA,ifelse(dt_sub$cut%in%c(4,5),"high","low"))
    }else if(input$gsea_km_cutoff_av_custom=="optimal"){
      if(input$gsea_km_type_av_custom=="OS"){
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
  #   average_table_custom()
  # })
  average_km_plot_custom = reactive({
    
    dt = average_table_group_custom()
    
    if(input$gsea_km_type_av_custom=="OS"){
      fit = survfit(Surv(OS.time,OS)~cut,data = dt)
    }else{
      fit = survfit(Surv(PFI.time,PFI)~cut,data = dt)
    }
    p = ggsurvplot(fit,data = dt,pval=TRUE,legend.title="Median",legend.labs=c("high","low"),palette = c("#DF8F44FF","#374E55FF"))
    p = p$plot
    
    p
  })
  
  output$average_km_plot_custom = renderPlot({
    average_km_plot_custom()
  })
  
  output$average_table_custom_output = renderDT({
    res = average_table_custom()
    #res = data.frame(res[,c("Row.names","value","purity","OS","OS.time","PFI","PFI.time","age","gender")])
    
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
    output[["download_average_table_custom"]]<-downloadHandler(
      filename = function() {
        paste0(input$cancer_type,"_joint_average_custom",'.csv', sep = '')
        
      },
      content = function(file){
        fwrite(average_table_custom(),file,row.names = T)
      }
    )
    
  })
  
  observe({
    output[["download_average_km_custom"]]<-downloadHandler(
      filename = function() {
        paste0(input$cancer_type,"_joint_average_custom",'.png', sep = '')
        
      },
      content = function(file){
        ggsave(file, average_km_plot_custom(),type="cairo-png",width = 12, height = 12,units="cm")
      }
    )
  })
  
  #-------------------------------------------------------------------------------
  ###################### EGAnet #####################################
  #-------------------------------------------------------------------------------
  output$subnet_ui_custom = renderUI({
    div(
        column(2,
               actionButton("search_subnet_custom", "Click to Run", style = "width: 50%; height: 38px; background-color: #FFD700; color: black; border: none;")
        )
      
      
      
    )
    

  })
  
  #----------------------------------------------------------------------------------------------------
  ######################################### EGA subnetwork #################################################
  #----------------------------------------------------------------------------------------------------
  subnetwork_res_custom = eventReactive(input$search_subnet_custom,{
    omics = data.frame(omics())
    omics = subset(omics,omics$geneName%in%input$gene_custom)
    rownames(omics) = omics[,1]
    omics = data.frame(omics[,-1])
    omics = data.frame(t(omics))
    dt_sub_ega = omics
    ega = EGA(dt_sub_ega,model = "glasso",plot.EGA = T)
    
  },ignoreNULL = F)
  
  output$subnetwork_custom = renderPlot({
    subnetwork_res_custom()
  })
  
  ###################### forest plot #########################
  forest_res_custom = eventReactive(input$search_subnet_custom,{
    omics = data.frame(omics())
    sel_genes =input$gene_custom
    omics = subset(omics,omics$geneName%in%sel_genes)
    rownames(omics) = omics[,1]
    omics = data.frame(omics[,-1])
    omics = data.frame(t(omics))
    
    dt_sub_ega = omics
    
    res = subnetwork_res_custom()$dim.variables
    res = na.omit(res)
    dt_sub_ega_t = data.frame(t(dt_sub_ega))
    res_m = merge(res,dt_sub_ega_t,by="row.names")
    
    
    for(i in 4:ncol(res_m)){
      res_m[,i] = log2(res_m[,i]+1)
    }
    res_mean = lapply(4:ncol(res_m),function(x){
      res =  res_m %>%
        group_by(dimension) %>%
        summarise_at(colnames(res_m)[x], funs(mean(., na.rm=TRUE)))
      res[,-1]
      
      
    })
    
    res_final = do.call(cbind,res_mean)
    res_final_t = t(res_final)
    
    pheno = data.frame(pheno())
    res_final_t = merge(pheno,res_final_t,by.x = "Row.names",by.y = "row.names")
    
    covariates = names(res_final_t)[-c(1:3)]
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
    univ_results$names=paste0("group",rownames(univ_results))
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
  
  output$subnetwork_forest_custom = renderPlot(forest_res_custom())
  
}

# Run the application 
shinyApp(ui = ui, server = server)
