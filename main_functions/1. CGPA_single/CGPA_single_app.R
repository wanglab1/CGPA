#------------------------------------------------------------------------------
packages = c("shiny","rstudioapi","tidyverse","bslib","shinyWidgets","DT","shinycustomloader","shinyBS","fst","data.table","shinydashboard","shinydashboardPlus",
             "survival","survminer","circlize","broom","EGAnet","ComplexHeatmap","qgraph","plyr","DBI","RSQLite") # rstudioapi, not needed

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

folder_main = "/CGPA_single"

setwd(folder_main)
source("functions.R")
source("Pan-cancer_summary.R")
source("KM_surv_adj.R")
source("Multivariable_Analysis_uni.R")
source("Multivariable_Analysis_multi.R")
source("Multivariable_Analysis_Top_prognostic_genes.R")
source("Gene-hallmark_interation.R")
source("lncRNA_top_prog_across_cancers.R")
source("lncRNA_TIDE.R")
source("lncRNA_TIDE_km.R")
source("lncRNA_exploration_coexpression.R")

dark = bs_theme(version = 3,bg="black",fg = "white",warning = "#FFD300" )%>%
  bs_add_rules(sass::sass_file("www/style/style.scss"))

ui <- fluidPage(
  
  theme = dark,
  useShinydashboard(),
  # Include custom CSS rules
  
  div(class="navbar1",
      navbarPage(title = "CANCER GENE PROGNOSIS ATLAS",
                 
                 fluid = TRUE, 
                 collapsible = TRUE,
                 #==============================================================================================================================================================
                 #################################### Tab 1: CGPA dashboard ##################################################################################################
                 #==============================================================================================================================================================                           
                 # --------------------------------------------------
                 # tab panel 1 - Main
                 #---------------------------------------------------
                 tabPanel("Pan-Cancer Summary",
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
                                              div(style="text-align: center;",imageOutput("forest")),
                                              br(),
                                              div(style="text-align: left;",
                                                  textOutput("summary_forest1"),
                                                  textOutput("summary_forest2"))
                                              
                                         )) 
                                       )
                                   )
                                   
                            ),
                            column(5, pan_ui("box2","Gene Expression Profile",600,"cir_bar")),
                            column(7,pan_ui("box2","Kanplan-Meier Plots (OS optimal cutoff)",530,"KM_general")),
                            column(5, pan_ui("box4","PPI network (STRING)",540,"string_net"))
                          )
                          
                 ),
                 tabPanel("Multivariable Analysis", 
                          br(),
                          tags$style(HTML("
                                            .tabbable > .nav > li > a {font-size:20px; font-family:'Arial' serif;}
                                            .navbar.navbar-default.navbar-static-top{font-size: 18px;}
                                            .navbar .navbar-header {font-size:35px;font-weight:bold;text-align:left;width:350px;}
                                            #header {text-align: center;color: #fdfdfd;text-shadow: 0 0 1px #000;font-weight: bold;}
                                            h3{font-family: 'Montserrat', sans-serif;z-index: 2;}
                                            .navbar-nav {float: right; }
                                            .navbar-nav > li:head > a {font-size: 20px; font-weight: bold;}
                                           ")),
                          
                          fluidRow(
                            column(2,
                                   selectizeInput("cancer_type","Cancer type",
                                                  choices=    c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA",
                                                                "GBM","HNSC","KICH","KIRC","KIRP","LGG","LIHC",
                                                                "LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ",
                                                                "SARC","SKCM","STAD","TGCT","THCA","UCEC","UCS","UVM",
                                                                "BRCA-BASAL","BRCA-NON-BASAL","SKCM_Metastasis"),multiple=T,selected="HNSC"),
                                   
                                   sel_input_ui("input_gene_dash","Select a gene"),
                                   sel_input_ui("survtype","OS or PFI")
                            ),
                            column(10,
                                   tabsetPanel(type="pills",id="mainnav", 
                                               #---------------------------------------------------
                                               ################ Tab 2: survival ##################
                                               #---------------------------------------------------
                                               tabPanel("MULTIVARIABLE ANALYSIS",
                                                        br(),
                                                        uiOutput("multi_ui"),
                                                        
                                                        fluidRow(
                                                          column(9,br(),fluidRow(column(12,uiOutput("survival_dash")),
                                                                                 DTOutput("testt")
                                                                                 
                                                          ))
                                                        
                                                         
                                                        )
                                                        
                                               ),
                                               
                                               #-------------------------------------------------------------------------------
                                               ##################### Tab3: Most prognostic gene ###############################
                                               #-------------------------------------------------------------------------------
                                               tabPanel("TOP PROGNOSTIC GENES",
                                                        br(),
                                                        fluidRow(
                                                          uiOutput("top_prog_ui_ui")
                                                        ),
                                                       # br(),
                                                        fluidRow(
                                                          column(9,
                                                                 h5(strong("Please select one cancer type for this part")),
                                                                 fluidRow(
                                                                   column(6,top_prog_ui(box_id="box1",title="Top prognostic genes (mRNA)",condition="search_top_within",output="most_sig_genes1",drop_down_info="Click to see subnetwork",subnet_action="subnet",output_subnet="subnetwork_res")),
                                                                   column(6, top_prog_ui(box_id="box1",title="Top prognostic lncRNAs",condition="search_top_within",output="most_sig_genes2",drop_down_info="Subnetwork is not avaliable for lncRNAs"))
                                                                 )
                                                          ),
                                                          column(3,
                                                                 p("NOTE: Adjusting for tumor purity is not avaliable for CHOL, ESCA, MESO, PCPG, SARC, TGCT, YHYM and UVM, adjusting for sex is not avaliable for CESC, OV, PRAD, TGCT, UCEC, UCS, BASAL-BRCA and NON-BASAL-BRCA'"), 
                                                                 
                                                                 )
                                                         
                                                        )
                                                        
                                                        
                                               )
                                   )
                            )
                            
                            
                          )
                          
                 ),
                 #==============================================================================================================================================================
                 #################################### Tab 3: GHI model  ##################################################################################################
                 #==============================================================================================================================================================
                 tabPanel("Gene-Hallmark Interaction",
                          br(),
                          tabsetPanel(type="pills",id="mainnav_ghi",
                                      # tabPanel("COMMON GENES DISCOVERY",
                                      #          br(),
                                      #          br(),
                                      #          fluidRow(
                                      #            column(3,
                                      #                   selectizeInput("cancer_type_ghi_common","Cancer type",
                                      #                                  choices=    c("LUAD","BLCA","BRCA","CESC","HNSC","KIRC","KIRP","LGG","LIHC","LUSC","OV","PRAD",
                                      #                                                "SKCM","STAD","THCA","UCEC"
                                      #                                  ),multiple=T,selected="HNSC"),
                                      #                   selectizeInput("hallmark_match_common","Choose hallmarks",choices = NULL,multiple=T),
                                      #                   numericInput("top_n_ghi","Top n numbers",100),
                                      #                   radioButtons("tumor_purity_ghi_common","Adjust for tumor purity?",c("Yes","No"),"Yes"),
                                      #                   actionButton("ghi_common_run","Run")
                                      #                   
                                      #            ),
                                      #            
                                      #            column(9,
                                      #                   p("Discover common genes that appear in at least 2 cancer types or 2 hallmarks for overall survival"),
                                      #                   DTOutput("ghi_common_tb", width = "60%")
                                      #                   )
                                      #          )
                                      #   
                                      #          ),
                                      tabPanel("SINGLE GENE",
                                               br(),
                                               fluidRow(
                                                 column(3,
                                                        selectizeInput("cancer_type_ghi","Cancer type",
                                                                       choices=    c("BLCA","BRCA","CESC","HNSC","KIRC",
                                                                                     "KIRP","LGG","LIHC","LUAD","LUSC","OV",
                                                                                     "PRAD","SKCM","STAD","THCA","UCEC"),multiple=F,selected="HNSC"),
                                                        sel_input_ui("gene_single_ghi","Choose a gene",choices=NULL,multiple=F),
                                                        
                                                        selectizeInput("survtype_ghi","OS or PFI",c("OS","PFI"),"OS"),
                                                        sel_input_ui("hallmark_match","Choose hallmarks",choices = NULL,multiple=T),
                                                        radioButtons("tumor_purity_ghi","Adjust for tumor purity?",c("Yes","No"),"Yes"),
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
                                                              # 
                                                              # conditionalPanel(
                                                              #   condition= paste0("input.","search_single_ghi",">0"),
                                                              #   p("Interaction model"),
                                                              #   div(withLoader(plotOutput("edge_plot",width="80%",height="400px"),type="html",loader="loader1"),align="center")
                                                              # ),
                                                              # conditionalPanel(
                                                              #   condition= paste0("input.","search_single_ghi","==0"),
                                                              #   br(),
                                                              #   p("Click 'Run' button to get/update the figure")
                                                              # )
                                                            )
                                                            
                                                        )
                                                 )
                                               )
                                               
                                               # uiOutput("single_ghi")
                                      ),
                                      tabPanel("MULTIPLE GENES",
                                               br(),
                                               uiOutput("multi_ghi")
                                      )
                                      
                          )   
                 ),
                 #=====================================================================================================
                 #################################### Tab4: progsplicing ##############################
                 #=====================================================================================================
                 tabPanel("ProgSplicing",
                          br(),
                          tabsetPanel(type="pills",id="mainnav_splice",
                                      fluidRow(
                                        column(3,
                                               
                                               sel_input_ui("input_gene_splice","Select a gene",choices = NULL,multiple=F),  
                                               selectizeInput("surv_type_splice","OS or PFI",c("OS","PFI"),"OS"),
                                               selectizeInput("cutoffs_splice","Cutoffs",c("Median",'Optimal'),"Optimal"),
                                               selectizeInput("database_splice","Database (From OncoSplicing)",c("SpliceSeq","SplAdder"),"SpliceSeq")
                                              # actionButton("Run_splice", "Run")
                                        ),
                                        column(8,
                                               div(id="box1",
                                                   box(
                                                     title="Prognostic Landscape",collapsible = TRUE,status="warning",
                                                     solidHeader = TRUE,
                                                     width = 12,height=800,
                                                     dropdownButton(
                                                       numericInput("text_size_splice","Adjust the size of text",11,min=1,max=30,step = 1),
                                                       downloadButton("download_splice", "Download"),
                                                       
                                                       circle = TRUE, status = "warning", icon = icon("gear"), width = "300px",
                                                       tooltip = tooltipOptions(title = "Click to see inputs !")
                                                     ),
                                                     
                                                     withLoader(plotOutput("splice_plot",height="800px"),type="html",loader="loader1")
                                                     
                                                   )
                                               ),
                                               div(id="box1",
                                                   box(
                                                     title="Splice Event Annotation (From OncoSplicing)",collapsible = TRUE,status="warning",
                                                     solidHeader = TRUE,
                                                     width = 12,height=600,
                                                     dropdownButton(
                                                       downloadButton("download_splice_tb", "Download"),
                                                       circle = TRUE, status = "warning", icon = icon("gear"), width = "300px",
                                                       tooltip = tooltipOptions(title = "Click to see inputs !")
                                                     ),
                                                     withLoader(DTOutput("splice_dt"),type="html",loader="loader1")
                                                   )
                                               )
                                        )
                                      ))

                          ),
                 #==============================================================================================================================================================
                 #################################### Tab 5: lncRNA toolbox ##################################################################################################
                 #==============================================================================================================================================================
                 tabPanel("LncRNA Exploration",
                          
                          br(),
                          tabsetPanel(type="pills",id="mainnav",
                                      #---------------------------------------------------------------------------
                                      ##################### Tab4: TIDE lncRNA and mRNA #########################
                                      #--------------------------------------------------------------------------
                                      tabPanel("TOP PROGNOSTIC GENES ACROSS CANCERS",
                                               br(),
                                               uiOutput("common_genes"),
                                               fluidRow(
                                                 column(3),
                                                 column(9, p("Click 'Run' to get/update the figure"))
                                               )
                                               
                                               
                                               
                                      ),
                                      tabPanel("TIDE - OS INTERACTION",
                                               br(),
                                               tabsetPanel(type="pills",id="mainnav_tide",
                                                           
                                                           tabPanel("Explore common genes",
                                                                    fluidRow(
                                                                      column(3),column(9,p("Click 'Run' button to get/update results"))
                                                                    ),
                                                                    uiOutput("tide_res")
                                                           ),
                                                           tabPanel("KM plot of interaction",
                                                                    fluidRow(
                                                                      column(3),column(9,p("Take a second for the figure"))
                                                                    ),
                                                                    uiOutput("tide_km")
                                                           )
                                               )
                                      ),
                                      
                                      tabPanel("COEXPRESSION ANALYSIS",
                                               br(),
                                               fluidRow(
                                                 column(3,
                                                        selectizeInput("input_gene_cor","select a lncRNA",choices = NULL,multiple=FALSE),
                                                 ),column(9,p("Click 'Run' button to get/update results"))
                                               ),
                                               uiOutput("corr_genes")
                                               
                                               
                                               
                                               
                                               
                                               
                                      )
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
      )
  )
)

server = function(input,output,session){
  observe({
    query = parseQueryString(session$clientData$url_search)
    cat('input_gene ', query$geneid, '\n')
    
    if (!is.null(query[['geneid']])) {
      updateTextInput(session, "input_gene", "",value = toupper(query[['geneid']]))
    }

  })
  
  
  
  #============================================================================================================================================
  ########################################################## Pan cancer summary  #########################################################
  #==========================================================================================================================================
  ############################### Forest+heatplot ############################################
  #forest = eventReactive(input$search,{
  forest = reactive({
    pan_server_fig(box_id="box1",file = "www/forest/",file_name = "_forest_bubble_fig.png",input_gene = input[["input_gene"]])
  })
  output$forest <- renderImage({
    forest()
  },deleteFile=FALSE)
  
  ########################### summarize OS ############################
  uni_cox1 <- reactive({
    tryCatch({
      
      file ="www/files/unicox/"
      res = fread(paste0(file,make.names(input$input_gene),"_","OS.csv"))
      res
    },error = function(e){
      res = data.frame(matrix(vector(), 0, 3,
                              dimnames=list(c(), c("Date", "File", "User"))),
                       stringsAsFactors=F)
      res
    })
    
  })
  
  
  summary_forest1 = reactive({
    pan_server_text(input_data = uni_cox1(),type="OS")
    
  })
  output$summary_forest1 = renderText({
    paste("OS:", summary_forest1())
  })
  
  ############### Summarize forest plot (PFI) ####################
  uni_cox2 <- reactive({
    tryCatch({
      file = "www/files/unicox/"
      res = fread(paste0(file,make.names(input$input_gene),"_","PFI.csv"))
      res
    },error = function(e){
      res = data.frame(matrix(vector(), 0, 3,
                              dimnames=list(c(), c("Date", "File", "User"))),
                       stringsAsFactors=F)
      res
    })
  })
  
  summary_forest2 = reactive({
    pan_server_text(input_data = uni_cox2(),type="PFI")
    
  })
  output$summary_forest2 = renderText({
    paste("PFI:", summary_forest2())
  })
  ############################# CIRCULATE BAR ##############################################
  cir_bar = reactive({
    pan_server_fig(box_id="box2",file = "www/circulate_bar/",file_name = "_cirbar.png",input_gene = input[["input_gene"]])
  })
  output$cir_bar <- renderImage({
    cir_bar()
  },deleteFile=FALSE)
  ################################## General KM ###################################
  KM_general = reactive({
    pan_server_fig(box_id="box2",file = "www/KM_home/",file_name = "_OS_KM.png",input_gene = input[["input_gene"]])
  })
  output$KM_general = renderImage({
    KM_general()
  },deleteFile=FALSE)
  ##################################### string network ######################################
  src= reactive({
    pan_server_fig(box_id="box4",input_gene = input[["input_gene"]])
  })
  
  output$string_net<-  renderUI({
    tags$img(src = src(),alt="string network",width=1000,height=400)
  }) 
  
  #============================================================================================================================================
  ########################################################## CGPA dash #########################################################
  #============================================================================================================================================
  ###########update genesearch bar for tab2 ######
  geneName = reactive({
    read_fst("www/geneName.fst")
  })
  geneName_mrna = reactive({
    read_fst("www/geneName_mrna.fst")
  })
  geneName_lnc = reactive({
    read_fst("www/geneName_lnc.fst")
  })
  
  cancer_type = reactive({
    c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA",
      "GBM","HNSC","KICH","KIRC","KIRP","LGG","LIHC",
      "LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ",
      "SARC","SKCM","STAD","TGCT","THCA","UCEC","UCS","UVM")
  })
  observe({
    sel_input_server(input,output,session,"input_gene_dash","Select a gene",geneName()$geneName,input$input_gene)  
  })
  
  observe({
    sel_input_server(input,output,session,"survtype","OS or PFI",c("OS","PFI"),"OS")
  })
  
  
  
  #-----------------------------------------------------------------------------------------------------------------
  ######################################### KM and adj KM #################################################
  #------------------------------------------------------------------------------------------------------------------
  output$multi_ui = renderUI({ # UI for KM and unadjusted KM
    div(
      fluidRow(
        column(2,  # Adjusted from 2 for better spacing, adapt as needed
               style = "display: flex; align-items: flex-start;",  # Align items in the center
               tags$div("KM or adj.KM", style = "margin-right: 5px; font-size: 16px; white-space: nowrap; line-height: 1.5;margin-top:5px;"),  # Label
               selectizeInput("km_type", NULL, choices = c("KM","Adjusted_KM"), multiple=F)
        ),
        column(2,  # Adjusted from 2 for better spacing, adapt as needed
               style = "display: flex; align-items: flex-start;",  # Align items in the center
               tags$div("cutoff for KM",style = "margin-right: 5px; font-size: 16px; white-space: nowrap; line-height: 1;margin-top:5px;"),  # Keep label styling simple # Label
               selectizeInput("cutoff_multi", NULL, choices = c("median", "optimal", "quartile"), multiple=F)
        ),
        column(3,
               conditionalPanel(
                 condition="input.km_type=='Adjusted_KM'",
                 style = "display: flex; align-items: flex-start;",  # Align items in the center
                 tags$div("Adjusted covariates",style = "margin-right: 5px; font-size: 16px; white-space: nowrap; line-height: 1;margin-top:5px;"),  # Keep label styling simple # Label
                 selectizeInput("adj_cov_multi",NULL,choices =c("tumor_purity","age","sex","CTL"),multiple=T,selected="tumor_purity"),
               )
               )

      ),
      fluidRow(
        p("NOTE: Adjusting for tumor purity is not avaliable for CHOL, ESCA, MESO, PCPG, SARC, TGCT, YHYM and UVM, adjusting for sex is not avaliable for CESC, OV, PRAD, TGCT, UCEC, UCS, BASAL-BRCA and NON-BASAL-BRCA'") 
        
      )
    )
    
    # div(
    #   sel_input_ui("km_type","KM or adj.KM",choices=c("KM","Adjusted_KM"),multiple=F),
    #   sel_input_ui("cutoff_multi","Cutoff for KM",choices=c("Quartile","Optimal","Median"),multiple=F),
    #   conditionalPanel(
    #     condition="input.km_type=='Adjusted_KM'",
    #     sel_input_ui("adj_cov_multi","Adjusted covariates",choices =c("tumor_purity","age","sex","CTL"),multiple=T,selected="tumor_purity"),
    #     p("NOTE: Adjusting for tumor purity is not avaliable for CHOL, ESCA, MESO, PCPG, SARC, TGCT, YHYM and UVM, adjusting for sex is not avaliable for CESC, OV, PRAD, TGCT, UCEC, UCS, BASAL-BRCA and NON-BASAL-BRCA'") 
    #   )
    # )
  })
  ############################# unadjusted KM and cox ##############################
  
  ########## KM plot #############
  output$survival_dash = renderUI({
    
    surv_list = lapply(1:length(input$cancer_type),function(x){
      div(
        survival_UI(headings=input$cancer_type[[x]],output_info=paste0("survival",x),
                    more_cox= paste0("morecox",x)),
        style="display: inline-block; ")
    })
    tagList(surv_list)
  })
  
  input_data_adjsurv <- reactive({
    lapply(1:length(input$cancer_type),function(x){
      survival_adj_data(input$cancer_type[[x]],input$adj_cov_multi,input$cutoff_multi,input$input_gene_dash,geneName())
    })
  })
  
  observe({
    # observeEvent(input$search_multi,{
    tryCatch({
      if(input$km_type=="KM"){
        tryCatch({
          lapply(1:length(input$cancer_type),function(x){
            output[[paste0("survival",x)]]<-renderPlot({
              survival_server1(input,output,session,geneName =geneName(),
                               cancer_multi=input$cancer_type[[x]])
            })
          }
          )
        },error = function(e){
          "Create an empty image with 'please input a valid geneName'"
        })
      }else{
        tryCatch({
          lapply(1:length(input$cancer_type),function(x){
            output[[paste0("survival",x)]]<-renderPlot({
              survival_adj_server1(input,output,session,inputdata=input_data_adjsurv()[[x]],survtype=input$survtype,cutoff=input$cutoff_multi,input_gene_surv=make.names(input$input_gene_dash),adj_cov=input$adj_cov_multi)
            })
          }
          )
        },error = function(e){
          "Create an empty image with 'please input a valid geneName'"
        })
      }
      
    },error=function(e){
      ""
    })
    
    
  })
  
  ########## coxph result ##############
  input_data_adjsurv_cox <- reactive({
    lapply(1:length(input$cancer_type),function(x){
      survival_adj_data(input$cancer_type[[x]],input$adj_cov_multi,"raw",input$input_gene_dash,geneName())
    })
  })
  
  observe({
    # observeEvent(input$search_multi,{
    tryCatch({
      if(input$km_type=="KM"){
        tryCatch({
          lapply(1:length(input$cancer_type),function(x){
            output[[paste0("morecox",x)]]<-renderDataTable({
              res = survival_server2(input,output,session,geneName =geneName(),cancer_multi=input$cancer_type[[x]])
              names(res) = c("Gene","HR (95% CI)","p_value","n(events)")
              DT:::datatable(
                data.frame(res),rownames = FALSE,colnames=FALSE,options = list(dom='t',scrollX=TRUE,ordering=F,
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
              #%>% formatStyle(columns=colnames(res),color='white',background = '#4B4B4B',target = 'row')
              
            })
          }
          )
        },error = function(e){
          "Create an empty image with 'please input a valid geneName'"
        })
      }else{
        tryCatch({
          lapply(1:length(input$cancer_type),function(x){
            output[[paste0("morecox",x)]]<-renderDT({
              
              res = survival_adj_server2(input,output,session,inputdata=input_data_adjsurv_cox()[[x]],survtype=input$survtype,input_gene_surv=make.names(input$input_gene_dash),adj_cov_input=input$adj_cov_multi)
              # datatable(data.frame(res))
              DT:::datatable(
                data.frame(res),rownames = FALSE,colnames=FALSE,options = list(dom='t',scrollX=TRUE,ordering=F,
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
          )
        },error = function(e){
          "Create an empty image with 'please input a valid geneName'"
        }) 
      }
    },error = function(e){
      "Please Wait..."
    })
    
    
  })
  
  

  
  #-------------------------------------------------------------------------------
  ####################### most prognostic genes tab1#################################
  #-------------------------------------------------------------------------------
  output$top_prog_ui_ui = renderUI({
    div(
      fluidRow(
        column(3,  # Adjusted from 2 for better spacing, adapt as needed
               style = "display: flex; align-items: flex-start;",  # Align items in the center
               tags$div("Adjusted covariates", style = "margin-right: 5px; font-size: 16px; white-space: nowrap; line-height: 1.5;margin-top:5px;"),  # Label
               selectizeInput("surv_cov_most", NULL, choices =c("None","Age + Sex","Age + Sex + tumor_Purity","Age + Sex + CTL"), selected = "None",multiple=F)
        ),
        column(2,  # Adjusted from 2 for consistency, adapt as needed
               actionButton("search_top_within", "Click to Run", style = "width: 50%; height: 38px; background-color: #FFD700; color: black; border: none;")
        )
      )

    )
    
  })
  ############## top prognostic table ###############
  most_sig_genes_tb1 = eventReactive(input$search_top_within,{
    if(length(input$cancer_type)>1){
      validate("Please select one cancer type at a time, and click 'Run' button")
    }
    tryCatch({
      
      top_prog_server(file= "www/mrna_cox_all/",adj_cov=input$surv_cov_most,cancer_type=input$cancer_type,surv_type=input$survtype)
      
    },error=function(e){
      data.frame("Cancer type not exist")
    })
    
  })
  
  output$most_sig_genes1 = renderDT({
    if(length(input$cancer_type)>1){
      validate("Please select one cancer type at a time, and click 'Run' button")
    }
    
    tryCatch({
      res = most_sig_genes_tb1()[-c(1:2),]
      # # change names
      geneName = geneName_mrna()
      geneName$make.names = make.names(geneName$geneName)
      res = merge(geneName,res,by.x = "make.names",by.y = "Genes")
      names(res)[2] = "Genes"
      res = res[,-1]
      res = res[order(res$p.value,decreasing = FALSE),]
      
      DT:::datatable(
        res,rownames = FALSE,colnames=FALSE,escape = TRUE,selection ="single",options = list(scrollX=TRUE,scrollY="350px",
                                                                                             initComplete = JS(
                                                                                               "function(settings, json) {",
                                                                                               "$(this.api().table().header()).css({'background-color': 'white', 'color': '#4B4B4B'});",
                                                                                               "}")
        )
        ,
        container = tags$table(
          class="stripe row-border hover",
          tags$thead(tags$tr(lapply(colnames(res), tags$th)))
        )
      ) %>% formatStyle(columns=colnames(res),color='#4B4B4B',background = 'white',target = 'row')  %>%
        formatSignif(columns=colnames(res)[-1], digits=3)
    },error = function(e){
      NULL
    })
    
    
  })
  
  most_sig_genes_tb2 = eventReactive(input$search_top_within,{
    if(length(input$cancer_type)>1){
      validate("Please select one cancer type at a time, and click 'Run' button")
    }
    tryCatch({
      top_prog_server(file= "www/lnc_cox_all/",adj_cov=input$surv_cov_most,cancer_type=input$cancer_type,surv_type=input$survtype)
      
    },error=function(e){
      data.frame("Cancer type not exist")
    })
    
  })
  
  output$most_sig_genes2 = renderDT({
    if(length(input$cancer_type)>1){
      validate("Please select one cancer type at a time, and click 'Run' button")
    }
    res = most_sig_genes_tb2()[-c(1:2),]
    # # change names
    geneName = geneName_lnc()
    geneName$make.names = make.names(geneName$geneName)
    res = merge(geneName,res,by.x = "make.names",by.y = "Genes")
    names(res)[2] = "Genes"
    res = res[,-1]
    res = res[order(res$p.value,decreasing = FALSE),]
    tryCatch({
      DT:::datatable(
        res
        ,rownames = FALSE,colnames=FALSE,escape = TRUE,selection ="single",options = list(scrollX=TRUE,scrollY="350px",
                                                                                          initComplete = JS(
                                                                                            "function(settings, json) {",
                                                                                            "$(this.api().table().header()).css({'background-color': 'white', 'color': '#4B4B4B'});",
                                                                                            "}")
        ),
        container = tags$table(
          class="stripe row-border hover",
          tags$thead(tags$tr(lapply(colnames(res), tags$th)))
        )
      )  %>% formatStyle(columns=colnames(res),color='#4B4B4B',background = 'white',target = 'row') %>%
        formatSignif(columns=colnames(res)[-1], digits=3)
    },error = function(e){
      NULL
    })
    
  })
  
  
  #################### subnetwork ##########################
  dt_net = reactive({
    if(length(input$cancer_type)>1){
      validate("Please select one cancer type at a time")
    }
    
    top_subnet(top_tb = most_sig_genes_tb1(),input$input_gene_dash,cancer_type=input$cancer_type)
  })
  
  subnetwork = eventReactive(input$subnet,{
    
    dt_sub_ega =  dt_net()
    ega = EGA(dt_sub_ega,model = "glasso",plot.EGA = T)
    
  },ignoreNULL = F)
  
  output$subnetwork_res = renderPlot({
    subnetwork()
  })
  
  #------------------------------------------------------------------------------------
  ############################### GHI single ##########################################
  #------------------------------------------------------------------------------------
  # output$single_ghi = renderUI({
  #     ghi_ui_single("box4",title="Gene-hallmark interaction model (z-score)",cancer_type_ghi="cancer_type_ghi",
  #                   gene_single_ghi="gene_single_ghi",survtype_ghi="survtype_ghi",hallmark_match="hallmark_match",tumor_purity_ghi="tumor_purity_ghi",
  #                   output_dropdown = "GHI_zscore",output_plot = "edge_plot",action_button="search_single_ghi")
  # })
  
  top_mut = reactive({
    fread("www/GHI/mutated_genes.csv")
  })
  
  top_cna = reactive({
    fread("www/GHI/CNA.csv")
    
  })  
  
  ghi_top_mut = reactive({
    res = top_mut()[[input$cancer_type_ghi]]
    res[res==""]=NA
    res = res[complete.cases(res)]
    res
  })
  ghi_top_cna = reactive({
    res = top_cna()[[input$cancer_type_ghi]]
    res[res==""]=NA
    res = res[complete.cases(res)]
    res
  })
  
  observe({
    tryCatch({
      #trigger = input$mainnav_ghi
      sel_input_server(input,output,session,"hallmark_match","Choose hallmarks for interaction models",
                       input_res =  list(
                         "Clinical covariates" = unique(c("age","sex","Tumor_burden")),
                         "Top mutated genes" = paste0(unique(ghi_top_mut()),"_mut"),
                         "Top CNA" = paste0(unique(ghi_top_cna()),"_cna"),
                         "GE" = unique(c("IFN.Hallmark_ssgsea","IFN.Hallmark_mean","ISG.RS_ssgsea",
                                         "ISG.RS_mean","Hypoxia_mean","Hypoxia_ssgsea","BUFFA_HYPOXIA_SCORE")),
                         "Others" = unique(c("CYT","CTL","AR"))
                       ),sel_input =  unique(c("age","sex","Tumor_burden",( paste0(ghi_top_mut(),"_mut")[1:3]),
                                               (paste0(ghi_top_cna(),"_cna")[1:3]),"IFN.Hallmark_ssgsea","ISG.RS_ssgsea","CYT","CTL","AR")))
    },error=function(e){
      ""
    })
    
  })
  
  observe({
    #trigger = input$mainnav_ghi
    sel_input_server(input,output,session,"gene_single_ghi","Choose a gene",input_res=geneName()$geneName,sel_input =input$input_gene)
  })
  
  # GE data
  dt_match = reactive({
    if(length(input$cancer_type_ghi)>1){
      validate("Please select one cancer type at a time and click 'Run' button")
    }
    
    if(length(which(geneName()$geneName==input$gene_single_ghi))>1){
      cov_list = unique(c(paste0(make.names(input$gene_single_ghi),".x"),make.names(input$hallmark_match)))
      dt = read_fst(paste0("www/GHI/processed_app/",input$cancer_type_ghi,".fst"),columns=c("Patient.ID","OS.time","OS",
                                                                                            "PFI.time","PFI",cov_list,"Tumor_purity"))
      names(dt)[which(names(dt)==c(paste0(make.names(input$gene_single_ghi),".x")))] = make.names(input$gene_single_ghi)
    }else{
      cov_list = unique(c(make.names(input$gene_single_ghi),make.names(input$hallmark_match)))
      dt = read_fst(paste0("www/GHI/processed_app/",input$cancer_type_ghi,".fst"),columns=c("Patient.ID","OS.time","OS",
                                                                                            "PFI.time","PFI",cov_list,"Tumor_purity"))
    }
    
    # if(input$cancer_type_ghi=="HNSC"){
    #   names(dt)[which(names(dt)=="cancer_subtype")] = "HPV_status"
    # }
    dt
  })
  
  ############ edge plot ##############
  res_edge = eventReactive(input$search_single_ghi,{
    if(length(input$cancer_type_ghi)>1){
      validate("Please select one cancer type at a time")
    }
    
    ghi_model(cov=input$hallmark_match,input_gene=input$gene_single_ghi,hallmark_dt = dt_match(),response_match=input$survtype_ghi,geneName(),
              tumor_purity_adj=input$tumor_purity_ghi)
  },ignoreNULL = F)
  
  output$GHI_zscore = renderDT({
    res = t(res_edge())
    DT:::datatable(
      data.frame(res),options = list(searchable =F)
    )    %>%
      formatSignif(columns = c("GHI.part","GHI.full"), digits = 3)
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
  
  #---------------------------------------------------------------------------------------------------------
  #################################### multiple genes GHI ############################################
  #------------------------------------------------------------------------------------------------------
  output$multi_ghi = renderUI({
    ghi_ui_multi("box4",title="Gene-hallmark interaction model (z-score)",cancer_type_ghi="cancer_type_ghi_multi",
                 gene_single_ghi="gene_multi_ghi",survtype_ghi="survtype_ghi_multi",hallmark_match="hallmark_match_multi",tumor_purity_ghi="tumor_purity_ghi_multi",
                 output_dropdown = "GHI_zscore_multi",output_plot = "edge_plot_multi",action_button="search_multi_ghi")
  })
  
  
  ghi_top_mut_multi = reactive({
    res = top_mut()[[input$cancer_type_ghi_multi]]
    res[res==""]=NA
    res = res[complete.cases(res)]
    res
  })
  ghi_top_cna_multi = reactive({
    res = top_cna()[[input$cancer_type_ghi_multi]]
    res[res==""]=NA
    res = res[complete.cases(res)]
    res
  })
  
  
  
  observe({
    tryCatch({
      trigger = input$mainnav_ghi
      sel_input_server(input,output,session,"hallmark_match_multi","Choose hallmarks for interaction models",
                       input_res =  list(
                         "Clinical covariates" = unique(c("age","sex","Tumor_burden")),
                         "Top mutated genes" = paste0(unique(ghi_top_mut_multi()),"_mut"),
                         "Top CNA" = paste0(unique(ghi_top_cna_multi()),"_cna"),
                         "GE" = unique(c("IFN.Hallmark_ssgsea","IFN.Hallmark_mean","ISG.RS_ssgsea",
                                         "ISG.RS_mean","Hypoxia_mean","Hypoxia_ssgsea","BUFFA_HYPOXIA_SCORE")),
                         "Others" = unique(c("CYT","CTL","AR"))
                       ),sel_input =  unique(c("age","sex","CTL","AR")))
    },error=function(e){
      ""
    })
    
  })
  
  observe({
    trigger = input$mainnav_ghi
    sel_input_server_fixed(input,output,session,"gene_multi_ghi","Choose genes (maximum allowance is 5)",input_res=geneName()$geneName,
                           sel_input = unique(c("FGFR2","TBC1D24","PNMA3","TLR1","HCN4")))
  })
  
  # GE data
  dt_match_multi = reactive({
    if(length(input$cancer_type_ghi_multi)>1){
      validate("Please select one cancer type at a time and click 'Run' button")
    }
    
    gene_multi_ghi = lapply(1:length(input$gene_multi_ghi),function(x){
      if(length(which(geneName()$geneName==input$gene_multi_ghi[x]))>1){
        gene_multi_ghi = paste0(make.names(input$gene_multi_ghi[x]),".x")
      }else{
        gene_multi_ghi =  input$gene_multi_ghi[x] 
      }
    })
    
    gene_multi_ghi = unlist(gene_multi_ghi)
    cov_list = unique(c(gene_multi_ghi,input$hallmark_match_multi))
    dt = read_fst(paste0("www/GHI/processed_app/",input$cancer_type_ghi_multi,".fst"),columns=c("Patient.ID","OS.time","OS",
                                                                                                "PFI.time","PFI",cov_list,"Tumor_purity"))
    names(dt) = c("Patient.ID","OS.time","OS",
                  "PFI.time","PFI",input$gene_multi_ghi,input$hallmark_match_multi,"Tumor_purity")
    
    
    dt
  })
  
  ############ edge plot ##############
  res_edge_multi = eventReactive(input$search_multi_ghi,{
    if(length(input$cancer_type_ghi_multi)>1){
      validate("Please select one cancer type at a time")
    }
    
    ghi_model_multi(cov=input$hallmark_match_multi,input_gene=input$gene_multi_ghi,hallmark_dt = dt_match_multi(),response_match=input$survtype_ghi_multi,geneName(),
                    tumor_purity_adj=input$tumor_purity_ghi_multi)
  },ignoreNULL = F)
  
  output$GHI_zscore_multi = renderDT({
    #dt_match_multi()
    res = t(res_edge_multi())
    DT:::datatable(
      data.frame(res),options = list(searchable =F,scroller=T)
    ) %>%
      formatSignif(columns = colnames(res), digits = 3)
  })
  
  output$edge_plot_multi = renderPlot({
    par(mar=c(0.5,1,5,1),cex=0.7)
    chordDiagram( res_edge_multi(), annotationTrack = "grid",  transparency = 0.6,
                  link.visible=abs(res_edge_multi())>1.96,
                  preAllocateTracks = list(track.height = 0.2))
    
    circos.track(track.index = 1, panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, cex=1.5,
                  facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
    }, bg.border = NA) # here set bg.border to NA is important
  })
  
  
  # #--------------------------------------------------------------------------------------
  #   ######################### Common GHI #############################################
  # #--------------------------------------------------------------------------------------
  #   ghi_top_mut_common = reactive({
  #     res = select(top_mut(),input$cancer_type_ghi_common)
  #     res = data.frame(res)
  #     res[res==""]=NA
  #     res =apply(res,1,as.list)
  #     res = unique(unlist(res,use.names = F))
  #     res = res[complete.cases(res)]
  #     res
  #   })
  #   ghi_top_cna_common = reactive({
  #     res = select(top_cna(),input$cancer_type_ghi_common)
  #     res = data.frame(res)
  #     res[res==""]=NA
  #     res =apply(res,1,as.list)
  #     res = unique(unlist(res,use.names = F))
  #     res = res[complete.cases(res)]
  #     res
  #   })
  #   
  #   observe({
  #     tryCatch({
  #       updateSelectizeInput(session,"hallmark_match_common","Choose hallmarks for interaction models",
  #                        choices  =  list(
  #                          "Clinical covariates" = unique(c("age","sex","Tumor_burden")),
  #                          "Top mutated genes" = paste0(unique(ghi_top_mut_common()),"_mut"),
  #                          "Top CNA" = paste0(unique(ghi_top_cna_common()),"_cna"),
  #                          "GE" = unique(c("IFN.Hallmark_ssgsea","IFN.Hallmark_mean","ISG.RS_ssgsea",
  #                                          "ISG.RS_mean","Hypoxia_mean","Hypoxia_ssgsea","BUFFA_HYPOXIA_SCORE")),
  #                          "Others" = unique(c("CYT","CTL","AR"))
  #                        ),selected =  unique(c("age","Tumor_burden","CYT","AR")),server = TRUE)
  #     },error=function(e){
  #       ""
  #     })
  #     
  #   })
  #   
  #   cancer_dt = reactive({
  #     if(input$tumor_purity_ghi_common=="Yes"){
  #       file_folder = "www/GHI/GHI_tumor_purity/"
  #     }else{
  #       file_folder = "www/GHI/GHI_unadj/"
  #     }
  #     lapply(1:length(input$cancer_type_ghi_common),function(x){
  #       read_fst(paste0(file_folder,input$cancer_type_ghi_common[x],".fst"))
  #     })
  #   })
  #   
  # 
  # ghi_common_table = eventReactive(input$ghi_common_run,{
  #   ghi_model_common(input,output,session,cancer_type=input$cancer_type_ghi_common,
  #                    hallmarks=input$hallmark_match_common,top_n = input$top_n_ghi,cancer_dt = cancer_dt())
  # },ignoreNULL = F)
  # 
  # 
  # output$ghi_common_tb = renderDT({
  #   res = data.frame(ghi_common_table())
  #   DT::datatable(
  #     res,options = list(scroller=T,scrollX=TRUE)
  #   )%>%
  #     formatSignif(columns = colnames(res)[-1], digits = 3)
  # })
  
  #========================================================================================================================================
  ############################################ LNCRNA Tool box #################################################################
  #========================================================================================================================================
  #--------------------------------------------------------------------------
  ############ Top prognostic genes across cancers ###################
  #--------------------------------------------------------------------------
  output$common_genes = renderUI({
    common_ui(gene_type = "lnc_mrna_a",cancer_type="cancer_surv_multiple",OS_PFI = "OS_PFI_a",adj_cov="surv_cov_a",
              rank_order = "rank_z_a",top_n = "Top_n_a",font_size = "font_size",action_button = "search_top_across",
              cancer_sel = cancer_type(),output_plot = "heat_top_across")
  })
  
  dt_top_cb = eventReactive(input$search_top_across,{
    common_dt_server(input,output,session,geneName=geneName())  
  })
  
  
  output$heat_top_across = renderPlot({
    data_wide_plot = dt_top_cb()
    col_fun = colorRamp2(c(-5, 0,2,3,4,5), c("#08589E", "white","#FF2400", "#CA3433","#8D021F","#420C09"))
    Heatmap(data_wide_plot,cluster_rows =F,cluster_columns = F,rect_gp = gpar(col = "white", lwd = 2), col = col_fun,column_names_rot = 45,
            heatmap_legend_param=list(title = "Z score"),column_names_gp = grid::gpar(fontsize = input$font_size),
            row_names_gp = grid::gpar(fontsize = input$font_size),cell_fun = function(j, i, x, y, width, height, fill) {
              grid.text(sprintf("%.2f", data_wide_plot[i, j]), x, y, gp = gpar(fontsize = input$font_size))
            })
  })
  
  #------------------------------------------------------------------------
  ################ TIDE interaction #############################
  #-----------------------------------------------------------------------
  ########## TIDE heatmap and pvalue table ##############
  output$tide_res = renderUI({
    tide_ui(gene_type="data_type_interact",cancer_type="cancer_type_interact",z_score="z_score_interact",
            top_n="top_n_interact",action_button="search_tide_heat",cancer_options=cancer_type(),plot_output="tide_heat",table_ouput="tide_table")
  })
  
  dt_tide = eventReactive(input$search_tide_heat,{
    tide_dt(input,output,session)
  })
  
  dt_heat_tide = reactive({
    dt_tide_heat(input,output,session,dt_tide=dt_tide(),geneName=geneName())
    
  })
  
  dt_tide_table = reactive({
    dt_tide_pval(input,output,session,dt_tide=dt_tide(),geneName=geneName())
    
  })
  
  output$tide_table = renderDT(
    # DT:::datatable(
    #   dt_tide_table(),options = list(scrollX=TRUE,scrollY="350px"
    #   )
    # )
    # 
    DT:::datatable(
      dt_tide_table(),options = list(dom='t',scrollX=TRUE,scrollY="350px",
                                     initComplete = JS(
                                       "function(settings, json) {",
                                       "$(this.api().table().header()).css({'background-color': '#4B4B4B', 'color': '#fff'});",
                                       "}")
      )
    )
    
  )
  
  output$tide_heat = renderPlot({
    
    dt_heat_tide = dt_heat_tide()
    
    col_fun = colorRamp2(c(-5, 0, 5), c("#08589E", "white", "#D7301F"))
    
    Heatmap(dt_heat_tide,show_row_dend=F,show_column_dend = F,rect_gp = gpar(col = "white", lwd = 2), col = col_fun)
  })   
  
  ########### TIDE KM plot ##############
  output$tide_km = renderUI({
    tide_km_ui("gene_tide","cutoff_km_tide","cutoff_km_cyt","cancer_type_km_tide",cancer_type(),"KM_tide_hi","KM_tide_lo")
  })
  
  observe({
    trigger = input$mainnav_tide
    update_ui_select(session,"gene_tide",geneID =  as.character(rownames(dt_tide_table())))
  })
  
  tide_km = reactive({
    tide_km_data(input,output,session)
  })
  
  output$KM_tide_hi = renderPlot({
    tryCatch({
      dt = tide_km_fit(input,output,session,tide_km(),hi_lo="high")
      fit = survfit(Surv(X1,X2)~X3,data=dt)
      p = ggsurvplot(fit,data = dt,pval = T,legend.labs=c("high","low"),palette = c("#DF8F44FF","#374E55FF"))
      p$plot+ggtitle(paste0(input$gene_tide,"-high"))
    },error = function(e){
      ggplot() + theme_void()+ggtitle("Not able to generate the figure")
    })
    
  })
  
  output$KM_tide_lo = renderPlot({
    tryCatch({
      dt = tide_km_fit(input,output,session,tide_km(),hi_lo="low")
      fit = survfit(Surv(X1,X2)~X3,data=dt)
      p = ggsurvplot(fit,data = dt,pval = T,legend.labs=c("high","low"),palette = c("#DF8F44FF","#374E55FF"))
      p$plot+ggtitle(paste0(input$gene_tide,"-low"))
    },error = function(e){
      ggplot() + theme_void()+ggtitle("Not able to generate the figure")
    })
    
  })
  #---------------------------------------------------------------------------------
  ################# correlation between mRNA and lncRNA ###########################
  #---------------------------------------------------------------------------------
  observe({
    updateSelectizeInput(session,"input_gene_cor","Select a lncRNA",choices=geneName_lnc()$geneName,selected = "HOTAIR",server=T)
  })

  output$corr_genes = renderUI({
    corr_ui("input_gene_cor","cancer_surv_cor",cancer_type(),"search_cor","most_cor_genes","most_cor_genes_fig")
  })


  dt_cor = eventReactive(input$search_cor,{
    dt = read_fst(paste0("www/correlation/lnc_mrna/",make.names(input$cancer_surv_cor),"_cor.fst"),columns =c("row.names",make.names(input$input_gene_cor)))
    names(dt)[2]="correlation"
    dt$correlation = round(dt$correlation,digits = 4)
    dt = dt %>%
      arrange(desc(abs(correlation)))
    dt
  })

  output$most_cor_genes = renderDT(
    DT:::datatable(
      dt_cor(),options = list(scrollX=TRUE,scrollY="350px",
                              initComplete = JS(
                                "function(settings, json) {",
                                "$(this.api().table().header()).css({'background-color': '#4B4B4B', 'color': '#fff'});",
                                "}")
      )
    )
  )

  dt_cor_plot = eventReactive(input$search_cor,{

    dt = read_fst(paste0("www/correlation/lnc_mrna/",make.names(input$cancer_surv_cor),"_cor.fst"),columns =c("row.names",make.names(input$input_gene_cor)))
    names(dt)[2]="correlation"
    dt$correlation = round(dt$correlation,digits = 4)
    dt= dt[order(abs(dt$correlation),decreasing = T),]
    dt_plot = dt[1:100,]
    dt_plot = rbind(1,dt_plot)

    rownames_dt_plot =dt_plot$row.names;dt_plot=data.frame(dt_plot[,-1]);rownames(dt_plot) = rownames_dt_plot
    dt_plot = as.matrix(dt_plot)
    output <- dt_plot %*% t(dt_plot)
    diag(output) <- 1
    output[-1,-1] = 0
    output[1,1] = 0
    rownames(output)[1] = input$input_gene_cor
    colnames(output)[1] = input$input_gene_cor
    output

  })
  output$most_cor_genes_fig = renderPlot({
    output_plt = dt_cor_plot()
    colors = c("#374E55FF",rep("#DF8F44FF",(nrow(output_plt)-1)))
    qgraph(output_plt,layout="spring",vsize=5,label.prop=0.9,label.cex=0.9,labels=colnames(output_plt),borders=F,posCol="#003399",esize=10,color = colors)
  }
  )
  #-------------------------------------------------------------------------------
  ############### splicing ###############
  #-------------------------------------------------------------------------------
  observe({
 #   tryCatch({
      trigger = input$mainnav_splice
      sel_input_server(input,output,session,"input_gene_splice","Choose a gene",input_res=geneName_mrna()$geneName,sel_input =input$input_gene)
      
    # },error=function(e){
    #   NULL
    # })
  })
    
  dt_splice_plt = reactive({
  #  tryCatch({
      validate(
        need(input$input_gene_splice, "Please select a gene to display the data")
      )
      file_splice = "M:/dept/Dept_BBSR/Projects/Wang_Xuefeng/CGPA/move_022924/CGPA_single/www/Splicing/" #### Need to change !!!
      if(input$database_splice=="SpliceSeq"){
        con <- dbConnect(RSQLite::SQLite(), paste0(file_splice,"spliceseq_full.db"))
        selected_genes <- input$input_gene_splice  # replace with your actual genes of interest
        result <- dbGetQuery(con, sprintf("SELECT rowid FROM GeneSymbol WHERE gene_symbol = ('%s')",selected_genes))
        
        id_list_str <- paste0("(", paste(result$rowid, collapse = ","), ")")
        query_str <- paste0("SELECT * FROM spliceseq WHERE rowid IN ", id_list_str)
        dt =  dbGetQuery(con, query_str)
        dbDisconnect(con) 
        
      }else{
        con <- dbConnect(RSQLite::SQLite(), paste0(file_splice,"spladder_all.db"))
        selected_genes <- input$input_gene_splice  # replace with your actual genes of interest
        result <- dbGetQuery(con, sprintf("SELECT rowid FROM GeneSymbol WHERE gene_symbol = ('%s')",selected_genes))
        
        id_list_str <- paste0("(", paste(result$rowid, collapse = ","), ")")
        query_str <- paste0("SELECT * FROM spladder WHERE rowid IN ", id_list_str)
        dt =  dbGetQuery(con, query_str)
        dbDisconnect(con) 
        
      }

      dt
    # },error=function(e){
    #   ""
    # })

  })
  
  dt_splice_tb = reactive({
    dt = data.frame(dt_splice_plt())
    if(input$database_splice=="SpliceSeq"){
      dt$Splice_Type = paste0(dt$Splice_Type," (",dt$Event,")")
      dt_sub = dt %>% 
        select(Splice_Event,Gene_Symbol,Chr_Strand,Start,End,Splice_Type,Upstream_Exon,
               Downstream_Exon,Alt_Exons,AltExons_IsoName,SpliceIn_IsoName,SpliceOut_IsoName,Splice_Novel)      
    }else{
      dt_sub = dt %>% 
        select(Splice_Event,Gene_Symbol,Chr_Strand,Alt_Region,Splice_Type,AltRegion_IsoName,SpliceIn_IsoName,SpliceOut_IsoName,Splice_Novel )  
    }

    dt_sub
  })
  
  splice_plt = reactive({
   # tryCatch({
      validate(
        need(input$input_gene_splice, "Please select a gene to display the data")
      )
      
      dt = data.frame(dt_splice_plt())
      if(input$surv_type_splice=="OS"&input$cutoffs_splice=="Median"){
        dt_gene=dt%>%select(Splice_Event,Splice_Type,HRmedOS,pvalHRmedOS,cancer) 
      }else if(input$surv_type_splice=="OS"&input$cutoffs_splice=="Optimal") {
        dt_gene=dt%>%select(Splice_Event,Splice_Type,HRfitOS,pvalHRfitOS,cancer) 
      }else if (input$surv_type_splice=="PFI"&input$cutoffs_splice=="Optimal"){
        dt_gene=dt%>%select(Splice_Event,Splice_Type,HRfitPFI,pvalHRfitPFI,cancer) 
      }else if(input$surv_type_splice=="PFI"&input$cutoffs_splice=="Median"){
        dt_gene=dt%>%select(Splice_Event,Splice_Type,HRmedPFI,pvalHRmedPFI,cancer) 
      }
      
      names(dt_gene)[3]="HR"
      names(dt_gene)[4]="pval"
      dt_gene=dt_gene%>%filter(!(is.na(pval)))
      dt_gene$HR=ifelse(dt_gene$HR>1,"Poor prognosis","Favorable prognosis")
      dt_gene$HR=replace(dt_gene$HR, dt_gene$pval>0.05,"NS (p > 0.05)")
      dt_gene$HR=factor(dt_gene$HR,levels=c("Favorable prognosis","Poor prognosis","NS (p > 0.05)"))
      
      dt_gene=dt_gene%>%mutate(size_label=case_when(pval<0.05 ~ "-log10(0.05)",
                                                    pval<0.01 ~ "-log10(0.01)",
                                                    pval>=0.05 ~ "NS"))
      
      
      # Calculate dynamic text size
     # num_rows <- length(unique(dt_gene$Splice_Event))
    #  text_size <- max(3, min(12, 30 / sqrt(num_rows)))  # Adjust '30' or other numbers as needed for better scaling
      
      
      p=dt_gene %>%
        ggplot(aes(Splice_Event, cancer, fill = HR)) +
        geom_point(shape = 20,aes(color=HR,size=-log10(pval))) +
        scale_color_manual(values=c("Poor prognosis"="red", "Favorable prognosis"="blue","NS (p > 0.05)"="grey")) + 
        scale_size_continuous(name= paste0("-log10 (",input$surv_type_splice," P-value)"),
                              range  = c(1, 15), 
                              breaks = c(-log10(0.5), -log10(0.1), -log10(0.05),-log10(0.01)),
                              labels=c("-log10(0.5)","-log10(0.1)","-log10(0.05)","-log10(0.01)")
        )+
        labs(x = NULL, y = NULL) +
        scale_y_discrete(position = "right")+
        theme_classic()+
        theme(axis.text.x = element_text(angle = 60, vjust = 0.8, hjust=0.1,size=12),
              axis.text.y = element_text(size=input$text_size_splice),
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 14),
              plot.title = element_text(size = 20)
              )+
        coord_flip()+theme(aspect.ratio=0.4) +
        ggtitle(paste0(input$database_splice,"_",input$input_gene_splice,"_",input$surv_type_splice,"_",input$cutoffs_splice))
      p
    # },error=function(e){
    #   plot(1, type="n", xlab="X-axis", ylab="Y-axis", main="")
    # })
   
  })
  
  output$splice_plot = renderPlot({
    tryCatch({
      p <- splice_plt()
      p
    },error = function(e) {
      plot(1, type = "n", xlab = "X-axis", ylab = "Y-axis", main = "Error generating plot")
    })
  })
  
  # output$splice_plot = renderPlot({
  #   tryCatch({
  #     # Retrieve the plot from the reactive expression
  #     p <- splice_plt()
  #     
  #     # Extract the number of items directly from the reactive context if needed
  #     # or set a default value if it's static or not dependent on user input.
  #     num_items <- length(unique(p$data$Splice_Event))  # Make sure this is fetching the updated data
  #     
  #     # Calculate dynamic heights based on items
  #     item_height <- 30  # You might need to adjust this based on your actual item size
  #     base_height <- 120  # Additional height for non-item elements (like titles, labels, etc.)
  #     total_height <- num_items * item_height + base_height
  #     
  #     # Set the height dynamically for the plot
  #     # The height value passed here will override the static height set in the UI
  #     p
  #     
  #   }, error = function(e) {
  #     plot(1, type = "n", xlab = "X-axis", ylab = "Y-axis", main = "Error generating plot")
  #   })
  #   
  # }, height = function() { 
  #   # Use the dynamic height calculation here as well
  #   num_items <- length(unique(p$data$Splice_Event))  # Fetch the updated number of items
  #   item_height <- 30
  #   base_height <- 120
  #   total_height <- num_items * item_height + base_height
  #   
  #   # Return the dynamic height for the plot
  #   total_height
  # })
  
  observe({
      output[["download_splice"]]<-downloadHandler(
        filename = function() {
          paste0(input$database_splice,"_",input$input_gene_splice,"_",input$surv_type_splice,"_",input$cutoffs_splice,'.png', sep = '')
          
        },
        content = function(file){
          ggsave(file, splice_plt(),type="cairo-png",width = 30, height = 25,units="cm")
        }
      )
  })
  
  output$splice_dt = renderDT({
    
   dt = data.frame(dt_splice_tb())
   dt = distinct(dt, Splice_Event, .keep_all = TRUE)
   dt = dt[order(dt$Splice_Event,decreasing = T),]
   
   
   DT:::datatable(
     dt,
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
       tags$thead(tags$tr(lapply(colnames(dt), tags$th)))

     )
   )
   
  })
  
  observe({
      output[["download_splice_tb"]]<-downloadHandler(
        filename = function() {
          paste0(input$database_splice,"_",input$input_gene_splice,'.csv', sep = '')
          
        },
        content = function(file){
          fwrite(dt_splice_tb(),file,row.names = T)
        }
      )
    
  })
  
  #-------------------------------------------------------------------------------
  ############### user mannual ###############
  #-------------------------------------------------------------------------------

  output$tab_instructions<-renderUI({
    list(
      h3("Pan-Cancer Summary"),
      div(img(align="center",src=paste0("images/","pan-cancer_summary.jpg"),height="400px"), style="text-align: center;"),
      br(),
      h3("Multivariable Analysis"),
      div(img(align="center",src=paste0("images/","multi1.jpg"),height="620px"), style="text-align: center;"),
      br(),
      h3("Top Prognostic Genes"),
      div(img(align="center",src=paste0("images/","multi2.jpg"),height="600px"), style="text-align: center;"),
      br(),
      h3("Gene-Hallmark Interaction"),
      div(img(align="center",src=paste0("images/","ghi.jpg"),height="580px"), style="text-align: center;"),
      br(),
      h3("ProgSplicing"),
      div(img(align="center",src=paste0("images/","splicing.jpg"),height="375px"), style="text-align: center;"),
      br(),
      h3("LncRNA Exploration - Top Prognostic Genes Across Cancers"),
      div(img(align="center",src=paste0("images/","lnc_top_prog_genes.jpg"),height="265px"), style="text-align: center;"),
      br(),
      h3("LncRNA Exploration - TIDE"),
      div(img(align="center",src=paste0("images/","TIDE.jpg"),height="650px"), style="text-align: center;"),
      br(),
      h3("LncRNA Exploration - Coexpression Analysis"),
      div(img(align="center",src=paste0("images/","coexpression.jpg"),height="300px"), style="text-align: center;"),
      br()

      
      
    )
  }) 
}
shinyApp(ui, server)

