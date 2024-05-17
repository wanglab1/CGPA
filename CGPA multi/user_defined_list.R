folder_main = "C:/Users/4467777/Desktop/TCGA_shiny/cgpa_multi_031322/"

options(shiny.maxRequestSize = 500*1024^2)  # 500MB in bytes

ui = fluidPage(
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
                                 ),
                                 tabPanel("PROGNOSTIC RANKING",
                                          fluidRow(column(9,
                                                          br(),
                                                          fluidRow(       
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
                                                            )),
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
                                                            ))
                                                            
                                                          )
                                          ),
                                          column(3,
                                                 br(),
                                                 uiOutput("gbm_surv_ui_custom")
                                          )
                                          
                                          )
                                 ),
                                 tabPanel("JOINT SIGNATURE",
                                          br(),
                                          fluidRow(
                                            column(12,
                                                   tabsetPanel(type="pills",id="mainnav_surv",
                                                               tabPanel("Average score",
                                                                        fluidRow(
                                                                          column(9,br(),div(id="box2",
                                                                                            box(
                                                                                              title="Average score KM plot",collapsible = TRUE, status="warning",
                                                                                              solidHeader = TRUE,
                                                                                              width = 12,height=500,
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
                                                                          column(3,br(),
                                                                                 uiOutput("average_km_ui_custom")
                                                                          )
                                                                        )
                                                                        
                                                               ),
                                                               tabPanel("ssGSEA",
                                                                        fluidRow(
                                                                          column(
                                                                            9,br(),div(id="box2",
                                                                                       box(
                                                                                         title="ssGSEA KM plot",collapsible = TRUE, status="warning",
                                                                                         solidHeader = TRUE,
                                                                                         width = 12,height=500,
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
                                                                          column(3,br(),
                                                                                 uiOutput("gsea_km_ui_custom")
                                                                          )
                                                                          
                                                                        )
                                                                        
                                                               )
                                                   ),
                                                   
                                            )
                                            
                                          )
                                 ),
                                 tabPanel("SUBNETWORK",
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
                                                   
                                            ),
                                            column(3,
                                                   br(),
                                                   uiOutput("subnet_ui_custom")
                                            )
                                          ))
                     ))
         )
)

)


server = function(input,output,session){
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
    data.table(merged_dt$data[1:500,1:100],options = list(scrollX = TRUE))
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
  
  div(
    selectizeInput("gbm_surv_type_custom","OS or PFI",result,result[1])
  )
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
output$average_km_ui_custom = renderUI({
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
  div(
    selectizeInput("gsea_km_type_av_custom","OS or PFI",result,result[1]),
    selectizeInput("gsea_km_cutoff_av_custom","Choose a cutoff",c("median","optimal","quartile"),"median"),
    actionButton("search_average_joint_custom", "Run")      
  )
  
})

output$gsea_km_ui_custom = renderUI({
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
  
  div(
    selectizeInput("gsea_km_type_custom","OS or PFI",result,result[1]),
    selectizeInput("gsea_km_cutoff_custom","Choose a cutoff",c("median","optimal","quartile"),"median"),
    actionButton("search_gsea_custom", "Run")
  )
  

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

output$gsea_km_plot_custom = renderPlot({
  
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
output$average_km_plot_custom = renderPlot({

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

#-------------------------------------------------------------------------------
###################### EGAnet #####################################
#-------------------------------------------------------------------------------
output$subnet_ui_custom = renderUI({
  div(
    p("Click 'Run' button to get the results"),
    actionButton("search_subnet_custom", "Run")
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
  sel_genes = input$gene_custom
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
shinyApp(ui, server)

