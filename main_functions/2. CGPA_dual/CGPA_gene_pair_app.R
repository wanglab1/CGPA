########################### Code for CGPA-dual ####################################
#------------------------------------------------------------------------------
############ CGPA dual ###################
#------------------------------------------------------------------------------
packages = c("shiny","bslib","shinyWidgets","DT","shinycustomloader","shinyBS","fst","plotly","broom","expss","forestmodel",
             "data.table","shinydashboard","shinydashboardPlus","survival","survminer","MoffittFunctions","tidyverse"
) # rstudioapi, not needed
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

folder_main = "/CGPA_dual"
setwd(folder_main)
source("functions.R")

dark = bs_theme(version = 3,bg="black",fg = "white",warning = "#FFD300" )%>%
  bs_add_rules(sass::sass_file("www/style/style.scss"))

ui <- fluidPage(
  
  theme = dark,
  useShinydashboardPlus(),
  tags$head(
    tags$style(HTML("
        /* Styles for the tabbox */
        .nav-tabs-custom {
          border: none !important;
        }
        .nav-tabs-custom > .nav-tabs {
          background: #242424;
        }
        .nav-tabs-custom > .nav-tabs > li > a {
          color: white;
        }
        .nav-tabs-custom > .nav-tabs > li.active > a {
          background: #b58900;
          color: white;
        }
        
        /* Styles for the box and content */
        .nav-tabs-custom > .tab-content {
          background: #242424;
          color: white;
        }

      "))
  ),
  div(class="navbar1",
      navbarPage(title = "CANCER GENE PROGNOSIS ATLAS",
                 fluid = TRUE, 
                 collapsible = TRUE,
                 
                 tabPanel("CGPA-DUAL",
                          br(),
                          fluidRow(
                            column(2,
                                   br(),br(),
                                   selectizeInput("cancer_type","Cancer type",
                                                  choices=    c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA",
                                                                "GBM","HNSC","KICH","KIRC","KIRP","LGG","LIHC",
                                                                "LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ",
                                                                "SARC","SKCM","STAD","TGCT","THCA","UCEC","UCS","UVM"),multiple=F,selected="HNSC"),
                                   sel_input_ui("survtype","OS or PFI"),
                                   selectizeInput("gene_pair","Select two genes",choices=NULL,multiple=T,options = list(maxItems = 2)),
                                   sel_input_ui("KM_cutoff","KM cutoff",choices=c("Quartile","Optimal","Median"),multiple=F,selected = "Median")
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
                 tabPanel("Help",
                          mainPanel(
                            column(2),
                            column(10,
                                   br(),br(),
                                   textOutput("method_summary"),
                                   br(),br(),
                                   uiOutput("tab_instructions")
                                   )
                        
                          )
                          )
      ))
)
server = function(input,output,session){
  query = reactive({
    query = parseQueryString(session$clientData$url_search)
    query
  })
  

  
  
  geneName = reactive({
    read_fst("www/geneName.fst")
  })
  geneName_mrna = reactive({
    read_fst("www/geneName_mrna.fst")
  })
  geneName_lnc = reactive({
    read_fst("www/geneName_lnc.fst")
  })
  
  observe({
    sel_input_server(input,output,session,"survtype","OS or PFI",c("OS","PFI"),"OS")
  })
  
  
  observe({
    query=query()
    #updateSelectizeInput(session,"geneset",choices = c("None",genesets_pca$genesets),selected  = "None")
    
    if (!is.null(query[['geneid']])) {
      genes = strsplit(query[['geneid']],",")[[1]]
      updateSelectizeInput(session,"gene_pair",choices = geneName()[,1],selected  = genes,server = T)
      
    }else{
      updateSelectizeInput(session,"gene_pair",choices = geneName()[,1],selected = c("CD3E","HAVCR2"),server = T)
      
    }
  })
  
  ############# change on server ###############
  folder_direct = reactive({
  # "/data/CGPA/TCGA_prognostic_app/CGPA_single/"
    "C:/Users/4467777/Desktop/TCGA_shiny/cgpa_single031222/"
    
  })
  ##############################################

#--------------------------------------------------------------------------------
  ###################### Univariate cox model ############################
#-------------------------------------------------------------------------------
  raw_dt = reactive({
    req(input$gene_pair, input$cancer_type)
    dt = lapply(1:2,function(x){
      if(length(which(geneName()$geneName==input$gene_pair[x]))>1){
        dt = read_fst(paste0(folder_direct(),"server_result/surv_data/",input$cancer_type,"_","raw",".fst"),
                      columns = c("OS","OS.time","PFI","PFI.time",paste0(make.names(input$gene_pair[x]),".x"))) 
        names(dt)[ncol(dt)] = make.names(input$gene_pair[x])
      }else{
        dt = read_fst(paste0(folder_direct(),"server_result/surv_data/",input$cancer_type,"_","raw",".fst"),
                      columns = c("OS","OS.time","PFI","PFI.time",make.names(input$gene_pair[x])))    
        names(dt)[ncol(dt)] = make.names(input$gene_pair[x])
      }
      dt
    })
    
    ###################### calculate ratio  ###############################
    dt = data.frame(do.call(cbind,dt))

    dt$ratio = dt[,make.names(input$gene_pair[1])]/dt[,make.names(input$gene_pair[2])]
    dt$ratio  = log2( dt$ratio +1)
    dt[,make.names(input$gene_pair)] = log2(dt[,make.names(input$gene_pair)]+1)
    dt$OS.time = dt$OS.time/30.417
    dt$PFI.time= dt$PFI.time/30.417
    dt = apply_labels(dt,ratio = paste0((input$gene_pair[1]),"/",(input$gene_pair[2])))
    
    
    ###################### calculate ratio group ###############################
    if(input$KM_cutoff=="Median"){
      dt$ratio_group =  findInterval(dt$ratio,median(dt$ratio,na.rm=T))   
      dt$ratio_group = ifelse(dt$ratio_group==1,"high","low")
      
    } else if(input$KM_cutoff=="Quartile"){
      dt$ratio_group =  findInterval(dt$ratio,quantile(dt$ratio,na.rm=T))   
      dt$ratio_group = ifelse(dt$ratio_group%in%c(2,3),NA,ifelse(dt$ratio_group%in%c(4,5),"high","low"))
      
    } else if(input$KM_cutoff=="Optimal"){
      res.cut <- surv_cutpoint(dt,time = "OS.time", event = "OS", variables="ratio")$cutpoint[[1]]
      dt$ratio_group  = ifelse(dt$ratio >=res.cut ,"high","low") 
    }
    
    dt = apply_labels(dt,ratio_group = paste0((input$gene_pair[1]),"/",(input$gene_pair[2])))
    
    dt
  })
  

  
  uni_cox_fig1 = reactive({
    
    req(input$gene_pair, input$survtype)
    
   if(is.null(input$run_dual)){
     uni_cox = lapply(1:2,function(x){
       formula =  as.formula(paste('Surv(OS.time, OS)~',make.names(input$gene_pair[x])))
       fit = coxph(formula,data=raw_dt())
       
       
       
       forest_model(fit)
     })
     
     
   }else{
     if(input$survtype=="OS"){
       uni_cox = lapply(1:2,function(x){
         formula =  as.formula(paste('Surv(OS.time, OS)~',make.names(input$gene_pair[x])))
         fit = coxph(formula,data=raw_dt())
         
         
         
         forest_model(fit)
       })
       
     }else if (input$survtype=="PFI"){
       uni_cox = lapply(1:2,function(x){
         formula =  as.formula(paste('Surv(PFI.time, PFI)~',make.names(input$gene_pair[x])))
         fit = coxph(formula,data=raw_dt())
         forest_model(fit)
         
       })
       
     }
   }
    

    
  })

  
output$uni_cox_plot1 = renderPlot({
  tryCatch({
    uni_cox_fig1()[[1]]
    
  },error=function(e){
    ""
  })
})
       
output$uni_cox_plot2 = renderPlot({
  tryCatch({
    uni_cox_fig1()[[2]]
  },error= function(e){
    ""
  })
})

#--------------------------------------------------------------------------------
###################### Univariate cox KM ############################
#-------------------------------------------------------------------------------            
cutoff_dt = reactive({
  dt = lapply(1:2,function(x){
    if(length(which(geneName()$geneName==input$gene_pair[x]))>1){
      dt = read_fst(paste0(folder_direct(),"server_result/surv_data/",input$cancer_type,"_",tolower(input$KM_cutoff),".fst"),
                    columns = c("OS","OS.time","PFI","PFI.time",paste0(make.names(input$gene_pair[x]),".x"))) 
      names(dt)[ncol(dt)] = make.names(input$gene_pair[x])
    }else{
      dt = read_fst(paste0(folder_direct(),"server_result/surv_data/",input$cancer_type,"_",tolower(input$KM_cutoff),".fst"),
                    columns = c("OS","OS.time","PFI","PFI.time",make.names(input$gene_pair[x])))    
      names(dt)[ncol(dt)] = make.names(input$gene_pair[x])
    }
    dt
  })
  
  dt = data.frame(do.call(cbind,dt))
  dt$OS.time = dt$OS.time/30.417
  dt$PFI.time= dt$PFI.time/30.417
  dt$gene_combine = paste0(dt[[make.names(input$gene_pair[1])]],"+",dt[[make.names(input$gene_pair[2])]])
  dt
})    
  
# KM plot
uni_km_plot = reactive({
  req(input$gene_pair, input$survtype)
  
  if(input$survtype=="OS"){
    km = lapply(1:2,function(x){
      fit <- eval(parse(text = paste0("survfit(Surv(OS.time,OS) ~ ", make.names(input$gene_pair[x]), ", data = cutoff_dt())")))
      p = ggsurvplot(fit,data = cutoff_dt(),pval=T,legend.title=input$gene_pair[x],legend.labs=c("high","low"),palette = c("#DF8F44FF","#374E55FF"),size=1)
      p
    })
  }else{
    km = lapply(1:2,function(x){
      fit <- eval(parse(text = paste0("survfit(Surv(PFI.time,PFI) ~ ", make.names(input$gene_pair[x]), ", data = cutoff_dt())")))
      p = ggsurvplot(fit,data = cutoff_dt(),pval=T,legend.title=input$gene_pair[x],legend.labs=c("high","low"),palette = c("#DF8F44FF","#374E55FF"),size=1)
      p
    })
  }

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
    if(input$survtype=="OS"){
      time = "OS.time";status = "OS"
    }else{
      time = "PFI.time";status="PFI"
    }
    OS_KM_fit_table <- purrr::map_dfr( make.names(input$gene_pair), run_pretty_km_output, model_data = cutoff_dt(), time_in = time, 
                                       event_in = status, event_level = '1') %>% 
      select(Group, Level, everything())
    data.frame(OS_KM_fit_table)
  },error = function(e){
    data.table()
  })

})

output$uni_km_table = renderDT({
  tryCatch({
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
  if(input$survtype=="OS"){
    formula =  as.formula(paste('Surv(OS.time, OS)~ratio'))
  }else{
    formula =  as.formula(paste('Surv(PFI.time, PFI)~ratio'))
  }
  fit = coxph(formula,data= raw_dt())
})

output$cox_ratio = renderPlot({
  tryCatch({
    forest_model(ratio_cox())
    
  },error = function(e){
    ""
  })
})
#------------------------------------------------------------------------------
#################### Gene ratio km ######################################
#------------------------------------------------------------------------------
ratio_km = reactive({
 # eventReactive(input$run_dual,{
  if(input$survtype=="OS"){
    formula =  as.formula(paste('Surv(OS.time, OS)~ratio_group'))
  }else{
    formula =  as.formula(paste('Surv(PFI.time, PFI)~ratio_group'))
  }
  fit = survfit(formula,data = raw_dt())
  fit$call$formula <- formula
  
  fit
})

output$km_ratio = renderPlot({
  tryCatch({
    p = ggsurvplot(ratio_km(),data = raw_dt(),
                   pval=T,legend.title=paste0(input$gene_pair[1],"/",input$gene_pair[2]),
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
    if(input$survtype=="OS"){
      time = "OS.time";status = "OS"
    }else{
      time = "PFI.time";status="PFI"
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
  if(input$survtype=="OS"){
    formula = as.formula(paste0("Surv(OS.time,OS)~",paste(make.names(input$gene_pair),collapse = "+")))
  }else{
    formula = as.formula(paste0("Surv(PFI.time,PFI)~",paste(make.names(input$gene_pair),collapse = "+")))
  }
  fit = coxph(formula,data = raw_dt())
  fit
})

output$multi_cox = renderPlot({
  tryCatch({
    forest_model(multi_cox())
    
  },error=function(e){
    ""
  })
})

#-------------------------------------------------------------------------------
################ KM with 4 groups ########################
#-------------------------------------------------------------------------------
########## KM plot ###########

cutoff_dt_multi = reactive({
  na.omit(cutoff_dt())
})

multi_km_plot = reactive({
  #eventReactive(input$run_dual,{
  if(input$survtype=="OS"){
    formula = as.formula("Surv(OS.time,OS)~gene_combine")
  }else{
    formula = as.formula("Surv(PFI.time,PFI)~gene_combine")
  }
  fit = survfit(formula,data = cutoff_dt_multi())
  fit$call$formula <- formula
  if(length(unique(cutoff_dt_multi()$gene_combine))==2){
    p = ggsurvplot(fit,data = cutoff_dt_multi(),pval=T,legend.title=paste0(input$gene_pair[1],"+",input$gene_pair[2]),
    legend.labs=c("high+high","low+low"),palette = c("#DF8F44FF","#374E55FF"))
  }else{
    p = ggsurvplot(fit,data = cutoff_dt_multi(),pval=T,legend.title=paste0(input$gene_pair[1],"+",input$gene_pair[2]),
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
    #eventReactive(input$run_dual,{
    if(input$survtype=="OS"){
      time = "OS.time";status = "OS"
    }else{
      time = "PFI.time";status="PFI"
    }
    OS_KM_fit_table <- purrr::map_dfr( "gene_combine", run_pretty_km_output, model_data = cutoff_dt_multi(), time_in = time,
                                       event_in = status, event_level = '1') %>%
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
    if(input$survtype=="OS"){
      formula = as.formula("Surv(OS.time,OS)~gene_combine")
    }else{
      formula = as.formula("Surv(PFI.time,PFI)~gene_combine")
    }
    res = pairwise_survdiff(formula,data=cutoff_dt_multi())
    res = data.frame(broom::tidy(res))
    res$p.value = ifelse(res$p.value>=0.001,round(res$p.value,3),
                         formatC(res$p.value,format = "e",digits=3))
    names(res)[1] = input$gene_pair[1]
    names(res)[2] = input$gene_pair[2]
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
    if(input$survtype=="OS"){
      formula = as.formula(paste0("Surv(OS.time,OS)~",paste(make.names(input$gene_pair),collapse = "*")))
    }else{
      formula = as.formula(paste0("Surv(PFI.time,PFI)~",paste(make.names(input$gene_pair),collapse = "*")))
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
km_dt_interaction = reactive({
  if(length(which(geneName()$geneName==input$gene_pair[1]))>1){
    dt = read_fst(paste0(folder_direct(),"server_result/surv_data/",input$cancer_type,"_",tolower(input$KM_cutoff),".fst"),
                  columns = c("OS","OS.time","PFI","PFI.time",paste0(make.names(input$gene_pair[1]),".x"))) 
    names(dt)[ncol(dt)] = make.names(input$gene_pair[1])
  }else{
    dt = read_fst(paste0(folder_direct(),"server_result/surv_data/",input$cancer_type,"_",tolower(input$KM_cutoff),".fst"),
                  columns = c("OS","OS.time","PFI","PFI.time",make.names(input$gene_pair[1])))    
    names(dt)[ncol(dt)] = make.names(input$gene_pair[1])
  }
  dt
  
  dt$OS.time = dt$OS.time/30.417
  dt$PFI.time= dt$PFI.time/30.417
  
  if(length(which(geneName()$geneName==input$gene_pair[2]))>1){
    dt_raw = read_fst(paste0(folder_direct(),"server_result/surv_data/",input$cancer_type,"_",tolower("raw"),".fst"),
                      columns = c(paste0(make.names(input$gene_pair[2]),".x"))) 
    names(dt_raw)[ncol(dt_raw)] = make.names(input$gene_pair[2])
  }else{
    dt_raw = read_fst(paste0(folder_direct(),"server_result/surv_data/",input$cancer_type,"_",tolower("raw"),".fst"),
                      columns = c(make.names(input$gene_pair[2])))    
    names(dt_raw)[ncol(dt_raw)] = make.names(input$gene_pair[2])
  }
  dt_raw
  
  dt = data.frame(cbind(dt,dt_raw))   
})

km_interaction_hi = reactive({
  dt = km_dt_interaction()
  dt1 = subset(dt,dt[[input$gene_pair[1]]]=="high")
  if(input$survtype=="OS"){
    time="OS.time";status="OS"
  }else{
    time="PFI.time";status="PFI"
  }
  if(input$KM_cutoff=="Median"){
    dt1$cut = findInterval(dt1[[input$gene_pair[2]]], median(dt1[[input$gene_pair[2]]],na.rm=T))
    dt1$cut = ifelse(dt1$cut==1,"high","low")
  }else if(input$KM_cutoff=="Optimal"){
    res.cut <- surv_cutpoint(dt1,time = time, event = status, variables=input$gene_pair[2])
    dt1 <- surv_categorize(res.cut)
    names(dt1)[3]="cut"
  }else{
    dt1$cut = findInterval(dt1[[input$gene_pair[2]]], quantile(dt1[[input$gene_pair[2]]],na.rm=T))
    dt1$cut = ifelse(dt1$cut%in%c(2,3),NA,ifelse(dt1$cut%in%c(4,5),"high","low"))
  }
  dt1
})

km_interaction_lo = reactive({
  dt = km_dt_interaction()
  dt1 = subset(dt,dt[[input$gene_pair[1]]]=="low")
  if(input$survtype=="OS"){
    time="OS.time";status="OS"
  }else{
    time="PFI.time";status="PFI"
  }
  
  if(input$KM_cutoff=="Median"){
    dt1$cut = findInterval(dt1[[input$gene_pair[2]]], median(dt1[[input$gene_pair[2]]],na.rm=T))
    dt1$cut = ifelse(dt1$cut==1,"high","low")
  }else if(input$KM_cutoff=="Optimal"){
    res.cut <- surv_cutpoint(dt1,time = time, event = status, variables=input$gene_pair[2])
    dt1 <- surv_categorize(res.cut)
    names(dt1)[3]="cut"
  }else{
    dt1$cut = findInterval(dt1[[input$gene_pair[2]]], quantile(dt1[[input$gene_pair[2]]],na.rm=T))
    dt1$cut = ifelse(dt1$cut%in%c(2,3),NA,ifelse(dt1$cut%in%c(4,5),"high","low"))
  }
  dt1
})

# KM plot 
km_interation_plot = reactive({
  #eventReactive(input$run_dual,{
  if(input$survtype=="OS"){
    formula = as.formula("Surv(OS.time,OS)~cut")
  }else{
    formula = as.formula("Surv(PFI.time,PFI)~cut")
  }
  
  fit1 = survfit(formula,data = km_interaction_hi())
  fit2 = survfit(formula,data = km_interaction_lo())
  fit1$call$formula <- formula
  fit2$call$formula <- formula
  
  p1 = ggsurvplot(fit1,data = km_interaction_hi(),pval=T,legend.title=input$gene_pair[2],legend.labs=c("high","low"),palette = c("#DF8F44FF","#374E55FF")) + ggtitle(paste0(input$gene_pair[1],":High"))
  p2 = ggsurvplot(fit2,data = km_interaction_lo(),pval=T,legend.title=input$gene_pair[2],legend.labs=c("high","low"),palette = c("#DF8F44FF","#374E55FF"))+ ggtitle(paste0(input$gene_pair[1],":Low"))
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
  test_plt = data.frame(test[,c(input$gene_pair[1],input$gene_pair[2])])
  names(test_plt) = c("gene1","gene2")
  for(i in 1:2){
    test_plt[,i] = log2(test_plt[,i]+1)
  }
  
  pearson_corr <- round(cor(test_plt$gene1, test_plt$gene2, method = "pearson"), 2)
  spearman_corr <- round(cor(test_plt$gene1, test_plt$gene2, method = "spearman"), 2)
  
  scatter_plot <- ggscatter(test_plt, x = "gene1", y = "gene2",
                            add = "reg.line",                                 # Add regression line
                            conf.int = TRUE,                                  # Add confidence interval
                            add.params = list(color = "blue",
                                              fill = "lightgray")
  )+
   xlab(input$gene_pair[1])+
   ylab(input$gene_pair[2])
  
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
#==============================================================================
# Instruction page 
#==============================================================================
output$method_summary<-renderText({
  "CGPA Dual is primarily designed for conducting survival analysis when you have two gene inputs. It provides a comprehensive set of tools for analyzing the survival data associated with these genes."
})
output$tab_instructions<-renderUI({
  list(
    h3("Key Features"),
    div(img(align="center",src=paste0("images/","general.jpg"),height="332px"), style="text-align: center;"),
    br(),
    h4("Univariable Model"),
     div(img(align="center",src=paste0("images/","uni_survival.jpg"),height="88px"), style="text-align: center;"),
    br(),
    h4("Two-gene Model"),
     div(img(align="center",src=paste0("images/","multi_surv.jpg"),height="400px"), style="text-align: center;"),
    br(),
    h4("Interaction Model"),
    div(img(align="center",src=paste0("images/","interaction.jpg"),height="395px"), style="text-align: center;"),
    br(),
     h4("Gene Ratio"),
     div(img(align="center",src=paste0("images/","gene_ratio.jpg"),height="80px"), style="text-align: center;")
  

  )
})

}

shinyApp(ui, server)
