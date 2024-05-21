####################### Tab survvival, KM plot ########################################



survival_UI<-function(headings,output_info,more_cox){
  #column(12,
  div(id="box4",
      box(
        title=headings,collapsible = TRUE, status="warning",
        solidHeader = TRUE,
        width = 12,height= 450,
        plotOutput(output_info,width="100%",height="300px"),
        #withLoader(plotOutput(output_info,width="100%",height="300px"),type="html",loader="loader1"),
        # span(withLoader(plotOutput(output_info,width="100%",height="300px"),type="html",loader="loader1"),style="text-align: center;"),
        br(),
        p("Cox model result"),
        div(
          DTOutput(more_cox,width = 500)
        )
        # div(style="border-color: #b58900;display:inline-block",
        #     actionButton(clickx, "======Cox model result======"), style="display:center-align"),
        # 
        # bsModal(bs_id, "", clickx, size = "large",
        #         # h4("Log-rank test result"),
        #         # DTOutput(more_logrank),
        #         h4("Cox proportional hazard result"),
        #         withLoader(dataTableOutput(more_cox),type="html",loader="loader1")
        #         # h4("Interpretation"),
        #         # textOutput(more_inter)
        # )
        
      )
      #         )
      
  ) 

  
  
}

folder_main = dirname(getSourceEditorContext()$path)
setwd(folder_main)

#folder_main = "C:/Users/4467777/Desktop/TCGA_shiny/cgpa_single031222/"

survival_server1 <-function(input,output,session,geneName,cancer_multi){
  
  
  tryCatch({
  if(length(which(geneName$geneName==input$input_gene_dash))>1){
    dt = read_fst(paste0("server_result/surv_data/",cancer_multi,"_",tolower(input$cutoff_multi),".fst"),
                  columns = c("OS","OS.time","PFI","PFI.time",paste0(make.names(input$input_gene_dash),".x"))) 
    names(dt)[ncol(dt)] = make.names(input$input_gene_dash)
  }else{
    dt = read_fst(paste0("server_result/surv_data/",cancer_multi,"_",tolower(input$cutoff_multi),".fst"),
                  columns = c("OS","OS.time","PFI","PFI.time",make.names(input$input_gene_dash)))  
    names(dt)[ncol(dt)] = make.names(input$input_gene_dash)
  }
  
  
  
  dt$OS.time = dt$OS.time/30.417
  dt$PFI.time= dt$PFI.time/30.417
  
 # dt = na.omit(data.frame(dt))
  
  if(input$survtype=="OS"){
    fit <- eval(parse(text = paste0("survfit(Surv(OS.time,OS) ~ ", make.names(input$input_gene_dash), ", data = dt)")))
  }else if(input$survtype=="PFI"){
    fit <- eval(parse(text = paste0("survfit(Surv(PFI.time,PFI) ~ ", make.names(input$input_gene_dash), ", data = dt)")))
  }
  p = ggsurvplot(fit,data = dt,pval=T,legend.title=input$input_gene_dash,legend.labs=c("high","low"),palette = c("#DF8F44FF","#374E55FF"),size=0.5)
  p$plot
  },error=function(e){
       ggplot() + theme_void()+ggtitle("Not able to generate a figure")
   })
}


survival_server2 <-function(input,output,session,geneName,cancer_multi){
  
  tryCatch({
    
    if(length(which(geneName$geneName==input$input_gene_dash))>1){
      dt = read_fst(paste0("server_result/surv_data/",cancer_multi,"_","raw",".fst"),
                    columns = c("OS","OS.time","PFI","PFI.time",paste0(make.names(input$input_gene_dash),".x"))) 
      names(dt)[ncol(dt)] = make.names(input$input_gene_dash)
    }else{
      dt = read_fst(paste0("server_result/surv_data/",cancer_multi,"_","raw",".fst"),
                    columns = c("OS","OS.time","PFI","PFI.time",make.names(input$input_gene_dash)))    
      names(dt)[ncol(dt)] = make.names(input$input_gene_dash)
    }
    
    
    
    dt$OS.time = dt$OS.time/30.417
    dt$PFI.time= dt$PFI.time/30.417
    dt[make.names(input$input_gene_dash)] = log2(dt[make.names(input$input_gene_dash)]+1)

   # dtt = na.omit(data.frame(dt))
  
    if(input$survtype=="OS"){

      OS_KM_fit_table <- run_pretty_model_output(
        make.names(input$input_gene_dash),
        #"EGFR",
        model_data = dt,  y_in = "OS.time", event_in = "OS", event_level = '1',
        p_digits = 5) 
      OS_KM_fit_table = data.frame(OS_KM_fit_table[,c(-2,-(ncol(OS_KM_fit_table)-1))])
    }else if(input$survtype=="PFI"){
      OS_KM_fit_table <- run_pretty_model_output(
        make.names(input$input_gene_dash),
        model_data = dt,  y_in = "PFI.time", event_in = "PFI", event_level = '1',
        p_digits = 5) 
      OS_KM_fit_table = data.frame(OS_KM_fit_table[,c(-2,-(ncol(OS_KM_fit_table)-1))])
      
    }
    return(OS_KM_fit_table)
  },error=function(e){
    ""
  })

  
}


