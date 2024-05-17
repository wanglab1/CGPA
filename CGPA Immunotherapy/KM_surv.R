#-------------------------------------------------------------------------------
# KM plot
#------------------------------------------------------------------------------
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

        
      )

  ) 
  
  
  
}

folder_main = "C:/Users/4467777/Desktop/TCGA_shiny/CGPA_ICI/ICI_shiny2/"
setwd(folder_main)

#folder_main = "C:/Users/4467777/Desktop/TCGA_shiny/cgpa_single031222/"

survival_server1 <-function(geneName,cutoff,data){
  #-------------------------------------------------------------------------------
  # KM plot
  # data: import dataset
  # geneName: user input geneNames
  # cutoff: cutoffs for KM
  # time: survival time
  # status: survival status
  # data = dt_study;geneName=sel_gene;time="PFS.time";status = "PFS"
  #-------------------------------------------------------------------------------
  geneName = make.names(geneName)
  tryCatch({
    if(cutoff=="optimal"){
      res.cut <- surv_cutpoint(data, time = "time", event = "status", variables = "gene")
      res.cat <- surv_categorize(res.cut)
      fit <- survfit(Surv(time,status)~gene,data = res.cat)
      p = ggsurvplot(fit,data=res.cat,pval = T,legend.title = geneName,legend.labs=c("high","low"),palette =c("#DF8F44FF", "#374E55FF"))
      p = p$plot+
        xlab("Time (month)")
      
    }else if(cutoff=="median"){
      dt = data
      dt$cut = findInterval(dt$gene, median(dt$gene,na.rm=T))
      dt$cut = ifelse(dt$cut==1,"high","low")
      res.cat = data.frame(cbind("time"=dt$time,"status"=dt$status,"cut"=dt$cut))
      res.cat$time = as.numeric(as.character(res.cat$time))
      res.cat$status = as.numeric(as.character(res.cat$status))
      fit = survfit(Surv(time,status) ~cut, data = res.cat)
      tryCatch({
        if(length(unique(na.omit(res.cat$cut)))==1){
          p= ggplot() +
            theme_void() +
            geom_text(aes(0,0,label=paste0("Too many zeros")))+
            xlab("Time (month)")
        }else{
          p = ggsurvplot(fit,data=res.cat,pval = T,legend.title = geneName,legend.labs=c("high","low"),palette = c("#DF8F44FF", "#374E55FF"))
          p = p$plot+
            xlab("Time (month)")
        }
      },error=function(e){
        df <- data.frame()
        p = ggplot(df) + geom_point() + xlim(0, 10) + ylim(0, 100)+theme_classic()+
          ggtitle("Too many zeros")
      })
      
    }else if(cutoff=="quartile"){
      dt = data
      dt$cut = findInterval(dt$gene, quantile(dt$gene,na.rm=T))
      dt$cut = ifelse(dt$cut%in%c(2,3),NA,ifelse(dt$cut%in%c(4,5),"high","low"))
      res.cat = data.frame(cbind("time"=dt$time,"status"=dt$status,"cut"=dt$cut))
      res.cat$time = as.numeric(as.character(res.cat$time))
      res.cat$status = as.numeric(as.character(res.cat$status))
      fit = survfit(Surv(time, status) ~cut, data = res.cat)
      tryCatch({
        if(length(unique(na.omit(res.cat$gene)))==1){
          
          p= ggplot() +
            theme_void() +
            geom_text(aes(0,0,label=paste0("Too many zeros")))+
            xlab("Time (month)")
          
        }else{
          
          p = ggsurvplot(fit,data=res.cat,pval = T,legend.title = geneName,legend.labs=c("high","low"),palette = c("#DF8F44FF", "#374E55FF"))
          p = p$plot+
            xlab("Time (month)")
          
        }
      },error = function(e){
        df <- data.frame()
        p = ggplot(df) + geom_point() + xlim(0, 10) + ylim(0, 100)+theme_classic()+
          ggtitle("Too many zeros")
      })
    }
    
  },error = function(e){
    fit = survfit(Surv(time,status) ~1, data = data)
    p = ggsurvplot(fit,data=data,palette = c("#999999"))
    p= p$plot+
      ggtitle(paste0("Not able to generate the plot"))+
      xlab("Time (month)")
  })
  
  
  return(p)
}



survival_server2 <-function(geneName,input_data){
  #-----------------------------------------------------------------------------
  # Univariable coxph result
  # geneName: user-select genenames
  # study_type: study_type; geneName: input gene
  # geneName = sel_gene;input_data = dt_study
  #----------------------------------------------------------------------------- 
  geneName = make.names(geneName)
  
  tryCatch({
    dt =  input_data
    # dtt = na.omit(data.frame(dt))
    
    KM_fit_table <- run_pretty_model_output(
      "gene",
      model_data = dt,  y_in = "time", event_in = "status", event_level = '1',
      p_digits = 5) 
    KM_fit_table = data.frame(KM_fit_table[,c(-2,-(ncol(KM_fit_table)-1))])
    KM_fit_table$Variable = geneName
    
    return(KM_fit_table)
  },error=function(e){
    ""
  })
  
  
}

