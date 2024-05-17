library(MoffittFunctions)
library(AdjKMCIF)
#folder_main =  "/data/CGPA/"
# source("F:/Projects/Kim JF/survFunc/surv_adjust/adjusted_km/adj_km_051321/4_adjusted_KM_final.R")
# source("F:/Projects/Kim JF/survFunc/surv_adjust/adjusted_km/adj_km_051321/1_surv_prob.R")
# source("F:/Projects/Kim JF/survFunc/surv_adjust/adjusted_km/adj_km_051321/2_adj_km.R")

survival_adj_UI<-function(headings,output_info,bs_id,clickx,more_cox,more_inter){
  #column(12,
  div(id="box4",
      box(
        title=headings,collapsible = TRUE, status="warning",
        solidHeader = TRUE,
        width = 12,height= 450,
        withLoader(plotOutput(output_info,width="80%",height="300px"),type="html",loader="loader1"),
       # span(withLoader(plotOutput(output_info,width="100%",height="300px"),type="html",loader="loader1"),style="text-align: center;"),
        br(),
       p("Cox model result"),
       div(
         DTOutput(more_cox,width = 500)
           )
        # div(style="border-color: #b58900;display:inline-block",
        #     actionButton(clickx, "Click to check Cox model result"), style="display:center-align"),
        # 
        # bsModal(bs_id, "", clickx, size = "large",
        #         # h4("Log-rank test result"),
        #         # DTOutput(more_logrank),
        #         h4("Cox proportional hazard result"),
        #         DTOutput(more_cox)
        #         # h4("Interpretation"),
        #         # textOutput(more_inter)
        # )
        
      )

  ) 
  
  
}

survival_adj_data = function(cancer_type,adj_cov,cutoff_adj,input_gene,gene_name_list){
  #-------------------------------------------------------
  # cancer_type: cancer types
  # adj_cov: adjusted covariates
  # cutoff_adj: OS or PFI or raw
  # gene_name_list: geneName() in shiny
  #-------------------------------------------------------
  
  

      if("CTL"%in%adj_cov){
        dt = read_fst(paste0("server_result/surv_data/",cancer_type,"_raw",".fst"),
                      columns = c("CD8A","CD8B","GZMA","GZMB","PRF1"))
        dt = log2(dt+1)
        CTL = apply(dt,1,mean,na.rm=T)
        
        if(length(which(gene_name_list$geneName==input_gene))>1){
          dt = read_fst(paste0("server_result/surv_data/",cancer_type,"_",tolower(cutoff_adj),".fst"),
                        columns = c("OS","OS.time","PFI","PFI.time","age","sex","purity",paste0(make.names(input_gene),".x")))
          names(dt)[ncol(dt)] = make.names(input_gene)

          names(dt)[7] = "tumor_purity"
        }else{
          dt = read_fst(paste0("server_result/surv_data/",cancer_type,"_",tolower(cutoff_adj),".fst"),
                        columns = c("OS","OS.time","PFI","PFI.time","age","sex","purity",make.names(input_gene)))
          names(dt)[ncol(dt)] = make.names(input_gene)

          names(dt)[7] = "tumor_purity"
        }  
        
        dt = data.frame(cbind(CTL,dt))
        
      }else{
        if(length(which(gene_name_list$geneName==input_gene))>1){
          dt = read_fst(paste0("server_result/surv_data/",cancer_type,"_",tolower(cutoff_adj),".fst"),
                        columns = c("OS","OS.time","PFI","PFI.time","age","sex","purity",paste0(make.names(input_gene),".x")))
          names(dt)[ncol(dt)] = make.names(input_gene)

          names(dt)[7] = "tumor_purity"
        }else{
          dt = read_fst(paste0("server_result/surv_data/",cancer_type,"_",tolower(cutoff_adj),".fst"),
                        columns = c("OS","OS.time","PFI","PFI.time","age","sex","purity",make.names(input_gene)))
          names(dt)[ncol(dt)] = make.names(input_gene)

          names(dt)[7] = "tumor_purity"
        }  
        dt = data.frame(dt)
      }
      
      
      dt$OS.time = dt$OS.time/30.417
      dt$PFI.time= dt$PFI.time/30.417
      dt
    #  na.omit(data.frame(dt))
  
}


survival_adj_server1 <-function(input,output,session,inputdata,survtype,cutoff,input_gene_surv,adj_cov){
  
  tryCatch({
    
    if(survtype=="OS"){

      dtt = na.omit(data.frame(inputdata))
      
      res = adjusted_KM(data=dtt,time='OS.time',status="OS",group=input_gene_surv,covlist=adj_cov,stratified_cox = "Yes",reference_group="G&B")


         p = ggplot(res,aes(x=time,y = prob, group =class))+
           geom_step(aes(color = class),size=0.7)+
           theme_classic()+
           ylim(c(0,1))+
           theme(legend.position="top")+
           scale_color_manual(values=c("#DF8F44FF", "#374E55FF"))+
           ylab("Adjusted survival probability")+
           xlab("Time")
         p

    }else if(survtype=="PFI"){
      
          dtt = na.omit(data.frame(inputdata))
          res = adjusted_KM(data=dtt,time='PFI.time',status="PFI",group=input_gene_surv,covlist=adj_cov,stratified_cox = "Yes",reference_group="G&B")
          

          

        
          p =  ggplot(res,aes(x=time,y = prob, group =class))+
            geom_step(aes(color = class),size=0.7)+
            theme_classic()+
            ylim(c(0,1))+
            theme(legend.position="top")+
            scale_color_manual(values=c("#DF8F44FF", "#374E55FF"))+
            ylab("Adjusted survival probability")+
            xlab("Time")
          p

      
    }
    
   },error=function(e){
     ggplot() + theme_void()+ggtitle("Not able to generate a figure")
  })
  
  
}


survival_adj_server2 <-function(input,output,session,inputdata,survtype,input_gene_surv,adj_cov_input){

  tryCatch({

    if(survtype=="OS"){


      dtt = inputdata
      
      dtt[make.names(input_gene_surv)] = log2(dtt[make.names(input_gene_surv)]+1)
      OS_KM_fit_table <- run_pretty_model_output(
        c(adj_cov_input,make.names(input_gene_surv)),
        model_data = dtt,  y_in = "OS.time", event_in = "OS", event_level = '1',
        p_digits = 5) 

        OS_KM_fit_table = data.frame(OS_KM_fit_table[,-c(ncol(OS_KM_fit_table),(ncol(OS_KM_fit_table)-1))])
        colnames(OS_KM_fit_table)[3]="HR (95% CI)"
        OS_KM_fit_table



    }else if(survtype=="PFI"){


      dtt = inputdata
      dtt[make.names(input_gene_surv)] = log2(dtt[make.names(input_gene_surv)]+1)
      

      OS_KM_fit_table <- run_pretty_model_output(
        c(adj_cov_input,make.names(input_gene_surv)),
        model_data = dtt,  y_in = "PFI.time", event_in = "PFI", event_level = '1',
        p_digits = 5) 

        OS_KM_fit_table = data.frame(OS_KM_fit_table[,-c(ncol(OS_KM_fit_table),(ncol(OS_KM_fit_table)-1))])
        colnames(OS_KM_fit_table)[3]="HR (95% CI)"
        OS_KM_fit_table
        


    }

   },error=function(e){
   ""
   })


}
