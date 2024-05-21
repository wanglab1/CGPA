#----------------------------------------------------------------------------------
# Adjusted KM
##----------------------------------------------------------------------------------
library(MoffittFunctions)
library(AdjKMCIF)

survival_adj_data = function(cutoff,data,pheno_study){
  #-------------------------------------------------------
  # Adjusted KM plot
  # data: import dataset
  # cutoff: cutoffs for KM
  # data = dt_study;pheno_study =pheno_study;cutoff="optimal"
  #-------------------------------------------------------
  pheno_study = pheno_study %>% select(-c(OS.time, OS, PFS.time, PFS))
  if(cutoff=="optimal"){
    res.cut <- surv_cutpoint(data, time = "time", event = "status", variables = "gene")
    res.cat <- surv_categorize(res.cut)
    names(res.cat)[3] = "cut"
    res.cat = data.frame(cbind(res.cat,"gene"=data[["gene"]],pheno_study))
  }else if(cutoff=="median"){
    dt = data
    dt$cut = findInterval(dt$gene, median(dt$gene,na.rm=T))
    dt$cut = ifelse(dt$cut==1,"high","low")
    res.cat = data.frame(cbind("time"=dt$time,"status"=dt$status,"cut"=dt$cut))
    res.cat$time = as.numeric(as.character(res.cat$time))
    res.cat$status = as.numeric(as.character(res.cat$status))
    res.cat = data.frame(cbind(res.cat,"gene"=data[["gene"]],pheno_study))
    
  }else if(cutoff=="quartile"){
    dt = data
    dt$cut = findInterval(dt$gene, quantile(dt$gene,na.rm=T))
    dt$cut = ifelse(dt$cut%in%c(2,3),NA,ifelse(dt$cut%in%c(4,5),"high","low"))
    res.cat = data.frame(cbind("time"=dt$time,"status"=dt$status,"cut"=dt$cut))
    res.cat$time = as.numeric(as.character(res.cat$time))
    res.cat$status = as.numeric(as.character(res.cat$status))
    res.cat = data.frame(cbind(res.cat,"gene"=data[["gene"]],pheno_study))
    
  }
  return(res.cat)
}


survival_adj_server1 <-function(inputdata,adj_cov){
  #-------------------------------------------------
  # Generate adjusted KM plot
  # input_gene_surv: the searched gene, i.e. CD8A
  # inputdata = adj_dt;adj_cov = c("FMOne_mutation_burden_per_MB","Neoantigen_burden_per_MB")
  #-------------------------------------------------
#  tryCatch({
    
    
    dtt = data.frame(inputdata[,c("time","status",adj_cov,"cut")])
    dtt = na.omit(data.frame(dtt))
    
    res = adjusted_KM(data=dtt,time='time',status="status",group="cut",covlist=adj_cov,stratified_cox = "Yes",reference_group="G&B")
    
    
    p = ggplot(res,aes(x=time,y = prob, group =class))+
      geom_step(aes(color = class),size=0.7)+
      theme_classic()+
      ylim(c(0,1))+
      theme(legend.position="top")+
      scale_color_manual(values=c("#DF8F44FF", "#374E55FF"))+
      ylab("Adjusted survival probability")+
      xlab("Time")
    p
    
  # },error=function(e){
  #   ggplot() + theme_void()+ggtitle("Not able to generate a figure")
  # })
  # 
  
}


survival_adj_server2 <-function(inputdata,adj_cov_input){
  
  #------------------------------------------------------------------------------
  # multivariable cox
  # adj_cov_input: adjusted covariates
  # adj_cov_input=c("FMOne_mutation_burden_per_MB","Neoantigen_burden_per_MB");inputdata = adj_dt
  #------------------------------------------------------------------------------
  tryCatch({
    
    dtt = inputdata
    
    KM_fit_table <- run_pretty_model_output(
      c(adj_cov_input,"gene"),
      model_data = dtt,  y_in = "time", event_in = "status", event_level = '1',
      p_digits = 5) 
    
    
    KM_fit_table = data.frame(KM_fit_table[,-c(ncol(KM_fit_table),(ncol(KM_fit_table)-1))])
    colnames(KM_fit_table)[3]="HR (95% CI)"
      KM_fit_table
    
  },error=function(e){
    ""
  })
  
  
}

