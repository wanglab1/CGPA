ghi_model = function(cov,input_gene,hallmark_dt,response_match,tumor_purity_adj){
  #-------------------------------------------------------------------------
  # This function is to return result for full and partial interaction model 
  # cov: covariates for the interaction model
  # input_gene: input geneName
  # hallmark_dt: dataset for interaction model
  # response_match: OS or PFS
  # tumor_purity_adj: whether to adjust for tumor purity
  #-------------------------------------------------------------------------
  cov = unique(c(make.names(cov)))
  ########### full interaction model ##############
  res = lapply(cov,function(x){
    tryCatch({
      if(response_match=="OS"){
        # if(length(which(geneName$geneName==input_gene))>1){
        #   
        #   gene_sel = paste0(make.names(input_gene),".x")
        #   
        #   if(tumor_purity_adj=="No"){
        #     formula_sex =  as.formula(paste('Surv(OS.time, OS)~',gene_sel,"+",x,"+",paste0(x,":",gene_sel)))
        #   }else{
        #     formula_sex =  as.formula(paste('Surv(OS.time, OS)~',"Tumor_purity+",gene_sel,"+",x,"+",paste0(x,":",gene_sel)))
        #   }
        
        # }else{
        if(tumor_purity_adj=="No"){
          formula_sex =  as.formula(paste('Surv(OS.time, OS)~',make.names(input_gene),"+",x,"+",paste0(x,":",make.names(input_gene))))
        }else{
          formula_sex =  as.formula(paste('Surv(OS.time, OS)~',"Tumor_purity+",make.names(input_gene),"+",x,"+",paste0(x,":",make.names(input_gene))))
          
        }
        #}
        
      }else{
        # if(length(which(geneName$geneName==input_gene))>1){
        #   
        #   gene_sel = paste0(make.names(input_gene),".x")
        #   if(tumor_purity_adj=="No"){
        #     formula_sex =  as.formula(paste('Surv(PFS.time, PFS)~',gene_sel,"+",x,"+",paste0(x,":",gene_sel)))
        #   }else{
        #     formula_sex =  as.formula(paste('Surv(PFS.time, PFS)~',"Tumor_purity+",gene_sel,"+",x,"+",paste0(x,":",gene_sel)))
        #   }
        
        #}else{
        if(tumor_purity_adj=="No"){
          formula_sex =  as.formula(paste('Surv(PFS.time, PFS)~',make.names(input_gene),"+",x,"+",paste0(x,":",make.names(input_gene))))
        }else{
          formula_sex =  as.formula(paste('Surv(PFS.time, PFS)~',"Tumor_purity+",make.names(input_gene),"+",x,"+",paste0(x,":",make.names(input_gene))))
        }
        # }
        
      }
      fit = coxph(formula_sex,data=hallmark_dt)
      surv_res = tidy(fit,exponentiate = T)
      #surv_res$term[2] = paste0(surv_res$term[2],"_","EGFR")    
      res = data.frame(surv_res$statistic)
      
    },error = function(e){
      ""
    })
    
  })
  #res
  if(tumor_purity_adj=="No"){
    res = do.call(cbind,res)
    names(res) = paste0(cov)
    res$group = c("Gene","Hallmark","GMI")
    res = subset(res,res$group=="GMI")
    res$GMI = "Full.model"
    res_full = data.frame(res)
  }else{
    res = do.call(cbind,res)
    names(res) = paste0(cov)
    res$group = c("Tumor_purity","Gene","Hallmark","GMI")
    res = subset(res,res$group=="GMI")
    res$GMI = "Full.model"
    res_full = data.frame(res)
  }
  
  ############## partial interaction model ##########################    
  res = lapply(cov,function(x){
    tryCatch({
      
      if(response_match=="OS"){
        # if(length(which(geneName$geneName==input_gene))>1){
        #   
        #   gene_sel = paste0(make.names(input_gene),".x")
        #   if(tumor_purity_adj=="No"){
        #     formula_sex =  as.formula(paste('Surv(OS.time, OS)~',paste0(x,"/",gene_sel)))
        #   }else{
        #     formula_sex =  as.formula(paste('Surv(OS.time, OS)~',paste0("Tumor_purity+",x,"/",gene_sel)))
        #   }
        # 
        # }else{
        if(tumor_purity_adj=="No"){
          formula_sex =  as.formula(paste('Surv(OS.time, OS)~',paste0(x,"/",make.names(input_gene))))
        }else{
          formula_sex =  as.formula(paste('Surv(OS.time, OS)~',paste0("Tumor_purity+",x,"/",make.names(input_gene))))
        }
        # }
        
      }else{
        # if(length(which(geneName$geneName==input_gene))>1){
        #   
        #   gene_sel = paste0(make.names(input_gene),".x")
        #   if(tumor_purity_adj=="No"){
        #     formula_sex =  as.formula(paste('Surv(PFS.time, PFS)~',paste0(x,"/",gene_sel)))
        #   }else{
        #     formula_sex =  as.formula(paste('Surv(PFS.time, PFS)~',paste0("Tumor_purity+",x,"/",gene_sel)))
        #   }
        # 
        # }else{
        if(tumor_purity_adj=="No"){
          formula_sex =  as.formula(paste('Surv(PFS.time, PFS)~',paste0(x,"/",make.names(input_gene))))
        }else{
          formula_sex =  as.formula(paste('Surv(PFS.time, PFS)~',paste0("Tumor_purity+",x,"/",make.names(input_gene))))
          
        }
        #}
        
      }
      
      fit = coxph(formula_sex,data=hallmark_dt)
      surv_res = tidy(fit,exponentiate = T)
      res = data.frame(surv_res$statistic)
      
      if(tumor_purity_adj=="No"){
        if(nrow(res)==3){
          res = data.frame(res[-2,] )
        }else{
          res
        }
      }else{
        if(nrow(res)==4){
          res = data.frame(res[-3,] )
        }else{
          res
        }
      }
    },error = function(e){
      ""
    })
    
  })
  #res
  if(tumor_purity_adj=="No"){
    res = do.call(cbind,res)
    names(res) = paste0(cov)
    res$group = c("Gene","GMI")
    res = subset(res,res$group=="GMI")
    res$GMI = "Partial.model"
    res_part= res
  }else{
    res = do.call(cbind,res)
    names(res) = paste0(cov)
    res$group = c("Tumor_purity","Gene","GMI")
    res = subset(res,res$group=="GMI")
    res$GMI = "Partial.model"
    res_part= res
  }
  
  
  
  res_plot = data.frame(rbind.fill(res_full,res_part))
  res_plot = res_plot[,c(ncol(res_plot),1:(ncol(res_plot)-2))]
  
  res_plot = data.matrix(res_plot)
  rownames(res_plot) = c("Full.model","Partial.model");res_plot = data.frame(res_plot[,-1])
  res_plot = as.matrix(res_plot)
  res_plot
}