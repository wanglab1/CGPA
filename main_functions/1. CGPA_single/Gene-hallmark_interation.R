###############################################################################
############ GHI model ###############
#############################################################################

ghi_ui_single = function(box_id,title,cancer_type_ghi,gene_single_ghi,survtype_ghi,hallmark_match,tumor_purity_ghi,
                         output_dropdown,output_plot,action_button){ 
  #-------------------------------------------------------------
  # UI for single gene GHI model 
  # box_id: box_id for each panel, different box has diff styles
  # title: titles for each box
  # output_dropdown: dropdown output table
  # output_plot: output figure name
  # output_tb: univariate cox model
  #-------------------------------------------------------------
  tagList(
    fluidRow(
      column(3,
             selectizeInput(cancer_type_ghi,"Cancer type",
                            choices=    c("BLCA","BRCA","CESC","HNSC","KIRC",
                                          "KIRP","LGG","LIHC","LUAD","LUSC","OV",
                                          "PRAD","SKCM","STAD","THCA","UCEC"),multiple=F,selected="HNSC"),
             sel_input_ui(gene_single_ghi,"Choose a gene",choices=NULL,multiple=F),
             
             selectizeInput(survtype_ghi,"OS or PFI",c("OS","PFI"),"OS"),
             sel_input_ui(hallmark_match,"Choose hallmarks",choices = NULL,multiple=T),
             radioButtons(tumor_purity_ghi,"Adjust for tumor purity?",c("Yes","No"),"Yes"),
             actionButton(action_button,"Run")
             ),
      column(9,
             div(id=box_id,
                 box(
                   title=title,collapsible = TRUE, status="warning",
                   solidHeader = TRUE,
                   width = 8,height=800,
                   dropdown(
                     p("GHI full model: Hazard ~ a x Gene + b x Hallmark + c x Gene X Hallmark"),
                     p("GHI partial model: Hazard ~ a x Gene + c x Gene x Hallmark"),
                     br(),
                     p("z score of the interaction term between the gene and hallmark"),
                     DTOutput(output_dropdown),
                     circle = TRUE, status = "warning", icon = icon("gear"), width = "600px",
                     tooltip = tooltipOptions(title = "Click to see more")
                   ),
                   tags$style(HTML('#sw-content-dropdown, .sw-dropdown-in {background-color: gray;}')),
                   #DTOutput("test_match")
                   
                   conditionalPanel(
                     condition= paste0("input.",action_button,">0"),
                     p("Interaction model"),
                     div(withLoader(plotOutput(output_plot,width="80%",height="400px"),type="html",loader="loader1"),align="center")
                   ),
                   conditionalPanel(
                     condition= paste0("input.",action_button,"==0"),
                     br(),
                     p("Click 'Run' button to get/update the figure")
                   )
                 )
                 
             )
             )
    )
  )

  
}

ghi_ui_multi = function(box_id,title,cancer_type_ghi,gene_single_ghi,survtype_ghi,hallmark_match,tumor_purity_ghi,
                         output_dropdown,output_plot,action_button){ 
  #-------------------------------------------------------------
  # UI for multi gene GHI model 
  # box_id: box_id for each panel, different box has diff styles
  # title: titles for each box
  # output_dropdown: dropdown output table
  # output_plot: output figure name
  # output_tb: univariate cox model
  #-------------------------------------------------------------
  tagList(
    fluidRow(
      column(3,
             selectizeInput(cancer_type_ghi,"Cancer type",
                            choices=    c("BLCA","BRCA","CESC","HNSC","KIRC",
                                          "KIRP","LGG","LIHC","LUAD","LUSC","OV",
                                          "PRAD","SKCM","STAD","THCA","UCEC"),multiple=F,selected="HNSC"),
             sel_input_ui(gene_single_ghi,"Choose genes (maximum allowance is 5)",choices=NULL,multiple=T),
             
             selectizeInput(survtype_ghi,"OS or PFI",c("OS","PFI"),"OS"),
             sel_input_ui(hallmark_match,"Choose hallmarks",choices = NULL,multiple=T),
             radioButtons(tumor_purity_ghi,"Adjust for tumor purity?",c("Yes","No"),"Yes"),
             actionButton(action_button,"Run")
      ),
      column(9,
             div(id=box_id,
                 box(
                   title=title,collapsible = TRUE, status="warning",
                   solidHeader = TRUE,
                   width = 8,height=800,
                   dropdown(
 
                     br(),
                     p("z score of the interaction term between the gene and hallmark (full model)"),
                     DTOutput(output_dropdown),
                     circle = TRUE, status = "warning", icon = icon("gear"), width = "600px",
                     tooltip = tooltipOptions(title = "Click to see more")
                   ),
                   tags$style(HTML('#sw-content-dropdown, .sw-dropdown-in {background-color: gray;}')),
                   #DTOutput("test_match")
                   
                   conditionalPanel(
                     condition= paste0("input.",action_button,">0"),
                     p("Interaction model"),
                     div(withLoader(plotOutput(output_plot,width="80%",height="400px"),type="html",loader="loader1"),align="center"),
                   ),
                   conditionalPanel(
                     condition= paste0("input.",action_button,"==0"),
                     br(),
                     p("Click 'Run' button to get/update the figure")
                   )
                 )
                 
             )
      )
    )
  )
  
  
}

###############################################################################
############ GHI model functions ###############
#############################################################################
ghi_model = function(cov,input_gene,hallmark_dt,response_match,geneName,tumor_purity_adj){
  #-------------------------------------------------------------------------
  # This function is to return result for full and partial interaction model 
  # cov: covariates for the interaction model
  # input_gene: input geneName
  # hallmark_dt: dataset for interaction model
  # response_match: OS or PFI
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
        #     formula_sex =  as.formula(paste('Surv(PFI.time, PFI)~',gene_sel,"+",x,"+",paste0(x,":",gene_sel)))
        #   }else{
        #     formula_sex =  as.formula(paste('Surv(PFI.time, PFI)~',"Tumor_purity+",gene_sel,"+",x,"+",paste0(x,":",gene_sel)))
        #   }

        #}else{
          if(tumor_purity_adj=="No"){
            formula_sex =  as.formula(paste('Surv(PFI.time, PFI)~',make.names(input_gene),"+",x,"+",paste0(x,":",make.names(input_gene))))
          }else{
            formula_sex =  as.formula(paste('Surv(PFI.time, PFI)~',"Tumor_purity+",make.names(input_gene),"+",x,"+",paste0(x,":",make.names(input_gene))))
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
    res$GMI = "GMI.full"
    res_full = data.frame(res)
  }else{
    res = do.call(cbind,res)
    names(res) = paste0(cov)
    res$group = c("Tumor_purity","Gene","Hallmark","GMI")
    res = subset(res,res$group=="GMI")
    res$GMI = "GMI.full"
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
        #     formula_sex =  as.formula(paste('Surv(PFI.time, PFI)~',paste0(x,"/",gene_sel)))
        #   }else{
        #     formula_sex =  as.formula(paste('Surv(PFI.time, PFI)~',paste0("Tumor_purity+",x,"/",gene_sel)))
        #   }
        # 
        # }else{
          if(tumor_purity_adj=="No"){
            formula_sex =  as.formula(paste('Surv(PFI.time, PFI)~',paste0(x,"/",make.names(input_gene))))
          }else{
            formula_sex =  as.formula(paste('Surv(PFI.time, PFI)~',paste0("Tumor_purity+",x,"/",make.names(input_gene))))
            
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
    res$GMI = "GMI.part"
    res_part= res
  }else{
    res = do.call(cbind,res)
    names(res) = paste0(cov)
    res$group = c("Tumor_purity","Gene","GMI")
    res = subset(res,res$group=="GMI")
    res$GMI = "GMI.part"
    res_part= res
  }
  

  
  res_plot = data.frame(rbind.fill(res_full,res_part))
  res_plot = res_plot[,c(ncol(res_plot),1:(ncol(res_plot)-2))]
  
  res_plot = data.matrix(res_plot)
  rownames(res_plot) = c("GHI.full","GHI.part");res_plot = data.frame(res_plot[,-1])
  res_plot = as.matrix(res_plot)
  res_plot
}

ghi_model_multi = function(cov,input_gene,hallmark_dt,response_match,geneName,tumor_purity_adj){
  #-------------------------------------------------------------------------
  # This function is to return result for multiple genes interaction model
  # cov: covariates for the interaction model
  # input_gene: input geneName
  # hallmark_dt: dataset for interaction model
  # response_match: OS or PFI
  # tumor_purity_adj: whether to adjust for tumor purity
  #-------------------------------------------------------------------------
  cov = unique(c(make.names(cov)))
  genes = make.names(input_gene)
  ########### full interaction model ##############
  res_multi = lapply(1:length(genes),function(y){
   res = lapply(cov,function(x){
    tryCatch({
      if(response_match=="OS"){
        # if(length(which(geneName$geneName==genes[y]))>1){
        #   
        #   if(tumor_purity_adj=="No"){
        #     formula_sex =  as.formula(paste('Surv(OS.time, OS)~',paste0(make.names(genes[y]),".x"),"+",x,"+",paste0(x,":",paste0(make.names(genes[y]),".x"))))
        #   }else{
        #     formula_sex =  as.formula(paste('Surv(OS.time, OS)~',"Tumor_purity+",paste0(make.names(genes[y]),".x"),"+",x,"+",paste0(x,":",paste0(make.names(genes[y]),".x"))))
        #   }
        #   
        #}else{
          if(tumor_purity_adj=="No"){
            formula_sex =  as.formula(paste('Surv(OS.time, OS)~',make.names(genes[y]),"+",x,"+",paste0(x,":",make.names(genes[y]))))
          }else{
            formula_sex =  as.formula(paste('Surv(OS.time, OS)~',"Tumor_purity+",make.names(genes[y]),"+",x,"+",paste0(x,":",make.names(genes[y]))))
            
          }
       # }
        
      }else{
        # if(length(which(geneName$geneName==input_gene))>1){
        #   
        #   if(tumor_purity_adj=="No"){
        #     formula_sex =  as.formula(paste('Surv(PFI.time, PFI)~',paste0(make.names(genes[y]),".x"),"+",x,"+",paste0(x,":",paste0(make.names(genes[y]),".x"))))
        #   }else{
        #     formula_sex =  as.formula(paste('Surv(PFI.time, PFI)~',"Tumor_purity+",paste0(make.names(genes[y]),".x"),"+",x,"+",paste0(x,":",paste0(make.names(genes[y]),".x"))))
        #   }
          
       # }else{
          if(tumor_purity_adj=="No"){
            formula_sex =  as.formula(paste('Surv(PFI.time, PFI)~',make.names(genes[y]),"+",x,"+",paste0(x,":",make.names(genes[y]))))
          }else{
            formula_sex =  as.formula(paste('Surv(PFI.time, PFI)~',"Tumor_purity+",make.names(genes[y]),"+",x,"+",paste0(x,":",make.names(genes[y]))))
          }
        #}
        
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
    res$gene = genes[y]
    res_full = data.frame(res)
  }else{
    res = do.call(cbind,res)
    names(res) = paste0(cov)
    res$group = c("Tumor_purity","Gene","Hallmark","GMI")
    res = subset(res,res$group=="GMI")
    res$gene = genes[y]
    res_full = data.frame(res)
  }
   
})
  res_multi= data.frame(do.call(rbind,res_multi))
  rownames(res_multi) = res_multi$gene
  res_multi = as.matrix(res_multi[,1:(ncol(res_multi)-2)])
  return(res_multi)
}

ghi_model_common = function(input,output,session,cancer_type,hallmarks,top_n,cancer_dt){
#-------------------------------------------------------------------------------
# Find common genes across hallmarks or cancer types for GHI
# cancer_type: inputcancer types
# hallmarks: age, sex, etc.
# top_n: top prognostic markers
# cancer_dt: input cancer datasets (list format, can be multiple)
#-------------------------------------------------------------------------------
  n_cancer_type = length(cancer_type); n_hallmark = length(hallmarks)
  #################### common genes across hallmarks in one cancer type ########################
  if(n_cancer_type==1&n_hallmark==1){
    cancer_dt = cancer_dt[[1]]
    sel_vars = c(paste0("geneName_",hallmarks),paste0(hallmarks))
    dt = cancer_dt[,sel_vars]
    dt = dt[order(abs(dt[[hallmarks]]),decreasing=T),]
    dt = dt[1:top_n,]
    names(dt)[1] ="geneName"
    dt = subset(dt,abs(dt[[hallmarks]])>1.96)
    dt_final = dt
  }else if(n_cancer_type==1&n_hallmark>=2){
    cancer_dt = cancer_dt[[1]]
    geneNames = lapply(1:length(hallmarks),function(x){
      sel_vars = c(paste0("geneName_",hallmarks[x]),paste0(hallmarks[x]))
      dt = cancer_dt[,sel_vars]
      dt = dt[order(abs(dt[[hallmarks[x]]]),decreasing=T),]
      dt = dt[1:top_n,]
      names(dt)[1] ="geneName"
      dt = subset(dt,abs(dt[[hallmarks[x]]])>1.96)
      dt$geneName
    })
    
    common = Reduce(intersect, geneNames)
    
    # If exists in at least 2 hallmarks, then keep
    comb_num = combn(length(geneNames),2)
    common_two = lapply(1:length(geneNames),function(x){
      intersect(geneNames[[comb_num[1,x]]],geneNames[[comb_num[2,x]]])
    })
    common_two = unique(unlist(common_two))
    common_final = unique(c(common,common_two))
    
    dt = lapply(1:length(hallmarks),function(x){
      sel_vars = c(paste0("geneName_",hallmarks[x]),paste0(hallmarks[x]))
      dt = cancer_dt[,sel_vars]
      names(dt)[1]="geneName"
      dt = subset(dt,dt$geneName%in%common_final)
      dt
    })
    
    dt_final = Reduce(function(x, y) merge(x, y, all=TRUE), dt)
    
  } else if(n_cancer_type>=2&n_hallmark ==1){
    #################### common genes across cancertypes in one hallmark########################
    geneNames = lapply(1:length(cancer_dt),function(x){
      tryCatch({
        sel_vars = c(paste0("geneName_",hallmarks),paste0(hallmarks))
        dt = cancer_dt[[x]][,sel_vars]
        dt = dt[order(abs(dt[[hallmarks]]),decreasing=T),]
        dt = dt[1:top_n,]
        names(dt)[1] ="geneName"
        dt = subset(dt,abs(dt[[hallmarks]])>1.96)
        dt$geneName    
      },error = function(e){
        ""
      })
      
    })
    
    
    common = Reduce(intersect, geneNames)
    
    # If exists in at least 2 hallmarks, then keep
    comb_num = combn(length(geneNames),2)
    common_two = lapply(1:length(geneNames),function(x){
      intersect(geneNames[[comb_num[1,x]]],geneNames[[comb_num[2,x]]])
    })
    common_two = unique(unlist(common_two))
    common_final = unique(c(common,common_two))
    
    dt = lapply(1:length(cancer_dt),function(x){
      tryCatch({
        sel_vars = c(paste0("geneName_",hallmarks),paste0(hallmarks))
        dt = cancer_dt[[x]][,sel_vars]
        dt = dt[order(abs(dt[[hallmarks]]),decreasing=T),]
        dt = dt[1:top_n,]
        names(dt)[1] ="geneName"
        names(dt)[2] = paste0(cancer_type[x],"_",names(dt)[2])
        dt = subset(dt,dt$geneName%in%common_final)
        dt
      },error = function(e){
        ""
      })
      
    })
    
    dt_final = Reduce(function(x, y) merge(x, y, all=TRUE), dt)
  }else{
    #################### common genes across cancertypes in multiple hallmarks ########################
    geneNames = lapply(1:length(cancer_dt),function(x){
      lapply(1:length(hallmarks),function(y){
        tryCatch({
          sel_vars = c(paste0("geneName_",hallmarks[[y]]),paste0(hallmarks[[y]]))
          dt = cancer_dt[[x]][,sel_vars]
          dt = dt[order(abs(dt[[hallmarks[[y]]]]),decreasing=T),]
          dt = dt[1:top_n,]
          names(dt)[1] ="geneName"
          dt = subset(dt,abs(dt[[hallmarks[[y]]]])>1.96)
          dt$geneName    
        },error = function(e){
          ""
        })
        
      })  
    })
    
    # expand the list of list into list
    geneNames = lapply(rapply(geneNames, enquote, how="unlist"), eval)
    
    common = Reduce(intersect, geneNames)
    
    # If exists in at least 2 hallmarks, then keep
    comb_num = combn(length(geneNames),2)
    common_two = lapply(1:length(geneNames),function(x){
      intersect(geneNames[[comb_num[1,x]]],geneNames[[comb_num[2,x]]])
    })
    common_two = unique(unlist(common_two))
    common_final = unique(c(common,common_two))
    
    dt = map(1:length(cancer_dt),function(x){
      map(1:length(hallmarks),function(y){
        tryCatch({
          sel_vars = c(paste0("geneName_",hallmarks[[y]]),paste0(hallmarks[[y]]))
          dt = cancer_dt[[x]][,sel_vars]
          dt = dt[order(abs(dt[[hallmarks[[y]]]]),decreasing=T),]
          dt = dt[1:top_n,]
          names(dt)[1] ="geneName"
          names(dt)[2] = paste0(cancer_type[x],"_",names(dt)[2])
          dt = subset(dt,dt$geneName%in%common_final)    
        },error = function(e){
          ""
        })
        
      })  
    })
    
    # expand the list of list into list
    dt = rrapply::rrapply(dt, is.data.frame, classes = 'data.frame', how = 'flatten')
    
    
    dt_final = Reduce(function(x, y) merge(x, y, all=TRUE), dt)
    
  }

  
  return(dt_final)
}