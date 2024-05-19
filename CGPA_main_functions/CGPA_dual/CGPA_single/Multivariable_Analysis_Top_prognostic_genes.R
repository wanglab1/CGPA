#==============================================================================
############# Top prognostic results #########################
# 06/30/2022
#==============================================================================
top_prog_ui = function(box_id,title,condition,output,drop_down_info,subnet_action=NULL,output_subnet=NULL){
  #----------------------------------------------------
  # This function is for top prognostic tab UI part
  # title: title of the panel
  # condition: input.search for conditional panel
  # output: output info
  #----------------------------------------------------
  div(id=box_id,
      box(
        title=title,collapsible = TRUE, status="warning",
        solidHeader = TRUE,
        width = 12,height="130%",

        conditionalPanel(
          condition = paste0("input.",condition,">0"),
          withLoader(DTOutput(output),type="html",loader="loader1"),
          br(),
          div(
            if(box_id=="box1"){
              tagList(
                p("Click to see the subnetwork of the selected gene among the top protein-coding genes"),
                actionButton(subnet_action, "Subnetwork"),
                bsModal("subnet_modal", "", subnet_action, size = "large",
                        withLoader(plotOutput(output_subnet,width="100%",height="800px"),type="html",loader="loader1")
                )
              )
              
            }else{
              div(
              p("Subnetwork only works for protein-coding genes"),
              
              br(),br(),br()
              )
            }
          )
          
        ),
        conditionalPanel(
          condition=paste0("input.",condition,"==0"),
          br(),
          p("Click 'Run' button to get/update the table")
        )
      ))
}

top_prog_server = function(file,adj_cov,cancer_type,surv_type){
  #------------------------------------------------------------
  # Function for top prognostic genes
  # adj_cov: adjusted covariates
  # This function is for top prognostic table
  # file: locations
  # cancer_type
  # surv_type: OS or PFI
  #------------------------------------------------------------
  file = file

  if(adj_cov=="None"){
    res = read_fst(paste0(file,cancer_type,"_",surv_type,"_","FALSE_FALSE.fst"))
    #res = read_fst(paste0(file,fst_file[grepl(cancer_type,fst_file)&grepl("unadj",fst_file)]))
    #res = data.frame(res[,-2])
    #names(res) = c("Genes","HR (95% CI)","P.Value")
  }else if(adj_cov=="Age + Sex"){
    res = read_fst(paste0(file,cancer_type,"_",surv_type,"_","TRUE_FALSE.fst"))
    
    # res = read_fst(paste0(file,fst_file[grepl(cancer_type,fst_file)&grepl("adj",fst_file)]))
  }else if(adj_cov=="Age + Sex + tumor_Purity"){
    res = read_fst(paste0(file,cancer_type,"_",surv_type,"_","TRUE_TRUE.fst"))
  } else if(adj_cov=="Age + Sex + CTL"){
    res = read_fst(paste0(file,cancer_type,"_",surv_type,"_","CTL.fst"))
  }
  
  res = res[order(res$p.value,decreasing = F),]
  res = res[which(res$HR<=20),]
  if(cancer_type%in%c("BRCA-BASAL","BRCA-NON-BASAL","SKCM_Metastasis")){
    res = res
  } else{
    if(adj_cov=="Age + Sex + CTL"){
      res = res
    }else{
      res = res[,-ncol(res)]
    }
    
  }
  res
}

top_subnet = function(top_tb,input_gene,cancer_type){
  #------------------------------------------------------
  # Function for the dataset of subnetwork
  # top_tb: top prognostic protein-coding genes table
  # input_gene: input gene
  #-------------------------------------------------------
    dt_top = top_tb
    dt_top = data.frame(dt_top[order(dt_top$p.value,decreasing=FALSE),])
    dt_top = dt_top[1:100,]
    sel_genes = dt_top$Genes
    
    
  
      if(input_gene!=""){
        dt = read_fst("www/pca/mRNA_clin.fst",columns = c("Row.names","type",c(make.names(input_gene),make.names(sel_genes))))
        dt_sub = subset(dt,dt$type==cancer_type)
        dt_sub_ega = dt_sub[,c(make.names(input_gene),make.names(sel_genes))]
        names(dt_sub_ega)[1]=paste0("--->",names(dt_sub_ega)[1],"<---")
      }else{
        dt = read_fst("www/pca/mRNA_clin.fst",columns = c("Row.names","type",make.names(sel_genes)))
        dt_sub = subset(dt,dt$type==cancer_type)
        dt_sub_ega = dt_sub[,c(make.names(sel_genes))]
      }
      
    
    dt_sub_ega
  
}
