#==================================================================================
############## LncRNA tab: top prognostic genes across cancers##################
#==================================================================================
common_ui = function(gene_type,cancer_type,OS_PFI,adj_cov,rank_order,top_n,font_size,action_button,cancer_sel,output_plot){
  #--------------------------------------------------
  # UI for single gene GHI model 
  # UI for top prognostic genes across cancers
  # gene_type:lncRNA or mRNA
  # cancer_type: cancer_type
  # OS_PFI: OS or PFI
  # adj_cov: adjusted covariates
  # rank_order
  # top_n: top n numbers
  # font_size: font size for the heatmap
  # action_button
  #--------------------------------------------------
  tagList(
    fluidRow(
    column(3,
           sel_input_ui(gene_type,"lncRNA or mRNA",choices = c("mRNA","lncRNA")),
           sel_input_ui(cancer_type,"Cancer types",choices = cancer_sel, selected=c("HNSC","CESC","ESCA","LUSC","BLCA"),multiple=T),
           sel_input_ui(OS_PFI,"OS or PFI",c("OS","PFI")),
           sel_input_ui(adj_cov,"Adjusted covariates",choices = c("None","Age + Sex","Age + Sex + tumor_Purity")),
           sel_input_ui(rank_order,"Rank by average z pos or neg",c("Positive","Negative")),
           numericInput(top_n,"Top n numbers",30),
           numericInput(font_size,"Adjust font size",value = 10),
           p("Click 'Run' button to get the figure"),
           actionButton(action_button, "Run")
           ),
    column(9,
           withLoader(plotOutput(output_plot,width="60%",height="800px"),type="html",loader="loader1")
           
           # conditionalPanel(
           #   condition= paste0("input.",action_button>0),
           #   withLoader(plotOutput(output_plot,width="60%",height="800px"),type="html",loader="loader1")
           # ),
           # conditionalPanel(
           #   condition= paste0("input.",action_button==0),
           #   br(),
           #   p("Click 'Run' button to get the figure")
           # )
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
             sel_input_ui(gene_single_ghi,"Choose genes",choices=NULL,multiple=T),
             
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

common_dt_server = function(input,output,session,geneName){
  #---------------------------------------------------------------
  # This function is to generate dataset for the common heatmap
  # geneName
  #--------------------------------------------------------------
  if(input$lnc_mrna_a=="lncRNA"){
    folder = "www/top_lnc_combined/"
  }else if(input$lnc_mrna_a=="mRNA"){
    folder = "www/top_mrna_combined/"
  }
  if(input$OS_PFI_a=="OS"){
    if(input$surv_cov_a=="None"){
      dt = read_fst(paste0(folder,"combined_OS_FALSE_FALSE.fst"))
    }else if(input$surv_cov_a=="Age + Sex"){
      dt = read_fst(paste0(folder,"combined_OS_TRUE_FALSE.fst"))
    }else{
      dt = read_fst(paste0(folder,"combined_OS_TRUE_TRUE.fst"))
    }
  }else if(input$OS_PFI_a=="PFI"){
    if(input$surv_cov_a=="None"){
      dt = read_fst(paste0(folder,"combined_PFI_FALSE_FALSE.fst"))
    }else if(input$surv_cov_a=="Age + Sex"){
      dt = read_fst(paste0(folder,"combined_PFI_TRUE_FALSE.fst"))
    }else{
      dt = read_fst(paste0(folder,"combined_PFI_TRUE_TRUE.fst"))
    }
  }
  combined.df_sub = subset(dt,dt$type%in%input$cancer_surv_multiple)
  combined.df_sub = data.frame(combined.df_sub[,c("Genes","z_score","type")])
  data_wide_raw = spread(combined.df_sub, type, z_score)
  # 
  combined.df_sub_1 = combined.df_sub
  #combined.df_sub_1$z_score = abs(combined.df_sub_1$z_score )
  data_wide = spread(combined.df_sub_1, type, z_score )
  data_wide$z_mean = rowMeans(data_wide[,-1],na.rm = T)
  data_wide[,-c(1,ncol(data_wide))] = data_wide_raw[,-c(1,ncol(data_wide))]
  # 
  if(input$rank_z_a=="Positive"){
    data_wide = subset(data_wide,data_wide$z_mean>=0)
    data_wide = data_wide[order(data_wide$z_mean,decreasing = T),]
    
    data_wide_plot = data_wide[1:input$Top_n_a,]
  }else if(input$rank_z_a=="Negative"){
    data_wide = subset(data_wide,data_wide$z_mean<0)
    data_wide = data_wide[order(data_wide$z_mean,decreasing = F),]
    
    data_wide_plot = data_wide[1:input$Top_n_a,]
    
  }
  # data_wide_plot$z_mean = sapply(1:nrow(data_wide_plot),function(x){
  #   # all(data_wide_plot[x,c(-1,-ncol(data_wide_plot))]<=0)
  #   ifelse( all(data_wide_plot[x,c(-1,-ncol(data_wide_plot))]<=0),-(data_wide_plot$z_mean[x]),data_wide_plot$z_mean[x])
  #   
  # })
  data_wide_plot$Genes = factor(data_wide_plot$Genes,levels = rev(data_wide_plot$Genes))
  rownames(data_wide_plot) = data_wide_plot$Genes
  names(data_wide_plot)[ncol(data_wide_plot)] = "z mean"
  data_wide_plot = data.frame(data_wide_plot[,-1])
  
  geneName = geneName
  geneName$make.names = make.names(geneName$geneName)
  data_wide_plot = merge(geneName,data_wide_plot,by.x = "make.names",by.y = "row.names")
  data_wide_plot = data_wide_plot[!duplicated(data_wide_plot$geneName),]
  
  rownames(data_wide_plot) = data_wide_plot$geneName
  data_wide_plot = data_wide_plot[,-c(1:2)]
  data_wide_plot
}
