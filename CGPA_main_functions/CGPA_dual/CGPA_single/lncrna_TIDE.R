#============================================================================
################# Functions for TIDE #####################
#============================================================================
tide_ui = function(gene_type,cancer_type,z_score,top_n,action_button,cancer_options,plot_output,table_ouput){
#----------------------------------------------------------
# This functions is for UI part for TIDE-Interaction 
# gene_type: lncRNA or mRNA
# cancer_type: cancer type
# z_score: postive or negative z score
# top_n: top n features
# action_button: action button
# plot_output: output figure
# table_ouput:output table
#-----------------------------------------------------------
  fluidRow(
    column(2,
           sel_input_ui(gene_type,"lncRNA or mRNA",choices = c("lncRNA","mRNA")),
           sel_input_ui(cancer_type,"Cancer types",choices = cancer_options, selected=c("SKCM","HNSC","UCEC","BRCA"),multiple=T),
           radioButtons(z_score,"Z score positive or negative",c("positive","negative")),
           
           numericInput(top_n,"show top features",value=150,min = 0),
           actionButton(action_button, "Run"),

    ),
    column(5,
           
           div(id="box1",
               box(
                 title="Heatmap of common genes ",collapsible = TRUE,status="warning",
                 solidHeader = TRUE,
                 width = 12,height=700,
                 (div(style='max-width: 100%; height: 500px; width: auto;overflow-x: scroll;overflow-y: scroll;text-align: center;',
                      withLoader(plotOutput(plot_output,height = "500px"),type="html",loader="loader1")))      
                 
                 
                 
                 
               )
           )
    ),
    column(5,
           div(id="box1",
               box(
                 title="p values of common genes ",collapsible = TRUE,status="warning",
                 solidHeader = TRUE,
                 width = 12,height=600,
                 withLoader(DTOutput(table_ouput),type="html",loader="loader1")
               )
           )
    )
  )
}

tide_dt = function(input,output,session){
  #-------------------------------------------------------------------------
  # This function is to generate data for tide 
  #---------------------------------------------------------------------------
  if(input$data_type_interact=="lncRNA"){
    tide = read_fst("www/TIDE/combined_TIDE_lncRNA.fst")
  }else if(input$data_type_interact=="mRNA"){
    tide = read_fst("www/TIDE/combined_TIDE_mRNA.fst")
  }
  combined.df_sub = subset(tide,tide$type%in%input$cancer_type_interact)
  if(input$z_score_interact=="positive"){
    combined.df_sub_sub = subset(combined.df_sub,combined.df_sub$z>=0)
    
  }else{
    combined.df_sub_sub = subset(combined.df_sub,combined.df_sub$z<0)
    
  }
  # Top N highest values by group
  data_new2 <- combined.df_sub_sub %>%      
    group_by(type) %>%
    slice(1:input$top_n_interact)
  
  L = split(data_new2$geneName,data_new2$type)
  names(L) = NULL
  common = Reduce(intersect,L) # Find common genes among cancer types
  
  # If exist in two cancer type, then keep
  comb_num = combn(length(L),2)
  common_two = lapply(1:length(L),function(x){
    intersect(L[[comb_num[1,x]]],L[[comb_num[2,x]]])
  })
  common_two = unique(unlist(common_two))
  common_final = unique(c(common,common_two))
  
  data_new3 = combined.df_sub[which(combined.df_sub$geneName%in%common_final),]
  data_new3 = data_new3[order(data_new3$z),]
  data_new3
}

dt_tide_heat = function(input,output,session,dt_tide,geneName){
  #-------------------------------------------------
  # this func is for dataset for TIDE heatmap
  # dt_tide: output from tide_dt() func
  # geneName
  #--------------------------------------------------

    data_new= dt_tide[,c("geneName","z","type")]
    data_new = spread(data_new,type,z)
    rownames(data_new) = data_new$geneName
    data_new = data.frame(data_new[,-1])
    
    geneName = geneName
    geneName$make.names = make.names(geneName$geneName)
    data_new = merge(geneName,data_new,by.x = "make.names",by.y = "row.names")
    data_new = data_new[!duplicated(data_new$geneName),]
    
    rownames(data_new) = data_new$geneName
    data_new[,-c(1:2)]
  
}

dt_tide_pval= function(input,output,session,dt_tide,geneName){
  #-------------------------------------------------
  # this func is for dataset for TIDE p value table
  # dt_tide: output from tide_dt() func
  # geneName
  #--------------------------------------------------
  data_new= dt_tide[,c("geneName","pval","type")]
  data_new$pval = round(data_new$pval,4)
  data_new = spread(data_new,type,pval)
  rownames(data_new) = data_new$geneName
  data_new = data.frame(data_new[,-1])
  data_new
  
  geneName = geneName
  geneName$make.names = make.names(geneName$geneName)
  data_new = merge(geneName,data_new,by.x = "make.names",by.y = "row.names")
  data_new = data_new[!duplicated(data_new$geneName),]
  
  rownames(data_new) = data_new$geneName
  data_new[,-c(1:2)]
}
