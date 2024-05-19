#============================================================================
######## TIDE interaction model KM plot ##############
#============================================================================
tide_km_ui = function(gene_tide,cutoff_km_tide,cutoff_km_cyt,cancer_type,cancer_options,plot_output1,plot_output2){
  #-----------------------------------------------
  # This func is for UI part for tide KM tab 
  # gene_tide: genes from tide heatmap
  # cutoff_km_tide: cutoff for genes
  # cutoff_km_cyt: cutoff for CYT
  # cancer_type
  # cancer_optioms
  # plot_output1
  # plot_output2
  #-----------------------------------------------
  fluidRow(
    column(2,
           selectizeInput(gene_tide,"select a gene",NULL,multiple=FALSE),
           sel_input_ui(cutoff_km_tide,"select a cutoff for the gene",c("median","quantile_75","quantile_25")),
           sel_input_ui(cutoff_km_cyt,"select a cutoff for CYT",c("median","quartile","optimal")),
           sel_input_ui(cancer_type,"Cancer type",cancer_options,selected="BRCA")
           
    ),
    column(3,
           div(id="box1",
               box(
                 title="KM for high group",collapsible = TRUE,status="warning",
                 solidHeader = TRUE,
                 width = 12,height=300,
                 # DTOutput("test"),
                 (div(style='max-width: 100%; height: 500px; width: auto;overflow-x: scroll;overflow-y: scroll;text-align: center;',
                  withLoader(plotOutput(plot_output1,width="80%",height="400px"),type="html",loader="loader1")))
               )
           )
    ),
    column(3,
           div(id="box1",
               box(
                 title="KM for low group",collapsible = TRUE,status="warning",
                 solidHeader = TRUE,
                 width = 12,height=300,
                 (div(style='max-width: 100%; height: 500px; width: auto;overflow-x: scroll;overflow-y: scroll;text-align: center;',
                      withLoader(plotOutput(plot_output2,width="80%",height="400px"),type="html",loader="loader1")))
                 
                 
               )
           ))
  )
}




tide_km_data = function(input,output,session){
  #------------------------------------------------------
  # This function is for KM plot dataset
  #--------------------------------------------------------
  if(input$data_type_interact=="lncRNA"){
    dt = read_fst("www/TIDE/lncRNA_clin_CYT.fst",
                  columns = c("OS","OS.time","type",make.names(input$gene_tide),"CYT"))
  }else if(input$data_type_interact=="mRNA"){
    dt = read_fst("www/TIDE/mRNA_clin_CYT.fst",
                  columns = c("OS","OS.time","type",make.names(input$gene_tide),"CYT"))
  }
  
  dt_sub = data.frame(subset(dt,dt$type==input$cancer_type_km_tide))
  names(dt_sub) = c("OS","OS.time","type","gene","CYT")
  if(input$cutoff_km_tide=="median"){
    dt_sub$cut = findInterval(dt_sub$gene, median(dt_sub$gene,na.rm=T))
    dt_sub$cut  = ifelse(dt_sub$cut==1,"high","low")      
  }else if(input$cutoff_km_tide=="quantile_75"){
    dt_sub$cut = findInterval(dt_sub$gene, quantile(dt_sub$gene,na.rm=T))
    dt_sub$cut = ifelse(dt_sub$cut%in%c(1,2,3),"low","high")
  }else if(input$cutoff_km_tide=="quantile_25"){
    dt_sub$cut = findInterval(dt_sub$gene, quantile(dt_sub$gene,na.rm=T))
    dt_sub$cut = ifelse(dt_sub$cut%in%c(1,2),"low","high") 
  }
  
  dt_sub
}

tide_km_fit = function(input,output,session,tide_km_data,hi_lo){
  #------------------------------------------------------------
  # this function is for surfit for hi and lo CTL KM plot
  # hi_lo: CTL high or low
  #------------------------------------------------------------
  dt = tide_km_data
  dt = subset(dt,dt$cut==hi_lo)
  if(input$cutoff_km_cyt=="optimal"){
    res.cut <- surv_cutpoint(dt,time = "OS.time", event = "OS", variables="CYT")
    res.cat <- surv_categorize(res.cut)
    names(res.cat) = c("X1","X2","X3")
  }else if(input$cutoff_km_cyt=="quartile"){
    dt$cut = findInterval(dt$CYT, quantile(dt$CYT,na.rm=T))
    dt$cut = ifelse(dt$cut%in%c(2,3),NA,ifelse(dt$cut%in%c(4,5),"high","low"))
    res.cat = data.frame(cbind(dt$OS.time,dt$OS,dt$cut))
    res.cat$X1 = as.numeric(as.character(res.cat$X1))
    res.cat$X2 = as.numeric(as.character(res.cat$X2))
  }else if(input$cutoff_km_cyt=="median"){
    dt$cut = findInterval(dt$CYT, median(dt$CYT,na.rm=T))
    dt$cut = ifelse(dt$cut==1,"high","low")
    res.cat = data.frame(cbind(dt$OS.time,dt$OS,dt$cut))
    res.cat$X1 = as.numeric(as.character(res.cat$X1))
    res.cat$X2 = as.numeric(as.character(res.cat$X2))
  }
  res.cat
}
