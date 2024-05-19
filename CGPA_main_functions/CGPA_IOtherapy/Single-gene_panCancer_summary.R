#===============================================================================
# CGPA ICI
#===============================================================================
#=================================================================================
####################### Pan cancer dashboard tab ################################
#=================================================================================
pan_ui = function(box_id,title,height,output_fig){ 
  #-------------------------------------------------------------
  # box_id: box_id for each panel, different box has diff styles
  # title: titles for each box
  # height: box height
  # output figure name
  #-------------------------------------------------------------
  div(id=box_id,
      box(
        title=title,collapsible = TRUE, 
        status="warning",
        solidHeader = TRUE,
        width = 12,height=height,
        
        if(box_id=="box4"){
          div(style='max-width: 100%; height: 500px;overflow-x: auto;',
              div(style="text-align: left;",uiOutput(output_fig))
              
          )
        }else{
          (div(style='width: auto;max-width: 100%; height: 500px; ',
               div(style="text-align: center;width: auto;",imageOutput(output_fig)),
          ))
          
        }
        
      )
  )
}

pan_ui2 = function(box_id,title,height,output_fig){ 
  #-------------------------------------------------------------
  # for summary page KM plot
  # box_id: box_id for each panel, different box has diff styles
  # title: titles for each box
  # height: box height
  # output figure name
  #-------------------------------------------------------------
  div(id=box_id,
      box(
        title=title,collapsible = TRUE, 
        status="warning",
        solidHeader = TRUE,
        dropdown(
          selectizeInput("KM_cutoff_pan","KM cutoff",choices =c("median","optimal","quartile"),selected="optimal"),
          circle = TRUE, status = "warning", icon = icon("gear"), width = "300px",
          tooltip = tooltipOptions(title = "Click to see inputs !")
        ),
        width = 12,height=height,
        div(style='max-width: 100%; height: 500px;overflow-x: auto;',div(style="text-align: left;",uiOutput(output_fig)) )
  
        
      )
  )
}

pan_server_fig = function(box_id,file=NULL,file_name=NULL,input_gene){
  #----------------------------------------------------------------------
  # functions to import images
  # box_id: to identify functions in different panel
  # file: file path
  # file_name: file names
  # input_gene: geneNames
  #----------------------------------------------------------------------
  if(box_id=="box4"){
    src=paste0("https://string-db.org/api/svg/network?identifiers=%0D",input_gene,"&limit=10&network_flavor=evidence&species=9606")
    src
  }else{
    tryCatch({
      file = file
      fig = paste0(file,make.names(input_gene),file_name)
      list(src=fig,alt="fig")
    },error=function(e){
      ""
    })
  }
  
}

#================================================================================
# forest plot 
#================================================================================
safe_coxph <- function(data, formula) {
  tryCatch(
    {
      tidy(coxph(formula, data = data), exponentiate = TRUE, conf.int = TRUE)
    },
    error = function(e) {
      data.frame(
        term = NA,
        estimate = NA,
        std.error = NA,
        statistic = NA,
        p.value = NA,
        conf.low = NA,
        conf.high = NA
      )
    }
  )
}

forest_table_server = function(input,output,session,input_dt,surv_type,sel_gene){
  #----------------------------------------------------------------------
  # functions to generate univariable cox
  # input_dt: dataset
  # surv_type: OS or PFI
  # sel_gene: searched gene by user
  #----------------------------------------------------------------------
  
 dt_sub = input_dt
  
 if(surv_type=="OS"){
   os_formula = as.formula(paste("Surv(OS.time,OS)~", paste0(paste(sel_gene, collapse=" + "))))
   
   os_res <- dt_sub %>%
     group_by(study) %>%
     do(safe_coxph(., os_formula))
   
   os_res$symbol = ifelse(os_res$estimate>=1,">=1","<1")
   os_res = data.frame(os_res[order(os_res$estimate),])
   names(os_res) = c("study","term","HR","std","statistics","pval","HR_lower","HR_upper","symbol")
   res = os_res
 }else{
   # PFS
   pfs_formula = as.formula(paste("Surv(PFS.time,PFS)~", paste0(paste(sel_gene, collapse=" + "))))
   
   
   # Apply the function to each group
   pfs_res <- dt_sub %>%
     group_by(study) %>%
     do(safe_coxph(., pfs_formula))
   
   pfs_res$symbol = ifelse(pfs_res$estimate>=1,">=1","<1")
   pfs_res$study = factor(pfs_res$study,levels =pfs_res$study)
   
   names(pfs_res) = c("study","term","HR","std","statistics","pval","HR_lower","HR_upper","symbol")
   res = pfs_res
   
 }
return(res)
  
}

forest_fig = function(input,output,session,res,surv_type){
  #----------------------------------------------------------------------
  # functions to generate forest plot
  # res: univariable cox result from forest_table_server()
  # surv_type: OS or PFI
  #----------------------------------------------------------------------
  
  
  res = res
  if(surv_type=="OS"){
    xlab_name = "OS (HR)"
  }else{
    xlab_name = "PFS/PFI (HR)"
  }
  p=ggplot(res, aes(y=study, x=HR, xmin=HR_lower, xmax=HR_upper,color = symbol))+
    #Add data points and color them black
    geom_pointrange(shape=15,size = 1.5,fatten =1.5)+
    #add the CI error bars
    geom_errorbarh(height=.1)+
    ylab('')+
    xlab(xlab_name)+
    #Add a vertical dashed line indicating an effect size of zero, for reference
    geom_vline(xintercept=1, color='black', linetype='dotted',size=1)+
    theme_classic()+
    #scale_color_brewer(palette="Set1")+
    scale_color_manual(values=c("#374E55FF", "#DF8F44FF"))+
    theme(legend.position='none',axis.line = element_line(colour = "gray"),
          axis.text.y = element_text(face="bold", color="black",size=10),
          axis.text.x = element_text(face="bold", color="black"),
          axis.ticks.x = element_blank(),
          axis.title.x = element_text(size = 20)
    )
  return(p)  
}

pan_server_text = function(input_data,type){
  #----------------------------------------------------------------------
  # functions to summarize forest plot
  # input_data: dataset for the selected gene (from pan_server_text_data())
  # type: OS or PFI
  #----------------------------------------------------------------------
  
  uni_cox1 = input_data
  names_sig = uni_cox1$study[which(as.numeric(uni_cox1$pval) < 0.05)]
  
  if(length(names_sig)==0){
    text= "The prognostic marker is not significant in any cancer types"
  } else if(length(names_sig)!=0){
    names_fav = uni_cox1$study[which(uni_cox1$symbol=="<1"&uni_cox1$pval<0.05)]
    names_unfav = uni_cox1$study[which(uni_cox1$symbol==">=1"&uni_cox1$pval<0.05)]
    if(length(names_unfav)==0){
      text= paste0("The prognostic marker is significant and favourable in ",paste(names_fav,collapse=" ")," (hazard ratio < 1, low risk of death)")
    }else if(length(names_fav)==0){
      text= paste0("The prognostic marker is significant and unfavourable in ",paste(names_unfav,collapse=" ")," (hazard ratio > 1, high risk of death)")
    }else{
      text= paste0("The prognostic marker is significant and favourable in ",paste(names_fav,collapse=" ")," (hazard ratio < 1, low risk of death)",
                   ", and unfavourable in ",paste(names_unfav,collapse=" ")," (hazard ratio > 1, high risk of death)"," for ",type)
      
    }
  }else{
    text = "Too many zeros in the selected gene"
  }
  text
}

#===============================================================================
# KM plt for significant studies
#===============================================================================
# Function to generate KM plot with optimal cutoff for each group
generate_km_plot <- function(data, variable_col, study_type, cutoff) {
  #---------------------------------------
  # time_col = "OS.time"
  # event_col = "OS"
  # variable_col = sel_gene
  # study_type = unique(data$study)
  # data = dt_sub_km
  # cutoff = "median" 
  #---------------------------------------
  
  tryCatch({
    if(cutoff == "optimal"){
      res.cut <- surv_cutpoint(data, time = "OS.time", event = "OS", variables = variable_col)
      res.cat <- surv_categorize(res.cut)
      names(res.cat)[3] = "gene"
      fit <- survfit(Surv(OS.time,OS)~gene,data = res.cat)
      p <- ggsurvplot(
        fit,
        data = res.cat,
        pval = TRUE,
        legend.title = variable_col,
        legend.labs = c("high", "low"),
        palette = c("#DF8F44FF", "#374E55FF")
      )
      
      p <- p$plot +
        ggtitle(study_type) +
        xlab("Time (month)") 
    }else if(cutoff=="median"){
      dt = data
      dt$cut = findInterval(dt[[variable_col]], median(dt[[variable_col]],na.rm=T))
      dt$cut = ifelse(dt$cut==1,"high","low")
      res.cat = data.frame(cbind("time"=dt$OS.time,"status"=dt$OS,"cut"=dt$cut))
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
          p = ggsurvplot(fit,data=res.cat,pval = T,legend.title = variable_col,legend.labs=c("high","low"),palette = c("#DF8F44FF", "#374E55FF"))
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
      dt$cut = findInterval(dt[[variable_col]], quantile(dt[[variable_col]],na.rm=T))
      dt$cut = ifelse(dt$cut%in%c(2,3),NA,ifelse(dt$cut%in%c(4,5),"high","low"))
      res.cat = data.frame(cbind("time"=dt$OS.time,"status"=dt$OS,"cut"=dt$cut))
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
          
          p = ggsurvplot(fit,data=res.cat,pval = T,legend.title = variable_col,legend.labs=c("high","low"),palette = c("#DF8F44FF", "#374E55FF"))
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
    fit = survfit(Surv(OS.time,OS) ~1, data = data)
    p = ggsurvplot(fit,data=data,palette = c("#999999"))
    p= p$plot+
      ggtitle(paste0("The gene is not significant", " in ", study_type))+
      xlab("Time (month)")
  })
  
  
  return(p)
}

#------------------------------------------------------------------------------
# circular barchart 
#------------------------------------------------------------------------------
cir_bar_func = function(gene,input_dt){
  #----------------------------------------------
  # gene: input gene
  # input_dt: input dataset 
  # dt_final = dt_sub
  # gene = "CD8A"
  #----------------------------------------------
  gene = make.names(gene)
  med_group = input_dt %>%
    group_by(study) %>% 
    summarise_at(vars(gene), median)
  med_group = data.frame(na.omit(med_group))
  med_group = med_group[order(med_group[,gene],decreasing = F),]
  med_group$study = factor(med_group$study,levels = med_group$study)
  med_group$id = as.numeric(med_group$study)
  
 # med_group$idd = as.factor(1:nrow(med_group))
  med_group$idd  = med_group$study # change for circula
  
  # Get the name and the y position of each label
  label_data_md <- med_group
  number_of_bar_md <- nrow(label_data_md)
  angle_md <- 90 - 360 * (label_data_md$id-0.5) /number_of_bar_md    
  label_data_md$hjust <- ifelse( angle_md < -90, 1, 0)
  label_data_md$angle<- ifelse(angle_md < -90, angle_md+180, angle_md)
  names(label_data_md)[3]="value"
  
  val = max(med_group[[gene]])
  p = ggplot(data = med_group,aes_string(x="idd",y=gene))+
    geom_bar(stat="identity", alpha=1,fill="#ff7f0e") +
   # ylim(-(val+5),val+5)+ # change for circular
    theme_minimal() +
   # coord_polar() + # change for circular
    theme(
      #legend.position = "none", # change for circular
     # axis.text = element_blank(), # change for circular
      axis.title = element_blank(), 
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=15),# REMOVE for circular
      axis.text.y = element_text(size=15), # REMOVE for circular
      panel.grid = element_blank(),
     # plot.margin = unit(c(3, 3, 3, 3), "cm")      
     plot.margin = unit(c(3, 6, 3, 6), "cm")      
     
    )
    
  #+
  # geom_text(data=label_data_md, aes(x=idd, y=value+3.5, label=study , hjust=hjust), color="black", fontface="bold",alpha=0.6, size=4, angle= label_data_md$angle, inherit.aes = FALSE ) # change for circular 
  
  p
}