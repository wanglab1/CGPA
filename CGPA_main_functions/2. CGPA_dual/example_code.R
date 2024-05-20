#===============================================================================
# Gene-dual
#===============================================================================
library(forestplot)
library(survival)
library(fst)
library(broom)
library(expss)
library(forestmodel)

genes = c("EGFR","CD8A")
geneName = read_fst("C:/Users/4467777/Desktop/TCGA_shiny/CGPA_dual/www/geneName.fst")
folder_direct = "C:/Users/4467777/Desktop/TCGA_shiny/cgpa_single031222/"


#------------------------------------------------------------------------------
########################## Univariate cox #################################
#------------------------------------------------------------------------------

dt = lapply(1:2,function(x){
  if(length(which(geneName$geneName==genes[x]))>1){
    dt = read_fst(paste0(folder_direct,"server_result/surv_data/","HNSC","_","raw",".fst"),
                  columns = c("OS","OS.time","PFI","PFI.time",paste0(make.names(genes[x]),".x"))) 
    names(dt)[ncol(dt)] = make.names(genes[x])
  }else{
    dt = read_fst(paste0(folder_direct,"server_result/surv_data/","HNSC","_","raw",".fst"),
                  columns = c("OS","OS.time","PFI","PFI.time",make.names(genes[x])))    
    names(dt)[ncol(dt)] = make.names(genes[x])
  }
  dt
})

dt = data.frame(do.call(cbind,dt))
dt[,genes] = log2(dt[,genes]+1)

uni_cox = lapply(1:2,function(x){
  fit <- eval(parse(text = paste0("coxph(Surv(OS.time,OS) ~ ", make.names(genes[x]), ", data = dt)")))
  
  fit = tidy(fit,exponentiate = T,conf.int=T)
  fit
})
uni_cox = do.call(rbind,uni_cox)
uni_cox$symbol = ifelse(uni_cox$estimate>1,">1","<=1")

p1=ggplot(uni_cox, aes(y=term, x=estimate, xmin=conf.low, xmax=conf.high,color = symbol))+
  #Add data points and color them black
  geom_pointrange(shape=15,size = 1.5,fatten =1.5)+
  #add the CI error bars
  geom_errorbarh(height=.1)+
  #Specify the limits of the x-axis and relabel it to something more meaningful
  #  scale_x_continuous(limits=c(-2,2), name='Standardized Mean Difference (d)')+
  #Give y-axis a meaningful label
  ylab('')+
  xlab("OS (HR)")+
  #Add a vertical dashed line indicating an effect size of zero, for reference
  geom_vline(xintercept=1, color='black', linetype='dotted',size=1)+
  theme_classic()+
  #scale_color_brewer(palette="Set1")+
  scale_color_manual(values=c("#374E55FF", "#DF8F44FF"))+
  theme(legend.position='none',axis.line = element_line(colour = "gray"),
        axis.text.y = element_text(face="bold", color="gray"),
        axis.text.x = element_text(face="bold", color="gray"),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 18)
        
        
  )
ggplotly(p1)

#------------------------------------------------------------------------------
########################## Univariate KM #################################
#------------------------------------------------------------------------------
dt = lapply(1:2,function(x){
  if(length(which(geneName$geneName==genes[x]))>1){
    dt = read_fst(paste0(folder_direct,"server_result/surv_data/","HNSC","_",tolower("Median"),".fst"),
                  columns = c("OS","OS.time","PFI","PFI.time",paste0(make.names(genes[x]),".x"))) 
    names(dt)[ncol(dt)] = make.names(genes[x])
  }else{
    dt = read_fst(paste0(folder_direct,"server_result/surv_data/","HNSC","_",tolower("Median"),".fst"),
                  columns = c("OS","OS.time","PFI","PFI.time",make.names(genes[x])))    
    names(dt)[ncol(dt)] = make.names(genes[x])
  }
  dt
})

dt = data.frame(do.call(cbind,dt))
dt

dt$OS.time = dt$OS.time/30.417
dt$PFI.time= dt$PFI.time/30.417

km = lapply(1:2,function(x){
  fit <- eval(parse(text = paste0("survfit(Surv(OS.time,OS) ~ ", make.names(genes[x]), ", data = dt)")))
  p = ggsurvplot(fit,data = dt,pval=T,legend.title=genes[x],legend.labs=c("high","low"),palette = c("#DF8F44FF","#374E55FF"),size=0.5)
  p
})
km
ggarrange(km[[1]]$plot,km[[2]]$plot,ncol=2)

# KM table
OS_KM_fit_table <- purrr::map_dfr( make.names(genes), run_pretty_km_output, model_data = dt, time_in = 'OS.time', 
                                  event_in = 'OS', event_level = '1') %>% 
                  select(Group, Level, everything())
OS_KM_fit_table  
  
#------------------------------------------------------------------------------
########################## Gene ratio cox #################################
#------------------------------------------------------------------------------
dt = lapply(1:2,function(x){
  if(length(which(geneName$geneName==genes[x]))>1){
    dt = read_fst(paste0(folder_direct,"server_result/surv_data/","HNSC","_","raw",".fst"),
                  columns = c("OS","OS.time","PFI","PFI.time",paste0(make.names(genes[x]),".x"))) 
    names(dt)[ncol(dt)] = make.names(genes[x])
  }else{
    dt = read_fst(paste0(folder_direct,"server_result/surv_data/","HNSC","_","raw",".fst"),
                  columns = c("OS","OS.time","PFI","PFI.time",make.names(genes[x])))    
    names(dt)[ncol(dt)] = make.names(genes[x])
  }
  dt
})


dt = data.frame(do.call(cbind,dt))
dt$ratio = dt[,genes[1]]/dt[,genes[2]]
dt$ratio = log2(dt$ratio+1)
dt[,genes] = log2(dt[,genes]+1)
dt = apply_labels(dt,ratio = paste0(make.names(genes[1]),"/",make.names(genes[2])))

fit = coxph(Surv(OS.time,OS)~ratio,data = dt)
summary(fit)
forest_model(fit)


#------------------------------------------------------------------------------
########################## Gene ratio KM #################################
#------------------------------------------------------------------------------
cutoff = "median"
dt = lapply(1:2,function(x){
  if(length(which(geneName$geneName==genes[x]))>1){
    dt = read_fst(paste0(folder_direct,"server_result/surv_data/","HNSC","_","raw",".fst"),
                  columns = c("OS","OS.time","PFI","PFI.time",paste0(make.names(genes[x]),".x"))) 
    names(dt)[ncol(dt)] = make.names(genes[x])
  }else{
    dt = read_fst(paste0(folder_direct,"server_result/surv_data/","HNSC","_","raw",".fst"),
                  columns = c("OS","OS.time","PFI","PFI.time",make.names(genes[x])))    
    names(dt)[ncol(dt)] = make.names(genes[x])
  }
  dt
})
dt = data.frame(do.call(cbind,dt))
dt$ratio = dt[,genes[1]]/dt[,genes[2]]
dt$ratio = log2(dt$ratio+1)
dt$OS.time = dt$OS.time/30.417
dt$PFI.time= dt$PFI.time/30.417

if(cutoff=="median"){
  dt$ratio_group =  findInterval(dt$ratio,median(dt$ratio,na.rm=T))   
  dt$ratio_group = ifelse(dt$ratio_group==1,"high","low")

} else if(cutoff=="quartile"){
  dt$ratio_group =  findInterval(dt$ratio,quantile(dt$ratio,na.rm=T))   
  dt$ratio_group = ifelse(dt$ratio_group%in%c(2,3),NA,ifelse(dt$ratio_group%in%c(4,5),"high","low"))
  
} else if(cutoff=="optimal"){
  res.cut <- surv_cutpoint(dt,time = "OS.time", event = "OS", variables="ratio")$cutpoint[[1]]
  dt$ratio_group  = ifelse(dt$ratio >=res.cut ,"high","low") 
}

dt = apply_labels(dt,ratio_group = paste0(make.names(genes[1]),"/",make.names(genes[2])))


fit <- eval(parse(text = paste0("survfit(Surv(OS.time,OS) ~ ", "ratio_group", ", data = dt)")))
  p = ggsurvplot(fit,data = dt,pval=T,legend.title=paste0(genes[1],"/",genes[2]),legend.labs=c("high","low"),palette = c("#DF8F44FF","#374E55FF"),size=0.5)
  p


# KM table
OS_KM_fit_table <- purrr::map_dfr("ratio_group", run_pretty_km_output, model_data = dt, time_in = 'OS.time', 
                                   event_in = 'OS', event_level = '1') %>% 
  select(Group, Level, everything())
OS_KM_fit_table  

#------------------------------------------------------------------------------
########################## multivariate cox #################################
#------------------------------------------------------------------------------
dt = lapply(1:2,function(x){
  if(length(which(geneName$geneName==genes[x]))>1){
    dt = read_fst(paste0(folder_direct,"server_result/surv_data/","HNSC","_","raw",".fst"),
                  columns = c("OS","OS.time","PFI","PFI.time",paste0(make.names(genes[x]),".x"))) 
    names(dt)[ncol(dt)] = make.names(genes[x])
  }else{
    dt = read_fst(paste0(folder_direct,"server_result/surv_data/","HNSC","_","raw",".fst"),
                  columns = c("OS","OS.time","PFI","PFI.time",make.names(genes[x])))    
    names(dt)[ncol(dt)] = make.names(genes[x])
  }
  dt
})
dt = data.frame(do.call(cbind,dt))
dt[,genes] = log2(dt[,genes]+1)
dt$OS.time = dt$OS.time/30.417
dt$PFI.time= dt$PFI.time/30.417

formula = as.formula(paste0("Surv(OS.time,OS)~",paste(genes,collapse = "+")))
fit = coxph(formula,data = dt)
forest_model(fit)

#-------------------------------------------------------------------------------
################ KM with 4 groups ########################
#-------------------------------------------------------------------------------
genes = c("CD3E","HAVCR2")
dt = lapply(1:2,function(x){
  if(length(which(geneName$geneName==genes[x]))>1){
    dt = read_fst(paste0(folder_direct,"server_result/surv_data/","HNSC","_",tolower("quartile"),".fst"),
                  columns = c("OS","OS.time","PFI","PFI.time",paste0(make.names(genes[x]),".x"))) 
    names(dt)[ncol(dt)] = make.names(genes[x])
  }else{
    dt = read_fst(paste0(folder_direct,"server_result/surv_data/","HNSC","_",tolower("quartile"),".fst"),
                  columns = c("OS","OS.time","PFI","PFI.time",make.names(genes[x])))    
    names(dt)[ncol(dt)] = make.names(genes[x])
  }
  dt
})

dt = data.frame(do.call(cbind,dt))
dt

dt$OS.time = dt$OS.time/30.417
dt$PFI.time= dt$PFI.time/30.417
dt$gene_combine = paste0(dt[[genes[1]]],"+",dt[[genes[2]]])
dt = na.omit(dt)

formula = as.formula("Surv(OS.time,OS)~gene_combine")
fit = survfit(formula,data = dt)
fit$call$formula <- formula
p = ggsurvplot(fit,data = dt,pval=T,legend.title=paste0(genes[1],"+",genes[2]),legend.labs=c("high+high","high+low","low+high","low+low"),palette = c("#DF8F44FF","gray","#4682B433","#374E55FF"))
p


# KM table

OS_KM_fit_table <- purrr::map_dfr("gene_combine", run_pretty_km_output, model_data = dt, time_in = 'OS.time', 
                                  event_in = 'OS', event_level = '1') %>% 
  select(Group, Level, everything())
OS_KM_fit_table  

# pairwise logrank
formula = as.formula("Surv(OS.time,OS)~gene_combine")
res = pairwise_survdiff(formula,data=dt)
res = tidy(res)
res$p.value = round(res$p.value,3)
res
data.table(data.frame(tidy(res)))

#-------------------------------------------------------------------------------
###################### Interaction model ####################################
#-------------------------------------------------------------------------------
######## Interaction cox #############
dt = lapply(1:2,function(x){
  if(length(which(geneName$geneName==genes[x]))>1){
    dt = read_fst(paste0(folder_direct,"server_result/surv_data/","HNSC","_","raw",".fst"),
                  columns = c("OS","OS.time","PFI","PFI.time",paste0(make.names(genes[x]),".x"))) 
    names(dt)[ncol(dt)] = make.names(genes[x])
  }else{
    dt = read_fst(paste0(folder_direct,"server_result/surv_data/","HNSC","_","raw",".fst"),
                  columns = c("OS","OS.time","PFI","PFI.time",make.names(genes[x])))    
    names(dt)[ncol(dt)] = make.names(genes[x])
  }
  dt
})
dt = data.frame(do.call(cbind,dt))
dt[,genes] = log2(dt[,genes]+1)
dt$OS.time = dt$OS.time/30.417
dt$PFI.time= dt$PFI.time/30.417

formula = as.formula(paste0("Surv(OS.time,OS)~",paste(genes,collapse = "*")))
cox_fit = coxph(formula,data=dt)
res = tidy(cox_fit,exponentiate = T,conf.int = T)
res = data.frame(res[,c("term","estimate","conf.low","conf.high","p.value")])

res$p.value = ifelse(res$p.value>=0.001,round(res$p.value,3),
                     formatC(res$p.value,format = "e",digits=3))

for(i in 2:4){
  res[,i] = round(res[,i],3)
}
res
names(res) = c("term","HR","conf.low", "conf.high", "p.value")
res

######### interaction KM #############
genes = c("CD3E","HAVCR2")
  if(length(which(geneName$geneName==genes[1]))>1){
    dt = read_fst(paste0(folder_direct,"server_result/surv_data/","HNSC","_",tolower("quartile"),".fst"),
                  columns = c("OS","OS.time","PFI","PFI.time",paste0(make.names(genes[1]),".x"))) 
    names(dt)[ncol(dt)] = make.names(genes[1])
  }else{
    dt = read_fst(paste0(folder_direct,"server_result/surv_data/","HNSC","_",tolower("quartile"),".fst"),
                  columns = c("OS","OS.time","PFI","PFI.time",make.names(genes[1])))    
    names(dt)[ncol(dt)] = make.names(genes[1])
  }
  dt

dt$OS.time = dt$OS.time/30.417
dt$PFI.time= dt$PFI.time/30.417

if(length(which(geneName$geneName==genes[2]))>1){
  dt_raw = read_fst(paste0(folder_direct,"server_result/surv_data/","HNSC","_",tolower("raw"),".fst"),
                columns = c(paste0(make.names(genes[2]),".x"))) 
  names(dt_raw)[ncol(dt_raw)] = make.names(genes[2])
}else{
  dt_raw = read_fst(paste0(folder_direct,"server_result/surv_data/","HNSC","_",tolower("raw"),".fst"),
                columns = c(make.names(genes[2])))    
  names(dt_raw)[ncol(dt_raw)] = make.names(genes[2])
}
dt_raw

dt = data.frame(cbind(dt,dt_raw))       


dt1 = subset(dt,dt[[genes[1]]]=="high")

cutoff="median"
if(cutoff=="median"){
  dt1$cut = findInterval(dt1[[genes[2]]], median(dt1[[genes[2]]],na.rm=T))
  dt1$cut = ifelse(dt1$cut==1,"high","low")
}else if(cutoff=="optimal"){
  res.cut <- surv_cutpoint(dt1,time = "OS.time", event = "OS", variables=genes[2])
  dt1 <- surv_categorize(res.cut)
  named(dt1)[3]="cut"
}else{
  dt1$cut = findInterval(dt1[[genes[2]]], quantile(dt1[[genes[2]]],na.rm=T))
  dt1$cut = ifelse(dt1$cut%in%c(2,3),NA,ifelse(dt1$cut%in%c(4,5),"high","low"))
}

dt2 = subset(dt,dt[[genes[1]]]=="low")

cutoff="median"
if(cutoff=="median"){
  dt2$cut = findInterval(dt2[[genes[2]]], median(dt2[[genes[2]]],na.rm=T))
  dt2$cut = ifelse(dt2$cut==1,"high","low")
}else if(cutoff=="optimal"){
  dt2 <- surv_cutpoint(dt2,time = "OS.time", event = "OS", variables=genes[2])
  dt2 <- surv_categorize(res.cut)
}else{
  dt2$cut = findInterval(dt2[[genes[2]]], quantile(dt2[[genes[2]]],na.rm=T))
  dt2$cut = ifelse(dt2$cut%in%c(2,3),NA,ifelse(dt2$cut%in%c(4,5),"high","low"))
}

fit1 = survfit(Surv(OS.time,OS)~cut,data = dt1)
fit2 = survfit(Surv(OS.time,OS)~cut,data = dt2)

p1 = ggsurvplot(fit1,data = dt1,pval=T,legend.title=genes[2],legend.labs=c("high","low"),palette = c("#DF8F44FF","#374E55FF")) + ggtitle(paste0(genes[1],":High"))
p2 = ggsurvplot(fit2,data = dt2,pval=T,legend.title=genes[2],legend.labs=c("high","low"),palette = c("#DF8F44FF","#374E55FF"))+ ggtitle(paste0(genes[1],":Low"))
ggarrange(p1$plot,p2$plot)
