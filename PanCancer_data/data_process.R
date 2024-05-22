#========================================
# Data process for mRNA PanCancer Atlas
#========================================

library(data.table)
library(readxl)
library(fst)
library(tidyverse)
#----------------------------------------------------------------------------------
########################## TCGA mRNA ###########################################
#---------------------------------------------------------------------------------
# Data can be downloaded from 
# https://api.gdc.cancer.gov/data/3586c0da-64d0-4b74-a449-5ff4d9136611
dt<- fread("EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv")
################ Remove normal samples # 20531x11070######################
names(dt)=substr(names(dt),1,16)
names_temp = substr(names(dt),14,16)
# 10,11,12,13,14 refer to normal samples
#which(grepl("10",names_temp))
#tumor_num = which(!(grepl("11",names_temp)))
# only keep 01 and 06
tumor_num1 = which((grepl("01",names_temp)))
tumor_num2 = which((grepl("06",names_temp)))
tumor_num = c(tumor_num1,tumor_num2)
# only keep 01 and 06 
dt= data.frame(dt)
#grepl("12",names_temp)
#grepl("13",names_temp)
#grepl("14",names_temp)
dt = data.frame(dt[,c(1,tumor_num)]) # 20531x10101

############## Exclude any gene with ? in genenames ##############
include_number = which(!(grepl("\\?",dt$gene_id)))
dt = data.frame(dt[include_number,]) # 20502x10333
dt$gene_id = sapply(1:nrow(dt),function(x){
  strsplit(dt$gene_id[x],"\\|")[[1]][1]
  
})
#subset duplicated rows
dt_sub_dup=dt[duplicated(dt$gene_id)|duplicated(dt$gene_id, fromLast=TRUE),]
dt_sub_dup=dt_sub_dup[order(dt_sub_dup$gene_id),]
dt_sub_dup$row_sum=apply(dt_sub_dup[,-c(1)],1,sum)
dt_sub_dup <- dt_sub_dup[order(dt_sub_dup$gene_id, -abs(dt_sub_dup$row_sum) ), ] ### sort first
dt_sub_dup <- dt_sub_dup[ !duplicated(dt_sub_dup$gene_id), ]  ### Keep highest
# Remove dplicates rows
dt=dt[!(duplicated(dt$gene_id) | duplicated(dt$gene_id, fromLast = TRUE)), ]
t=as.numeric(dim(dt_sub_dup)[2])
dt_sub_dup=as.data.frame(dt_sub_dup)
colnames(dt_sub_dup)=c(colnames(dt),"rowsum")
dt=data.frame(data.frame(rbind(dt,dt_sub_dup[,-t]))) # 20501x10101


geneName = data.frame(dt$gene_id)
write.fst(geneName)
names(geneName) = "geneName"
write.fst(geneName,"geneName_mrna.fst")

rownames(dt) = dt$gene_id
dt = data.frame(dt[,-1])
save(dt,file = "mRNA.RData")


rownames(dt) = dt$gene_id
dt = data.frame(dt[,-1])
save(dt,file = "mRNA.RData")

######################### merge with phenotype data ###########################
load("clin.RData")

dt = data.frame(t(dt),check.names  = F)
names_dup= substr(rownames(dt),1,12)
check_samples = names_dup[duplicated(names_dup)]
dup_samples_id = rownames(dt)[which(names_dup%in%check_samples)]
uniq_sample_id =  rownames(dt)[which(!(names_dup%in%check_samples))]
# Convert the vector to a data frame for easier manipulation
df <- data.frame(dup_samples_id, stringsAsFactors = FALSE)

# Extract relevant parts for sorting and filtering
df <- df %>%
  mutate(
    base = substr(dup_samples_id, 1, 12),  
    type = case_when(
      grepl("\\.01A$", dup_samples_id) ~ "01A",
      grepl("\\.01$", dup_samples_id) ~ "01",
      grepl("\\.06$", dup_samples_id) ~ "06",
      TRUE ~ "other"
    )
  )

# Function to prioritize and filter rows
prioritize_rows <- function(df) {
  df %>%
    arrange(base, match(type, c("01A", "01", "06", "other"))) %>%
    group_by(base) %>%
    slice(1) %>%
    ungroup()
}

df_filtered <- prioritize_rows(df)

dup_rows =  df_filtered$dup_samples_id

omics = data.frame(dt[rownames(dt) %in% c(uniq_sample_id,dup_rows), ],check.names=F)
rownames(omics) =  substr(rownames(omics),1,12)

mRNA_clin = merge(clin,omics,by.x="Row.names",by.y="row.names")
save(mRNA_clin,file = "mRNA.RData")

#########################  merge with purity data ######################### 
# Import all the purity data
setwd("../purity/")
temp = list.files(pattern="*.txt")
cancers = sapply(1:length(temp),function(x) {
  t = strsplit(temp[x],"-")[[1]][2]  
  strsplit(t,"\\.")[[1]][1]
})
myfiles = lapply(temp, read.delim)
purity = do.call(rbind,myfiles)
purity = data.frame(purity[,c("sampleid","purity")])
purity$sampleid_temp = substr(purity$sampleid,1,12)
purity$sampleid[which(duplicated(purity$sampleid_temp))]
purity$sampleid[which(duplicated(purity$sampleid_temp,fromLast=T))]
# remove duplicates 
purity = data.frame(purity[-which(duplicated(purity$sampleid_temp)),])
names_temp = substr(purity$sampleid,14,16)
purity$sampleid = purity$sampleid_temp
purity = data.frame(purity[,c("sampleid","purity")]) # n = 8566

load("mRNA_clin.RData") # n =10069
dt = merge(purity,dt,by.x="sampleid",by.y ="Row.names",all.y=TRUE)
names(dt)[1]="Row.names"
save(dt,file = "/mRNA_clin.RData")

