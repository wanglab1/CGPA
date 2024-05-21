# SKCM GSE78220
library(GEOquery)
library(readr)
library(data.table)
library(readxl)
library(dplyr)
library(IOBR)

setwd("~/Library/CloudStorage/OneDrive-MoffittCancerCenter/Moffitt_DaisyS/Data/9_SKCM/GSE78220")
#download the data 
myGSE<-"GSE78220"
#gpl <- getGEO("GPL11154")
#head(Table(gpl))
gset <- getGEO(myGSE, GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL11154", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
cli<- pData(gset)
write.csv(cli, "GSE78220_pheno.csv")
# read matrix
expr<-read.csv("GSE78220_PatientFPKM.csv")
expr[1:5,1:5]
dim(expr)

## convert
expr_t<-t(expr)
expr_t[1:5,1:5]
write.table(expr_t, "GSE78220_temp.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)
expr_tpm <- read_tsv("GSE78220_temp.tsv")
expr_tpm[1:3,1:3]
cli<-read.csv("GSE78220_pheno_clean.csv")
head(cli)
# merge
merged_data <- merge(cli,expr_tpm,by = "patient_id")
merged_data[1:3,1:3]
dim(merged_data)
# [1]    28 25281
write.csv(merged_data, "SKCM_GSE78220.csv")

