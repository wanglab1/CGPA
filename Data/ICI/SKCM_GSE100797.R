# SKCM GSE100797
library(GEOquery)
library(readr)
library(data.table)
library(readxl)
library(dplyr)
library(IOBR)
library(biomaRt)
setwd("~/Library/CloudStorage/OneDrive-MoffittCancerCenter/Moffitt_DaisyS/Data/9_SKCM/GSE100797")
#download the data 
myGSE<-"GSE100797"
gset <- getGEO(myGSE, GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL11154", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
cli<- pData(gset)
write.csv(cli, "GSE100797_pheno.csv")
# convert expr file ID to gene symbol
cli<-read.csv("GSE100797_pheno_clean.csv")
expr<-read.csv("GSE100797_ProcessedData.csv")
head(expr)
## convert
expr_t<-t(expr)
expr_t[1:5,1:5]
write.table(expr_t, "GSE100797_temp.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)
expr_tpm <- read_tsv("GSE100797_temp.tsv")
expr_tpm[1:3,1:3]
# merge
merged_data <- merge(cli,expr_tpm,by = "sampleID")
merged_data[1:3,1:3]
dim(merged_data)
# [1]     25 18430
write.csv(merged_data, "SKCM_GSE100797.csv")
