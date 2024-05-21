# BRCA GSE173839
library(GEOquery)
library(readr)
library(data.table)
library(readxl)
library(dplyr)
setwd("~/Library/CloudStorage/OneDrive-MoffittCancerCenter/Moffitt_DaisyS/Data/4_BRCA/GSE173839")
#download the data 
myGSE<-"GSE173839"
gpl <- getGEO("GPL20078")
gset <- getGEO(myGSE, GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL20078", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
cli<- pData(gset)
write.csv(cli_GSE173839, "GSE173839_pheno.csv")
expr<- exprs(gset)
dim(expr)
expr[1:5,1:5]
annotation_data <- Table(gpl)[, c("ID", "GeneName")]
# mapping the gene name
aligned_data <- inner_join(data.frame(ProbeID = rownames(expr), expr, check.names = FALSE),
                           annotation_data, by = c("ProbeID" = "ID"))
dim(aligned_data) 
aligned_data[, 1] <- aligned_data[, 107]
aligned_data[1:5,1:5]

## convert
aligned_data_t<-t(aligned_data)
aligned_data_t[1:5,1:5]
write.table(aligned_data_t, "aligned_data_temp.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)
expr <- read_tsv("aligned_data_temp.tsv")
expr[1:3,1:3]
colnames(expr)[1] <- "geo_accession"
cli<-read.csv("GSE173839_pheno.csv")
# merge
merged_data <- merge(cli,expr,by = "geo_accession")
merged_data[1:3,1:3]
dim(merged_data)
# 105 32153
write.csv(merged_data, "BRCA_GSE173839.csv")

