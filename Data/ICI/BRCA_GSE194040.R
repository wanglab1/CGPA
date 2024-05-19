# BRCA GSE194040
library(GEOquery)
library(readr)
library(data.table)
library(readxl)
library(dplyr)
setwd("~/Library/CloudStorage/OneDrive-MoffittCancerCenter/Moffitt_DaisyS/Data/4_BRCA/GSE194040")
#download the data 
myGSE<-"GSE194040"
gpl <- getGEO("GPL20078")
gset <- getGEO(myGSE, GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL20078", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
cli<- pData(gset)
write.csv(cli, "GSE194040_pheno.csv")
## get the matrix
expr <- read.table("GSE194040_ISPY2ResID_AgilentGeneExp_990_FrshFrzn_GPL20078_ProbeLevel_n654.txt.gz", sep = "\t")
dim(expr)
expr[1:5,1:5]
annotation_data <- Table(gpl)[, c("ID", "GeneName")]
# mapping the gene name
aligned_data <- inner_join(data.frame(ProbeID = rownames(expr), expr, check.names = FALSE),
                           annotation_data, by = c("ProbeID" = "ID"))
dim(aligned_data) 
aligned_data[, 1] <- aligned_data[, 656]
aligned_data[1:5,1:5]
## convert
aligned_data_t<-t(aligned_data)
aligned_data_t[1:5,1:5]
write.table(aligned_data_t, "aligned_data_temp.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)
expr <- read_tsv("aligned_data_temp.tsv")
expr[1:3,1:3]
colnames(expr)[1] <- "patient_id" 
cli<-read.csv("GSE194040_pheno_clean.csv")
# merge
merged_data <- merge(cli,expr,by = "patient_id")
merged_data[1:3,1:3]
colnames(merged_data)[1] <- "geo_accession"
dim(merged_data)
# 653 32155
write.csv(merged_data, "BRCA_GSE194040.csv")

