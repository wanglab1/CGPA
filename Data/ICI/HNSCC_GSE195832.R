# HNSCC GSE195832
library(GEOquery)
library(readr)
library(data.table)
library(readxl)
library(dplyr)
library(IOBR)
setwd("/Users/daisys/Library/CloudStorage/OneDrive-MoffittCancerCenter/Moffitt_DaisyS/Data/8_HNSCC/GSE195832")
#download the data 
myGSE<-"GSE195832"
#gpl <- getGEO("GPL24676")
#head(Table(gpl))
gset <- getGEO(myGSE, GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL24676", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
cli<- pData(gset)
write.csv(cli, "GSE195832_pheno.csv")
## counts to TPM
expr <- read_tsv("GSE195832_TJU_featurecounts_symbol.tsv.gz")
head(expr)
expr <- expr[!is.na(expr$symbol),]
expr <- expr %>%
  group_by(symbol) %>%
  summarise_if(is.numeric, mean)
rownames(expr) <- expr$symbol
rownames <- row.names(expr)
expr <- expr[, -1]
rownames(expr) <- rownames
expr<-count2tpm(countMat = expr, idType = "symbol", org="hsa", source = "local" )
dim(expr)
expr[1:5,1:5]
# [1] 49434    56
# add rowname
expr<-data.frame(GeneSymbol = rownames(expr), expr, check.names = FALSE)
expr[1:5,1:5]
## convert
expr_t<-t(expr)
expr_t[1:5,1:5]
write.table(expr_t, "GSE195832_temp.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)
expr_tpm <- read_tsv("GSE195832_temp.tsv")
expr_tpm[1:3,1:3]
cli<-read.csv("metadata_TJcohort.csv")
head(cli)
colnames(expr_tpm)[1] <- "Sample"
# merge
merged_data <- merge(cli,expr_tpm,by = "Sample")
merged_data[1:3,1:3]
dim(merged_data)
# [1]    56 49443
write.csv(merged_data, "HNSCC_GSE195832.csv")

