# CTCL GSE162137
library(GEOquery)
library(readr)
library(data.table)
library(readxl)
library(dplyr)
library(IOBR)

setwd("~/Library/CloudStorage/OneDrive-MoffittCancerCenter/Moffitt_DaisyS/Data/11_CTCL/GSE162137/Archive")
#download the data 
myGSE<-"GSE162137"
gset <- getGEO(myGSE, GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL18573", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
cli<- pData(gset)
write.csv(cli, "GSE162137_pheno.csv")
## counts to TPM
expr <- read.csv("GSE162137_CTCL_raw_reads.csv")
head(expr)
dim(expr)
expr <- expr[!is.na(expr$gene.ID),]
expr <- expr %>%
  group_by(gene.ID) %>%
  summarise_if(is.numeric, mean)
rownames(expr) <- expr$gene.ID
rownames <- row.names(expr)
expr <- expr[, -1]
rownames(expr) <- rownames
expr_tpm<-count2tpm(countMat = expr, idType = "Ensembl", org="hsa", source = "local" )
dim(expr_tpm)
expr_tpm[1:5,1:5]
#[1] 56492    64
# add rowname
expr_tpm<-data.frame(GeneSymbol = rownames(expr_tpm), expr_tpm, check.names = FALSE)
expr_tpm[1:5,1:5]


## convert
expr_t<-t(expr_tpm)
expr_t[1:5,1:5]
write.table(expr_t, "GSE162137_temp.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)
expr_tpm <- read_tsv("GSE162137_temp.tsv")
expr_tpm[1:3,1:3]
cli<-read.csv("GSE162137_pheno_clean.csv")
head(cli)
colnames(expr_tpm)[1] <- "sampleID"
# merge
merged_data <- merge(cli,expr_tpm,by = "sampleID")
merged_data[1:3,1:3]
dim(merged_data)
#[1]    64 56496
write.csv(merged_data, "CTCL_GSE162137.csv")

