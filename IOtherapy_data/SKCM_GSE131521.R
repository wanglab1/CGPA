# SKCM GSE131521
library(GEOquery)
library(readr)
library(data.table)
library(readxl)
library(dplyr)
library(IOBR)

setwd("~/Library/CloudStorage/OneDrive-MoffittCancerCenter/Moffitt_DaisyS/Data/9_SKCM/GSE131521/Archive")
#download the data 
myGSE<-"GSE131521"
gset <- getGEO(myGSE, GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL16791", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
cli<- pData(gset)
write.csv(cli, "GSE131521_pheno.csv")
## counts to TPM
expr <- read.table("GSE131521_counts.txt", header = TRUE, sep = "\t", comment.char = "!")
head(expr)
dim(expr)
expr_s<-expr[, c(1, 7:ncol(expr))]

expr_s <- expr_s[!is.na(expr$Geneid),]
expr_s <- expr_s %>%
  group_by(Geneid) %>%
  summarise_if(is.numeric, mean)
rownames(expr_s) <- expr_s$Geneid
rownames <- row.names(expr_s)
expr_s <- expr_s[, -1]
rownames(expr_s) <- rownames
expr_tpm<-count2tpm(countMat = expr_s, idType = "Ensembl", org="hsa", source = "local" )
dim(expr_tpm)
expr_tpm[1:5,1:5]
#[1] 49171    17
# add rowname
expr_tpm<-data.frame(GeneSymbol = rownames(expr_tpm), expr_tpm, check.names = FALSE)
expr_tpm[1:5,1:5]


## convert
expr_t<-t(expr_tpm)
expr_t[1:5,1:5]
write.table(expr_t, "GSE131521_temp.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)
expr_tpm <- read_tsv("GSE131521_temp.tsv")
expr_tpm[1:3,1:3]
cli<-read.csv("GSE131521_pheno_clean.csv")
head(cli)
colnames(expr_tpm)[1] <- "sampleID"
# merge
merged_data <- merge(cli,expr_tpm,by = "sampleID")
merged_data[1:3,1:3]
dim(merged_data)
# [1]    17 49178
write.csv(merged_data, "SKCM_GSE131521.csv")

