# SKCM GSE91061
library(GEOquery)
library(readr)
library(data.table)
library(readxl)
library(dplyr)
library(IOBR)
library(biomaRt)
setwd("~/Library/CloudStorage/OneDrive-MoffittCancerCenter/Moffitt_DaisyS/Data/9_SKCM/GSE91061")
#download the data 
myGSE<-"GSE91061"
gset <- getGEO(myGSE, GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL9052", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
cli<- pData(gset)
write.csv(cli, "GSE91061_pheno.csv")
# convert expr file ID to gene symbol
expr<-read.csv("GSE91061_expr.csv")
dim(expr)
ncbi_gene_ids <-as.character(expr[[1]])
gene_symbols <- select(org.Hs.eg.db, keys=ncbi_gene_ids, columns=c("SYMBOL"), keytype="ENTREZID")
dim(gene_symbols)
joined_data <- merge(expr, gene_symbols, by = "ENTREZID")
head(joined_data)
dim(joined_data)
joined_data[, 1] <- joined_data[, 111]
head(joined_data)
joined_data_new<-joined_data[, -111]
colnames(joined_data_new)[1] <- "GeneSymbol" 
head(joined_data_new)
## convert
joined_data_new_t<-t(joined_data_new)
joined_data_new_t[1:5,1:5]
write.table(joined_data_new_t, "GSE91061_temp.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)
expr_tpm <- read_tsv("GSE91061_temp.tsv")
expr_tpm[1:3,1:3]
cli_1<-read.csv("GSE91061_pheno_clean.csv")
head(cli_1)
cli_2<-read.csv("mmc2.csv")
head(cli_2)
# merge
merged_data <- merge(cli_1,cli_2,by = "Patient")
merged_data[1:3,1:3]
dim(merged_data)
# [1]    28 25281
write.csv(merged_data, "merged_cli.csv")
cli<-read.csv("merged_cli.csv")
colnames(expr_tpm)[1] <- "sampleID"
merged<- merge(cli,expr_tpm,by = "sampleID")
write.csv(merged, "SKCM_GSE91061.csv")
