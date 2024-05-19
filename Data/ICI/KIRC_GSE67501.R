## getting data from GEO
library(GEOquery)
library(readr)
library(biomaRt)
### KIRC
# GSE67501

setwd("/Users/daisys/Library/CloudStorage/OneDrive-MoffittCancerCenter/Moffitt_DaisyS/Data/10_KIRC/GSE67501")

myGSE<-"GSE67501"
gset <- getGEO(myGSE, GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL14951", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
cli_GSE67501<- pData(gset)
write.csv(cli_GSE67501, "GSE67501_pheno.csv")
expr<-exprs(gset)
head(expr)
write.csv(expr, "GSE67501_expr.csv")
########################################################
cli<-read.csv("GSE67501_pheno_clean.csv")
expr<-read.csv("GSE67501_expr.csv")
expr_t<-t(expr)
write.table(expr_t, "expr_t_tmp.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)
expr <- read_tsv("expr_t_tmp.tsv")
expr[1:3,1:3]
cli[1:3,1:3]
colnames(expr)[1] <- "geo_accession"
## orgnize to data matrix
merged_data <- merge(cli,expr,by = "geo_accession")
merged_data[1:3,1:3]
dim(merged_data)
# [1]   11 29384
write.csv(merged_data,"KIRC_GSE67501.csv")
