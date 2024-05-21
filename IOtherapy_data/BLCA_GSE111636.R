## getting data from GEO
library(GEOquery)
library(readr)
library(biomaRt)
library(annotate)
library(hta20transcriptcluster.db)
library(dplyr)

### BLCA
# GSE111636

setwd("/Users/daisys/Library/CloudStorage/OneDrive-MoffittCancerCenter/Moffitt_DaisyS/Data/1_BLCA.Kidney/1_Study/1_ICI_treatment/2_NoFollowUp/GSE111636_response_BLCA")

myGSE<-"GSE111636"
gset <- getGEO(myGSE, GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL17586", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
cli_GSE111636<- pData(gset)
ex_GSE111636 <- exprs(gset)

write.csv(cli_GSE111636, "GSE111636_pheno_data.csv")
write.csv(ex_GSE111636, "GSE111636_exp.csv")

## loaded the cleaned clinical data processed at background
cli_GSE111636<-read.csv("GSE111636_pheno_clean.csv")
ex_GSE111636<-read.csv("GSE111636_exp.csv")
ex_GSE111636_t<-t(ex_GSE111636)
write.table(ex_GSE111636_t, "ex_GSE111636_t.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)
ex_GSE111636 <- read_tsv("ex_GSE111636_t.tsv")
ex_GSE111636[1:3,1:3]
cli_GSE111636[1:3,1:3]
## orgnize to data matrix
merged_data <- merge(cli_GSE111636,ex_GSE111636,by = "geo_accession")
merged_data[1:3,1:3]
dim(merged_data)
# [1]    11 32650
write.csv(merged_data,"BLCA_GSE111636.csv")
