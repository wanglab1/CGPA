## getting data from GEO
library(GEOquery)
library(readr)
library(biomaRt)
### BLCA
# GSE248167

setwd("/Users/daisys/Library/CloudStorage/OneDrive-MoffittCancerCenter/Moffitt_DaisyS/Data/1_BLCA.Kidney/1_Study/1_ICI_treatment/2_NoFollowUp/GSE248167_BLCA")

myGSE<-"GSE248167"
gset <- getGEO(myGSE, GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL24676", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
cli_GSE248167<- pData(gset)
write.csv(cli_GSE248167, "GSE248167_pheno_data.csv")
########################################################
# laod expression 
## matrix 1
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE248167&format=file&file=GSE248167%5Fgene%5Ffpkm1%2Etxt%2Egz"
path <- paste(urld, "acc=GSE248167", "file=GSE248167_gene_fpkm1.txt.gz", sep="&");
ex_GSE248167_1 <- data.table::fread(path)
ex_GSE248167_1[, 1] <- ex_GSE248167_1[, 12]
ex_GSE248167_1[1:3,1:3]
write.csv(ex_GSE248167_1,"GSE248167_exp_1.csv")
# matrix 2
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE248167&format=file&file=GSE248167%5Fgene%5Ffpkm2%2Etxt%2Egz"
path <- paste(urld, "acc=GSE248167", "file=GSE248167_gene_fpkm2.txt.gz", sep="&");
ex_GSE248167_2 <- data.table::fread(path)
ex_GSE248167_2[, 1] <- ex_GSE248167_2[, 23]
ex_GSE248167_2[1:3,1:3]
write.csv(ex_GSE248167_2,"GSE248167_exp_2.csv")
## matrix 3
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE248167&format=file&file=GSE248167%5Fgene%5Ffpkm3%2Etxt%2Egz"
path <- paste(urld, "acc=GSE248167", "file=GSE248167_gene_fpkm3.txt.gz", sep="&");
ex_GSE248167_3 <- data.table::fread(path)
head(ex_GSE248167_3)
ex_GSE248167_3[, 1] <- ex_GSE248167_3[, 22]
ex_GSE248167_3[1:3,1:3]
write.csv(ex_GSE248167_3,"GSE248167_exp_3.csv")


## loaded and convert the data

ex_GSE248167_1<-read.csv("GSE248167_exp_1.csv")
ex_GSE248167_1<-t(ex_GSE248167_1)
write.table(ex_GSE248167_1, "BLCA_GSE248167_expr.1.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)

ex_GSE248167_2<-read.csv("GSE248167_exp_2.csv")
ex_GSE248167_2<-t(ex_GSE248167_2)
write.table(ex_GSE248167_2, "BLCA_GSE248167_expr.2.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)

ex_GSE248167_3<-read.csv("GSE248167_exp_3.csv")
ex_GSE248167_3<-t(ex_GSE248167_3)
write.table(ex_GSE248167_3, "BLCA_GSE248167_expr.2.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)
