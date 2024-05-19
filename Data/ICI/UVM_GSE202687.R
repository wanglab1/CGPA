## getting data from GEO
library(GEOquery)
library(readr)
library(biomaRt)
### UVM
# GSE202687

setwd("~/Library/CloudStorage/OneDrive-MoffittCancerCenter/Moffitt_DaisyS/Data/13_UVM/GSE202687")

myGSE<-"GSE202687"
gset <- getGEO(myGSE, GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL24676", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
cli_GSE202687<- pData(gset)
write.csv(cli_GSE202687, "GSE202687_pheno.csv")

########################################################
# laod expression 
expr<-read.csv("merged_raw_counts.csv")
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
colnames(expr)[1] <- "ensembl_gene_id"
gene_info <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), filters = 'ensembl_gene_id', values = expr$ensembl_gene_id, mart = ensembl)
data_matrix_with_symbols <- merge(expr, gene_info, by = "ensembl_gene_id", all.x = TRUE)
tail(data_matrix_with_symbols)
dim(data_matrix_with_symbols)
#[1] 63677    17
data_matrix_with_symbols[, 1] <- data_matrix_with_symbols[, 17]
data_matrix<-data_matrix_with_symbols[, -17]
write.csv(data_matrix,"GSE202687_TPM.csv")


## loaded the cleaned clinical data processed at background
cli<-read.csv("GSE202687_pheno_clean.csv")
expr<-read.csv("GSE202687_TPM.csv")
expr_t<-t(expr)
write.table(expr_t, "expr_t_tmp.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)
expr <- read_tsv("expr_t_tmp.tsv")
expr[1:3,1:3]
cli[1:3,1:3]
colnames(expr)[1] <- "sampleID"
## orgnize to data matrix
merged_data <- merge(cli,expr,by = "sampleID")
merged_data[1:3,1:3]
dim(merged_data)
# [1]   15 63695
write.csv(merged_data,"UVM_GSE202687.csv")
