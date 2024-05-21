## getting data from GEO
library(GEOquery)
library(readr)
library(biomaRt)
### NSCLC
# GSE135222

setwd("/Users/daisys/Library/CloudStorage/OneDrive-MoffittCancerCenter/Moffitt_DaisyS/Data/3_NSCLC/GSE135222")

myGSE<-"GSE135222"
gset <- getGEO(myGSE, GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL16791", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
cli_GSE135222<- pData(gset)
write.csv(cli_GSE135222, "GSE135222_pheno_data.csv")
########################################################
# laod expression 
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE135222&format=file&file=GSE135222%5FGEO%5FRNA%2Dseq%5Fomicslab%5Fexp%2Etsv%2Egz"
path <- paste(urld, "acc=GSE135222", "file=GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz", sep="&");
ex_GSE135222 <- data.table::fread(path)
colnames(ex_GSE135222)[1] <- "title" # change the first column name
ex_GSE135222[1:3,1:3]
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
colnames(ex_GSE135222)[1] <- "ensembl_gene_id"
gene_info <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                     +                    filters = 'ensembl_gene_id',
                     +                    values = ex_GSE135222$ensembl_gene_id,
                     +                    mart = ensembl)

data_matrix_with_symbols <- merge(ex_GSE135222, gene_info, by = "ensembl_gene_id", all.x = TRUE)
data_matrix_with_symbols[1:3,1:3]
# ensembl_gene_id NSCLC378 NSCLC1104
1 ENSG00000000003    16.86     34.47
2 ENSG00000000005     0.00      0.00
3 ENSG00000000419    79.04     62.69
write.csv(data_matrix_with_symbols,"GSE135222_exp.csv")
save(ex_GSE135222, cli_GSE135222, file = "GSE135222.RData")

## loaded the cleaned clinical data processed at background
cli_GSE135222<-read.csv("GSE135222_pheno_clean.csv")
ex_GSE135222<-read.csv("GSE135222_exp.csv")
ex_GSE135222_t<-t(ex_GSE135222)
write.table(ex_GSE135222_t, "ex_GSE135222_t.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)
ex_GSE135222 <- read_tsv("ex_GSE135222_t.tsv")
ex_GSE135222[1:3,1:3]
cli_GSE135222[1:3,1:3]
## orgnize to data matrix
merged_data <- merge(cli_GSE135222,ex_GSE135222,by = "title")
merged_data[1:3,1:3]
dim(merged_data)
# [1]    27 39674
write.csv(merged_data,"NSCLC_GSE135222.csv")
