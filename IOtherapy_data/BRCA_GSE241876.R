## getting data from GEO
library(GEOquery)
library(readr)
### BRCA_GSE241876
# GSE241876

setwd("~/Library/CloudStorage/OneDrive-MoffittCancerCenter/Moffitt_DaisyS/1_Projects/Biwei_CGPA_CancerRes_DS/Organized.data/BRCA/GSE241876")
myGSE<-"GSE241876"

gset <- getGEO(myGSE, GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL20301", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
cli_GSE241876<- pData(gset)
write.csv(cli_GSE241876, "GSE241876_pheno_data.csv")
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE241876&format=file&file=GSE241876%5FDeseqNormalizedCount%2Ecsv%2Egz"
path <- paste(urld, "acc=GSE241876", "file=GSE241876_DeseqNormalizedCount.csv.gz", sep="&");
ex_GSE241876 <- as.matrix(data.table::fread(path), rownames=1)
ex_GSE241876 <- data.table::fread(path)
ex_GSE241876[, 1] <- ex_GSE241876[, 22]
ex_GSE241876 <- ex_GSE241876[, -((ncol(ex_GSE241876)-2):ncol(ex_GSE241876))]
colnames(ex_GSE241876)[1] <- "Title" # change the first column name
save(ex_GSE241876, cli_GSE241876, file = "GSE241876.RData")

## loaded the cleaned clinical data processed at background
cli_GSE241876<-read.csv("GSE241876_pheno_clean.csv")
ex_GSE241876_t<-t(ex_GSE241876)
write.table(ex_GSE241876_t, "ex_GSE241876_t.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)
ex_GSE241876 <- read_tsv("ex_GSE241876_t.tsv")
ex_GSE241876[1:3,1:3]
cli_GSE241876[1:3,1:3]
## orgnize to data matrix
merged_data <- merge(cli_GSE241876,ex_GSE241876,by = "Title")
merged_data[1:3,1:3]
dim(merged_data)
# [1]    19 50660
write.csv(merged_data,"BRCA_GSE241876.csv")
