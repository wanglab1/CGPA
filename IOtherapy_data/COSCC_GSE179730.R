## getting data from GEO
library(GEOquery)
library(readr)
### GSE179730
# oral-cavity squamous cell carcinoma/COSCC
setwd("/Users/daisys/Library/CloudStorage/OneDrive-MoffittCancerCenter/Moffitt_DaisyS/Data/2_OralCavitySquamousCellCarcinoma.COSCC/GSE179730")

myGSE<-"GSE179730"
gset <- getGEO(myGSE, GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL24676	", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
cli_GSE179730<- pData(gset)
write.csv(cli_GSE179730, "GSE179730_pheno_data.csv")
### getting matrix 
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE179730&format=file&file=GSE179730%5FRNAseq%2DcombinedCPM%2Etxt%2Egz"
path <- paste(urld, "acc=GSE179730", "file=GSE179730_RNAseq-combinedCPM.txt.gz", sep="&");
ex_GSE179730 <- data.table::fread(path)
colnames(ex_GSE179730)[1] <- "Title" # change the first column name
save(ex_GSE179730, cli_GSE179730, file = "GSE179730.RData")

## loaded the cleaned clinical data processed at background
cli_GSE179730<-read.csv("GSE179730_pheno_clean.csv")
ex_GSE179730_t<-t(ex_GSE179730)
write.table(ex_GSE179730_t, "ex_GSE179730_t.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)
ex_GSE179730 <- read_tsv("ex_GSE179730_t.tsv")
ex_GSE179730[1:3,1:3]
cli_GSE179730[1:3,1:3]
## orgnize to data matrix
merged_data <- merge(cli_GSE179730,ex_GSE179730,by = "Title")
merged_data[1:3,1:3]
dim(merged_data)
# [1] 23 24867
write.csv(merged_data,"COSCC_GSE179730.csv")

