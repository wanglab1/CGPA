## getting data from GEO
library(GEOquery)
library(limma)
library(umap)
### HNSCC
setwd("/Users/daisys/Library/CloudStorage/OneDrive-MoffittCancerCenter/DaisyS/Data/XuefengW_lab/HNSCC")

myGSE<-"GSE159067"
gset <- getGEO(myGSE, GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL18573	", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

#warning: not every GEO data can get the expression by this funciton. 
# ex_159067 <- exprs(gset)
ex_159067 <- read.table("GSE159067_IHN_log2cpm_data.txt.gz", header = TRUE, sep = "\t", comment.char = "!")
cli_159067<- pData(gset)

# write.csv(ex_159067, "GSE12345_expression.csv")
# write.csv(cli_159067, "GSE12345_pheno_data.csv")

save(ex_159067, cli_159067, file = "GSE159067.RData")
####

ex_159067<-read.csv("GSE12345_expression.csv")
dim(ex_159067)
merged_data <- merge(cli_159067,ex_159067,by = "geo_accession")
write.csv(merged_data,"HNSCC_GEO159067.csv")
