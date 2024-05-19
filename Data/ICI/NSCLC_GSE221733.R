## getting data from GEO
library(GEOquery)
library(limma)
library(umap)
### NSCLC
# GSE221733

setwd("~/Library/CloudStorage/OneDrive-MoffittCancerCenter/Moffitt_DaisyS/1_Projects/Biwei_CGPA_CancerRes_DS/Organized.data/NSCLC")

myGSE<-"GSE221733"
gset <- getGEO(myGSE, GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL18573", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

cli_221733<- pData(gset)

write.csv(cli_221733, "GSE221733_pheno_data.csv")

#save(ex_248249, cli_248249, file = "GSE248249.RData")

setwd("~/Library/CloudStorage/OneDrive-MoffittCancerCenter/Moffitt_DaisyS/1_Projects/Biwei_CGPA_CancerRes_DS/Organized.data/NSCLC/GSE221733")
cli_221733<-read.csv("GSE221733_pheno_clean.csv")
ex_221733<-read.csv("GSE221733_4301_CTA_norm.csv")
merged_data <- merge(cli_221733,ex_221733,by = "Study")
write.csv(merged_data,"NSCLC_GEO221733.csv")
