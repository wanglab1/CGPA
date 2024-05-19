# SKCM GSE115821
library(GEOquery)
library(readr)
library(data.table)
library(readxl)
library(dplyr)
library(IOBR)
library(biomaRt)
setwd("~/Library/CloudStorage/OneDrive-MoffittCancerCenter/Moffitt_DaisyS/Data/9_SKCM/GSE115821")
#download the data 
myGSE<-"GSE115821"
gset <- getGEO(myGSE, GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL11154", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
cli<- pData(gset)
write.csv(cli, "GSE115821_pheno.csv")
# convert expr file ID to gene symbol
cli<-read.csv("GSE115821_pheno_clean.csv")
expr<-read.csv("GSE115821_MGH_counts.csv.gz")
expr[1:5,1:5]
dim(expr)
gene_lengths <- expr$Length
counts <- expr[, 7:ncol(expr)]
# Calculate RPK (Reads Per Kilobase)
rpk <- sweep(counts, 1, gene_lengths / 1000, "/")

# Calculate the scaling factor
scaling_factors <- colSums(rpk)

# Calculate TPM
tpm <- sweep(rpk, 2, scaling_factors, "/") * 1e6

# Convert to data frame for better handling
tpm <- as.data.frame(tpm)
tpm[1:5,1:5]
dim(tpm)
# [1] 75253    37
# Save the TPM data to a CSV file if needed
write.csv(tpm, "GSE115821_MGH_TPM.csv")
## convert
expr<-read.csv("GSE115821_MGH_TPM.csv")
cli<-read.csv("GSE115821_pheno_clean.csv")
expr_t<-t(expr)
expr_t[1:5,1:5]
write.table(expr_t, "GSE115821_temp.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)
expr_tpm <- read_tsv("GSE115821_temp.tsv")
expr_tpm[1:3,1:3]
# merge
merged_data <- merge(cli,expr_tpm,by = "sampleID")
merged_data[1:3,1:3]
dim(merged_data)
# [1]    23 21191
write.csv(merged_data, "SKCM_GSE115821.csv")
