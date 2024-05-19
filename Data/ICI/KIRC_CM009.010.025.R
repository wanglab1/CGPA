# KIRC CM009.010.025
library(readr)
library(data.table)
library(readxl)
setwd("/Users/daisys/Library/CloudStorage/OneDrive-MoffittCancerCenter/Moffitt_DaisyS/Data/10_KIRC/KIRC_CM010")
expr <- data.table(read_excel("KIRC_CM009.010.025.expr.cli.xlsx", sheet = "expr"))
cli<-data.table(read_excel("KIRC_CM009.010.025.expr.cli.xlsx", sheet = "cli"))
dim(expr) 
# [1] 43893   312
expr[1:3,1:3]
cli[1:3,1:3]
# colnames(expr)[1] <- "RNA_ID"
## convert
expr<-t(expr)
write.table(expr, "expr_temp.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)
expr <- read_tsv("expr_temp.tsv")

colnames(expr)[1] <- "RNA_ID"
expr[1:3,1:3]
# merge
merged_data <- merge(cli,expr,by = "RNA_ID")
merged_data[1:3,1:3]
dim(merged_data)
# 311 43913
write.csv(merged_data,"KIRC_CM009.010.025.csv")
