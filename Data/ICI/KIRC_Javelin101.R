library(readr)
setwd("~/OneDrive - Moffitt Cancer Center/Moffitt_DaisyS/Data/BLCA.Kidney/Study/Robert.J.Motzer_Renal_2020")
ex<-read.csv("Expr_matrix.csv")
ex[1:3,1:3]
ex_t<-t(ex)
write.table(ex_t, "ex_t.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)
ex <- read_tsv("ex_t.tsv")
ex[1:3,1:3]
cli<-read.csv("Cli_matrix.csv")
cli[1:3,1:3]
## orgnize to data matrix
merged_data <- merge(cli,ex,by = "ID")
merged_data[1:3,1:3]
dim(merged_data)
# [1]   726 22963
write.csv(merged_data,"KIRC_Javelin101.csv")
