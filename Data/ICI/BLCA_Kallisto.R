# Kallisto BLCA
library(readr)
setwd("~/OneDrive - Moffitt Cancer Center/Moffitt_DaisyS/Data/BLCA.Kidney/Study/ICI_treatment/FollowUp/Kallisto_BLCA_2017")
ex<-read.csv("data_kallisto_reformed.csv")
ex[1:3,1:3]
ex_t<-t(ex)
write.table(ex_t, "data_kallisto_reformed_t.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)
ex <- read_tsv("data_kallisto_reformed_t.tsv")
ex[1:3,1:3]
cli<-read.csv("data_clinical_cleaned.csv")
cli[1:3,1:3]
## orgnize to data matrix
merged_data <- merge(cli,ex,by = "patient_id")
merged_data[1:3,1:3]
dim(merged_data)
# [1]    25 35456
write.csv(merged_data,"Kallisto_BLCA.csv")
