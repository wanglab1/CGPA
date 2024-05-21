
## Melanoma
library(readr)
setwd("~/Library/CloudStorage/OneDrive-MoffittCancerCenter/DaisyS/Data/XuefengW_lab/Melanoma/harmonized_melanoma") #file path

melanoma_expr <- read_tsv("melanoma_all_TPM.tsv")

melanoma_expr_t<-t(melanoma_expr)

setwd("~/OneDrive - Moffitt Cancer Center/DaisyS/1_Projects/Biwei_CGPA_CancerRes_DS/Organized.data/Melanoma") # save path

write.table(melanoma_expr_t, "melanoma_all_TPM_t.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)

# this conversation can delete the V1,V2, etc
melanoma_expr <- read_tsv("melanoma_all_TPM_t.tsv")
melanoma_cli<- read.csv("RNA-CancerCell-MORRISON1-metadata.csv")

# write.csv(melanoma_expr,"melanoma_all_TPM_t_test.csv")
# melanoma_expr<-read.csv("melanoma_all_TPM_t_test.csv")
# melanoma_expr_reduce<-melanoma_expr[,-1]

melanoma_expr[1:3,1:3]
dim(melanoma_expr)
# [1]   442 14051
melanoma_cli[1:3,1:3]
dim(melanoma_cli)
# [1] 442  14

save(melanoma_expr, melanoma_cli, file = "melanoma_meta.RData")
merged_data <- merge(melanoma_cli,melanoma_expr,by = "PatientID")
write.csv(merged_data,"Melanoma_meta.csv")


