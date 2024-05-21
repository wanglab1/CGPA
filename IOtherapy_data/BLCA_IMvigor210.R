# BLCA
## IMvigor210

### BLCA
library(readr)

load("~/Library/CloudStorage/OneDrive-MoffittCancerCenter/DaisyS/Data/XuefengW_lab/BLCA/imvigor/IMvigor210.RData")

setwd("~/Library/CloudStorage/OneDrive-MoffittCancerCenter/DaisyS/1_Projects/Biwei_CGPA_CancerRes_DS/Organized.data/BLCA") # save path

write.table(eset, "IMvigor210_expression.tsv", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
write.csv(pdata, "IMvigor210_cli.csv")

BLCA_expr <- read_tsv("IMvigor210_expression.tsv")
BLCA_cli<-read.csv("IMvigor210_cli.csv")
head(BLCA_cli)

BLCA_expr[1:3,1:3]
BLCA_expr_t<-t(BLCA_expr)
BLCA_expr_t[1:3,1:3]
## conversion the format
write.table(BLCA_expr_t, "BLCA_TPM_t.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)
BLCA_expr_2 <- read_tsv("BLCA_TPM_t.tsv")
BLCA_expr_2[1:3,1:3]
dim(BLCA_expr_2)
# [1]   348 26131

merged_data <- merge(BLCA_cli,BLCA_expr_2,by = "ID")

write.csv(merged_data,"BLCA_IMvigor210.csv")
