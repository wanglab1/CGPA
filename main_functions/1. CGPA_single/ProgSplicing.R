#-------------------------------------------------
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(DBI)
library(RSQLite)
#-------------------------------------------------------------------------------
# Process data, process splseq, combine all spliceseq the cancer type together
#-------------------------------------------------------------------------------
file_folder = "/home/4467777/CGPA_2024/oncosplicing/splseq/"


sub_func = function(input_file){
  dt = fread(paste0(file_folder,input_file))
  dt=dt%>%select(Splice_Event,Gene_Symbol,Chr_Strand,Splice_Type,Upstream_Exon,Downstream_Exon,Alt_Exons,AltExons_IsoName,SpliceIn_IsoName,
                  SpliceOut_IsoName,Splice_Novel,HRmedOS,pvalHRmedOS,HRfitOS,pvalHRfitOS,
                  HRmedPFI,pvalHRmedPFI,HRfitPFI,pvalHRfitPFI
  )
  dt$cancer=paste0(sapply(strsplit(input_file, "_//s*"), tail, 1))
  dt$cancer=gsub(".csv","",dt$cancer)
  return(dt)
}

ACC= sub_func("SpliceSeq_info_ACC.csv")
BLCA= sub_func("SpliceSeq_info_BLCA.csv")
BRCA= sub_func("SpliceSeq_info_BRCA.csv")
CESC= sub_func("SpliceSeq_info_CESC.csv")
CHOL= sub_func("SpliceSeq_info_CHOL.csv")
COAD= sub_func("SpliceSeq_info_COAD.csv")
DLBC= sub_func("SpliceSeq_info_DLBC.csv")
ESCA= sub_func("SpliceSeq_info_ESCA.csv")
GBM= sub_func("SpliceSeq_info_GBM.csv")
HNSC= sub_func("SpliceSeq_info_HNSC.csv")
KICH= sub_func("SpliceSeq_info_KICH.csv")
KIRC= sub_func("SpliceSeq_info_KIRC.csv")
KIRP= sub_func("SpliceSeq_info_KIRP.csv")


LGG= sub_func("SpliceSeq_info_LGG.csv")
LIHC= sub_func("SpliceSeq_info_LIHC.csv")
LUAD= sub_func("SpliceSeq_info_LUAD.csv")
LUSC= sub_func("SpliceSeq_info_LUSC.csv")
MESO= sub_func("SpliceSeq_info_MESO.csv")
OV= sub_func("SpliceSeq_info_OV.csv")
PAAD= sub_func("SpliceSeq_info_PAAD.csv")


SARC= sub_func("SpliceSeq_info_SARC.csv")
SKCM= sub_func("SpliceSeq_info_SKCM.csv")
STAD= sub_func("SpliceSeq_info_STAD.csv")
UCS= sub_func("SpliceSeq_info_UCS.csv")
UVM= sub_func("SpliceSeq_info_UVM.csv")
READ= sub_func("SpliceSeq_info_READ.csv")

LAML = fread(paste0(file_folder,"SpliceSeq_info_LAML.csv"))
LAML=LAML%>%select(Splice_Event,Gene_Symbol,Chr_Strand,Splice_Type,Upstream_Exon,Downstream_Exon,Alt_Exons,AltExons_IsoName,SpliceIn_IsoName,
               SpliceOut_IsoName,Splice_Novel,HRmedOS,pvalHRmedOS,HRfitOS,pvalHRfitOS)
LAML$HRmedPFI = NA
LAML$pvalHRmedPFI = NA
LAML$HRfitPFI = NA
LAML$pvalHRfitPFI = NA
LAML$cancer = "LAML"
LAML=LAML%>%select(Splice_Event,Gene_Symbol,Chr_Strand,Splice_Type,Upstream_Exon,Downstream_Exon,Alt_Exons,AltExons_IsoName,SpliceIn_IsoName,
               SpliceOut_IsoName,Splice_Novel,HRmedOS,pvalHRmedOS,HRfitOS,pvalHRfitOS,
               HRmedPFI,pvalHRmedPFI,HRfitPFI,pvalHRfitPFI,cancer
)

names(LAML) = c("Splice_Event","Gene_Symbol","Chr_Strand","Splice_Type","Upstream_Exon","Downstream_Exon","Alt_Exons","AltExons_IsoName","SpliceIn_IsoName",
                "SpliceOut_IsoName","Splice_Novel","HRmedOS","pvalHRmedOS","HRfitOS","pvalHRfitOS",
                "HRmedPFI","pvalHRmedPFI","HRfitPFI","pvalHRfitPFI","cancer")

PCPG = fread(paste0(file_folder,"SpliceSeq_info_PCPG.csv"))
PCPG=PCPG%>%select(Splice_Event,Gene_Symbol,Chr_Strand,Splice_Type,Upstream_Exon,Downstream_Exon,Alt_Exons,AltExons_IsoName,SpliceIn_IsoName,
                   SpliceOut_IsoName,Splice_Novel,HRmedPFI,pvalHRmedPFI,HRfitPFI,pvalHRfitPFI)
PCPG$HRmedOS = NA
PCPG$pvalHRmedOS = NA
PCPG$HRfitOS = NA
PCPG$pvalHRfitOS = NA
PCPG$cancer = "PCPG"
PCPG= PCPG %>%select(Splice_Event,Gene_Symbol,Chr_Strand,Splice_Type,Upstream_Exon,Downstream_Exon,Alt_Exons,AltExons_IsoName,SpliceIn_IsoName,
                     SpliceOut_IsoName,Splice_Novel,HRmedOS,pvalHRmedOS,HRfitOS,pvalHRfitOS,
                     HRmedPFI,pvalHRmedPFI,HRfitPFI,pvalHRfitPFI,cancer)

names(PCPG) = c("Splice_Event","Gene_Symbol","Chr_Strand","Splice_Type","Upstream_Exon","Downstream_Exon","Alt_Exons","AltExons_IsoName","SpliceIn_IsoName",
                "SpliceOut_IsoName","Splice_Novel","HRmedOS","pvalHRmedOS","HRfitOS","pvalHRfitOS",
                "HRmedPFI","pvalHRmedPFI","HRfitPFI","pvalHRfitPFI","cancer")


PRAD = fread(paste0(file_folder,"SpliceSeq_info_PRAD.csv"))
PRAD=PRAD%>%select(Splice_Event,Gene_Symbol,Chr_Strand,Splice_Type,Upstream_Exon,Downstream_Exon,Alt_Exons,AltExons_IsoName,SpliceIn_IsoName,
                   SpliceOut_IsoName,Splice_Novel,HRmedPFI,pvalHRmedPFI,HRfitPFI,pvalHRfitPFI)
PRAD$HRmedOS = NA
PRAD$pvalHRmedOS = NA
PRAD$HRfitOS = NA
PRAD$pvalHRfitOS = NA
PRAD$cancer = "PRAD"
PRAD= PRAD %>%select(Splice_Event,Gene_Symbol,Chr_Strand,Splice_Type,Upstream_Exon,Downstream_Exon,Alt_Exons,AltExons_IsoName,SpliceIn_IsoName,
                     SpliceOut_IsoName,Splice_Novel,HRmedOS,pvalHRmedOS,HRfitOS,pvalHRfitOS,
                     HRmedPFI,pvalHRmedPFI,HRfitPFI,pvalHRfitPFI,cancer)

names(PRAD) = c("Splice_Event","Gene_Symbol","Chr_Strand","Splice_Type","Upstream_Exon","Downstream_Exon","Alt_Exons","AltExons_IsoName","SpliceIn_IsoName",
                "SpliceOut_IsoName","Splice_Novel","HRmedOS","pvalHRmedOS","HRfitOS","pvalHRfitOS",
                "HRmedPFI","pvalHRmedPFI","HRfitPFI","pvalHRfitPFI","cancer")


TGCT = fread(paste0(file_folder,"SpliceSeq_info_TGCT.csv"))
TGCT=TGCT%>%select(Splice_Event,Gene_Symbol,Chr_Strand,Splice_Type,Upstream_Exon,Downstream_Exon,Alt_Exons,AltExons_IsoName,SpliceIn_IsoName,
                   SpliceOut_IsoName,Splice_Novel,HRmedPFI,pvalHRmedPFI,HRfitPFI,pvalHRfitPFI)
TGCT$HRmedOS = NA
TGCT$pvalHRmedOS = NA
TGCT$HRfitOS = NA
TGCT$pvalHRfitOS = NA
TGCT$cancer = "TGCT"
TGCT= TGCT %>%select(Splice_Event,Gene_Symbol,Chr_Strand,Splice_Type,Upstream_Exon,Downstream_Exon,Alt_Exons,AltExons_IsoName,SpliceIn_IsoName,
                     SpliceOut_IsoName,Splice_Novel,HRmedOS,pvalHRmedOS,HRfitOS,pvalHRfitOS,
                     HRmedPFI,pvalHRmedPFI,HRfitPFI,pvalHRfitPFI,cancer)

names(TGCT) = c("Splice_Event","Gene_Symbol","Chr_Strand","Splice_Type","Upstream_Exon","Downstream_Exon","Alt_Exons","AltExons_IsoName","SpliceIn_IsoName",
                "SpliceOut_IsoName","Splice_Novel","HRmedOS","pvalHRmedOS","HRfitOS","pvalHRfitOS",
                "HRmedPFI","pvalHRmedPFI","HRfitPFI","pvalHRfitPFI","cancer")

THCA = fread(paste0(file_folder,"SpliceSeq_info_THCA.csv"))
THCA=THCA%>%select(Splice_Event,Gene_Symbol,Chr_Strand,Splice_Type,Upstream_Exon,Downstream_Exon,Alt_Exons,AltExons_IsoName,SpliceIn_IsoName,
                   SpliceOut_IsoName,Splice_Novel,HRmedPFI,pvalHRmedPFI,HRfitPFI,pvalHRfitPFI)
THCA$HRmedOS = NA
THCA$pvalHRmedOS = NA
THCA$HRfitOS = NA
THCA$pvalHRfitOS = NA
THCA$cancer = "THCA"
THCA= THCA %>%select(Splice_Event,Gene_Symbol,Chr_Strand,Splice_Type,Upstream_Exon,Downstream_Exon,Alt_Exons,AltExons_IsoName,SpliceIn_IsoName,
                     SpliceOut_IsoName,Splice_Novel,HRmedOS,pvalHRmedOS,HRfitOS,pvalHRfitOS,
                     HRmedPFI,pvalHRmedPFI,HRfitPFI,pvalHRfitPFI,cancer)

names(THCA) = c("Splice_Event","Gene_Symbol","Chr_Strand","Splice_Type","Upstream_Exon","Downstream_Exon","Alt_Exons","AltExons_IsoName","SpliceIn_IsoName",
                "SpliceOut_IsoName","Splice_Novel","HRmedOS","pvalHRmedOS","HRfitOS","pvalHRfitOS",
                "HRmedPFI","pvalHRmedPFI","HRfitPFI","pvalHRfitPFI","cancer")

THYM = fread(paste0(file_folder,"SpliceSeq_info_THYM.csv"))
THYM=THYM%>%select(Splice_Event,Gene_Symbol,Chr_Strand,Splice_Type,Upstream_Exon,Downstream_Exon,Alt_Exons,AltExons_IsoName,SpliceIn_IsoName,
                   SpliceOut_IsoName,Splice_Novel,HRmedPFI,pvalHRmedPFI,HRfitPFI,pvalHRfitPFI)
THYM$HRmedOS = NA
THYM$pvalHRmedOS = NA
THYM$HRfitOS = NA
THYM$pvalHRfitOS = NA
THYM$cancer = "THYM"
THYM= THYM %>%select(Splice_Event,Gene_Symbol,Chr_Strand,Splice_Type,Upstream_Exon,Downstream_Exon,Alt_Exons,AltExons_IsoName,SpliceIn_IsoName,
                     SpliceOut_IsoName,Splice_Novel,HRmedOS,pvalHRmedOS,HRfitOS,pvalHRfitOS,
                     HRmedPFI,pvalHRmedPFI,HRfitPFI,pvalHRfitPFI,cancer)
names(THYM) = c("Splice_Event","Gene_Symbol","Chr_Strand","Splice_Type","Upstream_Exon","Downstream_Exon","Alt_Exons","AltExons_IsoName","SpliceIn_IsoName",
                "SpliceOut_IsoName","Splice_Novel","HRmedOS","pvalHRmedOS","HRfitOS","pvalHRfitOS",
                "HRmedPFI","pvalHRmedPFI","HRfitPFI","pvalHRfitPFI","cancer")


UCEC = fread(paste0(file_folder,"spliceseq_info_UCEC.csv"))
UCEC=UCEC%>%select(Splice_Event,Gene_Symbol,Chr_Strand,Splice_Type,Upstream_Exon,Downstream_Exon,Alt_Exons,AltExons_IsoName,SpliceIn_IsoName,
                   SpliceOut_IsoName,Splice_Novel,HRmedOS,pvalHRmedOS,HRfitOS,pvalHRfitOS,HRmedPFI,pvalHRmedPFI,HRfitPFI,pvalHRfitPFI)

UCEC$cancer = "UCEC"
UCEC= UCEC %>%select(Splice_Event,Gene_Symbol,Chr_Strand,Splice_Type,Upstream_Exon,Downstream_Exon,Alt_Exons,AltExons_IsoName,SpliceIn_IsoName,
                     SpliceOut_IsoName,Splice_Novel,HRmedOS,pvalHRmedOS,HRfitOS,pvalHRfitOS,
                     HRmedPFI,pvalHRmedPFI,HRfitPFI,pvalHRfitPFI,cancer)
names(UCEC) = c("Splice_Event","Gene_Symbol","Chr_Strand","Splice_Type","Upstream_Exon","Downstream_Exon","Alt_Exons","AltExons_IsoName","SpliceIn_IsoName",
                "SpliceOut_IsoName","Splice_Novel","HRmedOS","pvalHRmedOS","HRfitOS","pvalHRfitOS",
                "HRmedPFI","pvalHRmedPFI","HRfitPFI","pvalHRfitPFI","cancer")

all = data.frame(rbind(ACC,BLCA,BRCA,CESC,CHOL,COAD,DLBC,ESCA,GBM,HNSC,KICH,KIRC,KIRP,LAML,LGG,LIHC,LUAD,LUSC,MESO,OV,PAAD,PCPG,PRAD,READ,SARC,SKCM,STAD,TGCT,THCA,THYM,UCEC,UCS,UVM))
write_fst(all,"/home/4467777/CGPA_2024/oncosplicing/splseq/spliceseq_all.fst")
#write_parquet(all, "/home/4467777/CGPA_2024/oncosplicing/splseq/spliceseq_all.parquet") # faster

# save databases into sql databases
all = read_fst("/home/4467777/CGPA_2024/oncosplicing/splseq/spliceseq_all.fst")
annot = read_fst("/home/4467777/CGPA_2024/oncosplicing/splseq/annot_splseq.fst")

all = merge(all,annot,by.y="Splice_event",by.x="Splice_Event",all.x=TRUE)

setwd("/home/4467777/CGPA_2024/oncosplicing/splseq/")


gene_symbol = data.frame(all$Gene_Symbol)
names(gene_symbol)  = "gene_symbol"
full_data = all
con <- dbConnect(RSQLite::SQLite(), "spliceseq_full.db") # create an empty data base

dbWriteTable(con, "GeneSymbol", gene_symbol, overwrite = TRUE, row.names = FALSE)
dbWriteTable(con, "spliceseq", full_data, overwrite = TRUE, row.names = FALSE)

dbDisconnect(con) # this database contains two datasets, genesymbol and spliceseq
#------------------------------------------------------------------------------
# Check how to get the results
#------------------------------------------------------------------------------
# Read only the 'gene_symbol' column

con <- dbConnect(RSQLite::SQLite(), "M:/dept/Dept_BBSR/Projects/Wang_Xuefeng/CGPA/move_022924/CGPA_single/www/Splicing/spliceseq_full.db")
selected_genes <- c("FN1")  # replace with your actual genes of interest
result <- dbGetQuery(con, sprintf("SELECT rowid FROM GeneSymbol WHERE gene_symbol = ('%s')",selected_genes))

id_list_str <- paste0("(", paste(result$rowid, collapse = ","), ")")
query_str <- paste0("SELECT * FROM spliceseq WHERE rowid IN ", id_list_str)
dt =  dbGetQuery(con, query_str)

# test (too slow)
dt = read_fst("M:/dept/Dept_BBSR/Projects/Wang_Xuefeng/CGPA/move_022924/CGPA_single/Splicing/spliceseq_all.fst")
# test end

#-------------------------------------------------------------------------------
# Generate the final figure #
#-------------------------------------------------------------------------------

dt_gene=dt%>%select(Splice_Event,Splice_Type,HRmedOS,pvalHRmedOS,cancer) 

names(dt_gene)[3]="HR_OS"
names(dt_gene)[4]="pval_OS"
OS_fit = dt_gene
OS_fit=OS_fit%>%filter(!(is.na(pval_OS)))
OS_fit$HR=ifelse(OS_fit$HR_OS>1,"Poor prognosis","Favorable prognosis")
OS_fit$HR=replace(OS_fit$HR, OS_fit$pval_OS>0.05,"NS (p > 0.05)")
OS_fit$HR=factor(OS_fit$HR,levels=c("Favorable prognosis","Poor prognosis","NS (p > 0.05)"))

OS_fit=OS_fit%>%mutate(size_label=case_when(pval_OS<0.05 ~ "-log10(0.05)",
                                            pval_OS<0.01 ~ "-log10(0.01)",
                                            pval_OS>=0.05 ~ "NS"))



p=OS_fit %>%
  ggplot(aes(Splice_Event, cancer, fill = HR)) +
  geom_point(shape = 20,aes(color=HR,size=-log10(pval_OS))) +
  scale_color_manual(values=c("Poor prognosis"="blue", "Favorable prognosis"="red","NS (p > 0.05)"="grey")) + 
  scale_size_continuous(name="-log10(OS P-value)",
                        range  = c(1, 15), 
                        breaks = c(-log10(0.5), -log10(0.1), -log10(0.05),-log10(0.01)),
                        labels=c("-log10(0.5)","-log10(0.1)","-log10(0.05)","-log10(0.01)")
  )+
  labs(x = NULL, y = NULL) +
  scale_y_discrete(position = "right")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 60, vjust = 0.8, hjust=0.1))+
  coord_flip()+theme(aspect.ratio=0.4) +
  ggtitle("OS splicing")

p


dt_table = dt %>% select(Splice_Event,Gene_Symbol,Chr_Strand,Splice_Type,Upstream_Exon,Downstream_Exon,Alt_Exons,AltExons_IsoName,SpliceIn_IsoName,SpliceOut_IsoName,Splice_Novel)
dt_table = distinct(dt_table, Splice_Event, .keep_all = TRUE)
dt_table = dt_table[order(dt_table$Splice_Event,decreasing = T),]
dt_table

#------------------------------------------------------------------------------
# Annotation file #
#------------------------------------------------------------------------------
gene_location<-function(gene_name){
  data<-paste0("http://47.98.127.64:8080/bedFiles?fileName=",gene_name,".bed")
  bed <- read.table(data  ,header = F,fill = T,stringsAsFactors = FALSE)
  bed=as.data.frame(bed)
  inds=which(bed$V1 == "track")
  bed_seq=bed[c((inds[2]+1):(inds[3]-1)),]
  bed_seq2=bed_seq[,c(1:4,6)]
  names(bed_seq2)=c("Chr","Start","End","Splice_event","Strand")
  bed_seq2$Event=gsub("^[^_]+_","",bed_seq2$Splice_event)
  bed_seq2$Event=gsub("_.*","",bed_seq2$Event)
  bed_seq2$Event<-replace(bed_seq2$Event, bed_seq2$Event=="AA", "Alternate Acceptor")
  bed_seq2$Event<-replace(bed_seq2$Event, bed_seq2$Event=="AT", "Alternate Terminator")
  bed_seq2$Event<-replace(bed_seq2$Event, bed_seq2$Event=="ES", "Exon Skip")
  bed_seq2$Event<-replace(bed_seq2$Event, bed_seq2$Event=="RI", "Retained Intron")
  bed_seq2$Event<-replace(bed_seq2$Event, bed_seq2$Event=="AP", "Alternate Promoter ")
  bed_seq2$Event<-replace(bed_seq2$Event, bed_seq2$Event=="AD", "Alternate Donor")
  bed_seq2$Event<-replace(bed_seq2$Event, bed_seq2$Event=="ME", "Mutually Exclusive Exons")
  bed_seq2 = data.frame(bed_seq2[,c("Splice_event","Event","Chr","Strand","Start","End")])
  
  return(bed_seq2)
}

genes = read_fst("C:/Users/4467777/Desktop/TCGA_shiny/cgpa_single031222/www/geneName_mrna.fst")
genes = genes$geneName

annot_tb = lapply(1:10,function(x){
  tryCatch({
    gene_location(genes[x])
    
  },error=function(e){
    NULL
  })

})
annot_tb_final= data.frame(do.call(rbind,annot_tb))

#===================================================================================
#--------------------------------------------------------------------------------
# Get spladder datasets
#-------------------------------------------------------------------------------
file_folder = "/home/4467777/CGPA_2024/oncosplicing/spladder/"

input_file = "SplAdder_info_ACC.csv.gz"

sub_func = function(input_file){
  dt = fread(paste0(file_folder,input_file))
  dt=dt%>%select(Splice_Event,Gene_Symbol,Chr_Strand,Splice_Type, Event_Region, Alt_Region,AltRegion_IsoName,SpliceIn_IsoName,SpliceOut_IsoName,
                 Splice_Novel,HRmedOS,pvalHRmedOS,HRfitOS,pvalHRfitOS,
                 HRmedPFI,pvalHRmedPFI,HRfitPFI,pvalHRfitPFI
  )
  dt$cancer=paste0(sapply(strsplit(input_file, "_"), tail, 1))
  dt$cancer=gsub(".csv.gz","",dt$cancer)
  return(dt)
}

files = list.files(file_folder)
res = sapply(1:length(files),function(x){
  tryCatch({
    sub_func(files[x])
  },error=function(e){
    NULL
  })
})

res = data.frame(do.call(rbind,res))

PCPG = fread(paste0(file_folder,"SplAdder_info_PCPG.csv.gz"))
PCPG=PCPG%>%select(Splice_Event,Gene_Symbol,Chr_Strand,Splice_Type, Event_Region, Alt_Region,AltRegion_IsoName,SpliceIn_IsoName,SpliceOut_IsoName,
                   Splice_Novel,HRmedPFI,pvalHRmedPFI,HRfitPFI,pvalHRfitPFI)
PCPG$HRmedOS = NA
PCPG$pvalHRmedOS = NA
PCPG$HRfitOS = NA
PCPG$pvalHRfitOS = NA
PCPG$cancer = "PCPG"
PCPG= PCPG %>%select(Splice_Event,Gene_Symbol,Chr_Strand,Splice_Type, Event_Region, Alt_Region,AltRegion_IsoName,SpliceIn_IsoName,SpliceOut_IsoName,
                     Splice_Novel,HRmedOS,pvalHRmedOS,HRfitOS,pvalHRfitOS,
                     HRmedPFI,pvalHRmedPFI,HRfitPFI,pvalHRfitPFI,cancer)

names(PCPG) = c("Splice_Event","Gene_Symbol","Chr_Strand","Splice_Type","Event_Region","Alt_Region","AltRegion_IsoName","SpliceIn_IsoName","SpliceOut_IsoName",
                "Splice_Novel","HRmedOS","pvalHRmedOS","HRfitOS","pvalHRfitOS",
                "HRmedPFI","pvalHRmedPFI","HRfitPFI","pvalHRfitPFI","cancer")


PRAD = fread(paste0(file_folder,"SplAdder_info_PRAD.csv.gz"))
PRAD=PRAD%>%select(Splice_Event,Gene_Symbol,Chr_Strand,Splice_Type, Event_Region, Alt_Region,AltRegion_IsoName,SpliceIn_IsoName,SpliceOut_IsoName,
                   Splice_Novel,HRmedPFI,pvalHRmedPFI,HRfitPFI,pvalHRfitPFI)
PRAD$HRmedOS = NA
PRAD$pvalHRmedOS = NA
PRAD$HRfitOS = NA
PRAD$pvalHRfitOS = NA
PRAD$cancer = "PRAD"
PRAD= PRAD %>%select(Splice_Event,Gene_Symbol,Chr_Strand,Splice_Type, Event_Region, Alt_Region,AltRegion_IsoName,SpliceIn_IsoName,SpliceOut_IsoName,
                     Splice_Novel,HRmedOS,pvalHRmedOS,HRfitOS,pvalHRfitOS,
                     HRmedPFI,pvalHRmedPFI,HRfitPFI,pvalHRfitPFI,cancer)

names(PRAD) = c("Splice_Event","Gene_Symbol","Chr_Strand","Splice_Type","Event_Region","Alt_Region","AltRegion_IsoName","SpliceIn_IsoName","SpliceOut_IsoName",
                "Splice_Novel","HRmedOS","pvalHRmedOS","HRfitOS","pvalHRfitOS",
                "HRmedPFI","pvalHRmedPFI","HRfitPFI","pvalHRfitPFI","cancer")


TGCT = fread(paste0(file_folder,"SplAdder_info_TGCT.csv.gz"))
TGCT=TGCT%>%select(Splice_Event,Gene_Symbol,Chr_Strand,Splice_Type, Event_Region, Alt_Region,AltRegion_IsoName,SpliceIn_IsoName,SpliceOut_IsoName,
                   Splice_Novel,HRmedPFI,pvalHRmedPFI,HRfitPFI,pvalHRfitPFI)
TGCT$HRmedOS = NA
TGCT$pvalHRmedOS = NA
TGCT$HRfitOS = NA
TGCT$pvalHRfitOS = NA
TGCT$cancer = "TGCT"
TGCT= TGCT %>%select(Splice_Event,Gene_Symbol,Chr_Strand,Splice_Type, Event_Region, Alt_Region,AltRegion_IsoName,SpliceIn_IsoName,SpliceOut_IsoName,
                     Splice_Novel,HRmedOS,pvalHRmedOS,HRfitOS,pvalHRfitOS,
                     HRmedPFI,pvalHRmedPFI,HRfitPFI,pvalHRfitPFI,cancer)

names(TGCT) = c("Splice_Event","Gene_Symbol","Chr_Strand","Splice_Type","Event_Region","Alt_Region","AltRegion_IsoName","SpliceIn_IsoName","SpliceOut_IsoName",
                "Splice_Novel","HRmedOS","pvalHRmedOS","HRfitOS","pvalHRfitOS",
                "HRmedPFI","pvalHRmedPFI","HRfitPFI","pvalHRfitPFI","cancer")

THCA = fread(paste0(file_folder,"SplAdder_info_THCA.csv.gz"))
THCA=THCA%>%select(Splice_Event,Gene_Symbol,Chr_Strand,Splice_Type, Event_Region, Alt_Region,AltRegion_IsoName,SpliceIn_IsoName,SpliceOut_IsoName,
                   Splice_Novel,HRmedPFI,pvalHRmedPFI,HRfitPFI,pvalHRfitPFI)
THCA$HRmedOS = NA
THCA$pvalHRmedOS = NA
THCA$HRfitOS = NA
THCA$pvalHRfitOS = NA
THCA$cancer = "THCA"
THCA= THCA %>%select(Splice_Event,Gene_Symbol,Chr_Strand,Splice_Type, Event_Region, Alt_Region,AltRegion_IsoName,SpliceIn_IsoName,SpliceOut_IsoName,
                     Splice_Novel,HRmedOS,pvalHRmedOS,HRfitOS,pvalHRfitOS,
                     HRmedPFI,pvalHRmedPFI,HRfitPFI,pvalHRfitPFI,cancer)

names(THCA) = c("Splice_Event","Gene_Symbol","Chr_Strand","Splice_Type","Event_Region","Alt_Region","AltRegion_IsoName","SpliceIn_IsoName","SpliceOut_IsoName",
                "Splice_Novel","HRmedOS","pvalHRmedOS","HRfitOS","pvalHRfitOS",
                "HRmedPFI","pvalHRmedPFI","HRfitPFI","pvalHRfitPFI","cancer")

THYM = fread(paste0(file_folder,"SplAdder_info_THYM.csv.gz"))
THYM=THYM%>%select(Splice_Event,Gene_Symbol,Chr_Strand,Splice_Type, Event_Region, Alt_Region,AltRegion_IsoName,SpliceIn_IsoName,SpliceOut_IsoName,
                   Splice_Novel,HRmedPFI,pvalHRmedPFI,HRfitPFI,pvalHRfitPFI)
THYM$HRmedOS = NA
THYM$pvalHRmedOS = NA
THYM$HRfitOS = NA
THYM$pvalHRfitOS = NA
THYM$cancer = "THYM"
THYM= THYM %>%select(Splice_Event,Gene_Symbol,Chr_Strand,Splice_Type, Event_Region, Alt_Region,AltRegion_IsoName,SpliceIn_IsoName,SpliceOut_IsoName,
                     Splice_Novel,HRmedOS,pvalHRmedOS,HRfitOS,pvalHRfitOS,
                     HRmedPFI,pvalHRmedPFI,HRfitPFI,pvalHRfitPFI,cancer)
names(THYM) = c("Splice_Event","Gene_Symbol","Chr_Strand","Splice_Type","Event_Region","Alt_Region","AltRegion_IsoName","SpliceIn_IsoName","SpliceOut_IsoName",
               "Splice_Novel","HRmedOS","pvalHRmedOS","HRfitOS","pvalHRfitOS",
               "HRmedPFI","pvalHRmedPFI","HRfitPFI","pvalHRfitPFI","cancer")

res_final = data.frame(rbind(res,PCPG,PRAD,TGCT,THCA,THYM))

write_fst(res_final,"/home/4467777/CGPA_2024/oncosplicing/spladder/spladder_all.fst")

#
setwd("/home/4467777/CGPA_2024/oncosplicing/spladder/")
res_final = read_fst("/home/4467777/CGPA_2024/oncosplicing/spladder/spladder_all.fst")
res_final$Splice_Type = ifelse(res_final$Splice_Type=="A3","Alternative 3' site",ifelse(res_final$Splice_Type=="A5","Alternative 5' site",
                                                                                  ifelse(res_final$Splice_Type=="ES","Exon skip",
                                                                                         ifelse(res_final$Splice_Type=="ME","Mutually exclusive exons",
                                                                                               "Intron retention" ))))
gene_symbol = data.frame(res_final$Gene_Symbol)
names(gene_symbol)  = "gene_symbol"
full_data = res_final
con <- dbConnect(RSQLite::SQLite(), "spladder_all.db") # create an empty data base

dbWriteTable(con, "GeneSymbol", gene_symbol, overwrite = TRUE, row.names = FALSE)
dbWriteTable(con, "spladder", full_data, overwrite = TRUE, row.names = FALSE)

dbDisconnect(con) 


# test
con <- dbConnect(RSQLite::SQLite(), "M:/dept/Dept_BBSR/Projects/Wang_Xuefeng/CGPA/move_022924/CGPA_single/www/Splicing/spladder_all.db")
selected_genes <- c("BTNL9")  # replace with your actual genes of interest
result <- dbGetQuery(con, sprintf("SELECT rowid FROM GeneSymbol WHERE gene_symbol = ('%s')",selected_genes))

id_list_str <- paste0("(", paste(result$rowid, collapse = ","), ")")
query_str <- paste0("SELECT * FROM spladder WHERE rowid IN ", id_list_str)
dt =  dbGetQuery(con, query_str)
head(dt)
names(dt)
dbDisconnect(con) 
