# Calculate the correlation coefficient between tisues in three different source datasets.
# Long deyu
# 2022.05.22

library(circlize)
library(ggplot2)
library(dplyr)
library(openxlsx)
library(writexl)
library(readxl)
library(Hmisc)
library(dplyr)
library(ggraph)
library(igraph)
library(tidyverse)
library(dendextend)
library(colormap)
library(kableExtra)
library(readr)
library(data.table)
library(reshape2)


ub_wer<-read.csv("human_ubiquitination_regulators.csv")

gtex<-fread("rna_tissue_gtex.tsv")
colnames(gtex)
colnames(ub_wer)
colnames(gtex)[2]<-"gene_symbol"
gtex<-gtex[,c(2,3,6)]
ub_gtex<-left_join(ub_wer,gtex,by="gene_symbol") 
colnames(ub_gtex)
ub_gtex<-na.omit(ub_gtex)
ub_gtex_kuan<-reshape2::dcast(ub_gtex,gene_symbol+Family+Type+ensembl_id+Species~Tissue,
                              value.var = "nTPM",
                              fun.aggregate = mean
)
# 
FANTOM<-fread("rna_tissue_fantom.tsv")
colnames(FANTOM)
colnames(FANTOM)<-c("Gene","Gene name","Tissue","TPM","pTPM","nTPM") # 参照链接：https://www.proteinatlas.org/about/assays+annotation#hpa_rna
colnames(ub_wer)
colnames(FANTOM)[2]<-"gene_symbol"
FANTOM<-FANTOM[,c(2,3,6)]

ub_FANTOM<-left_join(ub_wer,FANTOM,by="gene_symbol") 
colnames(ub_FANTOM)
ub_FANTOM<-na.omit(ub_FANTOM)
ub_FANTOM_kuan<-reshape2::dcast(ub_FANTOM,gene_symbol+Family+Type+ensembl_id+Species~Tissue,
                                value.var = "nTPM",
                                fun.aggregate = mean
)
# 
HPA<-fread("rna_tissue_hpa.tsv")
colnames(HPA)
colnames(ub_wer)
colnames(HPA)[2]<-"gene_symbol"
HPA<-HPA[,c(2,3,6)]

gtex_tissue<-c('skeletal muscle','tongue','bone marrow','testis','thymus',
               'tonsil','lymph node','parathyroid gland','heart muscle','liver',
               'ovary','thyroid gland','placenta','spleen','appendix',
               'cervix','lung','small intestine','stomach','duodenum',
               'rectum','colon','breast','endometrium','urinary bladder',
               'esophagus','smooth muscle','kidney','adipose tissue','pancreas',
               'skin','salivary gland','prostate','fallopian tube','gallbladder',
               'epididymis','seminal vesicle','adrenal gland')
gtex_tissue<-gtex_tissue[order(gtex_tissue)]

HPA_raw<-HPA 
for(i in 1:length(gtex_tissue)){
  # i=17
  print(gtex_tissue[i])
  if(i==1){
    total_hpa<-subset(HPA,HPA$Tissue==gtex_tissue[i])
  }else{
    temp<-subset(HPA,HPA$Tissue==gtex_tissue[i])
    total_hpa<-rbind(total_hpa,temp)
  }
}
HPA<-total_hpa # 
ub_HPA<-left_join(ub_wer,HPA,by="gene_symbol")
ub_HPA<-na.omit(ub_HPA)
colnames(ub_HPA)
ub_HPA_kuan<-reshape2::dcast(ub_HPA,gene_symbol+Family+Type+ensembl_id+Species~Tissue,
                             value.var = "nTPM",
                             fun.aggregate = mean
)


ub_gtex_kuan1<-ub_gtex_kuan[!duplicated(ub_gtex_kuan$gene_symbol),]
rownames(ub_gtex_kuan1)<-ub_gtex_kuan1$gene_symbol
ub_FANTOM_kuan1<-ub_FANTOM_kuan[!duplicated(ub_FANTOM_kuan$gene_symbol),]
rownames(ub_FANTOM_kuan1)<-ub_FANTOM_kuan1$gene_symbol
ub_HPA_kuan1<-ub_HPA_kuan[!duplicated(ub_HPA_kuan$gene_symbol),]
rownames(ub_HPA_kuan1)<-ub_HPA_kuan1$gene_symbol


gene_common<-intersect(ub_gtex_kuan1$gene_symbol,ub_FANTOM_kuan1$gene_symbol)
gene_common<-intersect(gene_common,ub_HPA_kuan1$gene_symbol)

ub_gtex_kuan2<-ub_gtex_kuan1[gene_common,]
ub_FANTOM_kuan2<-ub_FANTOM_kuan1[gene_common,]
ub_HPA_kuan2<-ub_HPA_kuan1[gene_common,]

colnames(ub_gtex_kuan2)<-paste(colnames(ub_gtex_kuan2),"gtex",sep = "_")
colnames(ub_gtex_kuan2)[1]<-"gene_symbol"
colnames(ub_FANTOM_kuan2)<-paste(colnames(ub_FANTOM_kuan2),"fantom",sep = "_")
colnames(ub_FANTOM_kuan2)[1]<-"gene_symbol"
colnames(ub_HPA_kuan2)<-paste(colnames(ub_HPA_kuan2),"hpa",sep = "_")
colnames(ub_HPA_kuan2)[1]<-"gene_symbol"

ub_FANTOM_kuan2<-ub_FANTOM_kuan2[,-c(2:5)]
ub_HPA_kuan2<-ub_HPA_kuan2[,-c(2:5)]

all_data<-inner_join(ub_gtex_kuan2,ub_FANTOM_kuan2,by="gene_symbol")
all_data<-inner_join(all_data,ub_HPA_kuan2,by="gene_symbol")
colnames(all_data)


all_data_corr<-rcorr(as.matrix(all_data[,6:length(all_data)]),type = "pearson")

all_data_corr_r<-all_data_corr$r
all_data_corr_p<-all_data_corr$P


all_r_chang<-as.data.frame(all_data_corr_r)
all_r_chang$tissue<-rownames(all_r_chang)
all_r_chang<-reshape2::melt(all_r_chang,id.vars="tissue")

all_p_chang<-as.data.frame(all_data_corr_p)
all_p_chang$tissue<-rownames(all_p_chang)
all_p_chang<-reshape2::melt(all_p_chang,id.vars="tissue")

identical(all_p_chang$tissue,all_r_chang$tissue)
identical(all_p_chang$variable,all_r_chang$variable)

all_corr<-cbind(all_r_chang,all_p_chang$value)
colnames(all_corr)<-c("tissue","tissue_target","r_value","p_value")


all_r_kuan<-all_data_corr_r[which(grepl("_hpa",rownames(all_data_corr_r))==FALSE),
                            which(grepl("_gtex",colnames(all_data_corr_r))==FALSE)]
all_p_kuan<-all_data_corr_p[which(grepl("_hpa",rownames(all_data_corr_p))==FALSE),
                            which(grepl("_gtex",colnames(all_data_corr_p))==FALSE)]

all_r_kuan[which(grepl("_fantom",rownames(all_r_kuan))==TRUE),which(grepl("_fantom",colnames(all_r_kuan))==TRUE)]<-0
all_p_kuan[which(grepl("_fantom",rownames(all_p_kuan))==TRUE),which(grepl("_fantom",colnames(all_p_kuan))==TRUE)]<-1


all_r_kuan1<-as.data.frame(all_r_kuan)
all_r_kuan1$tissue<-rownames(all_r_kuan1)
all_r_kuan2<-reshape2::melt(all_r_kuan1,id.vars="tissue")
all_p_kuan1<-as.data.frame(all_p_kuan)
all_p_kuan1$tissue<-rownames(all_p_kuan1)
all_p_kuan2<-reshape2::melt(all_p_kuan1,id.vars="tissue")

identical(all_r_kuan2$tissue,all_p_kuan2$tissue)
identical(all_r_kuan2$variable,all_p_kuan2$variable)

hebing_corr<-cbind(all_r_kuan2,all_p_kuan2$value)
colnames(hebing_corr)<-c("tissue","tissue_target","r_value","p_value")

hebing_corr1<-subset(hebing_corr,abs(hebing_corr$r_value)>0.5)
hebing_corr1<-subset(hebing_corr1,hebing_corr1$p_value<0.05)

hebing_corr2<-hebing_corr1
hebing_corr2$source_at<-gsub(".*_","",hebing_corr2$tissue)
hebing_corr2$target_at<-gsub(".*_","",hebing_corr2$tissue_target)
hebing_corr2$node1<-gsub("_.*","",hebing_corr2$tissue)
hebing_corr2$node2<-gsub("_.*","",hebing_corr2$tissue_target)

write.csv(hebing_corr2,"Correlation network diagram.csv",row.names = FALSE,quote = FALSE)
