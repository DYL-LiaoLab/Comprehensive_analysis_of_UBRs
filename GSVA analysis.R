# GSVA analysis.
# Long deyu
# 2022.09.08

library(msigdbr)
library(data.table)
library(GSVA)
library(limma)
library(dplyr)
library(reshape2)
library(ggplot2)
library(stringr)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 

all_tpm<-read.table(gzfile("TcgaTargetGtex_rsem_gene_tpm.gz"),
                    sep = "\t")
head(all_tpm[1:4,1:3])
colnames(all_tpm)<-all_tpm[1,]
all_tpm<-all_tpm[-1,]


head(all_tpm$sample)
head(gsub("\\..*","",all_tpm$sample))
rownames(all_tpm)<-gsub("\\..*","",all_tpm$sample)
all_tpm<-all_tpm[,-1]
all_tpm[1:3,1:4]


table(gsub("-.*","",colnames(all_tpm)))

length(which(grepl("TCGA",colnames(all_tpm))==TRUE | grepl("GTEX",colnames(all_tpm))==TRUE))
all_tpm<-all_tpm[,which(grepl("TCGA",colnames(all_tpm))==TRUE | grepl("GTEX",colnames(all_tpm))==TRUE)]

# all_tpm<-all_tpm[,which(grepl("GTEX",colnames(all_tpm))==FALSE)]

pwer_anno<-read.csv("human_ubiquitination_regulators.csv")
colnames(pwer_anno)

all_tpm$ensembl_id<-rownames(all_tpm)

pwer_tpm<-inner_join(pwer_anno,all_tpm,by="ensembl_id")
rownames(pwer_tpm)<-pwer_tpm$ensembl_id

pwer_tpm1<-pwer_tpm[,-c(1:5)]


msigdbr_show_species()

h <- msigdbr(species = "Homo sapiens", 
             category = "H") 
h <- select(h, gs_name, ensembl_gene) %>% #
  as.data.frame %>% 
  split(., .$gs_name) %>% 
  lapply(., function(x)(x$ensembl_gene)) #


gs <- lapply(h, unique)

head(gs)

save(gs, file = "hallmark.gs.RData")

(load("hallmark.gs.RData"))

pwer_tpm1<-as.data.frame(lapply(pwer_tpm1,as.numeric)) 
rownames(pwer_tpm1)<-pwer_tpm$ensembl_id
gsva_es<-gsva(as.matrix(pwer_tpm1),gs,min.sz>1)

write.csv(gsva_es,"gsva_logtpm+0.001.csv")


gene_pathway<-read.csv("gsva_logtpm+0.001.csv",row.names = 1)

gene_pathway1<-as.data.frame(t(gene_pathway))

gene_exp<-read.csv("pwer基因所有癌症表达-origin.csv",row.names = 1)
gene_exp1<-as.data.frame(t(gene_exp))

identical(rownames(gene_exp1),rownames(gene_pathway1))

pwer_all<-cbind(gene_exp1,gene_pathway1)
write.csv(pwer_all,file = "gene-pathway-corr-input.csv")
# cor(pwer_all,method = "pearson")
# pheatmap::pheatmap(cor(pwer_all,method = "pearson"))

library(corrplot)
library(Hmisc)
library(RColorBrewer)
library(circlize)
pwer_all_corr<-rcorr(as.matrix(pwer_all),type = "pearson")
dev.off()
corrplot(pwer_all_corr$r[1:20,881:901],
         p.mat = pwer_all_corr$P[1:20,881:901]
)

temp<-c()
for(i in 853:901){
  print(colnames(pwer_all_corr$r)[i])
  temp[i-852]<-length(which(abs(pwer_all_corr$r[1:852,i])>0.5 & pwer_all_corr$P[1:852,i]<0.01))
}

temp_h<-c()
for(i in 1:852){
  temp_h[i]<-length(which(abs(pwer_all_corr$r[i,853:901])>0.5 & pwer_all_corr$P[i,853:901]<0.01))
}

