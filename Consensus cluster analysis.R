# Consensus cluster analysis.
# Long deyu
# 2022.09.20


library(ConsensusClusterPlus)
library(ggplot2)
library(dplyr)
library(reshape2)
library(cola)

pwer_tpm<-read.csv("ubiquitination_regulators-TPM.csv",check.names = FALSE)
cancer_name<-read.csv("cancer_name.csv",check.names = FALSE)
phenotype<-read.table(gzfile("TcgaTargetGTEX_phenotype.txt"),
                      sep = "\t")
colnames(phenotype)<-phenotype[1,]
phenotype<-phenotype[-1,]
colnames(phenotype)<-c("sample","detailed_category",        
                       "primary disease or tissue","primary_site",            
                       "sample_type","gender","study")

table(phenotype$primary_site)
tcga_phen<-subset(phenotype,phenotype$study=="TCGA")
cancer_name1<-as.data.frame(t(cancer_name))
cancer_name1$detailed_category<-rownames(cancer_name1)
cancer_name1[12,2]<-"Pheochromocytoma & Paraganglioma"
cancer_name1[25,2]<-"Diffuse Large B-Cell Lymphoma"
tcga_phen1<-inner_join(tcga_phen,cancer_name1,by = "detailed_category")
colnames(tcga_phen1)[8]<-"name"
write.csv(tcga_phen1,"cancer-name.csv",row.names = FALSE)
rownames(pwer_tpm)<-pwer_tpm$gene_symbol

cancer_sample<-pwer_tpm[,c(1:5,which(as.numeric(substr(colnames(pwer_tpm),14,15))<10))]


cancer_sample_tpm<-2^as.data.frame(lapply(cancer_sample[,-c(1:5)],as.numeric))
colnames(cancer_sample_tpm)<-gsub("\\.","-",colnames(cancer_sample_tpm))
cancer_sample_tpm<-cbind(cancer_sample[,1:5],cancer_sample_tpm)

hub_gene<-read.csv("mcode_cluster_top5-网络图_边文件.csv")
hub_gene<-subset(hub_gene,hub_gene$cluster=="Cluster1")

hub_gene1<-unique(hub_gene$Source)

temp_hub<-cancer_sample_tpm[hub_gene1,]
temp_hub<-temp_hub[order(temp_hub$Type),] 
temp_hub<-temp_hub[,-(1:5)]

cancer_type<-unique(tcga_phen1$name)
cancer_type<-cancer_type[order(cancer_type)]


for(i in 1:length(cancer_type)){
  # i=17
  print(i)
  print(cancer_type[i])
  temp1<-tcga_phen1[which(tcga_phen1$name==cancer_type[i]),] # 
  patient_id<-temp1$sample
  temp2<-temp_hub[,intersect(patient_id,colnames(temp_hub))] 
  temp_hub1<-sweep(temp2,1,apply(temp2,1,median)) 
  temp_hub1<-as.matrix(temp_hub1)
  cluster_pam <- ConsensusClusterPlus(
    d = temp_hub1, 
    maxK = 10, 
    seed = 2022, 
    reps = 2000, 
    pItem = 0.8, 
    pFeature = 1, 
    clusterAlg = 'pam', 
    distance = 'spearman', 
    innerLinkage = 'complete', 
    finalLinkage = 'complete', 
    corUse = 'pairwise.complete.obs', 
    title = paste(cancer_type[i],"pam_spearman_2000",sep = "_"), plot = 'png', #聚类簇评估结果的输出格式
    verbose = TRUE,writeTable = TRUE
  )
  
}

