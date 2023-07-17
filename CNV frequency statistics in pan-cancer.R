# Statistics of CNV gain and loss frequency for each cancer type.
# Long deyu
# 2022.08.30


library(TCGAbiolinks)
library(DT)
library(SummarizedExperiment)
library(ggplot2)
library(dplyr)
library(circlize)
library(reshape2)
library(tidyr)
library(tidyverse)
library(pheatmap)
library(ComplexHeatmap)

gdc_project<-getGDCprojects()

gdc_project_tcga<-subset(gdc_project,grepl("TCGA",gdc_project$project_id))


for(y in 1:length(gdc_project_tcga$id)){
  # y=1
  print(y)
  print(gdc_project_tcga$id[y])
  cancer_gistic2<-read.csv(file = paste0("GDC_Hub_CNV_gistic/",gdc_project_tcga$project_id[y],".gistic.tsv"),
                           sep = "\t",check.names = FALSE)

  cancer_gistic2$`Gene Symbol`<-gsub("\\..*","",cancer_gistic2$`Gene Symbol`)
  gain<-c()
  loss<-c()
  for(m in 1:length(cancer_gistic2$`Gene Symbol`)){
    # m=1
    # 
    gain_number<-length(which(cancer_gistic2[m,]==1))
    loss_number<-length(which(cancer_gistic2[m,]==-1))
    # 
    gain[m]<-round(gain_number/(length(cancer_gistic2)-1),3)
    loss[m]<-round(loss_number/(length(cancer_gistic2)-1),3)
  }

  temp<-data.frame(Gene_Symbol=cancer_gistic2$`Gene Symbol`,
                   gain=gain,
                   loss=loss)
  colnames(temp)<-c("Gene_Symbol",
                    paste(gdc_project_tcga$project_id[y],"gain",sep = "_"),
                    paste(gdc_project_tcga$project_id[y],"loss",sep = "_"))
  write.csv(temp,
            file = paste("GDC_Hub_CNV_gistic/",gdc_project_tcga$project_id[y],"_cnv_frequency.csv",sep = ""),
            row.names = FALSE,quote = FALSE)
  print("done！！")
  if(y==1){
    cnv_frequency<-temp
  }else{
    cnv_frequency<-full_join(cnv_frequency,temp,by="Gene_Symbol")
  }
  print(y)
  print("done")
}


write.csv(cnv_frequency,
          file = paste("GDC_Hub_CNV_gistic/","cancer_cnv_frequency.csv",sep = ""),
          row.names = FALSE,quote = FALSE)
