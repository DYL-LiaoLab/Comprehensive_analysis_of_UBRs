# Analysis of somatic mutation in pan-cancer.
# Long deyu
# 2022.09.02

library(maftools)
library(ggplot2)
library(DESeq2)
library(tidyverse)
library(dplyr)
library(readr)
library(readxl)
library(ComplexHeatmap)
library(pheatmap)
library(circlize)
library(tidyr)
library(reshape2)
library(TCGAbiolinks)
library(DT)
library(SummarizedExperiment)
library(RColorBrewer)

gdc_project<-getGDCprojects()
gdc_project_tcga<-subset(gdc_project,grepl("TCGA",gdc_project$project_id))

pwer<-read.csv("human_ubiquitination_regulators.csv")
gene_list<-unique(pwer$gene_symbol)

vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
vc_cols=c("Thistle1","RoyalBlue2","Cyan1","Brown2","HotPink",
          "Coral1","Gold1","Purple")
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)



for(x in 1:length(gdc_project_tcga$project_id)){
  # x=18
  print(x)
  print(gdc_project_tcga$project_id[x])
  clinical<-read.csv(paste(gdc_project_tcga$project_id[x],"clinical_info.csv",sep = "."))
  synapse_maf<-read.maf(maf = "mc3.v0.2.8.PUBLIC.maf.gz",
                        clinicalData = clinical)
  # 
  synapse_maf_pwer<-subsetMaf(maf = synapse_maf,
                              tsb = clinical$Tumor_Sample_Barcode)
  
  plotmafSummary(maf = synapse_maf_pwer,
                 rmOutlier = TRUE,
                 addStat = "median",
                 dashboard = TRUE,
                 titvRaw = FALSE,
                 color = vc_cols
  )
  oncostrip(maf = synapse_maf_pwer,top=20)
  pdf(paste(gdc_project_tcga$project_id[x],"titv.pdf",sep = "."),width = 12,height = 8)
  titv(maf = synapse_maf_pwer,plot = TRUE,useSyn = TRUE)
  dev.off()
  png(paste(gdc_project_tcga$project_id[x],"titv.png",sep = "."),width = 1200,height = 800)
  titv(maf = synapse_maf_pwer,plot = TRUE,useSyn = TRUE)
  dev.off()
  
  pdf(paste(gdc_project_tcga$project_id[x],"mafsummary.pdf",sep = "."),width = 12,height = 8)
  plotmafSummary(maf = synapse_maf_pwer,
                 rmOutlier = TRUE,
                 addStat = "median",
                 dashboard = TRUE,
                 titvRaw = FALSE,
                 color = vc_cols
  )
  dev.off()
  png(paste(gdc_project_tcga$project_id[x],"mafsummary.png",sep = "."),width = 1200,height = 800)
  plotmafSummary(maf = synapse_maf_pwer,
                 rmOutlier = TRUE,
                 addStat = "median",
                 dashboard = TRUE,
                 titvRaw = FALSE,
                 color = vc_cols
  )
  dev.off()
  # 
  pdf(paste(gdc_project_tcga$project_id[x],".pdf",sep = "."),width = 12,height = 12)
  oncoplot(maf = synapse_maf_pwer,
           top = 20,
           colors = vc_cols,
           # fontSize = 1,
           legendFontSize = 1,
           annotationFontSize = 1,
           # clinicalFeatures = c(
           #                      "neoplasm_histologic_grade"
           #                      # "Variant_Classification"
           #                      )
           # removeNonMutated = TRUE
  )
  dev.off()
  
  png(paste(gdc_project_tcga$project_id[x],".png",sep = "."),width = 1200,height = 800)
  oncoplot(maf = synapse_maf_pwer,
           top = 20,
           colors = vc_cols,
           # fontSize = 1,
           legendFontSize = 1,
           annotationFontSize = 1,
           # clinicalFeatures = c(
           #                      "neoplasm_histologic_grade"
           #                      # "Variant_Classification"
           #                      )
           # removeNonMutated = TRUE
  )
  dev.off()
  gene.summary<-synapse_maf_pwer@gene.summary
  write.csv(gene.summary,
            file = paste(gdc_project_tcga$project_id[x],"gene_summary.csv",sep = "."),
            row.names = FALSE)

  sample_count<-length(synapse_maf_pwer@variants.per.sample$Tumor_Sample_Barcode)
  
  synapse_maf_pwer@summary
  
  synapse_maf_pwer@variants.per.sample

  sample_var<-synapse_maf_pwer@variant.type.summary
  write.csv(sample_var,
            file = paste(gdc_project_tcga$project_id[x],"sample_var.csv",sep = "."),
            row.names = FALSE)

  sample_summary<-synapse_maf_pwer@variant.classification.summary
  getSampleSummary(synapse_maf_pwer)
  write.csv(sample_summary,
            file = paste(gdc_project_tcga$project_id[x],"sample_summary.csv",sep = "."),
            row.names = FALSE)
  
  mutation_frequancy<-data.frame(Hugo_Symbol=gene.summary$Hugo_Symbol,
                                 cancer_type=round(gene.summary$AlteredSamples/sample_count,3)) # 保留3位小数
  write.csv(mutation_frequancy,
            file = paste(gdc_project_tcga$project_id[x],"mutation_frequancy.csv",sep = "."),
            row.names = FALSE)
  
  
  print(x)
  print("done")
  
  rm(synapse_maf)
  rm(synapse_maf_pwer)
  gc()
}









