# Identification of differentially expressed genes in pan-cancer.
# Long deyu
# 2022.09.04

library(readr)
library(readxl)
library(TCGAbiolinks)
library(DT)
library(SummarizedExperiment)
library(dplyr)
library(TCGAbiolinksGUI)

gdc_project<-getGDCprojects()
gdc_project_tcga<-subset(gdc_project,grepl("TCGA",gdc_project$project_id))

phenotype<-read.table(gzfile("TcgaTargetGTEX_phenotype.txt.gz"),
                      sep = "\t")
colnames(phenotype)<-phenotype[1,]
phenotype<-phenotype[-1,]
colnames(phenotype)<-c("sample","detailed_category",        
                       "primary disease or tissue","primary_site",            
                       "sample_type","gender","study")

table(phenotype$primary_site)

all_tpm<-read.table(gzfile("TcgaTargetGtex_rsem_gene_tpm.gz"),
                    sep = "\t")

all_tpm[1:3,1:4]
colnames(all_tpm)<-all_tpm[1,]
all_tpm<-all_tpm[-1,]
all_tpm[1:3,1:4]

head(all_tpm$sample)
head(gsub("\\..*","",all_tpm$sample))
rownames(all_tpm)<-gsub("\\..*","",all_tpm$sample)
all_tpm<-all_tpm[,-1]
all_tpm[1:3,1:4]

head(phenotype$sample)
head(gsub("-.*","",phenotype$sample))
table(gsub("-.*","",phenotype$sample))

which(gsub("-.*","",phenotype$sample)=="K" | gsub("-.*","",phenotype$sample)=="TARGET")
phenotype<-phenotype[-which(gsub("-.*","",phenotype$sample)=="K" | gsub("-.*","",phenotype$sample)=="TARGET"),]
table(gsub("-.*","",phenotype$sample))

phenotype_tcga<-phenotype[phenotype[,7]=="TCGA",]
# tcga_cancer_tpye<-unique(phenotype_tcga$detailed_category)
tcga_cancer_tpye<-phenotype_tcga[,c(2,4)]
tcga_cancer_tpye_redum<-distinct(tcga_cancer_tpye)
tcga_cancer_tpye_redum<-tcga_cancer_tpye_redum[-34,] 

phenotype_gtex<-phenotype[phenotype[,7]=="GTEX",]
gtex_tissue_type<-phenotype_gtex[,c(4)]
gtex_tissue_type_redum<-unique(gtex_tissue_type)
gtex_tissue_type_redum<-gtex_tissue_type_redum[-13]
phenotype_gtex<-phenotype_gtex[phenotype_gtex[,5]=="Normal Tissue",]

pair_cancer<-c() 
unpair_cancer<-c() # 
none_normal_cancer<-c() # 
for(m in 1:length(tcga_cancer_tpye_redum$detailed_category)){
  # m=5
  print(m)
  print(tcga_cancer_tpye_redum$detailed_category[m])
  temp_tissue<-tcga_cancer_tpye_redum$primary_site[m]
  jishu=0
  for(n in 1:length(gtex_tissue_type_redum)){
   if(grepl(toupper(temp_tissue),toupper(gtex_tissue_type_redum[n])) | grepl(toupper(gtex_tissue_type_redum[n]),toupper(temp_tissue))){
      gtex<-phenotype_gtex[phenotype_gtex[,4]==gtex_tissue_type_redum[n],]
      jishu=jishu+1
      pair_cancer<-c(pair_cancer,tcga_cancer_tpye_redum$detailed_category[m])
    }else{
      next
    }
  }
  print(jishu)
  
  if(jishu==0){
    unpair_cancer<-c(unpair_cancer,tcga_cancer_tpye_redum$detailed_category[m])
  }
  
  tcga<-phenotype_tcga[phenotype_tcga[,2]==tcga_cancer_tpye_redum$detailed_category[m],]
  tcga_tumor_sample<-which(as.numeric(substr(tcga$sample,14,15))<10)
  tcga_normal_sample<-which(as.numeric(substr(tcga$sample,14,15))>=10)
  tcga_tumor_sampleid<-tcga$sample[tcga_tumor_sample]
  tcga_normal_sampleid<-tcga$sample[tcga_normal_sample]

  if(jishu==0 & length(tcga_normal_sampleid)==0){
    none_normal_cancer<-c(none_normal_cancer,tcga_cancer_tpye_redum$detailed_category[m])
    sample_count<-c(length(tcga_tumor_sampleid),
                    length(tcga_normal_sampleid),
                    0)
    write.csv(sample_count,
              file = paste(tcga_cancer_tpye_redum$detailed_category[m],"sample_number.csv",sep = "_"),row.names = FALSE)
    next
    
  }
  
  if(jishu==1){
    gtex_normal<-gtex$sample
    all_normal<-c(tcga_normal_sampleid,gtex_normal)
    sample_count<-c(length(tcga_tumor_sampleid),
                    length(tcga_normal_sampleid),
                    length(gtex_normal))
    write.csv(sample_count,
              file = paste(tcga_cancer_tpye_redum$detailed_category[m],"sample_number.csv",sep = "_"),row.names = FALSE)

    normal_sample_tpm<-subset(all_tpm,select = all_normal)
    tumor_sample_tpm<-subset(all_tpm,select = tcga_tumor_sampleid)
    
    all_sample<-c(all_normal,tcga_tumor_sampleid)
    all_sample_tpm<-subset(all_tpm,select = all_sample)
    rownames(all_sample_tpm)<-rownames(all_tpm)
    write.csv(all_sample_tpm,
              file =paste(tcga_cancer_tpye_redum$detailed_category[m],"FPKM.csv",sep = "_"),
              quote = FALSE )
    
    gene_id<-rownames(all_tpm)
    normal_barcode<-colnames(normal_sample_tpm)
    tumor_barcode<-colnames(tumor_sample_tpm)
    
    normal_sample_tpm<-2^as.data.frame(lapply(normal_sample_tpm,as.numeric))
    tumor_sample_tpm<-2^as.data.frame(lapply(tumor_sample_tpm,as.numeric))
    # normal_sample_tpm[1:3,1:4]
    colnames(normal_sample_tpm)<-normal_barcode
    colnames(tumor_sample_tpm)<-tumor_barcode
    rownames(normal_sample_tpm)<-gene_id
    rownames(tumor_sample_tpm)<-gene_id
    write.csv(normal_sample_tpm,
              file = paste(tcga_cancer_tpye_redum$detailed_category[m],"normal_FPKM+1.csv",sep = "_"),
              quote = FALSE)
    write.csv(tumor_sample_tpm,
              file = paste(tcga_cancer_tpye_redum$detailed_category[m],"tumor_FPKM+1.csv",sep = "_"),
              quote = FALSE)
    
    identical(rownames(normal_sample_tpm),rownames(tumor_sample_tpm))
    
    log2FC<-c()
    fold_change<-c()
    p.value<-c()
    normal_mean<-c()
    tumor_mean<-c()
    for(x in 1:length(rownames(all_tpm))){
      print(x)
      print(rownames(all_tpm)[x])
      temp_normal<-as.vector(t(normal_sample_tpm[x,]))
      temp_tumor<-as.vector(t(tumor_sample_tpm[x,]))
      temp_normal<-as.numeric(temp_normal)
      temp_tumor<-as.numeric(temp_tumor)
      normal_mean[x]<-mean(temp_normal)
      tumor_mean[x]<-mean(temp_tumor)
      fold_change[x]<-round(tumor_mean[x]/normal_mean[x],4)
      log2FC<-log2(fold_change)
      buffer<-wilcox.test(temp_normal,temp_tumor,conf.level = 0.95)
      p.value[x]<-buffer$p.value
    }
    
    cancer_result<-data.frame(Ensembl_id=rownames(all_tpm),
                              normal_mean=normal_mean,
                              tumor_mean=tumor_mean,
                              Fold_Change=fold_change,
                              log2FC=log2FC,
                              p.value=p.value)
    cancer_result1<-cancer_result[order(cancer_result$p.value),]
    cancer_result1$p.adjust<-p.adjust(cancer_result1$p.value,method = "BH")
    write.csv(cancer_result1,
              file = paste(tcga_cancer_tpye_redum$detailed_category[m],"wilcox_FPKM+0.001.csv",sep = "_"),
              row.names = FALSE,quote = FALSE)
    
  }else{
    gtex_normal<-c()
    all_normal<-c(tcga_normal_sampleid,gtex_normal)
    sample_count<-c(length(tcga_tumor_sampleid),
                    length(tcga_normal_sampleid),
                    length(gtex_normal))
    write.csv(sample_count,
              file = paste(tcga_cancer_tpye_redum$detailed_category[m],"sample_number.csv",sep = "_"),row.names = FALSE)
    normal_sample_tpm<-subset(all_tpm,select = all_normal)
    tumor_sample_tpm<-subset(all_tpm,select = tcga_tumor_sampleid)
    all_sample<-c(all_normal,tcga_tumor_sampleid)
    all_sample_tpm<-subset(all_tpm,select = all_sample)
    rownames(all_sample_tpm)<-rownames(all_tpm)
    write.csv(all_sample_tpm,
              file =paste(tcga_cancer_tpye_redum$detailed_category[m],"FPKM.csv",sep = "_"),
              quote = FALSE )
    gene_id<-rownames(all_tpm)
    normal_barcode<-colnames(normal_sample_tpm)
    tumor_barcode<-colnames(tumor_sample_tpm)
    
    normal_sample_tpm<-2^as.data.frame(lapply(normal_sample_tpm,as.numeric))
    tumor_sample_tpm<-2^as.data.frame(lapply(tumor_sample_tpm,as.numeric))
    colnames(normal_sample_tpm)<-normal_barcode
    colnames(tumor_sample_tpm)<-tumor_barcode
    rownames(normal_sample_tpm)<-gene_id
    rownames(tumor_sample_tpm)<-gene_id

    write.csv(normal_sample_tpm,
              file = paste(tcga_cancer_tpye_redum$detailed_category[m],"normal_FPKM+1.csv",sep = "_"),
              quote = FALSE)
    write.csv(tumor_sample_tpm,
              file = paste(tcga_cancer_tpye_redum$detailed_category[m],"tumor_FPKM+1.csv",sep = "_"),
              quote = FALSE)
    
    identical(rownames(normal_sample_tpm),rownames(tumor_sample_tpm))
    
    log2FC<-c()
    fold_change<-c()
    p.value<-c()
    normal_mean<-c()
    tumor_mean<-c()
    for(x in 1:length(rownames(all_tpm))){
      print(x)
      print(rownames(all_tpm)[x])
      temp_normal<-as.vector(t(normal_sample_tpm[x,]))
      temp_tumor<-as.vector(t(tumor_sample_tpm[x,]))
      temp_normal<-as.numeric(temp_normal)
      temp_tumor<-as.numeric(temp_tumor)
      normal_mean[x]<-mean(temp_normal)
      tumor_mean[x]<-mean(temp_tumor)
      fold_change[x]<-round(tumor_mean[x]/normal_mean[x],4)
      log2FC<-log2(fold_change)
      buffer<-wilcox.test(temp_normal,temp_tumor,conf.level = 0.95)
      p.value[x]<-buffer$p.value
    }
    
    cancer_result<-data.frame(Ensembl_id=rownames(all_tpm),
                              normal_mean=normal_mean,
                              tumor_mean=tumor_mean,
                              Fold_Change=fold_change,
                              log2FC=log2FC,
                              p.value=p.value)
    cancer_result1<-cancer_result[order(cancer_result$p.value),]
    cancer_result1$p.adjust<-p.adjust(cancer_result1$p.value,method = "BH")
    write.csv(cancer_result1,
              file = paste(tcga_cancer_tpye_redum$detailed_category[m],"wilcox_FPKM+0.001.csv",sep = "_"),
              row.names = FALSE,quote = FALSE)
    
    # ##############################################################################
  }
}


all_sample_number<-read.csv(paste(tcga_cancer_tpye_redum$detailed_category[1],"sample_number.csv",sep = "_"))
colnames(all_sample_number)<-tcga_cancer_tpye_redum$detailed_category[1]
rownames(all_sample_number)<-c("TCGA-tumor","TCGA-normal","GTEX")
for(i in 2:length(tcga_cancer_tpye_redum$detailed_category)){
  temp<-read.csv(paste(tcga_cancer_tpye_redum$detailed_category[i],"sample_number.csv",sep = "_"))
  all_sample_number<-cbind(all_sample_number,temp)
  colnames(all_sample_number)[i]<-tcga_cancer_tpye_redum$detailed_category[i]
}

tcga_cancer_tpye_rmnone<-tcga_cancer_tpye_redum[-c(1,15,25,33),]
all_cancer_foldchange<-read.csv(paste(tcga_cancer_tpye_rmnone$detailed_category[1],"wilcox_FPKM+0.001.csv",sep = "_"))
all_cancer_padj<-read.csv(paste(tcga_cancer_tpye_rmnone$detailed_category[1],"wilcox_FPKM+0.001.csv",sep = "_"))
# log2FC
all_cancer_foldchange<-all_cancer_foldchange[,c(1,5)]
all_cancer_padj<-all_cancer_padj[,c(1,7)]
colnames(all_cancer_foldchange)[2]<-tcga_cancer_tpye_rmnone$detailed_category[1]
colnames(all_cancer_padj)[2]<-tcga_cancer_tpye_rmnone$detailed_category[1]
for(i in 2:length(tcga_cancer_tpye_rmnone$detailed_category)){
  print(i)
  temp<-read.csv(paste(tcga_cancer_tpye_rmnone$detailed_category[i],"wilcox_FPKM+0.001.csv",sep = "_"))
  temp1<-temp[,c(1,5)]
  temp2<-temp[,c(1,7)]
  all_cancer_foldchange<-full_join(all_cancer_foldchange,temp1,by="Ensembl_id")
  all_cancer_padj<-full_join(all_cancer_padj,temp2,by="Ensembl_id")
  colnames(all_cancer_foldchange)[i+1]<-tcga_cancer_tpye_rmnone$detailed_category[i]
  colnames(all_cancer_padj)[i+1]<-tcga_cancer_tpye_rmnone$detailed_category[i]
}

write.csv(all_cancer_foldchange,
          file = "all_cancer_log2FC.csv",quote = FALSE,row.names=FALSE)
write.csv(all_cancer_padj,
          file = "all_cancer_padj.csv",quote = FALSE,row.names=FALSE)


tcga_cancer_tpye_rmnone<-tcga_cancer_tpye_redum[-c(1,15,25,33),]
all_cancer_diffgene_FC<-read.csv(paste(tcga_cancer_tpye_rmnone$detailed_category[1],"wilcox_FPKM+0.001.csv",sep = "_"))

all_cancer_diffgene_FC<-subset(all_cancer_diffgene_FC,
                               abs(all_cancer_diffgene_FC$log2FC)>=1 & all_cancer_diffgene_FC$p.adjust <0.01)
write.csv(all_cancer_diffgene_FC,
          file = paste(tcga_cancer_tpye_rmnone$detailed_category[1],"wilcox_diff_tpm.csv",sep = "_"),
          quote = FALSE,row.names = FALSE)

all_cancer_diffgene_FC<-all_cancer_diffgene_FC[,c(1,5)]
colnames(all_cancer_diffgene_FC)[2]<-tcga_cancer_tpye_rmnone$detailed_category[1]
for(i in 2:length(tcga_cancer_tpye_rmnone$detailed_category)){
  print(i)
  temp<-read.csv(paste(tcga_cancer_tpye_rmnone$detailed_category[i],"wilcox_FPKM+0.001.csv",sep = "_"))

  temp<-subset(temp,
               abs(temp$log2FC)>=1 & temp$p.adjust<0.01)
  write.csv(temp,
            file = paste(tcga_cancer_tpye_rmnone$detailed_category[i],"wilcox_diff_tpm.csv",sep = "_"),
            quote = FALSE,row.names = FALSE)
  temp1<-temp[,c(1,5)]
  all_cancer_diffgene_FC<-full_join(all_cancer_diffgene_FC,temp1,by="Ensembl_id")
  colnames(all_cancer_diffgene_FC)[i+1]<-tcga_cancer_tpye_rmnone$detailed_category[i]
}

write.csv(tcga_cancer_tpye_redum,file = "tcga_cancer_tpye_redum.csv",
          row.names = FALSE,quote = FALSE)

write.csv(all_cancer_diffgene_FC,
          file = "all_cancer_DEGs.csv",
          row.names = FALSE,quote = FALSE)