# Tissue-specific analysis of ubiquitination regulators in mice.
# Long deyu
# 2022.05.17

library(geneRal) 
library(dplyr)
library(openxlsx)
library(writexl)
library(readxl)
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(MASS)
library(tidyr)
library(reshape2)
library(plyr)
library(ggrepel)
library(grid)
library(gridExtra)


samp<-read.csv("PRJNA375882_SraRunTable.txt")
samp<-samp[,c(27,31)]
fpkm<-read.csv("table6-41598_2017_4520_MOESM2_ESM.csv",skip = 1)
#
samp$name<-gsub("_._6w_\\d","",samp$Sample.Name)

tissue<-unique(gsub("_._.*","",colnames(fpkm)))
fpkm_mean<-fpkm[,1:3]
for(i in 4:length(tissue)){
  # i=14
  print(tissue[i])
  temp<-which(samp$name==tissue[i])
  buff<-samp[temp,]
  print(buff[1,2])
  # 提取对应的列
  temp_fpkm<-fpkm[,which(grepl(tissue[i],colnames(fpkm))==TRUE)]
  temp_mean<-apply(temp_fpkm,1,mean)
  fpkm_mean<-cbind(fpkm_mean,temp_mean)
  colnames(fpkm_mean)[length(fpkm_mean)]<-buff[1,2]
}


pwer<-read.csv("mouse_ubiquitination_regulators.csv")
# colnames(pwer)[2]<-"Type"
colnames(fpkm_mean)[1]<-"ensembl_id"


fpkm_pwer_mean<-inner_join(pwer,fpkm_mean,by="ensembl_id")
iEKPD2_tpm2<-subset(fpkm_pwer_mean,select = c(3,5,1,8:24))

table(pwer$Type)
gene_type<-table(pwer$Type)
gene_type<-as.data.frame(gene_type)
colnames(gene_type)<-c("Type","Count")

cols<-c(Writer="#A7A0FC",Reader="#F3E1A8", Multi="#C3FFA7",Eraser="#CF9FFB")

labs<-paste0(gene_type$Type,"(",
             round(gene_type$Count/sum(gene_type$Count)*100,2),"%)")


pie(gene_type$Count,labels = labs,col = rev(cols),
    border = "black",
    edges = 10000,
    init.angle = 90,
    clockwise = TRUE,
    lty = 1)



iEKPD2_tpm2<-na.omit(iEKPD2_tpm2)

iekpd_scatter<-melt(iEKPD2_tpm2,c("gene_symbol","Type","ensembl_id"))
colnames(iekpd_scatter)<-c("Gene_name","Type","Ensembl_Gene_ID","Tissue","TPM_value")


iekpd_scatter$TPM_value <-log(iekpd_scatter$TPM_value+1)


genemean<-aggregate(iekpd_scatter[,5],list(iekpd_scatter$Ensembl_Gene_ID),mean)
colnames(genemean)<-c("Ensembl_Gene_ID","genemean")

iekpd_scatter2<-join(iekpd_scatter,genemean,by="Ensembl_Gene_ID")

iekpd_scatter2$abs<-iekpd_scatter2$TPM_value - iekpd_scatter2$genemean

res_vec<-c()

for(tissue in unique(iekpd_scatter2$Tissue)){
  print(tissue)
  temp<-iekpd_scatter2[with(iekpd_scatter2,iekpd_scatter2$Tissue %in% tissue),]
  res<-rlm(temp$TPM_value ~0 +temp$genemean)
  res_vec=c(res$residuals,res_vec)
}
threshold<-2.5*sd(res_vec)


specific_genes<-vector() 
for (tissue in unique(iekpd_scatter2$Tissue)){ #for each tissue
  
  temp <- iekpd_scatter2[with(iekpd_scatter2, iekpd_scatter2$Tissue %in% tissue),] #extract the data for a specific tissue
  res<- rlm(temp$TPM_value ~0 + temp$genemean)#linear model for that tissue 
  temp$res<- res$residuals #add residual values to the matrix
  temp$diff<- abs(temp$res)-threshold #difference between gene's residual and threshold
  spec<-subset(temp, diff>0 & abs>0) #extract specific genes in each tissue
  specific_genes<- rbind(spec,specific_genes) #add these genes to the initial data
  pdf(file=paste(tissue,"specificityplot_logtpm.pdf",sep="_"),height=5,width=5)
  print(ggplot(temp, aes(x=genemean, y=TPM_value,label=Gene_name)) +
          geom_point(data=temp, col="black",size=0.5)+ #All data points will be black
          geom_point(data=subset(temp, diff>0 & abs>0),col="red",size=2)+ #Except the specific genes
          geom_text_repel(data=subset(temp, diff>0 & abs>0),segment.size  = 0.4,segment.color = "grey50",)+ #Add text to the specific genes
          geom_smooth(method=rlm, formula = y ~0 + x, size=0.5)+ #abline will be from rlm function that passes through 0,0
          xlab("mRNA mean abundance All Human Tissues")+
          ylab(paste("mRNA mean abundance",tissue,sep=" "))+
          theme(panel.background = element_blank(),
                panel.border=element_rect(fill=NA),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                strip.background=element_blank(),
                axis.text.x=element_text(colour="black"),
                axis.text.y=element_text(colour="black"),
                axis.ticks=element_line(colour="black"),
                plot.margin=unit(c(1,1,1,1),"line")))
  dev.off()
  png(filename = paste(tissue,"specificityplot_logtpm.png",sep="_"))
  print(ggplot(temp, aes(x=genemean, y=TPM_value,label=Gene_name)) +
          geom_point(data=temp, col="black",size=0.5)+ #All data points will be black
          geom_point(data=subset(temp, diff>0 & abs>0),col="red",size=2)+ #Except the specific genes
          geom_text_repel(data=subset(temp, diff>0 & abs>0),segment.size  = 0.4,segment.color = "grey50",)+ #Add text to the specific genes
          geom_smooth(method=rlm, formula = y ~0 + x, size=0.5)+ #abline will be from rlm function that passes through 0,0
          xlab("mRNA mean abundance All Mouse Tissues")+
          ylab(paste("mRNA mean abundance",tissue,sep=" "))+
          theme(panel.background = element_blank(),
                panel.border=element_rect(fill=NA),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                strip.background=element_blank(),
                axis.text.x=element_text(colour="black"),
                axis.text.y=element_text(colour="black"),
                axis.ticks=element_line(colour="black"),
                plot.margin=unit(c(1,1,1,1),"line")))
  dev.off()
}


write.table(specific_genes, file="mouse_tissue_specific_genes_FPKM.tsv",quote=FALSE, row.names=FALSE,sep="\t")
table(specific_genes$Tissue)


rownames(iEKPD2_tpm2)<-iEKPD2_tpm2$ensembl_id
iEKPD2_tpm2<-iEKPD2_tpm2[order(iEKPD2_tpm2$Type),]
temp1<-iEKPD2_tpm2[,-c(1:3)]

iEKPD2_tpm2_norm1<-t(apply(temp1,1,scale)) 
colnames(iEKPD2_tpm2_norm1)<-colnames(temp1)


annotation_row<-subset(iEKPD2_tpm2,select = 2)

ann_colors=list(Type=c(Writer="#A7A0FC",Reader="#F3E1A8", Multi="#C3FFA7",Eraser="#CF9FFB"))

iEKPD2_tpm2_norm1<-as.matrix(iEKPD2_tpm2_norm1)
library(circlize)

pheatmap(iEKPD2_tpm2_norm1,
         clustering_distance_rows = "pearson",
         # clustering_distance_rows ="spearman",
         # clustering_distance_rows ="euclidean",
         annotation_row = annotation_row,
         annotation_colors = ann_colors,
         cluster_rows = FALSE,
         show_rownames = FALSE,
         # col = colorRamp2(c(-3,-1.5,0,1.5,3), c("#2c7bb6","#abd9e9","floralwhite","#fdae61", "#d7191c"),space = "RGB"),
         # col = colorRamp2(c(-3.6,-1.8,0,1.8,3.6), c("#2c7bb6","#abd9e9","floralwhite","#fdae61", "#d7191c"),space = "RGB"),
         # col = colorRamp2(c(-4,-2,0,2,4), c("#2c7bb6","#abd9e9","floralwhite","#fdae61", "#d7191c"),space = "RGB"),
         col = colorRamp2(c(-4,-2,0,2,4), c("#2c7bb6","#abd9e9","White","#FF8C69", "#FF0000"),space = "RGB"),
         name = "z-scale FPKM",
         # split=annotation_row$Type
)


iekpd2_box<-cbind(iEKPD2_tpm2[,2],as.data.frame(iEKPD2_tpm2_norm1[,1:17]))
colnames(iekpd2_box)[1]<-"type"
iekpd2_box1<-melt(iekpd2_box,c("type"))
colnames(iekpd2_box1)<-c("type","Tissue","zscore_TPM")
colnames(temp1)
tissue<-sort(colnames(temp1))

library(ggpubr)
compare_means(zscore_TPM~Tissue,data = iekpd2_box1,
              ref.group = "testis",method = "t.test")

p<-ggplot(iekpd2_box1,aes(x=factor(Tissue,levels = tissue),y=zscore_TPM))+
  geom_flat_violin(aes(fill=Tissue),position = position_nudge(x=.05),
                   color="black",adjust=0.5,width = 1.1)+
  # geom_jitter(aes(color=Tissue),width = 0.12,size=0.01)+
  geom_boxplot(aes(fill=Tissue),width=0.4,
               # outlier.color = NA,
               outlier.size = 1.5,outlier.shape = 21,
               position = position_nudge(x=-.2))+
  stat_summary(fun= "mean",geom = "point",
               fill="white",shape=21,size=2.5,
               position = position_nudge(x=-.2))+
  stat_compare_means(method = "anova",label.y = -2.5)+
  stat_compare_means(label = "p.signif",method = "t.test",ref.group = "testis",
                     hide.ns = FALSE)+
  xlab("")+
  ylab("z-score PTM")+
  theme_bw()+
  theme(axis.text.x = element_text(face = "bold",color = "black",
                                   size = 14,angle = 90,hjust = 1,vjust = 0.3),
        axis.text.y = element_text(face = "bold",color = "black",size = 14),
        axis.ticks.length = unit(-0.3,"cm"),
        axis.title = element_text(size = 18,face = "bold"))+
  coord_fixed()+
  guides(fill=FALSE)
p

ggplot(iekpd2_box1,aes(x=factor(Tissue,levels = tissue),y=zscore_TPM))+
  geom_flat_violin(aes(fill=Tissue),position = position_nudge(x=.05),
                   color="black",adjust=0.5,width = 1.1)+
  # geom_jitter(aes(color=Tissue),width = 0.12,size=0.01)+
  geom_boxplot(aes(fill=Tissue),width=0.4,
               # outlier.color = NA,
               outlier.size = 1.5,outlier.shape = 21,
               position = position_nudge(x=-.2))+
  stat_summary(fun.y = "mean",geom = "point",
               fill="white",shape=21,size=2.5,
               position = position_nudge(x=-.2))+
  stat_compare_means(method = "anova",label.y = -2.5)+
  stat_compare_means(label = "p.signif",method = "t.test",ref.group = "Testis",
                     hide.ns = FALSE)+
  xlab("Tissue")+ylab("z-score PTM")+
  theme_bw()+
  theme(axis.text.x = element_text(face = "bold",color = "black",
                                   size = 14,angle = 90,hjust = 1,vjust = 0.3),
        axis.text.y = element_text(face = "bold",color = "black",size = 14),
        axis.ticks.length = unit(-0.3,"cm"),
        axis.title = element_text(size = 18,face = "bold"))+
  coord_fixed()+
  guides(fill=FALSE)


temp1<-subset(specific_genes,specific_genes$Type=='Writer')
temp2<-subset(specific_genes,specific_genes$Type=='Eraser')
temp3<-subset(specific_genes,specific_genes$Type=='Reader')
temp4<-subset(specific_genes,specific_genes$Type=='Multi')
Writer<-as.data.frame(table(temp1$Tissue))
colnames(Writer)<-c("tissue","Writer")
Eraser<-as.data.frame(table(temp2$Tissue))
colnames(Eraser)<-c("tissue","Eraser")
Reader<-as.data.frame(table(temp3$Tissue))
colnames(Reader)<-c("tissue","Reader")
Multi<-as.data.frame(table(temp4$Tissue))
colnames(Multi)<-c("tissue","Multi")

temp6<-join(Writer,Eraser,by="tissue")
temp5<-join(temp6,Reader,by="tissue")
temp7<-join(temp5,Multi,by="tissue")
rownames(temp7)<-as.character(temp7$tissue)
temp8<-temp7[,2:5]
gene_tissue_heatmap<-t(temp8)

gene_tissue_heatmap[is.na(gene_tissue_heatmap)]<-0

col_anno<-HeatmapAnnotation(Sum=anno_barplot(apply(gene_tissue_heatmap,2,sum)))
row_anno<-rowAnnotation(Sum=anno_barplot(apply(gene_tissue_heatmap,1,sum)))


pheatmap(gene_tissue_heatmap,
         cluster_rows = FALSE,cluster_cols = FALSE,
         cellwidth = 18,cellheight = 16,
         display_numbers = TRUE,number_format = "%.1f",number_color = "black",
         col = colorRamp2(c(0,8,16,24), c("floralwhite","LightSkyBlue","Tomato", "#d7191c"),space = "RGB"),
         # col = c("floralwhite","Bisque1","#2c7bb6", "#d7191c"),
         top_annotation=col_anno,left_annotation=row_anno)




library(FactoMineR)
library(factoextra)


iekpd2_pca<-t(na.omit(iEKPD2_tpm2_norm1))

iekpd_pca<-PCA(iekpd2_pca,graph = TRUE,scale.unit = FALSE)
var<-get_pca_var(iekpd_pca)
var
head(var$coord,4)

iekpd_ind<-get_pca_ind(iekpd_pca)

fviz_pca_ind(iekpd_pca,
             # col.ind = "cos2",
             # gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,
             col.ind = "blue",
             col.ind.sup = "red"
               # pointsize="cos2"
)+theme_bw()+
  geom_point(size=1.5)
dev.off()
