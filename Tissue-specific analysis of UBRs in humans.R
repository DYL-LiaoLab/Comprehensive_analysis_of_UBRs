# Tissue-specific analysis of ubiquitination regulators in humans(GTEx).
# Long deyu
# 2022.05.16


library(geneRal)
library(dplyr)
library(openxlsx)
library(writexl)
library(readr)
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
library(data.table)


ub_wer<-read.csv("human_ubiquitination_regulators.csv")

dir.create("GTEx")

table(ub_wer$Type)
gene_type<-table(ub_wer$Type)
gene_type<-as.data.frame(gene_type)
colnames(gene_type)<-c("Type","Count")
gene_type<-gene_type[order(gene_type$Type,decreasing = TRUE),]

cols_pie<-c(Writer="#A7A0FC",Reader="#F3E1A8", Multi="#C3FFA7",Eraser="#CF9FFB"  )

labs<-paste0(gene_type$Type,"(",round(gene_type$Count/sum(gene_type$Count)*100,2),"%)")

pie(gene_type$Count,labels = labs,col = cols_pie,
    border = "black",
    edges = 10000,
    init.angle = 90,
    clockwise = TRUE,
    lty = 1)



gtex<-fread("rna_tissue_gtex.tsv")
colnames(gtex)
colnames(ub_wer)
colnames(gtex)[2]<-"gene_symbol"
gtex<-gtex[,c(2,3,6)]


ub_gtex<-left_join(ub_wer,gtex,by="gene_symbol") 
ub_gtex_raw<-ub_gtex
ub_gtex<-na.omit(ub_gtex)  
length(unique(ub_gtex$ensembl_id))
length(unique(ub_gtex$gene_symbol))
summary(ub_gtex$nTPM)


ub_gtex1<-ub_gtex

ub_gtex1<-ub_gtex1[,c(3,5,1,6,7)]
colnames(ub_gtex1)<-c("Gene_name","Type","Ensembl_Gene_ID","Tissue","TPM_value")

ub_gtex1$TPM_value<-log(ub_gtex1$TPM_value+1)


genemean<-aggregate(ub_gtex1[,5],list(ub_gtex1$Ensembl_Gene_ID),mean)
colnames(genemean)<-c("Ensembl_Gene_ID","genemean")

ub_gtex2<-join(ub_gtex1,genemean,by="Ensembl_Gene_ID")

ub_gtex2$abs<-ub_gtex2$TPM_value - ub_gtex2$genemean


iekpd_scatter2<-ub_gtex2
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
  pdf(file=paste("GTEx/",tissue,"specificityplot_logtpm.pdf",sep="_"),height=5,width=5)
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
  
  png(filename = paste("GTEx/",tissue,"specificityplot_logtpm.png",sep="_"))
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
}

write.table(specific_genes, file="GTEx/tissue_specific_genes_gtex_log(tpm+1)_human.tsv",quote=FALSE, row.names=FALSE,sep="\t")


ub_gtex_kuan<-ub_gtex
colnames(ub_gtex_kuan)
# ub_gtex_kuan<-spread(ub_gtex_kuan,key = gene_symbol,value = nTPM)
ub_gtex_kuan<-reshape2::dcast(ub_gtex_kuan,gene_symbol+Family+Type+ensembl_id+Species~Tissue,
                              value.var = "nTPM",
                              fun.aggregate = mean
)
ub_gtex_kuan<-na.omit(ub_gtex_kuan)

ub_gtex_kuan<-ub_gtex_kuan[order(ub_gtex_kuan$Type),]
rownames(ub_gtex_kuan)<-ub_gtex_kuan$gene_symbol

temp1<-ub_gtex_kuan[,-c(1:5)]
ub_gtex_kuan_norm<-t(apply(temp1,1,scale)) 
colnames(ub_gtex_kuan_norm)<-colnames(temp1)

annotation_row<-as.data.frame(ub_gtex_kuan[,c(2,3)])
annotation_row<-subset(annotation_row,select = 2)
ann_colors=list(Type=c(Writer="#A7A0FC",Reader="#F3E1A8", Multi="#C3FFA7",Eraser="#CF9FFB"))
# colnames(annotation_row)[1]<-"Type"


pheatmap(ub_gtex_kuan_norm,
         # 标准化
         # scale = "row",
         # 聚类方法
         # clustering_distance_rows = "pearson",
         # clustering_distance_rows ="spearman",
         clustering_distance_rows ="euclidean",
         annotation_row = annotation_row,
         annotation_colors = ann_colors,
         cluster_rows = FALSE,
         # cellwidth = 12,
         # cellheight = 2,
         show_rownames = FALSE,
         # col = colorRamp2(c(-4,-2,0,2,4), c("#383899","#1E90FF","White","#FF8C69", "#FF0000"),space = "RGB"),
         # col = colorRamp2(c(-3.6,-1.8,0,1.8,3.6), c("#2c7bb6","#abd9e9","floralwhite","#fdae61", "#d7191c"),space = "RGB"),
         # col = colorRamp2(c(-4,-2,0,2,4), c("#2c7bb6","#abd9e9","floralwhite","#fdae61", "#d7191c"),space = "RGB"),
         # col = colorRamp2(c(-2,-1,0,1,2), c("#2c7bb6","#abd9e9","floralwhite","#fdae61", "#d7191c"),space = "RGB"),
         col = colorRamp2(c(-4,-2,0,2,4), c("#2c7bb6","#abd9e9","White","#FF8C69", "#FF0000"),space = "RGB"),
         name = "z-scale nTPM",
         # split=annotation_row$Type
)



ub_gtex_kuan_norm<-as.data.frame(ub_gtex_kuan_norm)
ub_gtex_box<-cbind(ub_gtex_kuan[,c(1:5)],ub_gtex_kuan_norm)

colnames(ub_gtex_box)
ub_gtex_box1<-reshape2::melt(ub_gtex_box,c("gene_symbol","Family","Type","ensembl_id","Species" ))
colnames(ub_gtex_box1)[6:7]<-c("Tissue","zscore_nTPM")
#
library(ggpubr)
compare_means(zscore_nTPM~Tissue,data = ub_gtex_box1,
              ref.group = "testis",method = "t.test")

p<-ggplot(ub_gtex_box1,aes(x=Tissue,y=zscore_nTPM))+
  # position = position_nudge(x=.05)
  geom_flat_violin(aes(fill=Tissue),position = position_nudge(x=.05),
                   # 
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
  stat_compare_means(label = "p.signif",method = "t.test",ref.group = "testis",
                     hide.ns = FALSE)+
  # stat_compare_means(label = "p.signif",method = "t.test",ref.group = ".all.",
  #                    hide.ns = FALSE)+
  xlab("")+
  ylab("z-score nTPM")+
  theme_bw()+
  theme(axis.text.x = element_text(face = "bold",color = "black",
                                   size = 14,angle = 45,hjust = 1,vjust = 1),
        axis.text.y = element_text(face = "bold",color = "black",size = 14),
        axis.ticks.length = unit(-0.3,"cm"),
        axis.title = element_text(size = 18,face = "bold"))+
  coord_fixed()+
  guides(fill=FALSE)

p


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
#
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
         col = colorRamp2(c(0,20,40,60), c("floralwhite","LightSkyBlue","Tomato", "#d7191c"),space = "RGB"),
         # col = c("floralwhite","Bisque1","#2c7bb6", "#d7191c"),
         top_annotation=col_anno,left_annotation=row_anno)

pheatmap(gene_tissue_heatmap,
         cluster_rows = FALSE,cluster_cols = FALSE,
         cellwidth = 18,cellheight = 16,
         display_numbers = TRUE,number_format = "%.1f",number_color = "black",
         col = colorRamp2(c(0,20,40,60), c("floralwhite","LightSkyBlue","Tomato", "#d7191c"),space = "RGB"),
         # col = c("floralwhite","Bisque1","#2c7bb6", "#d7191c"),
         top_annotation=col_anno,left_annotation=row_anno)



# install.packages("FactoMineR")
# install.packages("factoextra")
library(FactoMineR)
library(factoextra)


ub_gtex_kuan_norm1<-na.omit(ub_gtex_kuan_norm)
ub_tissue_pca<-t(ub_gtex_kuan_norm1)

ub_tissue_pca1<-PCA(ub_tissue_pca,graph = TRUE,scale.unit = FALSE)

var<-get_pca_var(ub_tissue_pca1)
var
head(var$coord,4)

ub_ind<-get_pca_ind(ub_tissue_pca1)
ub_ind

fviz_pca_ind(ub_tissue_pca1,
             # col.ind = "cos2",
             # gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,
             col.ind = "blue",
             col.ind.sup = "red",
             
             # addEllipses = TRUE
             # alpha.ind = 5
             # pointsize="cos2"
)+theme_bw()+
  geom_point(size=1.5)

