# Expression pattern and enrichment analysis of ubiquitination regulators at different developmental stages of testis
# Long deyu
# 2022.06.11

library(marray)
library(Seurat)
library(Mfuzz)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(NbClust)
library(vegan)
library(cluster)
library(factoextra)
library(pheatmap)
library(ComplexHeatmap)
library(cluster)

tissue_mean<-readxl::read_excel("Mean value and sample size of all human tissues.xlsx",
                                sheet = "average for all human tissues")
# 
ub_wer<-read.csv("human_ubiquitination_regulators.csv")

colnames(ub_wer)
colnames(tissue_mean)
colnames(tissue_mean)[1]<-"ensembl_id"
ub_wer_tissue_mean<-inner_join(ub_wer,tissue_mean,by="ensembl_id")


length(unique(ub_wer_tissue_mean$gene_symbol))
length(unique(ub_wer_tissue_mean$ensembl_id))

ub_wer_tissue_mean<-ub_wer_tissue_mean[order(ub_wer_tissue_mean$Type,ub_wer_tissue_mean$Family),]

tissue<-unique(gsub("\\..*","",colnames(ub_wer_tissue_mean)[6:139]))
tissue

i=7
print(tissue[i])

buffer<-ub_wer_tissue_mean[,c(1:5,which(grepl(tissue[i],colnames(ub_wer_tissue_mean))==TRUE))]
rownames(buffer)<-buffer$gene_symbol

temp1<-buffer[,-c(1:5)]
temp1<-na.omit(temp1)
temp1<-as.matrix(temp1)
temp1<-temp1+0.0001


colnames(temp1)<-gsub(".*\\.","",colnames(temp1))
colnames(temp1)[14:21]<-c("IF","TD","YT","OT","YA","YM","OM","SI")
colnames(temp1)

mfuzz_class<-new("ExpressionSet",exprs=temp1)
mfuzz_class<-filter.NA(mfuzz_class,thres = 0.25)
mfuzz_class<-fill.NA(mfuzz_class,mode = "mean")
mfuzz_class<-filter.std(mfuzz_class,min.std = 0)

mfuzz_class<-standardise(mfuzz_class) 

set.seed(1000)
wss <- (nrow(mfuzz_class@assayData$exprs)-1)*sum(apply(mfuzz_class@assayData$exprs,2,var)) # 计算离均差平方和
for (i in 2:15) wss[i] <- sum(kmeans(mfuzz_class@assayData$exprs, 
                                     centers=i,nstart = 100)$withinss) #计算不同聚类个数的组内平方和

plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")


set.seed(1000)
cluster_num<-4
mfuzz_cluster<-mfuzz(mfuzz_class,c=cluster_num,m=mestimate(mfuzz_class))
head(mfuzz_cluster$membership)

library(RColorBrewer)
# color.2 <- colorRampPalette(rev(c("#d50000", "#ffea00", "#00e676")))(10)
color.1<-rev(c("grey11","grey71","gray91","grey51","grey31"))

mfuzz.plot2(mfuzz_class,
            cl=mfuzz_cluster,
            mfrow = c(2,2),
            time.labels = colnames(temp1),
            # colo = rev(mycolor),
            colo = color.1,
            # colo=brewer.pal(9,"Greys")[1:9],
            # colo = rainbow(10, s = 3/4, v = 1, start = 1/6,alpha = 0.5),
            # ax.col = color.2,
            centre.col = '#FF7256',
            centre = TRUE,centre.lwd = 6
)

cluster_size<-mfuzz_cluster$size
names(cluster_size)<-1:cluster_num

head(mfuzz_cluster$cluster)
head(mfuzz_cluster$membership)

temp1_cluster<-mfuzz_cluster$cluster
temp1_standard<-mfuzz_class@assayData$exprs
temp1_standard_cluster<-cbind(temp1_standard[names(temp1_cluster),],temp1_cluster)
write.table(temp1_standard_cluster,"testis-mfuzz-norm.txt",sep = "\t",col.names = NA,quote = FALSE)

colnames(ub_wer)
temp1_standard_cluster<-as.data.frame(temp1_standard_cluster)
temp1_standard_cluster$gene_symbol<-rownames(temp1_standard_cluster)
temp1_standard_cluster1<-inner_join(temp1_standard_cluster,ub_wer,by="gene_symbol")
rownames(temp1_standard_cluster1)<-temp1_standard_cluster1$gene_symbol

retu<-as.data.frame(temp1_standard_cluster1)
retu<-retu[order(retu$temp1_cluster,retu$Type),]

annotation_row<-subset(retu,select=c(22,27))

ann_colors= list(
  temp1_cluster=c("Cyan1","Gold","Purple","Sienna1","OliveDrab1"),
  Type=c(Writer="#A7A0FC",Reader="#F3E1A8", Multi="#C3FFA7",Eraser="#CF9FFB"),
  State=c("pre-born"="lemonChiffon","post-born"="lavender")
)
retu<-retu[,1:21]
library(circlize)

p<-ComplexHeatmap::pheatmap(retu,
                            name = "z-scale RPKM",
                            # column_title="The expression level of phosphorylation regulator in the developmental stage of human temp1 ",
                            # col = colorRamp2(c(-2.4,-1.2,0,1.2,2.4), c("#2c7bb6","#abd9e9","floralwhite","#fdae61", "#d7191c"),space = "RGB"),
                            col = colorRamp2(c(-2,-1,0,1,2), c("#2c7bb6","#abd9e9","White","#FF8C69", "#FF0000"),space = "RGB"),
                            clustering_distance_cols = "pearson",
                            cluster_rows = FALSE,
                            cluster_cols = FALSE,
                            row_dend_side = "left", #Dendogram on the right side
                            annotation_row = annotation_row,
                            annotation_colors = ann_colors,
                            fontsize = 10,
                            show_rownames = FALSE,
                            split = annotation_row$temp1_cluster, #Splitting by Class
)



library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DO.db)

id_trans<-bitr(
  rownames(temp1_standard_cluster),
  fromType = "SYMBOL",
  toType = c("ENSEMBL","ENTREZID"),
  OrgDb = org.Hs.eg.db
)
colnames(id_trans)[1]<-"gene_symbol"
temp1_standard_cluster<-as.data.frame(temp1_standard_cluster)
# 
temp1_genelist<-data.frame(gene_symbol=rownames(temp1_standard_cluster),
                           cluster=temp1_standard_cluster$temp1_cluster)

temp1_genelist1<-inner_join(temp1_genelist,id_trans,by="gene_symbol")
temp1_genelist1<-temp1_genelist1[,-3]

temp1_genelist1<-temp1_genelist1[!duplicated(temp1_genelist1),]

which(duplicated(temp1_genelist1$gene_symbol))


ego_go<-enrichGO(gene = temp1_genelist1$ENTREZID,
                 OrgDb = org.Hs.eg.db,
                 ont = "ALL",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05,
                 readable = TRUE)

write.csv(data.frame(ego_go),"human_testis_ub_wer_enrichGO.csv",row.names = FALSE)


for(i in 1:cluster_num){
  # i=1
  temp<-subset(temp1_genelist1,temp1_genelist1$cluster==i)

  ego_go<-enrichGO(gene = temp$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "ALL",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)
  
  write.csv(data.frame(ego_go),
            file = paste("Testis_enrichGO_cluster",i,".csv",sep = ""),
            row.names = FALSE)

}




