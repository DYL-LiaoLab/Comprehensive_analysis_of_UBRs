# Expression pattern and enrichment analysis of ubiquitination regulators at different developmental stages of testis
# Long deyu
# 2022.06.11

library(R.utils)
R.utils::setOption("clusterProfiler.download.method","auto")
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
library(clusterProfiler)
library(enrichplot)

human_tissue<-read.csv("mouse_ubiquitination_regulators.csv",check.names=FALSE)

testis<-human_tissue[,c(1,which(grepl("Testis",colnames(human_tissue))==TRUE))]
rownames(testis)<-testis[,1]
testis<-testis[,-1]

testis<-as.matrix(testis)
testis<-testis+0.0001

colnames(testis)<-gsub(".*\\.","",colnames(testis))
colnames(testis)

mfuzz_class<-new("ExpressionSet",exprs=testis)

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


library(RColorBrewer)
color.1<-rev(c("grey11","grey31","grey51","grey71","gray91"))
mfuzz.plot2(mfuzz_class,
            cl=mfuzz_cluster,
            mfrow = c(2,2),
            time.labels = colnames(testis),
            # colo = color.2,
            colo = color.1,
            # colo = rainbow(30, s = 3/4, v = 1, start = 1/6),
            # ax.col = color.2,
            centre.col = "DarkOrange",
            centre = TRUE,centre.lwd = 6
)

cluster_size<-mfuzz_cluster$size
names(cluster_size)<-1:cluster_num

testis_cluster<-mfuzz_cluster$cluster
testis_cluster<-cbind(testis[names(testis_cluster),],testis_cluster)
write.table(testis_cluster,"testis-mfuzz.txt",sep = "\t",col.names = NA,quote = FALSE)

head(mfuzz_cluster$cluster)

head(mfuzz_cluster$membership)

testis_cluster<-mfuzz_cluster$cluster
testis_standard<-mfuzz_class@assayData$exprs
testis_standard_cluster<-cbind(testis_standard[names(testis_cluster),],testis_cluster)
write.table(testis_standard_cluster,"testis-mfuzz-norm.txt",sep = "\t",col.names = NA,quote = FALSE)


retu<-as.data.frame(testis_standard_cluster)
retu<-retu[order(retu$testis_cluster),]

annotation_row<-subset(retu,select=15)

ann_colors= list(
  testis_cluster=c("Cyan1","Gold","Purple"),
  Type=c(PK="#7570B3",PP="#E7298A",PPBD="#66A61E"),
  State=c("pre-born"="lemonChiffon","post-born"="lavender")
)
retu<-retu[,1:14]
library(circlize)

p<-pheatmap(retu,
            name = "zscore RPKM",
            # column_title="The expression level of phosphorylation regulator in the developmental stage of human testis ",
            # row_title= "",
            # col = colorRamp2(c(-2.4,-1.2,0,1.2,2.4), c("#2c7bb6","#abd9e9","floralwhite","#fdae61", "#d7191c"),space = "RGB"),
            col = colorRamp2(c(-2,-1,0,1,2), c("#2c7bb6","#abd9e9","White","#FF8C69", "#FF0000"),space = "RGB"),
            clustering_distance_cols = "pearson",
            cluster_rows = FALSE,
            cluster_cols = FALSE,
            row_dend_side = "left", #Dendogram on the right side
            annotation_row = annotation_row,
            annotation_colors = ann_colors,
            # annotation_col = anno_col,
            # cellwidth = 18,
            # cellheight = 1,
            fontsize = 10,
            show_rownames = FALSE,
            split = annotation_row$testis_cluster, #Splitting by Class
)
p

library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Mm.eg.db)

id_trans<-bitr(
  rownames(testis_standard_cluster),
  fromType = "ENSEMBL",
  toType = c("SYMBOL","ENTREZID"),
  OrgDb = org.Mm.eg.db
)
colnames(id_trans)[1]<-"ensembl_id"

testis_standard_cluster<-as.data.frame(testis_standard_cluster)
#
testis_genelist<-data.frame(ensembl_id=rownames(testis_standard_cluster),
                            cluster=testis_standard_cluster$testis_cluster)

testis_genelist1<-left_join(testis_genelist,id_trans,by="ensembl_id")
# 
which(duplicated(testis_genelist1$ensembl_id))

# 
ego_go<-enrichGO(gene = testis_genelist1$ENTREZID,
                 OrgDb = org.Mm.eg.db,
                 ont = "ALL",#
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05,
                 readable = TRUE)

write.csv(data.frame(ego_go),"Testis_all_gene_enrichGO.csv",row.names = FALSE)

kk_all<-enrichKEGG(gene = testis_genelist1$ENTREZID,
                   organism = 'mmu', # 
                   pvalueCutoff = 0.05)
head(kk_all)
write.csv(data.frame(kk_all),"Testis_all_gene_enrichKEGG.csv",row.names = FALSE)


for(i in 1:cluster_num){
  # i=1
  temp<-subset(testis_genelist1,testis_genelist1$cluster==i)
  ego_go<-enrichGO(gene = temp$ENTREZID,
                   OrgDb = org.Mm.eg.db,
                   ont = "ALL",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)
  
  write.csv(data.frame(ego_go),
            file = paste("Testis_enrichGO_cluster",i,".csv",sep = "_"),
            row.names = FALSE)
  kk_all<-enrichKEGG(gene = temp$ENTREZID,
                     organism = 'mmu',
                     pvalueCutoff = 0.05)
  # head(kk_all)
  write.csv(data.frame(kk_all),
            file = paste("Testis_enrichKEGG_cluster",i,".csv",sep = "_"),
            row.names = FALSE)
  
}

