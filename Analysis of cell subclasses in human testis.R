# Analysis of cell subclasses in human testis.
# Long deyu
# 2022.08.20

library(ggplot2)
library(Seurat)
library(dplyr)
library(patchwork)
library(feature)
library(tidyverse)
library(ComplexHeatmap)
library(pheatmap)
library(viridis)

scrna<-readRDS("scrna_Harmony.rds")
head(scrna@meta.data)
table(scrna@meta.data$CellType)

scrna_sub<-scrna[,scrna@meta.data$CellType %in% c("Spermatocyte","Spermatid")]
table(scrna_sub@meta.data$CellType)

scrna_sub<-NormalizeData(scrna_sub,normalization.method = "LogNormalize", scale.factor = 1e4)

scrna_sub<-FindVariableFeatures(scrna_sub,selection.method = "vst",nfeatures = 2000)

scale.genes<-rownames(scrna_sub)
scrna_sub<-ScaleData(scrna_sub,features = scale.genes) 

scrna_sub<-RunPCA(scrna_sub,features = VariableFeatures(scrna_sub)) 
p<-ElbowPlot(scrna_sub,ndims = 50,reduction = "pca")+theme_bw()
dir.create("05human_subcluster")
pdf(file = "05human_subcluster/ElbowPlot_before.pdf",width = 6,height = 5)
p
dev.off()

pc.num=1:19
scrna_sub<-FindNeighbors(scrna_sub,dims = pc.num)
scrna_sub<-FindClusters(scrna_sub,resolution = 0.3)

table(scrna_sub@meta.data$seurat_clusters)
metadata<-scrna_sub@meta.data
cell_clusters<-data.frame(cell_ID=rownames(metadata),cluster_id=metadata$seurat_clusters)

# tsne
scrna_sub=RunTSNE(scrna_sub,dims = pc.num)
embed_tsne<-Embeddings(scrna_sub,"tsne")
write.csv(embed_tsne,"05human_subcluster/embed_tsne_before.csv")
plot1 = DimPlot(scrna_sub, reduction = "tsne") 
plot1+theme_bw()
ggsave("05human_subcluster/tSNE_before.pdf", plot = plot1+theme_bw(), width = 8, height = 7)
ggsave("05human_subcluster/tSNE_before.png", plot = plot1+theme_bw(), width = 8, height = 7)
#UMAP
scrna_sub <- RunUMAP(scrna_sub, dims = pc.num)
embed_umap <- Embeddings(scrna_sub, 'umap')
write.csv(embed_umap,'05human_subcluster/embed_umap.csv') 
plot2 = DimPlot(scrna_sub, reduction = "umap") 
plot2+theme_bw()
ggsave("05human_subcluster/UMAP_before.pdf", plot = plot2+theme_bw(), width = 8, height = 7)
ggsave("05human_subcluster/UMAP_before.png", plot = plot2+theme_bw(), width = 8, height = 7)
#合并tSNE与UMAP
plotc <- plot1+plot2+ plot_layout(guides = 'collect')
ggsave("05human_subcluster/tSNE_UMAP_before.pdf", plot = plotc+theme_bw(), width = 10, height = 5)
ggsave("05human_subcluster/tSNE_UMAP_before.png", plot = plotc+theme_bw(), width = 10, height = 5)

library(harmony)
options(repr.plot.height = 2.5, repr.plot.width = 6)
scrna<-RunHarmony(scrna_sub,group.by.vars = "sample",
                  plot_convergence = TRUE, max.iter.harmony = 30)

scrna_sub@reductions$harmony
harmony_embeddings<-Embeddings(scrna_sub,"harmony")
harmony_embeddings[1:5, 1:5]

p<-ElbowPlot(scrna_sub,ndims = 50,reduction = "harmony")+theme_bw()
p
ggsave("05human_subcluster/harmony_scatter_int_after.pdf",plot = p,
       width = 6,height = 5)

pc.num=1:30
scrna_sub<-FindNeighbors(scrna_sub,dims = pc.num,reduction = "harmony")
scrna_sub<-FindClusters(scrna_sub,resolution = 0.38)
scrna_sub<-RunTSNE(scrna_sub,reduction = "harmony",dims = pc.num)
scrna_sub<-RunUMAP(scrna_sub,reduction = "harmony",dims = pc.num)

p<-DimPlot(scrna_sub,reduction = "harmony",group.by = "stage")+theme_bw()
p
ggsave("05human_subcluster/harmony_stage.pdf",plot = p,
       width = 6,height = 5)
p<-DimPlot(scrna_sub,reduction = "harmony",group.by = "sample")+theme_bw()
p
ggsave("05human_subcluster/harmony_sample.pdf",plot = p,
       width = 6,height = 5)
p<-DimPlot(scrna_sub,reduction = "tsne",group.by = "sample")+theme_bw()
p
ggsave("05human_subcluster/tsne_sample_harmony.pdf",plot = p,
       width = 6,height = 5)
p<-DimPlot(scrna_sub,reduction = "umap",group.by = "sample")+theme_bw()
p
ggsave("05human_subcluster/umap_sample_harmony.pdf",plot = p,
       width = 6,height = 5)

p<-DimPlot(scrna_sub,reduction = "tsne",group.by = "CellType",label = TRUE)+theme_bw()
p
ggsave("05human_subcluster/tsne_CellType_harmony.pdf",plot = p,
       width = 6,height = 5)
p<-DimPlot(scrna_sub,reduction = "umap",group.by = "CellType",label = TRUE)+theme_bw()
p
ggsave("05human_subcluster/umap_CellType_harmony.pdf",plot = p,
       width = 6,height = 5)

p<-DimPlot(scrna_sub,reduction = "tsne",group.by = "stage",label = TRUE)+theme_bw()
p
ggsave("05human_subcluster/tsne_stage_harmony.pdf",plot = p,
       width = 6,height = 5)
p<-DimPlot(scrna_sub,reduction = "umap",group.by = "stage",label = TRUE)+theme_bw()
p
ggsave("05human_subcluster/umap_stage_harmony.pdf",plot = p,
       width = 6,height = 5)

saveRDS(scrna_sub, "05human_subcluster/scrna_sub_Harmony.rds")

scrna_sub<-readRDS("05human_subcluster/scrna_sub_Harmony.rds")
diff.wilcox<-FindAllMarkers(scrna_sub)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(all.markers, "05human_subcluster/diff_genes_wilcox_Harmony.csv", row.names = F)
write.csv(top10, "05human_subcluster/top10_diff_genes_wilcox_Harmony.csv", row.names = F)

p<-DimPlot(scrna_sub,reduction = "umap",label = TRUE)+theme_bw()
p
ggsave("05human_subcluster/human_umap_cluster.pdf",plot = p,width = 7,height = 6)
p<-DimPlot(scrna_sub,reduction = "tsne",label = TRUE)+theme_bw()
p
ggsave("05human_subcluster/human_tsne_cluster.pdf",plot = p,width = 7,height = 6)

head(scrna_sub@meta.data)

new.cluster.ids<-c("S4","S2","S1","S3","S2","SPG7","Z","L","P & D")
names(new.cluster.ids)<-levels(scrna_sub)
scrna_sub<-RenameIdents(scrna_sub,new.cluster.ids)
scrna_sub$Cell.type<-scrna_sub$RNA_snn_res.0.38
scrna_sub$Cell.type<-recode(scrna_sub$Cell.type,
                            "0"="S4","1"="S2","2"="S1","3"="S3",
                            "4"="S2","5"="SPG7","6"="Z","7"="L",
                            "8"="P & D")
DimPlot(scrna_sub,reduction = "umap",label = TRUE,group.by = "Cell.type")+theme_bw()
p<-DimPlot(scrna_sub,reduction = "umap",label = TRUE)+theme_bw()
p
ggsave("05human_subcluster/human_umap_celltype.pdf",plot = p,width = 8,height = 6)
p<-DimPlot(scrna_sub,reduction = "umap",label = TRUE)+theme_bw()
p
ggsave("05human_subcluster/human_umap_celltype-5.pdf",plot = p,width = 6,height = 5)
p<-DimPlot(scrna_sub,reduction = "tsne",label = TRUE)+theme_bw()+
  theme(panel.background = element_rect(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none')

p
ggsave("05human_subcluster/human_tsne_celltype.pdf",plot = p,width = 6,height = 6)

FeaturePlot(scrna_sub,reduction = "tsne",label = TRUE,
            features = "SCML1",cols = c("grey61","Khaki1","OrangeRed1"))+theme_bw()+
  theme(panel.background = element_rect(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

markgene<-c("SCML1","DPH7","TDRG1","CCNA1","TJP3","SUN5","TEX29","CCIN","NFKBIB","TNP1","IQCF3",
            "PRM1","LELP1")

for(i in 1:length(markgene)){
  print(markgene[i])
  p<-FeaturePlot(scrna_sub,
                 features = markgene[i],order = TRUE,
                 reduction = "tsne",label = TRUE,slot ="scale.data",
                 keep.scale = "feature",
                 cols = c("grey61","Khaki1","OrangeRed1"))+theme_bw()+
    theme(panel.background = element_rect(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  p
  ggsave(paste0("05human_subcluster/",markgene[i],"-featureplot.pdf",sep=""),p,
         width = 4,height = 4)
  ggsave(paste0("05human_subcluster/",markgene[i],"-featureplot.jpg",sep=""),p,
         width = 5,height = 4)
}


library(FlexDotPlot)
dp<-DotPlot(scrna_sub,
            features = markgene)+
  RotatedAxis()+theme_bw()+
  theme(panel.background = element_rect(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")+
  scale_y_discrete(limits=rev(c("L","Z","P & D","SPG7","S1","S2","S3","S4")))+ # 调整坐标轴的顺序
  scale_x_discrete(limits=markgene)+
  scale_color_gradientn(colours = c("#2c7bb6","White","#FF8C69", "#FF0000"))
dp
ggsave("05human_subcluster/human_markgene_dopplot.pdf",plot = dp,width = 12,height = 6)

human_cluster<-read.csv('human_cell_cluster.csv')
human_c1<-subset(human_cluster,human_cluster$cluster==1)
human_c2<-subset(human_cluster,human_cluster$cluster==2)
human_c6<-subset(human_cluster,human_cluster$cluster==6)


# data_scale<-GetAssayData(scrna_sub,slot = "scale.data")
human_c1_c2<-subset(human_cluster,human_cluster$cluster==1 | human_cluster$cluster==2 | human_cluster$cluster==6 )
human_c1_c2<-human_c1_c2[order(human_c1_c2$cluster),]

scrna_sub_c1c2<-subset(scrna_sub,features=human_c1_c2$gene_symbol)

data.scale<-GetAssayData(scrna_sub_c1c2,slot = "scale.data") 
cluster_info<-scrna_sub_c1c2@meta.data
cell_order<-c("L","Z","P & D","SPG7","S1","S2","S3","S4")
for(i in 1:length(cluster_info$Cell.type)){
  if(cluster_info$Cell.type[i]=="L"){
    cluster_info$Cell.type1[i]<-"A1"
  }else if(cluster_info$Cell.type[i]=="Z"){
    cluster_info$Cell.type1[i]<-"A2"
  }else if(cluster_info$Cell.type[i]=="P & D"){
    cluster_info$Cell.type1[i]<-"A3"
  }else if(cluster_info$Cell.type[i]=="SPG7"){
    cluster_info$Cell.type1[i]<-"A4"
    # }else if(cluster_info$Cell.type[i]=="MI"){
    #   cluster_info$Cell.type1[i]<-"A5"
  }else{
    cluster_info$Cell.type1[i]<-as.vector(cluster_info$Cell.type[i])
  }
}


for(i in 1:length(cell_order)){
  # i=6
  print(cell_order[i])
  buff<-subset(cluster_info,cluster_info$Cell.type==cell_order[i])
  if(i==1){
    cluster_order_info<-buff
  }else{
    cluster_order_info<-rbind(cluster_order_info,buff)
  }
}


data.scale1<-data.scale[,rownames(cluster_order_info)]
data.scale2<-data.scale1[human_c1_c2$gene_symbol,] 

anno_col<-subset(cluster_order_info,select = c(17))
identical(colnames(data.scale2),rownames(cluster_order_info))
# 
anno_row<-subset(human_c1_c2,select = c(1,2))
rownames(anno_row)<-human_c1_c2$gene_symbol

ann_colors= list(
  cluster=c("Coral1","PeachPuff3","LightBlue1",
            "Plum","OliveDrab3","RosyBrown2"),
  Type=c(Writer="#A7A0FC",Reader="#F3E1A8", Multi="#C3FFA7",Eraser="#CF9FFB")
)

library(circlize)
library(RColorBrewer)

p<-ComplexHeatmap::pheatmap(data.scale2,
                            use_raster=FALSE,
                            cluster_rows = FALSE,
                            cluster_cols = FALSE,
                            annotation_col = anno_col,
                            annotation_colors = ann_colors,
                            show_colnames = FALSE,
                            show_rownames = FALSE,
                            annotation_row = anno_row,
                            color = colorRamp2(c(-2,-1,0,1,2),c("NavyBlue","LightSkyBlue","gray81","Salmon1","Red1"),space = "RGB"),
                            row_split = anno_row$cluster,
                            column_split = anno_col$Cell.type1
)
pdf("05human_subcluster/human_heatmap-7.pdf",width=7)
p
dev.off()
png("05human_subcluster/human_heatmap.png",height = 600,width = 600)
p
dev.off()

human_data.scale<-as.data.frame(data.scale2)
human_data.scale_c2<-human_data.scale[human_c2$gene_symbol,]
p<-ComplexHeatmap::pheatmap(human_data.scale_c2,
                            use_raster=FALSE,
                            cluster_rows = FALSE,
                            cluster_cols = FALSE,
                            annotation_col = anno_col,
                            annotation_colors = ann_colors,
                            show_colnames = FALSE,
                            show_rownames = FALSE,
                            # annotation_row = anno_row,
                            color = colorRamp2(c(-2,-1,0,1,2),c("NavyBlue","LightSkyBlue","gray81","Salmon1","Red1"),space = "RGB"),
                            # row_split = anno_row$cluster,
                            column_split = anno_col$Cell.type1
)
p


human_c3<-subset(human_cluster,human_cluster$cluster==6)
human_data.scale_c3<-human_data.scale[human_c3$gene_symbol,]
p<-ComplexHeatmap::pheatmap(human_data.scale_c3,
                            use_raster=FALSE,
                            cluster_rows = FALSE,
                            cluster_cols = FALSE,
                            annotation_col = anno_col,
                            annotation_colors = ann_colors,
                            show_colnames = FALSE,
                            show_rownames = FALSE,
                            # annotation_row = anno_row,
                            color = colorRamp2(c(-2,-1,0,1,2),c("NavyBlue","LightSkyBlue","gray81","Salmon1","Red1"),space = "RGB"),
                            # row_split = anno_row$cluster,
                            column_split = anno_col$Cell.type1
)
p


data.scale_all<-GetAssayData(scrna_sub,slot = "scale.data")
cluster_info<-scrna_sub@meta.data
cell_order<-c("L","Z","P & D","SPG7","S1","S2","S3","S4")
for(i in 1:length(cluster_info$Cell.type)){
  if(cluster_info$Cell.type[i]=="L"){
    cluster_info$Cell.type1[i]<-"A1"
  }else if(cluster_info$Cell.type[i]=="Z"){
    cluster_info$Cell.type1[i]<-"A2"
  }else if(cluster_info$Cell.type[i]=="P & D"){
    cluster_info$Cell.type1[i]<-"A3"
  }else if(cluster_info$Cell.type[i]=="SPG7"){
    cluster_info$Cell.type1[i]<-"A4"
    # }else if(cluster_info$Cell.type[i]=="MI"){
    #   cluster_info$Cell.type1[i]<-"A5"
  }else{
    cluster_info$Cell.type1[i]<-as.vector(cluster_info$Cell.type[i])
  }
}


for(i in 1:length(cell_order)){
  # i=6
  print(cell_order[i])
  buff<-subset(cluster_info,cluster_info$Cell.type==cell_order[i])
  if(i==1){
    cluster_order_info<-buff
  }else{
    cluster_order_info<-rbind(cluster_order_info,buff)
  }
}

data.scale_all<-as.data.frame(data.scale_all)
data.scale_all1<-data.scale_all[,rownames(cluster_order_info)]
# data.scale_all<-as.data.frame(data.scale_all)
data.scale2<-data.scale_all1[markgene,] #

p<-ComplexHeatmap::pheatmap(data.scale2,
                            use_raster=FALSE,
                            cluster_rows = FALSE,
                            cluster_cols = FALSE,
                            annotation_col = anno_col,
                            annotation_colors = ann_colors,
                            show_colnames = FALSE,
                            # show_rownames = FALSE,
                            # annotation_row = anno_row,
                            color = colorRamp2(c(-2,-1,0,1,2),c("NavyBlue","LightSkyBlue","gray81","Salmon1","Red1"),space = "RGB"),
                            # row_split = anno_row$cluster,
                            column_split = anno_col$Cell.type1
)
p


for(i in 1:length(scrna_sub@meta.data$Cell.type)){
  if(scrna_sub@meta.data$Cell.type[i]=="L"){
    scrna_sub@meta.data$Cell.type1[i]<-"A1"
  }else if(scrna_sub@meta.data$Cell.type[i]=="Z"){
    scrna_sub@meta.data$Cell.type1[i]<-"A2"
  }else if(scrna_sub@meta.data$Cell.type[i]=="P & D"){
    scrna_sub@meta.data$Cell.type1[i]<-"A3"
    # }else if(scrna_sub@meta.data$Cell.type[i]=="D"){
    #   scrna_sub@meta.data$Cell.type1[i]<-"A4"
  }else if(scrna_sub@meta.data$Cell.type[i]=="SPG7"){
    scrna_sub@meta.data$Cell.type1[i]<-"A4"
  }else{
    scrna_sub@meta.data$Cell.type1[i]<-as.vector(scrna_sub@meta.data$Cell.type[i])
  }
}


p<-DotPlot(scrna_sub,
           features =markgene,
           group.by = "Cell.type1",
           cols = c("DarkBlue","Tomato1"),
           dot.min = 0,
           dot.scale = 8,
           col.min = -1,col.max = 1,
)+theme_bw()+
  # RotatedAxis()+
  scale_fill_gradientn(colours = c("NavyBlue","LightSkyBlue","gray81","Salmon1","Red1"))
p

library(hrbrthemes)
library(patchwork)
library(GGally)
library(tidyverse)
library(dplyr)
library(reshape2)
library(ggplot2)
library(hexbin)

human_cluster<-read.csv('human_cell_cluster.csv')
human_mean<-AverageExpression(scrna_sub,group.by = "Cell.type",
                              slot = "scale.data")
human_mean<-as.data.frame(human_mean$RNA)

head(human_mean)

ub_human_mean<-human_mean[human_cluster$gene_symbol,]
ub_human_mean<-ub_human_mean[,c(7,6,8,5,3,2,4,1)]


for(i in 1:6){
  # i=6
  if(i==3 | i==4 | i==5){next}
  temp_cluster<-subset(human_cluster,human_cluster==i)
  temp_ub_mean<-ub_human_mean[temp_cluster$gene_symbol,]
  ggparcoord(temp_ub_mean,
             columns = 1:8,
             alphaLines = 0.3)+
    geom_boxplot()
  temp_ub_mean$gene<-rownames(temp_ub_mean)
  temp_ub_chang<-melt(temp_ub_mean,id.vars = "gene") 
  colnames(temp_ub_chang)
  
  p<-ggplot(temp_ub_chang,aes(variable,value))+
    # geom_violin(scale = "width")+
    stat_bin_hex()+
    scale_fill_gradient(low = "#B0E2FF",high = "#FF00FF")+
    geom_point(position = position_jitter(width = 0.36,height = 0),colour="grey71",shape=24,size=0.6)+
    theme_bw()+
    ggtitle(paste0("group",i))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    stat_summary(fun = "mean",geom = "point",shape=24,
                 size=2.5,fill="#7FFF00")+
    stat_summary(fun = "median",geom = "point",shape=25,
                 size=2.5,fill="#FFD700")
  p
  ggsave(paste0("05human_subcluster/","human_cluster",i,".pdf"),plot = p,
         width = 6,height = 5.75)
  ggsave(paste0("05human_subcluster/","human_cluster",i,".jpg"),plot = p,
         width = 6,height = 5.75)
}
