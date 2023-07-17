# Batch removal using harmony.
# Long deyu 
# 2022.08.16

library(harmony)
library(Seurat)
library(cowplot)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(reshape2)

scrna <- readRDS("scrna.rds")

head(scrna@meta.data)
cellinfo <- scrna@meta.data

scrna <- CreateSeuratObject(scrna@assays$RNA@counts, meta.data = cellinfo) # scrna@assays$RNA@counts为原始矩阵文件
scrna<-NormalizeData(scrna,normalization.method = "LogNormalize",
                     # scale.factor = 10000,
                     verbose = FALSE) # 
scrna<-FindVariableFeatures(scrna,selection.method = "vst",nfeatures = 3000) 

scale.genes<-rownames(scrna)
scrna<-ScaleData(scrna,
                 features = scale.genes, # 
                 verbose = FALSE) # 
# pca
scrna<-RunPCA(scrna,features = VariableFeatures(scrna),npcs = 50,verbose = FALSE)

p<-ElbowPlot(scrna,ndims = 50)+theme_bw()
p
ggsave("scatter_int_before.pdf",plot = p,
       width = 6,height = 5)
p<-DimPlot(scrna,reduction = "pca",group.by = "sample")+theme_bw()
p
ggsave("pca_sample_int_before.pdf",plot = p,
       width = 6,height = 5)
p<-DimPlot(scrna,reduction = "pca",group.by = "stage")+theme_bw()
p
ggsave("pca_stage_int_before.pdf",plot = p,
       width = 6,height = 5)

pc.num=1:30
scrna<-FindNeighbors(scrna,dims = pc.num)
scrna<-FindClusters(scrna,resolution = 0.5)
# tSNE
scrna<-RunTSNE(scrna,dims = pc.num)

scrna<-RunUMAP(scrna,dims = pc.num)
p<-DimPlot(scrna,reduction = "tsne",group.by = "stage")+theme_bw()
p
ggsave("tsne_stage_int_before.pdf",plot = p,
       width = 6,height = 5)
p<-DimPlot(scrna,reduction = "umap",group.by = "stage")+theme_bw()
p
ggsave("umap_stage_int_before.pdf",plot = p,
       width = 6,height = 5)
p<-DimPlot(scrna,reduction = "tsne",group.by = "sample")+theme_bw()
p
ggsave("tsne_sample_int_before.pdf",plot = p,
       width = 6,height = 5)
p<-DimPlot(scrna,reduction = "umap",group.by = "sample")+theme_bw()
p
ggsave("umap_sample_int_before.pdf",plot = p,
       width = 6,height = 5)

p<-DimPlot(scrna,reduction = "tsne",group.by = "CellType")+theme_bw()
p
ggsave("tsne_CellType_int_before.pdf",plot = p,
       width = 6,height = 5)
p<-DimPlot(scrna,reduction = "umap",group.by = "CellType")+theme_bw()
p
ggsave("umap_CellType_int_before.pdf",plot = p,
       width = 6,height = 5)
head(scrna@meta.data)


# Run Harmony

head(scrna@meta.data)
options(repr.plot.height = 2.5, repr.plot.width = 6)
scrna<-RunHarmony(scrna,group.by.vars = "sample",plot_convergence = TRUE, max.iter.harmony = 30)


scrna@reductions$harmony

harmony_embeddings<-Embeddings(scrna,"harmony")
harmony_embeddings[1:5, 1:5]

p<-ElbowPlot(scrna,ndims = 50,reduction = "harmony")+theme_bw()
p
ggsave("harmony_scatter_int_after.pdf",plot = p,
       width = 6,height = 5)

pc.num=1:30
scrna<-FindNeighbors(scrna,dims = pc.num,reduction = "harmony")
scrna<-FindClusters(scrna,resolution = 0.5)
scrna<-RunTSNE(scrna,reduction = "harmony",dims = pc.num)
scrna<-RunUMAP(scrna,reduction = "harmony",dims = pc.num)

p<-DimPlot(scrna,reduction = "harmony",group.by = "stage")+theme_bw()
p
ggsave("harmony_stage.pdf",plot = p,
       width = 6,height = 5)
p<-DimPlot(scrna,reduction = "harmony",group.by = "sample")+theme_bw()
p
ggsave("harmony_sample.pdf",plot = p,
       width = 6,height = 5)
p<-DimPlot(scrna,reduction = "tsne",group.by = "sample")+theme_bw()
p
ggsave("tsne_sample_harmony.pdf",plot = p,
       width = 6,height = 5)
p<-DimPlot(scrna,reduction = "umap",group.by = "sample")+theme_bw()
p
ggsave("umap_sample_harmony.pdf",plot = p,
       width = 8,height = 6)

p<-DimPlot(scrna,reduction = "tsne",group.by = "CellType",label = TRUE)+theme_bw()
p
ggsave("tsne_CellType_harmony.pdf",plot = p,
       width = 8,height = 6)
p<-DimPlot(scrna,reduction = "umap",group.by = "CellType",label = TRUE)+theme_bw()
p
ggsave("umap_CellType_harmony.pdf",plot = p,
       width = 8,height = 6)

p<-DimPlot(scrna,reduction = "tsne",group.by = "stage",label = TRUE)+theme_bw()
p
ggsave("tsne_stage_harmony.pdf",plot = p,
       width = 6,height = 5)
p<-DimPlot(scrna,reduction = "umap",group.by = "stage",label = TRUE)+theme_bw()
p
ggsave("umap_stage_harmony.pdf",plot = p,
       width = 6,height = 5)
# 
saveRDS(scrna, "scrna_Harmony.rds")


raw_counts<-GetAssayData(scrna,slot = "counts",assay = "RNA")
norm_counts<-GetAssayData(scrna,slot = "data",assay = "RNA")
scale_counts<-GetAssayData(scrna,slot = "scale.data",assay = "RNA")


pwer<-read.csv("human_ubiquitination_regulators.csv")


head(scrna@meta.data)

cell_type<-unique(scrna@meta.data$CellType)

cell_id2type<-data.frame(rownames(scrna@meta.data),scrna@meta.data$CellType)
colnames(cell_id2type)<-c("human_cell_id","cell_type")

for(i in 1:length(cell_type)){
  print(cell_type[i])
  a<-cell_type[i]
  temp<-subset(cell_id2type,cell_id2type$cell_type==a)
  temp1<-c()
  for(j in 1:length(temp$human_cell_id)){
    # print(i)
    for(k in 1:length(colnames(scale_counts))){
      if(temp$human_cell_id[j]==colnames(scale_counts)[k]){
        temp1<-c(temp1,k)
      }
    }
  }
  print(length(temp1))
  buffer<-subset(scale_counts,select = temp1)
  if(i==1){
    gene_mean<-data.frame(apply(buffer,1,mean))
    colnames(gene_mean)<-a
  }else{
    temp2<-apply(buffer,1,mean)
    gene_mean<-cbind(gene_mean,temp2)
    colnames(gene_mean)[i]<-a
  }
}


# head(scrna@meta.data)
# gene_mean_all<-AverageExpression(scrna,assays = "RNA" ,slot = "scale.data",
#                                  group.by = "CellType",verbose = TRUE)[[1]]

rownames(gene_mean)<-gsub("-","_",rownames(gene_mean))


gene_mean$gene_symbol<-rownames(gene_mean)
pwer_expre<-inner_join(pwer,gene_mean,by="gene_symbol")


head(cellinfo)
samp<-unique(cellinfo$sample)
type_cell<-unique(cellinfo$CellType)


bardata<-as.data.frame(matrix(data = 0,nrow = 8,ncol = 8))
rownames(bardata)<-samp
colnames(bardata)<-type_cell


for(i in 1:length(samp)){
  for(j in 1:length(type_cell)){
    buff<-0 # 
    for(k in 1:length(cellinfo$orig.ident)){
      
      if(cellinfo$sample[k]==samp[i] & cellinfo$CellType[k]==type_cell[j]){
        buff<-buff+1
      }
    }
    bardata[i,j]<-buff
  }
}

bardata$samp<-rownames(bardata)
bardata_long<-melt(bardata,id.vars = "samp")
bardata<-bardata[,-9]
samp_sum<-apply(bardata,1,sum)
head(bardata_long)
p<-ggplot(bardata_long,aes(x=samp,y=value,fill=variable))+
  geom_col(position = "fill")+
  scale_y_continuous(labels = scales::percent,expand = c(0,0))+
  scale_x_discrete(limits=c(samp[1:2],rev(samp[5:6]),samp[3:4],samp[7:8]))+
  scale_fill_manual(values = c("#00BE67","#FF61CC","#F8766D","#CD9600",
                                        "#7CAE00","#00BFC4","#00A9FF","#C77CFF"))+
                                          theme_bw()

p

