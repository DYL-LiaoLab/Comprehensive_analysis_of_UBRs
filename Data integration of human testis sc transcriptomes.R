# Integrative analysis of the human testis single-cell transcriptome.
# Long deyu
# 2022.08.16


library(Seurat)
library(harmony)
library(tidyverse)
library(patchwork)
library(dplyr)
library(reshape2)
library(readxl)


infant<-read.csv("GSE120506_infant_combined_UMI.txt",sep = "\t",check.names = FALSE,
                 header = TRUE,row.names = 1,stringsAsFactors = FALSE)
pubertal<-read.csv("GSE134144_Pubertal_combined_UMI.txt",sep = "\t",check.names = FALSE,
                   header = TRUE,row.names = 1,stringsAsFactors = FALSE)
adult<-read.csv("GSE112013_Adult_Combined_UMI_table.txt",sep = "\t",check.names = FALSE,
                header = TRUE,row.names = 1,stringsAsFactors = FALSE)

celltype<-read_xlsx("mmc2 (2).xlsx")
celltype$cell_id<-gsub("-P\\d*-","-",celltype$CellID)

table(gsub("-.*","",celltype$CellID))
infant_raw<-data.frame(cellid_raw=colnames(infant),
                       sample=gsub("-.*","",colnames(infant)),
                       cell_id=gsub("Infant\\d*-","Y1-",colnames(infant)))
pubertal_raw<-data.frame(cellid_raw=colnames(pubertal),
                         sample=gsub("-.*","",colnames(pubertal)),
                         cell_id=colnames(pubertal))
adult_raw<-data.frame(cellid_raw=colnames(adult),
                      sample=gsub("-.*","",colnames(adult)),
                      cell_id=gsub("Donor\\d*-","Y25-",colnames(adult)))

intersect(celltype$cell_id,colnames(pubertal))
pubertal<-subset(pubertal,select=intersect(celltype$cell_id,colnames(pubertal)))

head(colnames(infant))
colnames(infant)<-gsub("Infant\\d*-","Y1-",colnames(infant))
intersect(celltype$cell_id,colnames(infant))

table(gsub("-.*","",celltype$CellID))
head(colnames(adult))
# colnames(adult)<-gsub("Donor1-","Y17-",colnames(adult))
# colnames(adult)<-gsub("Donor2-","Y24-",colnames(adult))
# colnames(adult)<-gsub("Donor3-","Y25-",colnames(adult))
colnames(adult)<-gsub("Donor\\d*-","Y25-",colnames(adult))
intersect(celltype$cell_id,colnames(adult)) # 
adult<-subset(adult,select=intersect(celltype$cell_id,colnames(adult))) # 
# 
infant_sc<-CreateSeuratObject(infant,project = "infant")
pubertal_sc<-CreateSeuratObject(pubertal,project = "pubertal")
adult_sc<-CreateSeuratObject(adult,project = "adult")

# 
scrna<-merge(infant_sc,c(pubertal_sc,adult_sc))

# 
table(scrna$orig.ident)
summary(scrna)
#
human_metadata<-scrna@meta.data
human_metadata_raw<-human_metadata
human_metadata$cell_id<-rownames(human_metadata)
head(human_metadata)
human_metadata<-left_join(human_metadata,celltype,by="cell_id")

which(duplicated(human_metadata$cell_id))

cell_remove<-human_metadata$cell_id[which(duplicated(human_metadata$cell_id))]
cell_retain<-setdiff(human_metadata$cell_id,cell_remove) 

human_cell_retain<-human_metadata[1,]
for(i in 1:length(cell_retain)){
  # i=1
  print(i)
  buff<-subset(human_metadata,human_metadata$cell_id==cell_retain[i])
  human_cell_retain<-rbind(human_cell_retain,buff)
}
human_cell_retain<-human_cell_retain[-1,] # 


human_raw<-rbind(infant_raw,pubertal_raw)
human_raw<-rbind(human_raw,adult_raw)

human_cell<-inner_join(human_cell_retain,human_raw,by="cell_id")
which(duplicated(human_cell$cell_id))

cell_remove1<-human_cell$cell_id[which(duplicated(human_cell$cell_id))]
cell_retain1<-setdiff(human_cell$cell_id,cell_remove1)

human_cleancell<-human_cell[1,]
for(i in 1:length(cell_retain1)){
  # i=1
  print(i)
  buff<-subset(human_cell,human_cell$cell_id==cell_retain1[i])
  human_cleancell<-rbind(human_cleancell,buff)
}
human_cleancell<-human_cleancell[-1,]
human_cleancell<-human_cleancell[,c(1:5,10:12)]
rownames(human_cleancell)


scrna<-subset(scrna,cells=human_cleancell$cell_id)
rownames(human_cleancell)<-human_cleancell$cell_id

identical(rownames(scrna@meta.data),human_cleancell$cell_id)
head(scrna@meta.data)
table(scrna@meta.data$orig.ident)
identical(rownames(scrna@meta.data),rownames(human_cleancell))
human_cleancell$stage<-human_cleancell$orig.ident

head(scrna@meta.data)
head(human_cleancell)
scrna@meta.data<-human_cleancell
head(scrna@meta.data)

scrna[["percent.mt"]]<-PercentageFeatureSet(scrna,pattern = "^MT-")
head(scrna@meta.data)

violin<-VlnPlot(scrna,features = c("nFeature_RNA","nCount_RNA","percent.mt"),
                # cols = rainbow(col.num),
                pt.size = 0.01, #
                ncol = 1)+
  theme(axis.title = element_blank(),axis.title.x = element_blank(),
        axis.ticks = element_blank())
violin

ggsave("vlnplot_qc.pdf",plot=violin, width = 9, height = 8)


theme.set2 = theme(axis.title.x=element_blank())

plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt")
group = "orig.ident"

plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scrna, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=3)   
violin
ggsave("vlnplot_qc_style1.pdf",plot=violin, width = 9, height = 8)


# 
scrna <- NormalizeData(scrna) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData()
scrna <- RunPCA(scrna, verbose = F)
p<-ElbowPlot(scrna, ndims = 50)
ggsave("scatter.pdf",plot=p, width = 6, height = 5)

pc.num=1:30 
scrna <- scrna %>% RunTSNE(dims=pc.num) %>% RunUMAP(dims=pc.num)
scrna <- FindNeighbors(scrna, dims=pc.num) %>% FindClusters()
head(scrna@meta.data)
p<-DimPlot(scrna, group.by = "stage",label = T)+theme_bw()
p
ggsave("before_umap_stage.pdf",plot=p, width = 6, height = 5)
p<-DimPlot(scrna, group.by = "sample",label = T)+theme_bw()
p
ggsave("before_umap_sample.pdf",plot=p, width = 6, height = 5)
p<-DimPlot(scrna, group.by = "stage", 
           split.by = "stage", ncol = 2)+theme_bw()
p
ggsave("before_umap_stage_splite.pdf",plot=p, width = 12, height = 10)
p<-DimPlot(scrna, group.by = "sample", 
           split.by = "sample", ncol = 3)+theme_bw()
p
ggsave("before_umap_sample_splite.pdf",plot=p, width = 12, height = 10)
p<-DimPlot(scrna,group.by = "CellType",label = TRUE)+theme_bw()
p
ggsave("before_umap_CellType.pdf",plot=p, 
       width = 6, height = 5)
p<-DimPlot(scrna,group.by = "CellType",label = TRUE,reduction = "tsne")+theme_bw()
p
ggsave("before_tsne_CellType.pdf",plot=p, 
       width = 6, height = 5)
# 
saveRDS(scrna, "scrna.rds")

