# K-means clustering analysis of cells in human testis.
# Long deyu
# 2022.08.18

library(ComplexHeatmap)
library(openxlsx)
library(ggplot2)
library(dplyr)
library(reshape2)
library(tidyverse)
library(circlize)

# 
human_cell_mean<-read.csv("cell-ubiquitination-mean.csv",check.name=FALSE,row.names = 1)
rownames(human_cell_mean)<-human_cell_mean$gene_symbol
human_cell_mean<-human_cell_mean[,c(5,3,2,7,13,12,6,10,8:9,11)]
# colnames(human_cell_mean)[1]<-"Type"

pwer_cell_mean<-human_cell_mean
heatmap<-human_cell_mean

anno_row<-heatmap[,1:3]
rownames(anno_row)<-anno_row$gene_symbol
anno_row<-subset(anno_row,select = c(1,3))

summary(heatmap)
retu<-heatmap[,-c(1:3)]
pheatmap(retu,
         annotation_row = anno_row,
         # cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         show_colnames = TRUE,
         # col = colorRamp2(c(-3,-1.5,0,1.5,3), c("#2c7bb6","#abd9e9","floralwhite","#fdae61", "#d7191c"),space = "RGB"),
         # col = colorRamp2(c(-3.6,-1.8,0,1.8,3.6), c("#2c7bb6","#abd9e9","floralwhite","#fdae61", "#d7191c"),space = "RGB"),
         # col = colorRamp2(c(-4,-2,0,2,4), c("#2c7bb6","#abd9e9","floralwhite","#fdae61", "#d7191c"),space = "RGB"),
         col = colorRamp2(c(-1,-0.5,0,0.5,1), c("#2c7bb6","#abd9e9","floralwhite","#fdae61", "#d7191c"),space = "RGB"),
         name = "Expression levels",
         # column_split=anno_col$Tissue
         # split=anno_col$tissue
)

library(cluster)
library(factoextra)
fviz_nbclust(retu,kmeans,method = "wss")+
  geom_vline(xintercept = 4,linetype=2)
wss <- (nrow(retu)-1)*sum(apply(retu,2,var)) 
for (i in 2:15) wss[i] <- sum(kmeans(retu, 
                                     centers=i)$withinss) 
# 绘图
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

set.seed(1000)
cluster<-as.vector(kmeans(retu,6,nstart = 10000)$cluster) 
data_testis<-cbind(cluster,heatmap) 
data_testis<-as.data.frame(data_testis)
table(data_testis$cluster)

data_testis<-data_testis[order(data_testis$cluster,data_testis$Type),]
annotation_row<-data_testis[,c(1,2)]
colnames(annotation_row)<-c("Group","Type")
annotation_row$Group<-paste("group",annotation_row$Group,sep = "")

ann_colors= list(
  Group=c(group1="Coral1",group2="PeachPuff3",group3="LightBlue1",
          group4="Plum",group5="OliveDrab3",group6="RosyBrown2"),
  Type=c(Writer="#A7A0FC",Reader="#F3E1A8", Multi="#C3FFA7",Eraser="#CF9FFB")
)
library(circlize)
p<-ComplexHeatmap::pheatmap(data_testis[,-c(1:4)],
                            name = "z-scale expression",
                            # column_title="The expression level of phosphorylation regulator in the developmental stage of human testis ",
                            row_title= "",
                            # col = colorRamp2(c(-0.9,-0.45,0,0.45,0.9), c("#2c7bb6","#abd9e9","floralwhite","#fdae61", "#d7191c"),space = "RGB"),
                            # col = colorRamp2(c(-2,-1,0,1,2), c("#383899","#1E90FF","floralwhite","#FF8C69", "#FF0000"),space = "RGB"),
                            # col = colorRamp2(c(-1,-0.5,0,0.5,1), c("#383899","#1E90FF","floralwhite","#FF8C69", "#FF0000"),space = "RGB"),
                            col = colorRamp2(c(-1,-0.5,0,0.5,1), c("#2c7bb6","#abd9e9","White","#FF8C69", "#FF0000"),space = "RGB"),
                            clustering_distance_cols = "pearson",
                            cluster_rows = FALSE,
                            cluster_cols = FALSE,
                            row_dend_side = "left", #Dendogram on the right side
                            annotation_row = annotation_row,
                            annotation_colors = ann_colors,
                            fontsize = 10,
                            show_rownames = FALSE,
                            split = annotation_row$Group, #Splitting by Class
)
p

# 6. Violin Plot
########################################
colnames(data_testis)
data_testis<-data_testis[,-4]
datah7<- melt(data_testis,c("gene_symbol","cluster","Type"))
colnames(datah7)[2]<-"Group"
datah8<-datah7 %>% mutate(Group=recode(Group,"1"="group1","2"="group2",
                                       "3"="group3","4"="group4","5"="group5","6"="group6"))

p<-ggplot(datah8, aes(variable, value)) + 
  geom_violin(aes(fill = variable),scale="width",draw_quantiles = c(0.5))+
  # geom_violin(aes(fill = variable))+
  #stat_summary(fun.y=mean, geom="point", shape=19, size=1)+
  geom_jitter(height = 0, width = 0.1,size=0.2)+
  facet_wrap(~Group,ncol = 2)+
  # scale_fill_manual(values=c("#fc5185","#ff9e74","#e3c4a8","#3fc1c9","#364f6b"))+
  xlab("")+
  ylab("z-scale expression")+
  theme(panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold",size = rel(1.5)),
        strip.background = element_rect(fill = "lightblue",
                                        colour = "black",size = 1),
        legend.position = "none",
        axis.text.x=element_text(colour="black",angle=45,vjust = 1,hjust = 1,face = "bold",size = 12),
        # axis.text = element_text(angle = 30,hjust = 1,vjust = 1),
        axis.text.y=element_text(colour="black",face = "bold",size = 12),
        axis.title = element_text(angle = 30,hjust = 1,vjust = 1,colour="black",face = "bold",size = 14),
        axis.title.x = element_text(),
        axis.ticks=element_line(colour="black"),
        plot.margin=unit(c(1,1,1,1),"line"))

p


p<-ggplot(datah8, aes(variable, value)) + 
  geom_violin(aes(fill = variable),scale="width")+
  facet_wrap(~Group,ncol = 2)+
  ylab("zscore FPKM")+
  theme(panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold",size = rel(1.5)),
        strip.background = element_rect(fill = "lightblue",
                                        colour = "black",size = 1),
        legend.position = "none",
        axis.text.x=element_text(colour="black",angle=45,vjust = 1,hjust = 1,face = "bold",size = 12),
        axis.text.y=element_text(colour="black",size = 12),
        axis.title = element_text(colour="black",size = 14),
        axis.title.x = element_text(),
        axis.ticks=element_line(colour="black"),
        plot.margin=unit(c(1,1,1,1),"line"))
p


