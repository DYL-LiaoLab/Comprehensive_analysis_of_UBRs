# CNV plotting in pan-cancer.
# Long deyu
# 2022.08.31

library(dplyr)
library(ComplexHeatmap)
library(stringr)
library(circlize)
library(ggplot2)
library(reshape2)
library(circlize)


cnv_gain<-read.csv("cnv_gain.txt",sep = "\t")
cnv_loss<-read.csv("cnv_loss.txt",sep = "\t")

cnv_gain$sum<-apply(cnv_gain,1,sum)
cnv_loss$sum<-apply(cnv_loss,1,sum)

cnv_gain1<-cnv_gain[order(cnv_gain$sum,decreasing = TRUE),]
cnv_loss1<-cnv_loss[order(cnv_loss$sum,decreasing = TRUE),]
rownames(cnv_gain1)[1:10]
rownames(cnv_loss1)[1:10]

p1<-ComplexHeatmap::pheatmap(as.matrix(cnv_gain1[,-34]*100),
                             clustering_distance_rows = "pearson",
                             # clustering_distance_rows ="spearman",
                             # clustering_distance_rows ="euclidean",
                             cluster_rows = FALSE,
                             cellwidth = 18,
                             show_rownames = FALSE,
                             # col = colorRamp2(c(0,10,20,30,40), c("floralwhite","#abd9e9","#2c7bb6","#fdae61", "#d7191c"),space = "RGB"),
                             col = colorRamp2(c(0,7,14,21,28), c("floralwhite","LightSkyBlue","DeepSkyBlue","LightSalmon1", "OrangeRed1"),space = "RGB"),
                             name = "CNV Gain",
                             # split=annotation_row$Type
)
p1
p2<-ComplexHeatmap::pheatmap(as.matrix(cnv_loss1[,-34]*100),
                             clustering_distance_rows = "pearson",
                             # clustering_distance_rows ="spearman",
                             # clustering_distance_rows ="euclidean",
                             # annotation_row = anno_row1,
                             # annotation_colors = ann_colors,
                             cluster_rows = FALSE,
                             cellwidth = 18,
                             # cellheight = 2,
                             show_rownames = FALSE,
                             # col = colorRamp2(c(0,10,20,30,40), c("floralwhite","#abd9e9","#2c7bb6","#fdae61", "#d7191c"),space = "RGB"),
                             # col = colorRamp2(c(0,10,20,30,40), c("floralwhite","Coral1","Tomato1","OrangeRed1", "Red1"),space = "RGB"),
                             col = colorRamp2(c(0,7,14,21,28), c("floralwhite","LightSkyBlue","DeepSkyBlue","LightSalmon1", "OrangeRed1"),space = "RGB"),
                             name = "CNV Loss",
                             # split=annotation_row$Type
)
pdf(file = "CNV_gain_all.pdf",height = 8,width = 8)
p1
dev.off()
png(filename = "CNV_gain_all.png",width = 800,height = 800)
p1
dev.off()
pdf(file = "CNV_loss_all.pdf",height = 8,width = 8)
p2
dev.off()
png(filename = "CNV_loss_all.png",width = 800,height = 800)
p2
dev.off()

cnv_gain$type<-"gain"
cnv_gain$ensembl<-rownames(cnv_gain)
cnv_gain<-cnv_gain[,-34] 
gain_boxplot<-melt(cnv_gain,id.vars = c("type","ensembl")) 

colnames(gain_boxplot)

ggplot(gain_boxplot,aes(x=variable,y=value,fill=variable))+
  geom_violin(scale = "width",adjust=1,trim = TRUE)+
  geom_boxplot(width=0.1,outlier.colour = NA)+
  theme_bw()

ggplot(gain_boxplot,aes(x=variable,y=value))+
  geom_boxplot()


cnv_loss<-cnv_loss[,-34]
cnv_loss$type<-"loss"
cnv_loss$ensembl<-rownames(cnv_loss)
loss_boxplot<-melt(cnv_loss,id.vars = c("type","ensembl"))

ggplot(loss_boxplot,aes(x=variable,y=value,fill=variable))+
  geom_violin(scale = "width",adjust=1,trim = TRUE)+
  geom_boxplot(width=0.1,outlier.colour = NA)+
  theme_bw()

library(ggsignif)
library(ggpubr)
gain_loss_boxplot<-rbind(gain_boxplot,loss_boxplot) 
gain_loss_boxplot1<-na.omit(gain_loss_boxplot)
colnames(gain_loss_boxplot)

shunxu<-as.vector(unique(gain_loss_boxplot$variable))
shunxu<-shunxu[order(shunxu)]
p<-ggplot(data = gain_loss_boxplot,aes(x=variable,y=value,fill=type))+
  geom_violin(scale = "width",trim = TRUE,adjust=1,size=0.01)+
  geom_boxplot(width=0.2,position = position_dodge(0.9),outlier.colour = NA,size=0.01)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values = c("Coral1","DarkSlateGray2"))+
  # scale_y_continuous(expand = c(0,0))+
  stat_compare_means(aes(group=type),
                     method = "t.test",label = "p.signif")+
  scale_y_reverse()+scale_x_discrete(limits=shunxu)+
  coord_cartesian(ylim = c(0.4,0))
p

p<-ggplot(data = gain_loss_boxplot,aes(x=variable,y=value,fill=type))+
  geom_violin(scale = "width",trim = TRUE,adjust=1,size=0.1)+
  geom_boxplot(width=0.2,position = position_dodge(0.9),outlier.colour = NA,size=0.1)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values = c("Coral1","DarkSlateGray2"))+
  coord_cartesian(ylim = c(0,0.6))+scale_x_discrete(limits=shunxu)+
  facet_grid(type~.)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 1))


compaired<-list(c("gain","loss"))
p<-ggplot(data = gain_loss_boxplot,aes(x=type,y=value,fill=type))+
  geom_violin(scale = "width",trim = TRUE,adjust=1,size=0.01)+
  geom_boxplot(width=0.2,position = position_dodge(0.9),outlier.colour = NA,size=0.01)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values = c("Coral1","DarkSlateGray2"))+
  # scale_y_continuous(expand = c(0,0))+
  # stat_compare_means(aes(group=type),
  #                    method = "t.test",label = "p.signif")+
  geom_signif(comparisons = compaired,test = wilcox.test,vjust = 1.5)+
  # scale_y_reverse()+
  # coord_cartesian(ylim = c(0.6,0))+
  facet_wrap(~variable,ncol=9)
p 
