library(Seurat)
library(qs)
library(tidyverse)
library(ggsci)
library(ggplot2)
library(cowplot)
library(data.table)
library(ggplot2)
library(ggprism)
library(RColorBrewer)
library(ggpubr)
library(ggplot2)
library(cowplot)

## 载入数据
Epi = qread("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/Epi/epi.qs")
Mye = qread("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/Mye/mye.qs")
TNK = qread("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/TNK/tnk.qs")
Fib = qread("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/Fib/fib.qs")
End = qread("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/End/end.qs")

Epi = subset(Epi, dataset == 'dataset3')
Mye = subset(Mye, dataset == 'dataset3')
TNK = subset(TNK, dataset == 'dataset3')
Fib = subset(Fib, dataset == 'dataset3')
End = subset(End, dataset == 'dataset3')


# Angiogenic_EC
Angiogenic_EC = c('INSR','ITGA5','ICAM1', 'ICAM2', 'HLA-DPB1', 'CD40', 'SEMA4C', 'LAMB1', 'COL4A1', 'JAG1', 'JAG2', 'PDGFB', 'CCL14', 'DLL1', 'CXCL12'
,'CD34', 'NECTIN2', 'CD44', 'SELP', 'SELE', 'EPHA4', 'EPHA2', 'FN1', 'JAM2', 'HSPG2', 'GRN', 'NOTCH1', 'DLL4', 'CD36', 'TNFRSF1A', 'TNFRSF1B',
'NOTCH3', 'APP', 'PECAM1')

plot_data = cbind(End@meta.data[,'Tissue', drop=F], t(as.matrix(End[Angiogenic_EC,]@assays$RNA@data)))
colnames(plot_data) = colnames(plot_data) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

## 画箱线图
## 5个一组
p1 = ggplot(plot_data,aes(x=Tissue,y=TNFRSF1B),fill=Tissue)+
  scale_fill_manual(values=c('#20B2AA', '#FF1493'))+
  geom_boxplot(aes(fill=Tissue),show.legend = F)+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
  labs(x=NULL,
       y=NULL,
       title = "TNFRSF1B")+
  #theme_minimal() 
  #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
  coord_cartesian(ylim = c(0, 5))+
  stat_compare_means(symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif", method='wilcox.test', label.y = 5, label.x = 0.5)
p2 = ggplot(plot_data,aes(x=Tissue,y=NOTCH3),fill=Tissue)+
  scale_fill_manual(values=c('#20B2AA', '#FF1493'))+
  geom_boxplot(aes(fill=Tissue),show.legend = F)+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
  labs(x=NULL,
       y=NULL,
       title = "NOTCH3")+
  #theme_minimal() 
  #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
  coord_cartesian(ylim = c(0, 5))+
  stat_compare_means(symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif", method='wilcox.test', label.y = 5, label.x = 0.5)

p3 = ggplot(plot_data,aes(x=Tissue,y=APP),fill=Tissue)+
  scale_fill_manual(values=c('#20B2AA', '#FF1493'))+
  geom_boxplot(aes(fill=Tissue),show.legend = F)+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
  labs(x=NULL,
       y=NULL,
       title = "APP")+
  #theme_minimal() 
  #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
  coord_cartesian(ylim = c(0, 5))+
  stat_compare_means(symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif", method='wilcox.test', label.y = 5, label.x = 0.5)

p4 = ggplot(plot_data,aes(x=Tissue,y=PECAM1),fill=Tissue)+
  scale_fill_manual(values=c('#20B2AA', '#FF1493'))+
  geom_boxplot(aes(fill=Tissue),show.legend = F)+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
  labs(x=NULL,
       y=NULL,
       title = "PECAM1")+
  #theme_minimal() 
  #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
  coord_cartesian(ylim = c(0, 5))+
  stat_compare_means(symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif", method='wilcox.test', label.y = 5, label.x = 0.5)

p5 = ggplot(plot_data,aes(x=Tissue,y=TNFRSF1A),fill=Tissue)+
  scale_fill_manual(values=c('#20B2AA', '#FF1493'))+
  geom_boxplot(aes(fill=Tissue),show.legend = F)+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
  labs(x=NULL,
       y=NULL,
       title = "TNFRSF1A")+
  #theme_minimal() 
  #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
  coord_cartesian(ylim = c(0, 5))+
  stat_compare_means(symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif", method='wilcox.test', label.y = 5, label.x = 0.5)

plot_grid(p1, p2, p3, p4, p5, nrow=1)
ggsave("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/diff_compare/Angiogenic/TNFRSF1B-NOTCH3-APP-PECAM1-TNFRSF1A.pdf", width = 5.2, height =3, device = cairo_pdf, units = 'in')