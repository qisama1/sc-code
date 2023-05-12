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
library(tidyverse)
library(Seurat)
library(qs)


scRNA = qread("/public/home/yuwenqi/sc-data/selected/35/all3.qs")
meta = read.csv("/public/home/yuwenqi/sc-data/selected/35/plot_use.csv", row.names=1)

scRNA[['cluster']] =  meta['cluster']
scRNA[['type']] = meta['type']
scRNA$type <- factor(scRNA$type,levels = c("Tumor","Normal"))
p1 <- DimPlot(scRNA, label = T, pt.size = 1, group.by='cluster')+
  NoLegend()+labs(title = "Cluster") +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        #plot.title = element_blank()
  ) + scale_color_npg()

p2 <- DimPlot(scRNA, label = T, pt.size = 1, group.by='type')+
  NoLegend()+labs(title = "Tumor/Normal") +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        #plot.title = element_blank()
        legend.position = c(.75, .10),
        legend.background = element_rect(fill = "transparent")
  ) + scale_color_npg()

ggsave(filename = "/public/home/yuwenqi/sc-data/selected/35/pdac_major_umap.pdf", plot = p1 | p2, device = 'pdf', width = 8, height = 4)