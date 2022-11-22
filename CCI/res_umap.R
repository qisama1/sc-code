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
library(qs)

scRNA = qread("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/bc3_t.qs")
Epi = qread("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/Epi/epi.qs")
Mye = qread("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/Mye/mye.qs")
TNK = qread("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/TNK/tnk.qs")
Fib = qread("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/Fib/fib.qs")
End = qread("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/End/end.qs")

Epi = subset(Epi, Tissue == 'Tumor')
Mye = subset(Mye, Tissue == 'Tumor')
TNK = subset(TNK, Tissue == 'Tumor')
Fib = subset(Fib, Tissue == 'Tumor')
End = subset(End, Tissue == 'Tumor')

ligandSub = End
recvSub = TNK

p1 = FeaturePlot(scRNA, reduction = "umap", features = c('SELL', 'CD34'), label = T)
ggsave("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/cci_umap/CD8Teff/SELL-CD34.pdf", plot=p1, device='pdf', width=10, height=5)

p1 = FeaturePlot(ligandSub, reduction = "umap", features = c('TNFSF13B'), label = T, label.size = 4, repel = T)
ggsave("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/cci_umap/Angiogenic_EC/TNFSF13B_End.pdf", plot=p1, device='pdf', width=5, height=5)

p1 = FeaturePlot(recvSub, reduction = "umap", features = c('CD34'), label = T, label.size = 4, repel = T)
ggsave("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/cci_umap/CD8Teff/CD34_TNK.pdf", plot=p1, device='pdf', width=5, height=5)


# 补充画图
p1 = FeaturePlot(plot_clu, reduction = "umap", features = c('CXCR3', 'CD28'), label = T, label.size = 4, repel = T)
ggsave("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/cci_umap/sup/Treg.pdf", plot=p1, device='pdf', width=10, height=5)
