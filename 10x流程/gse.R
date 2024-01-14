library(Seurat)
library(tidyverse)
gse1 = readRDS("/public2022/humeiling/hcc/usedata/seurat/gse1.rds")
gse2 = readRDS("/public2022/humeiling/hcc/usedata/seurat/gse2.rds")
gse3 = readRDS("/public2022/humeiling/hcc/usedata/seurat/gse3.rds")
gse4 = readRDS("/public2022/humeiling/hcc/usedata/seurat/gse4.rds")

gse = merge(gse1, gse2)
gse = merge(gse, gse3)
gse = merge(gse, gse4)

scRNA = gse
maxGene=6000
pctMT=20
#scRNA=subset(scRNA, features=genes[0:19814, 'Gene_name'])

scRNA[["percent.mt"]]=PercentageFeatureSet(scRNA,pattern = "^MT-")#?????????????????????
scRNA = subset(scRNA, subset = nFeature_RNA < maxGene & percent.mt < pctMT)

scRNA <- NormalizeData(scRNA)
scRNA = FindVariableFeatures(scRNA, nfeatures = 3000)
scRNA = ScaleData(scRNA)
scRNA <- RunPCA(scRNA, verbose = F)
pc.num=1:30
scRNA <- RunUMAP(scRNA, dims=pc.num)
scRNA <- FindNeighbors(scRNA, dims=pc.num) %>% FindClusters(resolution=0.05)

setwd("/public/home/yuwenqi/works/humeiling/")
library(cowplot)
p3 <- DimPlot(scRNA, reduction = "umap", group.by = "orig.ident", pt.size=0.5)
p4 <- DimPlot(scRNA, reduction = "umap", group.by = "ident",   pt.size=0.5, label = TRUE,repel = TRUE)
fig_umap <- plot_grid(p3, p4, labels = c('orig.ident','ident'),align = "v",ncol = 2)

ggsave(filename = "./batch_plot.pdf", plot = fig_umap, device = 'pdf', width = 30, height = 15, units = 'cm')

saveRDS(scRNA, "/public/home/yuwenqi/works/humeiling/gse.rds")