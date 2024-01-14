library(Seurat)
library(harmony)
library(tidyverse)
library(patchwork)

rm(list = ls())
fs=list.files('./','^GSM')  
fs
library(stringr)
samples=str_split(fs,'_',simplify = T)[,1] 

lapply(unique(samples),function(x){
  #x=unique(samples)[1]
  y=fs[grepl(x,fs)]
  folder=paste0("5/", str_split(y[1],'_',simplify = T)[,1])
  dir.create(folder,recursive = T)
  
  file.rename(paste0("./",y[1]),file.path(folder,"barcodes.tsv.gz"))
  file.rename(paste0("./",y[2]),file.path(folder,"features.tsv.gz"))
  file.rename(paste0("./",y[3]),file.path(folder,"matrix.mtx.gz"))
})
samples=list.files('./','^GSM')

sceList = lapply(samples,function(pro){
  folder=file.path("5/",pro)
  CreateSeuratObject(counts = Read10X(folder),  #?????????????????????????????????
                     project = pro ,min.cells=50,min.features=800)
})
for (i in 1:length(sceList)) {
  sceList[[i]][["percent.mt"]]=PercentageFeatureSet(sceList[[i]],pattern = "^MT-")
}

scRNA <- merge(sceList[[1]], sceList[2:length(sceList)])

sce <- Read10X_h5(filename = "GSM4107899_LH16.3814_raw_gene_bc_matrices_h5.h5")
sce <- CreateSeuratObject(counts = sce)




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

library(cowplot)
p3 <- DimPlot(scRNA, reduction = "umap", group.by = "orig.ident", pt.size=0.5)
p4 <- DimPlot(scRNA, reduction = "umap", group.by = "ident",   pt.size=0.5, label = TRUE,repel = TRUE)
fig_umap <- plot_grid(p3, p4, labels = c('orig.ident','ident'),align = "v",ncol = 2)
ggsave(filename = "./batch_plot.pdf", plot = fig_umap, device = 'pdf', width = 30, height = 15, units = 'cm')
### harmony
scRNA <- readRDS("scRNA.rds")
cellinfo <- subset(scRNA@meta.data, select = c("orig.ident", "sub_cluster", 'dataset', 'Sample_Origin'))
scRNA <- CreateSeuratObject(scRNA@assays$RNA@counts, meta.data = cellinfo, min.cells=50, min.features=600)
### SCT
scRNA <- SCTransform(scRNA)
pc.num=1:30
### PCA
scRNA <- RunPCA(scRNA, npcs=pc.num, verbose=FALSE)
scRNA <- RunHarmony(scRNA, group.by.vars="orig.ident", max.iter.harmony = 8) 
scRNA <- FindNeighbors(scRNA, dims = 1:30)
scRNA <- FindClusters(scRNA, resolution = 0.2)

scRNA <- RunUMAP(scRNA, reduction="harmony", dims=pc.num)
library(cowplot)
p3 <- DimPlot(scRNA, reduction = "umap", group.by = "orig.ident", pt.size=0.5)
p4 <- DimPlot(scRNA, reduction = "umap", group.by = "seurat_clusters",   pt.size=0.5, label = TRUE,repel = TRUE)
#p5 <- DimPlot(scRNA, reduction = "umap", group.by = "type",   pt.size=0.5, label = TRUE,repel = TRUE)
fig_umap <- plot_grid(p3, p4, labels = c('orig.ident','ident'),align = "v",ncol = 2)
ggsave(filename = "./scanpy_batch2.pdf", plot = fig_umap, device = 'pdf', width = 45, height = 15, units = 'cm')

### Marker
library(future)
plan("multiprocess", workers = 16)
options(future.globals.maxSize = 200*1024^3)
markers <- FindAllMarkers(scRNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers, file = "./markers_batched2.csv", sep = ",", quote = F ,row.names = T)