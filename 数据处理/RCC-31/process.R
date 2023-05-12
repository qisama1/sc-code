library(Seurat)
library(qs)
b = readRDS("/public/home/yuwenqi/sc-data/selected/31/bcell.rds")
tnk = readRDS("/public/home/yuwenqi/sc-data/selected/31/tnk.rds")
mye = readRDS("/public/home/yuwenqi/sc-data/selected/31/myeloid.rds")
fib = readRDS("/public/home/yuwenqi/sc-data/selected/31/fib.rds")
epi = readRDS("/public/home/yuwenqi/sc-data/selected/31/epi.rds")
end = readRDS("/public/home/yuwenqi/sc-data/selected/31/end.rds")


new.cluster.ids <- c("B_1", "Follicular_B", "GC_B_1", "B_2", "del", "Naive B_1", "Naive B_2",
"GC_B_2", "B_3")
names(new.cluster.ids) <- levels(b)
b = RenameIdents(b, new.cluster.ids)
b[['ident']] = Idents(b)
b = subset(b, ident != 'del')
qsave(b, "/public/home/yuwenqi/sc-data/selected/31/bcell.qs")

new.cluster.ids <- c("CD8Tex_1", "NK_1", "CD8Tex_2", "CD8Tcm", "CD8Tem", "CD4Tn", "del", "Th", "CD4Tcm", "CD8Tpro", "Treg", "CD8Trm",
  "NK_2", "ILC_1", "NK_3", "ILC_2", "NK_4", "ISG16+T", "NK_5", "del", "del")
names(new.cluster.ids) <- levels(tnk)
tnk = RenameIdents(tnk, new.cluster.ids)
tnk[['ident']] = Idents(tnk)
tnk = subset(tnk, ident != 'del')
qsave(tnk, "/public/home/yuwenqi/sc-data/selected/31/tnk.qs")

new.cluster.ids <- c("Mono_1", "Macro_VEGFA", "Macro_APOE_SPP1", "Mono_2", "Macro_APOE", "cDC2", "del", "Mono_3", "Macro_FOLR2",
 "Macro_CXCL10", "del", "Macro_Cycing", "del", "del", 
"del", "del", "del", "cDC1", "del", "pDC1", "del", "del", "del")

names(new.cluster.ids) <- levels(mye)
mye = RenameIdents(mye, new.cluster.ids)
mye[['ident']] = Idents(mye)
mye = subset(mye, ident != 'del')
qsave(mye, "/public/home/yuwenqi/sc-data/selected/31/mye.qs")


new.cluster.ids <- c("NR4A1+Fib", "wound_myCAF", "del", "SMC", "del", "IL_iCAF", "IFNy_iCAF", "Parenchymal_CAF", "wound_myCAF", "del", "Cycling_CAF", "del", "del")
names(new.cluster.ids) <- levels(fib)
fib = RenameIdents(fib, new.cluster.ids)
fib[['ident']] = Idents(fib)
fib = subset(fib, ident != 'del')
qsave(fib, "/public/home/yuwenqi/sc-data/selected/31/fib.qs")

new.cluster.ids <- c("TNNT1+Epi", "TNNT1+Epi", "Tubule", "del", "Interferon_1", "del", "del", "Interferon_2", "loop_of_Henle",
"Epi_Cycling", "del", "pEMT", "Intercalated_cell_1", "Intercalated_cell_2", "del", "Intercalated_cell_3", "Podocyte", "Glandular", "Intercalated_cell_4")

names(new.cluster.ids) <- levels(epi)
epi = RenameIdents(epi, new.cluster.ids)
epi[['ident']] = Idents(epi)
epi = subset(epi, ident != 'del')
qsave(epi, "/public/home/yuwenqi/sc-data/selected/31/epi.qs")

new.cluster.ids <- c("capillary_1", "capillary_2", "angiogenic", "del", "del", "activated_pcv", "capillary_3", "del", "capillary_4",
"del")

names(new.cluster.ids) <- levels(end)
end = RenameIdents(end, new.cluster.ids)
end[['ident']] = Idents(end)
end = subset(end, ident != 'del')
qsave(end, "/public/home/yuwenqi/sc-data/selected/31/end.qs")

cellinfo <- subset(mye@meta.data, select = c("Sample", "Sample2", 'type', 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/31/mye_info.csv")
cellinfo <- subset(fib@meta.data, select = c("Sample", "Sample2", 'type', 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/31/fib_info.csv")
cellinfo <- subset(end@meta.data, select = c("Sample", "Sample2", 'type', 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/31/end_info.csv")
cellinfo <- subset(tnk@meta.data, select = c("Sample", "Sample2", 'type', 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/31/tnk_info.csv")
cellinfo <- subset(epi@meta.data, select = c("Sample", "Sample2", 'type', 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/31/epi_info.csv")
cellinfo <- subset(b@meta.data, select = c("Sample", "Sample2", 'type', 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/31/b_info.csv")

scRNA = merge(mye, end)
scRNA = merge(scRNA, epi)
scRNA = merge(scRNA, fib)
scRNA = merge(scRNA, b)
scRNA = merge(scRNA, tnk)

scRNA <- NormalizeData(scRNA) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData()
scRNA <- RunPCA(scRNA, verbose = F)
pc.num=1:30
scRNA <- scRNA %>% RunUMAP(dims=pc.num)
qsave(scRNA, "/public/home/yuwenqi/sc-data/selected/31/all.qs")