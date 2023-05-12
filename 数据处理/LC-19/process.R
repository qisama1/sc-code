library(Seurat)
library(qs)
b = readRDS("/public/home/yuwenqi/sc-data/selected/19/bcell.rds")
tnk = readRDS("/public/home/yuwenqi/sc-data/selected/19/tnk.rds")
mye = readRDS("/public/home/yuwenqi/sc-data/selected/19/myeloid.rds")
fib = readRDS("/public/home/yuwenqi/sc-data/selected/19/fib.rds")
epi = readRDS("/public/home/yuwenqi/sc-data/selected/19/epi.rds")
end = readRDS("/public/home/yuwenqi/sc-data/selected/19/end.rds")




new.cluster.ids <- c("Follicular_B_1", "Follicular_B_2", "Follicular_B_3", "Follicular_B_4", "MALT_B_1", "Follicular_B_5", "Follicular_B_6",
"Follicular_B_7", "MALT_B_2", "GC_B_1",
"Plasma_1", "del", "GC_B_2", "Plasma_2", "del", "del", "MALT_B_3")
names(new.cluster.ids) <- levels(b)
b = RenameIdents(b, new.cluster.ids)
b[['ident']] = Idents(b)
b = subset(b, ident != 'del')
qsave(b, "/public/home/yuwenqi/sc-data/selected/19/bcell.qs")

new.cluster.ids <- c("CD4Tem_1", "CD4Tcm_1", "CD4Tcm_2", "NK_1", "Th", "CD8Tem_1", "Treg", "CD4+Tn", "del", "NK_2", "CD8Teff_1", "del", "CD8Tex",
"Tfh", "NK_3", "ISG15+T", "CD8Teff_2", "CD8Tem_2", "CD8Tem_3", "CD8Tem_4", "del", "CD8Trm_1", "CD8Tpo", "del", "CD8Trm_2", "NK_4", "del", "del", "del", "del")

names(new.cluster.ids) <- levels(tnk)
tnk = RenameIdents(tnk, new.cluster.ids)
tnk[['ident']] = Idents(tnk)
tnk = subset(tnk, ident != 'del')
qsave(tnk, "/public/home/yuwenqi/sc-data/selected/19/tnk.qs")

new.cluster.ids <- c("Macro_alve_1", "Macro_FOLR2", "Macro_alve_2", "Mono_1", "Macro_APOE", "Macro_SPP1", "cDC2_1", "Macro_FTL", "Macro_cycling",
 "Mono_1", "Macro_alve_3", "del", "cDC2_2", "Macro_APOE_SPP1", 
"pDC1", "Macro_CXCL10", "mature_cDC_1", "cDC2_3", "mature_cDC_2", "cDC2_4", "Macro_alve_4")

names(new.cluster.ids) <- levels(mye)
mye = RenameIdents(mye, new.cluster.ids)
mye[['ident']] = Idents(mye)
mye = subset(mye, ident != 'del')
qsave(mye, "/public/home/yuwenqi/sc-data/selected/19/mye.qs")


new.cluster.ids <- c("wound_myCAF", "Parenchymal_CAF", "ecm_myCAF", "IFNy_iCAF", "del", "Pericytes", "SMC_1", "del", "del", "del", "Perivascular_CAF", "SMC_2", 
"del", "del", "del", "del")
names(new.cluster.ids) <- levels(fib)
fib = RenameIdents(fib, new.cluster.ids)
fib[['ident']] = Idents(fib)
fib = subset(fib, ident != 'del')
qsave(fib, "/public/home/yuwenqi/sc-data/selected/19/fib.qs")

new.cluster.ids <- c("Squamous_1", "Basal_1", "AT2_1", "Basal_2", "AT2_2", "del", "pEMT_1", "AT2_3", "Cycling",
"Club_cell_1", "pEMT_2", "Stress", "del", "Club_cell_2", "Cilliated_cells", "ISG15+Epi", "Interferon", "astrocyte", "Squamous_2", "AT2_4",
"AT2_5", "Mesenchymal_1", "AT1", "AT2_6", "Glandular", "Mesenchymal_2", "AT2_7", "Basal_3", "Glandular", "del", "del", "del")

names(new.cluster.ids) <- levels(epi)
epi = RenameIdents(epi, new.cluster.ids)
epi[['ident']] = Idents(epi)
epi = subset(epi, ident != 'del')
qsave(epi, "/public/home/yuwenqi/sc-data/selected/19/epi.qs")

new.cluster.ids <- c("capillary_1", "lymphatic", "capillary_2", "activated_PCV_1", "arterial", "angiogenic_1", "activated_PCV_2", "del", "capillary_3",
"del", "angiogenic_2", "del", "angiogenic_3")

names(new.cluster.ids) <- levels(end)
end = RenameIdents(end, new.cluster.ids)
end[['ident']] = Idents(end)
end = subset(end, ident != 'del')
qsave(end, "/public/home/yuwenqi/sc-data/selected/19/end.qs")

cellinfo <- subset(mye@meta.data, select = c("Sample", "Sample_Origin", 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/19/mye_info.csv")
cellinfo <- subset(fib@meta.data, select = c("Sample", "Sample_Origin", 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/19/fib_info.csv")
cellinfo <- subset(end@meta.data, select = c("Sample", "Sample_Origin", 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/19/end_info.csv")
cellinfo <- subset(tnk@meta.data, select = c("Sample", "Sample_Origin", 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/19/tnk_info.csv")
cellinfo <- subset(epi@meta.data, select = c("Sample", "Sample_Origin", 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/19/epi_info.csv")
cellinfo <- subset(b@meta.data, select = c("Sample", "Sample_Origin", 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/19/b_info.csv")

mye = qread("/public/home/yuwenqi/sc-data/selected/19/mye.qs")
end = qread("/public/home/yuwenqi/sc-data/selected/19/end.qs")
fib = qread("/public/home/yuwenqi/sc-data/selected/19/fib.qs")
tnk = qread("/public/home/yuwenqi/sc-data/selected/19/tnk.qs")
b = qread("/public/home/yuwenqi/sc-data/selected/19/bcell.qs")
epi = qread("/public/home/yuwenqi/sc-data/selected/19/epi.qs")

scRNA = merge(mye, end)
scRNA = merge(scRNA, epi)
scRNA = merge(scRNA, fib)
scRNA = merge(scRNA, b)
scRNA = merge(scRNA, tnk)

scRNA <- NormalizeData(scRNA) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData()
scRNA <- RunPCA(scRNA, verbose = F)
pc.num=1:30
scRNA <- scRNA %>% RunUMAP(dims=pc.num)
qsave(scRNA, "/public/home/yuwenqi/sc-data/selected/19/all.qs")