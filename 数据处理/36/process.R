b = readRDS("/public/home/yuwenqi/sc-data/selected/36/bcell.rds")
tnk = readRDS("/public/home/yuwenqi/sc-data/selected/36/tnk.rds")
mye = readRDS("/public/home/yuwenqi/sc-data/selected/36/myeloid.rds")
fib = readRDS("/public/home/yuwenqi/sc-data/selected/36/fib.rds")
epi = readRDS("/public/home/yuwenqi/sc-data/selected/36/epi.rds")



new.cluster.ids <- c("Follicular_B", "AIM2+B ", "Plasma_1", "Plasma_2", "Plasma_3", "del", "del", "del", "del", "del")
names(new.cluster.ids) <- levels(b)
b = RenameIdents(b, new.cluster.ids)
b[['ident']] = Idents(b)
b = subset(b, ident != 'del')
qsave(b, "/public/home/yuwenqi/sc-data/selected/36/bcell.qs")

new.cluster.ids <- c("CD4+Tn", "Th_1", "NK_1", "CD8+Tc", "CD8+Tem_1", "CD8+Tem_2", "CD8+Tem_3", "Treg", "NK_2", "del", "NK_3", "CD8+Tn", "Th_2", "CD8+Tem_4", "Th_3", "del", 
"NK_4", "NK_5", "NK_6")

names(new.cluster.ids) <- levels(tnk)
tnk = RenameIdents(tnk, new.cluster.ids)
tnk[['ident']] = Idents(tnk)
tnk = subset(tnk, ident != 'del')
qsave(tnk, "/public/home/yuwenqi/sc-data/selected/36/tnk.qs")

new.cluster.ids <- c("Granulocyte_1", "MMP9+TAM", "APOE+SPP1+TAM", "Monocyte_1", "Granulocyte_2", "Mast", "Mon-lineage", "cDC2", "del", "del", "del", "CXCL8+TAM", "del", "pDC1", "mature_cDC", 
"HSC", "APOE+TAM", "Granulocyte_3")

names(new.cluster.ids) <- levels(mye)
mye = RenameIdents(mye, new.cluster.ids)
mye[['ident']] = Idents(mye)
mye = subset(mye, ident != 'del')
qsave(mye, "/public/home/yuwenqi/sc-data/selected/36/mye.qs")


new.cluster.ids <- c("Fib_progenitors", "Fib_progenitors", "SMC_1", "Ecm_myCAF", "SMC_2", "SMC_3", "Wound_myCAF", "Epi_like_CAF", "Ecm_myCAF", "del", "Detox_iCAF", "del", 
"del", "Detox_iCAF")
names(new.cluster.ids) <- levels(fib)
fib = RenameIdents(fib, new.cluster.ids)
fib[['ident']] = Idents(fib)
fib = subset(fib, ident != 'del')
qsave(fib, "/public/home/yuwenqi/sc-data/selected/36/fib.qs")

new.cluster.ids <- c("CSC", "pEMT_Epi", "RPL+Epi", "Classical_Epi_1", "Acinar_cell_1", "Oxphos_Epi", "Cycling_Epi", "del", "Basal_Epi",
"del", "Endocrine_cell", "del", "Classical_Epi_2", "Acinar_cell_2")

names(new.cluster.ids) <- levels(epi)
epi = RenameIdents(epi, new.cluster.ids)
epi[['ident']] = Idents(epi)
epi = subset(epi, ident != 'del')
qsave(epi, "/public/home/yuwenqi/sc-data/selected/36/epi.qs")


cellinfo <- subset(mye@meta.data, select = c("orig.ident", "sample", 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/36/mye_info.csv")
cellinfo <- subset(fib@meta.data, select = c("orig.ident", "sample", 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/36/fib_info.csv")
cellinfo <- subset(tnk@meta.data, select = c("orig.ident", "sample", 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/36/tnk_info.csv")
cellinfo <- subset(epi@meta.data, select = c("orig.ident", "sample", 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/36/epi_info.csv")
cellinfo <- subset(b@meta.data, select = c("orig.ident", "sample", 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/36/b_info.csv")