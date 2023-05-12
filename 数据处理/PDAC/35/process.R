b = readRDS("/public/home/yuwenqi/sc-data/selected/35/bcell.rds")
tnk = readRDS("/public/home/yuwenqi/sc-data/selected/35/tnk.rds")
mye = readRDS("/public/home/yuwenqi/sc-data/selected/35/myeloid.rds")
fib = readRDS("/public/home/yuwenqi/sc-data/selected/35/fib.rds")
epi = readRDS("/public/home/yuwenqi/sc-data/selected/35/epi.rds")
end = readRDS("/public/home/yuwenqi/sc-data/selected/35/endo.rds")




new.cluster.ids <- c("Naïve_B", "Plasma_1", "GC B_1", "GC_B_2", "Follicular_B_cell", "del", "Plasma_2", "Plasma_3", "del")
names(new.cluster.ids) <- levels(b)
b = RenameIdents(b, new.cluster.ids)
b[['ident']] = Idents(b)
b = subset(b, ident != 'del')
qsave(b, "/public/home/yuwenqi/sc-data/selected/35/bcell.qs")

new.cluster.ids <- c("CD4+Tn", "CD8+Tpro", "Treg", "CD4+TIGIT+T", "CD8+Tem_1", "CD4+Tem_1", "Th", "CD8+Tem_2", "del", "del", "T_progenitors", "CD8+Tex")

names(new.cluster.ids) <- levels(tnk)
tnk = RenameIdents(tnk, new.cluster.ids)
tnk[['ident']] = Idents(tnk)
tnk = subset(tnk, ident != 'del')
qsave(tnk, "/public/home/yuwenqi/sc-data/selected/35/tnk.qs")

new.cluster.ids <- c("SPP1+TAM", "cDC2", "APOE+TAM", "Monocyte_1", "TRM", "Monocyte_2", "Cycing_TAM", "del", "CREM+TAM", "del", "ACP5+TAM", "Langerhans-like", "CXCL10+TAM", "del", 
"mature_cDC", "cDC1", "Mast", "APOE+SPP1+TAM", "pDC1", "del")

names(new.cluster.ids) <- levels(mye)
mye = RenameIdents(mye, new.cluster.ids)
mye[['ident']] = Idents(mye)
mye = subset(mye, ident != 'del')
qsave(mye, "/public/home/yuwenqi/sc-data/selected/35/mye.qs")


new.cluster.ids <- c("Detox_iCAF_1", "Ecm_myCAF_1", "Detox_iCAF_2", "Hypoxia-CAF", "Ecm_myCAF_2", "Detox_iCAF_3", "TGFβ_myCAF", "Wound_myCAF_1", "del", "Wound_myCAF_2", "IL_iCAF", "IFNγ_iCAF_1", 
"del", "IFNγ_iCAF_2", "Ecm_myCAF_3", "del", "del", "del", "del", "del")
names(new.cluster.ids) <- levels(fib)
fib = RenameIdents(fib, new.cluster.ids)
fib[['ident']] = Idents(fib)
fib = subset(fib, ident != 'del')
qsave(fib, "/public/home/yuwenqi/sc-data/selected/35/fib.qs")

new.cluster.ids <- c("Ductal_cell_A_1", "Ductal_cell_B_1", "Ductal_cell_A_2", "Ductal cell B_2", "Ductal_cell_B_3", "Ductal_cell_A_3", "Acinar_cell_1", "Ductal_cell A_4", "Interferon_Epi",
"Endocrine_cell", "Ductal_cell_B_4", "Acinar_cell_2", "Acinar_cell_3", "Basal", "Ductal_cell_A_5", "del", "Mesenchymal", "del", "del", "Ductal_cell_B_5")

names(new.cluster.ids) <- levels(epi)
epi = RenameIdents(epi, new.cluster.ids)
epi[['ident']] = Idents(epi)
epi = subset(epi, ident != 'del')
qsave(epi, "/public/home/yuwenqi/sc-data/selected/35/epi.qs")

new.cluster.ids <- c("Capillary_1", "PCV_1", "Angiogenic_1", "PCV_2", "Capillary_2", "Angiogenic_2", "Arterial_1", "Arterial_2", "del", "Venous", "Angiogenic_3", "del",
"Angiogenic_4", "del", "del", "del", "del")

names(new.cluster.ids) <- levels(end)
end = RenameIdents(end, new.cluster.ids)
end[['ident']] = Idents(end)
end = subset(end, ident != 'del')
qsave(end, "/public/home/yuwenqi/sc-data/selected/35/end.qs")

cellinfo <- subset(mye@meta.data, select = c("orig.ident", "sample", 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/35/mye_info.csv")
cellinfo <- subset(fib@meta.data, select = c("orig.ident", "sample", 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/35/fib_info.csv")
cellinfo <- subset(end@meta.data, select = c("orig.ident", "sample", 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/35/end_info.csv")
cellinfo <- subset(tnk@meta.data, select = c("orig.ident", "sample", 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/35/tnk_info.csv")
cellinfo <- subset(epi@meta.data, select = c("orig.ident", "sample", 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/35/epi_info.csv")
cellinfo <- subset(b@meta.data, select = c("orig.ident", "sample", 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/35/b_info.csv")