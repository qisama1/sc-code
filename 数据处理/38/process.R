new.cluster.ids <- c("Plasma", "B_cell", "del", "del", "del", "del")
names(new.cluster.ids) <- levels(b)
b = RenameIdents(b, new.cluster.ids)
b[['ident']] = Idents(b)
b = subset(b, ident != 'del')
qsave(b, "/public/home/yuwenqi/sc-data/selected/38/bcell.qs")

new.cluster.ids <- c("CD4+Tn", "CD8+Tem", "CD4+Tem", "NK_1", "Treg", "CD8+Tpro", "CD8+TRM", "NK_2", "ILC", "CD8+Tex", "NKT", "del")

names(new.cluster.ids) <- levels(tnk)
tnk = RenameIdents(tnk, new.cluster.ids)
tnk[['ident']] = Idents(tnk)
tnk = subset(tnk, ident != 'del')
qsave(tnk, "/public/home/yuwenqi/sc-data/selected/38/tnk.qs")

new.cluster.ids <- c("Monocyte_1", "APOE+SPP1+TAM_1", "cDC2_1", "del", "del", "FOLR2+TAM", "Monocyte_2", "TRM", "Cycing_TAM", "APOE+SPP1+TAM_2", "pDC1_1", "Monocyte_3", "Monocyte_4", "LYVE1+TAM", 
"del", "CXCL10+TAM", "cDC1_2", "pDC1_2", "mature_cDC")

names(new.cluster.ids) <- levels(mye)
mye = RenameIdents(mye, new.cluster.ids)
mye[['ident']] = Idents(mye)
mye = subset(mye, ident != 'del')
qsave(mye, "/public/home/yuwenqi/sc-data/selected/38/mye2.qs")

new.cluster.ids <- c("Stress_1", "Stress_2", "Ductal_cell_A", "Basal_1", "Basal_2", "del", "Endocrine_cell", "Mesenchymal")

names(new.cluster.ids) <- levels(epi)
epi = RenameIdents(epi, new.cluster.ids)
epi[['ident']] = Idents(epi)
epi = subset(epi, ident != 'del')
qsave(epi, "/public/home/yuwenqi/sc-data/selected/38/epi.qs")


cellinfo <- subset(mye@meta.data, select = c("orig.ident", "biosample_id", 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/38/mye_info2.csv")
cellinfo <- subset(tnk@meta.data, select = c("orig.ident", "biosample_id", 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/38/tnk_info.csv")
cellinfo <- subset(epi@meta.data, select = c("orig.ident", "biosample_id", 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/38/epi_info.csv")
cellinfo <- subset(b@meta.data, select = c("orig.ident", "biosample_id", 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/38/b_info.csv")