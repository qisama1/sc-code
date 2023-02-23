
tnk = readRDS("/public/home/yuwenqi/sc-data/selected/37/tnk.rds")
mye = readRDS("/public/home/yuwenqi/sc-data/selected/37/myeloid.rds")
fib = readRDS("/public/home/yuwenqi/sc-data/selected/37/fib.rds")
epi = readRDS("/public/home/yuwenqi/sc-data/selected/37/epi.rds")


new.cluster.ids <- c("Th", "del", "CD8+Tem", "Treg", "CD8+Tpro", "del", "del")

names(new.cluster.ids) <- levels(tnk)
tnk = RenameIdents(tnk, new.cluster.ids)
tnk[['ident']] = Idents(tnk)
tnk = subset(tnk, ident != 'del')
qsave(tnk, "/public/home/yuwenqi/sc-data/selected/37/tnk.qs")

new.cluster.ids <- c("FOLR2+TAM", "APOE+SPP1+TAM", "SLAMF8+TAM", "cDC2", "Monocyte_1", "Monocyte_2", "cDC1", "PCOLCE2+TAM", "del", "Cycing_TAM", "Monocyte_3")

names(new.cluster.ids) <- levels(mye)
mye = RenameIdents(mye, new.cluster.ids)
mye[['ident']] = Idents(mye)
mye = subset(mye, ident != 'del')
qsave(mye, "/public/home/yuwenqi/sc-data/selected/37/mye.qs")


new.cluster.ids <- c("Ecm_myCAF_1", "Ap_CAF", "IFNγ_iCAF_1", "del", "del", "del", "Epi_like_CAF", "Pericyte", "Ecm_myCAF_2", "del", "IFNγ_iCAF_2")
names(new.cluster.ids) <- levels(fib)
fib = RenameIdents(fib, new.cluster.ids)
fib[['ident']] = Idents(fib)
fib = subset(fib, ident != 'del')
qsave(fib, "/public/home/yuwenqi/sc-data/selected/37/fib.qs")

new.cluster.ids <- c("Acinar_cell_1", "Basal_1", "Classical", "REG4+Epi", "Ductal_cell_A", "Basal_2", "Mesenchymal_1", "Oxphos_Epi", "Ductal_cell_B",
"Mesenchymal_2", "Stress", "Ap_Epi", "Acinar_cell_2", "del", "del")

names(new.cluster.ids) <- levels(epi)
epi = RenameIdents(epi, new.cluster.ids)
epi[['ident']] = Idents(epi)
epi = subset(epi, ident != 'del')
qsave(epi, "/public/home/yuwenqi/sc-data/selected/37/epi.qs")


cellinfo <- subset(mye@meta.data, select = c("orig.ident", "sample", 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/37/mye_info.csv")
cellinfo <- subset(fib@meta.data, select = c("orig.ident", "sample", 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/37/fib_info.csv")
cellinfo <- subset(tnk@meta.data, select = c("orig.ident", "sample", 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/37/tnk_info.csv")
cellinfo <- subset(epi@meta.data, select = c("orig.ident", "sample", 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/37/epi_info.csv")
