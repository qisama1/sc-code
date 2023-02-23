
new.cluster.ids <- c("Follicular B_1", "Memory_B", "Plasma_1", "GC_B_1", "GC_B_2", "del", "Stress_B", "del", "Plasma_2", "Plasma_3", "Naïve_B", "Plasma_4", "Cycling_B", "del")
names(new.cluster.ids) <- levels(b)
b = RenameIdents(b, new.cluster.ids)
b[['ident']] = Idents(b)
b = subset(b, ident != 'del')
b = subset(b, ident != 'Myeloid cells')

new.cluster.ids <- c("CD4+Tn", "CD8+Tem", "CD8+Tex", "del", "Th_1", "NK_1", "Treg_1", "CD4+Tcm", "Tstr", "Th_2", "CD8+Tcm", "NK_2", "del", "Tpro", "NK_3", "Treg_2", "del", "del", "del", "del")

names(new.cluster.ids) <- levels(tnk)
tnk = RenameIdents(tnk, new.cluster.ids)
tnk[['ident']] = Idents(tnk)
tnk = subset(tnk, ident != 'del')


new.cluster.ids <- c("Monocyte_1", "Mature_cDC_1", "APOE+TAM", "THBS1+TAM", "APOE+SPP1+TAM", "Stress_TAM", "cDC2", "del", "del", "pDC1", "del", "cDC1", "Mature_cDC_2", "CXCL10+TAM", 
"Mature_cDC_3", "Monocyte_2", "del", "Cycing_TAM", "del", "del")

names(new.cluster.ids) <- levels(mye)
mye = RenameIdents(mye, new.cluster.ids)
mye[['ident']] = Idents(mye)
mye = subset(mye, ident != 'del')

new.cluster.ids <- c("IFNγ-iCAF_1", "Detox_iCAF_1", "Detox_iCAF_2", "del", "Ecm-myCAF_1", "Ecm-myCAF_2", "Ecm-myCAF_3", "del", "del", "del", "Detox_iCAF_3", "microglia", 
"del", "del")
names(new.cluster.ids) <- levels(fib)
fib = RenameIdents(fib, new.cluster.ids)
fib[['ident']] = Idents(fib)
fib = subset(fib, ident != 'del')


new.cluster.ids <- c("MAL+Epi", "Mucosal_1", "KRT17+Epi", "Stress_Epi_1", "Cycling Epi", "Mucosal_2", "CXCL17+Epi", "Ap_Epi_1", "IGF2BP2+Epi", "Mes_Epi", "Mucosal_3", "Mucosal_4",
"Mucosal_5", "Stress_Epi_2", "Stress_Epi_3", "Ap_Epi_2")

names(new.cluster.ids) <- levels(epi)
epi = RenameIdents(epi, new.cluster.ids)
epi[['ident']] = Idents(epi)
epi = subset(epi, ident != 'del')

new.cluster.ids <- c("PCV", "Lymphatic", "Angiogenic_1", "Aactivated_PCV", "Capillary", "Angiogenic_2", "Arterial", "del", "Pericyte", "del", "del", "Angiogenic_3",
"del", "Angiogenic_4", "del", "del")

names(new.cluster.ids) <- levels(end)
end = RenameIdents(end, new.cluster.ids)
end[['ident']] = Idents(end)
end = subset(end, ident != 'del')


cellinfo <- subset(mye@meta.data, select = c("orig.ident", "sample", 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/10/mye_info.csv")
cellinfo <- subset(fib@meta.data, select = c("orig.ident", "sample", 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/10/fib_info.csv")
cellinfo <- subset(end@meta.data, select = c("orig.ident", "sample", 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/10/end_info.csv")
cellinfo <- subset(tnk@meta.data, select = c("orig.ident", "sample", 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/10/tnk_info.csv")
cellinfo <- subset(epi@meta.data, select = c("orig.ident", "sample", 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/10/epi_info.csv")
cellinfo <- subset(b@meta.data, select = c("orig.ident", "sample", 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/10/b_info.csv")