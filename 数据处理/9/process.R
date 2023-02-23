new.cluster.ids <- c("Stress_Epi", "Cycling_Epi", "RPLP1+Epi", 
"Cycling_Epi", "NDRG1+Epi", "Mucosal_Epi", "Ap_Epi_1", "Mucosal_Epi", "MAFB+Epi", "Ap_Epi_2", "Mes_Epi", "A2M+Epi",
"Ap_Epi_3", "del", "del")

names(new.cluster.ids) <- levels(epi)
epi = RenameIdents(epi, new.cluster.ids)
epi[['ident']] = Idents(epi)
epi = subset(epi, ident != 'del')

new.cluster.ids <- c("Ecm_myCAF_1", "Detox_iCAF", "Ecm_myCAF_2", 
"IL_iCAF", "Wound_myCAF", "IFNγ_iCAF", "Adventitial_CAF", "Cycing_CAF", "Vascular_smooth_muscle", "del", "del")

names(new.cluster.ids) <- levels(fib)
fib = RenameIdents(fib, new.cluster.ids)
fib[['ident']] = Idents(fib)
fib = subset(fib, ident != 'del')

new.cluster.ids <- c("Naïve_B", "Follicular_B", "Plasma_1", 
"GC_B", "Cycling_B", "Plasma_2", "del")

names(new.cluster.ids) <- levels(b)
b = RenameIdents(b, new.cluster.ids)
b[['ident']] = Idents(b)
b = subset(b, ident != 'del')

new.cluster.ids <- c("End_progenitor", "Aactivated_PCV", "Angiogenic_1", 
"Arterial", "Angiogenic_2", "PCV", "Cycing_End")
names(new.cluster.ids) <- levels(end)
end = RenameIdents(end, new.cluster.ids)
end[['ident']] = Idents(end)
end = subset(end, ident != 'del')

new.cluster.ids <- c("Treg", "CD8+Tex", "CD8+Tem_1", 
"CD4+Tn", "NK_1", "Tfh", "CD8+Tc_1", "Th", "CD8+Tpro_1", "CD8+Tpro_2", "CD8+Tc_2", "CD8+Tem_2",
"Treg_2", "Treg_3", "NK_2", "CD8+Tem_3", 'del')

names(new.cluster.ids) <- levels(tnk)
tnk = RenameIdents(tnk, new.cluster.ids)
tnk[['ident']] = Idents(tnk)
tnk = subset(tnk, ident != 'del')


new.cluster.ids <- c("SPP1+APOE+TAM", "Mast", "CXCL10+TAM", 
"cDC2", "Monocyte_1", "Monocyte_2", "Monocyte_3", "SPP1+TAM", "mature_cDC", "pDC1", "Cycing_TAM", "cDC1")

names(new.cluster.ids) <- levels(mye)
mye = RenameIdents(mye, new.cluster.ids)
mye[['ident']] = Idents(mye)
mye = subset(mye, ident != 'del')

cellinfo <- subset(mye@meta.data, select = c("orig.ident", "sample", 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/9/mye_info2.csv")
cellinfo <- subset(fib@meta.data, select = c("orig.ident", "sample", 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/9/fib_info.csv")
cellinfo <- subset(end@meta.data, select = c("orig.ident", "sample", 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/9/end_info2.csv")
cellinfo <- subset(tnk@meta.data, select = c("orig.ident", "sample", 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/9/tnk_info.csv")
cellinfo <- subset(epi@meta.data, select = c("orig.ident", "sample", 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/9/epi_info.csv")
cellinfo <- subset(b@meta.data, select = c("orig.ident", "sample", 'ident'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/9/b_info2.csv")