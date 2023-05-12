meta = read.csv("/public2022/chenzixi/research/TCGA_data/new-RNAseq/clin.filter.csv", row.names = 1)

tpm = read.csv("/public2022/chenzixi/research/TCGA_data/new-RNAseq/all-tpm_unstranded.csv", row.names=1)
id_name = read.csv("/public2022/chenzixi/research/TCGA_data/new-RNAseq/id_name.csv", row.names=1)
tpm = t(tpm)
colnames(tpm) = id_name[colnames(tpm), 'Gene_name']
tpm = t(tpm)

lihc = read.csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/others/lihc_meta2.csv", row.names=1)
lihc_counts = tpm[, colnames(lihc)]
coad = read.csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/others/coad_meta2.csv", row.names=1)
coad_counts = tpm[, colnames(coad)]
brca = read.csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/others/brca_meta2.csv", row.names=1)
brca_counts = tpm[, colnames(brca)]



corr_res1 = cor(lihc_counts, method = 'spearman')
corr_res2 = cor(brca_counts, method = 'spearman')
corr_res3 = cor(coad_counts, method = 'spearman')
write.csv(corr_res1, "xx.csv")


plot_cells(cds,
           color_cells_by = "cell_type")


ggsave(filename = paste0("/public/home/yuwenqi/sc-data/selected/9/module/cci_an_sheet/monocle/epi/t.pdf"), plot = p, device = 'pdf', width = 6, height = 6) 