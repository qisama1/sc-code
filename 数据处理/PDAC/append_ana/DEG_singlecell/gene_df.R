scRNA = qread("/public/home/yuwenqi/sc-data/selected/35/all3.qs")
gene = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/DEG-singlecell/PDAC/cellchat/gene_needed.csv")
write.csv(as.matrix(scRNA[gene$gene, ]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/append_ana/DEG-singlecell/PDAC/cellchat/gene_df.csv")

gene = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/DEG-singlecell/PDAC/cellphonedb/gene_needed.csv")
write.csv(as.matrix(scRNA[gene$gene, ]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/append_ana/DEG-singlecell/PDAC/cellphonedb/gene_df.csv")
