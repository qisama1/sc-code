scRNA = qread("/public/home/yuwenqi/sc-data/selected/9/all.qs")
gene = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/DEG-singlecell/ESCA/cellchat/gene_needed.csv")
write.csv(as.matrix(scRNA[gene$gene, ]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/append_ana/DEG-singlecell/ESCA/cellchat/gene_df.csv")

gene = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/DEG-singlecell/ESCA/cellphonedb/gene_needed.csv")
write.csv(as.matrix(scRNA[gene$gene, ]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/append_ana/DEG-singlecell/ESCA/cellphonedb/gene_df.csv")
