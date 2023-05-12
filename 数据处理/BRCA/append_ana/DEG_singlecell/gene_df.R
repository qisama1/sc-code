scRNA = qread("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/bc3.qs")
gene = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/DEG-singlecell/BC/cellchat/gene_needed.csv")
write.csv(as.matrix(scRNA[gene$gene, ]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/append_ana/DEG-singlecell/BC/cellchat/gene_df.csv")

gene = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/DEG-singlecell/BC/cellphonedb/gene_needed.csv")
write.csv(as.matrix(scRNA[gene$gene, ]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/append_ana/DEG-singlecell/BC/cellphonedb/gene_df.csv")
