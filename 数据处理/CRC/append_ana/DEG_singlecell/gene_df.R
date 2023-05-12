scRNA = qread("/public/home/yuwenqi/sc-data/selected/CRC/workspace/CRC1/crc1_used.qs")
gene = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/DEG-singlecell/CRC/cellchat/gene_needed.csv")
write.csv(as.matrix(scRNA[gene$gene, ]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/append_ana/DEG-singlecell/CRC/cellchat/gene_df.csv")

gene = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/DEG-singlecell/CRC/cellphonedb/gene_needed.csv")
write.csv(as.matrix(scRNA[gene$gene, ]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/append_ana/DEG-singlecell/CRC/cellphonedb/gene_df.csv")
