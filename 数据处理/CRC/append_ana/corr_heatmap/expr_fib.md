# get expr
```R
scRNA = qread("/public/home/yuwenqi/sc-data/selected/CRC/workspace/CRC1/crc1_used.qs")
gene = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/corr_heatmap/caf_gene.csv")
write.csv(as.matrix(scRNA[gene$gene, ]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/append_ana/corr_heatmap/CRC/gene_df.csv")
```