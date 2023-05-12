# get expr
```R
scRNA = qread("/public/home/yuwenqi/sc-data/selected/9/all.qs")
gene = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/corr_heatmap/caf_gene.csv")
write.csv(as.matrix(scRNA[gene$gene, ]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/append_ana/corr_heatmap/ESCA/gene_df.csv")
```