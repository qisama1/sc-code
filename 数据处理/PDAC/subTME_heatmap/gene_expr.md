# gene expr

```R
scRNA = qread("/public/home/yuwenqi/sc-data/selected/35/all3.qs")
gene = read.csv("/public/home/yuwenqi/sc-data/selected/subTME_heatmap/gene_needed.csv")
write.csv(as.matrix(scRNA[gene$gene, ]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/subTME_heatmap/PDAC/gene_df.csv")
```